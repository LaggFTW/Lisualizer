#include <stdio.h>
#include <string.h>

#include <math.h>
#include <unistd.h>

#ifdef __gnu_linux__
#include <sys/time.h>
#endif

#ifdef __WINDOWS__
#include <windows.h>
#endif

#if defined(__APPLE__) && defined(__MACH__)
#include <sys/time.h>
#endif

#include <sndfile.h>

#include <portaudio.h>

#include "pa_ringbuffer.h"

#include "memutils.h"

/**
#ifdef __gnu_linux__
#include <pa_linux_alsa.h>
#endif

#ifdef __WINDOWS__
#include <pa_win_wasapi.h>
#endif
**/

//TODO: Figure out why gcc/g++ -libmp3lame flag gives errors
//#include <lame/lame.h>

#include <fftw3.h>

#ifndef GLFW_INCLUDE_GLU
    #define GLFW_INCLUDE_GLU
#endif

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#define ERRCODE_SUCCESS 0
#define ERRCODE_EOF 1
#define ERRCODE_FAILURE 2
#define ERRCODE_UNSUPPORTED -1

#define MILLION_INT 1000000

#define WINDOW_SIZE 4096

#define MOV_AVG_WINDOW_SIZE 4096

#define FRAME_RATE 60

//speed of sound, in meters per second
#define SOUND_SPEED 343.59

//distance between left and right audio receivers (in this case, ears), scaled by a damping factor
#define HEAD_WIDTH 0.35

//maximum interaural level difference, as a normalized ratio
#define MAX_ILD 1.6

//flag: 1 for a hanning window, 0 for a rectangular window
#define HANNING 1

//average value of squared hanning window
#define HANN_SCALING_FACTOR 0.375

typedef struct timeval timeval;

//user data struct; (mostly) useless for now
typedef struct {
    int poi;
} testdata;

/**
Parameters for processing; global for now
Will reorganize structure of code later anyways, so no point in trying to make
this structured now
**/
long window_increment; //used for static window shifting, if needed; for now, dynamic window shifts are being used
double *window[2];
double hann[WINDOW_SIZE];
fftw_complex *dft_window[2];
fftw_complex *dft_corre;
fftw_plan plan[2];
fftw_plan plancorre;
float samples[WINDOW_SIZE*2];
int transform_size;
int cutoff_1600; //index of dft coefficient corresponding to 1600 Hz

//separation cutoff data
fftw_complex *dft_bin[2];
int *freq_cutoffs;
int num_bins, num_cutoffs;

//visualization data
double *source_sines; //sine of the source angles
double *power_spectra; //current power spectrum per bin (separated "sound")
double *pspec_mov_avg; //moving average of previously obtained power spectra per bin

//audio stream data
SNDFILE* inputfile;
PaUtilRingBuffer io_buffer; //interface between audio input and audio output
PaUtilRingBuffer sample_buffer; //interface between audio input and audio processing
float io_buf_arr[WINDOW_SIZE*2];
float smp_buf_arr[WINDOW_SIZE*2];
int samplerate;
/**
End of parameters
**/


/**
GLFW callbacks
**/
static void error_callback(int error, const char* description)
{
    fputs(description, stderr);
}
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
}

/*
 * Returns: the timeval's toal value represented in microseconds.
 */
static long tval_to_usec(const timeval t){
    return (t.tv_sec * MILLION_INT) + t.tv_usec;
}

/*
 * Takes the difference between the timevals end and start, storing the value
 * of the result in elapsed.
 */
static void tval_sub(const timeval end, const timeval start, timeval *elapsed){
    elapsed->tv_sec = end.tv_sec - start.tv_sec;
    elapsed->tv_usec = end.tv_usec - start.tv_usec;
    /* check for carry */
    if (elapsed->tv_usec < 0){
        (elapsed->tv_sec)--;
        (elapsed->tv_usec) += MILLION_INT;
    }
}

static int pa_callback(const void *inputBuffer, void *outputBuffer, 
        unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo* timeInfo,
        PaStreamCallbackFlags statusFlags, void *userData){
    testdata *data = (testdata *)userData;
    float *in = (float *)inputBuffer;
    float *out = (float *)outputBuffer;
    float smp[2];
    int num_samples = framesPerBuffer * 2;
    memset(out, 0, num_samples * sizeof(float));
    for (int i = 0; i < num_samples; i++){
        if (PaUtil_GetRingBufferReadAvailable(&io_buffer)){
            PaUtil_ReadRingBuffer(&io_buffer, &(out[i]), 1);
        }
    }
    for (int i = 0; i < framesPerBuffer; i++){
        if ((PaUtil_GetRingBufferWriteAvailable(&io_buffer) >= 2) 
                && (PaUtil_GetRingBufferWriteAvailable(&io_buffer) >= 2)){
            if (sf_read_float(inputfile, smp, 2) < 2){
                //for processing purposes, if we encounter EOF, treat it as silence
                smp[0] = 0.0f; smp[1] = 0.0f;
                PaUtil_WriteRingBuffer(&sample_buffer, smp, 2);
                continue;
            }
            PaUtil_WriteRingBuffer(&io_buffer, smp, 2);
            PaUtil_WriteRingBuffer(&sample_buffer, smp, 2);
        }
    }
    return paContinue;
}

static int calc_cutoffs_2ch(fftw_complex *data[2], int *cutoffs, int start, int end, int cut_ind, int num_cutoffs){
    //notes on the input: cutoffs is an array of (num_cutoffs * 2) - 1 elements, which will store the
    //cutoffs indices based on the weighted average of values in the two channels of data,
    //and thus will allow the partitioning of the data into (num_cutoffs * 2) groups
    //this method also assumes that num_cutoffs is a power of 2
    double numer, denom;
    int cut = 0;
    numer = 0; denom = 0;
    for (int ch = 0; ch < 2; ch++){
        for (int i = start; i < end; i++){
            double norm = data[ch][i][0] * data[ch][i][0] + data[ch][i][1] * data[ch][i][1];
            numer += i * norm;
            denom += norm;
        }
    }
    cut = (denom > 0.0 || denom < 0.0)? (int)(numer/denom) : ((start + end)/2);
    cutoffs[cut_ind] = cut;
    num_cutoffs /= 2;
    if (num_cutoffs >= 1){
        int cut_before, cut_after;
        cut_before = cut_ind - num_cutoffs;
        cut_after = cut_ind + num_cutoffs;
        calc_cutoffs_2ch(data, cutoffs, start, cut, cut_before, num_cutoffs);
        calc_cutoffs_2ch(data, cutoffs, cut, end, cut_after, num_cutoffs);
    }
    return ERRCODE_SUCCESS;
}

//updates the moving average of power spectrum peaks, using a simple moving average approximation
static void mov_avg_update(){
    for (int i = 0; i < num_bins; i++){
        pspec_mov_avg[i] += (power_spectra[i] - pspec_mov_avg[i])/MOV_AVG_WINDOW_SIZE;
    }
}

static int init(const char *fname){
    SF_INFO fileinfo;

    //information for separation
    num_bins = 4;
    num_cutoffs = num_bins - 1;
    freq_cutoffs = (int *)(ucalloc(num_cutoffs, sizeof(int)));
    
    //information for visualization
    source_sines = (double *)(ucalloc(num_bins, sizeof(double)));
    power_spectra = (double *)(ucalloc(num_bins, sizeof(double)));
    pspec_mov_avg = (double *)(ucalloc(num_bins, sizeof(double)));
    
    /**
    Intial opening of file
    **/
    inputfile = sf_open(fname, SFM_READ, &fileinfo);
    if (!inputfile){
        printf("sndfile: file reading error\n");
        return ERRCODE_FAILURE;
    }
    if (fileinfo.channels != 2){
        return ERRCODE_UNSUPPORTED;
    }
    samplerate = fileinfo.samplerate;
    for (int i = 0; i < WINDOW_SIZE; i+=2){ //zero-initialize sample arrays and buffers
        io_buf_arr[i] = 0.0f;
        smp_buf_arr[i] = 0.0f;
        samples[i] = 0.0f;
    }
    //initialize ringbuffers
    if (PaUtil_InitializeRingBuffer(&io_buffer, sizeof(float), WINDOW_SIZE*2, io_buf_arr)){
        printf("io ringbuffer error\n");
    }
    if (PaUtil_InitializeRingBuffer(&sample_buffer, sizeof(float), WINDOW_SIZE*2, smp_buf_arr)){
        printf("sample ringbuffer error\n");
    }

    /**
    Hanning Window Calculation
    **/
    if (HANNING) {
        double pi = 3.141592653589793238463;
        long size = WINDOW_SIZE - 1;
        for (int i = 0; i < WINDOW_SIZE; i++){
            hann[i] = 0.5 * (1 - cos(2 * pi * i / size));
        }
    }
    
    /**
    Sliding Window Initialization
    **/
    transform_size = WINDOW_SIZE/2 + 1;
    cutoff_1600 = (int)((3200.0/samplerate) * (WINDOW_SIZE/2)); //cutoff = 1600 / nyquist freq * max transform index
    window_increment = samplerate / FRAME_RATE * 2; //for static window shift, if needed
    for (int i = 0; i < 2; i++){
        window[i] = (double *)(fftw_malloc(sizeof(double)*WINDOW_SIZE));
        dft_window[i] = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex)*transform_size));
        plan[i] = fftw_plan_dft_r2c_1d(WINDOW_SIZE, window[i], dft_window[i], FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
        
        dft_bin[i] = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex)*transform_size));
    }
    dft_corre = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex)*transform_size));
    plancorre = fftw_plan_dft_c2r_1d(WINDOW_SIZE, dft_corre, window[0], FFTW_ESTIMATE);
    return ERRCODE_SUCCESS;
}

/**
Sliding Window Analysis
**/
static int analyze(){
    //timing for synchronization
    timeval before, after, elapsed;
    unsigned long usec_remaining;
    gettimeofday(&before, NULL);
    
    /**
    Processing of current buffer waveform
    Current method uses a dynamic window, where the size is determined by
    how many samples are ready to be read from the ring buffer
    **/
    long dyn_win_shift, orig_win_offset;
    int i;
    //dyn_win_shift: the number of samples (frames * channels) to shift the window
    //orig_win_offset: the index of the first new sample
    dyn_win_shift = PaUtil_GetRingBufferReadAvailable(&sample_buffer);
    dyn_win_shift -= dyn_win_shift % 2; //ensure the number of samples shifted is a multiple of the number of channels (2)
    orig_win_offset = WINDOW_SIZE * 2 - dyn_win_shift;
    //bounds checking (pray that the below if block never executes)
    if (orig_win_offset < 0){
        orig_win_offset = 0;
        dyn_win_shift = WINDOW_SIZE * 2;
    }
    //shift over the old samples by dyn_win_shift
    for (int j = 0; j < orig_win_offset; j+=2){
        samples[j] = samples[j+dyn_win_shift];
        samples[j+1] = samples[j+1+dyn_win_shift];
    }
    //copy dyn_win_shift samples over to the space left over from the shift
    PaUtil_ReadRingBuffer(&sample_buffer, &(samples[orig_win_offset]), dyn_win_shift);
    i = 0;
    //split the channels, applying a windowing function if needed
    for (int k = 0; k < WINDOW_SIZE; i+=2, k++){
        if (HANNING){
            window[0][k] = samples[i] * hann[k];
            window[1][k] = samples[i+1] * hann[k];
        } else {
            window[0][k] = samples[i];
            window[1][k] = samples[i+1];
        }
    }
    //fourier analysis (separation, localization)
    fftw_execute(plan[0]);
    fftw_execute(plan[1]);
    //separation of the DFTs
    calc_cutoffs_2ch(dft_window, freq_cutoffs, 0, transform_size, num_cutoffs/2, num_bins/2);
    for (int i = 0; i < num_bins; i++){
        for (int k = 0; k < 2; k++){
            for (int j = 0; j < transform_size; j++){
                if ((i == 0 || j > freq_cutoffs[i-1]) && (i == num_cutoffs || j <= freq_cutoffs[i])){
                    memcpy(dft_bin[k][j], dft_window[k][j], sizeof(fftw_complex));
                } else {
                    dft_bin[k][j][0] = 0;
                    dft_bin[k][j][1] = 0;
                }
            }
        }
        //TODO: Add interaural level difference checking for bins representing audio above ~1600 Hz instead of cross-correlation
        //localization (take cross-correlation of the left and right using the DFTs already computed, compute a source angle)
        //negative angle: to the left, positive angle: to the right
        for (int j = 0; j < transform_size; j++){
            dft_corre[j][0] = (dft_bin[0][j][0] * dft_bin[1][j][0] + dft_bin[0][j][1] * dft_bin[1][j][1])/WINDOW_SIZE;
            dft_corre[j][1] = (dft_bin[0][j][0] * dft_bin[1][j][1] - dft_bin[0][j][1] * dft_bin[1][j][0])/WINDOW_SIZE;
        }
        fftw_execute(plancorre);
        {
            long max_index = 0;
            long max_delay, delay_value;
            double angle, factor, ratio;
            factor = SOUND_SPEED / (samplerate * HEAD_WIDTH);
            //maximum possible delay such that the angle calculated will be 90 degrees; 
            //used to distinguish between delay caused by sound coming in from an angle 
            //and delay caused by echoes
            max_delay = (long)((1 / factor)) + 2;
            if (max_delay > WINDOW_SIZE){
                max_delay = WINDOW_SIZE;
            }
            for (long j = 1; j < max_delay; j++){
                if (window[0][j] > window[0][max_index]){
                    max_index = j;
                }
                if (window[0][WINDOW_SIZE - j] > window[0][max_index]){
                    max_index = WINDOW_SIZE - j;
                }
            }
            delay_value = ((max_index > (WINDOW_SIZE / 2))? WINDOW_SIZE - max_index : -max_index);
            ratio = delay_value * factor;
            if (ratio > 1.0){
                ratio = 1.0;
            }
            if (ratio < -1.0){
                ratio = -1.0;
            }
            angle = asin(ratio);
            //printf("bin %d sample offset: %ld\n", i, value);
            //printf("if the value is negative, then the left precedes the right by the absolute value\n");
            //printf("otherwise, the right precedes the left by the value shown\n");
            printf("%lf, ", angle);//DEBUG
        }
    }
    printf("\n");//DEBUG
    
    //timing for synchronization
    gettimeofday(&after, NULL);
    tval_sub(after, before, &elapsed);
    usec_remaining = MILLION_INT/FRAME_RATE - tval_to_usec(elapsed);
    if (usec_remaining > 0){
        //printf("leeway ~= %ld\n", usec_remaining);//DEBUG
        usleep(usec_remaining);
    } else {
        printf("Warning: processing is taking too long\n.");
    }
}

/**
Sliding Window Analysis, no synchronization
**/
static void analyze_async(){
    /**
    Processing of current buffer waveform
    Current method uses a dynamic window, where the size is determined by
    how many samples are ready to be read from the ring buffer
    **/
    long dyn_win_shift, orig_win_offset;
    double window_sq; //stores value of WINDOW_SIZE ^ 2 for optimized analysis
    int i;
    //dyn_win_shift: the number of samples (frames * channels) to shift the window
    //orig_win_offset: the index of the first new sample
    dyn_win_shift = PaUtil_GetRingBufferReadAvailable(&sample_buffer);
    dyn_win_shift -= dyn_win_shift % 2; //ensure the number of samples shifted is a multiple of the number of channels (2)
    orig_win_offset = WINDOW_SIZE * 2 - dyn_win_shift;
    //bounds checking (pray that the below if block never executes)
    if (orig_win_offset < 0){
        orig_win_offset = 0;
        dyn_win_shift = WINDOW_SIZE * 2;
    }
    //shift over the old samples by dyn_win_shift
    for (int j = 0; j < orig_win_offset; j+=2){
        samples[j] = samples[j+dyn_win_shift];
        samples[j+1] = samples[j+1+dyn_win_shift];
    }
    //copy dyn_win_shift samples over to the space left over from the shift
    PaUtil_ReadRingBuffer(&sample_buffer, &(samples[orig_win_offset]), dyn_win_shift);
    i = 0;
    //split the channels, applying a windowing function if needed
    for (int k = 0; k < WINDOW_SIZE; i+=2, k++){
        if (HANNING){
            window[0][k] = samples[i] * hann[k];
            window[1][k] = samples[i+1] * hann[k];
        } else {
            window[0][k] = samples[i];
            window[1][k] = samples[i+1];
        }
    }
    //fourier analysis (separation, localization)
    fftw_execute(plan[0]);
    fftw_execute(plan[1]);
    //separation of the DFTs
    calc_cutoffs_2ch(dft_window, freq_cutoffs, 0, transform_size, num_cutoffs/2, num_bins/2);
    for (i = 0; i < num_bins; i++){
        double diff, intensity; //values to be used for ILD
        int loc_method = 0; //determines whether to use ITD or ILD
        power_spectra[i] = 0.0; //reset the current power spectrum values
        for (int k = 0; k < 2; k++){
            for (int j = 0; j < transform_size; j++){
                if ((i == 0 || j > freq_cutoffs[i-1]) && (i == num_cutoffs || j <= freq_cutoffs[i])){
                    memcpy(dft_bin[k][j], dft_window[k][j], sizeof(fftw_complex));
                } else {
                    dft_bin[k][j][0] = 0;
                    dft_bin[k][j][1] = 0;
                }
            }
        }
        if (i > 0){
            if (freq_cutoffs[i-1] >= cutoff_1600){
                loc_method = 1;
                diff = 0.0;
                intensity = 0.0;
            } else {
                if (i < num_cutoffs && (freq_cutoffs[i-1] + freq_cutoffs[i]) >= (cutoff_1600 * 2)){
                    loc_method = 1;
                    diff = 0.0;
                    intensity = 0.0;
                }
            }
        }
        /**
        This loop does two things
        1. Computes peak power of this bin
        2. Uses a sound localization method to find a source angle (done in two ways):
        Audio below ~1600 Hz: take cross-correlation of the left and right using the DFTs already computed, compute a source angle
            negative angle: to the left, positive angle: to the right
        Audio above ~1600 Hz: checks interaural level difference between left and right bins
        **/
        window_sq = WINDOW_SIZE * WINDOW_SIZE;
        for (int j = 0; j < transform_size; j++){
            double norm_l, norm_r, norm_a;
            norm_l = (dft_bin[0][j][0] * dft_bin[0][j][0] + dft_bin[0][j][1] * dft_bin[0][j][1])/(HANN_SCALING_FACTOR * window_sq);
            norm_r = (dft_bin[1][j][0] * dft_bin[1][j][0] + dft_bin[1][j][1] * dft_bin[1][j][1])/(HANN_SCALING_FACTOR * window_sq);
            norm_a = (norm_l + norm_r) / 2;
            //compute peak power
            if (norm_a > power_spectra[i]){
                power_spectra[i] = norm_a;
            }
            //prepare for sound localization
            if (loc_method){
                //audio above ~1600 Hz: use ILD
                diff += norm_r - norm_l;
                intensity += norm_a;
            } else {
                //audio below ~1600 Hz: use ITD, prepare for cross-correlation
                dft_corre[j][0] = (dft_bin[0][j][0] * dft_bin[1][j][0] + dft_bin[0][j][1] * dft_bin[1][j][1])/window_sq;
                dft_corre[j][1] = (dft_bin[0][j][0] * dft_bin[1][j][1] - dft_bin[0][j][1] * dft_bin[1][j][0])/window_sq;
            }
        }
        if (loc_method){
            //ILD
            double ratio;
            diff = (intensity > 0.0)? diff/intensity : 0.0;
            ratio = diff / MAX_ILD;
            if (ratio > 1.0){
                ratio = 1.0;
            }
            if (ratio < -1.0){
                ratio = -1.0;
            }
            source_sines[i] = ratio;
        } else {
            //ITD w/ cross-correlation
            long max_index = 0;
            long max_delay, delay_value;
            double factor, ratio;
            
            fftw_execute(plancorre); //finish computing cross-correlation
            
            factor = SOUND_SPEED / (samplerate * HEAD_WIDTH);
            //maximum possible delay such that the angle calculated will be 90 degrees; 
            //used to distinguish between delay caused by sound coming in from an angle 
            //and delay caused by echoes
            max_delay = (long)((1 / factor)) + 2;
            if (max_delay > WINDOW_SIZE){
                max_delay = WINDOW_SIZE;
            }
            for (long j = 1; j < max_delay; j++){
                if (window[0][j] > window[0][max_index]){
                    max_index = j;
                }
                if (window[0][WINDOW_SIZE - j] > window[0][max_index]){
                    max_index = WINDOW_SIZE - j;
                }
            }
            delay_value = ((max_index > (WINDOW_SIZE / 2))? WINDOW_SIZE - max_index : -max_index);
            ratio = delay_value * factor;
            if (ratio > 1.0){
                ratio = 1.0;
            }
            if (ratio < -1.0){
                ratio = -1.0;
            }
            //angle = asin(ratio); //we don't actually need the angle for rendering
            source_sines[i] = ratio;
        }
    }
}

static int cleanup(){
    for (int i = 0; i < 2; i++){
        fftw_destroy_plan(plan[i]);
        fftw_free(window[i]);
        fftw_free(dft_window[i]);
        
        fftw_free(dft_bin[i]);
    }
    free(freq_cutoffs);
    free(source_sines);
    free(power_spectra);
    free(pspec_mov_avg);
    fftw_destroy_plan(plancorre);
    fftw_free(dft_corre);
    fftw_cleanup();
    if (sf_close(inputfile)){
        printf("sndfile: error closing file\n");
    }
}

int main(int argc, char **argv){
    if (argc >= 2){
        PaStream *stream;
        PaError err;
        testdata data;
        double time_increment;
        
        /** GLFW initialization **/
        GLFWwindow* window;
        glfwSetErrorCallback(error_callback);
        if (!glfwInit())
            exit(EXIT_FAILURE);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
        glfwWindowHint(GLFW_SAMPLES, 4);
        window = glfwCreateWindow(720, 720, "Poi", NULL, NULL);
        if (!window)
        {
            glfwTerminate();
            exit(EXIT_FAILURE);
        }
        glfwMakeContextCurrent(window);
        glfwSwapInterval(1);
        glfwSetKeyCallback(window, key_callback);
        glewExperimental = GL_TRUE;
        if (glewInit() != GLEW_OK)
            exit(EXIT_FAILURE);
        glEnable(GL_FRAMEBUFFER_SRGB_EXT);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glShadeModel(GL_SMOOTH);
        glEnable(GL_MULTISAMPLE);
        /** **/
        
        err = Pa_Initialize(); //initialize portaudio
        if(err != paNoError){
            return 0;
        }
        init(argv[1]); //initialize parameters for playback + waveform processing
        //open default stereo output stream
        err = Pa_OpenDefaultStream(&stream, 0, 2, paFloat32, samplerate, paFramesPerBufferUnspecified, pa_callback, &data);
        if (err != paNoError){
            return 0;
        }
        //start stream
        err = Pa_StartStream(stream);
        if (err != paNoError){
            return 0;
        }
        /** original infinite loop
        //process waveform in parallel
        while (1){
            analyze();
        }
        **/
        
        ///** loop with GLFW window
        time_increment = 0.8/FRAME_RATE; //guarantees at least 60 FPS updating
        glfwSetTime(0.0);
        while (!glfwWindowShouldClose(window))
        {
            float ratio;
            int width, height;
            
            /** timing statistics
            timeval before, after, elapsed;
            unsigned long usec_remaining;
            gettimeofday(&before, NULL);
            **/
            
            if (glfwGetTime() >= time_increment){
                glfwSetTime(0.0);
                mov_avg_update();
                analyze_async();
            }
            glfwGetFramebufferSize(window, &width, &height);
            ratio = width / (float) height;
            glViewport(0, 0, width, height);
            glClear(GL_COLOR_BUFFER_BIT);
            glMatrixMode(GL_PROJECTION);
            glLoadIdentity();
            glOrtho(-ratio, ratio, -1.f, 1.f, 1.f, -1.f);
            glMatrixMode(GL_MODELVIEW);
            glLoadIdentity();
            
            //TODO: enable lighting, implement basic visualization (w/ shaders)
            
            //Test display of the calculated source locations
            glBegin(GL_TRIANGLES);
            for (int i = 0; i < num_bins; i++){
                float index, offset, rising_edge, intensity;
                index = i/((float)(num_cutoffs));
                offset = (float)(source_sines[i] * ratio);
                if (power_spectra[i] > 0.0){
                    intensity = 1.f + (log10f((float)(power_spectra[i])) / (1.f - log10f((float)(power_spectra[i]))));
                } else {
                    intensity = 0.f;
                }
                if (pspec_mov_avg[i] > 0.0){
                    rising_edge = (float)((power_spectra[i] / (1.25f * pspec_mov_avg[i])) * intensity);
                    rising_edge = (rising_edge > 1.f)? 1.f : ((rising_edge < 0.f)? 0.f : rising_edge);
                } else {
                    rising_edge = 0.f;
                }
                glColor4f(0.f + index, 0.f, 1.f - index, rising_edge);
                glVertex3f(-0.06f + offset, -1.f, 0.f);
                glVertex3f(0.06f + offset, -1.f, 0.f);
                glVertex3f(0.f + offset, 3.f * intensity - 1.f, 0.f);
            }
            glEnd();
            
            /**timing statistics
            gettimeofday(&after, NULL);
            tval_sub(after, before, &elapsed);
            usec_remaining = MILLION_INT/FRAME_RATE - tval_to_usec(elapsed);
            printf("leeway ~= %ld\n", usec_remaining);
            **/
            
            glfwSwapBuffers(window);
            glfwPollEvents();
        }
        //**/
        
        //sleep(1); //Pa_Sleep(1*1000);
        err = Pa_StopStream(stream);
        if (err != paNoError){
            return 0;
        }
        cleanup();
        
        /** close GLFW window **/
        glfwDestroyWindow(window);
        glfwTerminate();
        exit(EXIT_SUCCESS);
        /** **/
    }
    return 0;
}

