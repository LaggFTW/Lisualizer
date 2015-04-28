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

#ifndef GLFW_INCLUDE_GLU
    #define GLFW_INCLUDE_GLU
#endif

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "partitioner.h"

#include "visualizer.h"

#include "memutils.h"

//TODO: The code structure of this is still pretty disgusting, fix when implementing a better visualization

typedef struct timeval timeval;

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
    VisualizerData *data = (VisualizerData *)userData;
    float *in = (float *)inputBuffer;
    float *out = (float *)outputBuffer;
    float smp[2];
    int num_samples = framesPerBuffer * 2;
    memset(out, 0, num_samples * sizeof(float));
    for (int i = 0; i < num_samples; i++){
        if (PaUtil_GetRingBufferReadAvailable(&(data->io_buffer))){
            PaUtil_ReadRingBuffer(&(data->io_buffer), &(out[i]), 1);
        }
    }
    for (int i = 0; i < framesPerBuffer; i++){
        if ((PaUtil_GetRingBufferWriteAvailable(&(data->io_buffer)) >= 2) 
                && (PaUtil_GetRingBufferWriteAvailable(&(data->io_buffer)) >= 2)){
            if (sf_read_float(data->inputfile, smp, 2) < 2){
                //for processing purposes, if we encounter EOF, treat it as silence
                smp[0] = 0.0f; smp[1] = 0.0f;
                PaUtil_WriteRingBuffer(&(data->sample_buffer), smp, 2);
                continue;
            }
            PaUtil_WriteRingBuffer(&(data->io_buffer), smp, 2);
            PaUtil_WriteRingBuffer(&(data->sample_buffer), smp, 2);
        }
    }
    return paContinue;
}

//updates the moving average of power spectrum peaks, using a simple moving average approximation
static void mov_avg_update(VisualizerData *data){
    for (int i = 0; i < data->num_bins; i++){
        data->pspec_mov_avg[i] += (data->power_spectra[i] - data->pspec_mov_avg[i])/MOV_AVG_WINDOW_SIZE;
    }
}

static int init(const char *fname, VisualizerData *data){
    SF_INFO fileinfo;

    //information for separation
    data->num_bins = 4;
    data->num_cutoffs = data->num_bins - 1;
    data->freq_cutoffs = (int *)(ucalloc(data->num_cutoffs, sizeof(int)));
    
    //information for visualization
    data->source_sines = (double *)(ucalloc(data->num_bins, sizeof(double)));
    data->power_spectra = (double *)(ucalloc(data->num_bins, sizeof(double)));
    data->pspec_mov_avg = (double *)(ucalloc(data->num_bins, sizeof(double)));
    
    /**
    Intial opening of file
    **/
    data->inputfile = sf_open(fname, SFM_READ, &fileinfo);
    if (!data->inputfile){
        printf("sndfile: file reading error\n");
        return ERRCODE_FAILURE;
    }
    if (fileinfo.channels != 2){
        return ERRCODE_UNSUPPORTED;
    }
    data->samplerate = fileinfo.samplerate;
    for (int i = 0; i < WINDOW_SIZE; i+=2){ //zero-initialize sample arrays and buffers
        data->io_buf_arr[i] = 0.0f;
        data->smp_buf_arr[i] = 0.0f;
        data->samples[i] = 0.0f;
    }
    //initialize ringbuffers
    if (PaUtil_InitializeRingBuffer(&(data->io_buffer), sizeof(float), WINDOW_SIZE*2, data->io_buf_arr)){
        printf("io ringbuffer error\n");
    }
    if (PaUtil_InitializeRingBuffer(&(data->sample_buffer), sizeof(float), WINDOW_SIZE*2, data->smp_buf_arr)){
        printf("sample ringbuffer error\n");
    }

    /**
    Hanning Window Calculation
    **/
    if (HANNING) {
        double pi = 3.141592653589793238463;
        long size = WINDOW_SIZE - 1;
        for (int i = 0; i < WINDOW_SIZE; i++){
            data->hann[i] = 0.5 * (1 - cos(2 * pi * i / size));
        }
    }
    
    /**
    Sliding Window Initialization
    **/
    data->transform_size = WINDOW_SIZE/2 + 1;
    data->cutoff_1600 = (int)((3200.0/data->samplerate) * (WINDOW_SIZE/2)); //cutoff = 1600 / nyquist freq * max transform index
    data->window_increment = data->samplerate / FRAME_RATE * 2; //for static window shift, if needed
    for (int i = 0; i < 2; i++){
        data->window[i] = (double *)(fftw_malloc(sizeof(double)*WINDOW_SIZE));
        data->dft_window[i] = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex)*data->transform_size));
        data->plan[i] = fftw_plan_dft_r2c_1d(WINDOW_SIZE, data->window[i], data->dft_window[i], FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
        
        data->dft_bin[i] = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex)*data->transform_size));
    }
    data->dft_corre = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex)*data->transform_size));
    data->plancorre = fftw_plan_dft_c2r_1d(WINDOW_SIZE, data->dft_corre, data->window[0], FFTW_ESTIMATE);
    return ERRCODE_SUCCESS;
}

/**
Sliding Window Analysis, no synchronization
**/
static void analyze_async(VisualizerData *data){
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
    dyn_win_shift = PaUtil_GetRingBufferReadAvailable(&(data->sample_buffer));
    dyn_win_shift -= dyn_win_shift % 2; //ensure the number of samples shifted is a multiple of the number of channels (2)
    orig_win_offset = WINDOW_SIZE * 2 - dyn_win_shift;
    //bounds checking (pray that the below if block never executes)
    if (orig_win_offset < 0){
        orig_win_offset = 0;
        dyn_win_shift = WINDOW_SIZE * 2;
    }
    //shift over the old samples by dyn_win_shift
    for (int j = 0; j < orig_win_offset; j+=2){
        data->samples[j] = data->samples[j+dyn_win_shift];
        data->samples[j+1] = data->samples[j+1+dyn_win_shift];
    }
    //copy dyn_win_shift samples over to the space left over from the shift
    PaUtil_ReadRingBuffer(&(data->sample_buffer), &(data->samples[orig_win_offset]), dyn_win_shift);
    i = 0;
    //split the channels, applying a windowing function if needed
    for (int k = 0; k < WINDOW_SIZE; i+=2, k++){
        if (HANNING){
            data->window[0][k] = data->samples[i] * data->hann[k];
            data->window[1][k] = data->samples[i+1] * data->hann[k];
        } else {
            data->window[0][k] = data->samples[i];
            data->window[1][k] = data->samples[i+1];
        }
    }
    //fourier analysis (separation, localization)
    fftw_execute(data->plan[0]);
    fftw_execute(data->plan[1]);
    //separation of the DFTs
    calc_cutoffs_2ch(data->dft_window, data->freq_cutoffs, 0, data->transform_size, data->num_cutoffs/2, data->num_bins/2);
    for (i = 0; i < data->num_bins; i++){
        double diff, intensity; //values to be used for ILD
        int loc_method = 0; //determines whether to use ITD or ILD
        data->power_spectra[i] = 0.0; //reset the current power spectrum values
        for (int k = 0; k < 2; k++){
            for (int j = 0; j < data->transform_size; j++){
                if ((i == 0 || j > data->freq_cutoffs[i-1]) && (i == data->num_cutoffs || j <= data->freq_cutoffs[i])){
                    memcpy(data->dft_bin[k][j], data->dft_window[k][j], sizeof(fftw_complex));
                } else {
                    data->dft_bin[k][j][0] = 0;
                    data->dft_bin[k][j][1] = 0;
                }
            }
        }
        if (i > 0){
            if (data->freq_cutoffs[i-1] >= data->cutoff_1600){
                loc_method = 1;
                diff = 0.0;
                intensity = 0.0;
            } else {
                if (i < data->num_cutoffs && (data->freq_cutoffs[i-1] + data->freq_cutoffs[i]) >= (data->cutoff_1600 * 2)){
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
        for (int j = 0; j < data->transform_size; j++){
            double norm_l, norm_r, norm_a;
            fftw_complex **dft_bin = data->dft_bin;
            norm_l = (dft_bin[0][j][0] * dft_bin[0][j][0] + dft_bin[0][j][1] * dft_bin[0][j][1])/(HANN_SCALING_FACTOR * window_sq);
            norm_r = (dft_bin[1][j][0] * dft_bin[1][j][0] + dft_bin[1][j][1] * dft_bin[1][j][1])/(HANN_SCALING_FACTOR * window_sq);
            norm_a = (norm_l + norm_r) / 2;
            //compute peak power
            if (norm_a > data->power_spectra[i]){
                data->power_spectra[i] = norm_a;
            }
            //prepare for sound localization
            if (loc_method){
                //audio above ~1600 Hz: use ILD
                diff += norm_r - norm_l;
                intensity += norm_a;
            } else {
                //audio below ~1600 Hz: use ITD, prepare for cross-correlation
                data->dft_corre[j][0] = (dft_bin[0][j][0] * dft_bin[1][j][0] + dft_bin[0][j][1] * dft_bin[1][j][1])/window_sq;
                data->dft_corre[j][1] = (dft_bin[0][j][0] * dft_bin[1][j][1] - dft_bin[0][j][1] * dft_bin[1][j][0])/window_sq;
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
            data->source_sines[i] = ratio;
        } else {
            //ITD w/ cross-correlation
            long max_index = 0;
            long max_delay, delay_value;
            double factor, ratio;
            
            fftw_execute(data->plancorre); //finish computing cross-correlation
            
            factor = SOUND_SPEED / (data->samplerate * HEAD_WIDTH);
            //maximum possible delay such that the angle calculated will be 90 degrees; 
            //used to distinguish between delay caused by sound coming in from an angle 
            //and delay caused by echoes
            max_delay = (long)((1 / factor)) + 2;
            if (max_delay > WINDOW_SIZE){
                max_delay = WINDOW_SIZE;
            }
            for (long j = 1; j < max_delay; j++){
                if (data->window[0][j] > data->window[0][max_index]){
                    max_index = j;
                }
                if (data->window[0][WINDOW_SIZE - j] > data->window[0][max_index]){
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
            data->source_sines[i] = ratio;
        }
    }
}

static int cleanup(VisualizerData *data){
    for (int i = 0; i < 2; i++){
        fftw_destroy_plan(data->plan[i]);
        fftw_free(data->window[i]);
        fftw_free(data->dft_window[i]);
        
        fftw_free(data->dft_bin[i]);
    }
    free(data->freq_cutoffs);
    free(data->source_sines);
    free(data->power_spectra);
    free(data->pspec_mov_avg);
    fftw_destroy_plan(data->plancorre);
    fftw_free(data->dft_corre);
    fftw_cleanup();
    if (sf_close(data->inputfile)){
        printf("sndfile: error closing file\n");
    }
}

int visualizer_start(int argc, char **argv){
    if (argc >= 2){
        PaStream *stream;
        PaError err;
        VisualizerData data;
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
        init(argv[1], &data); //initialize parameters for playback + waveform processing
        //open default stereo output stream
        err = Pa_OpenDefaultStream(&stream, 0, 2, paFloat32, data.samplerate, paFramesPerBufferUnspecified, pa_callback, &data);
        if (err != paNoError){
            return 0;
        }
        //start stream
        err = Pa_StartStream(stream);
        if (err != paNoError){
            return 0;
        }
        
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
                mov_avg_update(&data);
                analyze_async(&data);
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
            for (int i = 0; i < data.num_bins; i++){
                float index, offset, rising_edge, intensity;
                index = i/((float)(data.num_cutoffs));
                offset = (float)(data.source_sines[i] * ratio);
                if (data.power_spectra[i] > 0.0){
                    intensity = 1.f + (log10f((float)(data.power_spectra[i])) / (1.f - log10f((float)(data.power_spectra[i]))));
                } else {
                    intensity = 0.f;
                }
                if (data.pspec_mov_avg[i] > 0.0){
                    rising_edge = (float)((data.power_spectra[i] / (1.25f * data.pspec_mov_avg[i])) * intensity);
                    rising_edge = (rising_edge > 1.f)? 1.f : rising_edge;
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
        cleanup(&data);
        
        /** close GLFW window **/
        glfwDestroyWindow(window);
        glfwTerminate();
        exit(EXIT_SUCCESS);
        /** **/
    }
    return 0;
}

