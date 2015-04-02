#include <cstdio>
#include <cstring>

#include <math.h>
#include <unistd.h>

#include <sndfile.hh>

/**
#include <portaudio.h>

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

#define ERRCODE_SUCCESS 0

#define SAMPLE_RATE 44100
#define WINDOW_SIZE 32768

#define FRAME_RATE 60

//speed of sound, in meters per second
#define SOUND_SPEED 343.59

//distance between left and right audio receivers (in this case, ears), scaled by a damping factor
#define HEAD_WIDTH 0.35

//flag: 1 for a hanning window, 0 for a rectangular window
#define HANNING 1

/**
static int pa_callback(const void *inputBuffer, void *outputBuffer, 
        unsigned long framesPerBuffer, const PaStreamCallbackTimeInfo* timeInfo,
        PaStreamCallbackFlags statusFlags, void *userData){
    return 0;
}
**/

static int read_file(const char *fname, short **samples_arr, int* sample_rate, long *num_samples, int *num_channels){
    long frame_ct; //number of samples
    SndfileHandle file; //object from libsndfile
    
	file = SndfileHandle(fname);
	*sample_rate = file.samplerate();
	*num_channels = file.channels();
	frame_ct = file.frames();
	*num_samples = frame_ct * (*num_channels); //total number of samples (each individual sample in every audio channel is included separately)
	
	printf("Opened file '%s'\n", fname);
	printf("    Sample rate : %d\n", *sample_rate);
	printf("    Channels    : %d\n", *num_channels);
	printf("    Frames      : %ld\n", frame_ct);
	
	*samples_arr = new short[*num_samples];
	file.read(*samples_arr, *num_samples);
	return ERRCODE_SUCCESS;
}

static int calc_cutoffs(fftw_complex *data, int *cutoffs, int start, int end, int cut_ind, int num_cutoffs){
    //notes on the input: cutoffs is an array of (num_cutoffs * 2) - 1 elements, which will store the
    //cutoffs indices based on the weighted average of values in data, and thus will allow the
    //partitioning of the data into (num_cutoffs * 2) groups
    //this method also assumes that num_cutoffs is a power of 2
    double numer, denom;
    int cut = 0;
    numer = 0; denom = 0;
    for (int i = start; i < end; i++){
        double norm = data[i][0] * data[i][0] + data[i][1] * data[i][1];
        numer += i * norm;
        denom += norm;
    }
    cut = (denom > 0.0 || denom < 0.0)? (int)(numer/denom) : ((start + end)/2);
    cutoffs[cut_ind] = cut;
    num_cutoffs /= 2;
    if (num_cutoffs >= 1){
        int cut_before, cut_after;
        cut_before = cut_ind - num_cutoffs;
        cut_after = cut_ind + num_cutoffs;
        calc_cutoffs(data, cutoffs, start, cut, cut_before, num_cutoffs);
        calc_cutoffs(data, cutoffs, cut, end, cut_after, num_cutoffs);
    }
    return ERRCODE_SUCCESS;
}

int main(int argc, char **argv){
    /**
    PaStream *stream;
    PaError err;
    paTestData data;
    err = Pa_Initialize();
    if(err != paNoError){
        return 0;
    }
    err = Pa_OpenDefaultStream(&stream, 2, 0, paFloat32, SAMPLE_RATE, &data);
    if (err != paNoError){
        return 0;
    }
    err = Pa_StartStream(stream);
    if (err != paNoError){
        return 0;
    }
    sleep(10); //Pa_Sleep(10*1000);
    err = Pa_StopStream(stream);
    if (err != paNoError){
        return 0;
    }
    **/
    if (argc >= 2){
        long num_samples, window_start, window_end, window_increment;
        double *window[2];
        double hann[WINDOW_SIZE];
        fftw_complex *dft_window[2];
        fftw_complex *dft_corre;
        fftw_plan plan[2];
        fftw_plan plancorre;
        short *samples;
        int transform_size;
        
        //separation cutoff data
        fftw_complex *dft_bin[2];
        int *freq_cutoffs[2];
        int num_bins, num_cutoffs;
        num_bins = 4;
        num_cutoffs = num_bins - 1;
        freq_cutoffs[0] = new int[num_cutoffs];
        freq_cutoffs[1] = new int[num_cutoffs];
        
        //read in the file
        {
            int channel_ct;
            int sample_rate;
            read_file(argv[1], &samples, &sample_rate, &num_samples, &channel_ct);
            if (channel_ct != 2 || sample_rate != SAMPLE_RATE){
                printf("Unsupported sample rate or number of channels; exiting\n");
                return 0;
            }
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
        Sliding Window Analysis
        **/
        transform_size = WINDOW_SIZE/2 + 1;
        for (int i = 0; i < 2; i++){
            window[i] = (double *)(fftw_malloc(sizeof(double)*WINDOW_SIZE));
            dft_window[i] = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex)*transform_size));
            plan[i] = fftw_plan_dft_r2c_1d(WINDOW_SIZE, window[i], dft_window[i], FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
            
            dft_bin[i] = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex)*transform_size));
        }
        dft_corre = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex)*transform_size));
        plancorre = fftw_plan_dft_c2r_1d(WINDOW_SIZE, dft_corre, window[0], FFTW_ESTIMATE);
        
        window_start = -WINDOW_SIZE+2;
        window_end = WINDOW_SIZE+2;
        window_increment = SAMPLE_RATE/FRAME_RATE * 2;
        //outer loop; shifts window position, will likely end up shifting it multiple frames a time (i.e. i+=(SAMPLE_RATE/FRAME_RATE)), 
        //as there is very little need to take the DFT every single time we get a new sample
        for (long w = 0; w < num_samples; w+=window_increment){
            //inner loop; updates window values
            int k = 0;
            for (long j = window_start; j < window_end; j+=2, k++){
                if (j < 0 || j >= num_samples){
                    window[0][k] = 0;
                    window[1][k] = 0;
                } else {
                    if (HANNING){
                        window[0][k] = samples[j] * hann[k];
                        window[1][k] = samples[j+1] * hann[k];
                    } else {
                        window[0][k] = samples[j];
                        window[1][k] = samples[j+1];
                    }
                }
            }
            //fourier analysis (separation, localization)
            fftw_execute(plan[0]);
            fftw_execute(plan[1]);
            //separation of the DFTs
            calc_cutoffs(dft_window[0], freq_cutoffs[0], 0, transform_size, num_cutoffs/2, num_bins/2);
            calc_cutoffs(dft_window[1], freq_cutoffs[1], 0, transform_size, num_cutoffs/2, num_bins/2);
            for (int i = 0; i < num_bins; i++){
                for (int k = 0; k < 2; k++){
                    for (int j = 0; j < transform_size; j++){
                        if ((i == 0 || j > freq_cutoffs[k][i-1]) && (i == num_cutoffs || j <= freq_cutoffs[k][i])){
                            memcpy(dft_bin[k][j], dft_window[k][j], sizeof(fftw_complex));
                        } else {
                            dft_bin[k][j][0] = 0;
                            dft_bin[k][j][1] = 0;
                        }
                    }
                }
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
                    factor = SOUND_SPEED / (SAMPLE_RATE * HEAD_WIDTH);
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
            window_start += window_increment; window_end += window_increment;
        }
        for (int i = 0; i < 2; i++){
            fftw_destroy_plan(plan[i]);
            fftw_free(window[i]);
            fftw_free(dft_window[i]);
            
            delete[] freq_cutoffs[i];
            fftw_free(dft_bin[i]);
        }
        fftw_destroy_plan(plancorre);
        fftw_free(dft_corre);
        fftw_cleanup();
        delete[] samples;
    }
    return 0;
}
