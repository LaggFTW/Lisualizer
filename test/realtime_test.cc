#include <cstdio>
#include <cstring>

#include <math.h>

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
#define WINDOW_SIZE 2048

#define FRAME_RATE 60

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
        long num_samples, window_start, window_end, window_increment, window_increment2;
        double *window[2];
        fftw_complex *dft_window[2];
        fftw_plan plan[2];
        short *samples;
        int transform_size;
        {
            int channel_ct;
            int sample_rate;
            //read in the file
            read_file(argv[1], &samples, &sample_rate, &num_samples, &channel_ct);
            if (channel_ct != 2 || sample_rate != SAMPLE_RATE){
                printf("Unsupported sample rate or number of channels; exiting\n");
                return 0;
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
        }
        window_start = -WINDOW_SIZE+2;
        window_end = WINDOW_SIZE+2;
        window_increment = SAMPLE_RATE/FRAME_RATE;
        window_increment2 = window_increment*2;
        //outer loop; shifts window position, will likely end up shifting it multiple frames a time (i.e. i+=(SAMPLE_RATE/FRAME_RATE)), 
        //as there is very little need to take the DFT every single time we get a new sample
        for (long i = 0; i < num_samples; i+=window_increment){
            //inner loop; updates window values
            int k = 0;
            for (int j = window_start; j < window_end; j+=2, k++){
                if (j < 0 || j >= num_samples){
                    window[0][k] = 0;
                    window[1][k] = 0;
                } else {
                    window[0][k] = samples[j];
                    window[1][k] = samples[j+1];
                }
            }
            //fourier analysis (separation, localization)
            fftw_execute(plan[0]);
            fftw_execute(plan[1]);
            //TODO: separation of the DFTs
            //TODO: localization (take cross-correlation of the left and right using the DFTs already computed, compute a source angle)
            //the above will require creating a few more arrays of doubles/fftw_complexes/fftw_plans, most likely
            window_start += window_increment2; window_end += window_increment2;
        }
        for (int i = 0; i < 2; i++){
            fftw_destroy_plan(plan[i]);
            fftw_free(window[i]);
            fftw_free(dft_window[i]);
        }
        fftw_cleanup();
        delete[] samples;
    }
    return 0;
}
