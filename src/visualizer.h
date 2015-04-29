#ifndef VISUALIZER_DATA_H
#define VISUALIZER_DATA_H

#include <portaudio.h>
#include <sndfile.h>
#include <fftw3.h>
#include "pa_ringbuffer.h"

#ifdef __cplusplus
extern "C" {
#endif

#define ERRCODE_SUCCESS 0
#define ERRCODE_EOF 1
#define ERRCODE_FAILURE 2
#define ERRCODE_UNSUPPORTED -1

#define MILLION_INT 1000000

#define WINDOW_SIZE 4096

#define MOV_AVG_WINDOW_SIZE 4096

#define FRAME_RATE 60

#define WALL_RES 64

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

typedef struct {
    /**
    Parameters for processing
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
} VisualizerData;

int visualizer_start(int argc, char **argv);

#ifdef __cplusplus
}
#endif

#endif
