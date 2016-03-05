#ifndef SEP_PROC_DATA_H
#define SEP_PROC_DATA_H

#include <fftw3.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    /**
    Parameters for waveform sound separation/partitioning processing
    Data includes necessary information for partitioning, which can then
    also be used as an input to compute other data for visualization purposes
    **/
    
    //TODO: determine the necessary parameters for this, and update this struct
    //temporary data/params below; will be heavily revised,
    //and will likely have many more parameters
    double *window[2];
    fftw_complex *dft_window[2];
    fftw_plan plan[2];
    int window_size;
    int transform_size;
    
    //separation cutoff data
    fftw_complex *dft_bin[2];
    int *freq_cutoffs;
    int num_bins, num_cutoffs;
} separator_data;

#ifdef __cplusplus
}
#endif

#endif
