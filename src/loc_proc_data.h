#ifndef LOC_PROC_DATA_H
#define LOC_PROC_DATA_H

#include <fftw3.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    /**
    Parameters for waveform localization processing
    Data includes necessary information for localization processing, which can
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
    int cutoff; //index of dft coefficient corresponding to the cutoff determining whether to use ITD or ILD
} localizer_data;

#ifdef __cplusplus
}
#endif

#endif
