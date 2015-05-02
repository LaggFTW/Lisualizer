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
    double head_width; //distance between left and right audio receivers (in this case, ears), scaled by a damping factor
    double max_damping; //maximum interaural level difference, as a normalized ratio
    int itd_cutoff; //index of highest frequency to attempt ITD on
    int ild_cutoff; //index of lowest frequency to attempt ILD on
} localizer_data;

#ifdef __cplusplus
}
#endif

#endif
