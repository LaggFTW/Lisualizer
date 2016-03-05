#ifndef PARTITIONER_H
#define PARTITIONER_H

#include <fftw3.h>
#include "sep_proc_data.h"

#ifdef __cplusplus
extern "C" {
#endif

void calc_cutoffs_2ch(fftw_complex *data[2], int *cutoffs, int start, int end, int cut_ind, int num_cutoffs);

#ifdef __cplusplus
}
#endif

#endif
