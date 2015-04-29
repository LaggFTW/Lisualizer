/**
Functions related to the partitioning of the waveform into different sounds
**/

#include "partitioner.h"

void calc_cutoffs_2ch(fftw_complex *data[2], int *cutoffs, int start, int end, int cut_ind, int num_cutoffs){
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
    //if there were no non-zero values in the range, simply split it evenly in half
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
}

/** (more methods likely to be added) **/

