#include <cstdio>
#include <cstring>
#include <iostream>
#include <fstream>

#include <math.h>

#include <sndfile.hh>

//TODO: Figure out why gcc/g++ -libmp3lame flag gives errors
//#include <lame/lame.h>

#include <fftw3.h>

#define ERRCODE_SUCCESS 0

//flag: 1 for a linear cross correlation, 0 for a circular cross correlation
#define LINEAR 0

//flag: 1 for a hanning window, 0 for a rectangular window
#define HANNING 0

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

static int write_file(const char *fname, int format, short *samples_arr, long num_samples, int num_channels, int sample_rate){
    SndfileHandle file = SndfileHandle(fname, SFM_WRITE, format, num_channels, sample_rate);
    file.write(samples_arr, num_samples);
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
    cut = (int)(numer/denom);
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

static int xcorrelation(short *samples_left, short *samples_right, long num_samples, long dft_size, long *output){
    double *samples_l, *samples_r, *samples_corre;
    fftw_complex *samples_dft_l, *samples_dft_r, *samples_dft_corre;
    fftw_plan forward_l, forward_r;
    fftw_plan backward;
    long transform_size;
    
    transform_size = dft_size/2 + 1;
    samples_l = (double *)(fftw_malloc(sizeof(double) * dft_size));
    samples_r = (double *)(fftw_malloc(sizeof(double) * dft_size));
    samples_corre = (double *)(fftw_malloc(sizeof(double) * dft_size));
    samples_dft_l = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex) * transform_size));
    samples_dft_r = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex) * transform_size));
    samples_dft_corre = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex) * transform_size));
    for (long i = 0; i < dft_size; i++){
        if (i < num_samples){
            double val = 1;
            if (HANNING){
                val = 0.5 * (1 - cos(2 * 3.141592653589793238463 * i / (num_samples - 1)));
            }
            samples_l[i] = samples_left[i] * val;
            samples_r[i] = samples_right[i] * val;
        } else {
            samples_l[i] = 0;
            samples_r[i] = 0;
        }
    }
    forward_l = fftw_plan_dft_r2c_1d(dft_size, samples_l, samples_dft_l, FFTW_ESTIMATE);
    forward_r = fftw_plan_dft_r2c_1d(dft_size, samples_r, samples_dft_r, FFTW_ESTIMATE);
    
    //run dft on the input channel samples
    fftw_execute(forward_l);
    fftw_execute(forward_r);
    
    //multiply the two fourier transforms
    for (long i = 0; i < transform_size; i++){
        samples_dft_corre[i][0] = (samples_dft_l[i][0] * samples_dft_r[i][0] + samples_dft_l[i][1] * samples_dft_r[i][1]) / dft_size;
        samples_dft_corre[i][1] = (samples_dft_l[i][0] * samples_dft_r[i][1] - samples_dft_l[i][1] * samples_dft_r[i][0]) / dft_size;
    }
    
    backward = fftw_plan_dft_c2r_1d(dft_size, samples_dft_corre, samples_corre, FFTW_ESTIMATE);
    
    //run idft for the cross-correlation
    fftw_execute(backward);
    
    //scale and copy the data into the output
    for (long i = 0; i < dft_size; i++){
        //printf("%lf\t", samples_corre[i]/dft_size);
        if (LINEAR){
            output[i] = (long)(samples_corre[((i+num_samples)%dft_size)]/dft_size);
        } else {
            output[i] = (long)(samples_corre[i]/dft_size);
        }
    }
    
    //free dynamic memory
    fftw_destroy_plan(forward_l);
    fftw_destroy_plan(forward_r);
    fftw_destroy_plan(backward);
    fftw_free(samples_l);
    fftw_free(samples_r);
    fftw_free(samples_corre);
    fftw_free(samples_dft_l);
    fftw_free(samples_dft_r);
    fftw_free(samples_dft_corre);
    return ERRCODE_SUCCESS;
}

int main(int argc, char** argv){
    if (argc >= 2){
        long num_samples;
        long samples_per_ch;
        short *samples;
        long *correlation;
        short **channels;
        int channel_ct;
        int sample_rate;
        
        //read in the file
        read_file(argv[1], &samples, &sample_rate, &num_samples, &channel_ct);
        
        //separate audio file into its channels
        samples_per_ch = num_samples/((long)(channel_ct));
        channels = new short *[channel_ct];
        correlation = new long[num_samples];
        for (int i = 0; i < channel_ct; i++){
            channels[i] = new short[samples_per_ch];
        }
        for (long i = 0; i < num_samples; i++){
            channels[(i%channel_ct)][(i/channel_ct)] = samples[i];
        }
        
        //cross-correlate
        if (LINEAR){
            xcorrelation(channels[0], channels[1], samples_per_ch, num_samples - 1, correlation);
        } else {
            xcorrelation(channels[0], channels[1], samples_per_ch, samples_per_ch, correlation);
        }
        
        /**
        //write to csv
        {
            std::ofstream fileout;
            fileout.open("./tmp/output.csv");
            for (long i = 0; i < num_samples; i++){
                char *s;
                int size;
                const char *s_format = "%ld,\n";
                size = snprintf(NULL, 0, s_format, correlation[i]) + 1;
                s = new char[size];
                snprintf(s, size, s_format, correlation[i]);
                fileout << s;
                delete[] s;
            }
            fileout.close();
        }
        **/
        
        {
            long max_index = 0;
            long value;
            for (long i = 1; i < num_samples; i++){
                if (correlation[i] > correlation[max_index]){
                    max_index = i;
                }
            }
            if (LINEAR){
                value = (samples_per_ch - max_index - 1);
            } else {
                value = ((max_index > (samples_per_ch / 2))? samples_per_ch - max_index : -max_index);
            }
            printf("sample offset: %ld\n", value);
            printf("if the value is negative, then the left precedes the right by the absolute value\n");
            printf("otherwise, the right precedes the left by the value shown\n");
        }
        
        //free dynamic memory
        for (int i = 0; i < channel_ct; i++){
            delete[] channels[i];
        }
        delete[] channels;
        delete[] correlation;
        delete[] samples;
    }
    return 0;
}

