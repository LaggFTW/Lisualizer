#include <cstdio>
#include <cstring>

#include <math.h>

#include <sndfile.hh>

//TODO: Figure out why gcc/g++ -libmp3lame flag gives errors
//#include <lame/lame.h>

#include <fftw3.h>

#define ERRCODE_SUCCESS 0

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

/** 
TODO: figure out frequency cutoffs (linearly is not gonna work, logarithmic would probably work, 
but something based on frequency statistics of the source waveform might work the best)
**/
static int freq_bin(short *samples_arr, long num_samples, int num_bins, int *freq_cutoffs, short ***freq_filtered){
    double *samples;
    double **samples_idfts;
    fftw_complex *samples_dft;
    fftw_complex **bin_dfts;
    fftw_plan forward;
    fftw_plan *backward;
    long transform_size, bin_width;
    
    transform_size = num_samples/2 + 1;
    samples = (double *)(fftw_malloc(sizeof(double) * num_samples));
    samples_dft = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex) * transform_size));
    for (long i = 0; i < num_samples; i++){
        samples[i] = samples_arr[i];
    }
    forward = fftw_plan_dft_r2c_1d(num_samples, samples, samples_dft, FFTW_ESTIMATE);
    
    //run dft on the input channel samples
    fftw_execute(forward);
    
    //set up bins
    if (num_bins > (transform_size / 2)){
        num_bins = transform_size;
    }
    bin_width = transform_size / num_bins;
    *freq_filtered = new short *[num_bins];
    samples_idfts = new double *[num_bins];
    bin_dfts = new fftw_complex *[num_bins];
    backward = new fftw_plan[num_bins];
    
    for (int i = 0; i < num_bins; i++){
        (*freq_filtered)[i] = new short[num_samples];
        samples_idfts[i] = (double *)(fftw_malloc(sizeof(double) * num_samples));
        bin_dfts[i] = (fftw_complex *)(fftw_malloc(sizeof(fftw_complex) * transform_size));
        for (long j = 0; j < transform_size; j++){
            (bin_dfts[i][j])[0] = 0;
            (bin_dfts[i][j])[1] = 0;
        }
        backward[i] = fftw_plan_dft_c2r_1d(num_samples, bin_dfts[i], samples_idfts[i], FFTW_ESTIMATE);
    }
    
    /**
    //split the frequency data up into these bins linearly
    for (long i = 0; i < transform_size; i++){
        int bin_index = i/bin_width;
        if (bin_index >= num_bins){
            bin_index = num_bins - 1;
        }
        memcpy(bin_dfts[bin_index][i], samples_dft[i], sizeof(fftw_complex));
    }
    **/
    
    /**
    //split the frequency data up into bins logarithmically
    //(well, it's really just a quad root to approximate log behavior)
    for (long i = 0; i < transform_size; i++){
        int bin_index = (int)(sqrt(sqrt(((double)(i))/transform_size))*num_bins);
        if (bin_index >= num_bins){
            bin_index = num_bins - 1;
        }
        memcpy(bin_dfts[bin_index][i], samples_dft[i], sizeof(fftw_complex));
    }
    **/
    
    ///**
    //split data into quartiles based on weighted averages of the values in the dft
    //(i.e. if the waveform is bass-heavy, then the splits would be mostly at lower frequencies,
    //while if the waveform is treble-heavy, the splits/cutoffs would be higher)
    //alternatively, approximate the split locations by taking the locations with the
    //largest values in the dft
    //ASSUMPTION: num_bins is a power of 2
    {
        int num_cutoffs = num_bins - 1;
        freq_cutoffs = new int[num_cutoffs];
        calc_cutoffs(samples_dft, freq_cutoffs, 0, transform_size, num_cutoffs/2, num_bins/2);
        {
            int bin_index = 0;
            for (long i = 0; i < transform_size; i++){
                while (bin_index < num_cutoffs && i > freq_cutoffs[bin_index]){
                    bin_index += 1;
                }
                memcpy(bin_dfts[bin_index][i], samples_dft[i], sizeof(fftw_complex));
            }
        }
        delete[] freq_cutoffs;
    }
    //**/
    
    //run idft on the bin_dfts, placing them into their respective samples_idfts bin
    for (int i = 0; i < num_bins; i++){
        fftw_execute(backward[i]);
    }
    
    //scale and copy the data from samples_idfts into the output: freq_filtered
    for (int i = 0; i < num_bins; i++){
        for (long j = 0; j < num_samples; j++){
            (*freq_filtered)[i][j] = (short)(samples_idfts[i][j]/num_samples);
        }
    }
    
    //free dynamic memory
    for (int i = 0; i < num_bins; i++){
        fftw_destroy_plan(backward[i]);
        fftw_free(samples_idfts[i]);
        fftw_free(bin_dfts[i]);
    }
    fftw_destroy_plan(forward);
    fftw_free(samples);
    fftw_free(samples_dft);
    delete[] backward;
    delete[] samples_idfts;
    delete[] bin_dfts;
    return ERRCODE_SUCCESS;
}

int main(int argc, char** argv){
    if (argc >= 2){
        long num_samples;
        long samples_per_ch;
        short ***freq_bin_samples;
        short *samples;
        short **channels;
        int channel_ct;
        int sample_rate;
        int num_bins = 4;
        
        //read in the file
        read_file(argv[1], &samples, &sample_rate, &num_samples, &channel_ct);
        
        //separate audio file into its channels
        samples_per_ch = num_samples/((long)(channel_ct));
        channels = new short *[channel_ct];
        freq_bin_samples = new short **[channel_ct];
        for (int i = 0; i < channel_ct; i++){
            channels[i] = new short[samples_per_ch];
        }
        for (long i = 0; i < num_samples; i++){
            channels[(i%channel_ct)][(i/channel_ct)] = samples[i];
        }
        
        //do basic frequency binning and analysis
        for (int i = 0; i < channel_ct; i++){
            freq_bin(channels[i], samples_per_ch, num_bins, NULL, &(freq_bin_samples[i]));
        }
        
        //combine channels, write to output
        for (int i = 0; i < num_bins; i++){
            char *fname;
            int size;
            const char *fname_format = "./audioout/bin_%d.wav";
            size = snprintf(NULL, 0, fname_format, i) + 1;
            fname = new char[size];
            snprintf(fname, size, fname_format, i);
            for (long j = 0; j < num_samples; j++){
                samples[j] = freq_bin_samples[(j%channel_ct)][i][(j/channel_ct)];
            }
            printf("Writing to file named '%s'\n", fname);
            write_file(fname, SF_FORMAT_WAV | SF_FORMAT_PCM_16, samples, num_samples, channel_ct, sample_rate);
        }
        
        //free dynamic memory
        for (int i = 0; i < channel_ct; i++){
            delete[] channels[i];
            for (int j = 0; j < num_bins; j++){
                delete[] freq_bin_samples[i][j];
            }
            delete[] freq_bin_samples[i];
        }
        delete[] channels;
        delete[] freq_bin_samples;
        delete[] samples;
    }
    return 0;
}

