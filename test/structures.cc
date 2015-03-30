#include <cstdio>
#include <cstring>

#include <math.h>

#include <sndfile.hh>

#include <portaudio.h>

//TODO: Figure out why gcc/g++ -libmp3lame flag gives errors
//#include <lame/lame.h>

#include <fftw3.h>

#define ERRCODE_SUCCESS 0

/**
 TODO:
 Actually fill in these classes, possibly adding other classes along the way
**/

class AudioTrackHandler{
    //class to store information regarding an audio track such as filename,
    //sample rate, codec type, array of frames/samples; also handles 
    //reading/writing of the audio file
}

class AudioStreamHandler{
    //class to store information regarding an audio stream such as stream id,
    //channels, sample rate, frame/sample buffer; also handles reading from the
    //stream via portaudio
}

class ProcessedWaveform{
    //class to store information regarding a waveform of a clip of audio such as
    //arrays of frames/samples organized by channel, arrays for each channel's
    //DFT, sample rate; also possibly handles processing of itself
}

/**
 TODO:
 Classes/structures for GUI system w/ wxWidgets
 Classes/structures to handle Visualization calculations + rendering w/ OpenGL
**/

