CC = gcc
INCLUDE_DIR = ../include
CFLAGS = -std=gnu99 -iquote $(INCLUDE_DIR)

L_DIR = ../lib

# change as necessary: -lsndfile -lfftw3 -lportaudio `pkg-config --cflags glfw3 glu` `pkg-config --libs glfw3 glu` -lGLEW -lm
LIBS = -lfftw3 -lm

O_DIR = ./obj/
ELF_DIR = ../bin/

# uncomment this once targets have been set up, and set up all and clean targets
# .PHONY: all clean

# gcc -std=gnu99 -o test.out realtime_playback_test2.c pa_ringbuffer.c memutils.c -lsndfile -lfftw3 -lportaudio `pkg-config --cflags glfw3 glu` `pkg-config --libs glfw3 glu` -lGLEW -lm

# gcc -std=gnu99 -o test.out main.c visualizer.c partitioner.c pa_ringbuffer.c memutils.c -lsndfile -lfftw3 -lportaudio `pkg-config --cflags glfw3 glu` `pkg-config --libs glfw3 glu` -lGLEW -lm

