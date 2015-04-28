#include <stdio.h>
#include "memutils.h"

/**
basic error-handling wrappers for dynamic allocation functions in stdlib.h
umalloc, ucalloc, urealloc: print error message and exit on failure
cbmalloc, cbcalloc, cbrealloc: call the callback function on failure
**/

typedef void (*cb_func_ptr)(void *);

void *umalloc(size_t size){
    void *to_ret = malloc(size);
    if (!to_ret){
        fprintf(stderr, "error: not enough space on heap; exiting\n");
        exit(EXIT_FAILURE);
    }
    return to_ret;
}

void *ucalloc(size_t nmemb, size_t size){
    void *to_ret = calloc(nmemb, size);
    if (!to_ret){
        fprintf(stderr, "error: not enough space on heap; exiting\n");
        exit(EXIT_FAILURE);
    }
    return to_ret;
}

void *urealloc(void *ptr, size_t size){
    void *to_ret = realloc(ptr, size);
    if (!to_ret){
        fprintf(stderr, "error: heap reallocation error; exiting\n");
        exit(EXIT_FAILURE);
    }
    return to_ret;
}

void *cbmalloc(size_t size, cb_func_ptr func, void *data){
    void *to_ret = malloc(size);
    if (!to_ret){
        func(data);
    }
    return to_ret;
}

void *cbcalloc(size_t nmemb, size_t size, cb_func_ptr func, void *data){
    void *to_ret = calloc(nmemb, size);
    if (!to_ret){
        func(data);
    }
    return to_ret;
}

void *cbrealloc(void *ptr, size_t size, cb_func_ptr func, void *data){
    void *to_ret = realloc(ptr, size);
    if (!to_ret){
        func(data);
    }
    return to_ret;
}

