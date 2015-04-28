#ifndef MEMUTILS_H
#define MEMUTILS_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void *umalloc(size_t size);
void *ucalloc(size_t nmemb, size_t size);
void *urealloc(void *ptr, size_t size);

void *cbmalloc(size_t size, void (*callback)(void *), void *data);
void *cbcalloc(size_t nmemb, size_t size, void (*callback)(void *), void *data);
void *cbrealloc(void *ptr, size_t size, void (*callback)(void *), void *data);

#ifdef __cplusplus
}
#endif

#endif
