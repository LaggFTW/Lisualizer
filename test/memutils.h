#ifndef MEM_UTILS_POI
#define MEM_UTILS_POI

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

void *umalloc(size_t);
void *ucalloc(size_t, size_t);
void *urealloc(void *, size_t);

void *cbmalloc(size_t size, void (*)(void *), void *data);
void *cbcalloc(size_t nmemb, size_t size, void (*)(void *), void *data);
void *cbrealloc(void *ptr, size_t size, void (*)(void *), void *data);

#ifdef __cplusplus
}
#endif

#endif
