#include "malloc.h"
#include <stdlib.h>
#include <stdio.h>

/* malloc監視 */
void *MyMalloc_Rep(size_t sz, const char *pcFileName, int nLine) {
    void *ptr;
    ptr = malloc(sz);
    fprintf(stderr, "malloc at %s %-4d : %p - %p, %ld byten", pcFileName, nLine, ptr, ptr + sz, sz);
    return ptr;
}

void MyFree_Rep(void *ptr, const char *pcFileName, int nLine) {
    fprintf(stderr, "freeing at %s %-4d : %pn", pcFileName, nLine, ptr);
    free(ptr);
}