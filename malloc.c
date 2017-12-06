#include "malloc.h"
#include "graph.h"
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

/* グラフ構造体が利用する配列のメモリを確保する関数 */
void resizeGraph(Graph *G, int n, int m) {

    int *temp;
    double *temp2;
	if (n < 1) {
		fprintf(stderr, "newGraph(): ERROR [2nd argument]\n");
		exit(EXIT_FAILURE);
	}
	if (m < 0) {
		fprintf(stderr, "newGraph(): ERROR [3rd argument]\n");
		exit(EXIT_FAILURE);
	}

	G->nmax = n; G->mmax = m;
	G->n = n++; G->m = m++;
	G->src = 0; G->sink = 0;
    
    if ((temp = (int *)realloc(G->tail, sizeof(int) * m)) == NULL) {
        fprintf(stderr, "resizeGraph(): ERROR [realloc(G->tail)]\n");
		exit(EXIT_FAILURE);
    } else {
        G->tail = temp;
        temp = NULL;
    }

	if ((temp = (int *)realloc(G->head, sizeof(int) * m)) == NULL) {
        fprintf(stderr, "resizeGraph(): ERROR [realloc(G->head)]\n");
		exit(EXIT_FAILURE);
    } else {
        G->head = temp;
        temp = NULL;
    }

	if ((temp2 = (double *)realloc(G->capa, sizeof(double) * m)) == NULL) {
        fprintf(stderr, "resizeGraph(): ERROR [realloc(G->capa)]\n");
		exit(EXIT_FAILURE);
    } else {
        G->capa = temp2;
        temp2 = NULL;
    }

    if ((temp = (int *)realloc(G->first, sizeof(int) * n)) == NULL) {
        fprintf(stderr, "resizeGraph(): ERROR [realloc(G->first)]\n");
		exit(EXIT_FAILURE);
    } else {
        G->first = temp;
        temp = NULL;
    }
	
	if ((temp = (int *)realloc(G->next, sizeof(int) * m)) == NULL) {
        fprintf(stderr, "resizeGraph(): ERROR [realloc(G->next)]\n");
		exit(EXIT_FAILURE);
    } else {
        G->next = temp;
        temp = NULL;
    }

	if ((temp = (int *)realloc(G->prev, sizeof(int) * m)) == NULL) {
        fprintf(stderr, "resizeGraph(): ERROR [realloc(G->prev)]\n");
		exit(EXIT_FAILURE);
    } else {
        G->prev = temp;
        temp = NULL;
    }

	G->revf = NULL; G->revn = NULL; G->revp = NULL;
}

