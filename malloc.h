#include <stdlib.h>
#include "graph.h"

#define MyMalloc(S) MyMalloc_Rep( (S), __FILE__, __LINE__ )
#define MyFree(P) MyFree_Rep( (P), __FILE__, __LINE__ )

void *MyMalloc_Rep( size_t sz, const char *pcFileName, int nLine );
void MyFree_Rep( void *ptr, const char *pcFileName, int nLine );
void resizeGraph(Graph *G, int n, int m);