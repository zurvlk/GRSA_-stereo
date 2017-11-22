#include "graph.h"
#include "bmp.h"

int function;
int nc2(int n);
double theta(double n, double T);
double p(int *label, int height, int width);
double energy(Graph *G, int *label, int *I, double T);
double pairwise(double i, double j, double T);
double data(int i, int label);
int gen_submodular_subsets(int label_size, int range_size, int **ls);
void set_edge(Graph *G, int height, int width, int *ls, int *label, int *I, double T);
int is_convex(int i, int j, double T);
int isin_array(int *ls, int target);

int make_label_index(Graph *G, int *label, int *label_index, int alpha, int beta);
void set_single_edges(Graph *G, int height, int width);
int cmparray(int *array1, int *array2, int size);
void cpyarray(int *terget, int *source, int size);