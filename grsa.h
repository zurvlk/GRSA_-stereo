#include "graph.h"
#include "bmp.h"

int function;
int nc2(int n);
double theta(double n, double T);
double p(int *label, int height, int width);
double energy(Graph *G, int *label, int *I, double T, int lamda);
double energy_str(Graph *G, int *label, double T, int lamda, int width, int *I_left, int *I_right);
double pairwise(double i, double j, double T, int lamda);
double data(int i, int label);
int gen_submodular_subsets(int label_size, int range_size, int **ls);
void set_edge(Graph *G, int height, int width, int *ls, int *label, int *I, double T, int lamda);
void set_edge_str(Graph *G, int height, int width, int *ls, int *label, double T, int lamda, int *I_left, int *I_right);
int is_convex(int i, int j, double T);
int isin_array(int *ls, int target);

int make_label_index(Graph *G, int *label, int *label_index, int alpha, int beta);
void set_single_edges(Graph *G, int height, int width);
int cmparray(int *array1, int *array2, int size);
void cpyarray(int *terget, int *source, int size);