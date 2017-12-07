#include "graph.h"
#include "bmp.h"

typedef struct __INPUT_BITMAPS__ {
    img raw_left;
    img raw_right;
    img truth;
    img output;
    int *left;
    int *right;
    int width;
    int height;
    int label_max;
    int scale;
}Image;

typedef struct __SUBMODULAR_SUBSETS__ {
    int **ls;
    int **pairs;
    int number;
    int T;
}Subsets;

int function;
void readStrBitmap(Image *image, char filename[], int scale);
int gen_submodular_subsets(int label_size, int range_size, Subsets *ss);
int nc2(int n);
double theta(double n, double T);
double p(int *label, int height, int width);
double energy(Graph *G, int *label, int *I, double T, int lamda);
double energy_str(Graph *G, int *label, double T, int lamda, Image image);
double pairwise(double i, double j, double T, int lamda);
double data(int i, int label);
int set_edge(Graph *G, int *ls, int *label, double T, int lamda, Image image);
int is_convex(int i, int j, double T);
int isin_array(int *ls, int target);

int make_label_index(Graph *G, int *label, int *label_index, int alpha, int beta);
void set_single_edges(Graph *G, int height, int width);
int cmparray(int *array1, int *array2, int size);
void cpyarray(int *terget, int *source, int size);
double err_rate(img output, Image image);
