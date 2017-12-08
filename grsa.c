#include "bmp.h"
#include "grsa.h"
#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <string.h>
// #include <math.h>

#define INF DBL_MAX
#define _CONSIDE_ALL_PAIRS_ 0

double dabs(double a, double b) {
    return a - b > 0 ? a - b : b - a;
}

double abss(double a) {
    return a > 0 ? a : - a;
}
int nc2(int n) {
    return n * (n - 1) / 2;
}

double fmin(double i, double j) {
    return i < j ? i : j;
}

double fmin3(double i, double j, double k) {
    return fmin(fmin(i, j), k);
}

double fmax(double i, double j) {
    return i > j ? i : j;
}
double theta(double n, double T) {


    if(!function) return fmin(n > 0 ? n : -n, T);
    else return fmin(n * n, T);
    // n = abss(n);
    // T = INF;
    // if (n <= 2) return 1.5 * n;
    // else if (n >= 2 && n <= 3) return 1.5 * 2;
    // else return 1.5 * (n - 1);
}

// リスト内にtargetが存在 1
int isin_array(int *ls, int target) {
    for (int i = 1; i <= ls[0]; i++) {
        if (ls[i] > target) break;
        if (ls[i] == target) return 1;
    }
    return 0;
}

int is_convex(int i, int j, double T) {
    if (theta(j, T) - theta(j - 1, T) >= (theta(j - 1, T) - theta(i, T)) / (j - 1 - i) &&
        i * (theta(i + 1, T) - theta(i, T)) >= theta(i, T) - theta(0, T) &&
        theta(i, T) - theta(0, T) >= 0) {
            return 1;
        }

    else return 0;
}

// 2つの配列の値がすべて同一 0
// 2つの配列の値が異なる 1
int cmparray(int *array1, int *array2, int size) {
    int i;
    for (i = 1; i <= size; i++) {
        if(array1[i] != array2[i]) return 1;
    }
    return 0;
}

void cpyarray(int *terget, int *source, int size) {
    int i;

    for (i = 1; i <= size; i++) {
        terget[i] = source[i];
    }
}

int mpair(int i, int j, int label_size) {
    int res = 0, c;
    int temp = i;
    if (i > j) {
        i = j;
        j = temp;
    }

    if(i == 0) return j;
    c = 0;
    while (c < i) {
        res += label_size - c;
        c++;
    }
    // printf("res: %d j: %d\n", res, j);
    res += j - 2 * c;
    return res;
}

int gen_submodular_subsets(int label_size, int range_size, Subsets *ss) {
    int i, j, k, l, m, n, prev, rs2, current_array, large_array, size;

    int ccvex = 0, label_max = label_size - 1;
    int convex[10][2];

    int **temp;

    // 候補区間抽出
    i = 0;
    j = 1;
    while (j < label_size) {
        j++;
        // printf("i: %d, j: %d isc: %f %f\n", i, j, theta(j, T) - theta(j - 1, T) , (theta(j - 1, T) - theta(i, T)) / (j - 1 - i));
        if (j > 1  && is_convex(i, j, ss->T)) {
            prev = 1;
            if (j == label_size) {
                convex[ccvex][0] = i;
                convex[ccvex][1] = j - 1;
                ccvex++;
            }
            // printf("%d - > %d\n", i, j);
        } else if (i != j - 1 && prev) {
            // printf("%d %d\n", i, j);
            convex[ccvex][0] = i;
            convex[ccvex][1] = j - 1;
            i = j - 1;
            ccvex++;
            prev = 0;
        } else {
            i = j - 1;
            prev = 0;
        }
        if (ccvex > 9) break;

    }
    for (i = 0; i < ccvex; i++) {
        printf("T : %.0d 候補区間: %d --> %d\n", ss->T, convex[i][0], convex[i][1]);
    }

    // 全てのラベルのペアを列挙
    if ((ss->pairs = (int **)malloc(sizeof(int*) * (nc2(label_size) + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main()->pairs]\n");
        exit(EXIT_FAILURE);
    }

    k = 1;
    for (i = 0; i < label_max; i++) {
        for (j = i + 1; j <= label_max; j++) {
            if ((ss->pairs[k] = (int*)malloc(sizeof(int) * 3)) == NULL) {
                fprintf(stderr, "Error!:malloc[main()-=>pairs]\n");
                exit(EXIT_FAILURE);
            }
            ss->pairs[k][0] = 0;
            ss->pairs[k][1] = i;
            ss->pairs[k][2] = j;
            k++;
        }
    }


    if ((ss->ls = (int **)malloc(sizeof(int*) * (nc2(label_size) + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main()->ls]\n");
        exit(EXIT_FAILURE);
    }

    ss->number = 0;
    if(label_size == range_size) {
        large_array = 1;
        ss->number = 1;
        if ((temp = (int **)realloc(ss->ls, 2 * sizeof(int *))) == NULL) {
            fprintf(stderr, "Error!:malloc[main()->ls]\n");
            exit(EXIT_FAILURE);
        } else {
            ss->ls = temp;
            temp = NULL;
        }

        if ((ss->ls[1] = (int *)malloc(sizeof(int) * (range_size + 1))) == NULL) {
            fprintf(stderr, "Error!:malloc[main()->ls]\n");
            exit(EXIT_FAILURE);
        }
        ss->ls[1][0] = range_size;
        for (i = 1; i <= range_size; i++) ss->ls[1][i] = i - 1;
        ss->number = 1;
    } else if (range_size == 2) {
        ss->ls = ss->pairs;
        for (i = 1; i <= nc2(label_size); i++) ss->ls[i][0] = 2;
        ss->number = nc2(label_size);

    } else {
        for (i = 0; i < ccvex; i++) {
            if (convex[i][0] == 0) {
                if (convex[i][1] > range_size) rs2 = range_size;
                else rs2 = convex[i][1];
                large_array = 1;
                j = rs2 - 1;
                do {
                    j += rs2 - 1;
                    large_array++;
                } while (j + rs2 - 1 < label_size);
                large_array++;

                printf("large_array %d\n", large_array);

                for (j = 1; j < large_array; j++) {
                    if ((ss->ls[j] = (int *)malloc(sizeof(int) * (rs2 + 1))) == NULL) {
                        fprintf(stderr, "Error!:malloc[main()->ls]\n");
                        exit(EXIT_FAILURE);
                    }
                    ss->number++;
                    ss->ls[j][0] = rs2;
                    if(j != 1) ss->ls[j][1] = ss->ls[j - 1][rs2];
                    else ss->ls[j][1] = 0;
                    for (k = 2; k <= rs2; k++) {
                        //segmentation fault
                        // printf("%d, %d\n", j, k);
                        ss->ls[j][k] = ss->ls[j][k - 1] + 1;
                        n = ss->ls[j][k];
                    }
                }
                size = label_size - n;
                if ((ss->ls[large_array] = (int *)malloc(sizeof(int) * size)) == NULL) {
                    fprintf(stderr, "Error!:malloc[main()->ls]\n");
                    exit(EXIT_FAILURE);
                }
                printf("size: %d\n", size);
                printf("n: %d\n", n);
                ss->ls[large_array][0] = size;
                ss->ls[large_array][1] = n;

                for (j = 2; ss->ls[large_array][j - 1] + 1 <= label_max; j++) {
                    ss->ls[large_array][j] = ss->ls[large_array][j - 1] + 1;
                    n = ss->ls[large_array][j];
                    printf("%d ", n);
                }
                ss->number++;
                printf("n: %d\n", n);

                for (current_array = 1; current_array <= large_array; current_array++) {
                    for (k = 0; k < ss->ls[current_array][0]; k++) {
                        for (l = k + 1; l <= ss->ls[current_array][0]; l++) {
                            m = mpair(ss->ls[current_array][k], ss->ls[current_array][l], label_size);
                            // printf("pair (%d, %d) == %d \n", ls[current_array][k], ls[current_array][l], m);

                            if (ss->pairs[m][1] == ss->ls[current_array][k] && ss->pairs[m][2] == ss->ls[current_array][l]) {
                                ss->pairs[m][0] = 1;
                            }
                        }
                    }
                }
            } else {
                large_array = 0;
                n = convex[i][1] / convex[i][0];
                if (n + 1 > range_size) n = range_size - 1;
                if (n >= 2) {
                    for (j = 0; j < convex[i][0]; j++) {
                        if (j + n * convex[i][0] <= label_max) large_array++;
                    }
                    for (j = 1; j <= large_array; j++) {
                        if ((ss->ls[ss->number + j] = (int *)malloc(sizeof(int) * (n + 2))) == NULL) {
                            fprintf(stderr, "Error!:malloc[main()->ls]\n");
                            exit(EXIT_FAILURE);
                        }
                        ss->ls[ss->number + j][0] = n + 1;
                        ss->ls[ss->number + j][1] = j - 1;
                        for (k = 2; k <= n + 1; k++) {
                            // printf("%d, %d\n",  ls[ss->number + j][k - 1] , convex[i][0]);
                            ss->ls[ss->number + j][k] = ss->ls[ss->number + j][k - 1] + convex[i][0];

                        }

                        for (k = 1; k < ss->ls[ss->number + j][0]; k++) {
                            for (l = k + 1; l <= ss->ls[ss->number + j][0]; l++) {
                                m = mpair(ss->ls[ss->number + j][k], ss->ls[ss->number + j][l], label_size);
                                // printf("p:%d %d\n", ls[ss->number + j][k], ls[ss->number + j][l]);
                                //ifいらないはず
                                if (ss->pairs[m][1] == ss->ls[ss->number + j][k] && ss->pairs[m][2] == ss->ls[ss->number + j][l]) {
                                    ss->pairs[m][0] = 1;
                                }
                            }
                        }
                    }
                    ss->number += large_array;
                }
            }
        }

    #if _CONSIDE_ALL_PAIRS_
        for (i = 1; i <= nc2(label_size); i++) {
            // printf("pairs[%d] = (%d, %d, %d)\n", i, pairs[i][0], pairs[i][1], pairs[i][2]);

            if (ss->pairs[i][0] == 0) {
                ss->number++;
                // printf("size: %d\n", size);
                if ((ss->ls[ss->number] = (int *)malloc(sizeof(int) * (3))) == NULL) {
                    fprintf(stderr, "Error!:malloc[main()->ls]\n");
                    exit(EXIT_FAILURE);
                }

                ss->ls[ss->number ][0] = 2;
                ss->ls[ss->number ][1] = ss->pairs[i][1];
                ss->ls[ss->number ][2] = ss->pairs[i][2];

            }
        }
    #endif

        if ((temp = (int **)realloc(ss->ls, ss->number * sizeof(int *))) == NULL) {
            fprintf(stderr, "Error!:malloc[main()->ls]\n");
            exit(EXIT_FAILURE);
        } else {
            ss->ls = temp;
            temp = NULL;
        }

    }
    printf("ss->number = %d\n", ss->number);
#if _OUTPUT_SUBMODULAR_SUBSETS_
    for (i = 1; i <= ss->number; i++) {
        for (j = 1; j <= ls[i][0]; j++) {
            printf("%d ", ls[i][j]);
        }
        printf("\n");
    }
#endif
    return ss->number;
}

void readStrBitmap(Image *image, char filename[], int scale) {
    int i, j, grids_node;
    char imgleft[100];
    char imgright[100];
    char imgtruth[100];

    strcpy(imgleft, filename);
    strcpy(imgright, filename);
    strcpy(imgtruth, filename);

    strcat(imgleft, "left.bmp");
    strcat(imgright, "right.bmp");
    strcat(imgtruth, "truth.bmp");

    ReadBmp(imgleft, &(image->raw_left));
    ReadBmp(imgright, &(image->raw_right));
    ReadBmp(imgtruth, &(image->truth));
    ReadBmp(imgtruth, &(image->output));

    image->width = image->raw_left.width;
    image->height = image->raw_left.height;
    image->scale = scale;
    grids_node = image->height * image->width;

    if(image->width != image->raw_right.width || image->height != image->raw_right.height) {
        fprintf(stderr, "Error %s と %s の解像度が異なります\n", imgleft, imgright);
        exit(EXIT_FAILURE);
    }

    Gray(&(image->raw_left), &(image->raw_left));
    Gray(&(image->raw_right), &(image->raw_right));
    Gray(&(image->truth), &(image->truth));


    if ((image->left = (int *)malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main()->image.left]\n");
        exit(EXIT_FAILURE);
    }
    if ((image->right = (int *)malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main()->image.right]\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i <  image->height; i++) {
        for (j = 0; j < image->width; j++) {
            image->left[i * image->width + j + 1] = image->raw_left.data[i][j].r / scale;
        }
    }
    for (i = 0; i <  image->height; i++) {
        for (j = 0; j < image->width; j++) {
            image->right[i * image->width + j + 1] = image->raw_right.data[i][j].r / scale;
        }
    }

}


double pairwise(double i, double j, double T, int lamda) {
    return lamda * theta(i - j, T);
}

double data(int i, int label) {
    return 1.0 * dabs(label, i);
}


double data_str(int i, int label, int width, int *I_left, int *I_right) {
    double data = 0;
    // return 1.0 * dabs(label, I_left[i]);
    //leftの中のものがどこにあるか
    if ((i - 1) / width == (i - label - 1) / width) {
        // data = (I_left[i] - I_right[i - label]) * (I_left[i] - I_right[i - label]);
        data = (I_left[i] - I_right[i - label]);
    }else data = INF;

    return abss(data);
    // return sqrt(data);

}



double between(double a, double b, double c){
    if((a <= b && b <= c) || (c <= b && b <= a))    return 0;
    else return 1;
}

double Dt(int x, Image *image, int i, int j) {
    double d_m, d_l, d_r, d_m2, d_l2, d_r2;
    double  I, I_1 = 0, I_2 = 0;

	if(x < 0 || x > image->label_max ) return INF;
    if (j - x < 0) return 10000;

    int flag = 1;
    while(flag){

	    d_m = image->raw_left.data[i][j].r;
	    d_m2 =  image->raw_right.data[i][j - x].r;

	    if((j - x - 1 >= 0) && (j - x + 1 <= image->raw_right.width - 1)){

	        d_l2 = (image->raw_right.data[i][j - x].r + image->raw_right.data[i][j - x - 1].r) / 2.0;
            d_r2 = (image->raw_right.data[i][j - x].r + image->raw_right.data[i][j - x + 1].r) / 2.0;
	        if(between(d_l2, d_m, d_m2) == 0)    break;
	        if(between(d_r2, d_m, d_m2) == 0)    break;
	        I_1 = fmin3(abss(d_m - d_l2), abss(d_m - d_m2), abss(d_m - d_r2));              /*特に制約なし*/

	    }else if(j - x - 1 < 0){

	        d_r2 = (image->raw_right.data[i][j - x].r + image->raw_right.data[i][j - x + 1].r) / 2.0;
	        if(between(d_m2, d_m, d_r2) == 0)    break;
	        I_1 = fmin(abss(d_m - d_m2), abss(d_m - d_r2));                                /*Rの左端が出る*/

	    }else{

	        d_l2 = (image->raw_right.data[i][j - x].r + image->raw_right.data[i][j - x - 1].r) / 2.0;
	        if(between(d_l2, d_m, d_m2) == 0)    break;
	        I_1 = fmin(abss(d_m - d_l2), abss(d_m - d_m2));                                /*Rの右端が出る*/
	    }

	    if((j != 0) && (j != image->raw_left.width - 1)){

	        d_l = (image->raw_left.data[i][j].r + image->raw_left.data[i][j - 1].r) / 2.0;
	        d_r = (image->raw_left.data[i][j].r + image->raw_left.data[i][j + 1].r) / 2.0;
	        if(between(d_l, d_m2, d_m) == 0)    break;
	        if(between(d_r, d_m2, d_m) == 0)    break;
	        I_2 = fmin3(abss(d_l - d_m2), abss(d_m - d_m2), abss(d_r - d_m2));              /*特に制約なし*/

	    }else if(j == 0){

	        d_r = (image->raw_left.data[i][j].r + image->raw_left.data[i][j + 1].r) / 2.0;
	        if(between(d_m, d_m2, d_r) == 0)    break;
	        I_2 = fmin(abss(d_m - d_m2), abss(d_r - d_m2));                                /*Lの左端が出る*/

	    }else{

	        d_l = (image->raw_left.data[i][j].r + image->raw_left.data[i][j - 1].r) / 2.0;
	        if(between(d_l, d_m2, d_m) == 0)    break;
	        I_2 = fmin(abss(d_l - d_m2), abss(d_m - d_m2));                                /*Lの右端が出る*/
	    }
	    flag = 0;
	}
    I = fmin(I_1 / (double)image->scale, I_2 / (double)image->scale);
    I = fmin(I, 20);
    return I;
    //I *= I;
    //return I * I;                                   /*Dを2乗する*/
    //return C_D * (I);
}

double energy(Graph *G, int *label, int *I, double T, int lamda) {
    int i;
    double energy = 0;
    //* Dterm
    for (i = 1; i <= G->n - 2; i++) {
        energy += data(I[i], label[i]);
    }
    // */
    // Vterm
    for (i = 1; i <= G->m - 2 * (G->n - 2); i++) {
        energy += pairwise(label[G->tail[i]], label[G->head[i]], T, lamda);
    }
    return energy;
}

double energy_str(Graph *G, int *label,  double T, int lamda, Image image) {
    int i;
    double energy = 0;
    //* Dterm
    //*
    for (i = 1; i <= G->n - 2; i++) {
        energy += data_str(i, label[i], image.width, image.left, image.right);
        // energy += data(I_left[i], label[i]);
    }
    /*/
    for(int i = 1; i < G->n - 1; i++){
        energy += Dt(label[i], &image, i / image.width, i % image.width);
    }
    // */
    // Vterm
    for (i = 1; i <= G->m - 2 * (G->n - 2); i++) {
        energy += pairwise(label[G->tail[i]], label[G->head[i]], T, lamda);
    }
    return energy;
}

double e_cost(int i, int j, double T) {
    double cost = 0.5 * (theta(i - j + 1, T)- 2 * theta(i - j, T) + theta(i - j - 1, T));
    return cost;
}

//枝を作る関数makeedge(グラフ,高さ,幅)
void set_single_edges(Graph *G, int height, int width) {
    int i, j, edge_count;
    int tail, head, source, sink;

    source = G->n - 1;
    sink = source + 1;
    setSource(G, source);
    setSink(G, sink);

    edge_count = 1;
    //点と点の間の枝（横）
    for (i = 1; i < height + 1; i++) {
        for (j = 1; j < width; j++) {
            tail =  (i - 1) * width + j;
            head =  tail + 1;
            setEdge(G, edge_count, tail, head, 0);
            edge_count++;
        }
    }

    //点と点の間の枝（縦）
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            tail = (i - 1) * width + j;
            head = tail + width;
            setEdge(G, edge_count, tail, head, 0);
            edge_count++;
        }
    }

    //sourceと点の間の枝
    for (i = 1; i < height * width + 1; i++){
       setEdge(G, edge_count, G->src, i, 0);
       edge_count++;
    }

    //点とsinkの間の枝
    for (i = 1; i < height * width + 1; i++){
       setEdge(G, edge_count, i, G->sink, 0);
       edge_count++;
    }
    return;
}

int make_label_index(Graph *G, int *label, int *label_index, int alpha, int beta) {
    int i, arraysize;
    for (i = 1; i <= G->n - 2; i++) label_index[i] = 0;
    arraysize = 1;
    for (i = 1; i <= G->n - 2; i++) {
        if (label[i] <= beta && label[i] >= alpha) {
            label_index[arraysize] = i;
            arraysize++;
        }
    }
    return arraysize;
}

double phi (int i, int j, int *ls, double T, int lamda) {
    double p = 0;
    if(j > i) return 0;
    if(1 < j && j <= i) {
        p = lamda * (theta(ls[i] - ls[j - 1], T) - theta(ls[i] - ls[j], T) - theta(ls[i - 1] - ls[j - 1], T) + theta(ls[i - 1] - ls[j], T));
        if(i == j) p *= 0.5;
    }
    return p;
}

double nn(int i, int label, int *ls, double T, int lamda) {
    double p = 0;
    p = pairwise(ls[i], label, T, lamda);
    return p;
}

double near_nodes(int i, int j, int height, int width, int *ls, int *label, int *isin, double T, int lamda) {
    int grids_node = height * width;

    double nnp_total = 0;

    if (i >= width + 1) {
        // 画素が一番上の行に存在しないとき(iの上が空白でないとき)
        if (!isin[i - width]){
            // iの上の点がLs内に含まれない
            nnp_total += pairwise(ls[j], label[i - width], T, lamda) ;
        }
    }

    if (i <= grids_node - width) {
        // 画素が一番下の行に存在しないとき(iの下が空白でないとき)
        if (!isin[i + width]){
            // iの下の点がLs内に含まれない
            nnp_total += pairwise(ls[j], label[i + width], T, lamda) ;
        }
    }

    if ((i % width) != 1) {
        // 画素が一番左の列に存在しないとき(iの左が空白でないとき)
        if (!isin[i - 1]){
            // iの左の点がLs内に含まれない
            nnp_total += pairwise(ls[j], label[i - 1], T, lamda) ;
        }
    }

    if ((i % width) != 0) {
        // 画素が一番右の列に存在しないとき(iの右が空白でないとき)
        if (!isin[i + 1]){
            // iの右の点がLs内に含まれない
            nnp_total += pairwise(ls[j], label[i + 1], T, lamda) ;
        }
    }
    return nnp_total;
}

// set_edge for grsa
int set_edge(Graph *G, int *ls, int *label, double T, int lamda, Image image) {
    int i, j, k, l;
    int tail, head, t_base, h_base, grids_node, source, sink, edge_count, current_edge;
    int s2i_begin, i2t_begin, depth_begin;
    int *isin;
    double *min;

    // 格子部分1階層分の点数合計
    grids_node = image.height * image.width;

    if (((min = (double *) malloc(sizeof(double) * G->n))) == NULL) {
        fprintf(stderr, "set_all_edge(): ERROR [min = malloc()]\n");
        exit (EXIT_FAILURE);
    }
    // min[i]の全てにINFを設定
    for (i = 0; i < G->n; i++) min[i] = INF;

    if (((isin = (int *) malloc(sizeof(int) * (grids_node + 1)))) == NULL) {
        fprintf(stderr, "set_all_edge(): ERROR [isin = malloc()]\n");
        exit (EXIT_FAILURE);
    }
    // min[i]の全てにINFを設定
    for (i = 1; i <= grids_node; i++) isin[i] = isin_array(ls, label[i]);

    for (i = 1; i < G->n; i++) G->capa[i] = 0;
    source = grids_node * ls[0] + 1;
    sink = source + 1;

    setSource(G, source);
    setSink(G, sink);

    edge_count = 1;
    // source->i1
    s2i_begin = edge_count;
    for (i = 1; i <= grids_node; i++) {
        setEdge(G, edge_count, source, i, 0);
        if(isin[i]) {
            G->capa[edge_count] = INF;
        }
        if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
        edge_count++;
    }

    // depth
    depth_begin = edge_count;
    for (i = 1; i <= grids_node; i++) {
        tail = i;
        head = i;
        for (j = 1; j < ls[0]; j++) {
            head = head + grids_node;
            setEdge(G, edge_count, tail, head, 0);
            if(isin[i]) {
                // G->capa[edge_count] = data(I_left[i],ls[j]) + nnp_4_grsa(i, j, height, width, ls, label, T, lamda);
                G->capa[edge_count] = data_str(i, ls[j], image.width, image.left, image.right) + near_nodes(i, j, image.height, image.width, ls, label, isin, T, lamda);
                // G->capa[edge_count] = Dt(ls[j], &image, i / image.width, i % image.width) + near_nodes(i, j, image.height, image.width, ls, label, isin, T, lamda);
            }

            if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
            edge_count++;
            tail = head;
        }
    }
    // ik->sink
    i2t_begin = edge_count;
    // rk = r(beta , label_size, grids_edge);
    for (i = 1; i <= grids_node; i++) {
        setEdge(G, edge_count, i + grids_node * (ls[0] - 1), sink, 0);
        if(isin[i]) {
            // G->capa[edge_count] = data(I_left[i], ls[ls[0]]) + nnp_4_grsa(i, ls[0], height, width, ls, label, T, lamda);
            G->capa[edge_count] = data_str(i, ls[ls[0]], image.width, image.left, image.right) + near_nodes(i, j, image.height, image.width , ls, label, isin, T, lamda);
            // G->capa[edge_count] = Dt(ls[ls[0]], &image, i / image.width, i % image.width) + near_nodes(i, j, image.height, image.width , ls, label, isin, T, lamda);
        }

        if (min[i] > G->capa[edge_count]) min[i] = G->capa[edge_count];
        edge_count++;
    }


    // reverce edge
    // depth
    for (i = 1; i <= grids_node; i++) {
        tail = i;
        head = i;
        for (j = 1; j < ls[0]; j++) {
            head = head + grids_node;
            setEdge(G, edge_count, head, tail, 0);
            if(isin[i]) {
                G->capa[edge_count] = INF;
            }
            edge_count++;
            tail = head;
        }
    }


    // new codes
    // horizonal
    for (i = 1; i <= image.height; i++) {
        for (j = 1; j < image.width; j++) {
            t_base =  (i - 1) * image.width + j;
            h_base =  t_base + 1;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(isin[t_base] && isin[h_base]) {
                        // head, tail in label_index
                        G->capa[edge_count] = phi(k + 1, l + 1, ls, T, lamda);
                    }
                    edge_count++;
                }
            }
        }
    }

    // vertical
    for (i = 1; i < image.height ; i++){
        for (j = 1; j < image.width + 1; j++) {
            t_base = (i - 1) * image.width + j;
            h_base = t_base + image.width;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(isin[t_base] && isin[h_base]) {
                        G->capa[edge_count] = phi(k + 1, l + 1, ls, T, lamda);
                        // head, tail in label_index
                    }
                    edge_count++;
                }
            }
        }
    }



    // new codes
    // horizonal
    for (i = 1; i <= image.height; i++) {
        for (j = 1; j < image.width; j++) {
            t_base =  (i - 1) * image.width + j;
            h_base =  t_base + 1;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    if(isin[t_base] && isin[h_base]) {
                        G->capa[edge_count] = phi(l + 1, k + 1, ls, T, lamda);
                    }
                    edge_count++;
                }
            }
        }
    }

    // vertical
    for (i = 1; i < image.height ; i++){
        for (j = 1; j < image.width + 1; j++) {
            t_base = (i - 1) * image.width + j;
            h_base = t_base + image.width;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    if(isin_array(ls, label[t_base]) && isin_array(ls, label[h_base])) {
                        G->capa[edge_count] = phi(l + 1, k + 1, ls, T, lamda);
                    }
                    edge_count++;
                }
            }
        }
    }


    //  s->tの一連の枝から定数値を引く処理
    for (i = s2i_begin; i <= grids_node; i++) {
        G->capa[i] -= min[i];
    }
    for (i = 1; i <= grids_node; i++) {
        G->capa[i + i2t_begin - 1] -= min[i];
    }
    current_edge = depth_begin;
    for (i = 1; i <= grids_node; i++) {
        for (j = 1; j < ls[0]; j++) {
            G->capa[current_edge] -=min[i];
            current_edge++;
        }
    }

    free(min);
    free(isin);

    // printf("total edge : %d\n", edge_count - 1);
    return edge_count - 1;
}

double err_rate(img output, Image image) {
    int i, error_count = 0;
    double err;
    //     // Gray(&truth, &truth);
    //     Gray(&(image.output), &(image.output));

    if (image.truth.data[0][0].r) {
        for(i = 1; i <= (image.output.height) * (image.output.width); i++) {
            if (abs(image.output.data[(i - 1) / image.output.width][(i - 1) % image.output.width].r - image.truth.data[(i - 1) / image.truth.width][(i - 1) % image.truth.width].r )
                >= image.scale + 1) {
                error_count++;
            }
        }
    } else {
        for(i = 1; i <= (image.output.height) * (image.output.width); i++) {
            if ((i - 1) / image.output.width >= image.scale && (i - 1) % image.output.width >= image.scale &&
                (i - 1) / image.output.width <= image.output.height - image.scale && (i - 1) % image.output.width <= image.output.width - image.scale) {
                if (abs(image.output.data[(i - 1) / image.output.width][(i - 1) % image.output.width].r - image.truth.data[(i - 1) / image.truth.width][(i - 1) % image.truth.width].r )
                    >= image.scale + 1) {
                    error_count++;
                }
            }
        }
    }

    err = 100 * error_count / (double)(image.truth.height * image.truth.width);
    return err;
}
