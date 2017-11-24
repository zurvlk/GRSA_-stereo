#include "bmp.h"
#include "grsa.h"
#include "graph.h"
#include <stdlib.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

#define INF DBL_MAX

double dabs(double a, double b) {
    return a - b > 0 ? a - b : b - a;
}

int nc2(int n) {
    return n * (n - 1) / 2;
}

double fmin(double i, double j) {
    return i < j ? i : j;
}

double fmax(double i, double j) {
    return i > j ? i : j;
}
double theta(double n, double T) {

    if(!function) return fmin(n > 0 ? n : -n, T);
    else return fmin(n * n, T);
    // return fmin(n >= 0 ? n : -n, fmax(n >= 0 ? n - 5 : -n - 5, T));

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


int gen_submodular_subsets(int label_size, int range_size, int **ls) {
    int i, j, k, l, size, large_array, total_ss_count;
    int label_max = label_size - 1;
    if (label_size != range_size) {
        large_array = 0;
        if (range_size > 2) {
            i = 0;
            do {
                i += range_size - 1;
                large_array++;
            } while (i + range_size < label_size);
            large_array++;
            size = label_max - i + 1;
            total_ss_count = large_array + nc2(label_size) - (large_array - 1) * nc2(range_size) - nc2(size);
        } else total_ss_count = nc2(label_size);
        
        
        // large_array = label_size / (range_size - 1) - 1;
        
        
        printf("size : %d large_array %d\n", size, large_array);
        if ((ls = (int **)malloc(sizeof(int*) * (total_ss_count + 1))) == NULL) {
            fprintf(stderr, "Error!:malloc[main()->ls]\n");
            exit(EXIT_FAILURE);
        }

        if (range_size > 2) {
            for (i = 1; i < large_array; i++) {
                if ((ls[i] = (int *)malloc(sizeof(int) * (range_size + 1))) == NULL) {
                    fprintf(stderr, "Error!:malloc[main()->ls]\n");
                    exit(EXIT_FAILURE);
                }
                ls[i][0] = range_size;
                if(i != 1) ls[i][1] = ls[i - 1][range_size];
                else ls[i][1] = 0;
                // printf("%d ", ls[i][1]);
                for (j = 2; j <= range_size; j++) {
                    ls[i][j] = ls[i][j - 1] + 1;
                    k = ls[i][j];
                    // ls[i][j] = (i - 1) * range_size + j;
                    // printf("%d ", ls[i][j]);
                }
                // printf("\n");
            }
        
            if ((ls[large_array] = (int *)malloc(sizeof(int) * size)) == NULL) {
                fprintf(stderr, "Error!:malloc[main()->ls]\n");
                exit(EXIT_FAILURE);
            }
            ls[large_array][0] = size;
            ls[large_array][1] = k;
            // printf("%d ", ls[large_array][1]);
            for (j = 2; ls[large_array][j - 1] + 1 <= label_max; j++) {
                    ls[large_array][j] = ls[large_array][j - 1] + 1;
                    // ls[i][j] = (i - 1) * range_size + j;
                    // printf("%d ", ls[large_array][j]);
            }
            // printf("\n\n");
            for (i = large_array + 1; i <= total_ss_count; i++) {
                if ((ls[i] = (int *)malloc(sizeof(int) * (3))) == NULL) {
                    fprintf(stderr, "Error!:malloc[main()->ls]\n");
                    exit(EXIT_FAILURE);
                }
                ls[i][0] = 2;
            }
            i = 0;
            j = range_size;
            k = large_array + 1;
            l = range_size - 1;
            while(i < label_max - (ls[large_array][0] - 1)) {
                ls[k][1] = i;
                ls[k][2] = j;
                // printf("%d %d \n", ls[k][1], ls[k][2]);
                k++;
                if (j == label_max) {
                    i++;
                    j = l + 1;
                    if (i == l - 1) l += range_size - 1;
                }
                else j++;
            }
        } else {

            for (i = 1; i <= total_ss_count; i++) {
                if ((ls[i] = (int *)malloc(sizeof(int) * (3))) == NULL) {
                    fprintf(stderr, "Error!:malloc[main()->ls]\n");
                    exit(EXIT_FAILURE);
                }
                ls[i][0] = 2;
            }

            i = 0;
            j = 1;
            k = large_array + 1;
            while(i < label_max) {
                ls[k][1] = i;
                ls[k][2] = j;
                // printf("%d %d \n", ls[k][1], ls[k][2]);
                k++;
                if (j == label_max) {
                    i++;
                    j = i + 1;
                }
                else j++;
            }
        }
    } else {
        large_array = 1;
        total_ss_count = 1;
        if ((ls = (int **)malloc(sizeof(int*) * (total_ss_count + 1))) == NULL) {
            fprintf(stderr, "Error!:malloc[main()->ls]\n");
            exit(EXIT_FAILURE);
        }

        if ((ls[1] = (int *)malloc(sizeof(int) * (range_size + 1))) == NULL) {
            fprintf(stderr, "Error!:malloc[main()->ls]\n");
            exit(EXIT_FAILURE);
        }
        ls[1][0] = range_size;
        for (i = 1; i <= range_size; i++) ls[1][i] = i - 1;
    }

    for (i = 1; i <= total_ss_count; i++) {
        for (j = 1; j <= ls[i][0]; j++) {
            printf("%d ", ls[i][j]);
        }
        printf("\n");
    }
    return total_ss_count;
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
        data = (I_left[i] - I_right[i - label]) * (I_left[i] - I_right[i - label]);
    }else data = INF;
    

    return sqrt(data);
    
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

double energy_str(Graph *G, int *label,  double T, int lamda, int width, int *I_left, int *I_right) {
    int i;
    double energy = 0;
    //* Dterm
    for (i = 1; i <= G->n - 2; i++) {
        energy += data_str(i, label[i], width, I_left, I_right);
        // energy += data(I_left[i], label[i]);
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

double phi (int i, int j, int *ls, double T) {
    double p = 0;
    if(j > i) return 0;
    if(1 < j && j <= i) {
        p = theta(ls[i] - ls[j - 1], T) - theta(ls[i] - ls[j], T) - theta(ls[i - 1] - ls[j - 1], T) + theta(ls[i - 1] - ls[j], T);
        if(i == j) p *= 0.5;
    }
    return p;
}

double nn(int i, int label, int *ls, double T, int lamda) {
    double p = 0;
    p = pairwise(ls[i], label, T, lamda);
    return p;
}

int nnp_4_grsa(int i, int j, int height, int width, int *ls, int *label, double T, int lamda) {
    int grids_node = height * width;

    double nnp_total = 0;

    if (i >= width + 1) {
        // 画素が一番上の行に存在しないとき(iの上が空白でないとき)
        if (!isin_array(ls, label[i - width])){
            // iの上の点がLs内に含まれない
            nnp_total += pairwise(ls[j], label[i - width], T, lamda) ;
        }
    }
    
    if (i <= grids_node - width) {
        // 画素が一番下の行に存在しないとき(iの下が空白でないとき)
        if (!isin_array(ls, label[i + width])){
            // iの下の点がLs内に含まれない
            nnp_total += pairwise(ls[j], label[i + width], T, lamda) ;
        }
    }

    if ((i % width) != 1) {
        // 画素が一番左の列に存在しないとき(iの左が空白でないとき)
        if (!isin_array(ls, label[i - 1])){
            // iの左の点がLs内に含まれない
            nnp_total += pairwise(ls[j], label[i - 1], T, lamda) ;
        }
    }

    if ((i % width) != 0) {
        // 画素が一番右の列に存在しないとき(iの右が空白でないとき)
        if (!isin_array(ls, label[i + 1])){
            // iの右の点がLs内に含まれない
            nnp_total += pairwise(ls[j], label[i + 1], T, lamda) ;
        }
    }
    return nnp_total;
}

// set_edge for grsa
void set_edge(Graph *G, int height, int width, int *ls, int *label, int *I, double T, int lamda) {
    int i, j, k, l;
    int tail, head, t_base, h_base, grids_node, source, sink, edge_count, current_edge;
    int s2i_begin, i2t_begin, depth_begin;
    double *min;


    if (((min = (double *) malloc(sizeof(double) * G->n))) == NULL) {
        fprintf(stderr, "set_all_edge(): ERROR [min = malloc()]\n");
        exit (EXIT_FAILURE);
    }
    // min[i]の全てにINFを設定
    for (i = 0; i < G->n; i++) min[i] = INF;

    // 格子部分1階層分の点数合計
    grids_node = height * width;

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
        if(isin_array(ls, label[i])) {
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
            if(isin_array(ls , label[i])) {
                G->capa[edge_count] = data(I[i],ls[j]) + nnp_4_grsa(i, j, height, width, ls, label, T, lamda);
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
        if(isin_array(ls , label[i])) {
            G->capa[edge_count] = data(I[i], ls[ls[0]]) + nnp_4_grsa(i, ls[0], height, width, ls, label, T, lamda);
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
            if(isin_array(ls, label[i])) {
                G->capa[edge_count] = INF;
            }
            edge_count++;
            tail = head;
        }
    }

    

    // new codes
    // horizonal
    for (i = 1; i <= height; i++) {
        for (j = 1; j < width; j++) {
            t_base =  (i - 1) * width + j;
            h_base =  t_base + 1;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(isin_array(ls, label[t_base]) && isin_array(ls, label[h_base])) {
                        // head, tail in label_index
                        G->capa[edge_count] = phi(k + 1, l + 1, ls, T);
                    }
                    edge_count++;
                }
            }
        }
    }
    
    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            t_base = (i - 1) * width + j;
            h_base = t_base + width;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(isin_array(ls, label[t_base]) && isin_array(ls, label[h_base])) {
                        G->capa[edge_count] = phi(k + 1, l + 1, ls, T);
                        // head, tail in label_index
                    }
                    edge_count++;
                }
            }
        }
    }




    // new codes
    // horizonal
    for (i = 1; i <= height; i++) {
        for (j = 1; j < width; j++) {
            t_base =  (i - 1) * width + j;
            h_base =  t_base + 1;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    if(isin_array(ls, label[t_base]) && isin_array(ls, label[h_base])) {
                        G->capa[edge_count] = phi(l + 1, k + 1, ls, T);
                    }
                    edge_count++;
                }
            }
        }
    }

    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            t_base = (i - 1) * width + j;
            h_base = t_base + width;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    if(isin_array(ls, label[t_base]) && isin_array(ls, label[h_base])) {
                        G->capa[edge_count] = phi(l + 1, k + 1, ls, T);
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
    
    // printf("total edge : %d\n", edge_count - 1);
    return;
}

void set_edge_str(Graph *G, int height, int width, int *ls, int *label, double T, int lamda, int *I_left, int *I_right) {
    int i, j, k, l;
    int tail, head, t_base, h_base, grids_node, source, sink, edge_count, current_edge;
    int s2i_begin, i2t_begin, depth_begin;
    double *min;


    if (((min = (double *) malloc(sizeof(double) * G->n))) == NULL) {
        fprintf(stderr, "set_all_edge(): ERROR [min = malloc()]\n");
        exit (EXIT_FAILURE);
    }
    // min[i]の全てにINFを設定
    for (i = 0; i < G->n; i++) min[i] = INF;

    // 格子部分1階層分の点数合計
    grids_node = height * width;

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
        if(isin_array(ls, label[i])) {
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
            if(isin_array(ls , label[i])) {
                // G->capa[edge_count] = data(I_left[i],ls[j]) + nnp_4_grsa(i, j, height, width, ls, label, T, lamda);
                G->capa[edge_count] = data_str(i, ls[j], width, I_left, I_right) + nnp_4_grsa(i, j, height, width, ls, label, T, lamda);
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
        if(isin_array(ls , label[i])) {
            // G->capa[edge_count] = data(I_left[i], ls[ls[0]]) + nnp_4_grsa(i, ls[0], height, width, ls, label, T, lamda);
            G->capa[edge_count] = data_str(i, ls[ls[0]], width, I_left, I_right) + nnp_4_grsa(i, ls[0], height, width, ls, label, T, lamda);
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
            if(isin_array(ls, label[i])) {
                G->capa[edge_count] = INF;
            }
            edge_count++;
            tail = head;
        }
    }


    // new codes
    // horizonal
    for (i = 1; i <= height; i++) {
        for (j = 1; j < width; j++) {
            t_base =  (i - 1) * width + j;
            h_base =  t_base + 1;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(isin_array(ls, label[t_base]) && isin_array(ls, label[h_base])) {
                        // head, tail in label_index
                        G->capa[edge_count] = phi(k + 1, l + 1, ls, T);
                    }
                    edge_count++;
                }
            }
        }
    }
    
    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            t_base = (i - 1) * width + j;
            h_base = t_base + width;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(isin_array(ls, label[t_base]) && isin_array(ls, label[h_base])) {
                        G->capa[edge_count] = phi(k + 1, l + 1, ls, T);
                        // head, tail in label_index
                    }
                    edge_count++;
                }
            }
        }
    }



    // new codes
    // horizonal
    for (i = 1; i <= height; i++) {
        for (j = 1; j < width; j++) {
            t_base =  (i - 1) * width + j;
            h_base =  t_base + 1;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    if(isin_array(ls, label[t_base]) && isin_array(ls, label[h_base])) {
                        G->capa[edge_count] = phi(l + 1, k + 1, ls, T);
                    }
                    edge_count++;
                }
            }
        }
    }

    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            t_base = (i - 1) * width + j;
            h_base = t_base + width;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    if(isin_array(ls, label[t_base]) && isin_array(ls, label[h_base])) {
                        G->capa[edge_count] = phi(l + 1, k + 1, ls, T);
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
    
    // printf("total edge : %d\n", edge_count - 1);
    return;
}
/**
void set_edge_str(Graph *G, int height, int width, int *ls, int *label, double T, int lamda, int *I_right, int *I_left) {
    int i, j, k, l;
    int tail, head, t_base, h_base, grids_node, source, sink, edge_count, current_edge;
    int s2i_begin, i2t_begin, depth_begin;
    double *min;


    if (((min = (double *) malloc(sizeof(double) * G->n))) == NULL) {
        fprintf(stderr, "set_all_edge(): ERROR [min = malloc()]\n");
        exit (EXIT_FAILURE);
    }
    // min[i]の全てにINFを設定
    for (i = 0; i < G->n; i++) min[i] = INF;

    // 格子部分1階層分の点数合計
    grids_node = height * width;

    for (i = 1; i < G->n; i++) G->capa[i] = 0;
    source = grids_node * ls[0] + 1;
    sink = source + 1;

    setSource(G, source);
    setSink(G, sink);

    edge_count = 1;
    // source->i1
    s2i_begin = edge_count;data_str(i, ls[ls[0]], width, I_left, I_right)
    for (i = 1; i <= grids_node; i++) {
        setEdge(G, edge_count, source, i, 0);
        if(isin_array(ls, label[i])) {
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
            if(isin_array(ls , label[i])) {
                G->capa[edge_count] = data(I_left[i],ls[j]) + nnp_4_grsa(i, j, height, width, ls, label, T, lamda);
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
        if(isin_array(ls , label[i])) {
            // G->capa[edge_count] = data_str(i, ls[ls[0]], width, I_left, I_right) + nnp_4_grsa(i, ls[0], height, width, ls, label, T, lamda);
            G->capa[edge_count] = data(I_left[i], ls[ls[0]]) + nnp_4_grsa(i, ls[0], height, width, ls, label, T, lamda);
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
            if(isin_array(ls, label[i])) {
                G->capa[edge_count] = INF;
            }
            edge_count++;
            tail = head;
        }
    }

    

    // new codes
    // horizonal
    for (i = 1; i <= height; i++) {
        for (j = 1; j < width; j++) {
            t_base =  (i - 1) * width + j;
            h_base =  t_base + 1;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(isin_array(ls, label[t_base]) && isin_array(ls, label[h_base])) {
                        // head, tail in label_index
                        G->capa[edge_count] = phi(k + 1, l + 1, ls, T);
                    }
                    edge_count++;
                }
            }
        }
    }
    
    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            t_base = (i - 1) * width + j;
            h_base = t_base + width;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, tail, head, 0);
                    if(isin_array(ls, label[t_base]) && isin_array(ls, label[h_base])) {
                        G->capa[edge_count] = phi(k + 1, l + 1, ls, T);
                        // head, tail in label_index
                    }
                    edge_count++;
                }
            }
        }
    }




    // new codes
    // horizonal
    for (i = 1; i <= height; i++) {
        for (j = 1; j < width; j++) {
            t_base =  (i - 1) * width + j;
            h_base =  t_base + 1;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    if(isin_array(ls, label[t_base]) && isin_array(ls, label[h_base])) {
                        G->capa[edge_count] = phi(l + 1, k + 1, ls, T);
                    }
                    edge_count++;
                }
            }
        }
    }

    // vertical
    for (i = 1; i < height ; i++){
        for (j = 1; j < width + 1; j++) {
            t_base = (i - 1) * width + j;
            h_base = t_base + width;
            for (k = 1; k < ls[0]; k++) {
                tail = t_base + k * grids_node;
                for (l = 1; l < ls[0]; l++) {
                    head = h_base + l * grids_node;
                    setEdge(G, edge_count, head, tail, 0);
                    if(isin_array(ls, label[t_base]) && isin_array(ls, label[h_base])) {
                        G->capa[edge_count] = phi(l + 1, k + 1, ls, T);
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
    
    // printf("total edge : %d\n", edge_count - 1);
    return;
}
**/