#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <float.h>
// #include <math.h>
#include "bmp.h"
#include "grsa.h"
#include "graph.h"
#include "ford_fulkerson.h"



#define INF DBL_MAX

#define _OUTPUT_T_ 0     // BK-maxflow後のtの状態をファイルに出力 0:出力しない 1:出力する
#define _OUTPUT_INFO_ 0     // デバッグ情報出力 0:出力しない 1:出力する
#define _OUTPUT_GRAPH_ 0    // グラフ情報出力  0:出力しない 1:出力する
#define _OUTPUT_PROGRESS_ 0 // 処理過程ファイル出力 0:出力しない 1:出力する
#define _RUN_FIRST_ONLY_ 0 // 1度目の移動で終了(デバッグ用)
#define _SHOW_EACH_ENERGY_ 0 // 各移動時にエネルギー表示
#define _OUTPUT_SUBMODULAR_SUBSETS_ 1

int main(int argc, char *argv[]) {
    int i, j, k, node, edge, grids_node, flag, height, width;
    int scale, grids_edge, count, last_move, ci, total_ss_count, label_max;
    int *t, *label, *newlabel, *label_index, *left, *right;
    int label_size = 16;
    int range_size = 4;
    int errlog = 0;
    int lambda = 1;
    double decreace, prev_energy, before_energy, new_energy, err,  T = INF;
    double *f;
    char output_file[100];
    char imgleft[100];
    char imgright[100];
    char imgtruth[100];
    clock_t start;
    // Ge:エネルギー計算用
    Graph G, Ge;
    Subsets ss;
    img *raw_left, *raw_right, *truth, *output;
    
#if _OUTPUT_PROGRESS_
    int l = 0;
    char pf[100];
    system("rm output/*.bmp &> /dev/null");
#endif

    dterm = 0;
    function = 1;
    if (argc != 2 && argc != 3 && argc != 5 && argc != 6 && argc != 7 && argc != 8) {
        printf("Usage: %s <input_file> <output_file(option)> <range_size(option)> <scale(option)> <lambda (option)>  <Dterm 0: Dt, 1:normal(option)> <T (option)>\n", argv[0]);
        return (EXIT_FAILURE);
    }

    // Read Bitmap
    raw_left = (img *)malloc(sizeof(img));
    raw_right = (img *)malloc(sizeof(img));
    output = (img *)malloc(sizeof(img));
    truth = (img *)malloc(sizeof(img));

    strcpy(imgleft, argv[1]);
    strcpy(imgright, argv[1]);
    strcpy(imgtruth, argv[1]);

    strcat(imgleft, "left.bmp");
    strcat(imgright, "right.bmp");
    strcat(imgtruth, "truth.bmp");
    
    ReadBmp(imgleft, raw_left);
    ReadBmp(imgright, raw_right);
    ReadBmp(imgtruth, truth);
    ReadBmp(imgtruth, output);

    if(raw_left->width != raw_right->width || raw_left->height != raw_right->height) {
        fprintf(stderr, "Error %s と %s の解像度が異なります\n", imgleft, imgright);
        exit(EXIT_FAILURE);
    }

    Gray(raw_left, raw_left);
    Gray(raw_right, raw_right);
    Gray(truth, truth);

    height = raw_left->height;
    width = raw_left->width;

    if (argc == 2) strcpy(output_file, "/dev/null");
    else strcpy(output_file, argv[2]);
    if (argc >= 5 && argc <= 8) {
        range_size = atoi(argv[3]);
        scale = atoi(argv[4]);
        lambda = atoi(argv[5]);
        if(argc >= 7)  dterm = atoi(argv[6]);
        if(argc == 8)  T = atof(argv[7]);
    }

    if(T < range_size) {
        fprintf(stderr, "error! T (%f) < range_size (%d)\n", T, range_size);
        exit (EXIT_FAILURE);
    }

    label_size = 256 / scale;
    label_max = label_size - 1;
    if (range_size < 2) {
        fprintf(stderr, "Error! Range size == %d \n", range_size);
        exit (EXIT_FAILURE);
    }
    if (label_size < range_size) {
        fprintf(stderr, "Error! label_size < range_size \n");
        exit (EXIT_FAILURE);
    }

    grids_node = height * width;

    if ((left = (int *)malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main() int *left]\n");
        exit(EXIT_FAILURE);
    }
    if ((right = (int *)malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main() int *right]\n");
        exit(EXIT_FAILURE);
    }

    for (i = 0; i <  height; i++) {
        for (j = 0; j < width; j++) {
            left[i * width + j + 1] = raw_left->data[i][j].r / scale;
        }
    }
    for (i = 0; i <  height; i++) {
        for (j = 0; j < width; j++) {
            right[i * width + j + 1] = raw_right->data[i][j].r / scale;
        }
    }

    // printf("T: %f\n",theta(T, INF));
    if(T == INF) T = theta(range_size, INF);
    else{
        i = theta(T, INF);;
        T = i;
    } 
    ss.T = T;
    // ccvex 凸区間の数(count of convex)
    total_ss_count = gen_submodular_subsets(label_size, range_size, &ss);
    
    printf("submodular subsets: \n");
    for (i = 1; i <= total_ss_count; i++) {
        if (ss.ls[i][0] != 1) {
            printf("%d, (%d) ",i, ss.ls[i][0]);
            for (j = 1; j <= ss.ls[i][0]; j++) {
                printf("%d ", ss.ls[i][j]);
            }
            printf(" end\n");
        }
    }

    printf("----------------------------------------------\n");
    printf("input_file: %s\n", argv[1]);
    printf("output_file: %s\n", output_file);
    printf("label_size: %d\n", label_size);
    printf("range_size: %d\n", range_size);
    printf("lambda: %d\n", lambda);
    printf("T: %.2f\n", T);
    printf("Data term: ");
    if(dterm == 0) printf("Dt\n");
    else printf("normal\n");
    if(theta(2, 2 * 2) > 2) printf("Vpq(fp, fq) = (fp - fq)^2\n");
    else printf("Vpq(fp, fq) = |fp - fq|\n");
    // エネルギー計算用一層グラフ作成
    node = grids_node + 2;
    edge = (height - 1) * width + height * (width - 1) + 2 * grids_node;
    newGraph(&Ge, node, edge);

    set_single_edges(&Ge, height, width);
    initAdjList(&Ge);

    if ((label = (int *) malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "main(): ERROR [label = malloc()]\n");
        return (EXIT_FAILURE);
    }
    if ((newlabel = (int *) malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "main(): ERROR [label = malloc()]\n");
        return (EXIT_FAILURE);
    }
    if ((label_index = (int *) malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "main(): ERROR [label_index = malloc()]\n");
        return (EXIT_FAILURE);
    }

    // 輝度から初期ラベル設定
    for (i = 1; i <= grids_node ; i++) label[i] = 5;
    cpyarray(newlabel, label, grids_node);
    prev_energy = energy_str(&Ge, label,  T, lambda, height, width, label_max, left, right, raw_left, raw_right);
    printf("Energy (before): %.0lf\n", prev_energy);

    #if _OUTPUT_T_
    fprintf(fp, "Energy (before): %lf\n", prev_energy);
    fprintf(fp, "position :\n");
    for (i = 1; i <= Ge.n - 2; i++) {
        // printf("t[%d] : %d\n", i, t[i]);
        fprintf(fp, "%d ", i);
        if(i % raw_left->width == 0) fprintf(fp, "\n");
        if(i % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
    }
    fprintf(fp, "init_label:\n");
    for (i = 1; i <= Ge.n - 2; i++) {
        // printf("t[%d] : %d\n", i, t[i]);
        fprintf(fp, "%d ", label[i]);
        if(i % raw_left->width == 0) fprintf(fp, "\n");
        if(i % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
    }

    #endif

    last_move = total_ss_count + 1;
    decreace = 0;
    flag = 0;
    ci = 0;
    start = clock();

    do {
        prev_energy = energy_str(&Ge, label,  T, lambda, height, width, label_max, left, right, raw_left, raw_right);
        for(i = 1; i <= total_ss_count; i++) {
            if (last_move == i) {
                flag = 1;
                break;
            }
            before_energy = energy_str(&Ge, label,  T, lambda, height, width, label_max, left, right, raw_left, raw_right);

            #if _OUTPUT_T_
            fprintf(fp, "\n-------------------------------------\n");
            fprintf(fp, "submodular subsets: ");
            for (j = 1; j <= ss.ls[i][0]; j++) {
                fprintf(fp, "%d ", ss.ls[i][j]);
            }
            fprintf(fp, "\n");

            for (j = 1; j <= Ge.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", isin_array(ss.ls[i], label[j]) ? 1 : 0);
                if(j % raw_left->width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }

            fprintf(fp, "label: \n");
            for (j = 1; j <= Ge.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", label[j]);
                if(j % raw_left->width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }
            #endif

            if (ss.ls[i][0] == 1) continue;
            // printf("submodular subsets: ");
            // printf("%d, (%d) ",i, ss.ls[i][0]);
            // for (j = 1; j <= ss.ls[i][0]; j++) {
            //     printf("%d ", ss.ls[i][j]);
            // }
            // printf(" end\n");
            node = height * width * ss.ls[i][0] + 2;
            grids_edge = (height - 1) * width + height * (width - 1);
            edge = 2 * grids_node * ss.ls[i][0] + 2 * grids_edge * (ss.ls[i][0] - 1) * ((ss.ls[i][0] - 1));

            newGraph(&G, node, edge);
            set_edge(&G, ss.ls[i], label, T, lambda, label_max, left, right, raw_left, raw_right);
            // resizeGraph(&G, node, edge);
            // set_edge(&G, image.height, image.width, ss.ls[i], label, image.left, T, lambda);
            initAdjList(&G);

            if ((f = (double *) malloc(sizeof(double) * (G.m + 1))) == NULL) {
                fprintf(stderr, "main(): ERROR [f = malloc()]\n");
                return (EXIT_FAILURE);
            }
            if ((t = (int *) malloc(sizeof(int) * (G.n + 1))) == NULL) {
                fprintf(stderr, "main(): ERROR [t = malloc()]\n");
                return (EXIT_FAILURE);
            }

            for (j = 0; j < G.m + 1 ; j++) f[j] = 0;
            for (j = 0; j < G.n + 1 ; j++) t[j] = 0;
            boykov_kolmogorov(G, f, t);
            ci++;

            for (j = 1; j <= Ge.n - 2; j++) {
                if (isin_array(ss.ls[i], label[j])) {
                    k = j;
                    count = 0;
                    while (k <= ss.ls[i][0] * grids_node && t[k] == 1) {
                        // if (k + grids_node > G.n) break;
                        k += grids_node;
                        count++;
                    }
                    newlabel[j] = ss.ls[i][count];
                } else newlabel[j] = label[j];
            }
            new_energy = energy_str(&Ge, newlabel,  T, lambda, height, width, label_max, left, right, raw_left, raw_right);
            // printf("Energy : %.0lf\n", new_energy);
            if (new_energy < before_energy) {
                last_move = i;
                cpyarray(label, newlabel, grids_node);
            } else if (new_energy > before_energy) {
                errlog = 1;
                printf("err %lf -> %lf\n", before_energy, new_energy);
                for (j = 1; j <= ss.ls[i][0]; j++) printf("%d ", ss.ls[i][j]);
                printf("\n");
                exit(EXIT_FAILURE);
            }

            #if _OUTPUT_T_
            fprintf(fp, "t: \n");
            for (j = 1; j <= G.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", t[j]);
                if(j % raw_left->width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }
            fprintf(fp, "label: \n");
            for (j = 1; j <= Ge.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", newlabel[j]);
                if(j % raw_left->width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }
            #endif

            #if _OUTPUT_PROGRESS_
            for (j = 0; j <  raw_left->height; j++) {
                for (k = 0; k < raw_left->width; k++) {
                    output->data[j][k].r = label[j * raw_left->width + k + 1] * scale;
                    output->data[j][k].g = output->data[j][k].r;
                    output->data[j][k].b = output->data[j][k].r;
                }
            }
            sprintf(pf, "output/left_%04d.bmp", l);
            WriteBmp(pf, &output);
            l++;
            #endif
            // showGraph(&G);
            free(f);
            free(t);
            delGraph(&G);

            #if _RUN_FIRST_ONLY_
            flag = 1;
            break;
            #endif
        }
        if (flag) break;
        decreace = prev_energy - energy_str(&Ge, label, T, lambda, height, width, label_max, left, right, raw_left, raw_right);

        #if _SHOW_EACH_ENERGY_
        printf("Energy : %.0lf\n", energy_str(&Ge, label, T, lambda, height, width, label_max, left, right, raw_left, raw_right));
        #endif
    } while (decreace > 0);

    if ((f = (double *) malloc(sizeof(double) * (G.m + 1))) == NULL) {
        fprintf(stderr, "main(): ERROR [f = malloc()]\n");
        return (EXIT_FAILURE);
    }
    if ((t = (int *) malloc(sizeof(int) * (G.n + 1))) == NULL) {
        fprintf(stderr, "main(): ERROR [t = malloc()]\n");
        return (EXIT_FAILURE);
    }

    int alpha, beta;
    node = grids_node + 2;
    edge = (height - 1) * width + height * (width - 1) + 2 * grids_node;
    
    newGraph(&G, node, edge);
    set_abswap_edge(&G, height, width);
    initAdjList(&G);

    for (alpha = 0; alpha <= label_max; alpha++) {
        for (beta = 0; beta <= label_max; beta++) {
            for (i = 0; i <= G.m; i++) {
                f[i] = 0;
                G.capa[i] = 0;
            }

            // capacity設定
            set_capacity_abswap(&G, label, alpha, beta, T, lambda, raw_left, raw_right, label_max);
            // capacity(&G, label, I, alpha);
            ci++;
            boykov_kolmogorov(G, f, t);

            // tを基にラベル更新

            for (i = 1; i <= G.n - 2; i++) {
                if (label[i] == alpha || label[i] == beta) {
                    if(t[i] == 1) {
                        newlabel[i] = beta;
                    } else {
                        newlabel[i] = alpha;
                    }
                } else newlabel[i] = label[i];
            }

            if (energy_str(&Ge, newlabel, T, lambda, height, width, label_max, left, right, raw_left, raw_right) < energy_str(&Ge, label, T, lambda, height, width, label_max, left, right, raw_left, raw_right)) {
                cpyarray(label, newlabel, grids_node);
                // last = k;
            } else if (energy_str(&Ge, newlabel, T, lambda, height, width, label_max, left, right, raw_left, raw_right) > energy_str(&Ge, label, T, lambda, height, width, label_max, left, right, raw_left, raw_right)) {
                printf("エネルギー増加\n");
                exit(EXIT_FAILURE);
            }
            #if _SHOW_EACH_ENERGY_
            printf("Energy : %.0lf\n", energy_str(&Ge, label, T, lambda, height, width, label_max, left, right, raw_left, raw_right));
            #endif
        }
    }

    #if _OUTPUT_T_
    fprintf(fp, "result:\n");
    for (i = 1; i <= Ge.n - 2; i++) {
        fprintf(fp, "%d ", label[i]);
        if(i % raw_left->width == 0) fprintf(fp, "\n");
    }
    #endif

    printf("Energy (after): %.0lf\n", energy_str(&Ge, label,  T, lambda, height, width, label_max, left, right, raw_left, raw_right));
    printf("Iteration: %d\n", ci);
    printf("Run time[%.2lf]\n", (double) (clock() - start) / CLOCKS_PER_SEC);


    // output to bitmap file
    for (i = 0; i <  raw_left->height; i++) {
        for (j = 0; j < raw_left->width; j++) {
            output->data[i][j].r = label[i * raw_left->width + j + 1] * scale;
            output->data[i][j].g = output->data[i][j].r;
            output->data[i][j].b = output->data[i][j].r;
        }
    }
    WriteBmp(output_file, output);

    #if _OUTPUT_T_
    fprintf(fp, "Energy (after): %lf\n", energy_str(&Ge, label,  T, lambda, height, width, label_max, left, right, raw_left, raw_right));
    fclose(fp);
    #endif

    if(errlog) printf("エネルギーが増大する移動が確認されました\n");

    if (strcmp(output_file, "/dev/null") != 0){
        ReadBmp(output_file, (output));
        // Gray(&truth, &truth);
        Gray(output, output);
        err = err_rate(output, truth, scale);
        printf("Error rate : %lf\n", err);
    }

    // free meory
    delGraph(&Ge);
    for (i = 0; i <= total_ss_count; i++) {
        free(ss.ls[i]);
    }

    free(ss.ls);

    if (ss.pairs != NULL) {
        for (i = 0; i <= nc2(label_size); i++) {
            free(ss.pairs[i]);
        }
        free(ss.pairs);
    }

    free(left);
    free(right);
    free(label);
    free(newlabel);
    free(raw_left);
    free(raw_right);
    free(output);
    free(truth);

    printf("----------------------------------------------\n");
    return(EXIT_SUCCESS);
}



