#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
// #include <math.h>
#include "bmp.h"
#include "grsa.h"
#include "graph.h"
#include "ford_fulkerson.h"
#include "malloc.h"
#include <float.h>
#define INF DBL_MAX

/*

*/

#define _OUTPUT_T_ 0     // BK-maxflow後のtの状態をファイルに出力 0:出力しない 1:出力する
#define _OUTPUT_INFO_ 0     // デバッグ情報出力 0:出力しない 1:出力する
#define _OUTPUT_GRAPH_ 0    // グラフ情報出力  0:出力しない 1:出力する
#define _OUTPUT_PROGRESS_ 0 // 処理過程ファイル出力 0:出力しない 1:出力する
#define _RUN_FIRST_ONLY_ 0 // 1度目の移動で終了(デバッグ用)
#define _SHOW_EACH_ENERGY_ 0 // 各移動時にエネルギー表示
#define _OUTPUT_SUBMODULAR_SUBSETS_ 0


void GenAllPairs(int **pairs, int label_size) {
    int i, j, k;
    int label_max = label_size - 1;
    if ((pairs = (int **)malloc(sizeof(int*) * (nc2(label_size) + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main()->pairs]\n");
        exit(EXIT_FAILURE);
    }
    printf("%d\n", nc2(label_size));

    k = 1;
    for (i = 0; i < label_max; i++) {
        for (j = i + 1; j <= label_max; j++) {
            if ((pairs[k] = (int*)malloc(sizeof(int) * 3)) == NULL) {
                fprintf(stderr, "Error!:malloc[main()->pairs]\n");
                exit(EXIT_FAILURE);
            }
            pairs[k][0] = 0;
            pairs[k][1] = i;
            pairs[k][2] = j;

            printf("(%d, %d): %d\n", i, j, k);
            k++;
        }
    }
}

int main(int argc, char *argv[]) {
    int i, j, k, l, m, n, current_array, node, edge, grids_node, flag, size, ccvex, prev, rs2;
    int scale, label_max, grids_edge, count, last_move, ci, large_array, total_ss_count;
    int *t, *label, *newlabel, *label_index;
    int convex[10][2];
    // I->入力画像の輝度, t->2値変数, label->ラベル付け
    int label_size = 16;
    int range_size = 4;
    int errlog = 0;
    int error_count = 0;
    int lamda = 1;
    double decreace, prev_energy, before_energy, new_energy, err,  T = INF;
    double *f;
    char output_file[100];
    clock_t start;
    img output;
    // Ge:エネルギー計算用
    Graph G, Ge;
    Input image;
    Subsets ss;


#if _OUTPUT_INFO_
    double maxflow;
#endif

#if _OUTPUT_T_
    FILE *fp;
    fp = fopen("log/t.txt", "w");
    if (fp == NULL) {
        fprintf(stderr, "cannot open file[t.txt]\n");
        exit (EXIT_FAILURE);
    }
#endif
#if _OUTPUT_PROGRESS_
    int l = 0;
    char pf[100];
    system("rm output/*.bmp &> /dev/null");
#endif
    function = 1;
    if (argc != 2 && argc != 3 && argc != 5 && argc != 6) {
        printf("Usage: %s <input_file> <output_file(option)> <range_size(option)> <scale(option)> <lamda (option)>\n", argv[0]);
        return 1;
    }

    if (argc == 2) strcpy(output_file, "/dev/null");
    else strcpy(output_file, argv[2]);
    if (argc == 5 || argc == 6) {
        range_size = atoi(argv[3]);
        scale = atoi(argv[4]);
        lamda = atoi(argv[5]);
        if(argc == 7) T = atof(argv[6]);
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
    
    T = theta(range_size, INF);
    range_size = 5;
    // ccvex 凸区間の数(count of convex)
    gen_submodular_subsets(label_size, range_size, &ss);
    readStrBitmap(&image, argv[1], scale);

    printf("----------------------------------------------\n");
    printf("input_file: %s\n", argv[1]);
    printf("output_file: %s\n", output_file);
    printf("label_size: %d\n", label_size);
    printf("range_size: %d\n", range_size);
    printf("lambda: %d\n", lamda);
    printf("T: %.2f\n", T);
    if(theta(2, 2 * 2) > 2) printf("Vpq(fp, fq) = (fp - fq)^2\n");
    else printf("Vpq(fp, fq) = |fp - fq|\n");

    grids_node = image.height * image.width;

    // エネルギー計算用一層グラフ作成
    node = grids_node + 2;
    edge = (image.height - 1) * image.width + image.height * (image.width - 1) + 2 * grids_node;
    newGraph(&Ge, node, edge);

    // set_single_edge(&Ge, image.height, image.width);
    set_single_edges(&Ge, image.height, image.width);
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
    for (i = 1; i <= grids_node ; i++) label[i] = 0;
    cpyarray(newlabel, label, grids_node);
    prev_energy = energy_str(&Ge, label, T, lamda, image.width, image.left, image.right);
    printf("Energy (before): %.0lf\n", prev_energy);


#if _OUTPUT_T_
    fprintf(fp, "Energy (before): %lf\n", prev_energy);
    fprintf(fp, "position :\n");
    for (i = 1; i <= Ge.n - 2; i++) {
        // printf("t[%d] : %d\n", i, t[i]);
        fprintf(fp, "%d ", i);
        if(i % image.width == 0) fprintf(fp, "\n");
        if(i % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
    }
    fprintf(fp, "init_label:\n");
    for (i = 1; i <= Ge.n - 2; i++) {
        // printf("t[%d] : %d\n", i, t[i]);
        fprintf(fp, "%d ", label[i]);
        if(i % image.width == 0) fprintf(fp, "\n");
        if(i % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
    }

#endif

    last_move = total_ss_count + 1;
    decreace = 0;
    flag = 0;
    ci = 0;
    start = clock();

    do {
        prev_energy = energy_str(&Ge, label, T, lamda, image.width, image.left, image.right);
        for(i = 1; i <= total_ss_count; i++) {
            if (last_move == i) {
                flag = 1;
                break;
            }
            before_energy = energy_str(&Ge, label, T, lamda, image.width, image.left, image.right);

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
                if(j % image.width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }

            fprintf(fp, "label: \n");
            for (j = 1; j <= Ge.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", label[j]);
                if(j % image.width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }
#endif


            node = image.height * image.width * ss.ls[i][0] + 2;
            grids_edge = (image.height - 1) * image.width + image.height * (image.width - 1);
            edge = 2 * grids_node * ss.ls[i][0] + 2 * grids_edge * (ss.ls[i][0] - 1) * ((ss.ls[i][0] - 1));

            newGraph(&G, node, edge);
            set_edge_str(&G, image.height, image.width, ss.ls[i], label, T, lamda, image.left, image.right);
            // set_edge(&G, image.height, image.width, ss.ls[i], label, image.left, T, lamda);
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

            new_energy = energy_str(&Ge, newlabel, T, lamda, image.width, image.left, image.right);

            if (new_energy <= before_energy) {
                last_move = i;
                cpyarray(label, newlabel, grids_node);
            } else if (new_energy > before_energy) {
                errlog = 1;

                printf("err %lf -> %lf\n", energy_str(&Ge,  label, T, lamda, image.width, image.left, image.right), energy_str(&Ge, newlabel, T, lamda, image.width, image.left, image.right));
                for (j = 1; j <= ss.ls[i][0]; j++) printf("%d ", ss.ls[i][j]);
                printf("\n");
            }

#if _OUTPUT_T_
            fprintf(fp, "t: \n");
            for (j = 1; j <= G.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", t[j]);
                if(j % image.width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }
            fprintf(fp, "label: \n");
            for (j = 1; j <= Ge.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", newlabel[j]);
                if(j % image.width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }
#endif

#if _OUTPUT_PROGRESS_
            for (j = 0; j <  image.height; j++) {
                for (k = 0; k < image.width; k++) {
                    image.output.data[j][k].r = label[j * image.width + k + 1] * scale;
                    image.output.data[j][k].g = image.output.data[j][k].r;
                    image.output.data[j][k].b = image.output.data[j][k].r;
                }
            }
            sprintf(pf, "output/left_%04d.bmp", l);
            WriteBmp(pf, &image.output);
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
        decreace = prev_energy - energy_str(&Ge, label, T, lamda, image.width, image.left, image.right);
#if _SHOW_EACH_ENERGY_
        printf("Energy : %.0lf\n", energy_str(&Ge, label, T, lamda, image.width, image.left, image.right));
#endif
    } while (decreace > 0);

#if _OUTPUT_T_
    fprintf(fp, "result:\n");
    for (i = 1; i <= Ge.n - 2; i++) {
        fprintf(fp, "%d ", label[i]);
        if(i % image.width == 0) fprintf(fp, "\n");
    }
#endif

    printf("Energy (after): %.0lf\n", energy_str(&Ge, label, T, lamda, image.width, image.left, image.right));
    printf("Italation: %d\n", ci);
    printf("Run time[%.2lf]\n", (double) (clock() - start) / CLOCKS_PER_SEC);


    // output to bitmap file
    for (i = 0; i <  image.height; i++) {
        for (j = 0; j < image.width; j++) {
            image.output.data[i][j].r = label[i * image.width + j + 1] * scale;
            image.output.data[i][j].g = image.output.data[i][j].r;
            image.output.data[i][j].b = image.output.data[i][j].r;
        }
    }
    WriteBmp(output_file, &image.output);

#if _OUTPUT_T_
    fprintf(fp, "Energy (after): %lf\n", energy_str(&Ge, label, T, lamda, image.width, image.left, image.right));
    fclose(fp);
#endif

    if(errlog) printf("エネルギーが増大する移動が確認されました\n");

    if (strcmp(output_file, "/dev/null") != 0){
        ReadBmp(output_file, &(image.output));
        // Gray(&truth, &truth);
        Gray(&(image.output), &(image.output));

        if (image.truth.data[0][0].r) {
            for(i = 1; i <= (image.output.height) * (image.output.width); i++) {
                if (abs(image.output.data[(i - 1) / image.output.width][(i - 1) % image.output.width].r - image.truth.data[(i - 1) / image.truth.width][(i - 1) % truth.width].r )
                    >= scale + 1) {
                    error_count++;
                }
            }
        } else {
            for(i = 1; i <= (image.output.height) * (image.output.width); i++) {
                if ((i - 1) / image.output.width >= scale && (i - 1) % image.output.width >= scale &&
                    (i - 1) / image.output.width <= image.output.height - scale && (i - 1) % image.output.width <= image.output.width - scale) {
                    if (abs(image.output.data[(i - 1) / image.output.width][(i - 1) % image.output.width].r - image.truth.data[(i - 1) / image.truth.width][(i - 1) % truth.width].r )
                        >= scale + 1) {
                        error_count++;
                    }
                }
            }
        }

        err = 100 * error_count / (double)(image.truth.height * image.truth.width);
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

    free(image.left);

    free(image.right);

    free(label);

    free(newlabel);

    printf("----------------------------------------------\n");
    return 0;

}
