#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "bmp.h"
#include "grsa.h"
#include "graph.h"
#include "ford_fulkerson.h"

#include <float.h>
#define INF DBL_MAX

/*

*/

#define _OUTPUT_T_ 0     // BK-maxflow後のtの状態をファイルに出力 0:出力しない 1:出力する
#define _OUTPUT_INFO_ 0     // デバッグ情報出力 0:出力しない 1:出力する
#define _OUTPUT_GRAPH_ 0    // グラフ情報出力  0:出力しない 1:出力する
#define _OUTPUT_PROGRESS_ 0 // 処理過程ファイル出力 0:出力しない 1:出力する
#define _RUN_FIRST_ONLY_ 0 // 1度目の移動で終了(デバッグ用)
#define _SHOW_EACH_ENERGY_ 1 // 各移動時にエネルギー表示
#define _OUTPUT_SUBMODULAR_SUBSETS_ 0


int main(int argc, char *argv[]) {
    int i, j, k, l, node, edge, grids_node, flag, size, ccvex, prev;
    int scale, label_max, grids_edge, count, last_move, ci, large_array, total_ss_count;
    int *I_left, *I_right, *t, *label, *newlabel, *label_index;
    int **ls;
    int convex[10][2];
    // I->入力画像の輝度, t->2値変数, label->ラベル付け
    int label_size = 16;
    int range_size = 4;
    int errlog = 0;
    int error_count = 0;
    int lamda = 1;
    double decreace, prev_energy, before_energy, new_energy, err, T = INF;
    double *f;
    char output_file[100];
    clock_t start;
    img left, right, output, truth;
    // Ge:エネルギー計算用
    Graph G, Ge;

    char imgleft[100];
    char imgright[100];
    char imgtruth[100];
    
    // スタックサイズ変更
    system("ulimit -s unlimited");

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

    strcpy(imgleft, argv[1]);
    strcpy(imgright, argv[1]);
    strcpy(imgtruth, argv[1]);

    strcat(imgleft, "left.bmp");
    strcat(imgright, "right.bmp");
    strcat(imgtruth, "truth.bmp");

    T = theta(range_size, INF);
    // T = range_size;

    
    // ccvex 凸区間の数(count of convex)
    ccvex = 0;

    // for (i = 0; i < label_size; i++) {
    //     printf("%lf ", theta(i ,T));
    // }
    // printf("\n");

    i = 0;
    j = 1;
    while (j < label_size) {
        j++;

        // printf("i: %d, j: %d isc: %f %f\n", i, j, theta(j, T) - theta(j - 1, T) , (theta(j - 1, T) - theta(i, T)) / (j - 1 - i));
        if (j > 1  && is_convex(i, j, T)) {
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
        printf("T : %.0lf 候補区間: %d --> %d\n", T, convex[i][0], convex[i][1]);
    }

    // exit(EXIT_SUCCESS);
    // generate submodular subsets
    // ls[i][0] == 劣モジュラ部分集合iの要素数
    // ls[i][1] ~ ls[i][range_size] 劣モジュラ部分集合

    // total_ss_count = gen_submodular_subsets(label_size, range_size, ls);
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

        // printf("size : %d large_array %d\n", size, large_array);
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

                for (j = 2; j <= range_size; j++) {
                    ls[i][j] = ls[i][j - 1] + 1;
                    k = ls[i][j];
                }
            }
        
            if ((ls[large_array] = (int *)malloc(sizeof(int) * size)) == NULL) {
                fprintf(stderr, "Error!:malloc[main()->ls]\n");
                exit(EXIT_FAILURE);
            }
            
            ls[large_array][0] = size;
            ls[large_array][1] = k;

            for (j = 2; ls[large_array][j - 1] + 1 <= label_max; j++) {
                    ls[large_array][j] = ls[large_array][j - 1] + 1;
            }

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

#if _OUTPUT_SUBMODULAR_SUBSETS_
    for (i = 1; i <= total_ss_count; i++) {
        for (j = 1; j <= ls[i][0]; j++) {
            printf("%d ", ls[i][j]);
        }
        printf("\n");
    }
#endif

   
    printf("----------------------------------------------\n");
    printf("input_file: %s\n", argv[1]);
    printf("output_file: %s\n", output_file);
    printf("label_size: %d\n", label_size);
    printf("range_size: %d\n", range_size);
    printf("lambda: %d\n", lamda);
    printf("T: %.2f\n", T);
    if(theta(2, 2 * 2) > 2) printf("Vpq(fp, fq) = (fp - fq)^2\n");
    else printf("Vpq(fp, fq) = |fp - fq|\n");

    ReadBmp(imgleft, &left);
    ReadBmp(imgright, &right);
    ReadBmp(imgtruth, &truth);
    ReadBmp(imgtruth, &output);

    Gray(&left, &left);
    Gray(&right, &right);
    Gray(&truth, &truth);

    grids_node = left.height * left.width;


    if ((I_left = (int *)malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main()->I_left]\n");
        exit(EXIT_FAILURE);
    }
    if ((I_right = (int *)malloc(sizeof(int) * (grids_node + 1))) == NULL) {
        fprintf(stderr, "Error!:malloc[main()->I_right]\n");
        exit(EXIT_FAILURE);
    }



    printf("height %ld, width %ld\n", left.height, left.width);
    for (i = 0; i <  left.height; i++) {
        for (j = 0; j < left.width; j++) {
            I_left[i * left.width + j + 1] = left.data[i][j].r / scale;
        }
    }
    for (i = 0; i <  right.height; i++) {
        for (j = 0; j < right.width; j++) {
            I_right[i * right.width + j + 1] = right.data[i][j].r / scale;
        }
    }

    // エネルギー計算用一層グラフ作成
    node = grids_node + 2;
    edge = (left.height - 1) * left.width + left.height * (left.width - 1) + 2 * grids_node;
    newGraph(&Ge, node, edge);

    // set_single_edge(&Ge, left.height, left.width);
    set_single_edges(&Ge, left.height, left.width);
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
    prev_energy = energy_str(&Ge, label, T, lamda, left.width, I_left, I_right);
    printf("Energy (before): %.0lf\n", prev_energy);


#if _OUTPUT_T_
    fprintf(fp, "Energy (before): %lf\n", prev_energy);
    fprintf(fp, "position :\n");
    for (i = 1; i <= Ge.n - 2; i++) {
        // printf("t[%d] : %d\n", i, t[i]);
        fprintf(fp, "%d ", i);
        if(i % left.width == 0) fprintf(fp, "\n");
        if(i % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
    }
    fprintf(fp, "init_label:\n");
    for (i = 1; i <= Ge.n - 2; i++) {
        // printf("t[%d] : %d\n", i, t[i]);
        fprintf(fp, "%d ", label[i]);
        if(i % left.width == 0) fprintf(fp, "\n");
        if(i % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
    }

#endif

    last_move = total_ss_count + 1;
    decreace = 0;
    flag = 0;
    ci = 0;
    start = clock();

    do {
        prev_energy = energy_str(&Ge, label, T, lamda, left.width, I_left, I_right);
        for(i = 1; i <= total_ss_count; i++) {
            if (last_move == i) {
                flag = 1;
                break;
            }
            before_energy = energy_str(&Ge, label, T, lamda, left.width, I_left, I_right);
            
#if _OUTPUT_T_
            fprintf(fp, "\n-------------------------------------\n");
            fprintf(fp, "submodular subsets: ");
            for (j = 1; j <= ls[i][0]; j++) {
                fprintf(fp, "%d ", ls[i][j]);
            }
            fprintf(fp, "\n");
            
            for (j = 1; j <= Ge.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", isin_array(ls[i], label[j]) ? 1 : 0);
                if(j % left.width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }
            
            fprintf(fp, "label: \n");
            for (j = 1; j <= Ge.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", label[j]);
                if(j % left.width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }      
#endif
         
            
            node = left.height * left.width * ls[i][0] + 2;
            grids_edge = (left.height - 1) * left.width + left.height * (left.width - 1);
            edge = 2 * grids_node * ls[i][0] + 2 * grids_edge * (ls[i][0] - 1) * ((ls[i][0] - 1)); 

            newGraph(&G, node, edge);
            set_edge_str(&G, left.height, left.width, ls[i], label, T, lamda, I_left, I_right);
            // set_edge(&G, left.height, left.width, ls[i], label, I_left, T, lamda);
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
                if (isin_array(ls[i], label[j])) {
                    k = j;
                    count = 0;
                    while (k <= ls[i][0] * grids_node && t[k] == 1) {
                        // if (k + grids_node > G.n) break;
                        k += grids_node;
                        count++;
                    }
                    newlabel[j] = ls[i][count];
                } else newlabel[j] = label[j];
            }

            new_energy = energy_str(&Ge, newlabel, T, lamda, left.width, I_left, I_right);

            if (new_energy <= before_energy) {
                last_move = i;
                cpyarray(label, newlabel, grids_node);
            } else if (new_energy > before_energy) {
                errlog = 1;
                
                printf("err %lf -> %lf\n", energy_str(&Ge,  label, T, lamda, left.width, I_left, I_right), energy_str(&Ge, newlabel, T, lamda, left.width, I_left, I_right));
                for (j = 1; j <= ls[i][0]; j++) printf("%d ", ls[i][j]);
                printf("\n");
            }

#if _OUTPUT_T_
            fprintf(fp, "t: \n");
            for (j = 1; j <= G.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", t[j]);
                if(j % left.width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }
            fprintf(fp, "label: \n");
            for (j = 1; j <= Ge.n - 2; j++) {
                // printf("t[%d] : %d\n", i, t[i]);
                fprintf(fp, "%d ", newlabel[j]);
                if(j % left.width == 0) fprintf(fp, "\n");
                if(j % (grids_node) == 0) fprintf(fp, "-------------------------------------\n");
            }
#endif

#if _OUTPUT_PROGRESS_
            for (j = 0; j <  left.height; j++) {
                for (k = 0; k < left.width; k++) {
                    output.data[j][k].r = label[j * left.width + k + 1] * scale;
                    output.data[j][k].g = output.data[j][k].r;
                    output.data[j][k].b = output.data[j][k].r;
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
        decreace = prev_energy - energy_str(&Ge, label, T, lamda, left.width, I_left, I_right);
#if _SHOW_EACH_ENERGY_
        printf("Energy : %.0lf\n", energy_str(&Ge, label, T, lamda, left.width, I_left, I_right));
#endif
    } while (decreace > 0);
    
#if _OUTPUT_T_
    fprintf(fp, "result:\n");
    for (i = 1; i <= Ge.n - 2; i++) {
        fprintf(fp, "%d ", label[i]);
        if(i % left.width == 0) fprintf(fp, "\n");
    }
#endif

    printf("Energy (after): %.0lf\n", energy_str(&Ge, label, T, lamda, left.width, I_left, I_right));
    printf("Italation: %d\n", ci);
    printf("Run time[%.2lf]\n", (double) (clock() - start) / CLOCKS_PER_SEC);
    
    
    // output to bitmap file
    for (i = 0; i <  left.height; i++) {
        for (j = 0; j < left.width; j++) {
            output.data[i][j].r = label[i * left.width + j + 1] * scale;
            output.data[i][j].g = output.data[i][j].r;
            output.data[i][j].b = output.data[i][j].r;
        }
    }
    WriteBmp(output_file, &output);
    
#if _OUTPUT_T_
    fprintf(fp, "Energy (after): %lf\n", energy_str(&Ge, label, T, lamda, left.width, I_left, I_right));
    fclose(fp);
#endif

    if(errlog) printf("エネルギーが増大する移動が確認されました\n");
    
    if (strcmp(output_file, "/dev/null") != 0){
        ReadBmp(output_file, &output);
        // Gray(&truth, &truth);
        Gray(&output, &output);

        if (truth.data[0][0].r) {
            for(i = 1; i <= (output.height) * (output.width); i++) {   
                if (abs(output.data[(i - 1) / output.width][(i - 1) % output.width].r - truth.data[(i - 1) / truth.width][(i - 1) % truth.width].r ) 
                    >= scale + 1) {
                    error_count++;
                }                            
            }
        } else {
            for(i = 1; i <= (output.height) * (output.width); i++) {
                if ((i - 1) / output.width >= scale && (i - 1) % output.width >= scale &&
                    (i - 1) / output.width <= output.height - scale && (i - 1) % output.width <= output.width - scale) {
                    if (abs(output.data[(i - 1) / output.width][(i - 1) % output.width].r - truth.data[(i - 1) / truth.width][(i - 1) % truth.width].r ) 
                        >= scale + 1) {
                        error_count++;
                    }
                }                            
            }
        }

        err = 100 * error_count / (double)(truth.height * truth.width);
        printf("Error rate : %lf\n", err);
    }

    // free meory
    delGraph(&Ge);
    for (i = 0; i <= total_ss_count; i++) {
        free(ls[i]);
    }
    free(ls);

    free(I_left);
    free(I_right);
    free(label);
    free(newlabel);
    free(label_index);
    printf("----------------------------------------------\n");
    return 0;
    
}
