#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RKF45.h"

#define SIGMA 10.0
#define RHO 28.0
#define BETA (8.0 / 3.0)

void lorenz_system(double t, const double *y, double *dydt) {
    dydt[0] = SIGMA * (y[1] - y[0]);
    dydt[1] = y[0] * (RHO - y[2]) - y[1];
    dydt[2] = y[0] * y[1] - BETA * y[2];
}

int main() {
    const double t_0 = 0.0;
    const double t_end = 50.0;
    const double y_0[] = {1.0, 1.0, 1.0};
    const size_t n = 3;
    const double tolerance = 1e-5;
    const double initial_step = 0.1;
    const double min_step = 1.0e-3;

    // 結果格納用
    double **result = malloc(100000 * sizeof(double*)); // 十分大きな配列を確保
    if (result == NULL) {
        fprintf(stderr, "Failed to allocate memory for result array.\n");
        return EXIT_FAILURE;
    }

    // 数値計算
    int status = RKF45(lorenz_system, t_0, y_0, n, result, tolerance, t_end, initial_step, min_step);

    if (status != RKF45_SUCCESS) {
        if (status == RKF45_ERR_MEMORY) {
            fprintf(stderr, "Memory allocation error occurred during integration.\n");
        } else if (status == RKF45_ERR_MIN_STEP) {
            fprintf(stderr, "Step size fell below the minimum allowed value.\n");
        } else {
            fprintf(stderr, "An unknown error occurred.\n");
        }
        // メモリ解放
        for (size_t i = 0; i < 100000; i++) {
            if (result[i] != NULL) {
                free(result[i]);
            }
        }
        free(result);
        return EXIT_FAILURE;
    }

    // CSVファイルに結果を保存
    FILE *file = fopen("lorenz_output.csv", "w");
    if (file == NULL) {
        fprintf(stderr, "Failed to open output file.\n");
        // メモリ解放
        for (size_t i = 0; i < 100000; i++) {
            if (result[i] != NULL) {
                free(result[i]);
            }
        }
        free(result);
        return EXIT_FAILURE;
    }

    fprintf(file, "t,x,y,z\n");
    for (size_t i = 0; result[i] != NULL; i++) {
        fprintf(file, "%f,%f,%f,%f\n", result[i][0], result[i][1], result[i][2], result[i][3]);
        free(result[i]);
    }

    fclose(file);
    free(result);

    printf("Lorenz system solved successfully. Results saved to lorenz_output.csv\n");
    return EXIT_SUCCESS;
}
