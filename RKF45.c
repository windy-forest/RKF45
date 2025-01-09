#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

// エラーコード
#define RKF45_SUCCESS 0
#define RKF45_ERR_MEMORY 1
#define RKF45_ERR_MIN_STEP 2

typedef void (*ODEFunction)(double t, const double *y, double *dydt);

int RKF45(
    ODEFunction func,
    double t_0,
    const double *y_0,
    size_t n,
    double **result,
    double tolerance,
    double t_end,
    double initial_step, // 初期ステップ幅
    double min_step      // 最小ステップ幅
) {
    const double a[6] = {0.0, 0.25, 3.0 / 8.0, 12.0 / 13.0, 1.0, 0.5};
    const double b[6][5] = {
        {0},
        {0.25},
        {3.0 / 32.0, 9.0 / 32.0},
        {1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0},
        {439.0 / 216.0, -8.0, 3680.0 / 513.0, -845.0 / 4104.0},
        {-8.0 / 27.0, 2.0, -3544.0 / 2565.0, 1859.0 / 4104.0, -11.0 / 40.0}
    };
    const double c[6] = {16.0 / 135.0, 0.0, 6656.0 / 12825.0, 28561.0 / 56430.0, -9.0 / 50.0, 2.0 / 55.0};
    const double c_hat[6] = {25.0 / 216.0, 0.0, 1408.0 / 2565.0, 2197.0 / 4104.0, -1.0 / 5.0};

    double t = t_0;
    double *y = malloc(n * sizeof(double));
    if (y == NULL) {
        return RKF45_ERR_MEMORY;
    }

    double *k[6];
    for (int i = 0; i < 6; i++) {
        k[i] = malloc(n * sizeof(double));
        if (k[i] == NULL) {
            for (int j = 0; j < i; j++) {
                free(k[j]);
            }
            free(y);
            return RKF45_ERR_MEMORY;
        }
    }

    for (size_t i = 0; i < n; i++) {
        y[i] = y_0[i];
    }

    size_t step = 0;
    double h = initial_step;

    while (t < t_end) {
        // 必要なメモリ領域を確保
        result[step] = malloc((n + 1) * sizeof(double));
        if (result[step] == NULL) {
            for (size_t i = 0; i < step; i++) {
                free(result[i]);
            }
            free(y);
            for (int i = 0; i < 6; i++) {
                free(k[i]);
            }
            return RKF45_ERR_MEMORY;
        }

        result[step][0] = t;
        for (size_t i = 0; i < n; i++) {
            result[step][i + 1] = y[i];
        }
        step++;

        // kを計算
        func(t, y, k[0]);
        for (int i = 1; i < 6; i++) {
            double y_temp[n];
            for (size_t j = 0; j < n; j++) {
                y_temp[j] = y[j];
                for (int l = 0; l < i; l++) {
                    y_temp[j] += h * b[i][l] * k[l][j];
                }
            }
            func(t + a[i] * h, y_temp, k[i]);
        }

        // 次の値を計算
        double y_next[n];
        double error = 0.0;
        for (size_t i = 0; i < n; i++) {
            y_next[i] = y[i];
            double y_hat = y[i];
            for (int j = 0; j < 6; j++) {
                y_next[i] += h * c[j] * k[j][i];
                y_hat += h * c_hat[j] * k[j][i];
            }
            error = fmax(error, fabs(y_next[i] - y_hat));
        }

        // 誤差を基準にステップ幅を調整
        if (error <= tolerance) {
            t += h;
            for (size_t i = 0; i < n; i++) {
                y[i] = y_next[i];
            }
        }

        double safety_factor = 0.9;
        h *= safety_factor * pow(tolerance / error, 0.25);

        // ステップ幅が最小ステップ幅を下回った場合の処理
        if (fabs(h) < min_step) {
            free(y);
            for (int i = 0; i < 6; i++) {
                free(k[i]);
            }
            return RKF45_ERR_MIN_STEP;
        }

        // ステップ幅が終了時刻を超えないように調整
        if (t + h > t_end) {
            h = t_end - t;
        }
    }

    free(y);
    for (int i = 0; i < 6; i++) {
        free(k[i]);
    }

    return RKF45_SUCCESS;
}
