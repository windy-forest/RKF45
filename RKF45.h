#ifndef RKF45_H
#define RKF45_H

#include <stddef.h>
#include <stdbool.h>

// エラーコード定義
#define RKF45_SUCCESS 0
#define RKF45_ERR_MEMORY 1
#define RKF45_ERR_MIN_STEP 2

// ODEシステムを定義するコールバック関数の型
typedef void (*ODEFunction)(double t, const double *y, double *dydt);

/**
 * Runge-Kutta-Fehlberg法 (RK45) を用いて常微分方程式を数値的に解きます。
 *
 * @param func 微分方程式を定義する関数
 * @param t_0 初期時刻
 * @param y_0 初期条件配列
 * @param n 方程式の次元
 * @param result 解を格納する2次元配列のポインタ（呼び出し元でメモリ確保が必要）
 * @param tolerance 許容誤差
 * @param t_end 終了時刻
 * @param initial_step 初期ステップ幅
 * @param min_step 最小ステップ幅
 * @return 成功時にRKF45_SUCCESS、エラー時に対応するエラーコードを返します。
 */
int RKF45(
    ODEFunction func,
    double t_0,
    const double *y_0,
    size_t n,
    double **result,
    double tolerance,
    double t_end,
    double initial_step,
    double min_step
);

typedef struct model{
    double* (*dydt)(double t, double *y); // 勾配を計算する関数
    double *t_span; // 初期時刻と最終時刻
    int size; // 変数の数
    double *y_0; // 初期値
    double t; // 時刻
    double h; // 現在の刻み幅
    double tolerance; // 1ステップあたりの収束判定
    double initialstep; // 最初の刻み幅
    double min_step;
    double max_step;
    int calc_step;
} Model;

void init_model(Model *model, FILE output, double *dydt, double *t_span, int size, double *y_0, double tolerance, double initialstep, double min_step, double max_step){
    model->dydt = dydt;
    model->t_span = t_span;
    model->size = size;
    model->y_0 = y_0;
    model->t = 0.0;
    model->calc_step = 0;
    model->tolerance = tolerance;
    model->initialstep = initialstep;
    model->min_step = min_step;
    model->max_step = max_step;
}

void rkf45(Model *model){
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

    // kを確保する
    double *k[6];
    for (int i = 0; i < 6; i++) {
        k[i] = malloc(model->size * sizeof(double));
        if (k[i] == NULL) {
            for (int j = 0; j < i; j++) {
                free(k[j]);
            }
            free(model->y_0);
            return RKF45_ERR_MEMORY;
        }
    }

    while (model->t < model->t_span[1]) {
        // 必要なメモリ領域を確保
        double *next_y = (double *)malloc(model->size * sizeof(double));
        if (next_y == NULL) {
            for (int i = 0; i < model->size; i++) {
                free(next_y);
            }
            free(next_y);
            for (int i = 0; i < 6; i++) {
                free(k[i]);
            }
            return RKF45_ERR_MEMORY;
        }

        double t = model->t_span[0];
        double *y = malloc(model->size * sizeof(double));
        if (y == NULL) {
            return RKF45_ERR_MEMORY;
        }

        // kを計算
        k[0] = model->dydt(model->t, model->y_0);
        for (int i = 1; i < 6; i++) {
            double y_temp[model->size];
            for (int j = 0; j < model->size; j++) {
                y_temp[j] = y[j];
                for (int l = 0; l < i; l++) {
                    y_temp[j] += model->h * b[i][l] * k[l][j];
                }
            }
            func(t + a[i] * model->h, y_temp, k[i]);
        }

        // 次の値を計算
        for (int i = 0; i < model->size; i++) {
            next_y[i] = y[i];
            for (int j = 0; j < 6; j++) {
            next_y[i] += model->h * c[j] * k[j][i];
            }
        }

        // 誤差を計算
        double error = 0.0;
        for (int i = 0; i < model->size; i++) {
            double local_error = 0.0;
            for (int j = 0; j < 6; j++) {
            local_error += model->h * (c[j] - c_hat[j]) * k[j][i];
            }
            error = fmax(error, fabs(local_error));
        }

        // ステップ幅を調整
        if (error > model->tolerance) {
            model->h *= 0.9 * pow(model->tolerance / error, 0.25);
            if (model->h < model->min_step) {
            for (int i = 0; i < 6; i++) {
                free(k[i]);
            }
            free(next_y);
            free(y);
            return RKF45_ERR_MIN_STEP;
            }
            continue;
        }

        // 時刻と状態を更新
        model->t += model->h;
        for (int i = 0; i < model->size; i++) {
            y[i] = next_y[i];
        }

        // ステップ幅を増加
        model->h *= 0.9 * pow(model->tolerance / error, 0.2);
        if (model->h > model->max_step) {
            model->h = model->max_step;
        }

        free(next_y);
        
        for (int i = 0; i < 6; i++) {
        free(k[i]);
        }
        free(y);
    }
    return RKF45_SUCCESS;
}

#endif // RKF45_H
