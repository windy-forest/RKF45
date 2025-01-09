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

#endif // RKF45_H
