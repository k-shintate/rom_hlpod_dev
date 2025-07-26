#ifndef INC_SVD_H_
#define INC_SVD_H_

#include <stddef.h>
#include "monolis_wrapper_scalapack_c.h"


/* ---- 増分 SVD 用の軽量行列（column-major） ---- */
typedef struct {
  int rows;     /* 行数 */
  int cols;     /* 列数 */
  double *a;    /* a[i + j*rows] = A(i,j) */
} CMat;

/* ---- 増分 SVD の保持状態 ---- */
typedef struct {
  int K;          /* 行数 */
  int N;          /* 現在の列数 */
  int r;          /* 現在のランク */
  double *U;      /* [K x r] column-major */
  double *V;      /* [N x r] column-major */
  double *sigma;  /* [r] 特異値 */
} IncSVD;

/* ---- API ---- */
/* 列ブロック B（K×Delta）を追加して U,Σ,V を更新（Brand 型） */
void incsvd_update(IncSVD* S, const CMat *B, int r_max, double tol);

/* Incremental SVD + 圧縮 NNLS 本体（ユーザ関数） */
void ddhr_lb_write_selected_elements_para_1line_init_with_first_block(
    MONOLIS_COM*   monolis_com,
    BBFE_DATA*     fe,
    BBFE_BC*       bc,
    HLPOD_VALUES*  hlpod_vals,
    HLPOD_DDHR*    hlpod_ddhr,
    HLPOD_MAT*     hlpod_mat,
    HLPOD_META*    hlpod_meta,
    const int      total_num_elem,
    const int      total_num_snapshot,
    const int      total_num_modes,
    const int      num_subdomains,
    const int      max_iter,   /* NNLS: 未使用 */
    const double   tol,        /* NNLS: 未使用 */
    const int      dof,
    const char*    directory);

void ddhr_lb_write_selected_elements_para_1line_incsvd_update(
    MONOLIS_COM*   monolis_com,
    BBFE_DATA*     fe,
    BBFE_BC*       bc,
    HLPOD_VALUES*  hlpod_vals,
    HLPOD_DDHR*    hlpod_ddhr,
    HLPOD_MAT*     hlpod_mat,
    HLPOD_META*    hlpod_meta,
    const int      total_num_elem,
    const int      total_num_snapshot,
    const int      total_num_modes,
    const int      num_subdomains,
    const int      max_iter,   /* NNLS: 未使用 */
    const double   tol,        /* NNLS: 未使用 */
    const int      dof,
    const char*    directory);

void ddhr_lb_write_selected_elements_para_1line_incremental_svd(
    MONOLIS_COM*   monolis_com,
    BBFE_DATA*     fe,
    BBFE_BC*       bc,
    HLPOD_VALUES*  hlpod_vals,
    HLPOD_DDHR*    hlpod_ddhr,
    HLPOD_MAT*     hlpod_mat,
    HLPOD_META*    hlpod_meta,
    const int      total_num_elem,
    const int      total_num_snapshot,
    const int      total_num_modes,
    const int      num_subdomains,
    const int      max_iter,   /* NNLS */
    const double   tol,        /* NNLS */
    const int      dof,
    const char*    directory);


#endif /* INC_SVD_H_ */
