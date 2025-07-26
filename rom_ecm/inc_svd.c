#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <stdio.h>

#include "DDHR_para_lb.h"   /* 依存型や関数がここに実装されている想定 */
#include "monolis_wrapper_scalapack_c.h"

#include "inc_svd.h"

/* ============================================================
 *  このファイル内で使うユーティリティ（CMat）
 * ============================================================ */
static inline CMat cm_alloc(int rows, int cols){
  CMat M = {rows, cols, (double*)calloc((size_t)rows*cols, sizeof(double))};
  return M;
}
static inline void cm_free(CMat *M){
  free(M->a); M->a = NULL; M->rows = M->cols = 0;
}
static inline double cm_get(const CMat *M, int i, int j){
  return M->a[(size_t)i + (size_t)j*M->rows];
}
static inline void cm_set(CMat *M, int i, int j, double v){
  M->a[(size_t)i + (size_t)j*M->rows] = v;
}
static inline void cm_copy(const CMat *A, CMat *B){
  assert(A->rows==B->rows && A->cols==B->cols);
  memcpy(B->a, A->a, (size_t)A->rows*A->cols*sizeof(double));
}

/* C = A * B  (A: m×k, B: k×n, C: m×n) column-major */
static CMat matmul_nn(const CMat *A, const CMat *B){
  assert(A->cols == B->rows);
  int m=A->rows, k=A->cols, n=B->cols;
  CMat C = cm_alloc(m,n);
  for(int j=0;j<n;++j){
    for(int p=0;p<k;++p){
      double bpj = cm_get(B,p,j);
      const double *ai = &A->a[(size_t)p*A->rows];
      double *ci = &C.a[(size_t)j*C.rows];
      for(int i=0;i<m;++i) ci[i] += ai[i]*bpj;
    }
  }
  return C;
}

/* C = A^T * B  (A: m×k, B: m×n -> C: k×n) */
static CMat matmul_tn(const CMat *A, const CMat *B){
  assert(A->rows == B->rows);
  int m=A->rows, k=A->cols, n=B->cols;
  CMat C = cm_alloc(k,n);
  for(int j=0;j<n;++j){
    for(int p=0;p<k;++p){
      double s=0.0;
      for(int i=0;i<m;++i) s += cm_get(A,i,p)*cm_get(B,i,j);
      cm_set(&C,p,j,s);
    }
  }
  return C;
}

/* C = A - B  */
static CMat mat_sub(const CMat *A, const CMat *B){
  assert(A->rows==B->rows && A->cols==B->cols);
  CMat C = cm_alloc(A->rows, A->cols);
  int n=A->rows*A->cols;
  for(int t=0;t<n;++t) C.a[t] = A->a[t] - B->a[t];
  return C;
}

/* C = [A B]  (横結合, 行同一) */
static CMat concat_cols(const CMat *A, const CMat *B){
  assert(A->rows == B->rows);
  CMat C = cm_alloc(A->rows, A->cols + B->cols);
  for(int j=0;j<A->cols;++j)
    memcpy(&C.a[(size_t)j*C.rows], &A->a[(size_t)j*A->rows], (size_t)A->rows*sizeof(double));
  for(int j=0;j<B->cols;++j)
    memcpy(&C.a[(size_t)(j+A->cols)*C.rows], &B->a[(size_t)j*B->rows], (size_t)B->rows*sizeof(double));
  return C;
}

/* C = diag(A, B)  (ブロック対角) */
static CMat block_diag(const CMat *A, const CMat *B){
  CMat C = cm_alloc(A->rows + B->rows, A->cols + B->cols);
  for(int j=0;j<A->cols;++j)
    for(int i=0;i<A->rows;++i)
      cm_set(&C,i,j, cm_get(A,i,j));
  for(int j=0;j<B->cols;++j)
    for(int i=0;i<B->rows;++i)
      cm_set(&C, A->rows+i, A->cols+j, cm_get(B,i,j));
  return C;
}

/* I_n */
static CMat eye(int n){
  CMat Imat = cm_alloc(n,n); /* 変数名 I は complex.h の I と衝突するため使用しない */
  for(int i=0;i<n;++i) cm_set(&Imat, i, i, 1.0);
  return Imat;
}

/* 0 行列 */
static CMat zeros(int m, int n){ return cm_alloc(m,n); }

/* diag(sigma) */
static CMat diag_from_vector(const double *sigma, int r){
  CMat S = cm_alloc(r,r);
  for(int i=0;i<r;++i) cm_set(&S, i, i, sigma[i]);
  return S;
}

/* 2x2 ブロック結合 */
static CMat block_2x2(const CMat *A11, const CMat *A12, const CMat *A21, const CMat *A22){
  assert(A11->rows == A12->rows);
  assert(A21->rows == A22->rows);
  assert(A11->cols == A21->cols);
  assert(A12->cols == A22->cols);
  int r1=A11->rows, c1=A11->cols, r2=A22->rows, c2=A22->cols;
  CMat K = cm_alloc(r1+r2, c1+c2);
  for(int j=0;j<c1;++j)
    for(int i=0;i<r1;++i)
      cm_set(&K,i,j, cm_get(A11,i,j));
  for(int j=0;j<c2;++j)
    for(int i=0;i<r1;++i)
      cm_set(&K,i, c1+j, cm_get(A12,i,j));
  for(int j=0;j<c1;++j)
    for(int i=0;i<r2;++i)
      cm_set(&K, r1+i, j, cm_get(A21,i,j));
  for(int j=0;j<c2;++j)
    for(int i=0;i<r2;++i)
      cm_set(&K, r1+i, c1+j, cm_get(A22,i,j));
  return K;
}

/* ============================================================
 *  QR（MGS2 再直交化）
 * ============================================================ */
static void qr_mgs2(const CMat *E, CMat *Q, CMat *R, double tol){
  int K=E->rows, D=E->cols;
  *Q = cm_alloc(K, D);
  *R = cm_alloc(D, D);
  int q=0;
  for(int j=0;j<D;++j){
    for(int i=0;i<K;++i) cm_set(Q, i, q, cm_get(E,i,j));
    for(int k=0;k<q;++k){
      double r = 0.0;
      for(int i=0;i<K;++i) r += cm_get(Q,i,k)*cm_get(Q,i,q);
      cm_set(R,k,j, r);
      for(int i=0;i<K;++i) cm_set(Q, i, q, cm_get(Q,i,q) - r*cm_get(Q,i,k));
    }
    for(int k=0;k<q;++k){
      double r = 0.0;
      for(int i=0;i<K;++i) r += cm_get(Q,i,k)*cm_get(Q,i,q);
      cm_set(R,k,j, cm_get(R,k,j)+r);
      for(int i=0;i<K;++i) cm_set(Q, i, q, cm_get(Q,i,q) - r*cm_get(Q,i,k));
    }
    double nrm=0.0;
    for(int i=0;i<K;++i){ double v=cm_get(Q,i,q); nrm += v*v; }
    nrm = sqrt(nrm);
    if(nrm > tol){
      cm_set(R,q,j, nrm);
      for(int i=0;i<K;++i) cm_set(Q,i,q, cm_get(Q,i,q)/nrm);
      ++q;
    }else{
      /* rank 落ち */
    }
  }
  Q->cols = q;
  R->rows = q;
}

/* ============================================================
 *  小 SVD： LAPACK あり→dgesvd / 無し→Jacobi 近似
 * ============================================================ */
#ifdef HAVE_LAPACK
  #include <lapacke.h>
#endif

/* 対称固有分解（Jacobi）: Gin = Q Λ Q^T */
static void jacobi_eigensymm(const CMat *Gin, double *lam, CMat *Q, int max_sweeps, double tol){
  int n = Gin->rows; 
  CMat G = cm_alloc(n,n); cm_copy(Gin,&G);
  *Q = eye(n);
  for(int sweep=0; sweep<max_sweeps; ++sweep){
    double off=0.0;
    for(int j=0;j<n;++j) for(int i=0;i<j;++i){ double g=cm_get(&G,i,j); off+=g*g; }
    if(sqrt(off) < tol) break;
    for(int p=0;p<n-1;++p){
      for(int q=p+1;q<n;++q){
        double gpq=cm_get(&G,p,q); if(fabs(gpq)<=0.0) continue;
        double gpp=cm_get(&G,p,p), gqq=cm_get(&G,q,q);
        double tau=(gqq-gpp)/(2.0*gpq);
        double t = (tau>=0)? 1.0/(tau+sqrt(1.0+tau*tau)) : -1.0/(-tau+sqrt(1.0+tau*tau));
        double c=1.0/sqrt(1.0+t*t), s=t*c;
        for(int k=0;k<n;++k){
          if(k==p||k==q) continue;
          double gpk=cm_get(&G,p,k), gqk=cm_get(&G,q,k);
          double npk=c*gpk - s*gqk, nqk=s*gpk + c*gqk;
          cm_set(&G,p,k,npk); cm_set(&G,k,p,npk);
          cm_set(&G,q,k,nqk); cm_set(&G,k,q,nqk);
        }
        double gpp_new=c*c*gpp - 2*c*s*gpq + s*s*gqq;
        double gqq_new=s*s*gpp + 2*c*s*gpq + c*c*gqq;
        cm_set(&G,p,p,gpp_new); cm_set(&G,q,q,gqq_new);
        cm_set(&G,p,q,0.0); cm_set(&G,q,p,0.0);
        for(int i=0;i<n;++i){
          double vip=cm_get(Q,i,p), viq=cm_get(Q,i,q);
          cm_set(Q,i,p, c*vip - s*viq);
          cm_set(Q,i,q, s*vip + c*viq);
        }
      }
    }
  }
  for(int i=0;i<n;++i) lam[i]=cm_get(&G,i,i);
  cm_free(&G);
}

/* 固有値降順ソート（Q の列も並べ替え） */
static void sort_eigs_desc(double *lam, CMat *Q){
  int n=Q->cols;
  for(int i=0;i<n-1;++i){
    int p=i; for(int j=i+1;j<n;++j) if(lam[j]>lam[p]) p=j;
    if(p!=i){
      double t=lam[i]; lam[i]=lam[p]; lam[p]=t;
      for(int r=0;r<Q->rows;++r){
        double tmp=cm_get(Q,r,i);
        cm_set(Q,r,i, cm_get(Q,r,p));
        cm_set(Q,r,p, tmp);
      }
    }
  }
}

/* Ksmall = Uhat * diag(S) * Vhat^T を求める */
static void svd_dense(const CMat *Ksmall, CMat *Uhat, double **sigma_out, CMat *Vhat){
  int m=Ksmall->rows, n=Ksmall->cols, r=(m<n?m:n);
#ifdef HAVE_LAPACK
  CMat A = cm_alloc(m,n); cm_copy(Ksmall,&A);
  *Uhat = cm_alloc(m,m);
  *Vhat = cm_alloc(n,n);
  double *S = (double*)calloc((size_t)r,sizeof(double));
  int info = LAPACKE_dgesvd(LAPACK_COL_MAJOR,'A','A', m,n, A.a,m, S, Uhat->a,m, Vhat->a,n, NULL);
  assert(info==0);
  *sigma_out=S; cm_free(&A);
#else
  /* Vhat と σ：K^T K の固有分解 */
  CMat G = matmul_tn(Ksmall, Ksmall);      /* n×n = K^T K */
  CMat Q; double *lam = (double*)calloc((size_t)n,sizeof(double));
  jacobi_eigensymm(&G, lam, &Q, 50, 1e-12); cm_free(&G);
  sort_eigs_desc(lam, &Q);
  double *S = (double*)calloc((size_t)r,sizeof(double));
  for(int i=0;i<r;++i){ S[i] = lam[i]>0.0? sqrt(lam[i]) : 0.0; }
  *Vhat = Q; /* n×n */
  /* Uhat = K * Vhat * Σ^{-1} （上位 m 列を構成） */
  CMat KV  = matmul_nn(Ksmall, Vhat);      /* m×n */
  *Uhat = cm_alloc(m,m);
  for(int j=0;j<m;++j){
    double inv = (j<r && S[j]>1e-15)? 1.0/S[j] : 0.0;
    for(int i=0;i<m;++i){
      double v = cm_get(&KV,i,j) * inv;
      cm_set(Uhat,i,j, v);
    }
  }
  cm_free(&KV); free(lam);
  *sigma_out = S;
#endif
}

/* ============================================================
 *  ランク選択 & トリム
 * ============================================================ */
static int select_rank(const double *sigma, int len, int r_max, double tol){
  int r = 0; double s0 = (len>0? sigma[0]:0.0);
  for(int i=0;i<len && i<r_max; ++i){
    if(sigma[i] >= tol*s0) ++r; else break;
  }
  if(r==0 && len>0) r=1;
  return r;
}

/* U_new: K×(r+q) / sigma_new: ≥min(r+q, r+Δ) / V_new: (N+Δ)×(r+Δ) */
static void truncate_uv(CMat *U_new, double **sigma_new, CMat *V_new, int r_keep){
  int m=U_new->rows, nV=V_new->rows;
  /* U */
  CMat U2 = cm_alloc(m, r_keep);
  for(int j=0;j<r_keep;++j)
    memcpy(&U2.a[(size_t)j*U2.rows], &U_new->a[(size_t)j*U_new->rows], (size_t)m*sizeof(double));
  cm_free(U_new); *U_new = U2;
  /* sigma */
  double *s2 = (double*)calloc((size_t)r_keep,sizeof(double));
  memcpy(s2, *sigma_new, (size_t)r_keep*sizeof(double));
  free(*sigma_new); *sigma_new = s2;
  /* V */
  CMat V2 = cm_alloc(nV, r_keep);
  for(int j=0;j<r_keep;++j)
    memcpy(&V2.a[(size_t)j*V2.rows], &V_new->a[(size_t)j*V_new->rows], (size_t)nV*sizeof(double));
  cm_free(V_new); *V_new = V2;
}

/* ============================================================
 *  増分 SVD 更新（Brand 型）
 * ============================================================ */
void incsvd_update(IncSVD* S, const CMat *B, int r_max, double tol){
  const int K = S->K;
  const int Delta = B->cols;
  const int r = S->r;

  /* (1) C = U^T B */
  CMat U = (CMat){K, r, S->U};
  CMat C = (r>0) ? matmul_tn(&U, B) : cm_alloc(0, Delta); /* r×Δ or 0×Δ */

  /* (2) E = B - U C, QR */
  CMat UC = (r>0) ? matmul_nn(&U, &C) : cm_alloc(K, Delta);
  CMat E  = mat_sub(B, &UC); cm_free(&UC);
  CMat Q,R; qr_mgs2(&E, &Q, &R, 1e-14); cm_free(&E);
  int q = Q.cols;

  /* (3) 小行列 SVD */
  CMat Sdiag = diag_from_vector(S->sigma, r);
  CMat K11 = Sdiag;         /* r×r */
  CMat K12 = C;             /* r×Δ */
  CMat K21 = zeros(q, r);
  CMat K22 = R;             /* q×Δ */
  CMat Ksmall = block_2x2(&K11,&K12,&K21,&K22);

  CMat Uhat, Vhat; double *sigma_new=NULL;
  svd_dense(&Ksmall, &Uhat, &sigma_new, &Vhat);
  cm_free(&K11); cm_free(&K21); cm_free(&K22); cm_free(&Ksmall);

  /* (4) U,V 更新 */
  CMat Ucm  = (CMat){K, r, S->U};
  CMat U_aug = (q>0) ? concat_cols(&Ucm, &Q) : Ucm;  /* コピー生成 */
  CMat U_new = matmul_nn(&U_aug, &Uhat);

  CMat Vcm  = (CMat){S->N, r, S->V};
  CMat I_delta = eye(Delta);
  CMat V_aug = block_diag(&Vcm, &I_delta);
  CMat V_new = matmul_nn(&V_aug, &Vhat);

  /* (5) ランク選択＆トリム */
  int len_sigma = (Uhat.rows < Vhat.cols) ? Uhat.rows : Vhat.cols;
  int r_keep = select_rank(sigma_new, len_sigma, r_max, tol);
  truncate_uv(&U_new, &sigma_new, &V_new, r_keep);

  /* (6) 状態に反映 */
  free(S->U); free(S->V); free(S->sigma);
  S->U = U_new.a; S->V = V_new.a; S->sigma = sigma_new;
  S->r = r_keep; S->N += Delta;

  /* 片付け */
  if(q>0){ cm_free(&Q); }
  cm_free(&C); cm_free(&Uhat); cm_free(&Vhat); cm_free(&I_delta); cm_free(&V_aug);
}

/* ============================================================
 *  初期ブロック：monolis で SVD → IncSVD に格納
 * ============================================================ */
static void init_with_first_block_by_monolis(
    IncSVD *S,            /* out: U:K×r, V:N0×r, sigma:r */
    double **matrixRM,    /* in : 全体 A (row-major) */
    int K, int N0,        /* A0 のサイズ */
    int r_max, double tol,
    int comm, int scalapack_comm)
{
    /* A0 を切り出し（row-major） */
    double **A0 = BB_std_calloc_2d_double(A0, K, N0);
    for (int j=0;j<K;++j){
        for (int i=0;i<N0;++i){
            A0[j][i] = matrixRM[j][i];
        }
    }

    /* monolis で SVD:  A0 = U * Σ * (V^T) */
    double **U = BB_std_calloc_2d_double(U, K, N0);
    double **V_or_Vt = BB_std_calloc_2d_double(V_or_Vt, N0, N0); /* 返却仕様に依存 */
    double *D = BB_std_calloc_1d_double(D, N0);       /* Σ を対角に格納 */

    monolis_scalapack_gesvd_R(K, N0, A0, U, D, V_or_Vt, comm, scalapack_comm);


    /* ランク決定：特異値に基づく */
    int r_full = (K<N0?K:N0);
    int r = 0; double s0 = D[0];
    for (int i=0;i<r_full && i<r_max; ++i){
        if (D[i] >= tol*s0) ++r; else break;
    }
    if (r==0 && r_full>0) r=1;

    /* IncSVD へコピー（column-major に詰め替え） */
    S->K = K; S->N = N0; S->r = r;
    S->U = (double*)malloc((size_t)K*r*sizeof(double));
    S->V = (double*)malloc((size_t)N0*r*sizeof(double));
    S->sigma = (double*)malloc((size_t)r*sizeof(double));

    for (int j=0;j<r;++j){
        S->sigma[j] = D[j];
        for (int i=0;i<K;++i) S->U[i + (size_t)j*K] = U[i][j];
    }

    for (int j=0;j<r;++j)
      for (int i=0;i<N0;++i)
        S->V[i + (size_t)j*N0] = V_or_Vt[i][j];   /* V(i,j) */

    /* 後片付け */
    BB_std_free_2d_double(A0, K, N0);
    BB_std_free_2d_double(U,  K, N0);
    BB_std_free_2d_double(V_or_Vt, N0, N0);
    BB_std_free_1d_double(D,  N0);
}

/* Row-major(double**) → column-major(CMat) へのブロックコピー */
static CMat copy_block_colmajor_from_rowmajor(double **matrixRM, int K, int /*N*/, int j0, int Delta){
  CMat B = cm_alloc(K, Delta);
  for(int j=0;j<Delta;++j){
    int jj = j0 + j;
    for(int i=0;i<K;++i){
      B.a[(size_t)i + (size_t)j*K] = matrixRM[i][jj];
    }
  }
  return B;
}

/* row-major 配列のアロケータ＆解放（NNLS 入力用） */
static double** alloc2d_rowmajor(int rows, int cols){
  double **M = (double**)malloc((size_t)rows*sizeof(double*));
  for(int i=0;i<rows;++i) M[i] = (double*)calloc((size_t)cols,sizeof(double));
  return M;
}
static void free2d_rowmajor(double **M, int rows){
  for(int i=0;i<rows;++i) free(M[i]);
  free(M);
}

/* 圧縮系の構築： G_k = Σ V^T,  b_k = U^T b */
static void build_compressed_system_from_incsvd(
    const IncSVD* S, const double* b,
    double*** G_k_out, double** b_k_out)
{
  int r=S->r, K=S->K, N=S->N;
  double **G_k = alloc2d_rowmajor(r, N);
  for(int i=0;i<r;++i){
    double sii = S->sigma[i];
    for(int j=0;j<N;++j){
      double Vij = S->V[(size_t)j + (size_t)i*N];  /* V: N×r */
      G_k[i][j] = sii * Vij;                       /* Σ V^T */
    }
  }
  double *b_k = (double*)calloc((size_t)r,sizeof(double));
  for(int i=0;i<r;++i){
    double ssum=0.0;
    for(int p=0;p<K;++p) ssum += S->U[(size_t)p + (size_t)i*K] * b[p];
    b_k[i] = ssum;                                 /* U^T b */
  }
  *G_k_out = G_k;
  *b_k_out = b_k;
}

/* ---- モジュール内グローバル（サブドメインごとの SVD 状態） ---- */
static IncSVD *g_svd = NULL;      /* 長さ = g_svd_count */
static int     g_svd_count = 0;   /* = num_subdomains を保持 */
static int     g_block_cols = 512;/* BLOCK 幅（初期化と更新で共通に使う） */

static void ensure_svd_states(int num_subdomains){
    if (g_svd && g_svd_count == num_subdomains) return;
    /* 既存を破棄 */
    if (g_svd){
        for (int i=0;i<g_svd_count;++i){
            free(g_svd[i].U); free(g_svd[i].V); free(g_svd[i].sigma);
        }
        free(g_svd);
    }
    g_svd = (IncSVD*)calloc((size_t)num_subdomains, sizeof(IncSVD));
    g_svd_count = num_subdomains;
}

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
    const char*    directory)
{
    ensure_svd_states(num_subdomains);

    const int r_max   = total_num_modes; /* ランク上限 */
    const double svd_tol = 1.0e-15;

    const int comm = monolis_mpi_get_self_comm();
    int scalapack_comm; monolis_scalapack_comm_initialize(comm, &scalapack_comm);

    for (int m = 0; m < num_subdomains; ++m) {

        const int K = hlpod_ddhr->num_modes_1stdd[m] * total_num_snapshot + 1;
        const int N = hlpod_ddhr->num_elems[m];

        /* A 全体を row-major で作る（最小限：初期ブロック分だけ使う） */
        double **matrix = BB_std_calloc_2d_double(matrix, K, N);
        for (int j=0; j<K; ++j){
            for (int i=0; i<N; ++i)
                matrix[j][i] = hlpod_ddhr->matrix[j][i][m];
        }

        /* 初期ブロック SVD → g_svd[m] に保存 */
        init_with_first_block_by_monolis(&g_svd[m], matrix, K, N, r_max, svd_tol, comm, scalapack_comm);

        /* 行列はここでは不要になったので開放 */
        BB_std_free_2d_double(matrix, K, N);
    }
}


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
    const char*    directory)
{
    /* 初期化がまだなら安全のため作る（U,V=NULL のままでも incsvd_update は呼ばない） */
    ensure_svd_states(num_subdomains);

    const int r_max   = total_num_modes;
    const double svd_tol = 1.0e-15;

    for (int m = 0; m < num_subdomains; ++m) {

        const int K = hlpod_ddhr->num_modes_1stdd[m] * total_num_snapshot + 1;
        const int N = hlpod_ddhr->num_elems[m];
        //const int N0 = (N < g_block_cols ? N : g_block_cols);

        /* g_svd[m] が未初期化ならスキップ（先に init を呼んでください） */
        if (g_svd[m].U == NULL || g_svd[m].V == NULL || g_svd[m].sigma == NULL){
            fprintf(stderr, "[incsvd_update] subdomain %d: SVD state not initialized. Call init first.\n", m);
            continue;
        }

        /* A を row-major で作成（更新では全体のうち N0 以降の列だけ使う） */
        double **matrix = BB_std_calloc_2d_double(matrix, K, N);
        for (int j=0; j<K; ++j){
            for (int i=0; i<N; ++i)
                matrix[j][i] = hlpod_ddhr->matrix[j][i][m];
        }

        /* 残り列を BLOCK ごとに追加 */
        for (int j0 = N; j0 < N; j0 += g_block_cols){
            int Delta = (j0 + g_block_cols <= N ? g_block_cols : (N - j0));
            CMat B = copy_block_colmajor_from_rowmajor(matrix, K, N, j0, Delta);
            incsvd_update(&g_svd[m], &B, r_max, svd_tol);
            cm_free(&B);
        }

        BB_std_free_2d_double(matrix, K, N);
    }
}

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
    const char*    directory)
{
    const int max_ITER = (max_iter>0 ? max_iter : 1000);
    const double TOL = (tol>0 ? tol : 1.0e-10);

    /* ここでは g_svd[m] が「全列取り込み済み」である前提（init→update を済ませる） */
    for (int m = 0; m < num_subdomains; ++m) {

        if (g_svd[m].U == NULL || g_svd[m].V == NULL || g_svd[m].sigma == NULL){
            fprintf(stderr, "[finalize] subdomain %d: SVD state not ready.\n", m);
            continue;
        }

        const int K = g_svd[m].K;
        const int N = g_svd[m].N;

        /* 右辺 b を構築 */
        double *RH = BB_std_calloc_1d_double(RH, K);
        for (int j=0; j<K; ++j) RH[j] = hlpod_ddhr->RH[j][m];

        /* 圧縮系 G_k=ΣV^T, b_k=U^T b を作る */
        double **G_k = NULL; double *b_k = NULL;
        build_compressed_system_from_incsvd(&g_svd[m], RH, &G_k, &b_k);


        /* NNLS を解く */
        double *ans_vec = BB_std_calloc_1d_double(ans_vec, N);
        double residual = 0.0;
        monolis_optimize_nnls_R_with_sparse_solution(
            G_k, b_k, ans_vec, g_svd[m].r, hlpod_ddhr->num_elems[m], max_ITER, TOL, &residual);

        g_svd[m].K = g_svd[m].N = g_svd[m].r = 0;
    }

    /* すべて終わったら配列も解放（再利用するなら解放しない） */
    free(g_svd); g_svd=NULL; g_svd_count=0;
}
