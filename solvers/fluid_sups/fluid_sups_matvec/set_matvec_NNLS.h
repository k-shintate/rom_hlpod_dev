#pragma once

//残差ベクトルのみをNNLSに使う
void ddhr_set_matvec_RH_for_NNLS_para_only_residuals(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HLPOD_MAT*     hlpod_mat,
        HLPOD_VALUES*     hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int		num_subdomains,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_modes,
        const double    dt,
		double       	t);

//残差ベクトルのみをNNLSに使う
void ddhr_set_matvec_residuals_for_NNLS_para_only_residuals(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
    	BBFE_BC*     	bc,
        HLPOD_MAT*      hlpod_mat,
        HLPOD_VALUES*   hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int 		num_subdomains,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_neib,			//1 + monolis_com->recv_n_neib
        const double    dt,
		double       	t);
