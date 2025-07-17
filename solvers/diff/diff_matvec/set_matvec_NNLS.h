#pragma once

//残差ベクトルのみをNNLSに使う
void ddhr_set_matvec_RH_for_NNLS_para_only_residuals(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        LPOD_MATRIX*     lpod_mat,
        POD_VALUES*     pod_vals,
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
        LPOD_MATRIX*    lpod_mat,
        POD_VALUES*     pod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int 		num_subdomains,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_neib,			//1 + monolis_com->recv_n_neib
        const double    dt,
		double       	t);

//残差ベクトルのみをNNLSに使う
void ddhr_set_matvec_only_residuals_for_NNLS2(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
    	BBFE_BC*     	bc,
        POD_MATRIX*     pod_mat,
        POD_VALUES*     pod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int 		num_subdomains,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_modes,
        const double    dt,
		double       	t);

//残差ベクトルのみをNNLSに使う
void ddhr_set_matvec_RH_for_NNLS2_only_residuals(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        POD_MATRIX*     pod_mat,
        POD_VALUES*     pod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int		num_subdomains,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_modes,
        const double    dt,
		double       	t);

void ddhr_set_matvec_RH_for_NNLS_para(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        LPOD_MATRIX*     lpod_mat,
        POD_VALUES*     pod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int		num_subdomains,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_modes,
        const double    dt,
		double       	t);
    
void ddhr_set_matvec_residuals_for_NNLS_para(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
    	BBFE_BC*     	bc,
        LPOD_MATRIX*    lpod_mat,
        POD_VALUES*     pod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int 		num_subdomains,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_neib,
        const double    dt,
		double       	t);

void ddhr_set_matvec_RH_for_NNLS2(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        POD_MATRIX*     pod_mat,
        POD_VALUES*     pod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int		num_subdomains,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_modes,
        const double    dt,
		double       	t);

void ddhr_set_matvec_residuals_for_NNLS2(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
    	BBFE_BC*     	bc,
        POD_MATRIX*     pod_mat,
        POD_VALUES*     pod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int 		num_subdomains,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_modes,
        const double    dt,
		double       	t);

void hr_set_matvec_for_NNLS(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        POD_MATRIX*     pod_mat,
        POD_VALUES*     pod_vals,
        HLPOD_HR*       hlpod_hr,
        const int       index_snap,
        const int       num_modes,
        const double    dt,
		double       	t);
    
void hr_set_matvec_RH_for_NNLS(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        POD_MATRIX*     pod_mat,
        POD_VALUES*     pod_vals,
        HLPOD_HR*       hlpod_hr,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_modes,
        const double    dt,
		double       	t);

void hr_set_matvec_residuals_for_NNLS(
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
    	BBFE_BC*     	bc,
        POD_MATRIX*     pod_mat,
        POD_VALUES*     pod_vals,
        HLPOD_HR*       hlpod_hr,
        const int       index_snap,
        const int       num_snapshot,
        const int       num_modes,
        const double    dt,
		double       	t);