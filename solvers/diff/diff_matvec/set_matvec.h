#pragma once

#include "diff_dataset.h"

#include "core_FOM.h"
#include "read_BB.h"

void set_element_reduced_rhs(
	double*			vec,
	double**		pod_basis,
	const int		index,
	const double	integ_val,
	const int		num_modes);

void set_element_reduced_rhs_Dbc(
	BBFE_BC*     	bc,
	double*			vec,
	double**		pod_basis,
	const int		index_i,
	const int		index_j,
	const double	integ_val,
	const int		num_modes);

void set_element_reduced_mat(
	double**		mat,
	double**		pod_basis,
	const int		index_i,
	const int		index_j,
	const double	integ_val,
	const int		num_modes);


void ddhr_set_reduced_mat_para_debug(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*    hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt);

void hr_set_reduced_mat(
		MONOLIS*     monolis,
		BBFE_DATA*     fe,
		BBFE_BASIS* basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*     hlpod_mat,
        HLPOD_HR*       hlpod_hr,
        const int num_modes,
		const double    dt);

void hr_set_reduced_vec(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HLPOD_HR*       hlpod_hr,
    	HLPOD_MAT*     hlpod_mat,
        const int		num_modes,
        const double    dt,
		double       	t);

void hr_set_D_bc(
		MONOLIS*     monolis,
		BBFE_DATA*     fe,
		BBFE_BASIS* basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*     hlpod_mat,
        HLPOD_HR*       hlpod_hr,
        const int num_modes,
		const double    dt);

void ddhr_set_reduced_mat2(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*     hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt);

void ddhr_set_D_bc2(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*     hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int		num_modes,
		const int 		num_subdomains,
		const double    dt);

void ddhr_set_reduced_vec3(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HLPOD_DDHR*     hlpod_ddhr,
    	HLPOD_MAT*     hlpod_mat,
        const int		num_modes,
		const int 		num_subdomains,
        const double    dt,
		double       	t);
	
void ddhr_set_reduced_mat3(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*     hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt);

void ddhr_set_D_bc3(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*     hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int		num_modes,
		const int 		num_subdomains,
		const double    dt);

void ddhr_set_reduced_vec3(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HLPOD_DDHR*     hlpod_ddhr,
    	HLPOD_MAT*     hlpod_mat,
        const int		num_modes,
		const int 		num_subdomains,
        const double    dt,
		double       	t);

/*
void ddhr_set_reduced_mat_para(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*    hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt);

void ddhr_set_D_bc_para(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*    hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int		num_modes,
		const int 		num_subdomains,
		const double    dt);

void ddhr_set_reduced_vec_para(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HLPOD_DDHR*     hlpod_ddhr,
    	HLPOD_MAT*    hlpod_mat,
        const int		num_modes,
		const int 		num_subdomains,
        const double    dt,
		double       	t);
*/

void ddhr_lb_set_reduced_mat_para(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*    hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt);

void ddhr_lb_read_reduced_mat(
    double**		reduced_mat,
	const int		m,
    const int 		n,
    const char*		directory);
	
void ddhr_lb_set_reduced_mat_para_save_memory(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
		HLPOD_VALUES*     hlpod_vals,
    	HLPOD_MAT*    hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int 		num_modes,
		const int 		num_subdomains,
		const double    dt);

void ddhr_lb_set_D_bc_para(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*    hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int		num_modes,
		const int 		num_subdomains,
		const double    dt);

void ddhr_lb_set_reduced_vec_para(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HLPOD_DDHR*     hlpod_ddhr,
    	HLPOD_MAT*     hlpod_mat,
        const int		num_modes,
		const int 		num_subdomains,
        const double    dt,
		double       	t);

/*for global mode + parallel computation*/
void set_reduced_mat_global_para(
		MONOLIS*		monolis,
		MONOLIS_COM*	monolis_com,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*     hlpod_mat,
        const int 		num_modes,
		const double    dt);

void set_D_bc_global_para(
		MONOLIS*     	monolis,
		MONOLIS_COM*	monolis_com,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*     hlpod_mat,
        const int 		num_modes,
		const double    dt);
	
void set_reduced_vec_global_para(
		MONOLIS*     	monolis,
		MONOLIS_COM*	monolis_com,
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
		HLPOD_VALUES*		hlpod_vals,
    	HLPOD_MAT*     hlpod_mat,
        const int		num_modes,
        const double    dt,
		double       	t);
