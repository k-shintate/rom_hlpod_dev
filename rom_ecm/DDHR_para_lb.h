#pragma once

#include "monolis_nnls_c.h"
#include "hrom_dataset.h"
#include "hlpod_core_fe.h"
#include "write_std.h"
#include "write_BB.h"

void ddhr_lb_get_selected_elements_internal_overlap(
	    HLPOD_DDHR*     hlpod_ddhr,
//		const int		num_parallel_subdomains,
		const char*		directory);

void ddhr_lb_read_selected_elements_para(
		const int 		num_subdomains,
		const char*		directory);

void ddhr_lb_write_selected_elements_para(
		MONOLIS_COM*  	monolis_com,
        BBFE_DATA*     	fe,
        BBFE_BC*     	bc,
        HLPOD_VALUES* 	hlpod_vals,
	    HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_MAT*    hlpod_mat,
		HLPOD_META*		hlpod_meta,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
		const int 		num_subdomains,
        const int       max_iter, //NNLS
        const double    tol,      //NNLS
		const char*		directory);

//1列のみ(残差に関する項のみ：任意列数に拡張する必要あり)
void ddhr_lb_write_selected_elements_para_1line(
		MONOLIS_COM*  	monolis_com,
        BBFE_DATA*     	fe,
        BBFE_BC*     	bc,
        HLPOD_VALUES* 	hlpod_vals,
	    HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_MAT*    hlpod_mat,
		HLPOD_META*		hlpod_meta,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
		const int 		num_subdomains,
        const int       max_iter, //NNLS
        const double    tol,      //NNLS
		const char*		directory);

void get_meta_neib(
		MONOLIS_COM*  	monolis_com,
		HLPOD_META*		hlpod_meta,
		const char*     directory);

void ddhr_lb_set_neib(
		MONOLIS_COM*  	monolis_com,
		//BBFE_DATA*     	fe,
		HLPOD_MAT* 	hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_META*		hlpod_meta,
		const int 		num_subdomains,
        const int       num_snapshots,
		const char*     directory);

/*
void ddhr_lb_set_element_para(
		BBFE_DATA*     	fe,
        HLPOD_VALUES* 	hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
		const int 		num_subdomains,
		const char*     directory);
*/

void ddhr_lb_set_element_para2(
		BBFE_DATA*     	fe,
        HLPOD_DDHR*     hlpod_ddhr,
		const int 		num_subdomains,
		const char*     directory);

void ddhr_lb_get_selected_elements_para_add(
	    HLPOD_DDHR*     hlpod_ddhr,
		const int		num_parallel_subdomains,
		const char*		directory);

void ddhr_hlpod_calc_block_mat_bcsr_pad(
		MONOLIS*     	monolis,
		MONOLIS_COM*  	monolis_com,
        HLPOD_VALUES* 	hlpod_vals,
		//LPOD_COM* 		lpod_com,
		HLPOD_MAT* 	hlpod_mat,
		//LPOD_PRM*		lpod_prm,
		HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_META*		hlpod_meta,
		const int 		num_bases,
		const int		num_2nddd,
		const char*		directory);

void ddhr_hlpod_calc_block_mat_bcsr(
		MONOLIS*     	monolis,
		MONOLIS_COM*  	monolis_com,
        HLPOD_VALUES* 	hlpod_vals,
		//LPOD_COM* 		lpod_com,
		HLPOD_MAT* 	hlpod_mat,
		//LPOD_PRM*		lpod_prm,
		HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_META*		hlpod_meta,
		const int 		num_bases,
		const int		num_2nddd,
		const char*		directory);

void ddhr_hlpod_WTf_to_monollis_rhs_bcsr(
		MONOLIS*		monolis,
		MONOLIS_COM*	monolis_com,
		HLPOD_MAT*    hlpod_mat,
		//LPOD_COM*    	lpod_com,
		HLPOD_DDHR*     hlpod_ddhr,
		const int		num_bases);
	
void ddhr_lb_get_selected_elements_para(
        BBFE_DATA*     	fe,
        BBFE_BC*     	bc,
        HLPOD_VALUES* 	hlpod_vals,
	    HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_MAT*    hlpod_mat,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
		const int 		num_subdomains,
        const int       max_iter, //NNLS
        const double    tol,      //NNLS
		const char*		directory);

void ddhr_lb_get_selected_elements_para2(
		MONOLIS_COM*  	monolis_com,
        BBFE_DATA*     	fe,
        BBFE_BC*     	bc,
        HLPOD_VALUES* 	hlpod_vals,
	    HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_MAT*    hlpod_mat,
		HLPOD_META*		hlpod_meta,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
		const int 		num_subdomains,
        const int       max_iter, //NNLS
        const double    tol,      //NNLS
		const char*		directory);

void ddhr_lb_set_selected_elems_para(
		BBFE_DATA*     	fe,
        HLPOD_DDHR*     hlpod_ddhr,
        const int		total_num_nodes,
		const int		num_subdomains,
		const char*     directory);

//level1領域の最大基底本数の共有
void get_neib_max_num_modes_pad(
		MONOLIS_COM*  	monolis_com,
        HLPOD_VALUES* 	hlpod_vals,
		HLPOD_MAT* 	    hlpod_mat,
		const int       np,
		const int       pre_num_modes);

//level1領域の選択された基底(p-adaptive)本数の共有
void get_neib_num_modes_pad(
		MONOLIS_COM*  	monolis_com,
        HLPOD_VALUES* 	hlpod_vals,
		HLPOD_MAT* 	    hlpod_mat,
		const int       np,
		const int       num_my_modes);

//for arbit dof ddecm
void get_neib_subdomain_id(
        MONOLIS_COM*  	monolis_com,
        HLPOD_MAT* 	    hlpod_mat,
        const int 		num_modes);

void set_max_num_modes(
	HLPOD_VALUES*		hlpod_vals,
    const int       num_modes,
	const int       num_1st_dd,	//並列計算領域数
	const char*     directory);

void get_neib_coordinates_pre(
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_MAT*	    hlpod_mat,
    const int       np,				//並列計算領域数
	const int       max_num_basis);

void get_neib_coordinates_pad(
	MONOLIS_COM*  	monolis_com,
    HLPOD_VALUES* 	hlpod_vals,
	HLPOD_MAT*	    hlpod_mat,
    const int      np,				//並列計算領域数
	const int       max_num_basis,	//level 1の基底本数 (並列計算領域が担当する基底本数の総和)
	const int 		num_subdomains,
	const int		max_num_bases);


void hlpod_hr_sys_set_bc_id(
    BBFE_DATA* 	fe,
    BBFE_BC*   	bc,
    HLPOD_DDHR* hlpod_ddhr,
    const int   num_dofs_on_node,
    HLPOD_MAT*	hlpod_mat);

void hlpod_hr_sys_manusol_set_bc(
    BBFE_DATA* 	fe,
    BBFE_BC*   	bc,
    const int   num_dofs_on_node,
    double     	t,
    double      (*func)(double, double, double, double), // scalar function(x, y, z, t)
    HLPOD_MAT*	hlpod_mat);
