#pragma once

#include "monolis_nnls_c.h"
#include "hrom_dataset.h"
#include "hlpod_core_fe.h"
#include "write_std.h"
#include "write_BB.h"


void ddhr_memory_allocation_para_online(
            HLPOD_DDHR*     hlpod_ddhr,
	    HLPOD_MAT*    hlpod_mat,
        const int       total_num_nodes);


void ddhr_memory_allocation_para(
	    HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_MAT*     	hlpod_mat,
        const int       total_num_nodes,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
		const int		num_subdomains);
/*
void HROM_set_ansvec_para(
		VALUES*        vals,
		HLPOD_DDHR*     hlpod_ddhr,
		const int       total_num_nodes);
*/

void ddhr_set_element_para(
        HLPOD_DDHR*     hlpod_ddhr,
		const int 		num_subdomains,
		const char*     directory);

void ddhr_get_selected_elements_para(
        BBFE_DATA*     	fe,
        BBFE_BC*     	bc,
	    HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_MAT*    hlpod_mat,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
		const int 		num_subdomains,
        const int       max_iter, //NNLS
        const double    tol,      //NNLS
		const char*		directory);

void ddhr_get_selected_elements_para_add(
	    HLPOD_DDHR*     hlpod_ddhr,
		const int		num_parallel_subdomains,
		const char*		directory);

/*for visualization*/
void ddhr_set_selected_elems_para(
		BBFE_DATA*     	fe,
        HLPOD_DDHR*     hlpod_ddhr,
        const int		total_num_nodes,
		const int		num_subdomains,
		const char*     directory);

void ddhr_monolis_set_matrix_para(
        MONOLIS*     	monolis,
        HLPOD_DDHR*     hlpod_ddhr,
        HLPOD_MAT*    hlpod_mat,
        LPOD_COM*    	lpod_com,
        const int 		max_num_basis,
        const int		num_subdomains);


void ddhr_to_monollis_rhs_para(
		MONOLIS*		monolis,
	//	MONOLIS_COM*	monolis_com,
		HLPOD_DDHR*     hlpod_ddrh,
		HLPOD_MAT*    hlpod_mat,
		const int 		k);

void lpod_pad_calc_block_solution_local_para(
		MONOLIS_COM*	monolis_com,
		BBFE_DATA* 		fe,
		HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_MAT*    hlpod_mat,
		BBFE_BC*      	bc,
		const int		num_2nddd);
//		LPOD_PRM*		lpod_prm);

void ddhr_to_monollis_rhs_para_pad(
		MONOLIS*		monolis,
		HLPOD_DDHR*     hlpod_ddrh,
		HLPOD_MAT*    hlpod_mat,
		const int		num_2nddd,
		const int		max_num_bases);

void lpod_pad_calc_block_solution_local_para_pad(
		MONOLIS_COM*	monolis_com,
		BBFE_DATA* 		fe,
		HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_MAT*    hlpod_mat,
		BBFE_BC*      	bc,
		const int		num_2nddd,
		const int		max_num_bases)
//		LPOD_PRM*		lpod_prm);
