#pragma once

#include "monolis_nnls_c.h"
#include "hrom_dataset.h"
#include "hlpod_core_fe.h"
#include "write_std.h"
#include "write_BB.h"

void memory_allocation_hr_sol_vec(
	    HR_VALUES*      hr_vals,
        const int       total_num_nodes,
        const int       dof);

void ddhr_memory_allocation_para_online(
        HLPOD_VALUES*   hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
	    HLPOD_MAT*      hlpod_mat,
        const int       total_num_nodes);

void ddhr_memory_allocation_para(
        HLPOD_VALUES*   hlpod_vals,
	    HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_MAT*     	hlpod_mat,
        const int       total_num_nodes,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
		const int		num_subdomains);

/*for visualization*/
void ddhr_set_selected_elems_para(
		BBFE_DATA*     	fe,
        HLPOD_DDHR*     hlpod_ddhr,
        const int		total_num_nodes,
		const int		num_subdomains,
		const char*     directory);

void ddhr_to_monollis_rhs_para(
		MONOLIS*		monolis,
		HLPOD_DDHR*     hlpod_ddrh,
		HLPOD_MAT*      hlpod_mat,
		const int 		k);

void ddhr_to_monollis_rhs_para_pad(
		MONOLIS*		monolis,
		HLPOD_DDHR*     hlpod_ddrh,
		HLPOD_MAT*      hlpod_mat,
		const int		num_2nddd,
		const int		max_num_bases);

void lpod_pad_calc_block_solution_local_para_pad(
	MONOLIS_COM*	monolis_com,
	BBFE_DATA* 		fe,
    HR_VALUES*      hr_vals,
	HLPOD_MAT*      hlpod_mat,
	const int		num_2nddd,
    const int 		dof);