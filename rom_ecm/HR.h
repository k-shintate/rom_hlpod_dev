#pragma once

#include "pod_dataset.h"
//#include "diff_dataset.h"
#include "hlpod_write_fe.h"

#include "monolis_nnls_c.h"

void hr_memory_allocation(
        const int       total_num_nodes,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
        HLPOD_HR*       hlpod_hr);
/*
void HROM_set_ansvec(
		VALUES*        vals,
	    HLPOD_HR*     	hlpod_hr,
		const int       total_num_nodes);
*/

void hr_get_selected_elements(
        BBFE_DATA*     	fe,
        BBFE_BC*     	bc,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
        const int       max_iter, //NNLS
        const double    tol,      //NNLS
        HLPOD_HR*       hlpod_hr,
		const char*		directory);

void hr_calc_solution(
        BBFE_DATA* 		fe,
        POD_MATRIX*     pod_mat,
        HLPOD_HR*       hlpod_hr,
        BBFE_BC*     	bc,
        int 			num_base,
        LPOD_PRM*		lpod_prm);


void hr_monolis_set_matrix(
        MONOLIS*     	monolis,
        HLPOD_HR*      hlpod_hr,
        const int 		num_basis);

void hr_to_monollis_rhs(
        MONOLIS*		monolis,
        HLPOD_HR*       hlpod_rh,
        const int 		k);

void hr_set_selected_elems(
		BBFE_DATA*     	fe,
        HLPOD_HR*       hlpod_hr,
        const int		total_num_nodes,
		const char*     directory);