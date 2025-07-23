#pragma once

#include "BBFE/sys/FE_dataset.h"

#include "fluid_sups_dataset.h"

#include "set_reduced_matvec.h"

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
    HLPOD_MAT*      hlpod_mat,
    HLPOD_DDHR*     hlpod_ddhr,
    const int		num_modes,
    const int 		num_subdomains,
    const double    dt);

void ddhr_lb_set_reduced_vec_para(
    MONOLIS*     	monolis,
    BBFE_DATA*     	fe,
    BBFE_BASIS*	 	basis,
    HR_VALUES*      hr_vals,
    HLPOD_VALUES*   hlpod_vals,
    HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_MAT*      hlpod_mat,
    const int		num_modes,
    const int 		num_subdomains,
    const double    dt,
    double       	t);
