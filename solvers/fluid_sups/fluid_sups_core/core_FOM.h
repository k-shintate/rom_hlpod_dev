#pragma once

#include "fluid_core.h"
#include "hlpod_dataset.h"
#include "fluid_sups_dataset.h"

void output_cavity_center_vx(
    VALUES*        vals,
    const char*    method,
    const char*    directory);

void initialize_velocity_pressure_karman_vortex(
	double** v,
	double* p,
	const int total_num_nodes);

void BBFE_fluid_sups_read_Dirichlet_bc_karman_vortex(
    BBFE_BC*     bc,
    const char*  filename,
    const char*  directory,
    const int    total_num_nodes,
    const int    block_size);

void output_result_file_karman_vortex(
    BBFE_DATA*     fe,
    VALUES*        vals,
    double         t,
    const char*    directory);

void output_result_file_karman_vortex_pressure(
    BBFE_DATA*     fe,
    VALUES*        vals,
    double         t,
    const char*    directory);

void output_result_file_karman_vortex_pressure_inf(
    BBFE_DATA*     fe,
    VALUES*        vals,
    double         t,
    const char*    directory);

void memory_allocation_nodal_values(
    VALUES*         vals,
    const int       total_num_nodes);

void assign_default_values(
    VALUES*         vals);

void print_all_values(
    VALUES*         vals);

void read_calc_conditions(
    VALUES*         vals,
    const char*     directory);

void output_result_file_vtk(
    BBFE_DATA*     fe,
    VALUES*        vals,
    const char*    filename,
    const char*    directory,
    double         t);

void output_files(
    FE_SYSTEM*  sys,
    int         file_num,
    double      t);

void set_element_mat_pred(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_pred(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_corr(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec_ppe(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void BBFE_fluid_sups_read_Dirichlet_bc(
    BBFE_BC*     bc,
    const char*  filename,
    const char*  directory,
    const int    total_num_nodes,
    const int    block_size);

void set_element_mat(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

void set_element_vec(
    MONOLIS*     monolis,
    BBFE_DATA*   fe,
    BBFE_BASIS*  basis,
    VALUES*      vals);

//右辺ベクトルのアップデートと解ベクトルの求解
void solver_fom(
    FE_SYSTEM   sys,
    double      t,
    const int   step);

void solver_fom_collect_snapmat(
    FE_SYSTEM   sys,
    double      t,
    const int   step);