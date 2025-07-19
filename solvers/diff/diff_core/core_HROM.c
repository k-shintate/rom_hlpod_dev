
#include "core_HROM.h"

const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND          = "cond.dat";
static const char* OUTPUT_FILENAME_VTK          = "result_%06d.vtk";
static const char* OUTPUT_FILENAME_ASCII_TEMP   = "temparature_%06d.dat";
static const char* OUTPUT_FILENAME_ASCII_SOURCE = "source_%06d.dat";

void HROM_set_ansvec(
		VALUES*         vals,
	    HLPOD_HR*     	hlpod_hr,
		const int       total_num_nodes)
{
	for(int i = 0; i < total_num_nodes; i++){
		hlpod_hr->HR_T[i] = vals->T[i];
	}
}

void HROM_set_ansvec_para(
		VALUES*         vals,
		HLPOD_DDHR*     hlpod_ddhr,
		const int       total_num_nodes)
{
	for(int i = 0; i < total_num_nodes; i++){
		hlpod_ddhr->HR_T[i] = vals->T[i];
	}
}

/*for Hyper-reduction*/
void HR_output_result_file_vtk(
		BBFE_DATA*     fe,
		VALUES*        vals,
		HLPOD_VALUES*    hlpod_vals,
		HLPOD_HR*      hlpod_hr,		
		const char*    filename,
		const char*    directory,
		double         t)
{
	FILE* fp;
	fp = ROM_BB_write_fopen(fp, filename, directory);

	switch( fe->local_num_nodes ) {
		case 4:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
			break;

		case 8:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);
			break;
	}

	fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
	BB_vtk_write_point_vals_scalar(fp, vals->T, fe->total_num_nodes, "fem-temperature");
	BB_vtk_write_point_vals_scalar(fp, hlpod_vals->sol_vec, fe->total_num_nodes, "pod-temperature");
	BB_vtk_write_point_vals_scalar(fp, hlpod_hr->HR_T, fe->total_num_nodes, "hr-temperature");
	// for manufactured solution
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, vals->T, hlpod_hr->HR_T);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-fem_abs_error");
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, hlpod_vals->sol_vec, hlpod_hr->HR_T);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-pod_abs_error");

	double* source;
	source = BB_std_calloc_1d_double(source, fe->total_num_nodes);
	manusol_set_source(fe, source, t);
	BB_vtk_write_point_vals_scalar(fp, source, fe->total_num_nodes, "source");
	BB_std_free_1d_double(source, fe->total_num_nodes);

	fclose(fp);

}
/********/

/*for Hyper-reduction*/
void HR_output_result_file_vtk_para(
		BBFE_DATA*     fe,
		VALUES*        vals,
		HLPOD_VALUES*    hlpod_vals,
		HLPOD_DDHR*     hlpod_ddhr,		
		const char*    filename,
		const char*    directory,
		double         t)
{
	FILE* fp;
	fp = ROM_BB_write_fopen(fp, filename, directory);

	switch( fe->local_num_nodes ) {
		case 4:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_TETRA);
			break;

		case 8:
			BBFE_sys_write_vtk_shape(fp, fe, TYPE_VTK_HEXAHEDRON);
			break;
	}

	fprintf(fp, "POINT_DATA %d\n", fe->total_num_nodes);
	BB_vtk_write_point_vals_scalar(fp, vals->T, fe->total_num_nodes, "fem-temperature");
	BB_vtk_write_point_vals_scalar(fp, hlpod_vals->sol_vec, fe->total_num_nodes, "pod-temperature");
	BB_vtk_write_point_vals_scalar(fp, hlpod_ddhr->HR_T, fe->total_num_nodes, "hr-temperature");
	// for manufactured solution
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, vals->T, hlpod_ddhr->HR_T);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-fem_abs_error");
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, hlpod_vals->sol_vec, hlpod_ddhr->HR_T);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-pod_abs_error");

	double* source;
	source = BB_std_calloc_1d_double(source, fe->total_num_nodes);
	manusol_set_source(fe, source, t);
	BB_vtk_write_point_vals_scalar(fp, source, fe->total_num_nodes, "source");
	BB_std_free_1d_double(source, fe->total_num_nodes);

	fclose(fp);

}
/********/



/*for Hyper-reduction*/
void HR_output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	const char* filename;
	char fname_vtk[BUFFER_SIZE];
	char fname_tem[BUFFER_SIZE];
	char fname_sou[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, OUTPUT_FILENAME_VTK, file_num);
	snprintf(fname_tem, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_TEMP, file_num);
	snprintf(fname_sou, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_SOURCE, file_num);
/*
	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);

	if(monolis_mpi_get_global_comm_size() == 1 &&sys->rom.hlpod_vals.num_2nd_subdomains == 1){
        HR_output_result_file_vtk(
                &(sys->fe),
                &(sys->vals),
                &(sys->rom.hlpod_vals),
                &(sys->hrom.hlpod_hr),
                filename,
                sys->cond.directory,
                t);
    }
    else{  
        HR_output_result_file_vtk_para(
                &(sys->fe),
                &(sys->vals),
                &(sys->rom.hlpod_vals),
                &(sys->hrom.hlpod_ddhr),
                filename,
                sys->cond.directory,
                t);
    }

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_tem);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			sys->rom.hlpod_vals.sol_vec,
			filename,
			sys->cond.directory);
*/

	/**** for manufactured solution ****/
/*
	double* source;
	source = BB_std_calloc_1d_double(source, sys->fe.total_num_nodes);
	manusol_set_source(&(sys->fe), source, t);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_sou);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			source,
			filename,
			sys->cond.directory);
*/
	double L2_error_fem_hr;
	double L2_error_pod_hr;
	if(monolis_mpi_get_global_comm_size() == 1){
		if(sys->rom.hlpod_vals.num_2nd_subdomains==1){
			L2_error_fem_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->monolis_com),
					t,
					sys->hrom.hlpod_hr.HR_T,
					sys->vals.T);
					
			L2_error_pod_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->monolis_com),
					t,
					sys->rom.hlpod_vals.sol_vec,
					sys->hrom.hlpod_hr.HR_T);
		}
		else{
			L2_error_fem_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->monolis_com),
					t,
					sys->vals.T,
					sys->hrom.hlpod_ddhr.HR_T);
					
			L2_error_pod_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->monolis_com),
					t,
					sys->rom.hlpod_vals.sol_vec,
					sys->hrom.hlpod_ddhr.HR_T);
		}
	}
	else{
		L2_error_fem_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
				&(sys->fe),
				&(sys->basis),
				&(sys->monolis_com),
				t,
				sys->vals.T,
				sys->hrom.hlpod_ddhr.HR_T);
				
		L2_error_pod_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
				&(sys->fe),
				&(sys->basis),
				&(sys->monolis_com),
				t,
				sys->rom.hlpod_vals.sol_vec,
				sys->hrom.hlpod_ddhr.HR_T);
	}

	printf("%s L2 error fem-hrom: %e\n", CODENAME, L2_error_fem_hr);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_fem_hrom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_fem_hr);
		fclose(fp);
	}

	printf("%s L2 error rom-hrom: %e\n", CODENAME, L2_error_pod_hr);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_rom_hrom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_pod_hr);
		fclose(fp);
	}

	//BB_std_free_1d_double(source, sys->fe.total_num_nodes);
	/***********************************/
}
/********/



void HROM_pre_offline(
		FE_SYSTEM* sys,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains)
{
	monolis_initialize(&(sys->monolis_hr0));
	monolis_initialize(&(sys->monolis_hr));

	double t = monolis_get_time_global_sync();
//ROM部分へ移行

	//for arbit dof ddecm
	get_neib_subdomain_id(
		&(sys->monolis_com),
		&(sys->rom.hlpod_mat),
		sys->rom.hlpod_vals.num_2nd_subdomains);		//num_2nd_subdomains

    printf("get_neib_subdomain_id done\n");
    //exit(1);

	set_max_num_modes(
		&(sys->rom.hlpod_vals),
		sys->rom.hlpod_vals.num_modes,
//		sys->rom.hlpod_vals.num_1st_subdomains,
        sys->rom.hlpod_vals.num_2nd_subdomains,
		sys->cond.directory);
/*
	//level2領域の最大基底本数の共有
    get_neib_max_num_modes_pad(
		&(sys->mono_com_rom_solv),
        &(sys->rom.hlpod_vals),
		&(sys->rom.hlpod_mat),
        1 + sys->mono_com_rom_solv.recv_n_neib,
		sys->rom.hlpod_vals.num_modes_max);
*/
    double t1 = monolis_get_time_global_sync();

    get_neib_num_modes_pad(
        &(sys->mono_com_rom),
        &(sys->rom.hlpod_vals),
        &(sys->rom.hlpod_mat),
        1 + sys->mono_com_rom_solv.recv_n_neib,
        sys->rom.hlpod_vals.num_modes);

    get_neib_coordinates_pre(
        &(sys->rom.hlpod_vals),
        &(sys->rom.hlpod_mat),
        1 + sys->mono_com_rom_solv.recv_n_neib,
        sys->rom.hlpod_vals.num_modes_max);

	//level2領域の最大基底本数の共有
    get_neib_max_num_modes_pad(
		&(sys->mono_com_rom),
        &(sys->rom.hlpod_vals),
		&(sys->rom.hlpod_mat),
        1 + sys->mono_com_rom_solv.recv_n_neib,
		sys->rom.hlpod_vals.num_modes_max);

	get_meta_neib(
		&(sys->mono_com_rom_solv),
		&(sys->rom.hlpod_meta),
		sys->cond.directory);

	ddhr_lb_set_neib(
		&(sys->mono_com_rom_solv),
		//&(sys->fe),
		&(sys->rom.hlpod_mat),
		&(sys->hrom.hlpod_ddhr),
		&(sys->rom.hlpod_meta),
		num_2nd_subdomains,
		sys->rom.hlpod_vals.num_snapshot,
		sys->cond.directory);

}


void HROM_pre_offline2(
		FE_SYSTEM* sys,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains)
{

	ddhr_lb_write_selected_elements_para_1line(
		&(sys->mono_com_rom_solv),
		&(sys->fe),
		&(sys->bc),
        &(sys->rom.hlpod_vals),
		&(sys->hrom.hlpod_ddhr),
		&(sys->rom.hlpod_mat),
		&(sys->rom.hlpod_meta),
		sys->fe.total_num_elems,
		sys->rom.hlpod_vals.num_snapshot,
		sys->rom.hlpod_vals.num_modes,
		num_2nd_subdomains,
		10000,
		1.0e-6,
		sys->cond.directory);

    ddhr_lb_get_selected_elements_internal_overlap(
            &(sys->hrom.hlpod_ddhr),
            sys->cond.directory);

    double t_tmp = monolis_get_time_global_sync();

}


void HROM_pre_online(
		FE_SYSTEM* sys,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains)
{
	monolis_initialize(&(sys->monolis_hr0));
	monolis_initialize(&(sys->monolis_hr));

	double t = monolis_get_time_global_sync();
//ROM部分へ移行


    //printf("monolis_mpi_get_global_comm_size() = %d\n", monolis_mpi_get_global_comm_size());
    //exit(1);

    //for arbit dof ddecm
    get_neib_subdomain_id_2nddd(
        &(sys->monolis_com),
        &(sys->rom.hlpod_mat),
        sys->rom.hlpod_vals.num_2nd_subdomains);              //num_2nd_subdomains

	ddhr_lb_read_selected_elements_para(
		num_2nd_subdomains,
		sys->cond.directory);

	ddhr_lb_get_selected_elements_para_add(
		&(sys->hrom.hlpod_ddhr),
		//sys->rom.hlpod_vals.num_1st_subdomains,
        monolis_mpi_get_global_comm_size(),
		sys->cond.directory);

	/*
	ddhr_lb_set_reduced_mat_para(
		&(sys->monolis_hr0),
		&(sys->fe),
		&(sys->basis),
		&(sys->bc),
		&(sys->rom.hlpod_mat),
		&(sys->hrom.hlpod_ddhr),
		sys->rom.hlpod_vals.num_modes,
		num_2nd_subdomains,
		sys->vals.dt);

	ddhr_lb_write_reduced_mat(
		sys->hrom.hlpod_ddhr.reduced_mat,
		sys->rom.hlpod_mat.n_neib_vec,	
		sys->rom.hlpod_mat.n_neib_vec,
		sys->cond.directory);
	
	t = monolis_get_time_global_sync();
//	exit(1);
*/
/*
	ddhr_lb_read_reduced_mat(
		sys->hrom.hlpod_ddhr.reduced_mat,
		sys->rom.hlpod_mat.n_neib_vec,	
		sys->rom.hlpod_mat.n_neib_vec,
		sys->cond.directory);
*/


	lpod_hdd_lb_set_hr_podbasis(
		&(sys->monolis_com),
		&(sys->rom.hlpod_vals),
		&(sys->rom.hlpod_mat),
        sys->fe.total_num_nodes);
//        sys->rom.hlpod_vals.num_modes);

	ddhr_lb_set_reduced_mat_para_save_memory(
		&(sys->monolis_hr0),
		&(sys->fe),
		&(sys->basis),
		&(sys->bc),
		&(sys->rom.hlpod_vals),
		&(sys->rom.hlpod_mat),
		&(sys->hrom.hlpod_ddhr),
		sys->rom.hlpod_vals.num_modes,
		num_2nd_subdomains,
		sys->vals.dt);

/*
	hlpod_set_comm_nzpattern_bcsr(
		&(sys->monolis_hr0),
		&(sys->mono_com_rom_solv),
		&(sys->rom.hlpod_mat),
		sys->rom.hlpod_vals.num_modes,
		sys->cond.directory);

    ROM_std_hlpod_set_nonzero_pattern_bcsr(
            &(sys->monolis_hr0),
            &(sys->rom.hlpod_mat),
            &(sys->rom.hlpod_meta),
            sys->rom.hlpod_vals.num_modes_pre
*/

    monolis_get_nonzero_pattern_by_nodal_graph_V_R(
        &(sys->monolis_hr0),
        sys->rom.hlpod_meta.num_meta_nodes,
        sys->rom.hlpod_meta.n_dof_list,
        sys->rom.hlpod_meta.index,
        sys->rom.hlpod_meta.item);

	ddhr_hlpod_calc_block_mat_bcsr_pad(
		&(sys->monolis_hr0),
		&(sys->mono_com_rom_solv),
		//&(sys->lpod_com),
        &(sys->rom.hlpod_vals),
		&(sys->rom.hlpod_mat),
//		&(sys->lpod_prm),
		&(sys->hrom.hlpod_ddhr),
		&(sys->rom.hlpod_meta),
		sys->rom.hlpod_vals.num_modes,
		sys->rom.hlpod_vals.num_2nd_subdomains,
		sys->cond.directory);
	
}


void HROM_hierarchical_parallel(
    FE_SYSTEM sys,
    const int step_HR,
    const int step_POD,
    const double t)
{
    if((step_HR-step_POD) == 1){
        HROM_set_ansvec_para(
            &(sys.vals),
            &(sys.hrom.hlpod_ddhr),
            sys.fe.total_num_nodes);
    }

    double t1 = monolis_get_time();

    //monolis_copy_mat_value_R(&(sys.monolis_hr0), &(sys.monolis_hr));
	monolis_copy_mat_value_matrix_R(&(sys.monolis_hr0), &(sys.monolis_hr));
	monolis_clear_mat_value_rhs_R(&(sys.monolis_hr));

	double set_bc1 = monolis_get_time();
	hlpod_hr_sys_manusol_set_bc(
        &(sys.fe),
        &(sys.bc),
        BLOCK_SIZE,
        t,
        manusol_get_sol,
        &(sys.rom.hlpod_mat));
	double set_bc2 = monolis_get_time();
	
    ddhr_lb_set_reduced_vec_para(
        &(sys.monolis_hr),
        &(sys.fe),
        &(sys.basis),
        &(sys.rom.hlpod_vals),
        &(sys.hrom.hlpod_ddhr),
        &(sys.rom.hlpod_mat),
        sys.rom.hlpod_vals.num_modes,
        sys.rom.hlpod_vals.num_2nd_subdomains,
        sys.vals.dt,
        t);
    
    ddhr_lb_set_D_bc_para(
        &(sys.monolis_hr),
        &(sys.fe),
        &(sys.basis),
        &(sys.bc),
        &(sys.rom.hlpod_mat),
        &(sys.hrom.hlpod_ddhr),
        sys.rom.hlpod_vals.num_modes,
        sys.rom.hlpod_vals.num_2nd_subdomains,
        sys.vals.dt);

    ddhr_to_monollis_rhs_para_pad(
        &(sys.monolis_hr),
        &(sys.hrom.hlpod_ddhr),
        &(sys.rom.hlpod_mat),
		sys.rom.hlpod_vals.num_2nd_subdomains,
        sys.rom.hlpod_vals.num_modes);

    double t2 = monolis_get_time();

    if(monolis_mpi_get_global_my_rank() == 0){
        FILE* fp;
        fp = ROM_BB_write_add_fopen(fp, "calctime/hr_time_calc_mat.txt", sys.cond.directory);
        fprintf(fp, "%e %e\n", t, t2-t1);
        fclose(fp);
    }


    BBFE_sys_monowrap_solve(
        &(sys.monolis_hr),
        &(sys.mono_com_rom_solv),
        sys.rom.hlpod_mat.mode_coef,
        MONOLIS_ITER_CG,
        MONOLIS_PREC_DIAG,
        sys.vals.mat_max_iter,
        sys.vals.mat_epsilon);

/*
	BBFE_sys_monowrap_solve_V(
		&(sys.monolis_hr),
		&(sys.mono_com_rom_solv),
		sys.rom.hlpod_mat.mode_coef,
		MONOLIS_ITER_CG,
		MONOLIS_PREC_DIAG,
		sys.vals.mat_max_iter,
		sys.vals.mat_epsilon,
		sys.rom.hlpod_mat.n_dof_list);
*/
    t1 = monolis_get_time();

    lpod_pad_calc_block_solution_local_para_pad(
        &(sys.monolis_com),
        &(sys.fe),
        &(sys.hrom.hlpod_ddhr),
        &(sys.rom.hlpod_mat),
        &(sys.bc),
        sys.rom.hlpod_vals.num_2nd_subdomains,
		sys.rom.hlpod_vals.num_modes);
//        &(sys.rom.hlpod_mat));

    t2 = monolis_get_time();

    if(monolis_mpi_get_global_my_rank() == 0){
        FILE* fp;
        fp = ROM_BB_write_add_fopen(fp, "calctime/hr_time_calc_sol.txt", sys.cond.directory);
        fprintf(fp, "%e %e\n", t, t2-t1);
        fclose(fp);
    }

	output_hr_monolis_solver_prm(&(sys.monolis_hr), sys.cond.directory, t);

}
