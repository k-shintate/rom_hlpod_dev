
#include "core_HROM.h"

const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND          = "cond.dat";
static const char* OUTPUT_FILENAME_VTK          = "result_%06d.vtk";
static const char* OUTPUT_FILENAME_ASCII_TEMP   = "temparature_%06d.dat";
static const char* OUTPUT_FILENAME_ASCII_SOURCE = "source_%06d.dat";

/*for Hyper-reduction*/
void HR_output_result_file_vtk(
		BBFE_DATA*     fe,
		VALUES*        vals,
        HR_VALUES*      hr_vals,
		HLPOD_VALUES*    hlpod_vals,
		HLPOD_HR*      hlpod_hr,		
		const char*    filename,
		const char*    directory,
		double         t)
{
	FILE* fp;
	fp = ROM_BB_write_fopen(fp, filename, directory);
/*
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
	BB_vtk_write_point_vals_scalar(fp, hr_vals->sol_vec, fe->total_num_nodes, "hr-temperature");
	// for manufactured solution
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, vals->T, hr_vals->sol_vec);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-fem_abs_error");
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, hlpod_vals->sol_vec, hr_vals->sol_vec);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-pod_abs_error");

	double* source;
	source = BB_std_calloc_1d_double(source, fe->total_num_nodes);
	manusol_set_source(fe, source, t);
	BB_vtk_write_point_vals_scalar(fp, source, fe->total_num_nodes, "source");
	BB_std_free_1d_double(source, fe->total_num_nodes);

	fclose(fp);
*/
}
/********/

/*for Hyper-reduction*/
void HR_output_result_file_vtk_para(
		BBFE_DATA*     fe,
		VALUES*        vals,
        HR_VALUES*      hr_vals,
		HLPOD_VALUES*    hlpod_vals,
		HLPOD_DDHR*     hlpod_ddhr,		
		const char*    filename,
		const char*    directory,
		double         t)
{
	FILE* fp;
	fp = ROM_BB_write_fopen(fp, filename, directory);
/*
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
	BB_vtk_write_point_vals_scalar(fp, hr_vals->sol_vec, fe->total_num_nodes, "hr-temperature");
	// for manufactured solution
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, vals->T, hr_vals->sol_vec);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-fem_abs_error");
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, hlpod_vals->sol_vec, hr_vals->sol_vec);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-pod_abs_error");

	double* source;
	source = BB_std_calloc_1d_double(source, fe->total_num_nodes);
	manusol_set_source(fe, source, t);
	BB_vtk_write_point_vals_scalar(fp, source, fe->total_num_nodes, "source");
	BB_std_free_1d_double(source, fe->total_num_nodes);

	fclose(fp);
*/
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

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);
/*
	if(monolis_mpi_get_global_comm_size() == 1 &&sys->rom.hlpod_vals.num_2nd_subdomains == 1){
        HR_output_result_file_vtk(
                &(sys->fe),
                &(sys->vals),
                &(sys->hrom.hr_vals),
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
                &(sys->hrom.hr_vals),
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

	double L2_error_fem_hr;
	double L2_error_pod_hr;
	if(monolis_mpi_get_global_comm_size() == 1){
		if(sys->rom.hlpod_vals.num_2nd_subdomains==1){
			L2_error_fem_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->mono_com),
					t,
					sys->hrom.hr_vals.sol_vec,
					sys->vals.T);
					
			L2_error_pod_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->mono_com),
					t,
					sys->rom.hlpod_vals.sol_vec,
					sys->hrom.hr_vals.sol_vec);
		}
		else{
			L2_error_fem_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->mono_com),
					t,
					sys->vals.T,
					sys->hrom.hr_vals.sol_vec);
					
			L2_error_pod_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->mono_com),
					t,
					sys->rom.hlpod_vals.sol_vec,
					sys->hrom.hr_vals.sol_vec);
		}
	}
	else{
		L2_error_fem_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
				&(sys->fe),
				&(sys->basis),
				&(sys->mono_com),
				t,
				sys->vals.T,
				sys->hrom.hr_vals.sol_vec);
				
		L2_error_pod_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
				&(sys->fe),
				&(sys->basis),
				&(sys->mono_com),
				t,
				sys->rom.hlpod_vals.sol_vec,
				sys->hrom.hr_vals.sol_vec);
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
*/

}
/********/



void HROM_pre_offline(
		FE_SYSTEM* sys,
        ROM*            rom,
        HROM*           hrom,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains)
{
	monolis_initialize(&(sys->monolis_hr0));
	monolis_initialize(&(sys->monolis_hr));

	double t = monolis_get_time_global_sync();
	//for arbit dof ddecm
	get_neib_subdomain_id(
		&(sys->mono_com),
		&(rom->hlpod_mat),
		rom->hlpod_vals.num_2nd_subdomains);		//num_2nd_subdomains

    printf("get_neib_subdomain_id done\n");

    double t1 = monolis_get_time_global_sync();

    get_neib_num_modes_pad(
        &(sys->mono_com_rom),
        &(rom->hlpod_vals),
        &(rom->hlpod_mat),
        1 + sys->mono_com.recv_n_neib,
        rom->hlpod_vals.num_modes);


    get_neib_coordinates_pre(
        &(rom->hlpod_vals),
        &(rom->hlpod_mat),
        1 + sys->mono_com.recv_n_neib,
        rom->hlpod_vals.num_modes_max);

    printf("%d %d\n", rom->hlpod_vals.num_modes_max, rom->hlpod_vals.num_modes_pre);

	//level2領域の最大基底本数の共有
    get_neib_max_num_modes_pad(
		&(sys->mono_com_rom),
        &(rom->hlpod_vals),
		&(rom->hlpod_mat),
        1 + sys->mono_com.recv_n_neib,
		rom->hlpod_vals.num_modes_max);

	get_meta_neib(
		&(sys->mono_com_rom_solv),
		&(rom->hlpod_meta),
		sys->cond.directory);

	ddhr_lb_set_neib(
		&(sys->mono_com_rom_solv),
		&(rom->hlpod_mat),
		&(hrom->hlpod_ddhr),
		&(rom->hlpod_meta),
		num_2nd_subdomains,
		rom->hlpod_vals.num_snapshot,
		sys->cond.directory);

}


void HROM_pre_offline2(
		FE_SYSTEM* sys,
        ROM*            rom,
        HROM*           hrom,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains)
{
	
        ddhr_lb_write_selected_elements_para_1line(
            &(sys->mono_com_rom_solv),
            &(sys->fe),
            &(sys->bc),
            &(rom->hlpod_vals),
            &(hrom->hlpod_ddhr),
            &(rom->hlpod_mat),
            &(rom->hlpod_meta),
            sys->fe.total_num_elems,
            rom->hlpod_vals.num_snapshot,
            rom->hlpod_vals.num_modes_pre,
            num_2nd_subdomains,
            10000,
            1.0e-8,
            sys->cond.directory);

        ddhr_lb_get_selected_elements_internal_overlap(
                &(hrom->hlpod_ddhr),
                sys->cond.directory);

        double t_tmp = monolis_get_time_global_sync();


}


void HROM_pre_online(
		FE_SYSTEM* sys,
        ROM*            rom,
        HROM*           hrom,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains)
{
	monolis_initialize(&(sys->monolis_hr));

	double t = monolis_get_time_global_sync();

        get_neib_subdomain_id_2nddd(
            &(sys->mono_com),
            &(rom->hlpod_mat),
            rom->hlpod_vals.num_2nd_subdomains);

        ddhr_lb_read_selected_elements_para(
            num_2nd_subdomains,
            sys->cond.directory);

        ddhr_lb_get_selected_elements_para_add(
            &(hrom->hlpod_ddhr),
            monolis_mpi_get_global_comm_size(),
            sys->cond.directory);

        lpod_hdd_lb_set_hr_podbasis(
            &(sys->mono_com),
            &(rom->hlpod_vals),
            &(rom->hlpod_mat),
            sys->fe.total_num_nodes);

        ddhr_lb_set_reduced_mat_para_save_memory(
            &(sys->monolis_hr0),
            &(sys->fe),
            &(sys->vals),
            &(sys->basis),
            &(sys->bc),
            &(rom->hlpod_vals),
            &(rom->hlpod_mat),
            &(hrom->hlpod_ddhr),
            rom->hlpod_vals.num_modes_pre,
            num_2nd_subdomains,
            sys->vals.dt);

        monolis_get_nonzero_pattern_by_nodal_graph_V_R(
            &(sys->monolis_hr0),
            rom->hlpod_meta.num_meta_nodes,
            rom->hlpod_meta.n_dof_list,
            rom->hlpod_meta.index,
            rom->hlpod_meta.item);

        ddhr_hlpod_calc_block_mat_bcsr_pad(
            &(sys->monolis_hr0),
            &(sys->mono_com_rom_solv),
            &(rom->hlpod_vals),
            &(rom->hlpod_mat),
            &(hrom->hlpod_ddhr),
            &(rom->hlpod_meta),
            rom->hlpod_vals.num_modes_pre,
            rom->hlpod_vals.num_2nd_subdomains,
            sys->cond.directory);
    
	
}

void HROM_hierarchical_parallel(
    FE_SYSTEM   sys,
    ROM*        rom,
    HROM*       hrom,
    const int   step_HR,
    const int   step_POD,
    const double t)
{
    double t1 = monolis_get_time();

	monolis_copy_mat_value_matrix_R(&(sys.monolis_hr0), &(sys.monolis_hr));
	monolis_clear_mat_value_rhs_R(&(sys.monolis_hr));

	double set_bc1 = monolis_get_time();
	/*
    hlpod_hr_sys_manusol_set_bc(
        &(sys.fe),
        &(sys.bc),
        BLOCK_SIZE,
        t,
        manusol_get_sol,
        &(rom->hlpod_mat));
    BBFE_sys_monowrap_set_Dirichlet_bc(
        &(sys->monolis),
        sys->fe.total_num_nodes,
        4,
        &(sys->bc),
        sys->monolis.mat.R.B);
        */
	double set_bc2 = monolis_get_time();
	
    ddhr_lb_set_reduced_vec_para(
        &(sys.monolis_hr),
        &(sys.fe),
        &(sys.vals),
        &(sys.basis),
        &(hrom->hr_vals),
        &(rom->hlpod_vals),
        &(hrom->hlpod_ddhr),
        &(rom->hlpod_mat),
        rom->hlpod_vals.num_modes_pre,
        rom->hlpod_vals.num_2nd_subdomains,
        sys.vals.dt,
        t);
    
    ddhr_lb_set_D_bc_para(
        &(sys.monolis_hr),
        &(sys.fe),
        &(sys.vals),
        &(sys.basis),
        &(sys.bc),
        &(rom->hlpod_mat),
        &(hrom->hlpod_ddhr),
        rom->hlpod_vals.num_modes_pre,
        rom->hlpod_vals.num_2nd_subdomains,
        sys.vals.dt);

    ddhr_to_monollis_rhs_para_pad(
        &(sys.monolis_hr),
        &(hrom->hlpod_ddhr),
        &(rom->hlpod_mat),
		rom->hlpod_vals.num_2nd_subdomains,
        rom->hlpod_vals.num_modes_pre);

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
        rom->hlpod_mat.mode_coef,
        MONOLIS_ITER_BICGSTAB,
        MONOLIS_PREC_DIAG,
        sys.vals.mat_max_iter,
        sys.vals.mat_epsilon);

    t1 = monolis_get_time();

    lpod_pad_calc_block_solution_local_para_pad(
        &(sys.mono_com),
        &(sys.fe),
        &(hrom->hr_vals),
        &(rom->hlpod_mat),
        rom->hlpod_vals.num_2nd_subdomains,
		4);

	ROM_sys_hlpod_fe_add_Dbc(
        hrom->hr_vals.sol_vec,
		&(sys.bc),
		sys.fe.total_num_nodes,
		1);

	monolis_mpi_update_R(&(sys.mono_com), sys.fe.total_num_nodes, 4, sys.rom_sups.hlpod_vals.sol_vec);

    BBFE_fluid_sups_renew_velocity(
        sys.vals_rom.v, 
        sys.rom_sups.hlpod_vals.sol_vec,
        sys.fe.total_num_nodes);

    BBFE_fluid_sups_renew_pressure(
        sys.vals_rom.p, 
        sys.rom_sups.hlpod_vals.sol_vec,
        sys.fe.total_num_nodes);
        
    //output_files(sys, step, t);

    t2 = monolis_get_time();

    if(monolis_mpi_get_global_my_rank() == 0){
        FILE* fp;
        fp = ROM_BB_write_add_fopen(fp, "calctime/hr_time_calc_sol.txt", sys.cond.directory);
        fprintf(fp, "%e %e\n", t, t2-t1);
        fclose(fp);
    }

	output_hr_monolis_solver_prm(&(sys.monolis_hr), sys.cond.directory, t);

}


void HROM_pre(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom)
{
    if(monolis_mpi_get_global_comm_size() == 1){
    }
    else{
        HROM_pre_offline(sys, rom, hrom, rom->hlpod_vals.num_modes_pre, rom->hlpod_vals.num_snapshot, rom->hlpod_vals.num_2nd_subdomains);
    }
}


void HROM_memory_allocation(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom)
{
        ddhr_lb_set_element_para2(
                &(sys->fe),
                &(hrom->hlpod_ddhr),
                rom->hlpod_vals.num_2nd_subdomains,
                sys->cond.directory);

        //基底本数の分布が決定されてからメモリ割り当て
        ddhr_memory_allocation_para(
                &(rom->hlpod_vals),
                &(hrom->hlpod_ddhr),
                &(rom->hlpod_mat),
                sys->fe.total_num_nodes,
                sys->fe.total_num_elems,
                rom->hlpod_vals.num_snapshot,
                rom->hlpod_vals.num_modes_pre,
                rom->hlpod_vals.num_2nd_subdomains);

}


void HROM_set_matvec(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom,
        int step,
        double t)
{
        get_neib_coordinates_pad(
                &(sys->mono_com_rom),
                &(rom->hlpod_vals),
                &(rom->hlpod_mat),
                1 + sys->mono_com_rom_solv.recv_n_neib,
                rom->hlpod_vals.num_modes_max,
                rom->hlpod_vals.num_2nd_subdomains,
                rom->hlpod_vals.num_modes_pre);

        ddhr_set_matvec_RH_for_NNLS_para_only_residuals(
                &(sys->fe),
                &(sys->vals),
                &(sys->basis),
                &(rom->hlpod_mat),
                &(rom->hlpod_vals),
                &(hrom->hlpod_ddhr),
                rom->hlpod_vals.num_2nd_subdomains,
                step -1 ,   //index 0 start
                rom->hlpod_vals.num_snapshot,
                rom->hlpod_vals.num_modes_pre,
                sys->vals.dt,
                t);

        ddhr_set_matvec_residuals_for_NNLS_para_only_residuals(
                &(sys->fe),
                &(sys->vals),
                &(sys->basis),
                &(sys->bc),
                &(rom->hlpod_mat),
                &(rom->hlpod_vals),
                &(hrom->hlpod_ddhr),
                rom->hlpod_vals.num_2nd_subdomains,
                step -1 ,   //index 0 start
                rom->hlpod_vals.num_snapshot,
                1 + sys->mono_com.recv_n_neib,
                sys->vals.dt,
                t);

}


void HROM_pre_offline2(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom)
{
    if(monolis_mpi_get_global_comm_size() == 1){
        HROM_pre_offline2(sys, rom, hrom, rom->hlpod_vals.num_modes_pre, rom->hlpod_vals.num_snapshot, rom->hlpod_vals.num_2nd_subdomains);
    }
    else{
        HROM_pre_offline2(sys, rom, hrom, rom->hlpod_vals.num_modes_pre, rom->hlpod_vals.num_snapshot, rom->hlpod_vals.num_2nd_subdomains);
    }
}


void HROM_memory_allocation_online(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom)
{
    ddhr_lb_set_element_para2(
        &(sys->fe),
        &(hrom->hlpod_ddhr),
        rom->hlpod_vals.num_2nd_subdomains,
        sys->cond.directory);

    ddhr_memory_allocation_para_online(
        &(rom->hlpod_vals),
        &(hrom->hlpod_ddhr),
        &(rom->hlpod_mat),
        sys->fe.total_num_nodes);

}


void HROM_pre_online(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom)
{
    if(monolis_mpi_get_global_comm_size() == 1){		
		HROM_pre_online(sys, rom, hrom, rom->hlpod_vals.num_modes_pre, rom->hlpod_vals.num_snapshot, rom->hlpod_vals.num_2nd_subdomains);
	}
	else{
		HROM_pre_online(sys, rom, hrom, rom->hlpod_vals.num_modes_pre, rom->hlpod_vals.num_snapshot, rom->hlpod_vals.num_2nd_subdomains);
	}
}