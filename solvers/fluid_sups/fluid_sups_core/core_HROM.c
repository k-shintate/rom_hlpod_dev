
#include "core_HROM.h"

const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND          = "cond.dat";
static const char* OUTPUT_FILENAME_VTK          = "result_%06d.vtk";
static const char* OUTPUT_FILENAME_ASCII_TEMP   = "temparature_%06d.dat";
static const char* OUTPUT_FILENAME_ASCII_SOURCE = "source_%06d.dat";

static const char* ROM_OUTPUT_FILENAME_VTK          = "rom_result_%06d.vtk";


void HROM_output_vtk_shape(
		BBFE_DATA*     fe,
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

	fclose(fp);
}

void HROM_add_output_vtk_pressure(
		BBFE_DATA*     fe,
		double*        fem_pressure,
		double*    	   rom_pressure,
        double*    	   hrom_pressure,
		const char*    filename,
		const char*    directory,
		double         t)
{
	FILE* fp;
	fp = ROM_BB_write_add_fopen(fp, filename, directory);

	/*for pressure*/
	BB_vtk_write_point_vals_scalar(fp, rom_pressure, fe->total_num_nodes, "ROM-Pressure");

    BB_vtk_write_point_vals_scalar(fp, hrom_pressure, fe->total_num_nodes, "HROM-Pressure");

	double* error_p;
	error_p = BB_std_calloc_1d_double(error_p, fe->total_num_nodes);
	BBFE_manusol_calc_nodal_error_scalar(fe, error_p, fem_pressure, rom_pressure);
	BB_vtk_write_point_vals_scalar(fp, error_p, fe->total_num_nodes, "abs_error-FEM_ROM-Pressure");

	BBFE_manusol_calc_nodal_error_scalar(fe, error_p, fem_pressure, hrom_pressure);
	BB_vtk_write_point_vals_scalar(fp, error_p, fe->total_num_nodes, "abs_error-FEM_HROM-Pressure");

	BBFE_manusol_calc_nodal_error_scalar(fe, error_p, rom_pressure, hrom_pressure);
	BB_vtk_write_point_vals_scalar(fp, error_p, fe->total_num_nodes, "abs_error-ROM_HROM-Pressure");

	BB_vtk_write_point_vals_scalar(fp, fem_pressure, fe->total_num_nodes, "FEM-Pressure");
	BB_std_free_1d_double(error_p, fe->total_num_nodes);

	fclose(fp);
}


void HROM_add_output_vtk_velocity(
		BBFE_DATA*     fe,
		double**       fem_velocity,
		double**	   rom_velocity,
        double**	   hrom_velocity,
		const char*    filename,
		const char*    directory,
		double         t)
{
	FILE* fp;
	fp = ROM_BB_write_add_fopen(fp, filename, directory);

	double** error_v;
	error_v = BB_std_calloc_2d_double(error_v, fe->total_num_nodes, 3);
	for (int i = 0; i < fe->total_num_nodes; i++){
		for (int j = 0; j < 3; j++){
			error_v[i][j] =	abs(fem_velocity[i][j] - rom_velocity[i][j]);
		}
	}
	
	/*for velocity*/
   	BB_vtk_write_point_vals_vector(fp, fem_velocity, fe->total_num_nodes, "FEM-Velocity");
	BB_vtk_write_point_vals_vector(fp, rom_velocity, fe->total_num_nodes, "ROM-Velocity");
    BB_vtk_write_point_vals_vector(fp, hrom_velocity, fe->total_num_nodes, "HROM-Velocity");

	BB_vtk_write_point_vals_vector(fp, error_v, fe->total_num_nodes, "abs_error-FEM_ROM-Velosity");

	for (int i = 0; i < fe->total_num_nodes; i++){
		for (int j = 0; j < 3; j++){
			error_v[i][j] =	abs(fem_velocity[i][j] - hrom_velocity[i][j]);
		}
	}

    BB_vtk_write_point_vals_vector(fp, error_v, fe->total_num_nodes, "abs_error-FEM_HROM-Velosity");

	for (int i = 0; i < fe->total_num_nodes; i++){
		for (int j = 0; j < 3; j++){
			error_v[i][j] =	abs(rom_velocity[i][j] - hrom_velocity[i][j]);
		}
	}

    BB_vtk_write_point_vals_vector(fp, error_v, fe->total_num_nodes, "abs_error-ROM_HROM-Velosity");

	fclose(fp);

	BB_std_free_2d_double(error_v, fe->total_num_nodes, 3);
}


void HROM_output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	const char* filename;
	char fname_vtk[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, ROM_OUTPUT_FILENAME_VTK, file_num);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);

	HROM_output_vtk_shape(
			&(sys->fe),
			filename,
			sys->cond.directory,
			t);
	HROM_add_output_vtk_pressure(
			&(sys->fe),
			sys->vals.p,
			sys->vals_rom.p,
            sys->vals_hrom.p,
			filename,
			sys->cond.directory,
			t);
	HROM_add_output_vtk_velocity(
			&(sys->fe),
			sys->vals.v,
			sys->vals_rom.v,
            sys->vals_hrom.v,
			filename,
			sys->cond.directory,
			t);

	double L2_error_p = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			&(sys->mono_com),
			t,
			(const double*)sys->vals.p,
			(const double*)sys->vals_rom.p);

	printf("%s L2 error pressure FEM-ROM: %e\n", CODENAME, L2_error_p);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_pressure_fem-rom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_p);
		fclose(fp);
	}

	double L2_error_v = ROM_sys_hlpod_fe_equivval_relative_L2_error_vector(
			&(sys->fe),
			&(sys->basis),
			&(sys->mono_com),
			t,
			(const double**)sys->vals.v,
			(const double**)sys->vals_rom.v);

	printf("%s L2 error velocity FEM-HROM: %e\n", CODENAME, L2_error_v);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_velocity_fem-rom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_v);
		fclose(fp);
	}

    L2_error_p = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			&(sys->mono_com),
			t,
			(const double*)sys->vals.p,
			(const double*)sys->vals_hrom.p);

	printf("%s L2 error pressure FEM-HROM: %e\n", CODENAME, L2_error_p);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_pressure_fem-hrom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_p);
		fclose(fp);
	}

	L2_error_v = ROM_sys_hlpod_fe_equivval_relative_L2_error_vector(
			&(sys->fe),
			&(sys->basis),
			&(sys->mono_com),
			t,
			(const double**)sys->vals.v,
			(const double**)sys->vals_hrom.v);

	printf("%s L2 error velocity FEM-HROM: %e\n", CODENAME, L2_error_v);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_velocity_fem-hrom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_v);
		fclose(fp);
	}

    L2_error_p = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			&(sys->mono_com),
			t,
			(const double*)sys->vals_rom.p,
			(const double*)sys->vals_hrom.p);

	printf("%s L2 error pressure ROM-HROM: %e\n", CODENAME, L2_error_p);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_pressure_rom-hrom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_p);
		fclose(fp);
	}

	L2_error_v = ROM_sys_hlpod_fe_equivval_relative_L2_error_vector(
			&(sys->fe),
			&(sys->basis),
			&(sys->mono_com),
			t,
			(const double**)sys->vals_rom.v,
			(const double**)sys->vals_hrom.v);

	printf("%s L2 error velocity ROM-HROM: %e\n", CODENAME, L2_error_v);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_velocity_rom-hrom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_v);
		fclose(fp);
	}
}


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
    ddhr_set_matvec_RH_for_NNLS_para_volume_const(
        &(sys->fe),
        &(sys->vals),
        &(sys->basis),
        &(rom->hlpod_mat),
        &(rom->hlpod_vals),
        &(hrom->hlpod_ddhr),
        rom->hlpod_vals.num_2nd_subdomains,
        0 ,   //index 0 start
        rom->hlpod_vals.num_snapshot,
        1 + sys->mono_com.recv_n_neib,
        sys->vals.dt,
        0);

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
        4,
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
        sys->fe.total_num_nodes,
        4);

    monolis_get_nonzero_pattern_by_nodal_graph_V_R(
        &(sys->monolis_hr0),
        rom->hlpod_meta.num_meta_nodes,
        rom->hlpod_meta.n_dof_list,
        rom->hlpod_meta.index,
        rom->hlpod_meta.item);
}

void HROM_hierarchical_parallel(
    FE_SYSTEM   sys,
    ROM*        rom,
    HROM*       hrom,
    const int   step_HR,
    const int   step_POD,
    const double t)
{
    monolis_initialize(&(sys.monolis_hr));
    monolis_copy_mat_nonzero_pattern_R(&(sys.monolis_hr0), &(sys.monolis_hr));

    ddhr_lb_set_reduced_mat_para_save_memory(
        &(sys.monolis_hr),
        &(sys.fe),
        &(sys.vals_hrom),
        &(sys.basis),
        &(sys.bc),
        &(rom->hlpod_vals),
        &(rom->hlpod_mat),
        &(hrom->hlpod_ddhr),
        rom->hlpod_vals.num_modes_pre,
        rom->hlpod_vals.num_2nd_subdomains,
        sys.vals.dt);

    ddhr_hlpod_calc_block_mat_bcsr_pad(
        &(sys.monolis_hr),
        &(sys.mono_com_rom_solv),
        &(rom->hlpod_vals),
        &(rom->hlpod_mat),
        &(hrom->hlpod_ddhr),
        &(rom->hlpod_meta),
        rom->hlpod_vals.num_modes_pre,
        rom->hlpod_vals.num_2nd_subdomains,
        sys.cond.directory);
	
    ddhr_lb_set_reduced_vec_para(
        &(sys.monolis_hr),
        &(sys.fe),
        &(sys.vals_hrom),
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
        &(sys.vals_hrom),
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

    monolis_show_iterlog (&(sys.monolis_hr), false);

    BBFE_sys_monowrap_solve(
        &(sys.monolis_hr),
        &(sys.mono_com_rom_solv),
        rom->hlpod_mat.mode_coef,
        MONOLIS_ITER_BICGSTAB_N128,
        MONOLIS_PREC_DIAG,
        sys.vals.mat_max_iter,
        sys.vals.mat_epsilon);

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
		4);

	monolis_mpi_update_R(&(sys.mono_com), sys.fe.total_num_nodes, 4, hrom->hr_vals.sol_vec);

    BBFE_fluid_sups_renew_velocity(
        sys.vals_hrom.v, 
        hrom->hr_vals.sol_vec,
        sys.fe.total_num_nodes);

    BBFE_fluid_sups_renew_pressure(
        sys.vals_hrom.p, 
        hrom->hr_vals.sol_vec,
        sys.fe.total_num_nodes);

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
    printf("num_modes_pre: %d\n", rom->hlpod_vals.num_modes_pre);
    printf("num_2nd_subdomains: %d\n", rom->hlpod_vals.num_2nd_subdomains);
    printf("num_modes_max: %d\n", rom->hlpod_vals.num_modes_max);
    printf("num_snapshot: %d\n", rom->hlpod_vals.num_snapshot);
    //exit(1);
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


void HROM_pre_offline3(
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



void HROM_std_hlpod_pre_lpod_para(
        MONOLIS*     monolis_rom0,
        MONOLIS_COM* monolis_com,
        MONOLIS_COM* mono_com_rom,
        MONOLIS_COM* mono_com_rom_solv,
        ROM*		 rom,
        const int    dof,
        const char*	 directory)
{
    ROM_std_hlpod_get_neib_vec(
            monolis_com,
            &(rom->hlpod_vals),
            &(rom->hlpod_mat),
            rom->hlpod_vals.num_modes,
            dof);

    ROM_std_hlpod_get_neib_num_modes_para_subd(
            mono_com_rom,
            &(rom->hlpod_vals),
            &(rom->hlpod_mat),
            1 + monolis_com->recv_n_neib,
            rom->hlpod_vals.num_modes);

    ROM_std_hlpod_get_neib_num_modes_mode_subd(
            mono_com_rom,
            mono_com_rom_solv,
            &(rom->hlpod_mat),
            &(rom->hlpod_meta),
            1 + monolis_com->recv_n_neib,
            directory);

    ROM_std_hlpod_get_n_dof_list(
            mono_com_rom_solv,
            &(rom->hlpod_mat),
            &(rom->hlpod_meta),
            rom->hlpod_vals.num_modes_pre);

    monolis_get_nonzero_pattern_by_nodal_graph_V_R(
            monolis_rom0,
            rom->hlpod_meta.num_meta_nodes,
            rom->hlpod_meta.n_dof_list,
            rom->hlpod_meta.index,
            rom->hlpod_meta.item);
}


void HROM_std_hlpod_online_pre(
        MONOLIS*     monolis_rom0,
        MONOLIS_COM* mono_com,
        MONOLIS_COM* mono_com_rom,
        MONOLIS_COM* mono_com_rom_solv,
        ROM* 		 rom_sups,
        const int 	 total_num_nodes,
        const int 	 dof,
        const char*  metagraph,
        const char*  directory)
{
    ROM_std_hlpod_read_metagraph(
			monolis_rom0,
			mono_com_rom_solv,
			rom_sups,
			metagraph,
			directory);

    if(monolis_mpi_get_global_comm_size() == 1){
    
        if(rom_sups->hlpod_vals.num_1st_subdomains==0){
            printf("\nError : num_1st_subdomains is not set\n");
            exit(1);
        }
        else{
            
        }
    }		
    else{
        if(rom_sups->hlpod_vals.bool_global_mode==false){

            HROM_std_hlpod_pre_lpod_para(
                    monolis_rom0,
                    mono_com,
                    mono_com_rom,
                    mono_com_rom_solv,
                    rom_sups,
                    dof,
                    directory);

        }
        else{				
            ROM_std_hlpod_update_global_modes(
                    mono_com,
                    &(rom_sups->hlpod_mat),
                    total_num_nodes,
                    mono_com->n_internal_vertex,
                    rom_sups->hlpod_vals.num_modes_pre,
                    4);
                
        }
    }
}

