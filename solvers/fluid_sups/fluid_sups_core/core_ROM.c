
#include "core_ROM.h"

/*for ROM*/
const char*     ROM_ID_FINISH_TIME      = "#rom_finish_time";
const double  POD_DVAL_FINISH_TIME      = 1.0;
const char* POD_ID_OUTPUT_INTERVAL      = "#rom_output_interval";
const int POD_DVAL_OUTPUT_INTERVAL      = 1;
const char* POD_ID_SNAPSHOT_INTERVAL    = "#snapshot_interval";
const int POD_DVAL_SNAPSHOT_INTERVAL    = 1;

const char*         ROM_ID_DENSITY  = "#density";
const char*       ROM_ID_VISCOSITY  = "#viscosity";

static const char* POD_INPUT_FILENAME_COND          = "rom_cond.dat";
static const char* ROM_OUTPUT_FILENAME_VTK          = "rom_result_%06d.vtk";

static const char* INPUT_FILENAME_DENSITY           = "density.dat";    //for parametric study
static const char* INPUT_FILENAME_VISCOSITY         = "viscosity.dat";  //for parametric study

static const char* INPUT_DIRECTORYNAME_1STDD        = "parted.0/";
static const char* INPUT_DIRECTORYNAME_2NDD         = "parted.1/";
static const char* INPUT_DIRECTORYNAME_METAGRAPH    = "metagraph_parted.0/";
static const char* INPUT_FILENAME_METAGRAPH         = "metagraph.dat";

const int BUFFER_SIZE = 10000;


void initialize_velocity_pressure(
	double**    v,
	double*     p,
	const int   total_num_nodes)
{
    for (int i = 0; i < total_num_nodes; i++) {
        for (int j = 0; j < 3; j++) {
            v[i][j] = 0.0;
        }
    }

    for (int i = 0; i < total_num_nodes; i++) {
        p[i] = 0.0;
    }
}

void ROM_offline_set_reynolds_num_cases(
    VALUES*         vals,
    const char*     directory)
{
	FILE* fp;
	char id[BUFFER_SIZE];
	int num_case_density;
	int num_case_viscosity;

	fp = ROM_BB_read_fopen(fp, INPUT_FILENAME_DENSITY, directory);

	fscanf(fp, "%s %d", id, &(vals->num_cases));
	num_case_density = vals->num_cases;
	vals->density_cases = BB_std_calloc_1d_double(vals->density_cases, vals->num_cases);
	for (int i = 0; i < vals->num_cases; i++) {
		fscanf(fp, "%lf", &(vals->density_cases[i]));
	}
	fclose(fp);

	fp = ROM_BB_read_fopen(fp, INPUT_FILENAME_VISCOSITY, directory);

	fscanf(fp, "%s %d", id, &(num_case_viscosity));
	if(num_case_density != num_case_viscosity){
		printf("Error: The number of cases in density.dat and viscosity.dat are different.\n");
		exit(1);
	}
	vals->viscosity_cases = BB_std_calloc_1d_double(vals->viscosity_cases, vals->num_cases);
	for (int i = 0; i < vals->num_cases; i++) {
		fscanf(fp, "%lf", &(vals->viscosity_cases[i]));
	}
	fclose(fp);

}

void ROM_offline_set_reynolds_number( 
        VALUES*     vals,
        const int   case_id)
{
    vals->density   = vals->density_cases[case_id];
    vals->viscosity = vals->viscosity_cases[case_id];

    printf("\n%s ---------- Calculation Conditions: Parametric Study Case %d ----------\n", 
            CODENAME, case_id);

    printf("%s %s: %e\n", CODENAME, ROM_ID_DENSITY,   vals->density);
    printf("%s %s: %e\n", CODENAME, ROM_ID_VISCOSITY, vals->viscosity);
}


const char* ROM_std_hlpod_get_parted_file_name(
    const int   solver_type)
{
    static char fname[BUFFER_SIZE];

    if(solver_type==1){
        return NULL;
    }
    else if(solver_type==2){
        snprintf(fname, BUFFER_SIZE,"%s", INPUT_DIRECTORYNAME_1STDD);
        return fname;
    }
    else if(solver_type==3){
        snprintf(fname, BUFFER_SIZE,"%s", INPUT_DIRECTORYNAME_2NDD);
        return fname;
    }
    else if(solver_type==4){
        snprintf(fname, BUFFER_SIZE,"%s", INPUT_DIRECTORYNAME_1STDD);
        return fname;
    }
    else{
        printf("Error: solver type is not set (solver type = %d)\n", solver_type);
        exit(1);
    }

}

const char* ROM_std_hlpod_get_metagraph_name(
    const int   solver_type)
{
    static char fname[BUFFER_SIZE];

    if(solver_type==1){
        return NULL;
    }
    else if(solver_type==2){
        snprintf(fname, BUFFER_SIZE,"%s/%s", INPUT_DIRECTORYNAME_1STDD, INPUT_FILENAME_METAGRAPH);
        return fname;
    }
    else if(solver_type==3){
        snprintf(fname, BUFFER_SIZE,"%s/%s", INPUT_DIRECTORYNAME_METAGRAPH ,INPUT_FILENAME_METAGRAPH);
        return fname;
    }
    else if(solver_type==4){
        return NULL;
    }
    else{
        printf("Error: solver type is not set (solver type = %d)\n", solver_type);
        exit(1);
    }

}


void set_target_parameter(
		VALUES*         vals,
		const char*     directory)
{
	FILE* fp;
	char fname[BUFFER_SIZE];
	char id[BUFFER_SIZE];
	int tmp;
	double target_vals;

	snprintf(fname, BUFFER_SIZE, "target_density.dat");
	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s %d", id, &(tmp));
	fscanf(fp, "%lf", &(target_vals));
	vals->density = target_vals;
	fclose(fp);

	snprintf(fname, BUFFER_SIZE, "target_viscosity.dat");
	fp = ROM_BB_read_fopen(fp, fname, directory);
	fscanf(fp, "%s %d", id, &(tmp));
	fscanf(fp, "%lf", &(target_vals));
	vals->viscosity = target_vals;
	fclose(fp);

	printf("\n%s ---------- Calculation condition - ROM %d  ----------\n", CODENAME);

	printf("%s %s: %e\n", CODENAME, ROM_ID_DENSITY,          vals->density);
	printf("%s %s: %e\n", CODENAME, ROM_ID_VISCOSITY,        vals->viscosity);

}


void ROM_set_param(
    ROM*            rom,
    const int       num_subdomains,
    const int       num_modes,
    const double    rom_epsilon,
    const int       solver_type)
{
    rom->hlpod_vals.num_1st_subdomains = num_subdomains;
    rom->hlpod_vals.num_2nd_subdomains = num_subdomains;
    rom->hlpod_vals.num_modes_pre = num_modes;
    rom->hlpod_vals.rom_epsilon = 1 - rom_epsilon;

    rom->hlpod_vals.bool_global_mode = false;
    
    if(solver_type == 4){
        rom->hlpod_vals.bool_global_mode = true;
    }

}

void ROM_offline_assign_default_values(
		VALUES*     vals)
{
	vals->snapshot_interval  = POD_DVAL_SNAPSHOT_INTERVAL;
}


void ROM_offline_print_all_values(
		VALUES*     vals)
{
	printf("\n%s ---------- Calculation condition POD----------\n", CODENAME);

	printf("%s %s: %d\n", CODENAME, POD_ID_SNAPSHOT_INTERVAL,  vals->snapshot_interval);

	printf("%s -------------------------------------------\n\n", CODENAME);
}


void ROM_offline_read_calc_conditions(
		VALUES*         vals,
		const char* 	directory)
{
	printf("\n");

	ROM_offline_assign_default_values(vals);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, POD_INPUT_FILENAME_COND);

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s Calc condition file \"%s\" is not found.\n", CODENAME, filename);
		printf("%s Default values are used in this calculation.\n", CODENAME);
	}
	else {
		printf("%s Reading POD conditon file \"%s\".\n", CODENAME, filename);
		int num;

		num = BB_std_read_file_get_val_int_p(
				&(vals->snapshot_interval), filename, POD_ID_SNAPSHOT_INTERVAL, BUFFER_SIZE, CODENAME);

		fclose(fp);
	}

	ROM_offline_print_all_values(vals);

	printf("\n");
}


void ROM_online_assign_default_values(
		VALUES*     vals)
{
	vals->rom_finish_time  = POD_DVAL_FINISH_TIME;
	vals->rom_output_interval  = POD_DVAL_OUTPUT_INTERVAL;
}


void ROM_online_print_all_values(
		VALUES*     vals)
{
	printf("\n%s ---------- Calculation condition POD----------\n", CODENAME);

	printf("%s %s: %e\n", CODENAME, ROM_ID_FINISH_TIME,  vals->rom_finish_time);
	printf("%s %s: %d\n", CODENAME, ROM_ID_FINISH_TIME,  vals->rom_output_interval);

	printf("%s -------------------------------------------\n\n", CODENAME);
}


void ROM_online_read_calc_conditions(
		VALUES*         vals,
		const char* 	directory)
{
	printf("\n");

	ROM_online_assign_default_values(vals);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, POD_INPUT_FILENAME_COND);

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s Calc condition file \"%s\" is not found.\n", CODENAME, filename);
		printf("%s Default values are used in this calculation.\n", CODENAME);
	}
	else {
		printf("%s Reading POD conditon file \"%s\".\n", CODENAME, filename);
		int num;

		num = BB_std_read_file_get_val_double_p(
				&(vals->rom_finish_time), filename, ROM_ID_FINISH_TIME, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_int_p(
				&(vals->rom_output_interval), filename, POD_ID_OUTPUT_INTERVAL, BUFFER_SIZE, CODENAME);

		fclose(fp);
	}

	ROM_online_print_all_values(vals);

	printf("\n");
}

void ROM_output_vtk_shape(
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


void ROM_add_output_vtk_pressure(
		BBFE_DATA*     fe,
		double*        fem_pressure,
		double*    	   rom_pressure,
		const char*    filename,
		const char*    directory,
		double         t)
{
	FILE* fp;
	fp = ROM_BB_write_add_fopen(fp, filename, directory);

	/*for pressure*/
	BB_vtk_write_point_vals_scalar(fp, rom_pressure, fe->total_num_nodes, "ROM-Pressure");

	double* error_p;
	error_p = BB_std_calloc_1d_double(error_p, fe->total_num_nodes);
	BBFE_manusol_calc_nodal_error_scalar(fe, error_p, fem_pressure, rom_pressure);
	BB_vtk_write_point_vals_scalar(fp, error_p, fe->total_num_nodes, "abs_error-FEM_ROM-Pressure");

	BB_vtk_write_point_vals_scalar(fp, fem_pressure, fe->total_num_nodes, "FEM-Pressure");
	BB_std_free_1d_double(error_p, fe->total_num_nodes);

	fclose(fp);
}


void ROM_add_output_vtk_velocity(
		BBFE_DATA*     fe,
		double**       fem_velocity,
		double**	   rom_velocity,
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
	BB_vtk_write_point_vals_vector(fp, rom_velocity, fe->total_num_nodes, "ROM-Velocity");
	BB_vtk_write_point_vals_vector(fp, fem_velocity, fe->total_num_nodes, "FEM-Velocity");
	BB_vtk_write_point_vals_vector(fp, error_v, fe->total_num_nodes, "abs_error-FEM_ROM-Velosity");

	fclose(fp);

	BB_std_free_2d_double(error_v, fe->total_num_nodes, 3);
}


void ROM_output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	const char* filename;
	char fname_vtk[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, ROM_OUTPUT_FILENAME_VTK, file_num);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);
/*
	ROM_output_vtk_shape(
			&(sys->fe),
			filename,
			sys->cond.directory,
			t);
	ROM_add_output_vtk_pressure(
			&(sys->fe),
			sys->vals.p,
			sys->vals_rom.p,
			filename,
			sys->cond.directory,
			t);
	ROM_add_output_vtk_velocity(
			&(sys->fe),
			sys->vals.v,
			sys->vals_rom.v,
			filename,
			sys->cond.directory,
			t);
	*/
	double L2_error_p = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
			&(sys->fe),
			&(sys->basis),
			&(sys->mono_com),
			t,
			(const double*)sys->vals_rom.p,
			(const double*)sys->vals.p);

	printf("%s L2 error: %e\n", CODENAME, L2_error_p);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_pressure.txt", sys->cond.directory);
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

	printf("%s L2 error: %e\n", CODENAME, L2_error_v);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_velocity.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_v);
		fclose(fp);
	}
}


void solver_rom(
    FE_SYSTEM* sys,
    const int step_HR,
    const int step_rom,
    const double t)
{
    monolis_clear_mat_value_R(&(sys->monolis));
	monolis_clear_mat_value_R(&(sys->monolis_rom));
	monolis_com_initialize_by_self(&(sys->mono_com0));

	if(monolis_mpi_get_global_comm_size() == 1){
	}
	else{
		monolis_copy_mat_value_R(&(sys->monolis_rom0), &(sys->monolis_rom));
	} 

    set_element_mat(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->vals_rom));
        
    set_element_vec(
        &(sys->monolis),
        &(sys->fe),
        &(sys->basis),
        &(sys->vals_rom));

    BBFE_sys_monowrap_set_Dirichlet_bc(
        &(sys->monolis),
        sys->fe.total_num_nodes,
        4,
        &(sys->bc),
        sys->monolis.mat.R.B);
	
	ROM_std_hlpod_calc_reduced_mat(
		&(sys->monolis),
        &(sys->monolis_rom),
        &(sys->mono_com),
		&(sys->mono_com0),
        &(sys->mono_com_rom_solv),
		&(sys->rom_sups),
        sys->fe.total_num_nodes,
		4);
	
    ROM_std_hlpod_solve_ROM(
        &(sys->monolis),
        &(sys->monolis_rom),
        &(sys->mono_com_rom_solv),
		&(sys->rom_sups),
        sys->fe.total_num_nodes,
		4,
		sys->vals.mat_max_iter,
		sys->vals.mat_epsilon,
		MONOLIS_ITER_BICGSTAB,
		MONOLIS_PREC_DIAG);

	ROM_sys_hlpod_fe_add_Dbc(
        sys->rom_sups.hlpod_vals.sol_vec,
		&(sys->bc),
		sys->fe.total_num_nodes,
		4);
	
	monolis_mpi_update_R(&(sys->mono_com), sys->fe.total_num_nodes, 4, sys->rom_sups.hlpod_vals.sol_vec);

    BBFE_fluid_sups_renew_velocity(
        sys->vals_rom.v, 
        sys->rom_sups.hlpod_vals.sol_vec,
        sys->fe.total_num_nodes);

    BBFE_fluid_sups_renew_pressure(
        sys->vals_rom.p, 
        sys->rom_sups.hlpod_vals.sol_vec,
        sys->fe.total_num_nodes);
        
    output_files(sys, step_rom, t);

}


void solver_rom_global_para(
	MONOLIS*     monolis,
	MONOLIS_COM* monolis_com,
	ROM*		 rom,
    FE_SYSTEM*   sys,
    const int    step_HR,
    const int    step_rom,
    const double t)
{
    int total_num_nodes = sys->fe.total_num_nodes;
    double** mat = BB_std_calloc_2d_double(mat, rom->hlpod_vals.num_modes, rom->hlpod_vals.num_modes);
	double* rhs = BB_std_calloc_1d_double(rhs, rom->hlpod_vals.num_modes);
    double* mode_coef = BB_std_calloc_1d_double(mode_coef, rom->hlpod_vals.num_modes);
    double* ansvec = BB_std_calloc_1d_double(ansvec, total_num_nodes* 4);

	set_D_bc_global_para(
			monolis,
			monolis_com,
			&(sys->fe),
			&(sys->basis),
			&(sys->vals_rom),
			&(sys->bc),
            rhs,
            rom->hlpod_mat.pod_modes,
			rom->hlpod_vals.num_modes);
	
	set_reduced_vec_global_para(
			monolis,
			monolis_com,
			&(sys->fe),
			&(sys->basis),
			&(sys->vals_rom),
            rhs,
            rom->hlpod_mat.pod_modes,
			rom->hlpod_vals.num_modes);

	set_reduced_mat_global_para(
			monolis,
			monolis_com,
			&(sys->fe),
			&(sys->basis),
			&(sys->vals_rom),
            mat,
            rom->hlpod_mat.pod_modes,
			rom->hlpod_vals.num_modes);
	
    allreduce_global_para(
            monolis_com,
            mat,
            rhs,
            rom->hlpod_vals.num_modes);

	ROM_BB_gauss_elimination(
			rom->hlpod_vals.num_modes,
			mat,
			rhs,
			mode_coef);
		
    ROM_std_hlpod_calc_sol_global_para(
            ansvec,
            rom->hlpod_mat.pod_modes,
            mode_coef,
			sys->fe.total_num_nodes,
			rom->hlpod_vals.num_modes,
			4);
	
	ROM_sys_hlpod_fe_add_Dbc(
            ansvec,
			&(sys->bc),
			sys->fe.total_num_nodes,
			4);
				
	/*解ベクトルのupdate*/
	monolis_mpi_update_R(monolis_com, sys->fe.total_num_nodes, 4, ansvec);

	BB_std_free_2d_double(mat, rom->hlpod_vals.num_modes, rom->hlpod_vals.num_modes);
	BB_std_free_1d_double(rhs, rom->hlpod_vals.num_modes);
	BB_std_free_1d_double(mode_coef, rom->hlpod_vals.num_modes);

    BBFE_fluid_sups_renew_velocity(
        sys->vals_rom.v, 
        ansvec,
        sys->fe.total_num_nodes);

    BBFE_fluid_sups_renew_pressure(
        sys->vals_rom.p, 
        ansvec,
        sys->fe.total_num_nodes);
    
    BB_std_free_1d_double(ansvec, total_num_nodes* 4);
}