
#include "core_ROM.h"

static const char* OPTION_NUM_MODES     = "-nm";
static const char* OPTION_NUM_1STDD     = "-nd";
static const char* OPTION_PADAPTIVE     = "-pa";
static const char* OPTION_SOLVER_TYPE   = "-st";
static const char* OPTION_HOT_START     = "-hs";

static const char* INPUT_DIRECTORYNAME_METAGRAPH = "metagraph_parted.0/";
static const char* INPUT_FILENAME_METAGRAPH      = "metagraph.dat";

static const char* INPUT_FILENAME_COND    = "cond.dat";
static const char* INPUT_FILENAME_D_BC_V  = "D_bc_v.dat";
static const char* INPUT_FILENAME_D_BC_P  = "D_bc_p.dat";

const int BUFFER_SIZE = 1024;

void ROM_read_args(
    int 		argc,
    char* 		argv[],
    ROM_PRM*    rom_prm)
{
	int num;
    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_NUM_1STDD);
    if(num == -1) {
		printf("\nargs error num_subdomains");
		exit(1);
    }
    else {
        rom_prm->num_subdomains = atoi(argv[num+1]);
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_NUM_MODES);
    if(num == -1) {
		printf("\nargs error num_modes");
		exit(1);
    }
    else {
		rom_prm->num_modes = atoi(argv[num+1]);
    }

    num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_PADAPTIVE);
    if(num == -1) {
		printf("\nargs error rom_epsilon");
		exit(1);
    }
    else {
        rom_prm->rom_epsilon = atof(argv[num+1]);
    }

	num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_SOLVER_TYPE);
    if(num == -1) {
		printf("\nargs error solver_type");
		exit(1);
    }
    else {
        rom_prm->solver_type = atof(argv[num+1]);
    }

	num = BB_std_read_args_return_char_num(
                argc, argv, OPTION_HOT_START);
    if(num == -1) {
        rom_prm->hot_start = 0;
    }
    else {
        rom_prm->hot_start = 1;
        //rom_prm->hot_start_time = atof(argv[num+1]);
        printf("Hot start is enabled.\n");
    }

	printf("num_subdomains = %d\n", rom_prm->num_subdomains);
	printf("num_modes = %d\n", rom_prm->num_modes);
	printf("rom_epsilon = %lf\n", rom_prm->rom_epsilon);
    printf("solver_type = %d\n", rom_prm->solver_type);
    printf("hot_start = %d\n", rom_prm->hot_start);

}

double hot_start_read_initialize_val(
    double*     int_val,
    const char* input_fname,
    const char* directory)
{
    int BUFFER_SIZE = 1024;
    int total_num_nodes;
    int ndof;
    double t = 0.0;
	FILE* fp;
	char fname[BUFFER_SIZE];
	char id[BUFFER_SIZE];

	fp = BBFE_sys_read_fopen(fp, "hot_start/start_time.dat", directory);
	fscanf(fp, "%s", id);
    fscanf(fp, "%lf", &(t));
    fclose(fp);

	fp = BBFE_sys_read_fopen(fp, input_fname, directory);
	fscanf(fp, "%s", id);
    fscanf(fp, "%d %d", &(total_num_nodes), &(ndof));
    for(int i = 0; i < total_num_nodes; i++) {
        for(int j = 0; j < ndof; j++) {
            fscanf(fp, "%lf", &(int_val[i * ndof + j]));
        }
    }
	fclose(fp);

    return t;
}

void hot_start_write_initialize_val(
    double*         int_val,
    const int       total_num_nodes,
    const int       ndof,
    const double    time,
    const char*     output_fname,
    const char*     directory)
{
	FILE* fp;
	char fname[BUFFER_SIZE];
	char id[BUFFER_SIZE];

	fp = BBFE_sys_write_fopen(fp, output_fname, directory);
	fprintf(fp, "initialization\n");
    fprintf(fp, "%d %d\n", total_num_nodes, ndof);
    for(int i = 0; i < total_num_nodes; i++) {
        for(int j = 0; j < ndof; j++) {
            fprintf(fp, "%e ", int_val[i * ndof + j]);
        }
        fprintf(fp, "\n");        
    }
	fclose(fp);

    if(monolis_mpi_get_global_my_rank()==0){
        fp = BBFE_sys_write_fopen(fp, "hot_start/start_time.dat", directory);
        fprintf(fp, "start_time\n");
        fprintf(fp, "%lf", time);
        fclose(fp);
    }
}

int main (
		int argc,
		char* argv[])
{
	printf("\n");

	FE_SYSTEM sys;

	monolis_global_initialize();

    double t1 = monolis_get_time();
	double FOM_t1 = monolis_get_time();

	sys.cond.directory = BBFE_fluid_get_directory_name(argc, argv, CODENAME);	

	read_calc_conditions(&(sys.vals), sys.cond.directory);

	BBFE_fluid_pre(
			&(sys.fe), &(sys.basis),
			argc, argv, sys.cond.directory,
			sys.vals.num_ip_each_axis);

	const char* filename;	
	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC_V);
	BBFE_fluid_sups_read_Dirichlet_bc(
			&(sys.bc),
			filename,
			sys.cond.directory,
			sys.fe.total_num_nodes,
			4);
    
    memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);

	BBFE_elemmat_set_Jacobi_mat(&(sys.fe), &(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(&(sys.fe), &(sys.basis));

	BBFE_sys_monowrap_init_monomat(&(sys.monolis) , &(sys.mono_com), &(sys.fe), 4, sys.cond.directory);

	//intialize for velocity and pressure
	initialize_velocity_pressure_karman_vortex(sys.vals.v, sys.vals.p, sys.fe.total_num_nodes);

double t = 0.0;
	int step = 0;
	int file_num = 0;
	while (t < sys.vals.finish_time) {
		t += sys.vals.dt;
		step += 1;

		printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);

		monolis_clear_mat_value_R(&(sys.monolis));

		printf("%s --- prediction step ---\n", CODENAME);

		set_element_mat(
				&(sys.monolis),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));
		set_element_vec(
				&(sys.monolis),
				&(sys.fe),
				&(sys.basis),
				&(sys.vals));

		BBFE_sys_monowrap_set_Dirichlet_bc(
				&(sys.monolis),
				sys.fe.total_num_nodes,
				4,
				&(sys.bc),
				sys.monolis.mat.R.B);

		BBFE_sys_monowrap_solve(
				&(sys.monolis),
				&(sys.mono_com),
				sys.monolis.mat.R.X,
				MONOLIS_ITER_BICGSTAB,
				MONOLIS_PREC_DIAG,
				sys.fe.total_num_nodes,
				sys.vals.mat_epsilon);

		BBFE_fluid_sups_renew_velocity(
				sys.vals.v, 
				sys.monolis.mat.R.X,
				sys.fe.total_num_nodes);
        
        BBFE_fluid_sups_renew_pressure(
				sys.vals.p, 
				sys.monolis.mat.R.X,
				sys.fe.total_num_nodes);

		/**********************************************/

		if(step%sys.vals.output_interval == 0) {

			BBFE_fluid_sups_renew_pressure(
				sys.vals.p, 
				sys.monolis.mat.R.X,
				sys.fe.total_num_nodes);

			output_files(&sys, file_num, t);
			file_num += 1;

			// cavity
			//output_result_file_cavity_center_vx(&(sys.vals), sys.cond.directory);
		}

	}

	BBFE_fluid_finalize(&(sys.fe), &(sys.basis));
	BBFE_sys_memory_free_Dirichlet_bc(&(sys.bc), sys.fe.total_num_nodes, 4);
	monolis_finalize(&(sys.monolis));

	double t2 = monolis_get_time();
	int myrank = monolis_mpi_get_global_my_rank();

	if(myrank == 0) {
		printf("** Total time: %f\n", t2 - t1);
	}

	monolis_global_finalize();

	printf("\n");

	return 0;

}

