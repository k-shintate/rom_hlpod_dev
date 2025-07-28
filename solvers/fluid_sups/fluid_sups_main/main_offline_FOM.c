
#include "core_ROM.h"

static const char* OPTION_NUM_MODES     = "-nm";
static const char* OPTION_NUM_1STDD     = "-nd";
static const char* OPTION_PADAPTIVE     = "-pa";
static const char* OPTION_SOLVER_TYPE   = "-st";

static const char* INPUT_DIRECTORYNAME_METAGRAPH = "metagraph_parted.0/";
static const char* INPUT_FILENAME_METAGRAPH      = "metagraph.dat";

static const char* INPUT_FILENAME_COND    = "cond.dat";
static const char* INPUT_FILENAME_D_BC_V  = "D_bc_v.dat";
static const char* INPUT_FILENAME_D_BC_P  = "D_bc_p.dat";


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

	printf("num_subdomains = %d\n", rom_prm->num_subdomains);
	printf("num_modes = %d\n", rom_prm->num_modes);
	printf("rom_epsilon = %lf\n", rom_prm->rom_epsilon);
    printf("solver_type = %d\n", rom_prm->solver_type);

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

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);
	
	filename = monolis_get_global_input_file_name(MONOLIS_DEFAULT_TOP_DIR, MONOLIS_DEFAULT_PART_DIR, INPUT_FILENAME_D_BC_V);
	BBFE_fluid_sups_read_Dirichlet_bc(
			&(sys.bc),
			filename,
			sys.cond.directory,
			sys.fe.total_num_nodes,
			4);

	BBFE_elemmat_set_Jacobi_mat(&(sys.fe), &(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(&(sys.fe), &(sys.basis));

	BBFE_sys_monowrap_init_monomat(&(sys.monolis) , &(sys.mono_com), &(sys.fe), 4, sys.cond.directory);

	//intialize for velocity and pressure
	initialize_velocity_pressure(sys.vals.v, sys.vals.p, sys.fe.total_num_nodes);

    /*for ROM *****************************************/
	
    /*for ROM input data*/
	ROM_read_args(argc, argv, &(sys.rom_prm_p));
	ROM_read_args(argc, argv, &(sys.rom_prm_v));

    ROM_set_param(
            &(sys.rom_p),
            sys.rom_prm_p.num_subdomains,
            sys.rom_prm_p.num_modes,
            sys.rom_prm_p.rom_epsilon,
            sys.rom_prm_p.solver_type);
    
    ROM_set_param(
            &(sys.rom_v),
            sys.rom_prm_v.num_subdomains,
            sys.rom_prm_v.num_modes,
            sys.rom_prm_v.rom_epsilon,
            sys.rom_prm_v.solver_type);
    
    ROM_set_param(
            &(sys.rom_sups),
            sys.rom_prm_v.num_subdomains,
            sys.rom_prm_p.num_modes + sys.rom_prm_v.num_modes,
            sys.rom_prm_v.rom_epsilon,
            sys.rom_prm_v.solver_type);
    
	ROM_sys_hlpod_fe_set_bc_id(
            (&sys.bc),
            sys.fe.total_num_nodes,
            4,
            &(sys.rom_sups.rom_bc));

    const char* parted_file_name;
    parted_file_name = ROM_std_hlpod_get_parted_file_name(sys.rom_prm_v.solver_type);

    const char* metagraph_name;
    metagraph_name = ROM_std_hlpod_get_metagraph_name(sys.rom_prm_v.solver_type);

	ROM_std_hlpod_pre(
            &(sys.rom_v),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            3,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);
    
	ROM_std_hlpod_pre(
            &(sys.rom_p),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            1,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);

	ROM_std_hlpod_pre(
            &(sys.rom_sups),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            4,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);

    /******************/

    /*for offline*/
    ROM_offline_read_calc_conditions(&(sys.vals), sys.cond.directory);
	ROM_offline_set_reynolds_num_cases(&(sys.vals), sys.cond.directory);

	ROM_std_hlpod_offline_memory_allocation_snapmat(
			&(sys.rom_v),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            sys.vals.finish_time,
            sys.vals.dt,
            sys.vals.snapshot_interval,
            sys.vals.num_cases,
			3);

	ROM_std_hlpod_offline_memory_allocation_snapmat(
			&(sys.rom_p),
			sys.fe.total_num_nodes,
            sys.mono_com.n_internal_vertex,
            sys.vals.finish_time,
            sys.vals.dt,
            sys.vals.snapshot_interval,
            sys.vals.num_cases,
			1);
    
    /******************/

	/**********************************************/

	/****************** solver ********************/
	double t = 0.0;
	int step = 0;
	int file_num = 0;
	int count = 0;  //for ROM

	for(int i = 0; i < sys.vals.num_cases; i++){
		ROM_offline_set_reynolds_number(&(sys.vals), i);
		
		t = 0.0; step = 0; file_num = 0;
		
		while (t < sys.vals.finish_time) {
			t += sys.vals.dt;
			step += 1;

			printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);

			if(step > (sys.rom_v.hlpod_vals.num_snapshot * sys.vals.snapshot_interval)){
				break;
			}
			solver_fom_collect_snapmat(sys, t, count);

			count ++;

			if(step%sys.vals.output_interval == 0) {
				output_files(&sys, file_num, t);
				file_num += 1;
                
                // cavity
    			//output_result_file_cavity_center_vx(&(sys.vals), sys.cond.directory);
			}
		}
	}
    /**********************************************/

	ROM_std_hlpod_set_pod_modes_diag(
		&(sys.rom_v),
		&(sys.rom_p),
		&(sys.rom_sups),
		sys.fe.total_num_nodes,
		sys.mono_com.n_internal_vertex,
		3,
		1,
		"pod_modes_v",
		"pod_modes_p",
		sys.cond.directory);

    /*for writing vtk*/
    ROM_std_hlpod_read_pod_modes_diag(
		&(sys.rom_v),
		&(sys.rom_p),
		&(sys.rom_sups),
		sys.fe.total_num_nodes,
		sys.mono_com.n_internal_vertex,
		3,
		1,
		"pod_modes_v",
		"pod_modes_p",
		sys.cond.directory);

	ROM_sys_hlpod_fe_write_pod_modes_vtk_diag(
		&(sys.fe),
		&(sys.rom_sups),
		sys.fe.total_num_nodes,
		10,
		10,
		3,
		1,
		"pod_modes_vtk/pod_modes_v.vtk",
		"pod_modes_vtk/pod_modes_p.vtk",
		sys.cond.directory);
    /***************/

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

