
#include "core_ROM.h"
#include "core_HROM.h"
#include "set_matvec_NNLS.h"

static const char* OPTION_NUM_MODES = "-nm";
static const char* OPTION_NUM_1STDD = "-nd";
static const char* OPTION_PADAPTIVE = "-pa";
static const char* OPTION_SOLVER_TYPE = "-st";

static const char* INPUT_DIRECTORYNAME_METAGRAPH = "metagraph_parted.0/";
static const char* INPUT_FILENAME_METAGRAPH = "metagraph.dat";

static const char* INPUT_FILENAME_COND    = "cond.dat";


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

	sys.cond.directory = BBFE_convdiff_get_directory_name(argc, argv, CODENAME);
	read_calc_conditions(&(sys.vals), sys.cond.directory);

	BBFE_convdiff_pre(
			&(sys.fe), &(sys.basis), (&sys.bc), (&sys.monolis), (&sys.monolis_com),
			argc, argv, sys.cond.directory,
			sys.vals.num_ip_each_axis,
			true);

	memory_allocation_nodal_values(
			&(sys.vals),
			sys.fe.total_num_nodes);
	manusol_set_init_value(&(sys.fe), sys.vals.T);

	FILE* fp;
	fp = BBFE_sys_write_fopen(fp, "l2_error.txt", sys.cond.directory);
	fclose(fp);

	BBFE_elemmat_set_Jacobi_mat(
			&(sys.fe),
			&(sys.basis));
	BBFE_elemmat_set_shapefunc_derivative(
			&(sys.fe),
			&(sys.basis));

	monolis_initialize(&(sys.monolis0));
    monolis_initialize(&(sys.monolis));

	monolis_com_initialize_by_parted_files(
			&(sys.monolis_com),
			monolis_mpi_get_global_comm(),
			MONOLIS_DEFAULT_TOP_DIR,
			MONOLIS_DEFAULT_PART_DIR,
			"node.dat");

	monolis_get_nonzero_pattern_by_simple_mesh_R(
			&(sys.monolis0),
			sys.fe.total_num_nodes,
			sys.fe.local_num_nodes,
			1,
			sys.fe.total_num_elems,
			sys.fe.conn);

	set_element_mat(
			&(sys.monolis0),
			&(sys.fe),
			&(sys.basis),
			&(sys.vals));

    /*for ROM *****************************************/
    
    /*for ROM input data*/
	ROM_read_args(argc, argv, &(sys.rom_prm));

    ROM_set_param(
            &(sys.rom),
            sys.rom_prm.num_subdomains,
            sys.rom_prm.num_modes,
            sys.rom_prm.rom_epsilon,
            sys.rom_prm.solver_type);
    
	ROM_sys_hlpod_fe_set_bc_id(
            (&sys.bc),
            sys.fe.total_num_nodes,
            1,
            &(sys.rom.rom_bc));

    const char* parted_file_name;
    parted_file_name = ROM_std_hlpod_get_parted_file_name(sys.rom_prm.solver_type);

    const char* metagraph_name;
    metagraph_name = ROM_std_hlpod_get_metagraph_name(sys.rom_prm.solver_type);

	ROM_std_hlpod_pre(
            &(sys.rom),
			sys.fe.total_num_nodes,
            sys.monolis_com.n_internal_vertex,
            1,
            metagraph_name,
            parted_file_name,
			sys.cond.directory);
    
    /******************/

    
    /*for offline phase*/
    ROM_offline_read_calc_conditions(&(sys.vals), sys.cond.directory);

	ROM_std_hlpod_offline_set_num_snapmat(
			&(sys.rom),
            sys.vals.finish_time,
            sys.vals.dt,
            sys.vals.snapshot_interval,
            1);
    /******************/

	/*for online phase*/
    read_calc_conditions(&(sys.vals_rom), sys.cond.directory);                  //set vals

    memory_allocation_nodal_values(&(sys.vals_rom), sys.fe.total_num_nodes);	//set vals

	ROM_std_hlpod_online_memory_allocation_ansvec(&(sys.rom.hlpod_vals), sys.fe.total_num_nodes, 1);

    ROM_online_read_calc_conditions(&(sys.vals), sys.cond.directory);

	monolis_initialize(&(sys.monolis_rom0));
    monolis_initialize(&(sys.monolis_rom));
	monolis_com_initialize_by_self(&(sys.mono_com_rom));
	monolis_com_initialize_by_self(&(sys.mono_com_rom_solv));

	ROM_std_hlpod_set_monolis_comm(
			&(sys.monolis_com),
			&(sys.mono_com_rom),
			&(sys.mono_com_rom_solv),
			INPUT_DIRECTORYNAME_METAGRAPH,
            INPUT_FILENAME_METAGRAPH,
            sys.rom_prm.solver_type,
			sys.cond.directory);

	if(monolis_mpi_get_global_my_rank() == 0){
		fp = ROM_BB_write_fopen(fp, "l2_error_rom.txt", sys.cond.directory);
		fclose(fp);
	}

	if(monolis_mpi_get_global_my_rank() == 0){
        ROM_std_hlpod_write_solver_prm_fopen("fem_solver_prm", sys.cond.directory);
		ROM_std_hlpod_write_solver_prm_fopen("pod_solver_prm", sys.cond.directory);
	}

    ROM_std_hlpod_read_pod_modes(
            &(sys.rom),
            sys.fe.total_num_nodes,
            sys.monolis_com.n_internal_vertex,
            1,
            "pod_modes",
            sys.cond.directory);

	ROM_std_hlpod_online_pre(
            &(sys.monolis_rom0),
            &(sys.monolis_com),
            &(sys.mono_com_rom),
            &(sys.mono_com_rom_solv),
            &(sys.rom),
            sys.fe.total_num_nodes,
            metagraph_name,
            sys.cond.directory);
    
    monolis_copy_mat_nonzero_pattern_R(&(sys.monolis_rom0), &(sys.monolis_rom));
    monolis_com_initialize_by_self(&(sys.mono_com0));
    /******************/

    /**************************************************/

    /*for Hyper-reduction*/
    HROM_pre(&sys, &(sys.rom), &(sys.hrom));
    HROM_memory_allocation(&sys, &(sys.rom), &(sys.hrom));
    /*********************/


    monolis_copy_mat_R(&(sys.monolis0), &(sys.monolis));
    ROM_BB_vec_copy(sys.vals.T, sys.vals_rom.T, sys.fe.total_num_nodes);    

    int file_num = 0;
    int step = 0;
    double t = 0;
	while (t < sys.vals.finish_time) {
		t += sys.vals.dt;
		step += 1;

		printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);
        double calctime_fem_t1 = monolis_get_time();
        solver_fom(sys, t, step);	
        double calctime_fem_t2 = monolis_get_time();
		/**********************************************/

        /****************** ROM solver ****************/
        double calctime_rom_t1 = monolis_get_time();
        solver_rom(&(sys), step, t);        
        double calctime_rom_t2 = monolis_get_time();
		/**********************************************/

        HROM_set_matvec(&(sys),&(sys.rom),&(sys.hrom),step,t);

        if(step%sys.vals.output_interval == 0) {
			ROM_output_files(&sys, file_num, t);
                        
			file_num += 1;
		}
    }

   HROM_pre_offline2(&sys, &(sys.rom), &(sys.hrom));

    if(monolis_mpi_get_global_my_rank() == 0){
        fp = BBFE_sys_write_fopen(fp, "calctime/time_offline.txt", sys.cond.directory);
        //fprintf(fp, "%e\n", calctime_offline);
        fclose(fp);
    }

	BBFE_convdiff_finalize(&(sys.fe), &(sys.basis), &(sys.bc));

	monolis_finalize(&(sys.monolis));
	monolis_finalize(&(sys.monolis0));

	monolis_global_finalize();

	double t2 = monolis_get_time_global_sync();
	printf("** Total time: %f\n", t2 - t1);

	printf("\n");

	return 0;
}
