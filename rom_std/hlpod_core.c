
#include "hlpod_core.h"

static const char* CODENAME = "ROM";

void ROM_std_hlpod_online_memory_allocation_ansvec(
        HLPOD_VALUES*	hlpod_vals,
        const int		total_num_nodes,
        const int	    dof)
{
    hlpod_vals->sol_vec = BB_std_calloc_1d_double(hlpod_vals->sol_vec, total_num_nodes * dof);
}


void ROM_std_hlpod_offline_set_num_snapmat(
        ROM*            rom,
        const int		total_num_nodes,
        const int		n_internal_vertex,
        const double    finish_time,
        const double    dt,
        const int       snapshot_interval,
        const int       num_case,
        const int	    dof)
{
    double quotient = finish_time / dt / snapshot_interval;
    
    if (fmod(quotient, 1.0) == 0.0) {
    	rom->hlpod_vals.num_snapshot = finish_time / dt / snapshot_interval * num_case;
    	printf("%s: num_snapshot = %d\n", CODENAME, rom->hlpod_vals.num_snapshot);
    }
    else{
        printf("Error: num_snapshot = %d is not integer\n");
        exit(1);
    }
}


void ROM_std_hlpod_offline_memory_allocation_snapmat(
        ROM*            rom,
        const int		total_num_nodes,
        const int		n_internal_vertex,
        const double    finish_time,
        const double    dt,
        const int       snapshot_interval,
        const int       num_case,
        const int	    dof)
{
    double quotient = finish_time / dt / snapshot_interval;
    
    if (fmod(quotient, 1.0) == 0.0) {
    	rom->hlpod_vals.num_snapshot = finish_time / dt / snapshot_interval * num_case;
    	printf("%s: num_snapshot = %d\n", CODENAME, rom->hlpod_vals.num_snapshot);
    }
    else{
        printf("Error: num_snapshot = %d is not integer\n");
        exit(1);
    }

    if (monolis_mpi_get_global_comm_size() == 1){
        rom->hlpod_mat.snapmat = BB_std_calloc_2d_double(rom->hlpod_mat.snapmat, total_num_nodes * dof, rom->hlpod_vals.num_snapshot);
    }
    else{
        rom->hlpod_mat.snapmat = BB_std_calloc_2d_double(rom->hlpod_mat.snapmat, n_internal_vertex * dof, rom->hlpod_vals.num_snapshot);
    }
}



void ROM_std_hlpod_set_pod_modes_diag(
        ROM* 		rom_v,
        ROM* 		rom_p,
        ROM* 		rom_sups,
        const int 	total_num_nodes,
        const int 	n_internal_vertex,
        const int 	ndof1,
        const int 	ndof2,
        const char* label1,
        const char* label2,
        const char* directory)
{
    if(monolis_mpi_get_global_comm_size() == 1){

        if(rom_sups->hlpod_vals.num_1st_subdomains == 0){
            printf("\nError : num_1st_subdomains is not set\n");
            exit(1);
        }
        else if(rom_sups->hlpod_vals.num_1st_subdomains==1){
            ROM_std_hlpod_set_podmodes(
                    &(rom_v->hlpod_vals),
                    &(rom_v->hlpod_mat),
                    total_num_nodes,
                    rom_v->hlpod_vals.num_modes_pre,
                    rom_v->hlpod_vals.num_snapshot,
                    rom_v->hlpod_vals.rom_epsilon,
                    ndof1,
                    label1,
                    directory);

            ROM_std_hlpod_set_podmodes(
                    &(rom_p->hlpod_vals),
                    &(rom_p->hlpod_mat),
                    total_num_nodes,
                    rom_p->hlpod_vals.num_modes_pre,
                    rom_p->hlpod_vals.num_snapshot,
                    rom_p->hlpod_vals.rom_epsilon,
                    ndof2,
                    label2,
                    directory);

            ROM_std_hlpod_free_podmodes(
                    &(rom_v->hlpod_mat),
                    total_num_nodes,
                    rom_v->hlpod_vals.num_modes_pre,
                    ndof1);
            
            ROM_std_hlpod_free_podmodes(
                    &(rom_p->hlpod_mat),
                    total_num_nodes,
                    rom_p->hlpod_vals.num_modes_pre,
                    ndof2);
        }
        else{
            ROM_std_hlpod_set_podmodes_local(
                    &(rom_v->hlpod_vals),
                    &(rom_v->hlpod_mat),
                    total_num_nodes,
                    rom_v->hlpod_vals.num_modes_pre,
                    rom_v->hlpod_vals.num_snapshot,
                    rom_v->hlpod_vals.num_1st_subdomains,
                    rom_v->hlpod_vals.rom_epsilon,
                    3,
                    label1,
                    directory);

            ROM_std_hlpod_set_podmodes_local(
                    &(rom_p->hlpod_vals),
                    &(rom_p->hlpod_mat),
                    total_num_nodes,
                    rom_p->hlpod_vals.num_modes_pre,
                    rom_p->hlpod_vals.num_snapshot,
                    rom_p->hlpod_vals.num_1st_subdomains,
                    rom_p->hlpod_vals.rom_epsilon,
                    1,
                    label2,
                    directory);

            ROM_std_hlpod_free_local_podmodes(
                    &(rom_v->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_v->hlpod_vals.num_modes_pre,
                    rom_v->hlpod_vals.num_2nd_subdomains,
                    ndof1);
            
            ROM_std_hlpod_free_local_podmodes(
                    &(rom_p->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_p->hlpod_vals.num_modes_pre,
                    rom_p->hlpod_vals.num_2nd_subdomains,
                    ndof2);
        }
    }		
    else{
        if(rom_sups->hlpod_vals.bool_global_mode==false){

            ROM_std_hlpod_set_podmodes_local_para(
                    &(rom_v->hlpod_vals),
                    &(rom_v->hlpod_mat),
                    &(rom_sups->hlpod_meta),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_v->hlpod_vals.num_modes_pre,
                    rom_v->hlpod_vals.num_snapshot,
                    rom_sups->hlpod_vals.num_2nd_subdomains,
                    rom_v->hlpod_vals.rom_epsilon,
                    ndof1,
                    label1,
                    directory);
            
            ROM_std_hlpod_set_podmodes_local_para(
                    &(rom_p->hlpod_vals),
                    &(rom_p->hlpod_mat),
                    &(rom_sups->hlpod_meta),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_p->hlpod_vals.num_modes_pre,
                    rom_p->hlpod_vals.num_snapshot,
                    rom_sups->hlpod_vals.num_2nd_subdomains,
                    rom_p->hlpod_vals.rom_epsilon,
                    ndof2,
                    label2,
                    directory);

            ROM_std_hlpod_free_local_podmodes_para(
                    &(rom_v->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_v->hlpod_vals.num_modes_pre,
                    rom_v->hlpod_vals.num_2nd_subdomains,
                    ndof1);
            
            ROM_std_hlpod_free_local_podmodes_para(
                    &(rom_p->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_p->hlpod_vals.num_modes_pre,
                    rom_p->hlpod_vals.num_2nd_subdomains,
                    ndof2);

        }
        else{

            ROM_std_hlpod_set_podmodes_global_para(
                    &(rom_v->hlpod_vals),
                    &(rom_v->hlpod_mat),
                    rom_v->hlpod_mat.snapmat,
                    total_num_nodes,
                    n_internal_vertex,
                    rom_v->hlpod_vals.num_modes_pre,
                    rom_v->hlpod_vals.num_snapshot,
                    rom_v->hlpod_vals.rom_epsilon,
                    ndof1,
                    label1,
                    directory);

            ROM_std_hlpod_set_podmodes_global_para(
                    &(rom_p->hlpod_vals),
                    &(rom_p->hlpod_mat),
                    rom_p->hlpod_mat.snapmat,
                    total_num_nodes,
                    n_internal_vertex,
                    rom_p->hlpod_vals.num_modes_pre,
                    rom_p->hlpod_vals.num_snapshot,
                    rom_p->hlpod_vals.rom_epsilon,
                    ndof2,
                    label2,
                    directory);

            ROM_std_hlpod_free_global_podmodes(
                    &(rom_v->hlpod_mat),
                    total_num_nodes,
                    rom_v->hlpod_vals.num_modes_pre,
                    ndof1);

            ROM_std_hlpod_free_global_podmodes(
                    &(rom_p->hlpod_mat),
                    total_num_nodes,
                    rom_p->hlpod_vals.num_modes_pre,
                    ndof2);
        }
    }
}

void ROM_std_hlpod_read_pod_modes_diag(
        ROM* 		rom_v,
        ROM* 		rom_p,
        ROM* 		rom_sups,
        const int 	total_num_nodes,
        const int 	n_internal_vertex,
        const int 	ndof1,
        const int 	ndof2,
        const char* label1,
        const char* label2,
        const char* directory)
{
    if(monolis_mpi_get_global_comm_size() == 1){

        if(rom_sups->hlpod_vals.num_1st_subdomains == 0){
            printf("\nError : num_1st_subdomains is not set\n");
            exit(1);
        }
        else if(rom_sups->hlpod_vals.num_1st_subdomains==1){
            ROM_std_hlpod_read_podmodes(
                    &(rom_v->hlpod_vals),
                    &(rom_v->hlpod_mat),
                    total_num_nodes,
                    rom_v->hlpod_vals.num_modes_pre,
                    rom_v->hlpod_vals.num_snapshot,
                    rom_v->hlpod_vals.rom_epsilon,
                    ndof1,
                    label1,
                    directory);

            ROM_std_hlpod_read_podmodes(
                    &(rom_p->hlpod_vals),
                    &(rom_p->hlpod_mat),
                    total_num_nodes,
                    rom_p->hlpod_vals.num_modes_pre,
                    rom_p->hlpod_vals.num_snapshot,
                    rom_p->hlpod_vals.rom_epsilon,
                    ndof2,
                    label2,
                    directory);

            ROM_std_hlpod_set_podmodes_diag(
                    &(rom_sups->hlpod_vals),
                    &(rom_sups->hlpod_mat),
                    rom_v->hlpod_mat.pod_modes,
                    rom_p->hlpod_mat.pod_modes,
                    total_num_nodes,
                    rom_v->hlpod_vals.num_modes,
                    rom_p->hlpod_vals.num_modes,
                    ndof1,
                    ndof2,
                    directory);
            
            rom_sups->hlpod_vals.num_modes_pre = rom_v->hlpod_vals.num_modes_pre + rom_v->hlpod_vals.num_modes_pre;

            ROM_std_hlpod_free_podmodes(
                    &(rom_v->hlpod_mat),
                    total_num_nodes,
                    rom_v->hlpod_vals.num_modes_pre,
                    ndof1);
            
            ROM_std_hlpod_free_podmodes(
                    &(rom_p->hlpod_mat),
                    total_num_nodes,
                    rom_p->hlpod_vals.num_modes_pre,
                    ndof2);
        }
        else{
            ROM_std_hlpod_read_podmodes_local(
                    &(rom_v->hlpod_vals),
                    &(rom_v->hlpod_mat),
                    total_num_nodes,
                    rom_v->hlpod_vals.num_modes_pre,
                    rom_v->hlpod_vals.num_snapshot,
                    rom_v->hlpod_vals.num_1st_subdomains,
                    3,
                    label1,
                    directory);

            ROM_std_hlpod_read_podmodes_local(
                    &(rom_p->hlpod_vals),
                    &(rom_p->hlpod_mat),
                    total_num_nodes,
                    rom_p->hlpod_vals.num_modes_pre,
                    rom_p->hlpod_vals.num_snapshot,
                    rom_p->hlpod_vals.num_1st_subdomains,
                    1,
                    label2,
                    directory);

            ROM_std_hlpod_set_podmodes_local_diag(
                    &(rom_sups->hlpod_vals),
                    &(rom_sups->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_v->hlpod_mat.pod_modes,
                    rom_p->hlpod_mat.pod_modes,
                    rom_v->hlpod_mat.num_modes_internal,
                    rom_p->hlpod_mat.num_modes_internal,
                    rom_p->hlpod_mat.n_internal_vertex_subd,
                    rom_p->hlpod_mat.node_id,
                    rom_sups->hlpod_vals.num_1st_subdomains,
                    rom_v->hlpod_vals.num_modes,
                    rom_p->hlpod_vals.num_modes,
                    ndof1,
                    ndof2,
                    directory);
            
            rom_sups->hlpod_vals.num_modes_pre = rom_v->hlpod_vals.num_modes_pre + rom_v->hlpod_vals.num_modes_pre;

            ROM_std_hlpod_free_local_podmodes(
                    &(rom_v->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_v->hlpod_vals.num_modes_pre,
                    rom_v->hlpod_vals.num_2nd_subdomains,
                    ndof1);
            
            ROM_std_hlpod_free_local_podmodes(
                    &(rom_p->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_p->hlpod_vals.num_modes_pre,
                    rom_p->hlpod_vals.num_2nd_subdomains,
                    ndof2);

        }
    }		
    else{
        if(rom_sups->hlpod_vals.bool_global_mode==false){
            ROM_std_hlpod_read_podmodes_local_para(
                    &(rom_v->hlpod_vals),
                    &(rom_v->hlpod_mat),
                    &(rom_sups->hlpod_meta),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_v->hlpod_vals.num_modes_pre,
                    rom_v->hlpod_vals.num_snapshot,
                    rom_sups->hlpod_vals.num_2nd_subdomains,
                    ndof1,
                    rom_v->hlpod_vals.rom_epsilon,
                    label1,
                    directory);
            
            ROM_std_hlpod_read_podmodes_local_para(
                    &(rom_p->hlpod_vals),
                    &(rom_p->hlpod_mat),
                    &(rom_sups->hlpod_meta),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_p->hlpod_vals.num_modes_pre,
                    rom_p->hlpod_vals.num_snapshot,
                    rom_sups->hlpod_vals.num_2nd_subdomains,
                    ndof2,
                    rom_p->hlpod_vals.rom_epsilon,
                    label2,
                    directory);
            
            ROM_std_hlpod_set_podmodes_local_para_diag(
                    &(rom_sups->hlpod_vals),
                    &(rom_sups->hlpod_mat),
                    &(rom_sups->hlpod_meta),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_v->hlpod_mat.pod_modes,
                    rom_p->hlpod_mat.pod_modes,
                    rom_v->hlpod_mat.num_modes_internal,
                    rom_p->hlpod_mat.num_modes_internal,
                    rom_v->hlpod_vals.num_modes,
                    rom_p->hlpod_vals.num_modes,
                    ndof1,
                    ndof2,
                    rom_sups->hlpod_vals.num_2nd_subdomains,
                    directory);
                    
            
            rom_sups->hlpod_vals.num_modes_pre = rom_v->hlpod_vals.num_modes_pre + rom_v->hlpod_vals.num_modes_pre;

            ROM_std_hlpod_free_local_podmodes_para(
                    &(rom_v->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_v->hlpod_vals.num_modes_pre,
                    rom_v->hlpod_vals.num_2nd_subdomains,
                    ndof1);
            
            ROM_std_hlpod_free_local_podmodes_para(
                    &(rom_p->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_p->hlpod_vals.num_modes_pre,
                    rom_p->hlpod_vals.num_2nd_subdomains,
                    ndof2);
        }
        else{
            ROM_std_hlpod_read_podmodes_global_para(
                    &(rom_v->hlpod_vals),
                    &(rom_v->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_v->hlpod_vals.num_modes_pre,
                    ndof1,
                    label1,
                    directory);

            ROM_std_hlpod_read_podmodes_global_para(
                    &(rom_p->hlpod_vals),
                    &(rom_p->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_p->hlpod_vals.num_modes_pre,
                    ndof2,
                    label2,
                    directory);
            
            ROM_std_hlpod_set_podmodes_global_para_diag(
                    &(rom_sups->hlpod_vals),
                    &(rom_sups->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom_v->hlpod_mat.pod_modes,
                    rom_p->hlpod_mat.pod_modes,
                    rom_v->hlpod_vals.num_modes_pre,
                    rom_p->hlpod_vals.num_modes_pre,
                    ndof1,
                    ndof2,
                    directory);
            
            rom_sups->hlpod_vals.num_modes_pre = rom_v->hlpod_vals.num_modes_pre + rom_v->hlpod_vals.num_modes_pre;

            ROM_std_hlpod_free_global_podmodes(
                    &(rom_v->hlpod_mat),
                    total_num_nodes,
                    rom_v->hlpod_vals.num_modes_pre,
                    ndof1);

            ROM_std_hlpod_free_global_podmodes(
                    &(rom_p->hlpod_mat),
                    total_num_nodes,
                    rom_p->hlpod_vals.num_modes_pre,
                    ndof2);
        }
    }
}

void ROM_std_hlpod_set_pod_modes(
        ROM* 		rom,
        const int 	total_num_nodes,
        const int 	n_internal_vertex,
        const int 	ndof,
        const char* label,
        const char* directory)
{
    if(monolis_mpi_get_global_comm_size() == 1){
        if(rom->hlpod_vals.num_1st_subdomains == 0){
            printf("\nError : num_1st_subdomains is not set\n");
            exit(1);
        }
        else if(rom->hlpod_vals.num_1st_subdomains==1){

        ROM_std_hlpod_set_podmodes(
                &(rom->hlpod_vals),
                &(rom->hlpod_mat),
                total_num_nodes,
                rom->hlpod_vals.num_modes_pre,
                rom->hlpod_vals.num_snapshot,
                rom->hlpod_vals.rom_epsilon,
                ndof,
                label,
                directory);
        }
        else{
            ROM_std_hlpod_set_podmodes_local(
                    &(rom->hlpod_vals),
                    &(rom->hlpod_mat),
                    total_num_nodes,
                    rom->hlpod_vals.num_modes_pre,
                    rom->hlpod_vals.num_snapshot,
                    rom->hlpod_vals.num_1st_subdomains,
                    rom->hlpod_vals.rom_epsilon,
                    ndof,
                    label,
                    directory);
        }
    }		
    else{
        if(rom->hlpod_vals.bool_global_mode==false){

            ROM_std_hlpod_set_podmodes_local_para(
                    &(rom->hlpod_vals),
                    &(rom->hlpod_mat),
                    &(rom->hlpod_meta),
                    total_num_nodes,
                    n_internal_vertex,
                    rom->hlpod_vals.num_modes_pre,
                    rom->hlpod_vals.num_snapshot,
                    rom->hlpod_vals.num_2nd_subdomains,
                    rom->hlpod_vals.rom_epsilon,
                    ndof,
                    label,
                    directory);
        
        }
        else{

            ROM_std_hlpod_set_podmodes_global_para(
                    &(rom->hlpod_vals),
                    &(rom->hlpod_mat),
                    rom->hlpod_mat.snapmat,
                    total_num_nodes,
                    n_internal_vertex,
                    rom->hlpod_vals.num_modes_pre,
                    rom->hlpod_vals.num_snapshot,
                    rom->hlpod_vals.rom_epsilon,
                    ndof,
                    label,
                    directory);

        }
    }
}


void ROM_std_hlpod_read_pod_modes(
        ROM* 		rom,
        const int 	total_num_nodes,
        const int 	n_internal_vertex,
        const int 	dof,
        const char* label,
        const char* directory)
{
    if(monolis_mpi_get_global_comm_size() == 1){

        if(rom->hlpod_vals.num_1st_subdomains == 0){
            printf("\nError : num_1st_subdomains is not set\n");
            exit(1);
        }
        else if(rom->hlpod_vals.num_1st_subdomains==1){
            ROM_std_hlpod_read_podmodes(
                    &(rom->hlpod_vals),
                    &(rom->hlpod_mat),
                    total_num_nodes,
                    rom->hlpod_vals.num_modes_pre,
                    rom->hlpod_vals.num_snapshot,
                    rom->hlpod_vals.rom_epsilon,
                    dof,
                    label,
                    directory);

        }
        else{
            ROM_std_hlpod_free_local_podmodes(
                    &(rom->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom->hlpod_vals.num_modes_pre,
                    rom->hlpod_vals.num_2nd_subdomains,
                    dof);

            ROM_std_hlpod_read_podmodes_local(
                    &(rom->hlpod_vals),
                    &(rom->hlpod_mat),
                    total_num_nodes,
                    rom->hlpod_vals.num_modes_pre,
                    rom->hlpod_vals.num_snapshot,
                    rom->hlpod_vals.num_1st_subdomains,
                    dof,
                    label,
                    directory);

        }
    }		
    else{
        if(rom->hlpod_vals.bool_global_mode==false){
            ROM_std_hlpod_free_local_podmodes_para(
                    &(rom->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom->hlpod_vals.num_modes_pre,
                    rom->hlpod_vals.num_2nd_subdomains,
                    dof);

            ROM_std_hlpod_read_podmodes_local_para(
                    &(rom->hlpod_vals),
                    &(rom->hlpod_mat),
                    &(rom->hlpod_meta),
                    total_num_nodes,
                    n_internal_vertex,
                    rom->hlpod_vals.num_modes_pre,
                    rom->hlpod_vals.num_snapshot,
                    rom->hlpod_vals.num_2nd_subdomains,
                    dof,
                    rom->hlpod_vals.rom_epsilon,
                    label,
                    directory);
        }
        else{
            ROM_std_hlpod_read_podmodes_global_para(
                    &(rom->hlpod_vals),
                    &(rom->hlpod_mat),
                    total_num_nodes,
                    n_internal_vertex,
                    rom->hlpod_vals.num_modes_pre,
                    dof,
                    label,
                    directory);

        }
    }
}


void ROM_std_hlpod_pre(
        ROM*         rom,
        const int    total_num_nodes,
        const int    n_internal_vertex,
        const int    dof,
        const char*  metagraph,
        const char*  label_pod_subd,
        const char*	 directory)
{            
    if(monolis_mpi_get_global_comm_size() == 1){
        if(rom->hlpod_vals.num_1st_subdomains==0){
            printf("\nError : num_1st_subdomains is not set\n");
            exit(1);
        }
        else if(rom->hlpod_vals.num_1st_subdomains==1){

        }
        else{
            ROM_std_hlpod_read_node_id_pod_subd(
                    &(rom->hlpod_mat),
                    &(rom->hlpod_meta),
                    total_num_nodes,
                    rom->hlpod_vals.num_2nd_subdomains,
                    label_pod_subd,
                    directory);
        }
    }
    else{
        if(rom->hlpod_vals.bool_global_mode==false){
            ROM_std_hlpod_get_subdomain_id(
                    &(rom->hlpod_vals),
                    &(rom->hlpod_meta),
                    monolis_mpi_get_global_my_rank(),
                    metagraph,
                    directory);

            ROM_std_hlpod_read_node_id_pod_subd_para(
                    &(rom->hlpod_mat),
                    &(rom->hlpod_meta),
                    total_num_nodes,
                    n_internal_vertex,
                    rom->hlpod_vals.num_2nd_subdomains,
                    label_pod_subd,
                    directory);
        }
    }
}

void ROM_std_hlpod_set_monolis_comm(
        MONOLIS_COM* monolis_com,
        MONOLIS_COM* mono_com_rom,
        MONOLIS_COM* mono_com_rom_solv,
        const char*	 metagraph_parted0,
        const char*  metagraph,
        const int    solver_type,
        const char*	 directory)
{
    if(monolis_mpi_get_global_comm_size() == 1){
        monolis_com_initialize_by_self(mono_com_rom_solv);
    }
    else{
        if(solver_type == 3){
            ROM_std_hlpod_set_comm_table_para_subd(
                    monolis_com,
                    mono_com_rom,
                    monolis_com->recv_n_neib);
            
            monolis_com_initialize_by_parted_files(
                    mono_com_rom_solv,
                    monolis_mpi_get_global_comm(),
                    directory,
                    metagraph_parted0,
                    metagraph);
        }
    }
}

void ROM_std_hlpod_read_metagraph(
        MONOLIS*     monolis_rom0,
        MONOLIS_COM* mono_com_rom_solv,
        ROM*		 rom,
        const char*  metagraph,
        const char*	 directory)
{
    if(monolis_mpi_get_global_comm_size() == 1){
        if(rom->hlpod_vals.num_1st_subdomains==0){
            printf("\nError : num_1st_subdomains is not set\n");
            exit(1);
        }
        else if(rom->hlpod_vals.num_1st_subdomains==1){
            monolis_com_initialize_by_self(mono_com_rom_solv);

            ROM_std_hlpod_set_nonzero_pattern(
                    monolis_rom0,
                    &(rom->hlpod_mat),
                    rom->hlpod_vals.num_modes_pre);
        }
        else{
            monolis_com_initialize_by_self(mono_com_rom_solv);

            ROM_std_hlpod_set_nonzero_pattern_bcsr(
                    monolis_rom0,
                    &(rom->hlpod_mat),
                    &(rom->hlpod_meta),
                    rom->hlpod_vals.num_modes_pre,
                    metagraph,
                    directory);
        }
    }
    else{
        if(rom->hlpod_vals.bool_global_mode==false){
            ROM_std_hlpod_set_nonzero_pattern_bcsr_para(
                    monolis_rom0,
                    mono_com_rom_solv,
                    &(rom->hlpod_mat),
                    &(rom->hlpod_meta),
                    metagraph,
                    directory);

            ROM_std_hlpod_get_meta_neib(
                    mono_com_rom_solv,
                    &(rom->hlpod_meta),
                    metagraph,
                    directory);

        }
        else{
        }
    }

}


void ROM_std_hlpod_pre_lpod_para(
        MONOLIS*     monolis_rom0,
        MONOLIS_COM* monolis_com,
        MONOLIS_COM* mono_com_rom,
        MONOLIS_COM* mono_com_rom_solv,
        ROM*		 rom,
        const char*	 directory)
{
    /*
    ROM_std_hlpod_get_neib_vec_save_memory(
            monolis_com,
            &(rom->hlpod_vals),
            rom->hlpod_vals.num_modes);
    */

    ROM_std_hlpod_get_neib_vec(
            monolis_com,
            &(rom->hlpod_vals),
            &(rom->hlpod_mat),
            rom->hlpod_vals.num_modes);

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


void ROM_std_hlpod_online_pre(
        MONOLIS*     monolis_rom0,
        MONOLIS_COM* mono_com,
        MONOLIS_COM* mono_com_rom,
        MONOLIS_COM* mono_com_rom_solv,
        ROM* 		 rom_sups,
        const int 	 total_num_nodes,
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

            ROM_std_hlpod_pre_lpod_para(
                    monolis_rom0,
                    mono_com,
                    mono_com_rom,
                    mono_com_rom_solv,
                    rom_sups,
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


void ROM_std_hlpod_calc_reduced_mat(
        MONOLIS*     monolis,
        MONOLIS*     monolis_rom,
        MONOLIS_COM* monolis_com,
        MONOLIS_COM* mono_com0,
        MONOLIS_COM* mono_com_rom,
        ROM*		 rom,
        const int	 total_num_nodes,
        const int	 dof)
{
    if(monolis_mpi_get_global_comm_size() == 1){     
        
        if(rom->hlpod_vals.num_1st_subdomains==0){
            printf("\nError : num_1st_subdomains is not set\n");
        }
        else if(rom->hlpod_vals.num_1st_subdomains==1){
            ROM_std_hlpod_calc_reduced_mat_seq(
                    monolis,
                    monolis_com,
                    &(rom->hlpod_mat),
                    total_num_nodes,
                    rom->hlpod_vals.num_modes,
                    dof);

           int connectivity[] = {0};
            monolis_add_matrix_to_sparse_matrix_R(
                    monolis_rom,					    
                    1,							    
                    connectivity,				
                    rom->hlpod_mat.VTKV);				
        }
        else{
            ROM_std_hlpod_calc_reduced_mat_seq_block(
                    monolis,
                    monolis_com,
                    &(rom->hlpod_mat),
                    total_num_nodes,
                    rom->hlpod_vals.num_modes_pre,
                    rom->hlpod_vals.num_2nd_subdomains,
                    dof);
            
            ROM_std_hlpod_set_reduced_mat(
                    monolis_rom,
                    &(rom->hlpod_mat),
                    &(rom->hlpod_meta),
                    total_num_nodes,
                    rom->hlpod_vals.num_modes_pre,
                    rom->hlpod_vals.num_2nd_subdomains,
                    dof);

        }
    }
    else{

        ROM_std_hlpod_calc_reduced_mat_save_memory(
                monolis,
                monolis_com,
                mono_com0,
                &(rom->hlpod_vals),
                &(rom->hlpod_mat),
                total_num_nodes,
                rom->hlpod_vals.n_neib_vec,
                rom->hlpod_vals.num_2nd_subdomains,
                rom->hlpod_vals.num_modes,
                dof);

        ROM_std_hlpod_set_reduced_mat_para(
                monolis_rom,
                mono_com_rom,
                &(rom->hlpod_vals),
                &(rom->hlpod_mat),
                &(rom->hlpod_meta),
                rom->hlpod_vals.num_modes_pre,
                rom->hlpod_vals.num_modes,
                rom->hlpod_vals.num_2nd_subdomains);
//        double t2 = monolis_get_time_global_sync();
//    exit(1);

    }
}



void ROM_std_hlpod_solve_ROM(
    MONOLIS*     monolis,
    MONOLIS*     monolis_rom,
    MONOLIS_COM* mono_com_rom,
    ROM*		 rom,
    const int	 total_num_nodes,
    const int	 dof,
    const int 	 mat_max_iter,
    const double mat_epsilon,
    const int 	 label_solver,
    const int	 label_prec)
{
    if(monolis_mpi_get_global_comm_size() == 1){
        ROM_std_hlpod_calc_reduced_rhs(
                monolis,
                &(rom->hlpod_mat),
                rom->hlpod_vals.num_modes_pre,
                rom->hlpod_vals.num_2nd_subdomains,
                dof);

        ROM_std_hlpod_reduced_rhs_to_monollis(
                monolis_rom,
                &(rom->hlpod_mat),
                rom->hlpod_vals.num_2nd_subdomains,
                rom->hlpod_vals.num_modes_pre);

        ROM_monowrap_solve(
                monolis_rom,
                mono_com_rom,
                rom->hlpod_mat.mode_coef,
                label_solver,
                label_prec,
                mat_max_iter,
                mat_epsilon);

        ROM_std_hlpod_calc_sol(
                &(rom->hlpod_vals),
                &(rom->hlpod_mat), 
                total_num_nodes,
                rom->hlpod_vals.num_modes_pre,
                rom->hlpod_vals.num_2nd_subdomains,
                dof);

    }
    else{
//        double t1 = monolis_get_time_global_sync();
//        exit(1);
        ROM_std_hlpod_calc_reduced_rhs(
                monolis,
                &(rom->hlpod_mat),
                rom->hlpod_vals.num_modes_pre,
                rom->hlpod_vals.num_2nd_subdomains,
                dof);

        ROM_std_hlpod_reduced_rhs_to_monollis(
                monolis_rom,
                &(rom->hlpod_mat),
                rom->hlpod_vals.num_2nd_subdomains,
                rom->hlpod_vals.num_modes_pre);

//double t = monolis_get_time_global_sync();
//exit(1);

        ROM_monowrap_solve(
                monolis_rom,
                mono_com_rom,
                rom->hlpod_mat.mode_coef,
                label_solver,
                label_prec,
                mat_max_iter,
                mat_epsilon);

        ROM_std_hlpod_calc_sol(
                &(rom->hlpod_vals),
                &(rom->hlpod_mat),
                total_num_nodes,
                rom->hlpod_vals.num_modes_pre,
                rom->hlpod_vals.num_2nd_subdomains,
                dof);

    }
}

