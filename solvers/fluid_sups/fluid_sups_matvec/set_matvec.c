

#include "set_matvec.h"

/*for global mode + parallel computation*/
void set_reduced_mat_global_para(
		MONOLIS*        monolis,
		MONOLIS_COM*	monolis_com,
		BBFE_DATA*      fe,
		BBFE_BASIS*     basis,
		VALUES*         vals,
        double**        mat,
        double**        pod_modes,
        const int 		num_modes)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double*** val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_3d_double(val_ip     , 4 , 4, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	double A[4][4];

    bool* is_internal_elem;
	is_internal_elem = BB_std_calloc_1d_bool(is_internal_elem , fe->total_num_elems);

	monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
		fe->local_num_nodes, fe->conn, is_internal_elem);

	for(int e=0; e<(fe->total_num_elems); e++) {
		if (!is_internal_elem[e]) continue;

		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {

					double tau = BBFE_elemmat_fluid_sups_coef(
							vals->density, vals->viscosity, v_ip[p], h_e, vals->dt);

					BBFE_elemmat_fluid_sups_mat(
							A, basis->N[p][i], basis->N[p][j], 
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], 
							v_ip[p], vals->density, vals->viscosity, tau, vals->dt);

					for(int a=0; a<4; a++){
						for(int b=0; b<4; b++) {
							val_ip[a][b][p] = A[a][b];
							A[a][b] = 0.0;
						}
					}
				}

				for(int a=0; a<4; a++){
					for(int b=0; b<4; b++) {
						double integ_val = BBFE_std_integ_calc(
								np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

                        set_element_reduced_mat(
                                mat,
                                pod_modes,
                                fe->conn[e][i]*4 + a,
                                fe->conn[e][j]*4 + b,
                                integ_val,	
                                num_modes);
                    }
				}
			}
		}
	}

	BB_std_free_3d_double(val_ip, 4 , 4, np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);
}


void set_D_bc_global_para(
		MONOLIS*        monolis,
        MONOLIS_COM*    monolis_com,
		BBFE_DATA*      fe,
		BBFE_BASIS*     basis,
		VALUES*         vals,
    	BBFE_BC*     	bc,
        double*         rhs,
        double**        pod_modes,
        const int 		num_modes)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double*** val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_3d_double(val_ip     , 4 , 4, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	double A[4][4];

    bool* is_internal_elem;
	is_internal_elem = BB_std_calloc_1d_bool(is_internal_elem , fe->total_num_elems);

	monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
		fe->local_num_nodes, fe->conn, is_internal_elem);

	for(int e=0; e<(fe->total_num_elems); e++) {
		if (!is_internal_elem[e]) continue;

		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
		}

		for(int i=0; i<nl; i++) {
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {

					double tau = BBFE_elemmat_fluid_sups_coef(
							vals->density, vals->viscosity, v_ip[p], h_e, vals->dt);

					BBFE_elemmat_fluid_sups_mat(
							A, basis->N[p][i], basis->N[p][j], 
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], 
							v_ip[p], vals->density, vals->viscosity, tau, vals->dt);

					for(int a=0; a<4; a++){
						for(int b=0; b<4; b++) {
							val_ip[a][b][p] = A[a][b];
							A[a][b] = 0.0;
						}
					}
				}

				for(int a=0; a<4; a++){
					for(int b=0; b<4; b++) {
						double integ_val = BBFE_std_integ_calc(
								np, val_ip[a][b], basis->integ_weight, Jacobian_ip);

							set_element_reduced_rhs_Dbc(
                                bc,
                                rhs,
                                pod_modes,
                                fe->conn[e][i]*4 + a,
                                fe->conn[e][j]*4 + b,
                                integ_val,
                                num_modes);

					}
				}
			}
		}
	}


	BB_std_free_3d_double(val_ip     , 4 , 4, np);
	BB_std_free_1d_double(Jacobian_ip, np);

	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);
}

void set_reduced_vec_global_para(
		MONOLIS*        monolis,        
		MONOLIS_COM*	monolis_com,
		BBFE_DATA*      fe,
		BBFE_BASIS*     basis,
		VALUES*         vals,
        double*         rhs,
        double**        pod_modes,
        const int		num_modes)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double** val_ip;
	double*  Jacobian_ip;
	val_ip      = BB_std_calloc_2d_double(val_ip, 4, np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_v;
	local_v = BB_std_calloc_2d_double(local_v, nl, 3);

	double** v_ip; 
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);

	bool* is_internal_elem;
	is_internal_elem = BB_std_calloc_1d_bool(is_internal_elem , fe->total_num_elems);

	monolis_get_bool_list_of_internal_simple_mesh(monolis_com, fe->total_num_nodes, fe->total_num_elems,
		fe->local_num_nodes, fe->conn, is_internal_elem);

	for(int e=0; e<(fe->total_num_elems); e++) {
        if (!is_internal_elem[e]) continue;
		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(
				np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_v, fe, vals->v, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(v_ip[p], nl, local_v, basis->N[p]);
		}

		for(int i=0; i<nl; i++) {
			double integ_val[4];

			for(int p=0; p<np; p++) {
				double tau = BBFE_elemmat_fluid_sups_coef(
						vals->density, vals->viscosity, v_ip[p], h_e, vals->dt);

				double vec[4];
				BBFE_elemmat_fluid_sups_vec(
						vec, basis->N[p][i], fe->geo[e][p].grad_N[i],
						v_ip[p], vals->density, tau, vals->dt);

				for(int d=0; d<4; d++) {
					val_ip[d][p] = vec[d];
				}
			}

			for(int d=0; d<4; d++) {
				integ_val[d] = BBFE_std_integ_calc(
						np, val_ip[d], basis->integ_weight, Jacobian_ip);

                set_element_reduced_rhs(
                        rhs,
                        pod_modes,
                        4 * fe->conn[e][i] + d,
                        integ_val[d],
                        num_modes);
			}
		}
	}

	BB_std_free_2d_double(val_ip, 4, np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);
}


void allreduce_global_para(
    MONOLIS_COM*	monolis_com,
    double**        mat,
    double*         rhs,
    const int		num_modes)
{
    double* vec;
	vec = BB_std_calloc_1d_double(vec, num_modes);

	for(int k2 = 0; k2 < num_modes; k2++){
        vec[k2] = rhs[k2];
	}

	monolis_allreduce_R(
		num_modes,
		vec,
		MONOLIS_MPI_SUM,
		monolis_com->comm);

	for(int k2 = 0; k2 < num_modes; k2++){
        rhs[k2] = vec[k2];
	}


	double* matrix_buffer = BB_std_calloc_1d_double(matrix_buffer, num_modes * num_modes);
	
	for(int k1 = 0; k1 < num_modes; k1++){
		for(int k2 = 0; k2 < num_modes; k2++){
            matrix_buffer[k1 * num_modes + k2] = mat[k1][k2];
		}
	}

	monolis_allreduce_R(
		num_modes * num_modes,
		matrix_buffer,
		MONOLIS_MPI_SUM,
		monolis_com->comm);

	for(int k1 = 0; k1 < num_modes; k1++){
		for(int k2 = 0; k2 < num_modes; k2++){
            mat[k1][k2] = matrix_buffer[k1 * num_modes + k2];
		}
	}
}


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
		const double    dt)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x   = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;  double* a_ip;  double** v_ip;  double* k_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);

    for(int i = 0; i < hlpod_vals->n_neib_vec; i++){
        for(int j = 0; j < hlpod_vals->n_neib_vec; j++){
            hlpod_ddhr->reduced_mat[i][j] = 0.0;
        }
        hlpod_ddhr->reduced_RH[i] = 0.0;
    }

	for(int m=0; m < hlpod_ddhr->ovl_num_selected_elems; m++) {
		int e = hlpod_ddhr->ovl_id_selected_elems[m];

		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);

		}

		for(int i=0; i<nl; i++) {       //六面体一次要素は8
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;
/*
					double tau = BBFE_elemmat_convdiff_stab_coef_ns(
							k_ip[p], v_ip[p], a_ip[p], h_e, vals->dt);
*/
/*
					val_ip[p] += BBFE_elemmat_convdiff_mat_conv(
							basis->N[p][i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p]);
*/
					val_ip[p] -= BBFE_elemmat_convdiff_mat_diff(
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], k_ip[p]);
/*
					val_ip[p] += BBFE_elemmat_convdiff_mat_stab_conv(
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p], tau);
*/
					//val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_mat_mass(
                    val_ip[p] += 1.0/ dt  * BBFE_elemmat_convdiff_mat_mass(
								basis->N[p][i], basis->N[p][j], a_ip[p]);
/*
					val_ip[p] += 1.0/(vals->dt) *  BBFE_elemmat_convdiff_mat_stab_mass(
								fe->geo[e][p].grad_N[i], basis->N[p][j], a_ip[p], v_ip[p], tau);
*/
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);
                
                //基底本数ループ
                int index_i = fe->conn[e][i];
                int index_j = fe->conn[e][j];

				if( bc->D_bc_exists[index_j]){
				}
				else{
					int subdomain_id_i = hlpod_mat->subdomain_id_in_nodes[index_i];
					int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i];
					int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i + 1];

					int subdomain_id_j = hlpod_mat->subdomain_id_in_nodes[index_j];
					int JS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j];
					int JE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j + 1];

					int subdomain_id = hlpod_mat->subdomain_id_in_nodes_2nddd[index_j]; 

				if(subdomain_id_j < num_subdomains){
	                for(int k1 = IS; k1 < IE; k1++){
                        //for(int k2 = 0; k2 < hlpod_vals->n_neib_vec; k2++){
						for(int k2 = JS; k2 < JE; k2++){
							//ver2 変更点
		                    //double val = hlpod_mat->neib_vec[index_i][k1] * integ_val * hlpod_mat->neib_vec[index_j][k2];
							double val = hlpod_mat->pod_basis_hr[index_i][k1] * integ_val * hlpod_mat->pod_basis_hr[index_j][k2];
                            hlpod_ddhr->reduced_mat[k1][k2] += hlpod_ddhr->ovl_elem_weight[m] *  val;
                        }
                    }
				}
/*				else{
				*/

				//printf("subdomain_id_j = %d num_neib_modes_sum = %d JS = %d JE = %d\n", subdomain_id, hlpod_mat->num_neib_modes_sum[subdomain_id-1], JS, JE);
				if(subdomain_id_j >= num_subdomains){
						for(int k1 = IS; k1 < IE; k1++){
							for(int k2 = JS - hlpod_mat->num_neib_modes_sum[subdomain_id-1]; k2 < JE - hlpod_mat->num_neib_modes_sum[subdomain_id-1]; k2++){
	//                        for(int k2 = 0; k2 < hlpod_vals->n_neib_vec; k2++){
								double val = hlpod_mat->pod_basis_hr[index_i][k1] * integ_val * hlpod_mat->pod_basis_hr[index_j][k2] ;
								//hlpod_ddhr->reduced_mat[k1][k2] += hlpod_ddhr->ovl_elem_weight_D_bc[m] *  val;
								hlpod_ddhr->reduced_mat[k1][k2 + hlpod_mat->num_neib_modes_sum[subdomain_id-1]] += hlpod_ddhr->ovl_elem_weight[m] *  val;
							}
						}
				}

				}
			}
		}
	}


    for(int m=0; m < hlpod_ddhr->ovl_num_selected_elems_D_bc; m++) {
		int e = hlpod_ddhr->ovl_id_selected_elems_D_bc[m];

		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);

		}

		for(int i=0; i<nl; i++) {       //六面体一次要素は8
			for(int j=0; j<nl; j++) {

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;
/*
					double tau = BBFE_elemmat_convdiff_stab_coef_ns(
							k_ip[p], v_ip[p], a_ip[p], h_e, vals->dt);
*/
/*
					val_ip[p] += BBFE_elemmat_convdiff_mat_conv(
							basis->N[p][i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p]);
*/
					val_ip[p] -= BBFE_elemmat_convdiff_mat_diff(
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], k_ip[p]);
/*
					val_ip[p] += BBFE_elemmat_convdiff_mat_stab_conv(
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p], tau);
*/
					//val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_mat_mass(
                    val_ip[p] += 1.0/ dt  * BBFE_elemmat_convdiff_mat_mass(
								basis->N[p][i], basis->N[p][j], a_ip[p]);
/*
					val_ip[p] += 1.0/(vals->dt) *  BBFE_elemmat_convdiff_mat_stab_mass(
								fe->geo[e][p].grad_N[i], basis->N[p][j], a_ip[p], v_ip[p], tau);
*/
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

                //基底本数ループ
                int index_i = fe->conn[e][i];
                int index_j = fe->conn[e][j];

                if( bc->D_bc_exists[index_j]) {
                }
                else{
					int subdomain_id_i = hlpod_mat->subdomain_id_in_nodes[index_i];
					int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i];
					int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_i + 1];

					int subdomain_id_j = hlpod_mat->subdomain_id_in_nodes[index_j];
					int JS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j];
					int JE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id_j + 1];

					int subdomain_id = hlpod_mat->subdomain_id_in_nodes_2nddd[index_j]; 

	                //for(int k1 = 0; k1 < hlpod_mat->num_basis; k1++){
					//	for(int k2 = 0; k2 < hlpod_vals->max_num_modes; k2++){
					if(subdomain_id_j < num_subdomains){
						for(int k1 = IS; k1 < IE; k1++){
							for(int k2 = JS; k2 < JE; k2++){
	//                        for(int k2 = 0; k2 < hlpod_vals->n_neib_vec; k2++){
								double val = hlpod_mat->pod_basis_hr[index_i][k1] * integ_val * hlpod_mat->pod_basis_hr[index_j][k2] ;
								//hlpod_ddhr->reduced_mat[k1][k2] += hlpod_ddhr->ovl_elem_weight_D_bc[m] *  val;
								hlpod_ddhr->reduced_mat[k1][k2] += hlpod_ddhr->ovl_elem_weight_D_bc[m] *  val;
							}
						}
					}

					if(subdomain_id_j >= num_subdomains){
						for(int k1 = IS; k1 < IE; k1++){
							for(int k2 = JS - hlpod_mat->num_neib_modes_sum[subdomain_id-1]; k2 < JE - hlpod_mat->num_neib_modes_sum[subdomain_id-1]; k2++){
	//                        for(int k2 = 0; k2 < hlpod_vals->n_neib_vec; k2++){
								double val = hlpod_mat->pod_basis_hr[index_i][k1] * integ_val * hlpod_mat->pod_basis_hr[index_j][k2] ;
								//hlpod_ddhr->reduced_mat[k1][k2] += hlpod_ddhr->ovl_elem_weight_D_bc[m] *  val;
								hlpod_ddhr->reduced_mat[k1][k2 + hlpod_mat->num_neib_modes_sum[subdomain_id-1]] += hlpod_ddhr->ovl_elem_weight_D_bc[m] *  val;
							}
						}
					}
					
                }
			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x,   fe->local_num_nodes, 3);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);

}

void ddhr_lb_set_D_bc_para(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS* 	basis,
    	BBFE_BC*     	bc,
    	HLPOD_MAT*    hlpod_mat,
        HLPOD_DDHR*     hlpod_ddhr,
        const int		num_modes,
		const int 		num_subdomains,
		const double    dt)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;  double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;
	local_x   = BB_std_calloc_2d_double(local_x, nl, 3);

	double** x_ip;  double* a_ip;  double** v_ip;  double* k_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);

    for(int m=0; m < hlpod_ddhr->ovl_num_selected_elems_D_bc; m++) {
        int e = hlpod_ddhr->ovl_id_selected_elems_D_bc[m];

		BBFE_elemmat_set_Jacobian_array(Jacobian_ip, np, e, fe);

		double vol = BBFE_std_integ_calc_volume(np, basis->integ_weight, Jacobian_ip);
		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x, e, 3);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);

		}

//		for(int i=0; i<nl; i++) {       //六面体一次要素は8
		for(int j=0; j<nl; j++) {       //六面体一次要素は8
			int index_j = fe->conn[e][j];
			//int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index_j];
			//if(subdomain_id < num_subdomains){
			
			if( bc->D_bc_exists[index_j]) {
	
//			for(int j=0; j<nl; j++) {
			for(int i=0; i<nl; i++) {
				
				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;
/*
					double tau = BBFE_elemmat_convdiff_stab_coef_ns(
							k_ip[p], v_ip[p], a_ip[p], h_e, vals->dt);
*/
/*
					val_ip[p] += BBFE_elemmat_convdiff_mat_conv(
							basis->N[p][i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p]);
*/
					val_ip[p] -= BBFE_elemmat_convdiff_mat_diff(
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], k_ip[p]);
/*
					val_ip[p] += BBFE_elemmat_convdiff_mat_stab_conv(
							fe->geo[e][p].grad_N[i], fe->geo[e][p].grad_N[j], a_ip[p], v_ip[p], tau);
*/
					//val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_mat_mass(
                    val_ip[p] += 1.0/ dt  * BBFE_elemmat_convdiff_mat_mass(
								basis->N[p][i], basis->N[p][j], a_ip[p]);
/*
					val_ip[p] += 1.0/(vals->dt) *  BBFE_elemmat_convdiff_mat_stab_mass(
								fe->geo[e][p].grad_N[i], basis->N[p][j], a_ip[p], v_ip[p], tau);
*/
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);
                
                //基底本数ループ
                int index_i = fe->conn[e][i];
                int index_j = fe->conn[e][j];

				int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index_i];

				int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id];
				int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id + 1];

                //if( bc->D_bc_exists[index_j]) {
					/*
                    for(int k1 = 0; k1 < hlpod_vals->num_modes; k1++){
						double val = hlpod_mat->neib_vec[index_i][k1] * integ_val * bc->imposed_D_val[index_j];
                        hlpod_ddhr->reduced_RH[k1] += - hlpod_ddhr->ovl_elem_weight_D_bc[m] * val;
                    }
					*/

					for(int k1 = IS; k1 < IE; k1++){
						double val = hlpod_mat->pod_modes[index_i][k1] * integ_val * bc->imposed_D_val[index_j];
                        hlpod_ddhr->reduced_RH[k1] += - hlpod_ddhr->ovl_elem_weight_D_bc[m] * val;
                    }
				//}
			}

			}
		}
	}


	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x,   fe->local_num_nodes, 3);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
}


void ddhr_lb_set_reduced_vec_para(
		MONOLIS*     	monolis,
		BBFE_DATA*     	fe,
		BBFE_BASIS*	 	basis,
        HR_VALUES*      hr_vals,
        HLPOD_VALUES*   hlpod_vals,
        HLPOD_DDHR*     hlpod_ddhr,
    	HLPOD_MAT*    hlpod_mat,
        const int		num_modes,
		const int 		num_subdomains,
        const double    dt,
		double       	t)
{
	int nl = fe->local_num_nodes;
	int np = basis->num_integ_points;

	double* val_ip;
	double* Jacobian_ip;
	val_ip      = BB_std_calloc_1d_double(val_ip     , np);
	Jacobian_ip = BB_std_calloc_1d_double(Jacobian_ip, np);

	double** local_x;  double* local_T;
	local_x = BB_std_calloc_2d_double(local_x, nl, 3);
	local_T = BB_std_calloc_1d_double(local_T, nl);

	double** x_ip;  double* a_ip;  double** v_ip;  double* k_ip;  double* T_ip;  double* f_ip;
	x_ip = BB_std_calloc_2d_double(x_ip, np, 3);
	v_ip = BB_std_calloc_2d_double(v_ip, np, 3);
	k_ip = BB_std_calloc_1d_double(k_ip, np);
	a_ip = BB_std_calloc_1d_double(a_ip, np);
	T_ip = BB_std_calloc_1d_double(T_ip, np);
	f_ip = BB_std_calloc_1d_double(f_ip, np);

    for(int k = 0; k < hlpod_vals->n_neib_vec; k++){
        hlpod_ddhr->reduced_RH[ k ] = 0.0;
    }


    for(int m = 0; m < hlpod_ddhr->ovl_num_selected_elems; m++) {
        int e = hlpod_ddhr->ovl_id_selected_elems[m];

		BBFE_elemmat_set_Jacobian_array(
				Jacobian_ip,
				basis->num_integ_points,
				e,
				fe);

		double vol = BBFE_std_integ_calc_volume(
				basis->num_integ_points,
				basis->integ_weight,
				Jacobian_ip);

		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x,   e, 3);
		BBFE_elemmat_set_local_array_scalar(local_T, fe, hr_vals->sol_vec, e);
/*
        for(int i=0; i<nl; i++) {
            int index = fe->conn[e][i];
			int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index];

            int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id];
            int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id + 1];

            for(int k = IS; k < IE; k++){
                hr_vals->sol_vec[index] = 0.0;
            }
            for(int k = IS; k < IE; k++){
                hr_vals->sol_vec[index] += hlpod_mat->pod_modes[index][k] 
                    * hlpod_mat->mode_coef[k - IS + num_modes*subdomain_id];
            }
        }
        */

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
			T_ip[p] = BBFE_std_mapping_scalar(nl, local_T, basis->N[p]);
			f_ip[p] = manusol_get_source(x_ip[p], t, a_ip[p], v_ip[p], k_ip[p]);
		}

		for(int i=0; i<nl; i++) {
            int index = fe->conn[e][i];
			int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index];

			if(subdomain_id < num_subdomains){
				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;
	/*
					double tau = BBFE_elemmat_convdiff_stab_coef_ns(
							k_ip[p], v_ip[p], a_ip[p], h_e, vals->dt);
	*/
					val_ip[p] += BBFE_elemmat_convdiff_vec_source(
							basis->N[p][i], f_ip[p]);
	/*
					val_ip[p] += BBFE_elemmat_convdiff_vec_stab_source(
							fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], tau, f_ip[p]);
	*/
				//	val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
					val_ip[p] += 1.0/ dt * BBFE_elemmat_convdiff_vec_mass(
							basis->N[p][i], T_ip[p], a_ip[p]);
	/*
					val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_stab_mass(
							fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], T_ip[p], tau);
	*/
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id];
				int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id + 1];

                //printf("subdomain_id = %d, IS = %d, IE = %d\n", subdomain_id, IS, IE);

				for(int k = IS; k < IE; k++){
					double val = integ_val * hlpod_mat->pod_modes[index][k];
					hlpod_ddhr->reduced_RH[k] += hlpod_ddhr->ovl_elem_weight[m] * val;
				}
			}

/*
            for(int k = 0; k < hlpod_vals->num_modes; k++){
				double val = integ_val * hlpod_mat->neib_vec[index][k];
                hlpod_ddhr->reduced_RH[k] += hlpod_ddhr->ovl_elem_weight[m] * val;
            }
*/            
		}
	}

    for(int m=0; m<(hlpod_ddhr->ovl_num_selected_elems_D_bc); m++) {
        int e = hlpod_ddhr->ovl_id_selected_elems_D_bc[m];


		BBFE_elemmat_set_Jacobian_array(
				Jacobian_ip,
				basis->num_integ_points,
				e,
				fe);

		double vol = BBFE_std_integ_calc_volume(
				basis->num_integ_points,
				basis->integ_weight,
				Jacobian_ip);

		double h_e = cbrt(vol);

		BBFE_elemmat_set_local_array_vector(local_x, fe, fe->x,   e, 3);
		BBFE_elemmat_set_local_array_scalar(local_T, fe, hr_vals->sol_vec, e);

		for(int p=0; p<np; p++) {
			BBFE_std_mapping_vector3d(x_ip[p], nl, local_x, basis->N[p]);
			manusol_get_conv_vel(v_ip[p], x_ip[p]);
			a_ip[p] = manusol_get_mass_coef(x_ip[p]);
			k_ip[p] = manusol_get_diff_coef(x_ip[p]);
			T_ip[p] = BBFE_std_mapping_scalar(nl, local_T, basis->N[p]);
			f_ip[p] = manusol_get_source(x_ip[p], t, a_ip[p], v_ip[p], k_ip[p]);
		}

		for(int i=0; i<nl; i++) {
			int index = fe->conn[e][i];
			int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index];
			if(subdomain_id < num_subdomains){

				for(int p=0; p<np; p++) {
					val_ip[p] = 0.0;
	/*
					double tau = BBFE_elemmat_convdiff_stab_coef_ns(
							k_ip[p], v_ip[p], a_ip[p], h_e, vals->dt);
	*/
					val_ip[p] += BBFE_elemmat_convdiff_vec_source(
							basis->N[p][i], f_ip[p]);
	/*
					val_ip[p] += BBFE_elemmat_convdiff_vec_stab_source(
							fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], tau, f_ip[p]);
	*/
				//	val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_mass(
					val_ip[p] += 1.0/ dt * BBFE_elemmat_convdiff_vec_mass(
							basis->N[p][i], T_ip[p], a_ip[p]);
	/*
					val_ip[p] += 1.0/(vals->dt) * BBFE_elemmat_convdiff_vec_stab_mass(
							fe->geo[e][p].grad_N[i], a_ip[p], v_ip[p], T_ip[p], tau);
	*/
				}

				double integ_val = BBFE_std_integ_calc(
						np, val_ip, basis->integ_weight, Jacobian_ip);

				int index = fe->conn[e][i];
				int subdomain_id = hlpod_mat->subdomain_id_in_nodes[index];

				//if(subdomain_id < num_subdomains){
					int IS = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id];
					int IE = hlpod_ddhr->num_neib_modes_1stdd_sum[subdomain_id + 1];

					for(int k = IS; k < IE; k++){
						double val = integ_val * hlpod_mat->pod_modes[index][k];
						hlpod_ddhr->reduced_RH[k] += hlpod_ddhr->ovl_elem_weight_D_bc[m] * val;
					}
				//}
	/*
				for(int k = 0; k < hlpod_vals->num_modes; k++){
					double val = integ_val * hlpod_mat->neib_vec[index][k];
					hlpod_ddhr->reduced_RH[k] += hlpod_ddhr->ovl_elem_weight_D_bc[m] * val;
				}
	*/
			}
		}
	}

	BB_std_free_1d_double(val_ip,      basis->num_integ_points);
	BB_std_free_1d_double(Jacobian_ip, basis->num_integ_points);

	BB_std_free_2d_double(local_x, fe->local_num_nodes, 3);
	BB_std_free_1d_double(local_T, fe->local_num_nodes);

	BB_std_free_2d_double(x_ip, np, 3);
	BB_std_free_2d_double(v_ip, np, 3);
	BB_std_free_1d_double(k_ip, np);
	BB_std_free_1d_double(a_ip, np);
	BB_std_free_1d_double(T_ip, np);
	BB_std_free_1d_double(f_ip, np);
}