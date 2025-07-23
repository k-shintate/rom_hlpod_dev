
/*
#include "std.h"
#include "DDHR.h"
#include "hlpod_read.h"
#include "hlpod_write.h"
*/


#include "DDHR.h"

static const int BUFFER_SIZE = 10000;
static const char* INPUT_FILENAME_ELEM_ID          = "elem.dat.id";

void manusol_get_conv_vel3(
		double v[3],
		double x[3])
{
	v[0] = 1.0 + x[0]*x[0];
	v[1] = 1.0 + x[1]*x[1];
	v[2] = 1.0 + x[2]*x[2];
}


double manusol_get_mass_coef3(
		double x[3])
{
	double val = 1.0;

	return val;
}


double manusol_get_diff_coef3(
		double x[3])
{
	//double val = (2.0 + sin(1.0*x[0]) * sin(0.5*x[1]) * sin(0.25*x[2]));
	double val = (1.75 + sin(1.75*x[0]) * sin(0.5*x[1]) * sin(0.25*x[2]));
	//double val = 1.0;

	return val;
}


double manusol_get_source3(
		double x[3],
		double t,
		double a,
		double v[3],
		double k)
{
	//double val = -sin( 0.25*x[0] ) * sin( 0.5*x[1] ) * sin( 1.0*x[2] ) *(-(1/16.0 + 1/4.0 +1.0)*sin( 1.0*t )-cos(1.0*t));
	double val = 0.0;

	return val;
}

void ddhr_memory_allocation(
        const int       total_num_nodes,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
		const int		num_subdomains,
        HLPOD_DDHR*     hlpod_ddhr)
{

    //hr_vals->sol_vec = BB_std_calloc_1d_double(hr_vals->sol_vec, total_num_nodes);

//for NNLS
    hlpod_ddhr->matrix = BB_std_calloc_3d_double(hlpod_ddhr->matrix, total_num_snapshot * total_num_modes * 2, total_num_elem, num_subdomains);
    hlpod_ddhr->RH = BB_std_calloc_2d_double(hlpod_ddhr->RH, total_num_snapshot * total_num_modes * 2, num_subdomains);

    hlpod_ddhr->reduced_mat = BB_std_calloc_2d_double(hlpod_ddhr->reduced_mat, total_num_modes * num_subdomains, total_num_modes * num_subdomains);
    hlpod_ddhr->reduced_RH = BB_std_calloc_1d_double(hlpod_ddhr->reduced_RH, total_num_modes * num_subdomains);

}
/********/

void ddhr_set_element(
        HLPOD_DDHR*     hlpod_ddhr,
		const int 		num_subdomains,
		const char*     directory)
{
	int id;
	int tmp;
    char fname[BUFFER_SIZE];
    char char_id[BUFFER_SIZE];
	FILE* fp;

	int* num_elems;
	num_elems = BB_std_calloc_1d_int(num_elems, num_subdomains);

    for(int m = 0; m < num_subdomains; m++){
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_ELEM_ID, m);

        fp = ROM_BB_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", char_id);
        fscanf(fp, "%d %d", &(num_elems[m]), &(tmp));
        fclose(fp);
	}

	int max_num_elem = ROM_BB_findMax(num_elems, num_subdomains);

	hlpod_ddhr->elem_id_local = BB_std_calloc_2d_int(hlpod_ddhr->elem_id_local, max_num_elem, num_subdomains);
	hlpod_ddhr->num_elems = BB_std_calloc_1d_int(hlpod_ddhr->num_elems, num_subdomains);

	//要素数の取得
    for(int m = 0; m < num_subdomains; m++){
		hlpod_ddhr->num_elems[m] = num_elems[m];
	}

    for(int m = 0; m < num_subdomains; m++){	
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_ELEM_ID, m);

        fp = ROM_BB_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", char_id);
        fscanf(fp, "%d %d", &(num_elems[m]), &(tmp));

        for(int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
            fscanf(fp, "%d", &(hlpod_ddhr->elem_id_local[i][m]));
        }

        fclose(fp);
	}

}

void ddhr_get_selected_elements(
        BBFE_DATA*     	fe,
        BBFE_BC*     	bc,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
		const int 		num_subdomains,
        const int       max_iter, //NNLS
        const double    tol,      //NNLS
        HLPOD_DDHR*       hlpod_ddhr,
		const char*		directory)
{
    int nl = fe->local_num_nodes;

    const int max_ITER = 400;
    const double TOL = 1.0e-6;

    double residual;

    int NNLS_row = total_num_snapshot * total_num_modes *2;

    double* ans_vec;
	double** matrix;
	double* RH;
    bool* bool_elem;  
	int* total_id_selected_elems;
    double* total_elem_weight;

/**/
	hlpod_ddhr->D_bc_exists = BB_std_calloc_2d_bool(hlpod_ddhr->D_bc_exists, fe->total_num_nodes, num_subdomains);

	hlpod_ddhr->id_selected_elems = BB_std_calloc_2d_int(hlpod_ddhr->id_selected_elems, max_ITER, num_subdomains);
    hlpod_ddhr->id_selected_elems_D_bc = BB_std_calloc_2d_int(hlpod_ddhr->id_selected_elems_D_bc, max_ITER, num_subdomains);
    
    hlpod_ddhr->elem_weight = BB_std_calloc_2d_double(hlpod_ddhr->elem_weight, max_ITER, num_subdomains);
    hlpod_ddhr->elem_weight_D_bc = BB_std_calloc_2d_double(hlpod_ddhr->elem_weight_D_bc, max_ITER, num_subdomains);

	hlpod_ddhr->num_selected_nodes_D_bc = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_nodes_D_bc, num_subdomains);
	hlpod_ddhr->D_bc_node_id = BB_std_calloc_2d_int(hlpod_ddhr->D_bc_node_id, max_ITER, num_subdomains);

	hlpod_ddhr->num_selected_elems = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_elems, num_subdomains);
	hlpod_ddhr->num_selected_elems_D_bc = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_elems_D_bc, num_subdomains);
/**/

for(int m = 0; m < num_subdomains; m++){
	printf("m = %d, num_elems[m] = %d ", m, hlpod_ddhr->num_elems[m]);

    ans_vec = BB_std_calloc_1d_double(ans_vec, hlpod_ddhr->num_elems[m]);
	matrix = BB_std_calloc_2d_double(matrix, NNLS_row, hlpod_ddhr->num_elems[m]);
	RH = BB_std_calloc_1d_double(RH, NNLS_row);

	for(int j = 0; j < NNLS_row; j++){
    	for (int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
			matrix[j][i] = hlpod_ddhr->matrix[j][i][m];
		}
		RH[j] = hlpod_ddhr->RH[j][m];
    }

	residual = 0.0;

    monolis_optimize_nnls_R_with_sparse_solution(
		matrix,
		RH,
        ans_vec, NNLS_row, hlpod_ddhr->num_elems[m], max_ITER, TOL, &residual);

    printf("\n\nmax_iter = %d, tol = %lf, residuals = %lf\n\n", max_ITER, TOL, residual);

    int index = 0;

    for(int i = 0; i < hlpod_ddhr->num_elems[m]; i++){
        if(ans_vec[i] != 0.0){
            index++;
        }
    }

    int total_num_selected_elems = index;

    printf("\n\nnum_selected_elems = %d\n\n", index);

	hr_write_NNLS_residual(residual, m, 0, directory);
	hr_write_NNLS_num_elems(total_num_selected_elems, m, 0, directory);

    bool_elem = BB_std_calloc_1d_bool(bool_elem, total_num_selected_elems);
    total_id_selected_elems = BB_std_calloc_1d_int(total_id_selected_elems, total_num_selected_elems);
    total_elem_weight = BB_std_calloc_1d_double(total_elem_weight, total_num_selected_elems);

    index = 0;
    for(int i = 0; i < hlpod_ddhr->num_elems[m]; i++){
        if(ans_vec[i] != 0.0){
            total_id_selected_elems[index] = hlpod_ddhr->elem_id_local[i][m];
            total_elem_weight[index] = ans_vec[i];
            index++;
        }
    }

    for(int h = 0; h < (total_num_selected_elems); h++) {
        int e = total_id_selected_elems[h];

		for(int i=0; i<nl; i++) {       //六面体一次要素は8
			for(int j=0; j<nl; j++) {
                //基底本数ループ
                int index_i = fe->conn[e][i];
                int index_j = fe->conn[e][j];

                if( bc->D_bc_exists[index_j]) {
                    bool_elem[h] = true;
					hlpod_ddhr->D_bc_exists[index_j][m] = true;
                }
            }
        }
    }

	index = 0;
	for(int i = 0; i < fe->total_num_nodes; i++) {       //六面体一次要素は8
        if( hlpod_ddhr->D_bc_exists[i][m]){
			index++;
			printf("i = %d index = %d\n", i, index);
		}
    }

	index = 0;
    for(int i = 0; i < fe->total_num_nodes; i++) {
        if( hlpod_ddhr->D_bc_exists[i][m]) {
			hlpod_ddhr->D_bc_node_id[index][m] = i;
			index++;
        }
    }

    index = 0;
    for(int h = 0; h < (total_num_selected_elems); h++) {
        if(bool_elem[h]){
            index++;
        } 
    }

	printf("\n\n num_elem_D_bc = %d \n\n", index);

    hlpod_ddhr->num_selected_elems[m] = total_num_selected_elems - index;
    hlpod_ddhr->num_selected_elems_D_bc[m] = index;

    int index1 = 0;
    int index2 = 0;
    for(int h=0; h<(total_num_selected_elems); h++) {
        int e = total_id_selected_elems[h];

        if(bool_elem[h]) {
            hlpod_ddhr->id_selected_elems_D_bc[index1][m] = total_id_selected_elems[h];
            hlpod_ddhr->elem_weight_D_bc[index1][m] = total_elem_weight[h];
            printf("D_bc = %d, elem_weight = %lf \n", hlpod_ddhr->id_selected_elems_D_bc[index1][m], hlpod_ddhr->elem_weight_D_bc[index1][m]);
            index1++;
        }
        else{
            hlpod_ddhr->id_selected_elems[index2][m] = total_id_selected_elems[h];
            hlpod_ddhr->elem_weight[index2][m] = total_elem_weight[h];
			printf("%d, elem_weight = %lf\n", hlpod_ddhr->id_selected_elems[index2][m], hlpod_ddhr->elem_weight[index2][m]);
            index2++;
        }
    }

	BB_std_free_1d_double(ans_vec, hlpod_ddhr->num_elems[m]);
	BB_std_free_2d_double(matrix, NNLS_row, hlpod_ddhr->num_elems[m]);
	BB_std_free_1d_double(RH, hlpod_ddhr->num_elems[m]);

    BB_std_free_1d_bool(bool_elem, total_num_selected_elems);
	BB_std_free_1d_int(total_id_selected_elems, total_num_selected_elems);
    BB_std_free_1d_double(total_elem_weight, total_num_selected_elems);

}

    /*input
    hlpod_ddhr->matrix, 
    hlpod_ddhr->RH,
    NNLS conditions
    */

    /*output
    hlpod_ddhr->num_selected_elem
    hlpod_ddhr->id_selected_elems
    hlpod_ddhr->elem_weight
    */

}


void ddhr_calc_solution(
	BBFE_DATA* 		fe,
    HR_VALUES*      hr_vals,
	HLPOD_MAT*      hlpod_mat,
	HLPOD_DDHR*     hlpod_ddhr,
	BBFE_BC*     	bc,
    int 			num_base,
	const int		num_subdomains,
	const int		dof)
{
	int nl = fe->total_num_nodes;
	int k = num_base * num_subdomains;

	double t1 = monolis_get_time();

	for(int j = 0; j < nl; j++){
		hr_vals->sol_vec[j] = 0.0;
	}
/*
	for(int i = 0; i < k; i++){
		for(int j = 0; j < nl; j++){
			hr_vals->sol_vec[j] += hlpod_mat->pod_modes[j][i] * hlpod_mat->mode_coef[i];
		}
	}
*/
	int index_row = 0;
	int index_column1 = 0;
	int index_column2 = 0;
	int sum = 0;
	for(int k = 0; k < num_subdomains; k++){
		for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
		//for(int i = 0; i < 20; i++){
			for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
				for(int l = 0; l < dof; l++){
					//index_row = hlpod_mat->node_id[j]*dof + l + sum;
					index_row = hlpod_mat->node_id[j+sum]*dof + l;
					hr_vals->sol_vec[index_row] += hlpod_mat->pod_modes[index_row][index_column1 + i] * hlpod_mat->mode_coef[index_column1 + i];
				}
			}
		}
		index_column1 += hlpod_mat->num_modes_internal[k];
		//index_column += 20;
		index_column2 += num_base;
		sum += hlpod_mat->n_internal_vertex_subd[k] * dof;
	}

	for(int i = 0; i < nl; i++) {
		if( bc->D_bc_exists[i] ) {
			hr_vals->sol_vec[i] = bc->imposed_D_val[i];
		}
	}

	double t2 = monolis_get_time();
	//lpod_prm->time_calc_sol = t2-t1;
}


void ddhr_monolis_set_matrix(
	MONOLIS*     	monolis,
	HLPOD_DDHR*      hlpod_ddhr,
    const int 		num_basis,
	const int		num_subdomains)
{
	const int k = num_basis * num_subdomains;
	int* index;
	int* item;
	int* connectivity;

	index = BB_std_calloc_1d_int(index, 1);
	item = BB_std_calloc_1d_int(item, 1);
	
	index[0] = 0;
	item[0] = 1;

	monolis_get_nonzero_pattern_by_nodal_graph_R(
		monolis,
		1,					//nnode:節点数
		k,					//ndof:節点あたりの自由度
		index,				//節点グラフを定義するindex配列
		item);				//設定グラフを定義するitem配列

	connectivity = BB_std_calloc_1d_int(connectivity, 1);
	connectivity[0] = 0;

	monolis_add_matrix_to_sparse_matrix_R(
		monolis,					//MONOLIS* mat,
  		1,							//int      n_base,
  		connectivity,				//int*     connectivity,
		hlpod_ddhr->reduced_mat);	//double** val

	BB_std_free_1d_int(index, 1);
	BB_std_free_1d_int(item, 1);
	BB_std_free_1d_int(connectivity, 1);
}

void ddhr_monolis_set_matrix2(
	MONOLIS*     	monolis,
	HLPOD_MAT*     hlpod_mat,
	HLPOD_DDHR*     hlpod_ddhr,
    HLPOD_META*     hlpod_meta,
    const int 		num_base,
	const int		num_2nddd)
{
	double** matrix;
	matrix = BB_std_calloc_2d_double(matrix, num_base*num_2nddd, num_base*num_2nddd);

	int Index_column1 = 0;
	int Index_column2 = 0;

	for(int k1 = 0; k1 < num_2nddd; k1++){
		for(int i1 = 0; i1 < hlpod_mat->num_modes_internal[k1]; i1++){
		int index_row = 0;
		int sum = 0;
		int index_column1 = 0;
		int index_column2 = 0;
			for(int k = 0; k < num_2nddd; k++){
				for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
				//for(int i = 0; i < 20; i++){
					for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
						index_row = hlpod_mat->node_id[j + sum];
						//hlpod_mat->UTAU[index_column2 + i][Index_column2 + i1] += hlpod_mat->pod_basis[index_row][index_column1 + i] * hlpod_mat->AU[index_row][Index_column1 + i1];
						matrix[index_column2 + i][Index_column2 + i1] = hlpod_ddhr->reduced_mat[index_column1 + i][Index_column1 + i1];
					}
				}
				index_column1 += hlpod_mat->num_modes_internal[k];
				index_column2 += num_base;
				sum += hlpod_mat->n_internal_vertex_subd[k];
			}
		}
		Index_column1 += hlpod_mat->num_modes_internal[k1];
		Index_column2 += num_base;
	}


	for(int k = 0; k < num_2nddd; k++){
		for(int m = 0; m < num_base; m++){
			for(int n = 0; n < num_base; n++){
        	    double val = matrix[k * num_base + m][k * num_base + n];
                if(m<hlpod_mat->num_modes_internal[k] && n<hlpod_mat->num_modes_internal[k]){
                    monolis_add_scalar_to_sparse_matrix_R(
                        monolis,
                        k,
                        k,
                        m,
                        n,
                        val);
                }
			}
		}
	}
	
	for(int k = 0; k < num_2nddd; k++){
	int iS = hlpod_meta->index[k];
	int iE = hlpod_meta->index[k + 1];

	for(int i = iS; i < iE; i++){
			for(int m = 0; m < num_base; m++){
				for(int n = 0; n < num_base; n++){
					double val = matrix[k * num_base + m][hlpod_meta->item[i] * num_base + n];

                    if(m<hlpod_mat->num_modes_internal[k]&& n<hlpod_mat->num_modes_internal[hlpod_meta->item[i]]){
                        monolis_add_scalar_to_sparse_matrix_R(
                            monolis,
                            k,
                            hlpod_meta->item[i],
                            m,
                            n,
                            val);
                    }
				}
			}
		}
	}

	BB_std_free_2d_double(matrix, num_base*num_2nddd, num_base*num_2nddd);

}

void ddhr_to_monollis_rhs(
	MONOLIS*		monolis,
	HLPOD_MAT*     hlpod_mat,
    HLPOD_DDHR*     hlpod_ddhr,
	const int 		num_base,
	const int		num_subdomains)
{
	int index_column1 = 0;
	int index_column2 = 0;

	for(int k = 0; k < num_subdomains; k++){
		for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
			monolis->mat.R.B[i + index_column2] = hlpod_ddhr->reduced_RH[i + index_column1];
		}
		index_column1 += hlpod_mat->num_modes_internal[k];
        index_column2 += hlpod_mat->num_modes_internal[k];
	}

}
