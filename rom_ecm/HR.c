
/*
#include "hlpod_write.h"
#include "hlpod_write_fe.h"
#include "HR.h"
*/

#include "HR.h"

static const int BUFFER_SIZE = 10000;
static const char* INPUT_FILENAME_ELEM_ID          = "elem.dat.id";
static const char* OUTPUT_FILENAME_ECM_ELEM_VTK = "ECM_elem.vtk";


void hr_manusol_set_bc(
		BBFE_DATA* 	fe,
		BBFE_BC*   	bc,
		HLPOD_HR*   hlpod_hr,
		double     	t,
		double      (*func)(double, double, double, double) // scalar function(x, y, z, t)
		)
{
    for(int i=0; i < fe->total_num_nodes; i++) {
		bc->imposed_D_val[i] = 0.0;
    }

    for(int i = 0; i < hlpod_hr->num_selected_nodes_D_bc; i++) {
        int index = hlpod_hr->D_bc_node_id[i];
		bc->imposed_D_val[index] = func(fe->x[index][0], fe->x[index][1], fe->x[index][2], t);
    }
}

void qr_decomposition(double **A, int m, int n, double **Q, double **R) {
    int i, j, k;
    double *u = (double *)malloc(m * sizeof(double));
    double **V = (double **)malloc(m * sizeof(double *));
    for (i = 0; i < m; ++i) {
        V[i] = (double *)malloc(m * sizeof(double));
    }

    // Initialize Q as identity matrix
    for (i = 0; i < m; ++i) {
        for (j = 0; j < m; ++j) {
            Q[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    // Perform QR decomposition
    for (k = 0; k < n; ++k) {
        // Construct u vector
        for (i = 0; i < m; ++i) {
            u[i] = A[i][k];
        }
        for (j = 0; j < k; ++j) {
            double dot_product = 0.0;
            for (i = 0; i < m; ++i) {
                dot_product += Q[i][j] * A[i][k];
            }
            for (i = 0; i < m; ++i) {
                u[i] -= dot_product * Q[i][j];
            }
        }

        // Normalize u vector
        double norm_u = 0.0;
        for (i = 0; i < m; ++i) {
            norm_u += u[i] * u[i];
        }
        norm_u = sqrt(norm_u);
        for (i = 0; i < m; ++i) {
            u[i] /= norm_u;
        }

        // Update R matrix
        for (i = 0; i < m; ++i) {
            R[i][k] = u[i];
        }

        // Update Q matrix
        for (j = k; j < m; ++j) {
            double dot_product = 0.0;
            for (i = 0; i < m; ++i) {
                dot_product += Q[i][j] * u[i];
            }
            for (i = 0; i < m; ++i) {
                Q[i][j] -= 2.0 * dot_product * u[i];
            }
        }
    }

    // Free dynamically allocated memory
    free(u);
    for (i = 0; i < m; ++i) {
        free(V[i]);
    }
    free(V);
}

void sort_diagonal_indices(double **R, int m, int n, int *indices) {
    // Step 1: Get absolute values of diagonal elements of R and store with their indices
    typedef struct {
        double value;
        int index;
    } DiagonalPair;

    DiagonalPair *pairs = (DiagonalPair *)malloc(m * sizeof(DiagonalPair));

    for (int i = 0; i < m; ++i) {
        pairs[i].value = fabs(R[i][i]);
        pairs[i].index = i;
    }

    // Step 2: Sort pairs based on the absolute values (descending order)
    for (int i = 0; i < m - 1; ++i) {
        for (int j = i + 1; j < m; ++j) {
            if (pairs[i].value < pairs[j].value) {
                // Swap pairs
                DiagonalPair temp = pairs[i];
                pairs[i] = pairs[j];
                pairs[j] = temp;
            }
        }
    }

    // Step 3: Extract indices in sorted order
    for (int i = 0; i < m; ++i) {
        indices[i] = pairs[i].index;
    }

    // Free allocated memory
    free(pairs);
}


void hr_memory_allocation(
        const int       total_num_nodes,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
        HLPOD_HR*       hlpod_hr)
{
    //hr_vals->sol_vec = BB_std_calloc_1d_double(hr_vals->sol_vec, total_num_nodes);

    hlpod_hr->matrix = BB_std_calloc_2d_double(hlpod_hr->matrix, total_num_snapshot * total_num_modes, total_num_elem);
    hlpod_hr->RH = BB_std_calloc_1d_double(hlpod_hr->RH, total_num_snapshot * total_num_modes);

    hlpod_hr->reduced_mat = BB_std_calloc_2d_double(hlpod_hr->reduced_mat, total_num_modes, total_num_modes);
    hlpod_hr->reduced_RH = BB_std_calloc_1d_double(hlpod_hr->reduced_RH, total_num_modes);

}

void hr_memory_allocation_online(
        const int       total_num_nodes,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
        HLPOD_HR*       hlpod_hr)
{
    //hr_vals->sol_vec = BB_std_calloc_1d_double(hr_vals->sol_vec, total_num_nodes);

    hlpod_hr->reduced_mat = BB_std_calloc_2d_double(hlpod_hr->reduced_mat, total_num_modes, total_num_modes);
    hlpod_hr->reduced_RH = BB_std_calloc_1d_double(hlpod_hr->reduced_RH, total_num_modes);
}


void hr_get_selected_elements(
        BBFE_DATA*     	fe,
        BBFE_BC*     	bc,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
        const int       max_iter, //NNLS
        const double    tol,      //NNLS
        HLPOD_HR*       hlpod_hr,
		const char*		directory)
{
    int nl = fe->local_num_nodes;

    const int max_ITER = 2000;
    const double TOL = 1.0e-9;

    double* ans_vec;
    ans_vec = BB_std_calloc_1d_double(ans_vec, total_num_elem);

    double residual;

    int NNLS_row = total_num_snapshot * total_num_modes;

    for(int i = 0; i < NNLS_row; i++){
        //printf("matrix[%d][0] = %lf\n", i, hlpod_hr->matrix[i][0]);
        //printf("RH[%d][0] = %lf\n", i, hlpod_hr->RH[i]);
    }

    monolis_optimize_nnls_R_with_sparse_solution(
        hlpod_hr->matrix, 
        hlpod_hr->RH, 
        ans_vec, NNLS_row, total_num_elem, max_ITER, TOL, &residual);

/*
    double** Q;
    Q = BB_std_calloc_2d_double(Q, total_num_snapshot * total_num_modes * 2, total_num_elem);
    double** R;
    R = BB_std_calloc_2d_double(R, total_num_snapshot * total_num_modes * 2, total_num_elem);

    int NNLS_row = total_num_snapshot * total_num_modes;

    qr_decomposition(
        hlpod_hr->matrix,
        total_num_snapshot * total_num_modes * 2,
        total_num_elem,
        Q,
        R);

    printf("QR decomposition");

    int* diagonal_indces;
    diagonal_indces = BB_std_calloc_1d_int(diagonal_indces, total_num_elem);
    
    sort_diagonal_indices(
        R,
        total_num_snapshot * total_num_modes,
        total_num_elem,
        diagonal_indces);
    
    for(int i = 0; i < total_num_snapshot * total_num_modes * 2; i++){
        //printf("%d, %d\n", i, diagonal_indces[i]);
    }

    double** NNLS_matrix;
    NNLS_matrix = BB_std_calloc_2d_double(NNLS_matrix, NNLS_row, total_num_elem);
    double* NNLS_RH;
    NNLS_RH = BB_std_calloc_1d_double(NNLS_RH, NNLS_row);

    for(int i = 0; i < NNLS_row; i++){
        for(int j = 0; j < total_num_elem; j++){
            NNLS_matrix[i][j] = hlpod_hr->matrix[diagonal_indces[i]][j];
        }
        NNLS_RH[i] = hlpod_hr->RH[diagonal_indces[i]];
    }

        monolis_optimize_nnls_R_with_sparse_solution(
        NNLS_matrix, 
        NNLS_RH, 
        ans_vec, NNLS_row, total_num_elem, max_ITER, TOL, &residual);
*/    

    printf("\n\nmax_iter = %d, tol = %lf, residuals = %lf\n\n", max_ITER, TOL, residual);

    int index = 0;

    for(int i = 0; i < total_num_elem; i++){
        if(ans_vec[i] != 0.0){
            index++;
        }
    }

    int total_num_selected_elems = index;

    printf("\n\n num_selected_elems = %d\n\n", index);

	hr_write_NNLS_residual(residual, 0, 0, directory);
	hr_write_NNLS_num_elems(total_num_selected_elems, 0, 0, directory);

    int* total_id_selected_elems;
    double* total_elem_weight;

    total_id_selected_elems = BB_std_calloc_1d_int(total_id_selected_elems, total_num_selected_elems);
    total_elem_weight = BB_std_calloc_1d_double(total_elem_weight, total_num_selected_elems);

    bool* bool_elem;  
    bool_elem = BB_std_calloc_1d_bool(bool_elem, total_num_selected_elems);

	hlpod_hr->D_bc_exists = BB_std_calloc_1d_bool(hlpod_hr->D_bc_exists, fe->total_num_nodes);

    index = 0;
    for(int i = 0; i < total_num_elem; i++){
        if(ans_vec[i] != 0.0){
            total_id_selected_elems[index] = i;
            total_elem_weight[index] = ans_vec[i];
            index++;
        }
    }

    for(int h=0; h<(total_num_selected_elems); h++) {
        int e = total_id_selected_elems[h];

		for(int i=0; i<nl; i++) {       //六面体一次要素は8
			for(int j=0; j<nl; j++) {
                int index_i = fe->conn[e][i];
                int index_j = fe->conn[e][j];

                if( bc->D_bc_exists[index_j]) {
                    bool_elem[h] = true;
					hlpod_hr->D_bc_exists[index_j] = true;
                }
            }
        }
    }

	index = 0;
	for(int i = 0; i < fe->total_num_nodes; i++) {       //六面体一次要素は8
        if( hlpod_hr->D_bc_exists[i]){
			index++;
			printf("i = %d index = %d\n", i, index);
		}
    }

	//D_bcを付与する節点
	hlpod_hr->num_selected_nodes_D_bc = index;
	hlpod_hr->D_bc_node_id = BB_std_calloc_1d_int(hlpod_hr->D_bc_node_id, hlpod_hr->num_selected_nodes_D_bc);

	index = 0;
    for(int i = 0; i < fe->total_num_nodes; i++) {
        if( hlpod_hr->D_bc_exists[i]) {
			hlpod_hr->D_bc_node_id[index] = i;
			index++;
        }
    }

    index = 0;
    for(int h=0; h<(total_num_selected_elems); h++) {
        if(bool_elem[h]){
            index++;
        } 
    }

	printf("\n\n num_elem_D_bc = %d \n\n", index);

    //index = D_bcが付与された要素数
    hlpod_hr->num_selected_elems = total_num_selected_elems - index;
    hlpod_hr->num_selected_elems_D_bc = index;

    hlpod_hr->id_selected_elems = BB_std_calloc_1d_int(hlpod_hr->id_selected_elems, hlpod_hr->num_selected_elems);
    hlpod_hr->id_selected_elems_D_bc = BB_std_calloc_1d_int(hlpod_hr->id_selected_elems_D_bc, hlpod_hr->num_selected_elems_D_bc);
    
    hlpod_hr->elem_weight = BB_std_calloc_1d_double(hlpod_hr->elem_weight, hlpod_hr->num_selected_elems);
    hlpod_hr->elem_weight_D_bc = BB_std_calloc_1d_double(hlpod_hr->elem_weight_D_bc, hlpod_hr->num_selected_elems_D_bc);

    int index1 = 0;
    int index2 = 0;
    for(int h=0; h<(total_num_selected_elems); h++) {
        int e = total_id_selected_elems[h];

        if(bool_elem[h]) {
            hlpod_hr->id_selected_elems_D_bc[index1] = total_id_selected_elems[h];
            hlpod_hr->elem_weight_D_bc[index1] = total_elem_weight[h];
            index1++;
        }
        else{
            hlpod_hr->id_selected_elems[index2] = total_id_selected_elems[h];
            hlpod_hr->elem_weight[index2] = total_elem_weight[h];
            index2++;
        }
    }

    FILE* fp1;
    FILE* fp2;
    char fname1[BUFFER_SIZE];
    char fname2[BUFFER_SIZE];

    snprintf(fname1, BUFFER_SIZE,"DDECM/lb_selected_elem_D_bc.%d.txt", 0);
    snprintf(fname2, BUFFER_SIZE,"DDECM/lb_selected_elem.%d.txt", 0);

    fp1 = ROM_BB_write_fopen(fp1, fname1, directory);
    fp2 = ROM_BB_write_fopen(fp2, fname2, directory);

    fprintf(fp1, "%d\n", index1);
    fprintf(fp2, "%d\n", index2);

    int index_1 = 0;
    int index_2 = 0;

    for(int h=0; h<(total_num_selected_elems); h++) {
        if(bool_elem[h]) {
            fprintf(fp1, "%d %.30e\n", hlpod_hr->id_selected_elems_D_bc[index_1], hlpod_hr->elem_weight_D_bc[index_1]);
            index_1++;
        }
        else{
            fprintf(fp2, "%d %.30e\n", hlpod_hr->id_selected_elems[index_2], hlpod_hr->elem_weight[index_2]);

            index_2++;
        }
    }
    
    fclose(fp1);
    fclose(fp2);

    /*input
    hlpod_hr->matrix, 
    hlpod_hr->RH,
    NNLS conditions
    */

    /*output
    hlpod_hr->num_selected_elem
    hlpod_hr->id_selected_elems
    hlpod_hr->elem_weight
    */

}


void hr_calc_solution(
	BBFE_DATA* 		fe,
	HLPOD_MAT*      hlpod_mat,
	double*         HR_T,
	BBFE_BC*     	bc,
    int 			num_base)
{
	int nl = fe->total_num_nodes;
	int k = num_base;

	for(int j = 0; j < nl; j++){
		HR_T[j] = 0.0;
	}

	for(int i = 0; i < k; i++){
		for(int j = 0; j < nl; j++){
			HR_T[j] += hlpod_mat->pod_modes[j][i] * hlpod_mat->mode_coef[i];
		}
	}
}


void hr_monolis_set_matrix(
	MONOLIS*     	monolis,
	HLPOD_HR*      hlpod_hr,
    const int 		num_basis)
{
	const int k = num_basis;
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
		hlpod_hr->reduced_mat);		//double** val

	BB_std_free_1d_int(index, 1);
	BB_std_free_1d_int(item, 1);
	BB_std_free_1d_int(connectivity, 1);
}

void hr_to_monollis_rhs(
	MONOLIS*		monolis,
    HLPOD_HR*       hlpod_rh,
	const int 		k)
{
	for(int i = 0; i < k; i++){
        //monolis->mat.R.B[i] = 0.0;
		monolis->mat.R.B[i] = hlpod_rh->reduced_RH[i];
	}
}

void hr_set_selected_elems(
		BBFE_DATA*     	fe,
        HLPOD_HR*       hlpod_hr,
        const int		total_num_nodes,
		const char*     directory)
{
	int nl = fe->local_num_nodes;
	const int myrank = monolis_mpi_get_global_my_rank();
	double* ECM_elem;		//非ゼロ要素可視化用ベクトル
	double* ECM_elem_weight;//非ゼロ要素可視化用重みベクトル
	double* ECM_wireframe;	//wireframe可視化用ゼロベクトル

	ECM_elem = BB_std_calloc_1d_double(ECM_elem, total_num_nodes);
	ECM_elem_weight = BB_std_calloc_1d_double(ECM_elem_weight, total_num_nodes);
	ECM_wireframe = BB_std_calloc_1d_double(ECM_wireframe, total_num_nodes);

    for(int h=0; h<(hlpod_hr->num_selected_elems); h++) {
        int e = hlpod_hr->id_selected_elems[h];
        for(int i=0; i<nl; i++) {
            int index = fe->conn[e][i];
            
            ECM_elem[index] = myrank + 1;
            ECM_elem_weight[index] += hlpod_hr->elem_weight[h];
        }
    }

	for(int h=0; h<(hlpod_hr->num_selected_elems_D_bc); h++) {
        int e = hlpod_hr->id_selected_elems_D_bc[h];
        for(int i=0; i<nl; i++) {
            int index = fe->conn[e][i];
            ECM_elem[index] = myrank + 1;
            ECM_elem_weight[index] += hlpod_hr->elem_weight_D_bc[h];
        }
    }
	

	const char* filename;

	char fname[BUFFER_SIZE];

	snprintf(fname, BUFFER_SIZE,"DDECM/%s", OUTPUT_FILENAME_ECM_ELEM_VTK);
	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname);

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

	BB_vtk_write_point_vals_scalar(fp, ECM_elem, fe->total_num_nodes, "ECM-elem");
	BB_vtk_write_point_vals_scalar(fp, ECM_elem_weight, fe->total_num_nodes, "ECM-elem_weight");
	BB_vtk_write_point_vals_scalar(fp, ECM_wireframe, fe->total_num_nodes, "ECM-wireframe");

	BB_std_free_1d_double(ECM_elem, fe->total_num_nodes);

	fclose(fp);

}



void hr_lb_read_selected_elements(
    HLPOD_HR*     hlpod_hr,
	const int num_subdomains,
	const char* directory)
{
	const int myrank = monolis_mpi_get_global_my_rank();

	FILE* fp;
	char fname[BUFFER_SIZE];
	char id[BUFFER_SIZE];
	int ndof;


	FILE* fp1;
	FILE* fp2;
	FILE* fp3;
	FILE* fp4;
	char fname1[BUFFER_SIZE];
	char fname2[BUFFER_SIZE];

	char fname3[BUFFER_SIZE];
	char fname4[BUFFER_SIZE];

//	snprintf(fname3, BUFFER_SIZE, "DDECM/selected_elem_D_bc.%d.txt", monolis_mpi_get_global_my_rank());
//	snprintf(fname4, BUFFER_SIZE, "DDECM/selected_elem.%d.txt", monolis_mpi_get_global_my_rank());

//	fp3 = BBFE_sys_write_fopen(fp3, fname3, directory);
//	fp4 = BBFE_sys_write_fopen(fp4, fname4, directory);

	int Index1 = 0;
	int Index2 = 0;
	int tmp;
	double val;
	int index1 = 0;
	int index2 = 0;
	int num_selected_elems;
	int num_selected_elems_D_bc;

	for (int m = 0; m < num_subdomains; m++) {
		snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", m);
		snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", m);

		fp1 = BBFE_sys_read_fopen(fp1, fname1, directory);
		fp2 = BBFE_sys_read_fopen(fp2, fname2, directory);

		fscanf(fp1, "%d", &(num_selected_elems));
		fscanf(fp2, "%d", &(num_selected_elems_D_bc));
		Index1 += num_selected_elems;
		Index2 += num_selected_elems_D_bc;

		fclose(fp1);
		fclose(fp2);
	}

//    hlpod_hr->elem_weight = BB_std_calloc_1d_double(f_ip, np);
    hlpod_hr->num_selected_elems = Index1;
    hlpod_hr->num_selected_elems_D_bc  = Index2;

    printf("Index1 = %d, Index2 = %d\n", Index1, Index2);

    hlpod_hr->elem_weight = BB_std_calloc_1d_double(hlpod_hr->elem_weight, Index1);
    hlpod_hr->id_selected_elems = BB_std_calloc_1d_int(hlpod_hr->id_selected_elems, Index1);
    hlpod_hr->elem_weight_D_bc = BB_std_calloc_1d_double(hlpod_hr->elem_weight_D_bc, Index2);
    hlpod_hr->id_selected_elems_D_bc = BB_std_calloc_1d_int(hlpod_hr->id_selected_elems_D_bc, Index2);

	for (int m = 0; m < num_subdomains; m++) {
		snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", m);

		fp1 = BBFE_sys_read_fopen(fp1, fname1, directory);

		fscanf(fp1, "%d", &(num_selected_elems));
        for(int i = 0; i < num_selected_elems; i++){
            fscanf(fp1, "%d %lf", &(hlpod_hr->id_selected_elems[index1]), &(hlpod_hr->elem_weight[index1]));
            printf("%d %lf\n", hlpod_hr->id_selected_elems[index1], hlpod_hr->elem_weight[index1]);
            index1++;
        }

		fclose(fp1);
	}

    for (int m = 0; m < num_subdomains; m++) {
		snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", m);

		fp2 = BBFE_sys_read_fopen(fp2, fname2, directory);

		fscanf(fp2, "%d", &(num_selected_elems_D_bc));

        for(int i = 0; i < num_selected_elems_D_bc; i++){
            fscanf(fp2, "%d %lf", &(hlpod_hr->id_selected_elems_D_bc[index2]), &(hlpod_hr->elem_weight_D_bc[index2]));
            
            printf("%d %lf\n", hlpod_hr->id_selected_elems_D_bc[index2], hlpod_hr->elem_weight_D_bc[index2]);
            index2++;
        }

		fclose(fp2);
	}

}
