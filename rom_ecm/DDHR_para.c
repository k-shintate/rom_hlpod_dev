//DDHROMに関して、オーバーラップ要素を含んで計算する方式
//内部要素の総和が分割前の要素数の総和になることを利用

/*
#include "rom_dataset.h"
#include "diff_dataset.h"
#include "hlpod_write_fe.h"
#include "hlpod_write.h"
#include "std.h"
#include "hlpod_read.h"

#include "monolis_nnls_c.h"
*/


#include "DDHR_para.h"

static const int BUFFER_SIZE = 10000;
static const char* INPUT_FILENAME_ELEM_ID          = "elem.dat.id";
static const char* OUTPUT_FILENAME_ECM_ELEM_VTK = "ECM_elem.vtk";

void ddhr_memory_allocation_para_online(
        HLPOD_VALUES*   hlpod_vals,
	    HLPOD_DDHR*     hlpod_ddhr,
	    HLPOD_MAT*      hlpod_mat,
        const int       total_num_nodes)
{
	hlpod_ddhr->HR_T = BB_std_calloc_1d_double(hlpod_ddhr->HR_T, total_num_nodes);
    hlpod_ddhr->reduced_mat = BB_std_calloc_2d_double(hlpod_ddhr->reduced_mat, hlpod_vals->n_neib_vec, hlpod_vals->n_neib_vec);
    hlpod_ddhr->reduced_RH = BB_std_calloc_1d_double(hlpod_ddhr->reduced_RH, hlpod_vals->n_neib_vec);
}

void ddhr_memory_allocation_para(
        HLPOD_VALUES*   hlpod_vals,
	    HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_MAT*      hlpod_mat,
        const int       total_num_nodes,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
		const int		num_subdomains)
{
	int max_num_elem = 0;
	if(num_subdomains==1){
		max_num_elem = total_num_elem;
	}
	else{
		max_num_elem = ROM_BB_findMax(hlpod_ddhr->num_elems, num_subdomains);
	}

    hlpod_ddhr->HR_T = BB_std_calloc_1d_double(hlpod_ddhr->HR_T, total_num_nodes);

//for NNLS
    hlpod_ddhr->matrix = BB_std_calloc_3d_double(hlpod_ddhr->matrix, total_num_snapshot*hlpod_vals->n_neib_vec, max_num_elem, num_subdomains);
    hlpod_ddhr->RH = BB_std_calloc_2d_double(hlpod_ddhr->RH, total_num_snapshot*hlpod_vals->n_neib_vec, num_subdomains);

    hlpod_ddhr->reduced_mat = BB_std_calloc_2d_double(hlpod_ddhr->reduced_mat, hlpod_vals->n_neib_vec, hlpod_vals->n_neib_vec);
    hlpod_ddhr->reduced_RH = BB_std_calloc_1d_double(hlpod_ddhr->reduced_RH, hlpod_vals->n_neib_vec);

}

/*
void HROM_set_ansvec_para(
		VALUES*        vals,
		HLPOD_DDHR*     hlpod_ddhr,
		const int       total_num_nodes)
{
	for(int i = 0; i < total_num_nodes; i++){
		hlpod_ddhr->HR_T[i] = vals->T[i];
	}
}
*/


void ddhr_set_element_para(
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
        snprintf(fname, BUFFER_SIZE, "parted.0/elem.dat.n_internal.%d", monolis_mpi_get_global_my_rank());

        fp = ROM_BB_read_fopen(fp, fname, directory);
        fscanf(fp, "%s %d", char_id, &(tmp));
        fscanf(fp, "%d", &(num_elems[m]));
		printf("num_elems = %d\n", num_elems[m]);
        fclose(fp);
	}

	int max_num_elem = ROM_BB_findMax(num_elems, num_subdomains);

	hlpod_ddhr->elem_id_local = BB_std_calloc_2d_int(hlpod_ddhr->elem_id_local, max_num_elem, num_subdomains);
	hlpod_ddhr->num_elems = BB_std_calloc_1d_int(hlpod_ddhr->num_elems, num_subdomains);

    for(int m = 0; m < num_subdomains; m++){
		hlpod_ddhr->num_elems[m] = num_elems[m];
	}

    for(int m = 0; m < num_subdomains; m++){	
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_ELEM_ID, monolis_mpi_get_global_my_rank());

        fp = ROM_BB_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", char_id);
        fscanf(fp, "%d %d", &(num_elems[m]), &(tmp));

        for(int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
            fscanf(fp, "%d", &(hlpod_ddhr->elem_id_local[i][m]));
        }

        fclose(fp);
	}


	//DDECM_paraで追加した内容
	hlpod_ddhr->total_num_elems = BB_std_calloc_1d_int(hlpod_ddhr->total_num_elems, num_subdomains);	//ovl要素も含んだ全要素数

	for(int m = 0; m < num_subdomains; m++){	
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_ELEM_ID, monolis_mpi_get_global_my_rank());

        fp = ROM_BB_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", char_id);
        fscanf(fp, "%d %d", &(hlpod_ddhr->total_num_elems[m]), &(tmp));

		if(monolis_mpi_get_global_my_rank()==0){
		}
		else{
			hlpod_ddhr->ovl_elem_global_id = BB_std_calloc_2d_int(hlpod_ddhr->ovl_elem_global_id, hlpod_ddhr->total_num_elems[m] - hlpod_ddhr->num_elems[m], num_subdomains);

			for(int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
				fscanf(fp, "%d", &(tmp));
			}

			for(int i = 0; i < hlpod_ddhr->total_num_elems[m] - hlpod_ddhr->num_elems[m]; i++) {
				fscanf(fp, "%d", &(hlpod_ddhr->ovl_elem_global_id[i][m]));
			}
		}
        fclose(fp);
	}

}


void ddhr_get_selected_elements_para(
        BBFE_DATA*     	fe,
        BBFE_BC*     	bc,
        HLPOD_VALUES*   hlpod_vals,
	    HLPOD_DDHR*     hlpod_ddhr,
		HLPOD_MAT*    hlpod_mat,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
		const int 		num_subdomains,
        const int       max_iter, //NNLS
        const double    tol,      //NNLS
		const char*		directory)
{
    int nl = fe->local_num_nodes;

    const int max_ITER = 400;
    const double TOL = 1.0e-6;

    double residual;

    int NNLS_row = total_num_snapshot*hlpod_vals->n_neib_vec*2;	//2は残差ベクトル＋右辺ベクトルを採用しているため
	const int myrank = monolis_mpi_get_global_my_rank();

    double* ans_vec;
	double** matrix;
	double* RH;
	bool** bool_elem;  
	int* total_id_selected_elems;
    double* total_elem_weight;

	int* total_num_selected_elems;
/**/
	hlpod_ddhr->D_bc_exists = BB_std_calloc_2d_bool(hlpod_ddhr->D_bc_exists, fe->total_num_nodes, num_subdomains);

	hlpod_ddhr->id_selected_elems = BB_std_calloc_2d_int(hlpod_ddhr->id_selected_elems, max_ITER, num_subdomains);
    hlpod_ddhr->id_selected_elems_D_bc = BB_std_calloc_2d_int(hlpod_ddhr->id_selected_elems_D_bc, max_ITER, num_subdomains);
    
    hlpod_ddhr->elem_weight = BB_std_calloc_2d_double(hlpod_ddhr->elem_weight, max_ITER, num_subdomains);
    hlpod_ddhr->elem_weight_D_bc = BB_std_calloc_2d_double(hlpod_ddhr->elem_weight_D_bc, max_ITER, num_subdomains);

	hlpod_ddhr->num_selected_elems = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_elems, num_subdomains);
	hlpod_ddhr->num_selected_elems_D_bc = BB_std_calloc_1d_int(hlpod_ddhr->num_selected_elems_D_bc, num_subdomains);
/**/
	//lbから追加
	bool_elem = BB_std_calloc_2d_bool(bool_elem, max_ITER, num_subdomains);
	total_num_selected_elems = BB_std_calloc_1d_int(total_num_selected_elems, num_subdomains);

	FILE* fp1;
	FILE* fp2;
	char fname1[BUFFER_SIZE];
	char fname2[BUFFER_SIZE];

	snprintf(fname1, BUFFER_SIZE,"DDECM/selected_elem_D_bc.%d.txt", monolis_mpi_get_global_my_rank());
	snprintf(fname2, BUFFER_SIZE,"DDECM/selected_elem.%d.txt", monolis_mpi_get_global_my_rank());

	fp1 = ROM_BB_write_fopen(fp1, fname1, directory);
	fp2 = ROM_BB_write_fopen(fp2, fname2, directory);

	int Index1 = 0;
	int Index2 = 0;

	for(int m = 0; m < num_subdomains; m++){
		printf("m = %d, num_elems[m] = %d ", m, hlpod_ddhr->num_elems[m]);

		ans_vec = BB_std_calloc_1d_double(ans_vec, hlpod_ddhr->num_elems[m]);
		matrix = BB_std_calloc_2d_double(matrix, NNLS_row, hlpod_ddhr->num_elems[m]);
		RH = BB_std_calloc_1d_double(RH, NNLS_row);

		/*
		if(myrank==0){
			for(int i = 5000; i < NNLS_row; i++){
				printf("i = %d, %lf\n", i, hlpod_ddhr->RH[i][m]);
			}
		}
		*/

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

		total_num_selected_elems[m] = index;

		printf("\n\nnum_selected_elems = %d\n\n", index);

		hr_write_NNLS_residual(residual, myrank, m, directory);
		hr_write_NNLS_num_elems(total_num_selected_elems[m], myrank, m, directory);

		total_id_selected_elems = BB_std_calloc_1d_int(total_id_selected_elems, total_num_selected_elems[m]);
		total_elem_weight = BB_std_calloc_1d_double(total_elem_weight, total_num_selected_elems[m]);

		index = 0;
		for(int i = 0; i < hlpod_ddhr->num_elems[m]; i++){
			if(ans_vec[i] != 0.0){
				total_id_selected_elems[index] = i;
				total_elem_weight[index] = ans_vec[i];
				index++;
			}
		}

		for(int h = 0; h < (total_num_selected_elems[m]); h++) {
			int e = total_id_selected_elems[h];

			for(int i=0; i<nl; i++) {       //六面体一次要素は8
				for(int j=0; j<nl; j++) {
					//基底本数ループ
					int index_i = fe->conn[e][i];
					int index_j = fe->conn[e][j];

					if(bc->D_bc_exists[index_j]) {
						bool_elem[h][m] = true;
						hlpod_ddhr->D_bc_exists[index_j][m] = true;
					}
				}
			}
		}

		index = 0;
		for(int i = 0; i < fe->total_num_nodes; i++) {       //六面体一次要素は8
			if( hlpod_ddhr->D_bc_exists[i][m]){
				index++;
			}
		}

		index = 0;
		for(int h = 0; h < (total_num_selected_elems[m]); h++) {
			if(bool_elem[h][m]){
				index++;
			} 
		}

		printf("\n\n num_elem_D_bc = %d \n\n", index);

		//index = D_bcが付与された要素数
		hlpod_ddhr->num_selected_elems[m] = total_num_selected_elems[m] - index;
		hlpod_ddhr->num_selected_elems_D_bc[m] = index;

		int index1 = 0;
		int index2 = 0;
		for(int h=0; h<(total_num_selected_elems[m]); h++) {
			int e = total_id_selected_elems[h];

			if(bool_elem[h][m]) {
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

			Index1 += index1;
			Index2 += index2;

		}

		BB_std_free_1d_double(ans_vec, hlpod_ddhr->num_elems[m]);
		BB_std_free_2d_double(matrix, NNLS_row, hlpod_ddhr->num_elems[m]);
		BB_std_free_1d_double(RH, NNLS_row);

		BB_std_free_1d_int(total_id_selected_elems, total_num_selected_elems[m]);
		BB_std_free_1d_double(total_elem_weight, total_num_selected_elems[m]);
	}


	fprintf(fp1, "%d\n", Index1);
	fprintf(fp2, "%d\n", Index2);

	for(int m = 0; m < num_subdomains; m++){
		int index1 = 0;
		int index2 = 0;
		for(int h=0; h<(total_num_selected_elems[m]); h++) {
			if(bool_elem[h][m]) {
				fprintf(fp1, "%d %.15f\n", hlpod_ddhr->elem_id_local[hlpod_ddhr->id_selected_elems_D_bc[index1][m]][m], hlpod_ddhr->elem_weight_D_bc[index1][m]);
				index1++;
			}
			else{
				fprintf(fp2, "%d %.15f\n", hlpod_ddhr->elem_id_local[hlpod_ddhr->id_selected_elems[index2][m]][m], hlpod_ddhr->elem_weight[index2][m]);
				index2++;
			}
		}
	}

	fclose(fp1);
	fclose(fp2);

	/*lbから追加*/
	BB_std_free_2d_bool(bool_elem, max_ITER, num_subdomains);	
	/**/

    BB_std_free_3d_double(hlpod_ddhr->matrix, total_num_snapshot*hlpod_vals->n_neib_vec*2, total_num_elem, num_subdomains);
    BB_std_free_2d_double(hlpod_ddhr->RH, total_num_snapshot*hlpod_vals->n_neib_vec*2, num_subdomains);
}
/********/


void ddhr_get_selected_elements_para_add(
	    HLPOD_DDHR*     hlpod_ddhr,
		const int		num_parallel_subdomains,
		const char*		directory)
{
	double t = monolis_get_time_global_sync();

	const int myrank = monolis_mpi_get_global_my_rank();
	FILE* fp1;
	FILE* fp2;
	char fname1[BUFFER_SIZE];
	char fname2[BUFFER_SIZE];

	int val;

	int*	ovl_selected_elems;
	int*	ovl_selected_elems_D_bc;
	double*	ovl_selected_elems_weight;
	double*	ovl_selected_elems_weight_D_bc;

	int num_selected_elems = 0;
	int num_selected_elems_D_bc = 0;

	for(int m = 0; m < num_parallel_subdomains; m++){
		snprintf(fname1, BUFFER_SIZE,"DDECM/selected_elem.%d.txt", m);
		fp1 = ROM_BB_read_fopen(fp1, fname1, directory);

        fscanf(fp1, "%d", &(val));
		num_selected_elems += val;
		fclose(fp1);
	}

	for(int m = 0; m < num_parallel_subdomains; m++){
		snprintf(fname2, BUFFER_SIZE,"DDECM/selected_elem_D_bc.%d.txt", m);
		fp2 = ROM_BB_read_fopen(fp2, fname2, directory);

        fscanf(fp2, "%d", &(val));
		num_selected_elems_D_bc += val;
		fclose(fp2);
	}

	ovl_selected_elems = BB_std_calloc_1d_int(ovl_selected_elems, num_selected_elems);
	ovl_selected_elems_weight = BB_std_calloc_1d_double(ovl_selected_elems_weight, num_selected_elems);
	ovl_selected_elems_D_bc = BB_std_calloc_1d_int(ovl_selected_elems_D_bc, num_selected_elems_D_bc);
	ovl_selected_elems_weight_D_bc = BB_std_calloc_1d_double(ovl_selected_elems_weight, num_selected_elems_D_bc);

	int index = 0;

	for(int m = 0; m < num_parallel_subdomains; m++){
		snprintf(fname1, BUFFER_SIZE,"DDECM/selected_elem.%d.txt", m);
		fp1 = ROM_BB_read_fopen(fp1, fname1, directory);

        fscanf(fp1, "%d", &(val));

		for(int j = 0; j < val; j++){
			fscanf(fp1, "%d %lf", &(ovl_selected_elems[j+index]), &(ovl_selected_elems_weight[j+index]));
			printf("%d %d %lf\n", myrank, (ovl_selected_elems[j+index]), (ovl_selected_elems_weight[j+index]));
		}

		index += val;

		fclose(fp1);
	}

	index = 0;
	for(int m = 0; m < num_parallel_subdomains; m++){
		snprintf(fname2, BUFFER_SIZE,"DDECM/selected_elem_D_bc.%d.txt", m);
		fp2 = ROM_BB_read_fopen(fp2, fname2, directory);

        fscanf(fp2, "%d", &(val));

		for(int j = 0; j < val; j++){
			fscanf(fp2, "%d %lf", &(ovl_selected_elems_D_bc[j+index]), &(ovl_selected_elems_weight_D_bc[j+index]));
			printf("%d %d %lf\n", myrank,(ovl_selected_elems_D_bc[j+index]), (ovl_selected_elems_weight_D_bc[j+index]));
		}

		index += val;

		fclose(fp2);
	}

	bool*	bool_ovl_selected_elems;
	bool*	bool_ovl_selected_elems_D_bc;

	bool_ovl_selected_elems = BB_std_calloc_1d_bool(bool_ovl_selected_elems, num_selected_elems);
	bool_ovl_selected_elems_D_bc = BB_std_calloc_1d_bool(bool_ovl_selected_elems_D_bc, num_selected_elems_D_bc);

	int*	ovl_elem_local_id;
	int*	ovl_elem_local_id_D_bc;

	ovl_elem_local_id = BB_std_calloc_1d_int(ovl_elem_local_id, num_selected_elems);
	ovl_elem_local_id_D_bc = BB_std_calloc_1d_int(ovl_elem_local_id_D_bc, num_selected_elems_D_bc);


    int index1 = 0;
    int index2 = 0;

	if(monolis_mpi_get_global_my_rank()==0){
	}
	else{
		for(int i = 0; i < num_selected_elems; i++){
			for(int j = 0; j < hlpod_ddhr->total_num_elems[0] - hlpod_ddhr->num_elems[0]; j++){
				if(ovl_selected_elems[i] == hlpod_ddhr->ovl_elem_global_id[j][0]){
					bool_ovl_selected_elems[i] = true;
					ovl_elem_local_id[index1] = j + hlpod_ddhr->num_elems[0];
					index1++;
				}
			}
		}
		for(int i = 0; i < num_selected_elems_D_bc; i++){
			for(int j = 0; j < hlpod_ddhr->total_num_elems[0] - hlpod_ddhr->num_elems[0]; j++){
				if(ovl_selected_elems_D_bc[i] == hlpod_ddhr->ovl_elem_global_id[j][0]){
					bool_ovl_selected_elems_D_bc[i] = true;
					ovl_elem_local_id_D_bc[index2] = j + hlpod_ddhr->num_elems[0];
					index2++;
				}
			}
		}
	}

	printf("myrank = %d index1 = %d, index2 = %d\n", monolis_mpi_get_global_my_rank(), index1, index2);

	hlpod_ddhr->ovl_id_selected_elems = BB_std_calloc_1d_int(hlpod_ddhr->ovl_id_selected_elems, index1);
	hlpod_ddhr->ovl_elem_weight = BB_std_calloc_1d_double(hlpod_ddhr->ovl_elem_weight, index1);
	hlpod_ddhr->ovl_id_selected_elems_D_bc = BB_std_calloc_1d_int(hlpod_ddhr->ovl_id_selected_elems_D_bc, index2);
	hlpod_ddhr->ovl_elem_weight_D_bc = BB_std_calloc_1d_double(hlpod_ddhr->ovl_elem_weight_D_bc, index2);

	hlpod_ddhr->ovl_num_selected_elems = index1;
	hlpod_ddhr->ovl_num_selected_elems_D_bc = index2;

	index1 = 0;
    index2 = 0;

	if(monolis_mpi_get_global_my_rank()==0){
	}
	else{
		for(int i = 0; i < num_selected_elems; i++){
			if(bool_ovl_selected_elems[i]){
				hlpod_ddhr->ovl_id_selected_elems[index1] = ovl_elem_local_id[index1];
				hlpod_ddhr->ovl_elem_weight[index1] = ovl_selected_elems_weight[i];
				index1++;
			}
		}

		for(int i = 0; i < num_selected_elems_D_bc; i++){
			if(bool_ovl_selected_elems_D_bc[i]){
				hlpod_ddhr->ovl_id_selected_elems_D_bc[index2] = ovl_elem_local_id_D_bc[index2];
				hlpod_ddhr->ovl_elem_weight_D_bc[index2] = ovl_selected_elems_weight_D_bc[i];
				index2++;
			}
		}
	}

	//BB_std_free_2d_int(hlpod_ddhr->ovl_elem_global_id, hlpod_ddhr->total_num_elems[0], 1);

/*
	BB_std_free_1d_int(ovl_selected_elems, num_selected_elems);
	BB_std_free_1d_double(ovl_selected_elems_weight, num_selected_elems);
	BB_std_free_1d_int(ovl_selected_elems_D_bc, num_selected_elems_D_bc);
	BB_std_free_1d_double(ovl_selected_elems_weight, num_selected_elems_D_bc);

	BB_std_free_1d_bool(bool_ovl_selected_elems, num_selected_elems);
	BB_std_free_1d_bool(bool_ovl_selected_elems_D_bc, num_selected_elems_D_bc);
	BB_std_free_1d_int(ovl_elem_local_id, num_selected_elems);
	BB_std_free_1d_int(ovl_elem_local_id_D_bc, num_selected_elems_D_bc);
*/

}
/********/

/*for visualization*/
void ddhr_set_selected_elems_para(
		BBFE_DATA*     	fe,
        HLPOD_DDHR*     hlpod_ddhr,
        const int		total_num_nodes,
		const int		num_subdomains,
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

	for(int n=0; n < num_subdomains; n++) {
		for(int m=0; m < hlpod_ddhr->num_selected_elems[n]; m++) {
			int e = hlpod_ddhr->id_selected_elems[m][n];
			for(int i=0; i<nl; i++) {
				int index = fe->conn[e][i];
				
				ECM_elem[index] = myrank + 1;
				ECM_elem_weight[index] += hlpod_ddhr->elem_weight[m][n];
			}
		}

		for(int m=0; m<(hlpod_ddhr->num_selected_elems_D_bc[n]); m++) {
			int e = hlpod_ddhr->id_selected_elems_D_bc[m][n];
			for(int i=0; i<nl; i++) {
				int index = fe->conn[e][i];
				ECM_elem[index] = myrank + 1;
				ECM_elem_weight[index] += hlpod_ddhr->elem_weight_D_bc[m][n];
			}
		}
	}


    for(int m = 0; m < hlpod_ddhr->ovl_num_selected_elems; m++) {
		int e = hlpod_ddhr->ovl_id_selected_elems[m];
		for(int i=0; i<nl; i++) {
			int index = fe->conn[e][i];
			ECM_elem[index] = myrank + 1;
			ECM_elem_weight[index] += hlpod_ddhr->ovl_elem_weight[m];
		}
	}

    for(int m=0; m<(hlpod_ddhr->ovl_num_selected_elems_D_bc); m++) {
        int e = hlpod_ddhr->ovl_id_selected_elems_D_bc[m];
		for(int i=0; i<nl; i++) {
				int index = fe->conn[e][i];
				ECM_elem[index] = myrank + 1;
				ECM_elem_weight[index] += hlpod_ddhr->ovl_elem_weight_D_bc[m];
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
	BB_std_free_1d_double(ECM_elem_weight, fe->total_num_nodes);
	BB_std_free_1d_double(ECM_wireframe, fe->total_num_nodes);

	fclose(fp);

}
/********/

void ddhr_monolis_set_matrix_para(
	MONOLIS*     	monolis,
    HLPOD_VALUES*   hlpod_vals,
	HLPOD_DDHR*     hlpod_ddhr,
	HLPOD_MAT*    hlpod_mat,
	//LPOD_COM*    	lpod_com,
    const int 		max_num_basis,
	const int		num_subdomains)
{
	const int M = max_num_basis;
	const int num_basis = hlpod_vals->num_modes;		//自領域

	const int n_neib_vec = hlpod_vals->n_neib_vec;	//自領域+隣接領域
	//const int M_all = M * (1 + hlpod_vals->n_neib_vec);	//自領域+隣接領域

	const int n_neib = hlpod_vals->n_neib_vec;
	//const int n_internal_vertex = hlpod_mat->n_internal_vertex;

	//printf("M = %d, M_all = %d, num_basis = %d, n_neib_vec = %d", M, M_all, num_basis, n_neib_vec);

	/*add_matrix*/
	//hlpod_mat->L = BB_std_calloc_2d_double(hlpod_mat->L, num_basis , n_neib_vec);
	double** L_in;
	L_in = BB_std_calloc_2d_double(L_in, M , M);

	double t1 = monolis_get_time();

	for(int j = 0; j < M; j++){
		for(int i = 0; i < M; i++){
			L_in[j][i] = 0.0;
		}
	}
	
	for(int j = 0; j < hlpod_vals->num_modes; j++){
		for(int i = 0; i < hlpod_vals->num_modes; i++){
			L_in[i][j] = hlpod_ddhr->reduced_mat[i][j];
		}
	}

	for(int i = hlpod_vals->num_modes; i < M; i++){
		L_in[i][i] = 1.0;
	}

	int* connectivity;
	int* connectivity1;
	int* connectivity2;

	connectivity = BB_std_calloc_1d_int(connectivity, 1);
	connectivity[0] = 0;

	connectivity1 = BB_std_calloc_1d_int(connectivity1, 1);
	connectivity1[0] = 0;
	connectivity2 = BB_std_calloc_1d_int(connectivity2, 1);

	monolis_add_matrix_to_sparse_matrix_R(
		monolis,					//MONOLIS* mat,
  		1,							//int      n_base,
  		connectivity,				//int*     connectivity,
		L_in);						//double** val

	for(int j = 0; j < M; j++){
		for(int i = 0; i < M; i++){
			L_in[j][i] = 0.0;
		}
	}
	int iS = 0;
	int iE = hlpod_mat->num_modes_1stdd_neib[0];
	int index = 0;
	for(int k = 0; k < n_neib; k++){

		iS += hlpod_mat->num_modes_1stdd_neib[k];
		iE += hlpod_mat->num_modes_1stdd_neib[k + 1];
		index = 0;
		for(int j = iS; j < iE; j++){
			for(int i = 0; i < num_basis; i++){
				L_in[i][index] = hlpod_ddhr->reduced_mat[i][j];
			}
			index++;
		}

		connectivity2[0] = k + 1;
		monolis_add_matrix_to_sparse_matrix_offdiag_R(
			monolis,				//MONOLIS* mat,
  			1,						//int      n_base1,
			1,						//int      n_base1,
  			connectivity1,			//int*     connectivity1,
			connectivity2,			//int*     connectivity2,
			L_in);

		for(int j = 0; j < M; j++){
			for(int i = 0; i < M; i++){
				L_in[j][i] = 0.0;
			}
		}
	}

	double t2 = monolis_get_time();
	//lpod_prm->time_add_matrix = t2-t1;


	/*線形の場合*/
  	//hlpod_mat->mode_coef = BB_std_calloc_1d_double(hlpod_mat->mode_coef, M_all);
  	//hlpod_mat->WTf = BB_std_calloc_1d_double(hlpod_mat->WTf, M_all);

	//BB_std_free_2d_double(hlpod_mat->L, num_basis , n_neib_vec);
	BB_std_free_2d_double(L_in, M, M);

	BB_std_free_1d_int(connectivity1, 1);
	BB_std_free_1d_int(connectivity2, 1);
	BB_std_free_1d_int(connectivity, 1);
}


void ddhr_to_monollis_rhs_para(
	MONOLIS*		monolis,
    HLPOD_DDHR*     hlpod_ddrh,
	HLPOD_MAT*    hlpod_mat,
	const int 		k)
{
	for(int i = 0; i < k; i++){
        monolis->mat.R.B[i] = 0.0;
		monolis->mat.R.B[i] = hlpod_ddrh->reduced_RH[i];
	}
}

void lpod_pad_calc_block_solution_local_para(
	MONOLIS_COM*	monolis_com,
	BBFE_DATA* 		fe,
	HLPOD_DDHR*     hlpod_ddhr,
	HLPOD_MAT*    hlpod_mat,
	BBFE_BC*      	bc,
	const int		num_2nddd)
	//LPOD_PRM*		lpod_prm)
{
//	const int n_internal_vertex = hlpod_mat->n_internal_vertex;
    const int n_internal_vertex = monolis_com->n_internal_vertex;

	double t1 = monolis_get_time();

	for(int j = 0; j < fe->total_num_nodes; j++){
		hlpod_ddhr->HR_T[j] = 0.0;
	}

	int index_row = 0;
	int index_column = 0;
	int sum = 0;
	
	for(int k = 0; k < num_2nddd; k++){
		for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
			for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
				index_row = hlpod_mat->node_id[j + sum];
				hlpod_ddhr->HR_T[index_row] += hlpod_mat->pod_modes[index_row][index_column + i] * hlpod_mat->mode_coef[index_column + i];
			}
		}
		index_column += hlpod_mat->num_modes_internal[k];
		sum += hlpod_mat->n_internal_vertex_subd[k];
	}

	for(int i = 0; i < bc->num_D_bcs; i++) {
        int index = 0;
		if(index < n_internal_vertex){
			hlpod_ddhr->HR_T[index] = bc->imposed_D_val[index];
		}
    }

	//解ベクトルのupdate
	monolis_mpi_update_R(monolis_com, fe->total_num_nodes, 1, hlpod_ddhr->HR_T);

	double t2 = monolis_get_time();
}


void ddhr_to_monollis_rhs_para_pad(
	MONOLIS*		monolis,
    HLPOD_DDHR*     hlpod_ddrh,
	HLPOD_MAT*    hlpod_mat,
	const int		num_2nddd,
	const int		max_num_bases)
{
	int index_row = 0;
	int index_column = 0;
	int sum = 0;
	int index = 0;

	for(int k = 0; k < num_2nddd; k++){
		for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
			monolis->mat.R.B[index + i] = hlpod_ddrh->reduced_RH[index_column + i];

            if(monolis_mpi_get_global_my_rank() == 1){
                printf("k = %d, index = %d, index_column = %d, num_modes_internal[k] = %d, reduced_RH[index_column + i] = %lf\n",
                       k, index, index_column, hlpod_mat->num_modes_internal[k], hlpod_ddrh->reduced_RH[index_column + i]);
            }
		}
		index_column += hlpod_mat->num_modes_internal[k];
		index += hlpod_mat->num_modes_internal[k];		
		//index += max_num_bases - hlpod_mat->num_modes_internal[k];
		sum += hlpod_mat->n_internal_vertex_subd[k];
	}

}

void lpod_pad_calc_block_solution_local_para_pad(
	MONOLIS_COM*	monolis_com,
	BBFE_DATA* 		fe,
	HLPOD_DDHR*     hlpod_ddhr,
	HLPOD_MAT*      hlpod_mat,
	BBFE_BC*      	bc,
	const int		num_2nddd,
	const int		max_num_bases)
//	LPOD_PRM*		lpod_prm)
{
	//const int n_internal_vertex = hlpod_mat->n_internal_vertex;
    const int n_internal_vertex = monolis_com->n_internal_vertex;

	double t1 = monolis_get_time();

	for(int j = 0; j < fe->total_num_nodes; j++){
		hlpod_ddhr->HR_T[j] = 0.0;
	}

	int index_row = 0;
	int index_column = 0;
	int sum = 0;
	int index = 0;

	for(int k = 0; k < num_2nddd; k++){
		for(int i = 0; i < hlpod_mat->num_modes_internal[k]; i++){
			for(int j = 0; j < hlpod_mat->n_internal_vertex_subd[k]; j++){
				index_row = hlpod_mat->node_id[j + sum];
				hlpod_ddhr->HR_T[index_row] += hlpod_mat->pod_modes[index_row][index_column + i] * hlpod_mat->mode_coef[index + i];
			}
		}
		index_column += hlpod_mat->num_modes_internal[k];
		index += hlpod_mat->num_modes_internal[k];
		//index += max_num_bases - hlpod_mat->num_modes_internal[k];
		sum += hlpod_mat->n_internal_vertex_subd[k];
	}

	for(int i = 0; i < bc->num_D_bcs; i++) {
        //int index = lpod_prm->D_bc_node_id[i];
        int index = 0;
		if(index < n_internal_vertex){
			hlpod_ddhr->HR_T[index] = bc->imposed_D_val[index];
		}
    }

	//解ベクトルのupdate
	monolis_mpi_update_R(monolis_com, fe->total_num_nodes, 1, hlpod_ddhr->HR_T);

	double t2 = monolis_get_time();
	//lpod_prm->time_calc_sol = t2-t1;
}
