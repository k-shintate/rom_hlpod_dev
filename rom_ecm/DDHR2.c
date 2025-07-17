//DDHROMに関して、オーバーラップ要素を含んで計算する方式
//内部要素の総和が分割前の要素数の総和になることを利用

/*
#include "std.h"
#include "pod_dataset.h"
#include "hlpod_write_fe.h"
#include "hlpod_read.h"
#include "hlpod_write.h"
*/

#include "hlpod_core_fe.h"

static const int BUFFER_SIZE = 10000;
static const char* INPUT_FILENAME_ELEM_ID          = "elem.dat.id";
static const char* OUTPUT_FILENAME_ECM_ELEM_VTK = "ECM_elem.vtk";

void ddhr_memory_allocation2(
        const int       total_num_nodes,
        const int       total_num_elem,
        const int       total_num_snapshot,
        const int       total_num_modes,
		const int		num_subdomains,
        HLPOD_DDHR*     hlpod_ddhr)
{
	printf("\n%d\n", num_subdomains);

	int max_num_elem = findMax(hlpod_ddhr->num_elems, num_subdomains);

    hlpod_ddhr->HR_T = BB_std_calloc_1d_double(hlpod_ddhr->HR_T, total_num_nodes);

//for NNLS
    //hlpod_ddhr->matrix = BB_std_calloc_3d_double(hlpod_ddhr->matrix, total_num_snapshot*total_num_modes*num_subdomains*2, max_num_elem, num_subdomains);
    //hlpod_ddhr->RH = BB_std_calloc_2d_double(hlpod_ddhr->RH, total_num_snapshot*total_num_modes*num_subdomains*2, num_subdomains);

    hlpod_ddhr->reduced_mat = BB_std_calloc_2d_double(hlpod_ddhr->reduced_mat, total_num_modes*num_subdomains, total_num_modes*num_subdomains);
    hlpod_ddhr->reduced_RH = BB_std_calloc_1d_double(hlpod_ddhr->reduced_RH, total_num_modes*num_subdomains);

}
/********/

void ddhr_set_element2(
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
        snprintf(fname, BUFFER_SIZE, "parted.0/elem.dat.n_internal.%d", m);

        fp = BBFE_sys_hlpod_read_fopen(fp, fname, directory);
        fscanf(fp, "%s %d", char_id, &(tmp));
        fscanf(fp, "%d", &(num_elems[m]));
		printf("num_elems = %d\n", num_elems[m]);
        fclose(fp);
	}

	int max_num_elem = findMax(num_elems, num_subdomains);

	hlpod_ddhr->elem_id_local = BB_std_calloc_2d_int(hlpod_ddhr->elem_id_local, max_num_elem, num_subdomains);
	hlpod_ddhr->num_elems = BB_std_calloc_1d_int(hlpod_ddhr->num_elems, num_subdomains);

    for(int m = 0; m < num_subdomains; m++){
		hlpod_ddhr->num_elems[m] = num_elems[m];
	}

    for(int m = 0; m < num_subdomains; m++){	
        snprintf(fname, BUFFER_SIZE, "parted.0/%s.%d", INPUT_FILENAME_ELEM_ID, m);

        fp = BBFE_sys_hlpod_read_fopen(fp, fname, directory);
        fscanf(fp, "%s", char_id);
        fscanf(fp, "%d %d", &(num_elems[m]), &(tmp));

        for(int i = 0; i < hlpod_ddhr->num_elems[m]; i++) {
            fscanf(fp, "%d", &(hlpod_ddhr->elem_id_local[i][m]));
        }

        fclose(fp);
	}

}

void ddhr_lb_read_selected_elements(
    HLPOD_DDHR*     hlpod_ddhr,
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

//	fp3 = BBFE_sys_hlpod_write_fopen(fp3, fname3, directory);
//	fp4 = BBFE_sys_hlpod_write_fopen(fp4, fname4, directory);

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

		fp1 = BBFE_sys_hlpod_read_fopen(fp1, fname1, directory);
		fp2 = BBFE_sys_hlpod_read_fopen(fp2, fname2, directory);

		fscanf(fp1, "%d", &(num_selected_elems));
		fscanf(fp2, "%d", &(num_selected_elems_D_bc));
		Index1 += num_selected_elems;
		Index2 += num_selected_elems_D_bc;

		fclose(fp1);
		fclose(fp2);
	}

//    hlpod_ddhr->elem_weight = BB_std_calloc_1d_double(f_ip, np);
    hlpod_ddhr->ovl_num_selected_elems = Index1;
    hlpod_ddhr->ovl_num_selected_elems_D_bc  = Index2;

    hlpod_ddhr->ovl_elem_weight = BB_std_calloc_1d_double(hlpod_ddhr->ovl_elem_weight, Index1);
    hlpod_ddhr->ovl_id_selected_elems = BB_std_calloc_1d_int(hlpod_ddhr->ovl_id_selected_elems, Index1);
    hlpod_ddhr->ovl_elem_weight_D_bc = BB_std_calloc_1d_double(hlpod_ddhr->ovl_elem_weight_D_bc, Index2);
    hlpod_ddhr->ovl_id_selected_elems_D_bc = BB_std_calloc_1d_int(hlpod_ddhr->ovl_id_selected_elems_D_bc, Index2);

	for (int m = 0; m < num_subdomains; m++) {
		snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", m);

		fp1 = BBFE_sys_hlpod_read_fopen(fp1, fname1, directory);

		fscanf(fp1, "%d", &(num_selected_elems));
        for(int i = 0; i < num_selected_elems; i++){
            fscanf(fp1, "%d %lf", &(hlpod_ddhr->ovl_id_selected_elems[index1]), &(hlpod_ddhr->ovl_elem_weight[index1]));
            index1++;
        }

		fclose(fp1);
	}

    for (int m = 0; m < num_subdomains; m++) {
		snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", m);

		fp2 = BBFE_sys_hlpod_read_fopen(fp2, fname2, directory);

		fscanf(fp2, "%d", &(num_selected_elems_D_bc));

        for(int i = 0; i < num_selected_elems_D_bc; i++){
            fscanf(fp2, "%d %lf", &(hlpod_ddhr->ovl_id_selected_elems_D_bc[index2]), &(hlpod_ddhr->ovl_elem_weight_D_bc[index2]));
            
            printf("%d %lf\n", hlpod_ddhr->ovl_id_selected_elems_D_bc[index2], hlpod_ddhr->ovl_elem_weight_D_bc[index2]);
            index2++;
        }

		fclose(fp2);
	}


/*
	fprintf(fp3, "%d\n", Index1);
	fprintf(fp4, "%d\n", Index2);

	for (int m = 0; m < num_subdomains; m++) {
		snprintf(fname1, BUFFER_SIZE, "DDECM/lb_selected_elem_D_bc.%d.txt", subdomain_id[m]);
		snprintf(fname2, BUFFER_SIZE, "DDECM/lb_selected_elem.%d.txt", subdomain_id[m]);

		fp1 = BBFE_sys_hlpod_read_fopen(fp1, fname1, directory);
		fp2 = BBFE_sys_hlpod_read_fopen(fp2, fname2, directory);

		fscanf(fp1, "%d", &(num_selected_elems));
		fscanf(fp2, "%d", &(num_selected_elems_D_bc));

		for (int i = 0; i < num_selected_elems; i++) {
			fscanf(fp1, "%d %lf", &(tmp), &(val));
			fprintf(fp3, "%d %.30e\n", tmp, val);
			index1++;
		}

		for (int i = 0; i < num_selected_elems_D_bc; i++) {
			fscanf(fp2, "%d %lf", &(tmp), &(val));
			fprintf(fp4, "%d %.30e\n", tmp, val);
			index2++;
		}

		fclose(fp1);
		fclose(fp2);
	}

	fclose(fp3);
	fclose(fp4);
*/
//	t = monolis_get_time_global_sync();
}


//level1領域の選択された基底(p-adaptive)本数の共有
void get_neib_subdomain_id_nonpara(
    POD_MATRIX* 	pod_mat,
	HLPOD_DDHR* 	hlpod_ddhr,
    const int       num_subdomains)
{
/*
	lpod_mat->num_neib_modes = BB_std_calloc_1d_int(lpod_mat->num_neib_modes, np);
	lpod_mat->num_neib_modes[0] = num_my_modes;
	monolis_mpi_update_I(monolis_com, np, 1, lpod_mat->num_neib_modes);
	lpod_mat->num_neib_modes_sum = BB_std_calloc_1d_int(lpod_mat->num_neib_modes_sum, np);	
	lpod_mat->num_neib_modes_sum[0] = num_my_modes;
	for(int i = 1; i < np; i++){
		lpod_mat->num_neib_modes_sum[i] = lpod_mat->num_neib_modes_sum[i-1] + lpod_mat->num_neib_modes[i];
	}
	lpod_mat->num_neib_modes = BB_std_calloc_1d_int(lpod_mat->num_neib_modes, np);
	lpod_mat->num_neib_modes[0] = num_my_modes;
	monolis_mpi_update_I(monolis_com, np, 1, lpod_mat->num_neib_modes);
*/

	hlpod_ddhr->num_neib_modes_1stdd_sum = BB_std_calloc_1d_int(hlpod_ddhr->num_neib_modes_1stdd_sum, num_subdomains +1);
	hlpod_ddhr->num_neib_modes_1stdd_sum[0] = 0;
	for(int i = 1; i < num_subdomains + 1; i++){
		hlpod_ddhr->num_neib_modes_1stdd_sum[i] = hlpod_ddhr->num_neib_modes_1stdd_sum[i-1] + pod_mat->num_modes_internal[i-1];
	}
}


void ddhr_get_selected_elements2(
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

    int NNLS_row = total_num_snapshot*total_num_modes*num_subdomains*2;	//2は残差ベクトル＋右辺ベクトルを採用しているため

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

        hr_write_NNLS_residual(residual, 0, m, directory);
        hr_write_NNLS_num_elems(total_num_selected_elems, 0, m, directory);

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

        //index = D_bcが付与された要素数

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

	int max_num_elem = findMax(hlpod_ddhr->num_elems, num_subdomains);
	BB_std_free_3d_double(hlpod_ddhr->matrix, NNLS_row, max_num_elem, num_subdomains);
    BB_std_free_2d_double(hlpod_ddhr->RH, NNLS_row, num_subdomains);


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
/********/


void ddhr_set_selected_elems(
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
				
				ECM_elem[index] = n + 1;
				ECM_elem_weight[index] += hlpod_ddhr->elem_weight[m][n];
			}
		}

		for(int m=0; m<(hlpod_ddhr->num_selected_elems_D_bc[n]); m++) {
			int e = hlpod_ddhr->id_selected_elems_D_bc[m][n];
			for(int i=0; i<nl; i++) {
				int index = fe->conn[e][i];
				ECM_elem[index] = n + 1;
				ECM_elem_weight[index] += hlpod_ddhr->elem_weight_D_bc[m][n];
			}
		}
	}

	const char* filename;

	char fname[BUFFER_SIZE];

	snprintf(fname, BUFFER_SIZE,"DDECM/%s", OUTPUT_FILENAME_ECM_ELEM_VTK);
	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname);

	FILE* fp;
	fp = BBFE_sys_hlpod_write_fopen(fp, filename, directory);

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