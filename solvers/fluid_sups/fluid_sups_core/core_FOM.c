
#include "core_FOM.h"

const char* ID_NUM_IP_EACH_AXIS = "#num_ip_each_axis";
const int DVAL_NUM_IP_EACH_AXIS = 2;
const char*     ID_MAT_EPSILON  = "#mat_epsilon";
const double  DVAL_MAT_EPSILON  = 1.0e-8;
const char*    ID_MAT_MAX_ITER  = "#mat_max_iter";
const int    DVAL_MAT_MAX_ITER  = 10000;
const char*              ID_DT  = "#time_spacing";
const double           DVAL_DT  = 0.01;
const char*     ID_FINISH_TIME  = "#finish_time";
const double  DVAL_FINISH_TIME  = 1.0;
const char* ID_OUTPUT_INTERVAL  = "#output_interval";
const int DVAL_OUTPUT_INTERVAL  = 1;
const char*         ID_DENSITY  = "#density";
const double      DVAL_DENSITY  = 1000.0;
const char*       ID_VISCOSITY  = "#viscosity";
const double    DVAL_VISCOSITY  = 1.0;

const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND    = "cond.dat";
static const char* INPUT_FILENAME_D_BC_V  = "D_bc_v.dat";
static const char* INPUT_FILENAME_D_BC_P  = "D_bc_p.dat";

static const char* OUTPUT_FILENAME_VTK    = "result_%d_%06d.vtk";
static const char* OUTPUT_FILENAME_CAVITY = "cavity_Re%e.txt";

static const char* OUTPUT_FILENAME_KARMAN_VORTEX_P = "karman_vortex_p_val.dat";
static const char* OUTPUT_FILENAME_KARMAN_VORTEX_N = "karman_vortex_n_val.dat";

static const char* OUTPUT_FILENAME_KARMAN_VORTEX_CP = "karman_vortex_Cp_%d.dat";
static const char* OUTPUT_FILENAME_KARMAN_VORTEX_PINF = "karman_vortex_pinf_%d.dat";

double epsilon = 1.0e-4;

// メッシュは各方向41節点40要素 or 51節点50要素 or 101節点100要素 の六面体一次要素限定
void output_cavity_center_vx(
		VALUES*        vals,
		const char*    method,
		const char*    directory)
{
	double reynolds = vals->density / vals->viscosity;
	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, OUTPUT_FILENAME_CAVITY, method, reynolds);

	FILE* fp;
	fp = ROM_BB_write_fopen(fp, filename, directory);

	for(int i=0; i<51; i++){
	 	double z = 0.02*i;
	 	int n = i*51*51 + 25*51 + 25;
	 	double vx = vals->v[n][0];

		fprintf(fp, "%lf %lf\n", z, vx);
	}

	/*
	for(int i=0; i<101; i++){
		double z = 0.01*i;
		int n = i*101*101 + 50*101 + 50;
		double vx = vals->v[n][0];

		fprintf(fp, "%lf %lf\n", z, vx);
	}
	*/

	fclose(fp);
}


void initialize_velocity_pressure_karman_vortex(
	double** v,
	double* p,
	const int total_num_nodes)
{
    // Initialize velocity array
    for (int i = 0; i < total_num_nodes; i++) {
    //    for (int j = 0; j < 3; j++) {
            v[i][0] = 1.0;
    //    }
    }

    // Initialize pressure array
    for (int i = 0; i < total_num_nodes; i++) {
        p[i] = 0.0;
    }
}

void BBFE_fluid_sups_read_Dirichlet_bc_karman_vortex(
		BBFE_BC*     bc,
		const char*  filename,
		const char*  directory,
		const int    total_num_nodes,
		const int    block_size)
{
	bc->total_num_nodes = total_num_nodes;
	bc->block_size      = block_size;

	srand((unsigned)time(NULL));

	BBFE_sys_memory_allocation_Dirichlet_bc(bc, total_num_nodes, bc->block_size);
	int n = total_num_nodes * bc->block_size;

	for(int i=0; i<n; i++) {
		bc->D_bc_exists[i]   = false;
		bc->imposed_D_val[i] = 0.0;
	}

	FILE* fp;
	fp = BBFE_sys_read_fopen_without_error(fp, filename, directory);
	if( fp == NULL ) {
		printf("%s WARNING: Dirichlet B.C. file, \"%s\", is not found.\n",
				CODENAME, filename);
		return;
	}

	int tmp;
	BB_std_scan_line(&fp, BUFFER_SIZE,
			"%d %d", &(bc->num_D_bcs), &(tmp));
	printf("%s Num. Dirichlet B.C.: %d, Num. block size: %d\n", CODENAME, bc->num_D_bcs, tmp);

	for(int i=0; i<(bc->num_D_bcs); i++) {
		int node_id;  int block_id;  double val;
		BB_std_scan_line(&fp, BUFFER_SIZE,
				"%d %d %lf", &node_id, &block_id, &val);

		int index = (bc->block_size)*node_id + block_id;
		bc->D_bc_exists[ index ]   = true;
		

        if (val == 1.0){
            int r = rand() % 9 + 1;  // 1〜9 の整数を生成
            bc->imposed_D_val[index] = val * (1 + (double)r * epsilon);
        }
	else{
            bc->imposed_D_val[ index ] = val;
        }

	}

	fclose(fp);
}


void output_result_file_karman_vortex(
        BBFE_DATA*     fe,
		VALUES*        vals,
        double         t,
		const char*    directory)
{
	double reynolds = vals->density / vals->viscosity;
	char filename[BUFFER_SIZE];

	snprintf(filename, BUFFER_SIZE, OUTPUT_FILENAME_KARMAN_VORTEX_P);

	FILE* fp;
	fp = BBFE_sys_write_add_fopen(fp, filename, directory);

    double val = 0.02;

    double lowerx = 1.25 - val;
    double upperx = 1.25 + val;

    double lowery = 0.5 - val;
    double uppery = 0.5 + val;

    for (int i = 0; i < fe->total_num_nodes; i++) {
        if (fe->x[i][0] > lowerx && fe->x[i][0] < upperx 
            && fe->x[i][1] > lowery && fe->x[i][1] < uppery) {
            fprintf(fp, "%lf %d %lf %lf %lf %lf %lf\n",
                    t,
                    i,
                    fe->x[i][0],
                    fe->x[i][1],
                    vals->v[i][0],
                    vals->v[i][1],
                    vals->v[i][2]);
        }
    }

	fclose(fp);

	snprintf(filename, BUFFER_SIZE, OUTPUT_FILENAME_KARMAN_VORTEX_N);
	fp = BBFE_sys_write_add_fopen(fp, filename, directory);

    lowery = -0.5 - val;
    uppery = -0.5 + val;

    for (int i = 0; i < fe->total_num_nodes; i++) {
        if (fe->x[i][0] > lowerx && fe->x[i][0] < upperx 
            && fe->x[i][1] > lowery && fe->x[i][1] < uppery) {
            fprintf(fp, "%lf %d %lf %lf %lf %lf %lf\n",
                    t,
                    i,
                    fe->x[i][0],
                    fe->x[i][1],
                    vals->v[i][0],
                    vals->v[i][1],
                    vals->v[i][2]);
        }
    }

fclose(fp);


}

void output_result_file_karman_vortex_pressure(
        BBFE_DATA*     fe,
		VALUES*        vals,
        double         t,
		const char*    directory)
{
	double reynolds = vals->density / vals->viscosity;
	char filename[BUFFER_SIZE];

	snprintf(filename, BUFFER_SIZE, OUTPUT_FILENAME_KARMAN_VORTEX_CP, monolis_mpi_get_global_my_rank());

	FILE* fp;
	fp = BBFE_sys_write_add_fopen(fp, filename, directory);

    double val = 0.000001;

    double lowerr = 0.5 - val;
    double upperr = 0.5 + val;



    for (int i = 0; i < fe->total_num_nodes; i++) {
        double X = fe->x[i][0];
        double Y = fe->x[i][1];
        double r = sqrt(X*X + Y*Y);

        if (r > lowerr && r < upperr) {
            fprintf(fp, "%lf %d %lf %lf %lf\n",
                    t,
                    i,
                    fe->x[i][0],
                    fe->x[i][1],
                    vals->p[i]);
        }
    }

	fclose(fp);
}


void output_result_file_karman_vortex_pressure_inf(
        BBFE_DATA*     fe,
                VALUES*        vals,
        double         t,
                const char*    directory)
{
        double reynolds = vals->density / vals->viscosity;
        char filename[BUFFER_SIZE];

        snprintf(filename, BUFFER_SIZE, OUTPUT_FILENAME_KARMAN_VORTEX_PINF, monolis_mpi_get_global_my_rank());

        FILE* fp;
        fp = BBFE_sys_write_add_fopen(fp, filename, directory);

    double val = 0.1;

    double lowerx = 21.2 - val;
    double upperx = 21.2 + val;

    double lowery = 8.4 - val;
    double uppery = 8.4 + val;

    for (int i = 0; i < fe->total_num_nodes; i++) {
        double X = fe->x[i][0];
        double Y = fe->x[i][1];

        if (X > lowerx && X < upperx &&
			Y > lowery && Y < uppery) {
            fprintf(fp, "%lf %d %lf %lf %lf\n",
                    t,
                    i,
                    fe->x[i][0],
                    fe->x[i][1],
                    vals->p[i]);
        }
    }

        fclose(fp);
}

void memory_allocation_nodal_values(
		VALUES*         vals,
		const int       total_num_nodes)
{
	vals->v = BB_std_calloc_2d_double(vals->v, total_num_nodes, 3);
	vals->p = BB_std_calloc_1d_double(vals->p, total_num_nodes);
}


void assign_default_values(
		VALUES*     vals)
{
	vals->num_ip_each_axis = DVAL_NUM_IP_EACH_AXIS;
	vals->mat_epsilon      = DVAL_MAT_EPSILON;
	vals->mat_max_iter     = DVAL_MAT_MAX_ITER;

	vals->dt               = DVAL_DT;
	vals->finish_time      = DVAL_FINISH_TIME;
	vals->output_interval  = DVAL_OUTPUT_INTERVAL;

	vals->density          = DVAL_DENSITY;
	vals->viscosity        = DVAL_VISCOSITY;
}


void print_all_values(
		VALUES*  vals)
{
	printf("\n%s ---------- Calculation condition ----------\n", CODENAME);

	printf("%s %s: %d\n", CODENAME, ID_NUM_IP_EACH_AXIS, vals->num_ip_each_axis);
	printf("%s %s: %e\n", CODENAME, ID_MAT_EPSILON,      vals->mat_epsilon);
	printf("%s %s: %d\n", CODENAME, ID_MAT_MAX_ITER,     vals->mat_max_iter);

	printf("%s %s: %e\n", CODENAME, ID_DT,               vals->dt);
	printf("%s %s: %e\n", CODENAME, ID_FINISH_TIME,      vals->finish_time);
	printf("%s %s: %d\n", CODENAME, ID_OUTPUT_INTERVAL,  vals->output_interval);

	printf("%s %s: %e\n", CODENAME, ID_DENSITY,          vals->density);
	printf("%s %s: %e\n", CODENAME, ID_VISCOSITY,        vals->viscosity);
	printf("%s -------------------------------------------\n\n", CODENAME);
}


void read_calc_conditions(
		VALUES*     vals,
		const char* directory)
{
	printf("\n");

	assign_default_values(vals);

	char filename[BUFFER_SIZE];
	snprintf(filename, BUFFER_SIZE, "%s/%s", directory, INPUT_FILENAME_COND);

	FILE* fp;
	fp = fopen(filename, "r");
	if( fp == NULL ) {
		printf("%s Calc condition file \"%s\" is not found.\n", CODENAME, filename);
		printf("%s Default values are used in this calculation.\n", CODENAME);
	}
	else {
		printf("%s Reading conditon file \"%s\".\n", CODENAME, filename);
		int num;
		num = BB_std_read_file_get_val_int_p(
				&(vals->num_ip_each_axis), filename, ID_NUM_IP_EACH_AXIS, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->mat_epsilon), filename, ID_MAT_EPSILON, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_int_p(
				&(vals->mat_max_iter), filename, ID_MAT_MAX_ITER, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->dt), filename, ID_DT, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->finish_time), filename, ID_FINISH_TIME, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_int_p(
				&(vals->output_interval), filename, ID_OUTPUT_INTERVAL, BUFFER_SIZE, CODENAME);

		num = BB_std_read_file_get_val_double_p(
				&(vals->density), filename, ID_DENSITY, BUFFER_SIZE, CODENAME);
		num = BB_std_read_file_get_val_double_p(
				&(vals->viscosity), filename, ID_VISCOSITY, BUFFER_SIZE, CODENAME);


		fclose(fp);
	}

	print_all_values(vals);


	printf("\n");
}


void output_result_file_vtk(
		BBFE_DATA*       fe,
		VALUES*        vals,
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
	BB_vtk_write_point_vals_vector(fp, vals->v, fe->total_num_nodes, "Velocity");
	BB_vtk_write_point_vals_scalar(fp, vals->p, fe->total_num_nodes, "Pressure");

	fclose(fp);

}


void output_files(
		FE_SYSTEM* sys,
		const int file_num,
		double t)
{
	int myrank = monolis_mpi_get_global_my_rank();
	char fname_vtk[BUFFER_SIZE];
    
	const char* filename;
	snprintf(fname_vtk, BUFFER_SIZE, OUTPUT_FILENAME_VTK, file_num, myrank);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);
/*
	output_result_file_vtk(
			&(sys->fe), &(sys->vals), filename, sys->cond.directory, t);
*/
}

void BBFE_fluid_sups_read_Dirichlet_bc(
		BBFE_BC*     bc,
		const char*  filename,
		const char*  directory,
		const int    total_num_nodes,
		const int    block_size)
{
	bc->total_num_nodes = total_num_nodes;
	bc->block_size      = block_size;

	BBFE_sys_memory_allocation_Dirichlet_bc(bc, total_num_nodes, bc->block_size);
	int n = total_num_nodes * bc->block_size;

	for(int i=0; i<n; i++) {
		bc->D_bc_exists[i]   = false;
		bc->imposed_D_val[i] = 0.0;
	}

	FILE* fp;
	fp = ROM_BB_read_fopen_without_error(fp, filename, directory);
	if( fp == NULL ) {
		printf("%s WARNING: Dirichlet B.C. file, \"%s\", is not found.\n",
				CODENAME, filename);
		return;
	}

	int tmp;
	BB_std_scan_line(&fp, BUFFER_SIZE,
			"%d %d", &(bc->num_D_bcs), &(tmp));
	printf("%s Num. Dirichlet B.C.: %d, Num. block size: %d\n", CODENAME, bc->num_D_bcs, tmp);

	for(int i=0; i<(bc->num_D_bcs); i++) {
		int node_id;  int block_id;  double val;
		BB_std_scan_line(&fp, BUFFER_SIZE,
				"%d %d %lf", &node_id, &block_id, &val);

		int index = (bc->block_size)*node_id + block_id;
		bc->D_bc_exists[ index ]   = true;
		bc->imposed_D_val[ index ] = val;
	}

	fclose(fp);
}

void set_element_mat(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
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

	for(int e=0; e<(fe->total_num_elems); e++) {
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

						monolis_add_scalar_to_sparse_matrix_R(
								monolis, fe->conn[e][i], fe->conn[e][j], a, b, integ_val);
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


void set_element_vec(
		MONOLIS*     monolis,
		BBFE_DATA*   fe,
		BBFE_BASIS*  basis,
		VALUES*      vals)
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

	for(int e=0; e<(fe->total_num_elems); e++) {
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

				monolis->mat.R.B[ 4*fe->conn[e][i] + d ] += integ_val[d];
			}
		}
	}

	BB_std_free_2d_double(val_ip, 4, np);
	BB_std_free_1d_double(Jacobian_ip, np);
	BB_std_free_2d_double(local_v, nl, 3);
	BB_std_free_2d_double(v_ip, np, 3);
}

//右辺ベクトルのアップデートと解ベクトルの求解
void solver_fom(
    FE_SYSTEM   sys,
    double      t,
    const int   step)
{
		printf("\n%s ----------------- step %d ----------------\n", CODENAME, step);

		monolis_clear_mat_value_R(&(sys.monolis));

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

		monolis_show_timelog (&(sys.monolis), true);
		monolis_show_iterlog (&(sys.monolis), true);
		BBFE_sys_monowrap_solve(
				&(sys.monolis),
				&(sys.mono_com),
				sys.monolis.mat.R.X,
				MONOLIS_ITER_BICGSTAB_N128,
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
			output_files(&sys, step, t);

			//for cavity
			//output_result_file_cavity_center_vx(&(sys.vals), sys.cond.directory);
		}
}



//右辺ベクトルのアップデートと解ベクトルの求解
void solver_fom_collect_snapmat(
    FE_SYSTEM sys,
    double t,
    const int step)
{
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

		monolis_show_timelog (&(sys.monolis), true);
		monolis_show_iterlog (&(sys.monolis), true);
		BBFE_sys_monowrap_solve(
				&(sys.monolis),
				&(sys.mono_com),
				sys.monolis.mat.R.X,
				MONOLIS_ITER_BICGSTAB_N128,
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
			
			output_files(&sys, step, t);
			//for cavity
			//output_result_file_cavity_center_vx(&(sys.vals), sys.cond.directory);
		}

		double* vec;
		vec = BB_std_calloc_1d_double(vec, 3*sys.fe.total_num_nodes);

		ROM_BB_vec_copy_2d_to_1d(
			sys.vals.v,
			vec,
			sys.fe.total_num_nodes);

		if(step%sys.vals.snapshot_interval == 0) {
			printf("set modes p: %d\n", (int)(step/sys.vals.snapshot_interval));

			if(monolis_mpi_get_global_comm_size() == 1){
				ROM_std_hlpod_set_snapmat_nobc(
						sys.vals.p,
						&(sys.rom_p.hlpod_mat),
						sys.fe.total_num_nodes,
                        1,
						(int)(step/sys.vals.snapshot_interval));
			}
			else{
				ROM_std_hlpod_set_snapmat_nobc(
						sys.vals.p,
						&(sys.rom_p.hlpod_mat),
						sys.mono_com.n_internal_vertex,
                        1,
						(int)(step/sys.vals.snapshot_interval));
			}

		}

		if(step%sys.vals.snapshot_interval == 0) {
			printf("set modes v: %d\n", (int)(step/sys.vals.snapshot_interval));

			if(monolis_mpi_get_global_comm_size() == 1){
				ROM_sys_hlpod_fe_set_snap_mat_para(
						vec,
						&(sys.rom_v.hlpod_mat),
						&(sys.bc),
						&(sys.rom_sups.rom_bc),	//要変更
						sys.fe.total_num_nodes,
						3,
						((int)step/sys.vals.snapshot_interval));
			}
			else{
				ROM_sys_hlpod_fe_set_snap_mat_para(
						vec,
						&(sys.rom_v.hlpod_mat),
						&(sys.bc),
						&(sys.rom_sups.rom_bc),
                        sys.mono_com.n_internal_vertex,
						3,
						((int)step/sys.vals.snapshot_interval));
			}

		}

		BB_std_free_1d_double(vec, 3*sys.fe.total_num_nodes);

}
