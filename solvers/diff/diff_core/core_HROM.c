
#include "core_FOM.h"
#include "core_ROM.h"
#include "core_HROM.h"

#include "core_FOM.h"
#include "hlpod_dataset.h"
#include "diff_dataset.h"

#include "HR.h"
#include "DDHR.h"
#include "DDHR2.h"
#include "DDHR_para.h"
#include "DDHR_para_lb.h"
//#include "set_matvec.h"
#include "monowrap.h"

const int BUFFER_SIZE = 10000;

static const char* INPUT_FILENAME_COND          = "cond.dat";
static const char* OUTPUT_FILENAME_VTK          = "result_%06d.vtk";
static const char* OUTPUT_FILENAME_ASCII_TEMP   = "temparature_%06d.dat";
static const char* OUTPUT_FILENAME_ASCII_SOURCE = "source_%06d.dat";

void HROM_set_ansvec(
		VALUES*         vals,
	    HLPOD_HR*     	hlpod_hr,
		const int       total_num_nodes)
{
	for(int i = 0; i < total_num_nodes; i++){
		hlpod_hr->HR_T[i] = vals->T[i];
	}
}

void HROM_set_ansvec_para(
		VALUES*         vals,
		HLPOD_DDHR*     hlpod_ddhr,
		const int       total_num_nodes)
{
	for(int i = 0; i < total_num_nodes; i++){
		hlpod_ddhr->HR_T[i] = vals->T[i];
	}
}

/*for Hyper-reduction*/
void HR_output_result_file_vtk(
		BBFE_DATA*     fe,
		VALUES*        vals,
		HLPOD_VALUES*    hlpod_vals,
		HLPOD_HR*      hlpod_hr,		
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
	BB_vtk_write_point_vals_scalar(fp, vals->T, fe->total_num_nodes, "fem-temperature");
	BB_vtk_write_point_vals_scalar(fp, hlpod_vals->sol_vec, fe->total_num_nodes, "pod-temperature");
	BB_vtk_write_point_vals_scalar(fp, hlpod_hr->HR_T, fe->total_num_nodes, "hr-temperature");
	// for manufactured solution
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, vals->T, hlpod_hr->HR_T);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-fem_abs_error");
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, hlpod_vals->sol_vec, hlpod_hr->HR_T);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-pod_abs_error");

	double* source;
	source = BB_std_calloc_1d_double(source, fe->total_num_nodes);
	manusol_set_source(fe, source, t);
	BB_vtk_write_point_vals_scalar(fp, source, fe->total_num_nodes, "source");
	BB_std_free_1d_double(source, fe->total_num_nodes);

	fclose(fp);

}
/********/

/*for Hyper-reduction*/
void HR_output_result_file_vtk_para(
		BBFE_DATA*     fe,
		VALUES*        vals,
		HLPOD_VALUES*    hlpod_vals,
		HLPOD_DDHR*     hlpod_ddhr,		
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
	BB_vtk_write_point_vals_scalar(fp, vals->T, fe->total_num_nodes, "fem-temperature");
	BB_vtk_write_point_vals_scalar(fp, hlpod_vals->sol_vec, fe->total_num_nodes, "pod-temperature");
	BB_vtk_write_point_vals_scalar(fp, hlpod_ddhr->HR_T, fe->total_num_nodes, "hr-temperature");
	// for manufactured solution
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, vals->T, hlpod_ddhr->HR_T);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-fem_abs_error");
	BBFE_manusol_calc_nodal_error_scalar(
			fe, vals->error, hlpod_vals->sol_vec, hlpod_ddhr->HR_T);
	BB_vtk_write_point_vals_scalar(fp, vals->error   , fe->total_num_nodes, "hr-pod_abs_error");

	double* source;
	source = BB_std_calloc_1d_double(source, fe->total_num_nodes);
	manusol_set_source(fe, source, t);
	BB_vtk_write_point_vals_scalar(fp, source, fe->total_num_nodes, "source");
	BB_std_free_1d_double(source, fe->total_num_nodes);

	fclose(fp);

}
/********/



/*for Hyper-reduction*/
void HR_output_files(
		FE_SYSTEM* sys,
		int file_num,
		double t)
{
	const char* filename;
	char fname_vtk[BUFFER_SIZE];
	char fname_tem[BUFFER_SIZE];
	char fname_sou[BUFFER_SIZE];
	snprintf(fname_vtk, BUFFER_SIZE, OUTPUT_FILENAME_VTK, file_num);
	snprintf(fname_tem, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_TEMP, file_num);
	snprintf(fname_sou, BUFFER_SIZE, OUTPUT_FILENAME_ASCII_SOURCE, file_num);
/*
	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_vtk);

	if(monolis_mpi_get_global_comm_size() == 1 &&sys->rom.hlpod_vals.num_2nd_subdomains == 1){
        HR_output_result_file_vtk(
                &(sys->fe),
                &(sys->vals),
                &(sys->rom.hlpod_vals),
                &(sys->hrom.hlpod_hr),
                filename,
                sys->cond.directory,
                t);
    }
    else{  
        HR_output_result_file_vtk_para(
                &(sys->fe),
                &(sys->vals),
                &(sys->rom.hlpod_vals),
                &(sys->hrom.hlpod_ddhr),
                filename,
                sys->cond.directory,
                t);
    }

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_tem);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			sys->rom.hlpod_vals.sol_vec,
			filename,
			sys->cond.directory);
*/

	/**** for manufactured solution ****/
/*
	double* source;
	source = BB_std_calloc_1d_double(source, sys->fe.total_num_nodes);
	manusol_set_source(&(sys->fe), source, t);

	filename = monolis_get_global_output_file_name(MONOLIS_DEFAULT_TOP_DIR, "./", fname_sou);
	BBFE_write_ascii_nodal_vals_scalar(
			&(sys->fe),
			source,
			filename,
			sys->cond.directory);
*/
	double L2_error_fem_hr;
	double L2_error_pod_hr;
	if(monolis_mpi_get_global_comm_size() == 1){
		if(sys->rom.hlpod_vals.num_2nd_subdomains==1){
			L2_error_fem_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->monolis_com),
					t,
					sys->hrom.hlpod_hr.HR_T,
					sys->vals.T);
					
			L2_error_pod_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->monolis_com),
					t,
					sys->rom.hlpod_vals.sol_vec,
					sys->hrom.hlpod_hr.HR_T);
		}
		else{
			L2_error_fem_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->monolis_com),
					t,
					sys->vals.T,
					sys->hrom.hlpod_ddhr.HR_T);
					
			L2_error_pod_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
					&(sys->fe),
					&(sys->basis),
					&(sys->monolis_com),
					t,
					sys->rom.hlpod_vals.sol_vec,
					sys->hrom.hlpod_ddhr.HR_T);
		}
	}
	else{
		L2_error_fem_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
				&(sys->fe),
				&(sys->basis),
				&(sys->monolis_com),
				t,
				sys->vals.T,
				sys->hrom.hlpod_ddhr.HR_T);
				
		L2_error_pod_hr = ROM_sys_hlpod_fe_equivval_relative_L2_error_scalar(
				&(sys->fe),
				&(sys->basis),
				&(sys->monolis_com),
				t,
				sys->rom.hlpod_vals.sol_vec,
				sys->hrom.hlpod_ddhr.HR_T);
	}

	printf("%s L2 error fem-hrom: %e\n", CODENAME, L2_error_fem_hr);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_fem_hrom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_fem_hr);
		fclose(fp);
	}

	printf("%s L2 error rom-hrom: %e\n", CODENAME, L2_error_pod_hr);

	if(monolis_mpi_get_global_my_rank() == 0){
		FILE* fp;
		fp = ROM_BB_write_add_fopen(fp, "l2_error_rom_hrom.txt", sys->cond.directory);
		fprintf(fp, "%e %e\n", t, L2_error_pod_hr);
		fclose(fp);
	}

	//BB_std_free_1d_double(source, sys->fe.total_num_nodes);
	/***********************************/
}
/********/
