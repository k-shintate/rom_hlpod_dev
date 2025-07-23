#pragma once

#include "core_FOM.h"
#include "hlpod_dataset.h"
#include "fluid_sups_dataset.h"

#include "BBFE/manusol/manusol.h"

#include "ecm_write.h"
//#include "hlpod_dataset.h"
//#include "diff_dataset.h"

#include "HR.h"
#include "DDHR.h"
#include "DDHR2.h"
#include "DDHR_para.h"
#include "DDHR_para_lb.h"
#include "set_matvec.h"
#include "set_matvec_NNLS.h"
#include "set_modes.h"
#include "monowrap.h"

#include "core_FOM.h"
#include "core_ROM.h"

/*for Hyper-reduction*/
void HR_output_files(
		FE_SYSTEM*      sys,
		int             file_num,
		double          t);

void HROM_pre_offline(
		FE_SYSTEM* sys,
        ROM*            rom,
        HROM*           hrom,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains);

void HROM_pre_offline2(
		FE_SYSTEM* 		sys,
        ROM*            rom,
        HROM*           hrom,
		const int 		num_modes,
		const int 		num_snapshot,
		const int 		num_2nd_subdomains);

void HROM_pre_online(
		FE_SYSTEM* sys,
        ROM*            rom,
        HROM*           hrom,
		const int num_modes,
		const int num_snapshot,
		const int num_2nd_subdomains);

void HROM_nonparallel(
		FE_SYSTEM       sys,
        ROM*            rom,
        HROM*           hrom,
		const int       step_HR,
		const int       step_POD,
		const double    t);

void HROM_hierarchical_parallel(
		FE_SYSTEM 		sys,
        ROM*            rom,
        HROM*           hrom,
		const int 		step_HR,
		const int 		step_POD,
		const double 	t);

void HROM_pre(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom);

void HROM_memory_allocation(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom);

void HROM_set_matvec(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom,
        int step,
        double t);

void HROM_pre_offline2(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom);

void HROM_memory_allocation_online(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom);

void HROM_pre_online(
        FE_SYSTEM* sys,
        ROM*        rom,
        HROM*       hrom);