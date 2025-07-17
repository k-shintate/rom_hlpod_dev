#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double ROM_BB_sum(
    const int start,
    const int end,
    double* x);

void ROM_BB_vec_copy(
		double*  vec1,	//input
		double*  vec2,	//output
		const int num);

void ROM_BB_vec_copy_2d_to_1d(
		double**  in,
		double*   out,
		const int num);

int ROM_BB_findMax(
    int array[],
    const int size);

void ROM_BB_bubble_sort_with_id(
    int*    vec,
    int*    id,
    int     total_num);

int ROM_BB_binarySearch(
    int*    a,
    int     x, 
    const int N);

int ROM_BB_estimate_num_pod_modes(
    double* V,
    const int max_num_modes,
    const double rom_epsilon);

int*** ROM_BB_calloc_3d_int(
	int***  array,
	const int  size1,
	const int  size2,
	const int  size3);

int ROM_BB_gauss_elimination(
        int n,
        double **A,
        double *b,
        double *x);