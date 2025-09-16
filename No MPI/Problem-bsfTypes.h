/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: SAlEM - Stochastic Along Edges Movement method (No MPI)
Module: Problem-bsfTypes.h (Predefined BSF Problem Types)
Prefix: PT_bsf
Authors: Alexander E. Zhulev & Leonid B. Sokolinsky
This source code is a part of BSF Skeleton
==============================================================================*/
#pragma once
#include "Problem-Types.h"		// Problem Types 
//=========================== BSF Types =========================
struct PT_bsf_parameter_T {		// Type of Parameter for workers (current approximation)
	PT_vector_T v_cur;				// Current vertex
};

struct PT_bsf_mapElem_T {		// Type of map-list elements
	int workerNo;
};

struct PT_bsf_reduceElem_T {	// Type of reduce-list elements for Job 0 (default)	
	PT_vector_T v_nex;	// Next vertex
	double objF_nex;	// F(v_nex)
	double objF_grd;	// Value of objective function after one unit movement
	#ifdef PP_MIN_OF_DEGREE
	int numOfEdgeCombinations;		// Number of edge combinations in the next vertex
	#endif // PP_MIN_OF_DEGREE
};

struct PT_bsf_reduceElem_T_1 {	// Type of reduce-list elements for Job 1
	// Not used
};

struct PT_bsf_reduceElem_T_2 {	// Type of reduce-list elements for Job 2
	// Not used
};

struct PT_bsf_reduceElem_T_3 {	// Type of reduce-list elements for Job 3
	// Not used
};