/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: SAlEM - Stochastic Along Edges Movement method (No MPI)
Module: Problem-Types.h (BSF Types)
Prefix: PT
Authors: Alexander E. Zhulev & Leonid B. Sokolinsky
This source code has been produced with using BSF-skeleton
==============================================================================*/			
#pragma once
#include "Problem-Include.h"		// Problem "Include" Files
#include "Problem-Parameters.h"		// Problem Parameters 
//=========================== Problem Types =========================
typedef double PT_matrix_T[PP_MM][PP_N];
typedef double PT_vector_T[PP_N];
typedef double PT_column_T[PP_MM];

typedef char PT_MPS_name_T[PP_MPS_NAME_LENGTH];

struct PT_MPS_row_T {
	char type;
	PT_MPS_name_T name;
	double RHS_value;
};

struct PT_MPS_column_T {
	PT_MPS_name_T varName;
	PT_MPS_name_T rowName;
	double value;
	int j;
};

struct PT_MPS_upBound_T {
	int varIndex;
	double value;
};

struct PT_MPS_fxVariable_T {
	int varIndex;
	double value;
};