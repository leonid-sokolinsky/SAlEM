/*=============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: AlEM - Along Edges Movement Method
Module: _Problems100_1000-0.h (LP problems of dimensions 100...1000 without random inequalities)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
Start vertex *_v.mtx for these problems was calculated by VeSP https://github.com/leonid-sokolinsky/VeSP
LP problems were obtained using LPP-Generator https://github.com/leonid-sokolinsky/LPP-Generator
LP problems are available in https://github.com/leonid-sokolinsky/Set-of-LP-Problems/tree/main/Rnd-LP
===============================================================================*/
#pragma once

//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
#define PP_EPS_RELATIVE_ERROR			1E-11			// Used if defined PP_CHECK_MAX_OBJ_VALUE 
#define PP_MAX_PSEUDOPROJECTING_ITER	100000000		// Maximum acceptable number of iterations in Pseudoprojection on flat

//============================== Problem Parameters ============================
// PP_OBJECTIVE_VECTOR_LENGTH - direct dependence on dimension PD_n.
// P_EPS_ZERO - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH.
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Accuracy for comparison with zero
#define PP_EPS_PROJECTION			(PP_EPS_ZERO*10)		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_PROJECTION*100)	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7					// Length of Objective Vector
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT			100000000	// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//==============================================================================

/*============================== rnd100-0 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd100-0"
#define PP_KK	100		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	101		// Number of equations (number of rows in *.mtx)
#define PP_N	201		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 1009900 // =200*(n-1)*(2+n)/2+100
//-------------------------- Compilation Modes ---------------------------------
#define PP_MAXPROJECTION
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_ITER_COUNT			10				// Each PP_ITER_COUNT-th iteration to be outputted inside PC_bsf_MapF(*)
#define PP_PROJECTION_COUNT		1000000			// Each PP_PROJECTION_COUNT iteration to be outputted inside Flat_MaxProjection(*)
//------------------------------------------------------------------------------
// Elapsed time: 13.555384
// Number of iterations: 5
// Computed objective value: 1009900
// Maximal objective value:  1009900
// Relative error = 0
//------------------------------------------------------------------------------

/*============================== rnd150-0 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd150-0"
#define PP_KK	150		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	151		// Number of equations (number of rows in *.mtx)
#define PP_N	301		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 2264900
//-----------------------------------------------------------------------------

/*============================== rnd200-0 LP problem ==========================*
#define PP_PROBLEM_NAME	"rnd200-0"
#define PP_KK	200		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	201		// Number of equations (number of rows in *.mtx)
#define PP_N	401		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 4019900
//-------------------------- Compilation Modes ---------------------------------
#define PP_MAXPROJECTION
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_ITER_COUNT			10				// Each PP_ITER_COUNT-th iteration to be outputted inside PC_bsf_MapF(*)
#define PP_PROJECTION_COUNT		1000000			// Each PP_PROJECTION_COUNT iteration to be outputted inside Flat_MaxProjection(*)
//-----------------------------------------------------------------------------

/*============================== rnd250-0 LP problem ==========================*
#define PP_PROBLEM_NAME	"rnd250-0"
#define PP_KK	250		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	251		// Number of equations (number of rows in *.mtx)
#define PP_N	501		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 6274900
//-----------------------------------------------------------------------------

/*============================== rnd400-0 LP problem ==========================*
#define PP_PROBLEM_NAME	"rnd400-0"
#define PP_KK	400		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	401		// Number of equations (number of rows in *.mtx)
#define PP_N	801		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 16039900
//-----------------------------------------------------------------------------

/*============================== rnd600-0 LP problem ==========================*
#define PP_PROBLEM_NAME	"rnd600-0"
#define PP_KK	600		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	601		// Number of equations (number of rows in *.mtx)
#define PP_N	1201	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 36059900
//-----------------------------------------------------------------------------

/*============================== rnd800-0 LP problem ==========================*
#define PP_PROBLEM_NAME	"rnd800-0"
#define PP_KK	800		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	801		// Number of equations (number of rows in *.mtx)
#define PP_N	1601	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 64079900
//-----------------------------------------------------------------------------

/*============================== tcube1K LP problem =========================*/
#define PP_PROBLEM_NAME	"tcube1K" // Truncated hypercube
#define PP_KK	1000		// Maximal number of edges that include surface point (compilator limit: 2 147 483 647)
#define PP_M	1001		// Number of equations (number of rows in *.mtx)
#define PP_N	2001	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 100099900
//-----------------------------------------------------------------------------

/*=============================================================================*/