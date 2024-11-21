/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: SAlEM - Stochastic Along Edges Movement method
Module: _Problems200_1000-0.h (LP problems of dimensions 200...1000 without random inequalities)
Prefix: PP
Authors: Alexander E. Zhulev & Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
LP problems were obtained using BSF-LPP-Generator.
==============================================================================*/
#pragma once

//============================== Problem Parameters =============================
// PP_OBJECTIVE_VECTOR_LENGTH - direct dependence on dimension PD_n.
// P_EPS_ZERO - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH.
// PP_EPS_PROJECTION_ROUND - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH. 
//						This parameter affects terminate condition when 
//						calculating pseudoprojection.
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	(PP_EPS_ZERO*10)// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10)// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO		// Precision of rounding pseudoprojecting vectors
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+10			// Length of Objective Vector
//-------------------------------------------------------------------------------
#define PP_MAX_PROJECTING_ITER	1E+7	// Maximum acceptable number of iterations in PseudoprojectionOnFace()
#define PP_PROBE_LENGTH			0.001		// Length of probe shift
//=============================================================================

/*============================== rnd200-0 LP problem =========================*/
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd200-0"
#define PP_M	201		// Number of equations (number of rows in *.mtx)
#define PP_N	401		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE	3010000
//-------------------------------------------------------------------------------
#define PP_REAL_TIME 100 // This parameter limits the calculation time (compilator limit: 2 147 483 647)
//-------------------------------------------------------------------------------

/*============================== rnd400-0 LP problem =========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd400-0"
#define PP_M	401		// Number of equations (number of rows in *.mtx)
#define PP_N	801		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 12060000
//-----------------------------------------------------------------------------

/*============================== rnd600-0 LP problem =========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd600-0"
#define PP_M	601		// Number of equations (number of rows in *.mtx)
#define PP_N	1201	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 27030000
//-----------------------------------------------------------------------------

/*============================== rnd800-0 LP problem =========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd800-0"
#define PP_M	801		// Number of equations (number of rows in *.mtx)
#define PP_N	1601	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 48040000
//-----------------------------------------------------------------------------

/*============================== rnd1000-0 LP problem =========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd1000-0"
#define PP_M	1001		// Number of equations (number of rows in *.mtx)
#define PP_N	2001	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 75050000
//-----------------------------------------------------------------------------

/*=============================================================================*/