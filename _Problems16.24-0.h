/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: SAlEM - Stochastic Along Edges Movement method
Module: _Problems16.24-0.h (LP problems of dimensions 16...24 without random inequalities)
Prefix: PP
Authors: Alexander E. Zhulev & Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
LP problems were obtained using BSF-LPP-Generator.
Initial surface points for these problems were calculated using BSF-Apex-Quest.
==============================================================================*/
#pragma once

//============================== Problem Parameters =============================
// PP_OBJECTIVE_VECTOR_LENGTH - direct dependence on dimension PD_n.
// P_EPS_ZERO - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH.
// PP_EPS_PROJECTION_ROUND - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH. 
//						This parameter affects terminate condition when 
//						calculating pseudoprojection.
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9				// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO			// Precision for point to be in halfspace
#define PP_EPS_MOVING_ON_POLYTOPE	(PP_EPS_ZERO/100)	// Precision for moving on polytope (affects Shift = 0)
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO			// Precision of rounding pseudoprojecting vectors
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+10				// Length of Objective Vector
//-------------------------------------------------------------------------------
#define PP_MAX_PROJECTING_ITER	1E+7	// Maximum acceptable number of iterations in PseudoprojectionOnFace()
#define PP_PROBE_LENGTH			0.1		// Length of probe shift
//=============================================================================

/*============================== rnd16-0 LP problem ===========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd16-0"
#define PP_M	17		// Number of equations (number of rows in *.mtx)
#define PP_N	33		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 27100
//-----------------------------------------------------------------------------

/*============================== rnd17-0 LP problem ===========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd17-0"
#define PP_M	18		// Number of equations (number of rows in *.mtx)
#define PP_N	35		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 30500
//-----------------------------------------------------------------------------

/*============================== rnd18-0 LP problem ===========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd18-0"
#define PP_M	19		// Number of equations (number of rows in *.mtx)
#define PP_N	37		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 34100
//-----------------------------------------------------------------------------

/*============================== rnd19-0 LP problem ===========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd19-0"
#define PP_M	20		// Number of equations (number of rows in *.mtx)
#define PP_N	39		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 37900
//-----------------------------------------------------------------------------

/*============================== rnd20-0 LP problem =========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd20-0"
#define PP_M	21		// Number of equations (number of rows in *.mtx)
#define PP_N	41		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 41900
//-----------------------------------------------------------------------------

/*============================== rnd22-0 LP problem =========================*/
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd22-0"
#define PP_M	23		// Number of equations (number of rows in *.mtx)
#define PP_N	45		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 50500
//-----------------------------------------------------------------------------

/*============================== rnd24-0 LP problem =========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd24-0"
#define PP_M	25		// Number of equations (number of rows in *.mtx)
#define PP_N	49		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 59900
//-----------------------------------------------------------------------------

/*=============================================================================*/