/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: BIP (Block-lterative Projection) method (No MPI)
Module: _Problems-NetLib-LP.h (Problems from the NETLIB LP Test Problem Set)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
LP problems were obtained using BSF-LPP-Generator.
==============================================================================*/
#pragma once

#define PP_SIMPLE_CONVERSION

/*============================== adlittle LP problem ===============*
// Distance to polytope: 9.93450598263140739058175e-06
#define PP_PROBLEM_NAME		"adlittle"
#define PP_M 56		// Number of equations (number of rows in *.mtx)
#define PP_N 138	// Number of variables (number of cols in *.mtx)

//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		1E-5	// Accuracy of belonging to hyperplane
#define PP_EPS_COS					1E-9	// Precision for cos == 1
#define PP_EPS_MOVING				(PP_EPS_ZERO/100)	// Precision for mooving
#define PP_EPS_POINT_INSIDE_CONE	(PP_EPS_ZERO/100)	// Precision for point to be inside cone
//------------------------------------------------------------------/**/

/*============================== afiro LP problem ==================*/
#define PP_PROBLEM_NAME	"afiro"
#define PP_M 27		// Number of equations (number of rows in *.mtx)
#define PP_N 51		// Number of variables (number of cols in *.mtx)
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for MakeHyperplaneList()
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO		// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7			// Length of Objective Vector
//-----------------------------------------------------------------------
#define PP_MAX_OBJ_VALUE 			464.7531
#define PP_PROBE_LENGTH				0.001			// Length of probe shift
//-------------------------------------------------------------------------------

/*============================== beaconfd LP problem ==================*
// Distance to polytope: 9.99994230306940206655895e-06
#define PP_PROBLEM_NAME		"beaconfd"
#define PP_M 173	// Number of equations (number of rows in *.mtx)
#define PP_N 295	// Number of variables (number of cols in *.mtx)
#define PP_SIMPLE_CONVERSION
//---------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		1E-5	// Accuracy of belonging to hyperplane
#define PP_EPS_COS					1E-2	// Precision for cos == 1
#define PP_EPS_MOVING				(PP_EPS_ZERO/100)	// Precision for mooving
#define PP_EPS_POINT_INSIDE_CONE	(PP_EPS_ZERO)	// Precision for point to be inside cone
//------------------------------------------------------------------/**/

/*============================== blend LP problem ==================*
// Distance to polytope: 9.99548788756706405066784e-06
#define PP_PROBLEM_NAME		"blend"
#define PP_M 74		// Number of equations (number of rows in *.mtx)
#define PP_N 114	// Number of variables (number of cols in *.mtx)
//#undef	PP_COS_MODE
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-7	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		1E-5	// Accuracy of belonging to hyperplane
#define PP_EPS_COS					1E-10	// Precision for cos == 1
#define PP_EPS_MOVING				(PP_EPS_ZERO/100)	// Precision for mooving
#define PP_EPS_POINT_INSIDE_CONE	(PP_EPS_ZERO/100)	// Precision for point to be inside cone
//------------------------------------------------------------------/**/

/*============================== fit1d LP problem ==================*
// Zero point is feasible
#define PP_PROBLEM_NAME		"fit1d"
#define PP_M 24		// Number of equations (number of rows in *.mtx)
#define PP_N 1049	// Number of variables (number of cols in *.mtx)
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-7	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		1E-5	// Accuracy of belonging to hyperplane
#define PP_EPS_COS					1E-4	// Precision for cos == 1
#define PP_EPS_MOVING				(PP_EPS_ZERO/100)	// Precision for mooving
#define PP_EPS_POINT_INSIDE_CONE	(PP_EPS_ZERO/100)	// Precision for point to be inside cone
//------------------------------------------------------------------/**/

/*============================== israel LP problem ==================*
// Elapsed time : 167108 sec
// Number of iterations : 1594230556
#define PP_PROBLEM_NAME		"israel"
#define PP_M 174	// Number of equations (number of rows in *.mtx)
#define PP_N 316	// Number of variables (number of cols in *.mtx)
//#undef	PP_COS_MODE
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-6	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		2E-1	// Accuracy of belonging to hyperplane
#define PP_EPS_COS					0		// Precision for cos == 1
#define PP_EPS_MOVING				(PP_EPS_ZERO/100)	// Precision for mooving
#define PP_EPS_POINT_INSIDE_CONE	(PP_EPS_ZERO/100)	// Precision for point to be inside cone
//------------------------------------------------------------------/**/

/*============================== kb2 LP problem ==============================*
// Zero point is feasible
#define PP_PROBLEM_NAME		"kb2"
#define PP_M 43	// Number of equations (number of rows in *.mtx)
#define PP_N 68	// Number of variables (number of cols in *.mtx)
#define PP_COS_MODE
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		1E-5	// Accuracy of belonging to hyperplane
#define PP_EPS_COS					1E-4	// Precision for cos == 1
#define PP_EPS_MOVING				(PP_EPS_ZERO/100)	// Precision for mooving
#define PP_EPS_POINT_INSIDE_CONE	(PP_EPS_ZERO/100)	// Precision for point to be inside cone
//-------------------------------------------------------------------------/**/

/*============================== recipe LP problem ============================*
// Distance to polytope: 9.11030970442950212567909e-07
#define PP_PROBLEM_NAME		"recipe"
#define PP_M 91		// Number of equations (number of rows in *.mtx)
#define PP_N 204	// Number of variables (number of cols in *.mtx)
//#undef	PP_COS_MODE
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		1E-6	// Accuracy of belonging to hyperplane
#define PP_EPS_COS					1E-2	// Precision for cos == 1
#define PP_EPS_MOVING				(PP_EPS_ZERO/100)	// Precision for mooving
#define PP_EPS_POINT_INSIDE_CONE	(PP_EPS_ZERO/100)	// Precision for point to be inside cone
//-------------------------------------------------------------------------/**/

/*============================== sc50a LP problem ==============================*
// Distance to polytope: 9.98581046916478632388689e-08
#define PP_PROBLEM_NAME		"sc50a"
#define PP_M 50	// Number of equations (number of rows in *.mtx)
#define PP_N 78	// Number of variables (number of cols in *.mtx)
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		1E-7	// Accuracy of belonging to hyperplane
#define PP_EPS_COS					1E-3	// Precision for cos == 1
#define PP_EPS_MOVING				(PP_EPS_ZERO/100)	// Precision for mooving
#define PP_EPS_POINT_INSIDE_CONE	(PP_EPS_ZERO/100)	// Precision for point to be inside cone
//----------------------------------------------------------------------------/**/

/*============================== sc50b LP problem ==============================*
// Distance to polytope: 9.95749256169206670294823e-08
#define PP_PROBLEM_NAME		"sc50b"
#define PP_M 50		// Number of equations (number of rows in *.mtx)
#define PP_N 78	// Number of variables (number of cols in *.mtx)
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		1E-7	// Accuracy of belonging to polytope
#define PP_EPS_COS					1E-2	// Precision for cos == 1
#define PP_EPS_MOVING				(PP_EPS_ZERO/100)	// Precision for mooving
#define PP_EPS_POINT_INSIDE_CONE	(PP_EPS_ZERO/100)	// Precision for point to be inside cone
//----------------------------------------------------------------------------/**/

/*============================== sc105 LP problem ==============================*
// Distance to polytope: 9.96588878084987247957247e-08
#define PP_PROBLEM_NAME		"sc105"
#define PP_M 105	// Number of equations (number of rows in *.mtx)
#define PP_N 163	// Number of variables (number of cols in *.mtx)
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		1E-7	// Accuracy of belonging to polytope
#define PP_EPS_COS					1E-2	// Precision for cos == 1
#define PP_EPS_MOVING				(PP_EPS_ZERO/100)	// Precision for mooving
#define PP_EPS_POINT_INSIDE_CONE	(PP_EPS_ZERO/100)	// Precision for point to be inside cone
//----------------------------------------------------------------------------

/*============================== share2b LP problem ==============================*
#define PP_PROBLEM_NAME		"share2b"	
#define PP_M 96		// Number of equations (number of rows in *.mtx)
#define PP_N 162	// Number of variables (number of cols in *.mtx)
//#undef	PP_COS_MODE
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for point to be in halfspace
#define PP_EPS_ON_HYPERPLANE		1E-5	// Accuracy of belonging to polytope
#define PP_EPS_COS					1E-3	// Precision for cos == 1
#define PP_EPS_MOVING				(PP_EPS_ZERO/100)	// Precision for mooving
#define PP_EPS_POINT_INSIDE_CONE	(PP_EPS_ZERO/100)	// Precision for point to be inside cone
//----------------------------------------------------------------------------/**/