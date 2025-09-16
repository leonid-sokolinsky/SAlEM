/*==========================================================================
Project: LiFe - New Linear Programming Solvers
Theme: SAlEM - Stochastic Along Edges Movement method (No MPI)
Module: _Problems-NetLib-LP.h (Problems from the NETLIB LP Test Problem Set)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
Start vertex *_v.mtx for these problems was calculated by VeSP https://github.com/leonid-sokolinsky/VeSP
LP problems are available in https://github.com/leonid-sokolinsky/Set-of-LP-Problems/tree/main/NetLib-LP
============================================================================*/
#pragma once

#define PP_MPS_FORMAT

//=========================== Problem Parameters ===============================
// PP_OBJECTIVE_VECTOR_LENGTH - direct dependence on dimension PD_n.
// P_EPS_ZERO - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH.
//------------------------------------------------------------------------------
#define PP_MAX_PSEUDOPROJECTING_ITER	100000000		// Maximum acceptable number of iterations in Pseudoprojection on flat

/*============================== adlittle LP problem ===========================*
// Number of equations: 15
// Subspace dimension: 82
#define PP_PROBLEM_NAME		"adlittle"
#define PP_M 56	// Number of constraints in mps-file
#define PP_N 97	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE	-225494.96316238038228101176621492
#define PP_BATCH_SIZE		110
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-10	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		2E-5	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
#define PP_EPS_RELATIVE_ERROR		1E-8	// Termination criteria 
#define PP_BATCH_SIZE				110		// Ìaximum number of combinations to process in PC_bsf_MapF
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT			100000000	// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR					1E+3	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
#define PP_RANDOM_LAUNCH_VECTOR 0
//------------------------------------------------------------------------------
// Elapsed time: 0.9869802
// Number of iterations: 3
// Computed objective value: -225494.963164742744993418
// Maximal objective value:  -225494.963162380387075245
// Relative error = 1.05e-11
// Distance to polytope: 2.9723925e-09
//------------------------------------------------------------------------------

/*============================== afiro LP problem ==============================*
// Number of equations : 8
// Subspace dimension : 24
#define PP_PROBLEM_NAME	"afiro"
#define PP_M 27			// Number of constraints in mps-file
#define PP_N 32			// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 464.75314285714285714285714285714
//------------------------------------------------------------------------------
#define PP_EPS_ZERO						1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE			1E-8	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH		1E+6	// Length of Objective Vector
#define PP_EPS_RELATIVE_ERROR			1E-8	// Termination criteria 
#define PP_BATCH_SIZE					50		// Ìaximum number of combinations to process in PC_bsf_MapF
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT				100000000	// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR					1E+6		// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- Compilation Modes ---------------------------------
#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0.0149885
// Number of iterations: 2
// Computed objective value: 464.753142857169962098851
// Maximal objective value:  464.753142857142847788054
// Relative error = 5.83e-14
// Distance to polytope: 1.2876901e-11
//------------------------------------------------------------------------------

/*============================== beaconfd LP problem ===========================*
// Number of equations: 140
// Subspace dimension: 122
#define PP_PROBLEM_NAME		"beaconfd"
#define PP_M 173	// Number of constraints in mps-file
#define PP_N 262	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE -33592.4858072
//------------------------------------------------------------------------------
#define PP_EPS_ZERO						1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE			1.2E-6	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH		1E+6	// Length of Objective Vector
#define PP_EPS_RELATIVE_ERROR			1E-8	// Termination criteria 
#define PP_BATCH_SIZE					300		// Ìaximum number of combinations to process in PC_bsf_MapF
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT				100000000				// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR					1E+3					// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 36.162442
// Number of iterations: 2
// Computed objective value: -33652.428323099302360788
// Maximal objective value:  -33592.4858071999988169409
// Relative error = 0.00178
// Distance to polytope: 4.2262245e-07
//------------------------------------------------------------------------------

/*============================== blend LP problem ==============================*
// Number of equations: 43
// Subspace dimension: 40
#define PP_PROBLEM_NAME		"blend"
#define PP_M 74			// Number of constraints in mps-file
#define PP_N 83			// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 30.812149845828220173774356124984	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-5	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
#define PP_BATCH_SIZE				10000		// Ìaximum number of combinations to process in PC_bsf_MapF
#define PP_EPS_RELATIVE_ERROR		1E-8	// Termination criteria 
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT				10000000			// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//-------------------------- Compilation Modes ---------------------------------
//#define PP_MAXPROJECTION
//#define PP_GRADIENT
//#define PP_MEDIAN_LAUNCH_VECTOR
//#define PP_RANDOM_LAUNCH_VECTOR 283
//------------------------------------------------------------------------------
// Elapsed time: 18.54392
// Number of iterations: 1
// Computed objective value: 0
// Maximal objective value:  30.8121498458282196963864
// Relative error = 1
// Distance to polytope: 0
//------------------------------------------------------------------------------

/*============================== fit1d LP problem ==============================*
// Number of equations : 1
// Subspace dimension : 1025
#define PP_PROBLEM_NAME		"fit1d"
#define PP_M 24		// Number of constraints 
#define PP_N 1026	// Number of variables
#define PP_MAX_OBJ_VALUE 9146.3780924209269467749025024617	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-10	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-7	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
#define PP_BATCH_SIZE				30		// Ìaximum number of combinations to process in PC_bsf_MapF
#define PP_EPS_RELATIVE_ERROR		1E-11
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT				1000000				// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR					1E+5					// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//#define PP_MAXPROJECTION
//#define PP_MEDIAN_LAUNCH_VECTOR
//#define PP_RANDOM_LAUNCH_VECTOR			283
//------------------------------------------------------------------------------
// ?
//------------------------------------------------------------------------------

/*============================== grow7 LP problem ============================*
// Number of equations: 140
// Subspace dimension: 161
#define PP_PROBLEM_NAME		"grow7"
#define PP_M 140	// Number of equations (after conversion to standard form)
#define PP_N 301	// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 47787811.814711502616766956242865	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO						1E-11	// Accuracy for comparison with zero
//#define PP_EPS_PROJECTION				1E-10	// Precision to calculate projection
#define PP_EPS_ON_HYPERPLANE			1E-2	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH		1E+7	// Length of Objective Vector
#define PP_EPS_RELATIVE_ERROR			1E-11				 
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT				50000000				// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR					1E+1					// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- Compilation Modes ---------------------------------
//#define PP_MAXPROJECTION
//#define PP_MEDIAN_LAUNCH_VECTOR
//------------------------------------------------------------------------------
// Elapsed time: 12.21888
// Number of iterations: 11
// Computed objective value: 47787811.8147406205534935
// Maximal objective value:  47787811.8147115036845207
// Relative error = 6.09e-13
// Distance to polytope: 5.3484396e-06
//------------------------------------------------------------------------------

/*============================== israel LP problem =========================*
// Number of equations: 0
#define PP_PROBLEM_NAME		"israel"
#define PP_M 174	// Number of constraints in mps-file
#define PP_N 142	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 896644.82186304572966200464196045	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-8	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		2E-3	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
#define PP_EPS_RELATIVE_ERROR		1E-9	// Termination criteria 
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT			50000000				// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR				1E+3					// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- Compilation Modes ---------------------------------
#define PP_RANDOM_LAUNCH_VECTOR 11
//------------------------------------------------------------------------------
// Elapsed time: 14.481201
// Number of iterations: 5
// Computed objective value: 896644.928349998430348933
// Maximal objective value:  896644.821863045683130622
// Relative error = 1.19e-07
// Distance to polytope: 0.00028764043
//------------------------------------------------------------------------------

/*============================== kb2 LP problem ================================*
// Number of equations: 16
// Subspace dimension: 25
#define PP_PROBLEM_NAME		"kb2"
#define PP_M 43	// Number of equations (after conversion to standard form)
#define PP_N 41	// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 1749.9001299062057129526866493726
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-10	// Accuracy for comparison with zero
//#define PP_EPS_PROJECTION			1E-7	// Accuracy of calculating pseudoprojection
#define PP_EPS_ON_HYPERPLANE		2E-5	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
#define PP_EPS_RELATIVE_ERROR		1E-8	// Termination criteria 
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT				100000000			// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR					1E+5					// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- Compilation Modes ---------------------------------
//#define PP_MAXPROJECTION
//------------------------------------------------------------------------------
// Elapsed time: 0.0109538
// Number of iterations: 16
// Computed objective value: 1749.90012192081326247717
// Maximal objective value:  1749.90012990620562050026
// Relative error = 4.56e-09
// Distance to polytope: 3.4918643e-06
//------------------------------------------------------------------------------

/*============================== recipe LP problem =============================*
// Number of equations: 67
// Subspace dimension: 92 
#define PP_PROBLEM_NAME		"recipe"
#define PP_M 91	// Number of constraints in mps-file
#define PP_N 180	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 266.616 // Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-8	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
#define PP_BATCH_SIZE				10		// Ìaximum number of combinations to process in PC_bsf_MapF
#define PP_EPS_RELATIVE_ERROR		1E-8	// Termination criteria 
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT			100000000	// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//-------------------------- Compilation Modes ---------------------------------
//------------------------------------------------------------------------------
// Elapsed time: 0.1891606
// Number of iterations: 1
// Computed objective value: 266.615999996835569163522
// Maximal objective value:  266.615999999999985448085
// Relative error = 1.19e-11
// Distance to polytope: 9.9162087e-13
//------------------------------------------------------------------------------

/*============================== sc105 LP problem ==============================*
// Number of equations: 45
// Subspace dimension: 58
#define PP_PROBLEM_NAME		"sc105"
#define PP_M 104	// Number of constraints in mps-file
#define PP_N 103	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 52.202061211707248062628010857689 // Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-6	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
#define PP_BATCH_SIZE				4000		// Ìaximum number of combinations to process in PC_bsf_MapF
#define PP_EPS_RELATIVE_ERROR		1E-8	// Termination criteria
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT			100000000	// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR				1E+7					// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 27.492867
// Number of iterations: 2
// Computed objective value: 52.2020612131678589662442
// Maximal objective value:  52.2020612117072460023337
// Relative error = 2.8e-11
// Distance to polytope: 1.2629975e-09
//------------------------------------------------------------------------------

/*============================== sc50a LP problem ==============================*
// Number of equations: 20
// Subspace dimension: 28
#define PP_PROBLEM_NAME		"sc50a"
#define PP_M 49	// Number of constraints
#define PP_N 48	// Number of variables
#define PP_MAX_OBJ_VALUE 64.575077058564509026860413914575	// Exact maximum value of objective function
//-------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-4	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
#define PP_EPS_RELATIVE_ERROR		1E-8	// Termination criteria 
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT				100000000			// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR				1E+7					// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//----------------------------------------------------------------------------
// Number of iterations: 4
// Computed objective value: 64.5750770586504359016544
// Maximal objective value:  64.5750770585645028631916
// Relative error = 1.33e-12
//----------------------------------------------------------------------------

/*============================== sc50b LP problem ============================*
// Number of equations: 20
// Subspace dimension: 28
#define PP_PROBLEM_NAME		"sc50b"
#define PP_M 48	// Number of constraints
#define PP_N 48	// Number of variables
#define PP_MAX_OBJ_VALUE 70	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-12	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-7	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
#define PP_EPS_RELATIVE_ERROR		1E-8	// Termination criteria 
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR				1E+7					// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT				100000000			// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//--------------------------------------------------------------------------
// Elapsed time: 3.747997
// Number of iterations: 1
// Computed objective value: 70.0000000000000142108547
// Maximal objective value:  70
// Relative error = 2.03e-16
//--------------------------------------------------------------------------

/*============================== share2b LP problem ============================*
// Number of equations: 13
// Subspace dimension: 66
#define PP_PROBLEM_NAME		"share2b"
#define PP_M 96		// Number of constraints in *.mps
#define PP_N 79	// Number of variables in *.mps
#define PP_MAX_OBJ_VALUE 415.732240741419486545199108738 // Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-4	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
#define PP_EPS_RELATIVE_ERROR		1E-4	// Termination criteria 
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR				1E+6	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT				10000000			// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0.0594301
// Number of iterations: 6
// Computed objective value: 415.732255369600750327663
// Maximal objective value:  415.732240741419502683129
// Relative error = 3.52e-08
// Distance to polytope: 1.8331818e-07
//------------------------------------------------------------------------------

/*============================== scagr7 LP problem =============================*
// Number of equations : 84
// Subspace dimension : 56
#define PP_PROBLEM_NAME	"scagr7"
#define PP_M 129		// Number of constraints in mps-file
#define PP_N 140		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 2331389.824330984	// Exact maximum value of objective function
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		2E-5	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
#define PP_EPS_RELATIVE_ERROR		1E-8	// Termination criteria 
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR				1E+6	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT				10000000			// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//-------------------------- Compilation Modes ---------------------------------
//#define PP_RANDOM_LAUNCH_VECTOR 113
#define PP_MEDIAN_LAUNCH_VECTOR
//------------------------------------------------------------------------------
// Elapsed time: 12.0677
// Number of iterations: 19
// Computed objective value: 2331389.82462500780820847
// Maximal objective value:  2331389.8243309841491282
// Relative error = 1.26e-10
// Distance to polytope: 1.5948626e-07
//------------------------------------------------------------------------------

/*============================== stocfor1 LP problem ============================*/
// Number of equations: 63
// Subspace dimension: 48
#define PP_PROBLEM_NAME		"stocfor1"	
#define PP_M 117	// Number of constraints in mps-file
#define PP_N 111	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 41131.976219436406065682760731514 // Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-3	// Precision for point to be in halfspace
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
#define PP_BATCH_SIZE				10000		// Ìaximum number of combinations to process in PC_bsf_MapF
#define PP_EPS_RELATIVE_ERROR		1E-4	// Termination criteria 
//------------------------------ ifdef PP_SAVE_ITER_RESULT ---------------------
#define PP_SCALE_FACTOR				1E+6	// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT			100000000			// Each PP_PROJECTION_COUNT-th iteration to be outputted inside Flat_MaxProjection(*) or Flat_BipProjection(*)
//-------------------------- Compilation Modes ---------------------------------
//#define PP_RANDOM_LAUNCH_VECTOR
//------------------------------------------------------------------------------
// Elapsed time: 41.078658
// Number of iterations: 1
// Computed objective value: 25468.6697284531473997049
// Maximal objective value:  41131.9762194364084280096
// Relative error = 0.381
// Distance to polytope: 1.4961e-12
//------------------------------------------------------------------------------

//==============================================================================*/