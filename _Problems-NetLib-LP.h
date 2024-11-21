/*==========================================================================
Project: LiFe - New Linear Programming Solvers
Theme: BIP (Block-lterative Projection) method (No MPI)
Module: _Problems-NetLib-LP.h (Problems from the NETLIB LP Test Problem Set)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
============================================================================*/
#pragma once

#define PP_MPS_FORMAT

/*============================== adlittle LP problem =======================*
// Number of equations : 15
// Subspace dimension : 82
#define PP_PROBLEM_NAME		"adlittle"
#define PP_M 56	// Number of constraints in mps-file
#define PP_N 97	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 		-225494.96316238038228101176621492
//--------------------------------------------------------------------------
#define PP_REAL_TIME 1000 // This parameter limits the calculation time (compilator limit: 2 147 483 647)
//--------------------------------------------------------------------------

/*============================== afiro LP problem ==========================*
// Number of equations : 8
// Subspace dimension : 24
#define PP_PROBLEM_NAME	"afiro"
#define PP_M 27		// Number of constraints in mps-file
#define PP_N 32		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 464.75314285714285714285714285714
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-7			// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	(PP_EPS_ZERO*10)// Precision for MakeHyperplaneList()
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10)// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO		// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7			// Length of Objective Vector
#define PP_PROBE_LENGTH				0.001			// Length of probe shift
#define PP_REAL_TIME				300				// This parameter limits the calculation time (compilator limit: 2 147 483 647)
#define PP_MAX_PSEUDOPROJECTING_ITER 10000			// Maximum number of iterations to interrupt calculation of pseudoprojection on flat
//--------------------------------------------------------------------------
// Elapsed time: 1.2636478
// Number of iterations: 1
// Computed objective value: 464.7531543530281
// Maximal objective value:  464.7531428571428
// Relative error = 2.47e-08
//--------------------------------------------------------------------------

/*============================== beaconfd LP problem =======================*/
// Number of equations: 140
// Subspace dimension: 122
#define PP_PROBLEM_NAME		"beaconfd"
#define PP_M 173	// Number of constraints in mps-file
#define PP_N 262	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE -33592.4858072
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-7			// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	(PP_EPS_ZERO*10)// Precision for MakeHyperplaneList()
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10)// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO		// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7			// Length of Objective Vector
#define PP_PROBE_LENGTH				0.001			// Length of probe shift
#define PP_REAL_TIME				600			// This parameter limits the calculation time (compilator limit: 2 147 483 647)
#define PP_MAX_PSEUDOPROJECTING_ITER 100000			// Maximum number of iterations to interrupt calculation of pseudoprojection on flat
//--------------------------------------------------------------------------
// Elapsed time: 9608.0186
// Number of iterations: 0
// Computed objective value: -33837.75644768299
// Maximal objective value:  -33592.4858072
// Relative error = 0.0073
// //--------------------------------------------------------------------------

/*============================== blend LP problem ==========================*
// Number of equations: 43
// Subspace dimension: 40
#define PP_PROBLEM_NAME		"blend"
#define PP_M 74	// Number of constraints in mps-file
#define PP_N 83		// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 30.812149845828220173774356124984	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-3			// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	(PP_EPS_ZERO*10)// Precision for MakeHyperplaneList()
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10)// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO		// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7			// Length of Objective Vector
#define PP_PROBE_LENGTH				0.001			// Length of probe shift
#define PP_REAL_TIME				10000			// This parameter limits the calculation time (compilator limit: 2 147 483 647)
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------

/*============================== fit1d LP problem ==========================*
// Number of equations : 1
// Subspace dimension : 1025
#define PP_PROBLEM_NAME		"fit1d"
#define PP_M 24	// Number of equations (after conversion to standard form)
#define PP_N 1026	// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 9146.3780924209269467749025024617	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	(PP_EPS_ZERO*10)// Precision for MakeHyperplaneList()
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10)// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO		// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7			// Length of Objective Vector
#define PP_PROBE_LENGTH				0.001			// Length of probe shift
#define PP_REAL_TIME				10		// This parameter limits the calculation time (compilator limit: 2 147 483 647)
//--------------------------------------------------------------------------

/*============================== israel LP problem =========================*
// Number of equations: 0
#define PP_PROBLEM_NAME		"israel"
#define PP_M 174	// Number of constraints in mps-file
#define PP_N 142	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 896644.82186304572966200464196045	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	(PP_EPS_ZERO*10)// Precision for MakeHyperplaneList()
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10)// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO		// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7			// Length of Objective Vector
#define PP_PROBE_LENGTH				0.001			// Length of probe shift
#define PP_REAL_TIME				1		// This parameter limits the calculation time (compilator limit: 2 147 483 647)
//--------------------------------------------------------------------------

/*============================== kb2 LP problem ============================*
// Number of equations: 16
// Subspace dimension: 25
#define PP_PROBLEM_NAME		"kb2"
#define PP_M 43	// Number of equations (after conversion to standard form)
#define PP_N 41	// Number of variables in mps-file (after conversion to standard form)
#define PP_MAX_OBJ_VALUE 1749.9001299062057129526866493726
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-7			// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	(PP_EPS_ZERO*10)// Precision for MakeHyperplaneList()
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10)// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO		// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7			// Length of Objective Vector
#define PP_PROBE_LENGTH				0.001			// Length of probe shift
#define PP_REAL_TIME				3000			// This parameter limits the calculation time (compilator limit: 2 147 483 647)
//--------------------------------------------------------------------------

/*============================== recipe LP problem =========================*
// Number of equations: 79
// Subspace dimension: 101
#define PP_PROBLEM_NAME		"recipe"
#define PP_M 91		// Number of constraints in mps-file
#define PP_N 180	// Number of variables in mps-file
#define PP_MAX_OBJ_VALUE 266.616 // Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9			// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	(PP_EPS_ZERO*10)// Precision for MakeHyperplaneList()
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10)// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO		// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7			// Length of Objective Vector
#define PP_PROBE_LENGTH				0.001			// Length of probe shift
#define PP_REAL_TIME				1				// This parameter limits the calculation time (compilator limit: 2 147 483 647)
//--------------------------------------------------------------------------

/*============================== sc50a LP problem ==========================*
// Number of equations: 20
// Subspace dimension: 28
#define PP_PROBLEM_NAME		"sc50a"
#define PP_M 49	// Number of constraints
#define PP_N 48	// Number of variables
#define PP_MAX_OBJ_VALUE 64.575077058564509026860413914575	// Exact maximum value of objective function
//-------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-7			// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	(PP_EPS_ZERO*10)// Precision for MakeHyperplaneList()
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10)// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO		// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7			// Length of Objective Vector
#define PP_PROBE_LENGTH				0.001			// Length of probe shift
#define PP_REAL_TIME				1000			// This parameter limits the calculation time (compilator limit: 2 147 483 647)
//----------------------------------------------------------------------------
// Elapsed time: 17601
// Number of iterations: 3
// Computed objective value: 64.57507375411537
// Maximal objective value:  64.5750770585645
// Relative error = 5.12e-08
//----------------------------------------------------------------------------

/*============================== sc50b LP problem ============================*
// Number of equations: 20
// Subspace dimension: 28
#define PP_PROBLEM_NAME		"sc50b"
#define PP_M 48	// Number of constraints
#define PP_N 48	// Number of variables
#define PP_MAX_OBJ_VALUE 70	// Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-7			// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	(PP_EPS_ZERO*10)// Precision for MakeHyperplaneList()
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10)// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO		// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7			// Length of Objective Vector
#define PP_PROBE_LENGTH				0.001			// Length of probe shift
#define PP_REAL_TIME				3000			// This parameter limits the calculation time (compilator limit: 2 147 483 647)
//--------------------------------------------------------------------------
// Elapsed time: 2468.7376
// Number of iterations: 16
// Computed objective value: 68.62744550653315
// Maximal objective value:  70
// Relative error = 0.0196
//--------------------------------------------------------------------------

/*============================== sc105 LP problem ==========================*
#define PP_PROBLEM_NAME		"sc105"
#define PP_M 104	// Number of constraints in *.mps
#define PP_N 103	// Number of variables in *.mps
#define PP_MAX_OBJ_VALUE 52.202061211707248062628010857689 // Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_REAL_TIME 1000 // This parameter limits the calculation time (compilator limit: 2 147 483 647)
//----------------------------------------------------------------------------

/*============================== share2b LP problem ==============================*
#define PP_PROBLEM_NAME		"share2b"	
#define PP_M 96		// Number of constraints in *.mps
#define PP_N 162	// Number of variables in *.mps
#define PP_MAX_OBJ_VALUE 415.732240741419486545199108738 // Exact maximum value of objective function
//--------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9				// Precision for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	(PP_EPS_ZERO*10)	// Precision for MakeHyperplaneList()
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10)	// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO			// Precision of rounding vector r
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7				// Length of Objective Vector
#define PP_PROBE_LENGTH				0.0001				// Length of probe shift
#define PP_REAL_TIME				1000				// This parameter limits the calculation time (compilator limit: 2 147 483 647)
//----------------------------------------------------------------------------/**/