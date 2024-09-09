/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: SAlEM - Stochastic Along Edges Movement method (No MPI)
Module: _Problems-Klee-Minty-Cube.h (Problems from the LP-Klee-Minty-Cube Set)
Prefix: PP
Authors: Alexander E. Zhulev & Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
==============================================================================*/
#pragma once

//============================== Problem Parameters ======================
// PP_OBJECTIVE_VECTOR_LENGTH - direct dependence on dimension PD_n.
// P_EPS_ZERO - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH.
// PP_EPS_PROJECTION_ROUND - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH. 
//						This parameter affects terminate condition when 
//						calculating pseudoprojection.
//-----------------------------------------------------------------------
#define PP_PROBE_LENGTH			0.1		// Length of probe shift
//=============================================================================

/*============================== Klee-Minty5 LP problem ==============================*
// Starting point:	0 ... 0
// Exact solution:	0 ... 0	3125
#define PP_PROBLEM_NAME	"Klee-Minty5"
#define PP_D 5			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 3125
//---------------------------------- Method parameters ------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO		// Precision for MakeHyperplaneList()
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO		// Precision of rounding pseudoprojecting vectors
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+9			// Length of Objective Vector
//-----------------------------------------------------------------------------
// _________________________________________________ 1 _____________________________________________________
// Edge code: 4.   Hyperplanes: {5, 6, 7, 8}.      Shift = 4.9999952       F(x) = 80
// New vertex:                    5               0               0               0               0
// _________________________________________________ 2 _____________________________________________________
// Edge code: 4.   Hyperplanes: {0, 5, 6, 7}.      Shift = 5.0000012       F(x) = 120
// New vertex:                    5               5               0               0               0
// _________________________________________________ 3 _____________________________________________________
// Edge code: 4.   Hyperplanes: {0, 1, 5, 6}.      Shift = 65.000004       F(x) = 380
// New vertex:                    5               5              65               0               0
// _________________________________________________ 4 _____________________________________________________
// Edge code: 4.   Hyperplanes: {0, 1, 2, 5}.      Shift = 245.00002       F(x) = 870
// New vertex:                    5               5              65             245               0
// _________________________________________________ 5 _____________________________________________________
// Edge code: 4.   Hyperplanes: {0, 1, 2, 3}.      Shift = 1385    F(x) = 2255
// New vertex:                    5               5              65             245            1385
// _________________________________________________ 6 _____________________________________________________
// Edge code: 3.   Hyperplanes: {0, 1, 2, 4}.      Shift = 1010.1609       F(x) = 2745
// New vertex:                    5               5              65               0            2365
// _________________________________________________ 7 _____________________________________________________
// Edge code: 2.   Hyperplanes: {0, 1, 4, 6}.      Shift = 524.04676       F(x) = 3005
// New vertex:                    5               5               0               0            2885
// _________________________________________________ 8 _____________________________________________________
// Edge code: 1.   Hyperplanes: {0, 4, 6, 7}.      Shift = 80.156112       F(x) = 3045
// New vertex:                    5               0               0               0            2965
// _________________________________________________ 9 _____________________________________________________
// Edge code: 0.   Hyperplanes: {4, 6, 7, 8}.      Shift = 160.07813       F(x) = 3125
// New vertex:                    0               0               0               0            3125
//=============================================================================

/*============================== Klee-Minty6 LP problem ==============================*/
// Starting point:	0 ... 0
// Exact solution:	0 ... 0	15625
#define PP_PROBLEM_NAME	"Klee-Minty6"
#define PP_D 6			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 15625
//---------------------------------- Method parameters ------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	1E-8		// Precision for MakeHyperplaneList()
#define PP_EPS_PROJECTION_ROUND		1E-6			// Precision of rounding pseudoprojecting vectors
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+13			// Length of Objective Vector
//-----------------------------------------------------------------------------
// _________________________________________________ 1 _____________________________________________________
// Edge code: 5.   Hyperplanes: {0, 6, 7, 8, 9}.   Shift = 5       F(x) = 160
// New vertex:                    5               0               0               0               0               0
// _________________________________________________ 2 _____________________________________________________
// Edge code: 5.   Hyperplanes: {0, 1, 6, 7, 8}.   Shift = 5       F(x) = 240
// New vertex:                    5               5               0               0               0               0
// _________________________________________________ 3 _____________________________________________________
// Edge code: 5.   Hyperplanes: {0, 1, 2, 6, 7}.   Shift = 65      F(x) = 760
// New vertex:                    5               5              65               0               0               0
// _________________________________________________ 4 _____________________________________________________
// Edge code: 5.   Hyperplanes: {0, 1, 2, 3, 6}.   Shift = 245     F(x) = 1740
// New vertex:                    5               5              65             245               0               0
// _________________________________________________ 5 _____________________________________________________
// Edge code: 5.   Hyperplanes: {0, 1, 2, 3, 4}.   Shift = 1385    F(x) = 4510
// New vertex:                    5               5              65             245            1385               0
// _________________________________________________ 6 _____________________________________________________
// Edge code: 5.   Hyperplanes: {0, 1, 2, 3, 4}.   Shift = 6605    F(x) = 11115
// New vertex:                    5               5              65             245            1385            6605
// _________________________________________________ 7 _____________________________________________________
// Edge code: 4.   Hyperplanes: {0, 1, 2, 3, 7}.   Shift = 5710.5013       F(x) = 13885
// New vertex:                    5               5              65             245               0           12145
// _________________________________________________ 8 _____________________________________________________
// Edge code: 3.   Hyperplanes: {0, 1, 2, 7, 8}.   Shift = 1975.2531       F(x) = 14865
// New vertex:                    5               5              65               0               0           14105
// _________________________________________________ 9 _____________________________________________________
// Edge code: 2.   Hyperplanes: {0, 1, 7, 8, 9}.   Shift = 1042.0293       F(x) = 15385
// New vertex:                    5               5               0               0               0           15145
// _________________________________________________ 10 _____________________________________________________
// Edge code: 1.   Hyperplanes: {0, 7, 8, 9, 10}.  Shift = 160.07811       F(x) = 15465
// New vertex:                    5               0               0               0               0           15305
// _________________________________________________ 11 _____________________________________________________
// Edge code: 0.   Hyperplanes: {7, 8, 9, 10, 11}. Shift = 320.03906       F(x) = 15625
// New vertex:                    0               0               0               0               0           15625
//=============================================================================

/*============================== Klee-Minty7 LP problem =======================*
// Starting point:	0 ... 0
// Exact solution:	0 ... 0	78125
#define PP_PROBLEM_NAME	"Klee-Minty7"
#define PP_D 7			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 78125
//---------------------------------- Method parameters ------------------------
#define PP_EPS_ZERO					1E-9		// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	1E-6		// Precision for MakeHyperplaneList()
#define PP_EPS_PROJECTION_ROUND		1E-5		// Precision of rounding pseudoprojecting vectors
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+15		// Length of Objective Vector
// _________________________________________________ 1 _____________________________________________________
// Edge code: 6.   Hyperplanes: {0, 7, 8, 9, 10, 11}.      Shift = 5       F(x) = 320
// New vertex:                    5               0               0               0               0               0               0
// _________________________________________________ 2 _____________________________________________________
// Edge code: 6.   Hyperplanes: {0, 1, 7, 8, 9, 10}.       Shift = 5       F(x) = 480
// New vertex:                    5               5               0               0               0               0               0
// _________________________________________________ 3 _____________________________________________________
// Edge code: 6.   Hyperplanes: {0, 1, 2, 7, 8, 9}.        Shift = 65      F(x) = 1520
// New vertex:                    5               5              65               0               0               0               0
// _________________________________________________ 4 _____________________________________________________
// Edge code: 6.   Hyperplanes: {0, 1, 2, 3, 7, 8}.        Shift = 245     F(x) = 3480
// New vertex:                    5               5              65             245               0               0               0
// _________________________________________________ 5 _____________________________________________________
// Edge code: 6.   Hyperplanes: {0, 1, 2, 3, 4, 7}.        Shift = 1385    F(x) = 9020
// New vertex:                    5               5              65             245            1385               0               0
// _________________________________________________ 6 _____________________________________________________
// Edge code: 6.   Hyperplanes: {0, 1, 2, 3, 4, 5}.        Shift = 6605    F(x) = 22230
// New vertex:                    5               5              65             245            1385            6605               0
// _________________________________________________ 7 _____________________________________________________
// Edge code: 6.   Hyperplanes: {0, 1, 2, 3, 4, 5}.        Shift = 33665   F(x) = 55895
// New vertex:                    5               5              65             245            1385            6605           33665
// _________________________________________________ 8 _____________________________________________________
// Edge code: 5.   Hyperplanes: {0, 1, 2, 3, 4, 8}.        Shift = 27233.113       F(x) = 69105
// New vertex:                    5               5              65             245            1385               0           60085
// _________________________________________________ 9 _____________________________________________________
// Edge code: 4.   Hyperplanes: {0, 1, 2, 3, 8, 9}.        Shift = 11166.227       F(x) = 74645
// New vertex:                    5               5              65             245               0               0           71165
// _________________________________________________ 10 _____________________________________________________
// Edge code: 3.   Hyperplanes: {0, 1, 2, 8, 9, 10}.       Shift = 3927.6488       F(x) = 76605
// New vertex:                    5               5              65               0               0               0           75085
// _________________________________________________ 11 _____________________________________________________
// Edge code: 2.   Hyperplanes: {0, 1, 8, 9, 10, 11}.      Shift = 2081.0154       F(x) = 77645
// New vertex:                    5               5               0               0               0               0           77165
// _________________________________________________ 12 _____________________________________________________
// Edge code: 1.   Hyperplanes: {0, 8, 9, 10, 11, 12}.     Shift = 320.03906       F(x) = 77805
// New vertex:                    5               0               0               0               0               0           77485
// _________________________________________________ 13 _____________________________________________________
// Edge code: 0.   Hyperplanes: {8, 9, 10, 11, 12, 13}.    Shift = 640.01953       F(x) = 78125
// New vertex:                    0               0               0               0               0               0           78125
//=============================================================================

/*============================== Klee-Minty8 LP problem =======================*
// Starting point:	0 ... 0
// Exact solution:	0 ... 0	390625
// Computed objective value: 389984.9953464441
// Relative error = 0.00164
#define PP_PROBLEM_NAME	"Klee-Minty8"
#define PP_D 8			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 390625
//---------------------------------- Method parameters ------------------------
#define PP_EPS_ZERO					1E-9			// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	1E-4		// Precision for MakeHyperplaneList()
#define PP_EPS_PROJECTION_ROUND		1E-3			// Precision of rounding pseudoprojecting vectors
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+17			// Length of Objective Vector
//=============================================================================

/*============================== Klee-Minty9 LP problem =======================*
// Not solved!
// Start point:	   0 ... 0
// Exact solution: 0 ... 0	1953125
#define PP_PROBLEM_NAME	"Klee-Minty9"
#define PP_D 9			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 1953125
//---------------------------------- Method parameters ------------------------
#define PP_EPS_ZERO					1E-9	// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	1E-4	// Precision for MakeHyperplaneList()
#define PP_EPS_PROJECTION_ROUND		1E-3	// Precision of rounding pseudoprojecting vectors
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+18	// Length of Objective Vector
//=============================================================================

/*=============================================================================*/