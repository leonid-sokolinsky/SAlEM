/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: SAlEM - Stochastic Along Edges Movement method
Module: _Problems-Klee-Minty-Cube.h (Problems from the LP-Klee-Minty-Cube Set)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
Start vertex *_v.mtx for these problems was calculated by VeSP https://github.com/leonid-sokolinsky/VeSP
LP problems were obtained using Klee-Minty-Cube-Generator https://github.com/leonid-sokolinsky/Klee-Minty-Cube-Generator
LP problems are available in https://github.com/leonid-sokolinsky/Set-of-LP-Problems/tree/main/Klee-Minty-Cube
================================================================================*/
#pragma once

//#define PP_GRADIENT

//============================== Problem Parameters ============================
// PP_OBJECTIVE_VECTOR_LENGTH - direct dependence on dimension PD_n.
// P_EPS_ZERO - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH.
//------------------------------------------------------------------------------
#ifndef PP_GRADIENT
#define PP_EPS_ZERO					1E-11			// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10)// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION			PP_EPS_ZERO		// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7			// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
#define PP_MAXPROJECTION // It can cause a stuck in the loop. In this case, you should increase PP_EPS_PROJECTION or decrease PP_OBJECTIVE_VECTOR_LENGTH.
//------------------------------------------------------------------------------
#endif
#define PP_EPS_RELATIVE_ERROR		1E-3			// Used if defined PP_CHECK_MAX_OBJ_VALUE 
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT				100000000			// Each PP_PROJECTION_COUNT iteration to be outputted inside Flat_MaxProjection(*)

//------------------------------------------------------------------------------

/*============================== Klee-Minty5 LP problem ========================*
#define PP_PROBLEM_NAME	"Klee-Minty5"
#define PP_D 5			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 3125
//------------------------------------------------------------------------------
#ifdef PP_GRADIENT
#define PP_EPS_ZERO					1E-10				// Accuracy for comparison with zero
#define PP_EPS_PROJECTION			(PP_EPS_ZERO*10)	// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_PROJECTION*1000)	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7				// Length of Objective Vector
#endif
//------------------------------------------------------------------------------

/*============================== Klee-Minty6 LP problem ========================*
#define PP_PROBLEM_NAME	"Klee-Minty6"
#define PP_D 6			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 15625
//------------------------------------------------------------------------------
#ifdef PP_GRADIENT
#define PP_EPS_ZERO					1E-10				// Accuracy for comparison with zero
#define PP_EPS_PROJECTION			(PP_EPS_ZERO*10)	// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10000)	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7				// Length of Objective Vector
#endif
//------------------------------------------------------------------------------

/*============================== Klee-Minty7 LP problem ========================*
#define PP_PROBLEM_NAME	"Klee-Minty7"
#define PP_D 7			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 78125
//------------------------------------------------------------------------------
#ifdef PP_GRADIENT
#define PP_EPS_ZERO					1E-11				// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10000)	// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION			(PP_EPS_ZERO*10)	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7				// Length of Objective Vector
#endif
//------------------------------------------------------------------------------

/*============================== Klee-Minty8 LP problem ========================*
#define PP_PROBLEM_NAME	"Klee-Minty8"
#define PP_D 8			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 390625
//------------------------------------------------------------------------------
#ifdef PP_GRADIENT
#define PP_EPS_ZERO					1E-11				// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1e-06				// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION			(PP_EPS_ZERO*10)	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7				// Length of Objective Vector
#endif
//------------------------------------------------------------------------------

/*============================== Klee-Minty9 LP problem ========================*/
#define PP_PROBLEM_NAME	"Klee-Minty9"
#define PP_D 9			// Space dimension
#define PP_M PP_D		// Number of equations (number of rows in *.mtx)
#define PP_N (2*PP_D)	// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 1953125
//------------------------------------------------------------------------------
#ifdef PP_GRADIENT
#define PP_EPS_ZERO					1E-11				// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1e-06				// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION			(PP_EPS_ZERO*10)	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7				// Length of Objective Vector
#endif
//------------------------------------------------------------------------------

/*==============================================================================*/