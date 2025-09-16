/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: SAlEM - Stochastic Along Edges Movement method
Module: _Problems05-1.h (LP problems of dimension 5 with 1 random inequality: LP-Rnd-Problems Set)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
Start vertex *_v.mtx for these problems was calculated by VeSP https://github.com/leonid-sokolinsky/VeSP
LP problems were obtained using LPP-Generator https://github.com/leonid-sokolinsky/LPP-Generator
LP problems are available in https://github.com/leonid-sokolinsky/Set-of-LP-Problems/tree/main/Rnd-LP
================================================================================*/
#pragma once

//=========================== Problem Parameters ===============================
// PP_OBJECTIVE_VECTOR_LENGTH - direct dependence on dimension PD_n.
// P_EPS_ZERO - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH.
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
#define PP_M						6		// Number of equations (number of rows in *.mtx)
#define PP_N						11		// Number of variables (number of cols in *.mtx)
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT				100000				// Each PP_PROJECTION_COUNT iteration to be outputted inside Flat_MaxProjection(*)

//------------------------------------------------------------------------------
#define PP_EPS_RELATIVE_ERROR		1E-3			// Used if defined PP_CHECK_MAX_OBJ_VALUE 

/*============================== rnd5-0 LP problem =============================*
#define PP_PROBLEM_NAME	"rnd5-0"
#define PP_MAX_OBJ_VALUE 		2900
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11				// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*100)	// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION			PP_EPS_ZERO			// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7				// Length of Objective Vector
//------------------------------------------------------------------------------

/*============================== rnd5-1-1 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd5-1-1"
#define PP_MAX_OBJ_VALUE 2584.34948970919685962144
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11				// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*100)	// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION			PP_EPS_ZERO			// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7				// Length of Objective Vector
//------------------------------------------------------------------------------

/*============================== rnd5-1-2 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd5-1-2"
#define PP_MAX_OBJ_VALUE 2657.52561253995372680947
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11				// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*100)	// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION			(PP_EPS_ZERO*10)	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6				// Length of Objective Vector
//------------------------------------------------------------------------------

/*============================== rnd5-1-3 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd5-1-3"
#define PP_MAX_OBJ_VALUE 2424.91915381191074629896
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11				// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*100)	// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION			PP_EPS_ZERO			// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7				// Length of Objective Vector
//------------------------------------------------------------------------------

/*============================== rnd5-1-4 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd5-1-4"
#define PP_MAX_OBJ_VALUE 2300.11275869818382489029
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11				// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*100)	// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION			(PP_EPS_ZERO*10)	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6				// Length of Objective Vector
//------------------------------------------------------------------------------

/*============================== rnd5-1-5 LP problem ===========================*
#define PP_PROBLEM_NAME	"rnd5-1-5"
#define PP_MAX_OBJ_VALUE 2626.47323620733004645444
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11				// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*100)	// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION			(PP_EPS_ZERO*10)	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6				// Length of Objective Vector
//------------------------------------------------------------------------------

/*============================== rnd5-1-6 LP problem ===========================*/
#define PP_PROBLEM_NAME	"rnd5-1-6"
#define PP_MAX_OBJ_VALUE 2675.35199418665024495567
//------------------------------------------------------------------------------
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7				// Length of Objective Vector
#define PP_EPS_ON_HYPERPLANE		(PP_EPS_ZERO*10)	// Accuracy of belonging to hyperplane
#define PP_EPS_PROJECTION			(PP_EPS_ZERO*10)	// Accuracy of belonging to hyperplane
#define PP_EPS_ZERO					1E-11				// Accuracy for comparison with zero
//-------------------------- Compilation Modes ---------------------------------
#define PP_MAXPROJECTION // It can cause a stuck in the loop. In this case, you should increase PP_EPS_PROJECTION or decrease PP_OBJECTIVE_VECTOR_LENGTH.
//------------------------------------------------------------------------------

/*===============================================================================*/