/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: SAlEM - Stochastic Along Edges Movement method
Module: _Problems05-1.h (LP problems of dimension 5 with 1 random inequality: LP-Rnd-Problems Set)
Prefix: PP
Authors: Alexander E. Zhulev & Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
LP problems were obtained using BSF-LPP-Generator.
Initial surface points for these problems were calculated using BSF-Apex-Quest.
==============================================================================*/
#pragma once

//=========================== Problem Parameters ========================
// PP_OBJECTIVE_VECTOR_LENGTH - direct dependence on dimension PD_n.
// P_EPS_ZERO - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH.
// PP_EPS_PROJECTION_ROUND - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH. 
//						This parameter affects terminate condition when 
//						calculating pseudoprojection.
//-----------------------------------------------------------------------
#define PP_EPS_ZERO					1E-10		// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	1E-9		// Precision for point to be in halfspace
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO	// Precision of rounding pseudoprojecting vectors
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+9		// Length of Objective Vector
//-----------------------------------------------------------------------
#define PP_M						6		// Number of equations (number of rows in *.mtx)
#define PP_N						11		// Number of variables (number of cols in *.mtx)
#define PP_PROBE_LENGTH				0.1		// Length of probe shift
//-------------------------------------------------------------------------------

/*============================== rnd5-0 LP problem ==============================*
// Exact solution:	100   200   200   200   200
#define PP_PROBLEM_NAME	"rnd5-0"
#define PP_MAX_OBJ_VALUE 		2900
//-------------------------------------------------------------------------------

/*============================== rnd5-1-1 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd5-1-1"
#define PP_MAX_OBJ_VALUE 2584.349489709197
//-------------------------------------------------------------------------------

/*============================== rnd5-1-2 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd5-1-2"
#define PP_MAX_OBJ_VALUE 2657.525612539948
//-------------------------------------------------------------------------------

/*============================== rnd5-1-3 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd5-1-3"
#define PP_MAX_OBJ_VALUE 2424.919153811911
//-------------------------------------------------------------------------------

/*============================== rnd5-1-4 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd5-1-4"
#define PP_MAX_OBJ_VALUE 2300.112758698184
//-------------------------------------------------------------------------------

/*============================== rnd5-1-5 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd5-1-5"
#define PP_MAX_OBJ_VALUE 2626.473236207309
//-------------------------------------------------------------------------------

/*============================== rnd5-1-6 LP problem ============================*/
#define PP_PROBLEM_NAME	"rnd5-1-6"
#define PP_MAX_OBJ_VALUE 2675.35199418665
//-------------------------------------------------------------------------------

/*===============================================================================*/