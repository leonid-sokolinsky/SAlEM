/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: SAlEM - Stochastic Along Edges Movement method
Module: _Problems15-1.h (LP problems of dimension 15 with 1 random inequality)
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
// PP_EPS_BIPPROJECTION_ROUND - inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH. 
//						This parameter affects terminate condition when 
//						calculating pseudoprojection.
//-------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-9				// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	PP_EPS_ZERO			// Precision for point to be in halfspace
//#define PP_EPS_MOVING_ON_POLYTOPE	(PP_EPS_ZERO/100)	// Precision for moving on polytope (affects Shift = 0)
#define PP_EPS_BIPPROJECTION_ROUND		PP_EPS_ZERO			// Precision of rounding pseudoprojecting vectors
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+10				// Length of Objective Vector
//-------------------------------------------------------------------------------
#define PP_M	16		// Number of equations (number of rows in *.mtx)
#define PP_N	31		// Number of variables (number of cols in *.mtx)
//=============================================================================

/*============================== rnd15-0 LP problem ===========================*
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd15-0"
#define PP_MAX_OBJ_VALUE 23900
//-----------------------------------------------------------------------------

/*============================== rnd15-1-1 LP problem =========================*
#define PP_PROBLEM_NAME	"rnd15-1-1"
#define PP_MAX_OBJ_VALUE 21054.89779933161
//-----------------------------------------------------------------------------

/*============================== rnd15-1-2 LP problem =========================*
#define PP_PROBLEM_NAME	"rnd15-1-2"
#define PP_MAX_OBJ_VALUE 21052.43702378703		
//-----------------------------------------------------------------------------

/*============================== rnd15-1-3 LP problem =========================*/
#define PP_PROBLEM_NAME	"rnd15-1-3"
#define PP_MAX_OBJ_VALUE 22224.70857556401		
//-----------------------------------------------------------------------------

/*============================== rnd15-1-4 LP problem =========================*
#define PP_PROBLEM_NAME	"rnd15-1-4"
#define PP_MAX_OBJ_VALUE 19041.38418137156		
//-----------------------------------------------------------------------------

/*============================== rnd15-1-5 LP problem =========================*
#define PP_PROBLEM_NAME	"rnd15-1-5"
#define PP_MAX_OBJ_VALUE 19048.83466963216		
//-----------------------------------------------------------------------------

/*============================== rnd15-1-6 LP problem =========================*
#define PP_PROBLEM_NAME	"rnd15-1-6"
#define PP_MAX_OBJ_VALUE 21538.96480089556		
//-----------------------------------------------------------------------------

/*=============================================================================*/