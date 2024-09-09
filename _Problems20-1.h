/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: SAlEM - Stochastic Along Edges Movement method
Module: _Problems20-1.h (LP problems of dimension 20 with 1 random inequality)
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
#define PP_EPS_ZERO					1E-9				// Accuracy for comparison with zero
#define PP_EPS_POINT_IN_HALFSPACE	1E-8				// Precision for point to be in halfspace
//#define PP_EPS_MOVING_ON_POLYTOPE	(PP_EPS_ZERO/100)	// Precision for moving on polytope (affects Shift = 0)
#define PP_EPS_PROJECTION_ROUND		PP_EPS_ZERO			// Precision of rounding pseudoprojecting vectors
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+9				// Length of Objective Vector
//-------------------------------------------------------------------------------
#define PP_M	21		// Number of equations (number of rows in *.mtx)
#define PP_N	41		// Number of variables (number of cols in *.mtx)

#define PP_PROBE_LENGTH			1					// Length of probe shift
#define PP_START_SHIFT_LENGTH	100
//=============================================================================

/*============================== rnd20-0 LP problem =========================*/
// Solution:	100  200  ...  200
#define PP_PROBLEM_NAME	"rnd20-0"
#define PP_MAX_OBJ_VALUE 41900
//-----------------------------------------------------------------------------

/*============================== rnd20-1-1 LP problem =========================*
#define PP_PROBLEM_NAME	"rnd20-1-1"
#define PP_MAX_OBJ_VALUE 29920.01741110406
//-----------------------------------------------------------------------------

/*============================== rnd20-1-2 LP problem =========================*
#define PP_PROBLEM_NAME	"rnd20-1-2"
#define PP_MAX_OBJ_VALUE 37836.93296989189		
//-----------------------------------------------------------------------------

/*============================== rnd20-1-3 LP problem =========================*
#define PP_PROBLEM_NAME	"rnd20-1-3"
#define PP_MAX_OBJ_VALUE 40150.02344132274		
//-----------------------------------------------------------------------------

/*============================== rnd20-1-4 LP problem =========================*
#define PP_PROBLEM_NAME	"rnd20-1-4"
#define PP_MAX_OBJ_VALUE 35207.9205455114		
//-----------------------------------------------------------------------------

/*============================== rnd20-1-5 LP problem =========================*
#define PP_PROBLEM_NAME	"rnd20-1-5"
#define PP_MAX_OBJ_VALUE 33275.82728933131		
//-----------------------------------------------------------------------------

/*=============================================================================*/