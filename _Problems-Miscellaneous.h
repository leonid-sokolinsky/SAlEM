/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: SAlEM - Stochastic Along Edges Movement method
Module: _Problems-Miscellaneous.h (Miscellaneous LP problems)
Prefix: PP
Author: Leonid B. Sokolinsky
This include file is part of Problem-Parameters.h
Start vertex *_v.mtx for these problems was calculated by VeSP https://github.com/leonid-sokolinsky/VeSP
LP problems are available in https://github.com/leonid-sokolinsky/Set-of-LP-Problems/tree/main/Miscellaneous-LP
================================================================================*/
#pragma once

//=========================== problem Parameters ===============================
// PP_OBJECTIVE_VECTOR_LENGTH - direct dependence on dimension PD_n.
// PP_EPS_PROJECTION -	inverse dependence on PP_OBJECTIVE_VECTOR_LENGTH. 
//						This parameter affects terminate condition when 
//						calculating pseudoprojection.
#define PP_BATCH_SIZE					10
#define PP_EPS_RELATIVE_ERROR			1E-5			// Used if defined PP_CHECK_MAX_OBJ_VALUE 
#define PP_SCALE_FACTOR					1				// #ifdef PP_SAVE_LOCAL_RESULT; makes 9 digits before the decimal point of PP_MAX_OBJ_VALUE
#define PP_MAX_PSEUDOPROJECTING_ITER	100000000		// Maximum acceptable number of iterations in Pseudoprojection on flat
//------------------------------ ifdef PP_DEBUG --------------------------------
#define PP_PROJECTION_COUNT				100000				// Each PP_PROJECTION_COUNT iteration to be outputted inside Flat_MaxProjection(*)
//==============================================================================

/*============================== nguyen5 LP problem ============================*
#define PP_PROBLEM_NAME	"nguyen5"
#define PP_MPS_FORMAT
#define PP_M 4		// Number of constrains
#define PP_N 5		// Number of variables
#define PP_MAX_OBJ_VALUE 21.4549732313098893143888
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Precision for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 4
// Computed objective value: 21.4549732313098786562477
// Maximal objective value:  21.4549732313097933911195
// Relative error = 3.97e-15
// Distance to polytope: 1.7781997e-15
//------------------------------------------------------------------------------

/*============================== wiki LP problem ============================*
#define PP_PROBLEM_NAME	"wiki"
#define PP_MPS_FORMAT
#define PP_M 2		// Number of constrains
#define PP_N 3		// Number of variables
#define PP_MAX_OBJ_VALUE 20
// https://en.wikipedia.org/wiki/Simplex_algorithm
//------------------------------------------------------------------------------
// https://sagecell.sagemath.org/
// p = Polyhedron(ieqs = [[10,-3,-2,-1],[15,-2,-5,-3],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
// p.plot()
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Precision for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 2
// Computed objective value: 20.0000000000000213162821
// Maximal objective value:  20
// Relative error = 1.07e-15
// Distance to polytope: 7.7804045e-15
//------------------------------------------------------------------------------

/*============================== simpleCube LP problem =========================*
#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"simpleCube"
#ifdef PP_MPS_FORMAT
#define PP_M 3		// Number of constrains
#define PP_N 3		// Number of variables
#else
#define PP_M 3	// Number of rows in *.mtx
#define PP_N 6	// Number of cols in *.mtx
#endif
#define PP_MAX_OBJ_VALUE 		60000
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Precision for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7	// Length of Objective Vector
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 3
// Computed objective value: 60000
// Maximal objective value:  60000
// Relative error = 0
// Distance to polytope: 0
//------------------------------------------------------------------------------

/*============================== square3D LP problem =========================*
#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"square3D"
#define PP_M 4		// Number of constrains
#define PP_N 3		// Number of variables
#define PP_MAX_OBJ_VALUE 		50000
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Precision for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 3
// Computed objective value: 50000
// Maximal objective value:  50000
// Relative error = 0
//------------------------------------------------------------------------------

/*============================== square4D LP problem =========================*
#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"square4D"
#define PP_M 7		// Number of constrains
#define PP_N 4		// Number of variables
#define PP_MAX_OBJ_VALUE 		36200
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Precision for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_BASIC_VECTORS_ONLY
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 2
// Computed objective value: 36200
// Maximal objective value:  36200
// Relative error = 0
// Distance to polytope: 0
//------------------------------------------------------------------------------

/*============================== pyramid LP problem =========================*
#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"pyramid"
#define PP_M 3		// Number of constrains
#define PP_N 3		// Number of variables
#define PP_MAX_OBJ_VALUE 		60000
//------------------------------------------------------------------------------
// https://sagecell.sagemath.org/
// p = Polyhedron(ieqs = [[0,1,0,0],[0,0,1,0,],[0,-1,0,1],[0,0,-1,1],[200,0,0,-1,]])
// p.plot()
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6					// Length of Objective Vector
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 3
// Computed objective value: 50000
// Maximal objective value:  50000
// Relative error = 0
//------------------------------------------------------------------------------

/*============================== featheredCube LP problem ======================*
#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"featheredCube"
#ifdef PP_MPS_FORMAT
#define PP_M 15		// Number of constrains
#define PP_N 3		// Number of variables
#else
#define PP_M 3	// Number of rows in *.mtx
#define PP_N 6	// Number of cols in *.mtx
#endif
#define PP_MAX_OBJ_VALUE 		60000
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 3
// Computed objective value: 60000
// Maximal objective value:  60000
// Relative error = 0
// Distance to polytope: 0
//------------------------------------------------------------------------------

/*============================== cubeAndHyperplane LP problem ===================*
#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"cubeAndHyperplane"
#define PP_M 4		// Number of constrains
#define PP_N 4		// Number of variables
#define PP_MAX_OBJ_VALUE 		90000
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero
#define PP_EPS_PROJECTION			1E-10		// Accuracy of belonging to hyperplane
#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 2
// Computed objective value: 89999.9999999999854480848
// Maximal objective value:  90000
// Relative error = 1.62e-16
// Distance to polytope: 5.6843419e-14
//------------------------------------------------------------------------------

/*============================== simple1FxVar LP problem =======================*
// Simple LP problem & x_1=150
#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"simple1FxVar"
#define PP_M 4		// Number of constraints
#define PP_N 3		// Number of variables
#define PP_MAX_OBJ_VALUE 52500
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 1
// Computed objective value: 52500
// Maximal objective value:  52500
// Relative error = 0
// Distance to polytope: 0
//------------------------------------------------------------------------------

/*============================== simple1 LP problem ============================*
#define PP_PROBLEM_NAME	"simple1"
#define PP_MPS_FORMAT
#ifdef PP_MPS_FORMAT
#define PP_M 4		// Number of constrains
#define PP_N 3		// Number of variables
#else
#define PP_M 4		// Number of rows in *.mtx
#define PP_N 7		// Nnumber of cols in *.mtx
#endif
#define PP_MAX_OBJ_VALUE 		55000
//------------------------------------------------------------------------------
// https://sagecell.sagemath.org/
// p = Polyhedron(ieqs = [[200,-1,0,0],[200,0,-1,0],[200,0,0,-1],[500,-1,-1,-1],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
// p.plot()
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 3
// Computed objective value: 55000
// Maximal objective value:  55000
// Relative error = 0
// Distance to polytope: 0
//------------------------------------------------------------------------------

/*============================== simple1.1 LP problem ==========================*
// Simple LP problem with alternating objective function
#define PP_PROBLEM_NAME	"simple1.1"
#define PP_M 4		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 		40000
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 2
// Computed objective value: 40000
// Maximal objective value:  40000
// Relative error = 0
// Distance to polytope: 0
//------------------------------------------------------------------------------

/*============================== simple2 LP problem ============================*
// Simple LP problem & x_3=200; x_2>=110; x_0<=190
#define PP_PROBLEM_NAME	"simple2"
#define PP_MPS_FORMAT
#ifdef PP_MPS_FORMAT
#define PP_M 5		// Number of constrains
#define PP_N 4		// Number of variables
#else
#define PP_M 5		// Number of rows in *.mtx
#define PP_N 8		// Nnumber of cols in *.mtx
#endif
#define PP_MAX_OBJ_VALUE 		63500
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 1
// Computed objective value: 63500
// Maximal objective value:  63500
// Relative error = 0
// Distance to polytope: 0
//------------------------------------------------------------------------------

/*============================== simple2' LP problem ===========================*
// Simple LP problem & x_3=200; 2*x_3=400; x_2>=110; x_0<=190
#define PP_MPS_FORMAT
#define PP_PROBLEM_NAME	"simple2'"
#define PP_M 6		// Number of constraints
#define PP_N 4		// Number of variables
#define PP_MAX_OBJ_VALUE 	63500
//------------------------------------------------------------------------------
#define PP_EPS_ZERO				1E-14					// Accuracy for comparison with zero
#define PP_EPS_ON_HYPERPLANE	1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7				// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
#define PP_BASIC_VECTORS_ONLY
//#define PP_MAXPROJECTION
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 2
// Computed objective value: 63500
// Distance to polytope: 0
//------------------------------------------------------------------------------

/*============================== simple3 LP problem ============================*
#define PP_PROBLEM_NAME	"simple3"
#define PP_MPS_FORMAT
#ifdef PP_MPS_FORMAT
#define PP_M 5		// Number of constrains
#define PP_N 5		// Number of variables
#else
#define PP_M 6		// Number of rows in *.mtx
#define PP_N 8		// Nnumber of cols in *.mtx
#endif
#define PP_MAX_OBJ_VALUE 		55000
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 1
// Computed objective value: 55000.0000000000072759576
// Maximal objective value:  55000
// Relative error = 1.32e-16
//------------------------------------------------------------------------------

/*============================== simple1min LP problem =========================*
#define PP_PROBLEM_NAME	"simple1min"
#define PP_M 5		// Number of equations (number of rows in *.mtx)
#define PP_N 8		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 		-5000
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 2
// Computed objective value: -4999.9999999999990905053
// Maximal objective value:  -5000
// Relative error = 1.82e-16
//------------------------------------------------------------------------------

/*============================== simple_zcv LP problem =========================*
#define PP_PROBLEM_NAME	"simple_zcv"
#define PP_M 4		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 		50000
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 3
// Computed objective value: 50000
// Maximal objective value:  50000
// Relative error = 0
//------------------------------------------------------------------------------

/*============================== simple_lcv LP problem =========================*
#define PP_PROBLEM_NAME	"simple_lcv"
#define PP_M 4		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 		50000.2
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 3
// Computed objective value: 50000.199999999997089617
// Maximal objective value:  50000.199999999997089617
// Relative error = 0
//------------------------------------------------------------------------------

/*============================== simple_lcv_neg LP problem =====================*
#define PP_PROBLEM_NAME	"simple_lcv_neg"
#define PP_M 4		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 		49998
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 2
// Computed objective value: 49998
// Maximal objective value:  49998
// Relative error = 0
// Distance to polytope: 0
//------------------------------------------------------------------------------

/*============================== angle03 LP problem ============================*
#define PP_PROBLEM_NAME	"angle03"
#define PP_M 3		// Number of equations (number of rows in *.mtx)
#define PP_N 6		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 		3000
//------------------------------------------------------------------------------
// https://sagecell.sagemath.org/
// p = Polyhedron(ieqs = [[55,-0.5,-5,0],[48,-4,-0.8,0],[10,0,0,-1],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
// p.plot()
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 1
// Computed objective value: 3000
// Maximal objective value:  3000
// Relative error = 0
//------------------------------------------------------------------------------

/*============================== angle04 LP problem ============================*
#define PP_PROBLEM_NAME	"angle04"
#define PP_M 3		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 		3300
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 1
// Computed objective value: 3300
// Maximal objective value:  3300
// Relative error = 0
//------------------------------------------------------------------------------

/*============================== cone3-0 LP problem ============================*
#define PP_PROBLEM_NAME	"cone3-0"
#define PP_M 11		// Number of equations (number of rows in *.mtx)
#define PP_N 14		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 115
//------------------------------------------------------------------------------
// https://sagecell.sagemath.org/
// p = Polyhedron(ieqs = [[2500,-50,-44,-35],[2500,-50,-28,-42],[2500,-50,-31,-37],[2500,-50,-30,-39],[2500,-50,-43,-26],[2500,-50,-29,-45],[2500,-50,-41,-44],[2500,-50,-26,-27],[50,-1,0,0],[25,0,-1,0],[25,0,0,-1],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
// p.plot()
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-6	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+7					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 3
// Computed objective value: 114.999999999994329868969
// Maximal objective value:  115
// Relative error = 4.93e-14
// Distance to polytope: 3.907985e-14
//------------------------------------------------------------------------------

/*============================== rnd3_3_Olkhovsky56 LP problem =================*
#define PP_PROBLEM_NAME	"rnd3_3_Olkhovsky56"
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 9		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 452.604395604395620011928
//------------------------------------------------------------------------------
// https://sagecell.sagemath.org/
// p = Polyhedron(ieqs = [[200,-1,0,0],[200,0,-1,0],[200,0,0,-1],[25987,-83,-23,-91],[42580,-32,-121,-138],[62862,-166,-64,-169],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
// p.plot()
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 3
// Computed objective value: 452.604395604395620011928
// Maximal objective value:  452.604395604395620011928
// Relative error = 0
// Distance to polytope: 0
//------------------------------------------------------------------------------

/*============================== rnd3_3_Olkhovsky289 LP problem ================*
#define PP_PROBLEM_NAME	"rnd3_3_Olkhovsky289"
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 9		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 431.817389789161893531855
//------------------------------------------------------------------------------
// https://sagecell.sagemath.org/
// p = Polyhedron(ieqs = [[200,-1,0,0],[200,0,-1,0],[200,0,0,-1],[63439,-193,-121,-156],[43414,-112,-59,-150],[38965,-187,-12,-96],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
// p.plot()
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Precision for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-8	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_MAXPROJECTION 
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 5
// Computed objective value: 431.817389789161893531855
// Maximal objective value:  431.817389789161893531855
// Relative error = 0
// Distance to polytope: 0
//------------------------------------------------------------------------------

/*============================== rnd3_3_Olkhovsky331 LP problem ================*
#define PP_PROBLEM_NAME	"rnd3_3_Olkhovsky331"
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 9		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 441.843642611683833365532
//------------------------------------------------------------------------------
// https://sagecell.sagemath.org/
// p = Polyhedron(ieqs = [[200,-1,0,0],[200,0,-1,0],[200,0,0,-1],[31312,-64,-90,-63],[8324,-25,-30,-4],[56572,-168,-18,-173],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
// p.plot()
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 5
// Computed objective value: 441.843642611665870845172
// Maximal objective value:  441.843642611683833365532
// Relative error = 4.07e-14
// Distance to polytope: 0
//------------------------------------------------------------------------------

/*============================== rnd3_3_Olkhovsky336 LP problem ================*
#define PP_PROBLEM_NAME	"rnd3_3_Olkhovsky336"
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 9		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 436.613602166716816554981
//------------------------------------------------------------------------------
// https://sagecell.sagemath.org/
// p = Polyhedron(ieqs = [[200,-1,0,0],[200,0,-1,0],[200,0,0,-1],[25928,-32,-67,-92],[47149,-181,-53,-105],[62600,-24,-198,-193],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
// p.plot()
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 4
// Computed objective value: 436.6136021667168733984
// Maximal objective value:  436.613602166716816554981
// Relative error = 1.3e-16
// Distance to polytope: 6.7414062e-14
//------------------------------------------------------------------------------

/*============================== rnd3-10 LP problem ============================*
#define PP_PROBLEM_NAME	"rnd3-10"
#define PP_M 13		// Number of equations (number of rows in *.mtx)
#define PP_N 16		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 397.938973646311296761269
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11					// Precision for comparison with zero

#define PP_EPS_ON_HYPERPLANE		1E-9	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6					// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 5
// Computed objective value: 397.938973646311410448106
// Maximal objective value:  397.938973646311296761269
// Relative error = 2.86e-16
// Distance to polytope: 3.7706941e-13
//------------------------------------------------------------------------------

/*============================== rnd5-100 LP problem ===========================*/
#define PP_PROBLEM_NAME	"rnd5-100"
#define PP_M 105		// Number of equations (number of rows in *.mtx)
#define PP_N 110		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 617.504337008263632924354
//------------------------------------------------------------------------------
#define PP_EPS_ZERO					1E-11	// Precision for comparison with zero
#define PP_EPS_ON_HYPERPLANE		1E-8	// Accuracy of belonging to hyperplane
#define PP_OBJECTIVE_VECTOR_LENGTH	1E+6	// Length of Objective Vector
//-------------------------- Compilation Modes ---------------------------------
//#define PP_GRADIENT
//------------------------------------------------------------------------------
// Elapsed time: 0
// Number of iterations: 11
// Computed objective value: 617.504337008263746611192
// Maximal objective value:  617.504337008263632924354
// Relative error = 1.84e-16
// Distance to polytope: 1.0168552e-13
//------------------------------------------------------------------------------

/*==============================================================================*/