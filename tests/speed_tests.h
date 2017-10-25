/*===========================================================================
 * Name         : speed_tests.h includes function to test
 *                the time-performance of basic operations.
 * Version      : 1.0.0, 7 Oct 2017
 *
 * Copyright (c) 2017 Ioannis Karagiannis
 * All rights reserved

 * This file is part of the LinearAlgebra library.

 * LinearAlgebra is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.

 * You are free to use this library under the terms of the GNU General
 * Public License, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with LinearAlgebra.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact info: https://www.linkedin.com/in/ioannis-karagiannis-7129394a/
 * 				ioanniskaragiannis1987@gmail.com
==============================================================================*/

#ifndef SPEED_TESTS_H_
#define SPEED_TESTS_H_

#include <chrono>  // C++11 feature
#include "../include/mat.h"

using namespace std::chrono;


// Below the addition performance for different-size square matrices is presented.
//    ____________________________________
//    |         |       ADDITION         |
//    | SIZE    |        (time)          |
//    |-----------------------------------
//    |   4x4   |       0.002 [ms]       |
//    ------------------------------------
//    |  10x10  |       0.002 [ms]       |
//    ------------------------------------
//    |  50x50  |       0.027 [ms]       |
//    ------------------------------------
//    | 100x100 |       0.141 [ms]       |
//    |-----------------------------------
//    | 250*250 |       1.203 [ms]       |
//    |-----------------------------------
//    | 500*500 |       4.644 [ms]       |
//    ------------------------------------
//    |1000*1000|      16.786 [ms]       |
//    ------------------------------------
//    |4000*4000|     244.062 [ms]        |
//    ------------------------------------
//
// MACHINE DETAILS: Intel® Core™ i3 CPU M 330 @ 2.13GHz × 4 (64-bit), 8GB RAM.


inline void test_matrix_addition_performance ( size_t size ){

	algebra::mat a = algebra::rand(size, size);
	algebra::mat b = algebra::rand(size, size);
	algebra::mat c;
	size_t iteration_nbr = 1;

	if(size > 1000){
		iteration_nbr = 2;
	}else{
		iteration_nbr = 8;
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for(size_t i = 0; i < iteration_nbr; i++){
		c = a + b;
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>( t2 - t1 ).count();
	double duration_f = (double) duration; duration_f = duration_f/iteration_nbr;
	printf("addition of two squared (%zux%zu) matrices lasted %.3f [ms]\n", a.rows(), a.cols(), duration_f/1000.0);
}

// Below the transposition performance for different-size square matrices is presented.
//    ____________________________________
//    |         |     TRANSPOSITION      |
//    | SIZE    |        (time)          |
//    |-----------------------------------
//    |   4x4   |       0.001 [ms]       |
//    ------------------------------------
//    |  10x10  |       0.003 [ms]       |
//    ------------------------------------
//    |  50x50  |       0.017 [ms]       |
//    ------------------------------------
//    | 100x100 |       0.053 [ms]       |
//    |-----------------------------------
//    | 250*250 |       0.533 [ms]       |
//    |-----------------------------------
//    | 500*500 |       3.238 [ms]       |
//    ------------------------------------
//    |1000*1000|      19.628 [ms]       |
//    ------------------------------------
//    |4000*4000|     384.778 [ms]       |
//    ------------------------------------
//
// MACHINE DETAILS: Intel® Core™ i3 CPU M 330 @ 2.13GHz × 4 (64-bit), 8GB RAM.
inline void test_matrix_transpose_performance ( size_t size ){

	algebra::mat a = algebra::rand(size, size);
	algebra::mat trans;
	size_t iteration_nbr = 1;

	if(size > 1000){
		iteration_nbr = 2;
	}else{
		iteration_nbr = 8;
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for(size_t i = 0; i < iteration_nbr; i++){
		trans = algebra::transpose(a);
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>( t2 - t1 ).count();
	double duration_f = (double) duration; duration_f = duration_f/iteration_nbr;
	printf("transposition of a squared (%zux%zu) matrix lasted %.3f [ms]\n", trans.rows(), trans.cols(), duration_f/1000.0);

}

// IT ASSUMES THAT MATRICES WILL BE SQUARE. THAT'S WHY ONLY SIZE IS DETERMINED.
inline void test_matrix_normal_multiplication_performance ( size_t size ){

	algebra::mat a = algebra::rand(size, size);
	algebra::mat b = algebra::rand(size, size);
	algebra::mat c;
	size_t iteration_nbr = 1;

	if(size > 1000){
		iteration_nbr = 2;
	}else{
		iteration_nbr = 8;
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for(size_t i = 0; i < iteration_nbr; i++){
		c = a * b;
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>( t2 - t1 ).count();
	double duration_f = (double) duration; duration_f = duration_f/iteration_nbr;
	printf("normal-multiplication of two squared (%zux%zu) matrices lasted %.3f [ms]\n", a.rows(), a.cols(), duration_f/1000.0);

}

//  The table below presents the performance of the two multiplication methods
//  for several sizes of square matrices.
//    ______________________________________________________________
//    |         |  NORMAL MULTIPLICATION | STRASSEN MULTIPLICATION |
//    | SIZE    |        (time)          |        (time)           |
//    |-------------------------------------------------------------
//    |   4x4   |        0.002 [ms]      |        0.236 [ms]       |
//    --------------------------------------------------------------
//    |  10x10  |        0.007 [ms]      |        6.070 [ms]       |
//    |-------------------------------------------------------------
//    |  50x50  |        0.298 [ms]      |       65.528 [ms]       |
//    --------------------------------------------------------------
//    | 100x100 |        2.155 [ms]      |      112.273 [ms]       |
//    --------------------------------------------------------------
//    | 250x250 |       26.642 [ms]      |     241.910 [ms]        |
//    --------------------------------------------------------------
//    | 500x500 |      201.615 [ms]      |     733.256 [ms]        |
//    |-------------------------------------------------------------
//    |1000*1000|        1.891 [s]       |        2.921 [s]        |
//    |-------------------------------------------------------------
//    |4000*4000|      115.820 [s]       |       91.807 [s]        |
//    --------------------------------------------------------------

// MACHINE DETAILS: Intel® Core™ i3 CPU M 330 @ 2.13GHz × 4 (64-bit), 8GB RAM.
inline void test_matrix_strassen_multiplication_performance ( size_t size ){

	algebra::mat a = algebra::rand(size, size);
	algebra::mat b = algebra::rand(size, size);
	algebra::mat c;
	size_t iteration_nbr = 1;

	if(size > 1000){
		iteration_nbr = 2;
	}else{
		iteration_nbr = 8;
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for(size_t i = 0; i < iteration_nbr; i++){
		c = algebra::strassen(a, b);
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>( t2 - t1 ).count();
	double duration_f = (double) duration; duration_f = duration_f/iteration_nbr;
	printf("strassen-multiplication of two squared (%zux%zu) matrices lasted %.3f [ms]\n", a.rows(), a.cols(), duration_f/1000.0);

}

// Below the inversion performance for different-size square matrices is presented.
//    ____________________________________
//    |         |       INVERSION        |
//    | SIZE    |        (time)          |
//    |-----------------------------------
//    |   4x4   |       0.002 [ms]       |
//    ------------------------------------
//    |  10x10  |       0.011 [ms]       |
//    ------------------------------------
//    |  50x50  |       0.425 [ms]       |
//    ------------------------------------
//    | 100x100 |       2.999 [ms]       |
//    |-----------------------------------
//    | 250*250 |      44.178 [ms]       |
//    |-----------------------------------
//    | 500*500 |     471.896 [ms]       |
//    ------------------------------------
//    |1000*1000|       9.898 [s]        |
//    ------------------------------------
//    |4000*4000|      14.443 [min]       |
//    ------------------------------------
//
// MACHINE DETAILS: Intel® Core™ i3 CPU M 330 @ 2.13GHz × 4 (64-bit), 8GB RAM.

inline void test_matrix_inversion_performance ( size_t size ){

	algebra::mat a = algebra::rand(size, size);
	algebra::mat b = a;
	size_t iteration_nbr = 1;

	if(size > 1000){
		iteration_nbr = 2;
	}else{
		iteration_nbr = 8;
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for(size_t i = 0; i < iteration_nbr; i++){
		b = algebra::inv(a);
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>( t2 - t1 ).count();
	double duration_f = (double) duration; duration_f = duration_f/iteration_nbr;
	printf("inversion of a squared (%zux%zu) matrix lasted %.3f [ms]\n", a.rows(), a.cols(), duration_f/1000.0);
}

// Uncomment them all if you want to test all the operations at once.
inline void test_speed_of_basic_operations(size_t size){
	test_matrix_addition_performance ( size );
	test_matrix_transpose_performance ( size );
	test_matrix_normal_multiplication_performance ( size );
	test_matrix_strassen_multiplication_performance ( size );
	test_matrix_inversion_performance( size );
}

#endif /* SPEED_TESTS_H_ */
