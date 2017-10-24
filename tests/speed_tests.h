/*===========================================================================
 * Name         : speed_tests.h includes function to test
 * 				  the time-performance of basic operations.
 * Version      : 1.0.0, 7 Oct 2017
 *
 * Copyright (c) 2017 Ioannis Karagiannis
 * All rights reserved

 * This file is part of myLinearAlgebra library.

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
//    |  10x10  |       0.02  [ms]       |
//    ------------------------------------
//    |  50x50  |        0.4  [ms]       |
//    ------------------------------------
//    | 100x100 |        0.5  [ms]       |
//    |-----------------------------------
//    | 250*250 |        3.3  [ms]       |
//    |-----------------------------------
//    | 500*500 |       11.6  [ms]       |
//    ------------------------------------
//    |1000*1000|       39.2  [ms]       |
//    ------------------------------------
//    |4000*4000|       0.43  [s]        |
//    ------------------------------------
//
// MACHINE DETAILS: Intel® Core™ i3 CPU M 330 @ 2.13GHz × 4 (64-bit), 8GB RAM.


inline void test_matrix_addition_performance ( size_t size ){

	algebra::mat a = algebra::rand(size, size);
	algebra::mat b = algebra::rand(size, size);
	algebra::mat c;
	size_t iteration_nbr = 1;

	if(size > 500){
		iteration_nbr = 1;
	}else{
		iteration_nbr = 4;
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
//    |  10x10  |       0.01  [ms]       |
//    ------------------------------------
//    |  50x50  |        0.2  [ms]       |
//    ------------------------------------
//    | 100x100 |        0.5  [ms]       |
//    |-----------------------------------
//    | 250*250 |        3.5  [ms]       |
//    |-----------------------------------
//    | 500*500 |       11.9  [ms]       |
//    ------------------------------------
//    |1000*1000|         47  [ms]       |
//    ------------------------------------
//    |4000*4000|        0.8  [s]        |
//    ------------------------------------
//
// MACHINE DETAILS: Intel® Core™ i3 CPU M 330 @ 2.13GHz × 4 (64-bit), 8GB RAM.
inline void test_matrix_transpose_performance ( size_t size ){

	algebra::mat a = algebra::rand(size, size);
	algebra::mat trans;
	size_t iteration_nbr = 1;

	if(size > 500){
		iteration_nbr = 1;
	}else{
		iteration_nbr = 4;
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for(size_t i = 0; i < iteration_nbr; i++){
		trans = algebra::transpose<double>(a);
	}
	high_resolution_clock::time_point t2 = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>( t2 - t1 ).count();
	double duration_f = (double) duration; duration_f = duration_f/iteration_nbr;
	printf("transposed of a squared (%zux%zu) matrix lasted %.3f [ms]\n", trans.rows(), trans.cols(), duration_f/1000.0);

}

// IT ASSUMES THAT MATRICES WILL BE SQUARE. THAT'S WHY ONLY SIZE IS DETERMINED.
inline void test_matrix_normal_multiplication_performance ( size_t size ){

	algebra::mat a = algebra::rand(size, size);
	algebra::mat b = algebra::rand(size, size);
	algebra::mat c;
	size_t iteration_nbr = 1;

	if(size > 500){
		iteration_nbr = 1;
	}else{
		iteration_nbr = 4;
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
//    |  50x50  |        3.14  [ms]      |      352  [ms]          |
//    --------------------------------------------------------------
//    | 250x250 |        0.33  [s]       |     1.53  [s]           |
//    --------------------------------------------------------------
//    | 512x512 |         2.6  [s]       |        5  [s]           |
//    |-------------------------------------------------------------
//    |1000*1000|          21  [s]       |       22  [s]           |
//    |-------------------------------------------------------------
//    |2000*2000|         170  [s]       |      129  [s]           |
//    --------------------------------------------------------------
//    |4000*4000|        22.3  [min]     |     14.3  [min]         |
//    --------------------------------------------------------------

// MACHINE DETAILS: Intel® Core™ i3 CPU M 330 @ 2.13GHz × 4 (64-bit), 8GB RAM.
inline void test_matrix_strassen_multiplication_performance ( size_t size ){

	algebra::mat a = algebra::rand(size, size);
	algebra::mat b = algebra::rand(size, size);
	algebra::mat c;
	size_t iteration_nbr = 1;

	if(size > 500){
		iteration_nbr = 1;
	}else{
		iteration_nbr = 4;
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for(size_t i = 0; i < iteration_nbr; i++){
		c = algebra::strassen<double>(a, b);
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
//    |  10x10  |        0.1  [ms]       |
//    ------------------------------------
//    |  50x50  |        6.4  [ms]       |
//    ------------------------------------
//    | 100x100 |         49  [ms]       |
//    |-----------------------------------
//    | 250*250 |        0.7  [s]        |
//    |-----------------------------------
//    | 500*500 |        5.6  [s]        |
//    ------------------------------------
//    |1000*1000|         45  [s]        |
//    ------------------------------------
//    |4000*4000|      50.16 [min]       |
//    ------------------------------------
//
// MACHINE DETAILS: Intel® Core™ i3 CPU M 330 @ 2.13GHz × 4 (64-bit), 8GB RAM.

inline void test_matrix_inversion_performance ( size_t size ){

	algebra::mat a = algebra::rand(size, size);
	algebra::mat b = a;
	size_t iteration_nbr = 1;

	if(size > 500){
		iteration_nbr = 1;
	}else{
		iteration_nbr = 4;
	}

	high_resolution_clock::time_point t1 = high_resolution_clock::now();
	for(size_t i = 0; i < iteration_nbr; i++){
		b = algebra::inv<double>(a);
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
