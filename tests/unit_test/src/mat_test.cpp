/*====================================================================================================
 * Name         : mat_test.cpp implements a unit-test for the 'mat' class of the LinearAlgebra library.
 * Version      : 1.0.0, 23 Sep 2017
 *
 * Copyright (c) 2017 Ioannis Karagiannis
 * All rights reserved

 * This file is part of myLinearAlgebra library.

 * myLinearAlgebra is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.

 * You are free to use this library under the terms of the GNU General
 * Public License, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.

 * You should have received a copy of the GNU General Public License
 * along with myLinearAlgebra.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact info: https://www.linkedin.com/in/ioannis-karagiannis-7129394a/
 * 				ioanniskaragiannis1987@gmail.com
=====================================================================================================*/

#include "../include/catch.hpp"
#include <../../../include/base.h>

namespace algebra {

// *********************** TEST MEMBER FUNCTIONS *************************
// ***********************************************************************


TEST_CASE( " Test mat constructor." ){
	SECTION(" Test default constructor."){
		// Clear log file.
		clear_file(LOG_ERROR_FILE);
		clear_file(LOG_FILE);
		clear_file(WARNING_FILE);
		log_error(" ");
		log_error("===========================================================================================================");
		log_error("================================= NEW UNIT-TEST FOR MATRIX CLASS STARTED ==================================");
		log_error("===========================================================================================================");
		log_error(" ");
		mat m;
		INFO("Unit test failed in algebra::mat()");  // Only appears on a FAIL
		REQUIRE(m.size() == 0);
		REQUIRE(m.rows() == 0);
		REQUIRE(m.cols() == 0);
	}
	SECTION(" Test user-defined constructor."){
		mat m(2,2);
		INFO("Unit test failed in algebra::mat(size_t r, size_t c)");  // Only appears on a FAIL
		REQUIRE(m.rows() == 2);
		REQUIRE(m.cols() == 2);
		REQUIRE(m.size() == 4);
	}
	SECTION(" Test boundary conditions."){
		REQUIRE_THROWS(mat(MAX_ACCEPTABLE_VECTOR_SIZE + 1, MAX_ACCEPTABLE_VECTOR_SIZE + 1));
	}
}

TEST_CASE( " Test 'mat::set_size(size_t r, size_t c)' function." ){
	mat m(2,2);
	SECTION(" Test for normal conditions"){
		m.set_size(3, 2);
		REQUIRE(m.rows() == 3);
		REQUIRE(m.cols() == 2);
		REQUIRE(m.size() == 6);
	}
	SECTION(" Test for boundary conditions."){
		m.set_size(0, 1);
		REQUIRE(m.rows() == 0);
		REQUIRE(m.cols() == 0);
		REQUIRE(m.size() == 0);
		REQUIRE_THROWS( m.set_size(MAX_ACCEPTABLE_VECTOR_SIZE + 1, MAX_ACCEPTABLE_VECTOR_SIZE) );
		REQUIRE_THROWS( m.set_size(4, -5) );
	}
}

TEST_CASE( " Test 'mat::get(size_t r, size_t c)' function." ){
	mat m;
	SECTION(" Test for normal conditions"){
		m.set_size(2, 2);
		REQUIRE(m.get(0,0) == 0); REQUIRE(m.get(0,1) == 0);
		REQUIRE(m.get(1,0) == 0); REQUIRE(m.get(1,1) == 0);
	}
	SECTION(" Test for boundary conditions."){
		m.set_size(2, 2);
		REQUIRE_THROWS(m.get(2, 1));
		REQUIRE_THROWS(m.get(-1, 0));
		m.set_size(2, 0);
		REQUIRE_THROWS(m.get(0, 0)); // NULL MATRIX
	}
}


TEST_CASE( " Test 'mat::set(size_t r, size_t c, float value)' function." ){
	mat m(2,2);
	SECTION(" Test for normal conditions"){
		m.set(0,1, 3.510);
		REQUIRE(m.get(0,1) == Approx(3.51));
	}
	SECTION(" Test for boundary conditions."){
		REQUIRE_THROWS(m.set(2, 1, 3.4));
		REQUIRE_THROWS(m.set(0, -1, -5.7));
	}
}

TEST_CASE( " Test 'mat::set_row(size_t r, const vec& v1)' function." ){
	mat m(3,3);
	SECTION(" Test for normal conditions"){
		vec v; v = "[1 2 3]";
		m.set_row(1, v);
		REQUIRE(m.get(1,0) == 1);
		REQUIRE(m.get(1,1) == 2);
		REQUIRE(m.get(1,2) == 3);
	}
	SECTION(" Test for boundary conditions."){
		vec v; v = "[1 2 3 4]";
		REQUIRE_THROWS(m.set_row(1, v)); // dimension mismatch
		REQUIRE_THROWS(m.set_row(3, v)); // index out of bounds
		v = "[1 2]";
		REQUIRE_THROWS(m.set_row(1, v)); // dimension mismatch
	}
}

TEST_CASE( " Test 'mat::set_col(size_t r, const vec& v1)' function." ){
	mat m(3,3);
	SECTION(" Test for normal conditions"){
		vec v; v = "[1 2 3]";
		m.set_col(1, v);
		REQUIRE(m.get(0,1) == 1);
		REQUIRE(m.get(1,1) == 2);
		REQUIRE(m.get(2,1) == 3);
	}
	SECTION(" Test for boundary conditions."){
		vec v; v = "[1 2 3 4]";
		REQUIRE_THROWS(m.set_col(1, v)); // dimension mismatch
		REQUIRE_THROWS(m.set_col(3, v)); // index out of bounds
		v = "[1 2]";
		REQUIRE_THROWS(m.set_col(1, v)); // dimension mismatch
	}
}

TEST_CASE( " Test 'mat::set_rows(size_t r0, const mat& m1)' function." ){
	mat m(3,3), k; vec r0,r1;
	SECTION(" Test for normal conditions"){
		// Example 1:
		//     |0 0 0|      |1 2 3|							 |0 0 0|
		// m = |0 0 0|, k = |4 5 6|, m.set_rows(1, k) => m = |1 2 3|
		//     |0 0 0|										 |4 5 6|
		k.set_size(2,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]";
		k.set_row(0, r0); k.set_row(1, r1);
		m.set_rows(1, k);
		REQUIRE(m.get(0,0) == 0); REQUIRE(m.get(0,1) == 0); REQUIRE(m.get(0,2) == 0);
		REQUIRE(m.get(1,0) == 1); REQUIRE(m.get(1,1) == 2); REQUIRE(m.get(1,2) == 3);
		REQUIRE(m.get(2,0) == 4); REQUIRE(m.get(2,1) == 5); REQUIRE(m.get(2,2) == 6);

		// Example 2:
		//     |0 0 0|      |1 2|							 |0 0 0|
		// m = |0 0 0|, k = |3 4|, m.set_rows(1, k) =>   m = |1 2 0|
		//     |0 0 0|										 |3 4 0|
		m.set_size(3,3);
		k.set_size(2,2);
		r0 = "[1 2]"; r1 = "[3 4]";
		k.set_row(0, r0); k.set_row(1, r1);
		m.set_rows(1, k);
		REQUIRE(m.get(0,0) == 0); REQUIRE(m.get(0,1) == 0); REQUIRE(m.get(0,2) == 0);
		REQUIRE(m.get(1,0) == 1); REQUIRE(m.get(1,1) == 2); REQUIRE(m.get(1,2) == 0);
		REQUIRE(m.get(2,0) == 3); REQUIRE(m.get(2,1) == 4); REQUIRE(m.get(2,2) == 0);

		// Example 3:
		//     |0 0 0|          							 |0 0 0|
		// m = |0 0 0|, k = |1 2 3|, m.set_rows(1, k) => m = |1 2 3|
		//     |0 0 0|										 |0 0 0|
		m.set_size(3,3);
		k.set_size(1,3);
		r0 = "[1 2 3]";
		k.set_row(0, r0);
		m.set_rows(1, k);
		REQUIRE(m.get(0,0) == 0); REQUIRE(m.get(0,1) == 0); REQUIRE(m.get(0,2) == 0);
		REQUIRE(m.get(1,0) == 1); REQUIRE(m.get(1,1) == 2); REQUIRE(m.get(1,2) == 3);
		REQUIRE(m.get(2,0) == 0); REQUIRE(m.get(2,1) == 0); REQUIRE(m.get(2,2) == 0);
	}
	SECTION(" Test for boundary conditions."){
		m.set_size(3,3);
		k.set_size(2,4);
		r0 = "[1 2 3 4]"; r1 = "[5 6 7 8]";
		k.set_row(0, r0); k.set_row(1, r1);
		REQUIRE_THROWS(m.set_rows(1, k)); // dimension mismatch
		REQUIRE_THROWS(m.set_rows(3, k)); // index out of bounds
	}
}

TEST_CASE( " Test 'mat::set_cols(size_t r0, const mat& m1)' function." ){
	mat m(3,3), k; vec c0, c1;
	SECTION(" Test for normal conditions"){
		// Example 1:
		//     |0 0 0|      |1 4|							|0 1 4|
		// m = |0 0 0|, k = |2 5|, m.set_cols(1, k) =>  m = |0 2 5|
		//     |0 0 0|		|3 6|						    |0 3 6|
		k.set_size(3,2);
		c0 = "[1 2 3]"; c1 = "[4 5 6]";
		k.set_col(0, c0); k.set_col(1, c1);
		m.set_cols(1, k);
		REQUIRE(m.get(0,0) == 0); REQUIRE(m.get(0,1) == 1); REQUIRE(m.get(0,2) == 4);
		REQUIRE(m.get(1,0) == 0); REQUIRE(m.get(1,1) == 2); REQUIRE(m.get(1,2) == 5);
		REQUIRE(m.get(2,0) == 0); REQUIRE(m.get(2,1) == 3); REQUIRE(m.get(2,2) == 6);

		// Example 2:
		//     |0 0 0|      |1 2|							|0 1 2|
		// m = |0 0 0|, k = |3 4|, m.set_cols(1, k) =>  m = |0 3 4|
		//     |0 0 0|		    						    |0 0 0|
		m.set_size(3,3);
		k.set_size(2,2);
		c0 = "[1 3]"; c1 = "[2 4]";
		k.set_col(0, c0); k.set_col(1, c1);
		m.set_cols(1, k);
		REQUIRE(m.get(0,0) == 0); REQUIRE(m.get(0,1) == 1); REQUIRE(m.get(0,2) == 2);
		REQUIRE(m.get(1,0) == 0); REQUIRE(m.get(1,1) == 3); REQUIRE(m.get(1,2) == 4);
		REQUIRE(m.get(2,0) == 0); REQUIRE(m.get(2,1) == 0); REQUIRE(m.get(2,2) == 0);

		// Example 3:
		//     |0 0 0|      |1|    						   |0 1 0|
		// m = |0 0 0|, k = |2|  , m.set_cols(1, k) => m = |0 2 0|
		//     |0 0 0|		|3|    						   |0 3 0|
		m.set_size(3,3);
		k.set_size(3,1);
		c0 = "[1 2 3]";
		k.set_col(0, c0);
		m.set_cols(1, k);
		REQUIRE(m.get(0,0) == 0); REQUIRE(m.get(0,1) == 1); REQUIRE(m.get(0,2) == 0);
		REQUIRE(m.get(1,0) == 0); REQUIRE(m.get(1,1) == 2); REQUIRE(m.get(1,2) == 0);
		REQUIRE(m.get(2,0) == 0); REQUIRE(m.get(2,1) == 3); REQUIRE(m.get(2,2) == 0);
	}
	SECTION(" Test for boundary conditions."){
		k.set_size(4,2);
		m.set_size(3,3);
		c0 = "[1 2 3 4]"; c1 = "[5 6 7 8]";
		k.set_col(0, c0); k.set_col(1, c1);
		REQUIRE_THROWS(m.set_cols(1, k)); // dimension mismatch
		REQUIRE_THROWS(m.set_cols(3, k)); // index out of bounds
	}
}

TEST_CASE( " Test 'mat::set_submatrix(size_t r0, size_t c0, const mat& m)' function." ){
	mat m(3,3), k; vec c0, c1;
	SECTION(" Test for normal conditions"){
		// Example 1:
		//     |0 0 0|      |1 2|							        |0 0 0|
		// m = |0 0 0|, k = |3 4|, m.set_submatrix(1, 1, k) =>  m = |0 1 2|
		//     |0 0 0|		    						            |0 3 4|
		m.set_size(3,3);
		k.set_size(2,2);
		c0 = "[1 3]"; c1 = "[2 4]";
		k.set_col(0, c0); k.set_col(1, c1);
		m.set_submatrix(1, 1,  k);
		REQUIRE(m.get(0,0) == 0); REQUIRE(m.get(0,1) == 0); REQUIRE(m.get(0,2) == 0);
		REQUIRE(m.get(1,0) == 0); REQUIRE(m.get(1,1) == 1); REQUIRE(m.get(1,2) == 2);
		REQUIRE(m.get(2,0) == 0); REQUIRE(m.get(2,1) == 3); REQUIRE(m.get(2,2) == 4);

		// Example 2:
		//     |0 0 0|      |1 2|							        |1 2 0|
		// m = |0 0 0|, k = |3 4|, m.set_submatrix(0, 0, k) =>  m = |3 4 0|
		//     |0 0 0|		    						            |0 0 0|
		m.set_size(3,3);
		k.set_size(2,2);
		c0 = "[1 3]"; c1 = "[2 4]";
		k.set_col(0, c0); k.set_col(1, c1);
		m.set_submatrix(0, 0,  k);
		REQUIRE(m.get(0,0) == 1); REQUIRE(m.get(0,1) == 2); REQUIRE(m.get(0,2) == 0);
		REQUIRE(m.get(1,0) == 3); REQUIRE(m.get(1,1) == 4); REQUIRE(m.get(1,2) == 0);
		REQUIRE(m.get(2,0) == 0); REQUIRE(m.get(2,1) == 0); REQUIRE(m.get(2,2) == 0);
	}
	SECTION(" Test for boundary conditions."){
		m.set_size(3,3);
		k.set_size(2,2);
		c0 = "[1 3]"; c1 = "[2 4]";
		k.set_col(0, c0); k.set_col(1, c1);
		REQUIRE_THROWS(m.set_submatrix(0, 2,  k)); // k can't fit to m if we start from m(0,2) element
		REQUIRE_THROWS(m.set_submatrix(2, 2,  k)); // k can't fit to m if we start from m(2,2) element
	}
}

TEST_CASE( " Test 'mat::get_col(size_t c1) const' function." ){
	mat m(3,3); vec v, r0, r1, r2;
	SECTION(" Test for normal conditions"){
		// Example 1:
		//     |1 2 3|
		// m = |4 5 6|, m.get_col(0) = [1 4 7],  m.get_col(1) = [2 5 8], m.get_col(2) = [3 6 9]
		//     |7 8 9|
		m.set_size(3,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]"; r2 = "[7 8 9]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);

		v = m.get_col(0);
		REQUIRE(v.get(0) == 1); REQUIRE(v.get(1) == 4); REQUIRE(v.get(2) == 7);
		v = m.get_col(1);
		REQUIRE(v.get(0) == 2); REQUIRE(v.get(1) == 5); REQUIRE(v.get(2) == 8);
		v = m.get_col(2);
		REQUIRE(v.get(0) == 3); REQUIRE(v.get(1) == 6); REQUIRE(v.get(2) == 9);
	}
	SECTION(" Test for boundary conditions."){
		m.set_size(3,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]"; r2 = "[7 8 9]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);
		REQUIRE_THROWS(m.get_col(-1)); // index out of bounds
		REQUIRE_THROWS(m.get_col(3)); // index out of bounds
	}
}

TEST_CASE( " Test 'mat::get_cols(size_t c1, size_t c2) const' function." ){
	mat m(3,3), p; vec r0, r1, r2;
	SECTION(" Test for normal conditions"){
		// Example 1:
		//     |1 2 3|						  |2 3|
		// m = |4 5 6|, p = m.get_cols(1,2) = |5 6|
		//     |7 8 9|						  |8 9|
		m.set_size(3,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]"; r2 = "[7 8 9]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);

		p = m.get_cols(1,2);
		REQUIRE(p.rows() == m.rows());
		REQUIRE(p.cols() == 2);
		REQUIRE(p.get(0, 0) == 2); REQUIRE(p.get(0, 1) == 3);
		REQUIRE(p.get(1, 0) == 5); REQUIRE(p.get(1, 1) == 6);
		REQUIRE(p.get(2, 0) == 8); REQUIRE(p.get(2, 1) == 9);
	}
	SECTION(" Test for boundary conditions."){
		m.set_size(3,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]"; r2 = "[7 8 9]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);
		REQUIRE_THROWS(m.get_cols(-1, 2)); // index out of bounds
		REQUIRE_THROWS(m.get_cols(1,3)); // index out of bounds
	}
}

TEST_CASE( " Test 'mat::get_row(size_t r1) const' function." ){
	mat m(3,3); vec v, r0, r1, r2;
	SECTION(" Test for normal conditions"){
		// Example 1:
		//     |1 2 3|
		// m = |4 5 6|, m.get_row(0) = [1 2 3],  m.get_row(1) = [4 5 6], m.get_row(2) = [7 8 9]
		//     |7 8 9|
		m.set_size(3,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]"; r2 = "[7 8 9]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);

		v = m.get_row(0);
		REQUIRE(v.get(0) == 1); REQUIRE(v.get(1) == 2); REQUIRE(v.get(2) == 3);
		v = m.get_row(1);
		REQUIRE(v.get(0) == 4); REQUIRE(v.get(1) == 5); REQUIRE(v.get(2) == 6);
		v = m.get_row(2);
		REQUIRE(v.get(0) == 7); REQUIRE(v.get(1) == 8); REQUIRE(v.get(2) == 9);
	}
	SECTION(" Test for boundary conditions."){
		m.set_size(3,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]"; r2 = "[7 8 9]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);
		REQUIRE_THROWS(m.get_row(-1)); // index out of bounds
		REQUIRE_THROWS(m.get_row(3)); // index out of bounds
	}
}

TEST_CASE( " Test 'mat::get_rows(size_t r1, size_t r2) const' function." ){
	mat m(3,3), p; vec r0, r1, r2;
	SECTION(" Test for normal conditions"){
		// Example 1:
		//     |1 2 3|
		// m = |4 5 6|, p = m.get_rows(1,2) = |4 5 6|
		//     |7 8 9|						  |7 8 9|
		m.set_size(3,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]"; r2 = "[7 8 9]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);

		p = m.get_rows(1,2);
		REQUIRE(p.rows() == 2);
		REQUIRE(p.cols() == m.cols());
		REQUIRE(p.get(0, 0) == 4); REQUIRE(p.get(0, 1) == 5); REQUIRE(p.get(0, 2) == 6);
		REQUIRE(p.get(1, 0) == 7); REQUIRE(p.get(1, 1) == 8); REQUIRE(p.get(1, 2) == 9);
	}
	SECTION(" Test for boundary conditions."){
		m.set_size(3,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]"; r2 = "[7 8 9]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);
		REQUIRE_THROWS(m.get_rows(-1, 2)); // index out of bounds
		REQUIRE_THROWS(m.get_rows(1,3)); // index out of bounds
	}
}


TEST_CASE( " Test 'mat::get(size_t r1, size_t r2, size_t c1, size_t c2) const' function." ){
	mat m(3,3), p; vec r0, r1, r2;
	SECTION(" Test for normal conditions"){
		// Example 1:
		//     |1 2 3|
		// m = |4 5 6|, p = m.get(1,2,1,2)  = |5 6|
		//     |7 8 9|						  |8 9|
		m.set_size(3,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]"; r2 = "[7 8 9]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);

		p = m.get(1,2,1,2);
		REQUIRE(p.rows() == 2);
		REQUIRE(p.cols() == 2);
		REQUIRE(p.get(0, 0) == 5); REQUIRE(p.get(0, 1) == 6);
		REQUIRE(p.get(1, 0) == 8); REQUIRE(p.get(1, 1) == 9);
	}
	SECTION(" Test for boundary conditions."){
		m.set_size(3,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]"; r2 = "[7 8 9]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);
		REQUIRE_THROWS(m.get(-1, 2, 1, 2)); // index out of bounds
		REQUIRE_THROWS(m.get(1,3,1,2)); // index out of bounds
		REQUIRE_THROWS(m.get(1,0,1,2)); // r2 < r1
		REQUIRE_THROWS(m.get(1,2,2,1)); // c2 < c1
	}
}

TEST_CASE( " Test 'mat::get(const vec& r, const vec& c) const' function." ){
	mat m, p; vec r0, r1, r2, r, c;
	SECTION(" Test for normal conditions"){
		// Example 1:
		//     |1 2 3|
		// m = |4 5 6|, r = [0 2], c = |0 2| == > m.get(r,c) = |1 3|
		//     |7 8 9|						                   |7 9|
		m.set_size(3,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]"; r2 = "[7 8 9]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);

		r = "[0 2]"; c = "[0 2]";
		p = m.get(r, c);
		REQUIRE(p.rows() == r.size());
		REQUIRE(p.cols() == c.size());
		REQUIRE(p.get(0, 0) == 1); REQUIRE(p.get(0, 1) == 3);
		REQUIRE(p.get(1, 0) == 7); REQUIRE(p.get(1, 1) == 9);

		// Example 2:
		//     |1 2 3|
		// m = |4 5 6|, r = [0 2], c = |0 1| == > m.get(r,c) = |1 2|
		//     |7 8 9|						                   |7 8|
		m.set_size(3,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]"; r2 = "[7 8 9]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);

		r = "[0 2]"; c = "[0 1]";
		p = m.get(r, c);
		REQUIRE(p.rows() == r.size());
		REQUIRE(p.cols() == c.size());
		REQUIRE(p.get(0, 0) == 1); REQUIRE(p.get(0, 1) == 2);
		REQUIRE(p.get(1, 0) == 7); REQUIRE(p.get(1, 1) == 8);

		// Example 3:
		//     |1 2 3|											   |1 2 3|
		// m = |4 5 6|, r = [0 1 2], c = |0 1 2| == > m.get(r,c) = |4 5 6|
		//     |7 8 9|						                       |7 8 9|
		m.set_size(3,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]"; r2 = "[7 8 9]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);

		r = "[0 1 2]"; c = "[0 1 2]";
		p = m.get(r, c);
		REQUIRE(p.rows() == r.size());
		REQUIRE(p.cols() == c.size());
		REQUIRE(p.get(0, 0) == 1); REQUIRE(p.get(0, 1) == 2); REQUIRE(p.get(0,2) == 3);
		REQUIRE(p.get(1, 0) == 4); REQUIRE(p.get(1, 1) == 5); REQUIRE(p.get(1,2) == 6);
		REQUIRE(p.get(2, 0) == 7); REQUIRE(p.get(2, 1) == 8); REQUIRE(p.get(2,2) == 9);
	}
	SECTION(" Test for boundary conditions."){
		m.set_size(3,3);
		r0 = "[1 2 3]"; r1 = "[4 5 6]"; r2 = "[7 8 9]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);

		// Example 1:
		r = "[0 1 2 3]"; c = "[1 2]";
		REQUIRE_THROWS(m.get(r,c)); // row-index vector exceeds matrix row dimension
		// Example 2:
		r = "[0 1 2]"; c = "[1 2 3 4]";
		REQUIRE_THROWS(m.get(r,c)); // col-index vector exceeds matrix col dimension
		// Example 3:
		r = "[]"; c = "[1 2]";
		REQUIRE( (m.get(r,c)).size() == 0); // index out of bounds
	}
}

TEST_CASE( " Test 'mat::zeros()' " ){
	mat m(2,2);
	m.set(0, 0, 44); m.set(0,1,55); m.set(1,1, 32);
	m.zeros();
	REQUIRE(m.get(0,0) == 0); REQUIRE(m.get(0,1) == 0);
	REQUIRE(m.get(1,0) == 0); REQUIRE(m.get(1,1) == 0);
}

TEST_CASE( " Test 'mat::clear()' " ){
	mat m(2,2);
	m.set(0, 0, 44); m.set(0,1,55); m.set(1,1, 32);
	m.clear();
	REQUIRE(m.get(0,0) == 0); REQUIRE(m.get(0,1) == 0);
	REQUIRE(m.get(1,0) == 0); REQUIRE(m.get(1,1) == 0);
}

TEST_CASE( " Test 'mat::ones()' " ){
	mat m(2,2);
	m.set(0, 0, 44); m.set(0,1,55); m.set(1,1, 32);
	m.ones();
	REQUIRE(m.get(0,0) == 1); REQUIRE(m.get(0,1) == 1);
	REQUIRE(m.get(1,0) == 1); REQUIRE(m.get(1,1) == 1);
}

TEST_CASE( " Test 'mat::swap_rows(size_t i, size_t j)' " ){
	imat m, tmp;
	SECTION( "Test normal conditions" ){
		m = rand_i(3,3);
		tmp = m;
		m.swap_rows(0,1);
		REQUIRE(m(0,0) == tmp(1,0)); REQUIRE(m(0,1) == tmp(1,1)); REQUIRE(m(0,2) == tmp(1,2));
		REQUIRE(m(1,0) == tmp(0,0)); REQUIRE(m(1,1) == tmp(0,1)); REQUIRE(m(1,2) == tmp(0,2));
		REQUIRE(m(2,0) == tmp(2,0)); REQUIRE(m(2,1) == tmp(2,1)); REQUIRE(m(2,2) == tmp(2,2));
	}
	SECTION( "Test boundary conditions." ){
		m = rand_i(4,4);
		REQUIRE_THROWS( m.swap_rows(0, 4) ); // j > cols
		REQUIRE_THROWS( m.swap_rows(-1, 2) );// i < 0
	}

}

TEST_CASE( " Test 'mat::swap_cols(size_t i, size_t j)' " ){
	imat m, tmp;
	SECTION( "Test normal conditions" ){
		m = rand_i(3,3);
		tmp = m;
		m.swap_cols(0,2);
		REQUIRE(m(0,0) == tmp(0,2)); REQUIRE(m(0,1) == tmp(0,1)); REQUIRE(m(0,2) == tmp(0,0));
		REQUIRE(m(1,0) == tmp(1,2)); REQUIRE(m(1,1) == tmp(1,1)); REQUIRE(m(1,2) == tmp(1,0));
		REQUIRE(m(2,0) == tmp(2,2)); REQUIRE(m(2,1) == tmp(2,1)); REQUIRE(m(2,2) == tmp(2,0));
	}
	SECTION( "Test boundary conditions." ){
		m = rand_i(4,4);
		REQUIRE_THROWS( m.swap_cols(0, 4) ); // j > cols
		REQUIRE_THROWS( m.swap_cols(-1, 2) );// i < 0
	}

}

// *********************** TEST OVERLOADED OPERATORS *************************
// ***************************************************************************

TEST_CASE( " Test 'mat::operator=(const char* a)' " ){
	mat m;
	SECTION("Test normal conditions"){
		// Example 1 (Square Matrix)
		//	   							|1 2 3|
		// m = [1 2 3;4 5 6;7 8 9]	=	|4 5 6|
		//     							|7 8 9|

		m = "[1 2 3;4 5 6;7 8 9]";
		REQUIRE(m.rows() == 3); REQUIRE(m.cols() == 3);
		REQUIRE(m.get(0,0) == 1); REQUIRE(m.get(0,1) == 2); REQUIRE(m.get(0,2) == 3);
		REQUIRE(m.get(1,0) == 4); REQUIRE(m.get(1,1) == 5); REQUIRE(m.get(1,2) == 6);
		REQUIRE(m.get(2,0) == 7); REQUIRE(m.get(2,1) == 8); REQUIRE(m.get(2,2) == 9);

		// Example 2 (Non-Square Matrix)
		//
		// m = [1 2 3;4 5 6]	=	|1 2 3|
		//     						|4 5 6|
		m = "[1 2 3;4 5 6]";
		REQUIRE(m.rows() == 2); REQUIRE(m.cols() == 3);
		REQUIRE(m.get(0,0) == 1); REQUIRE(m.get(0,1) == 2); REQUIRE(m.get(0,2) == 3);
		REQUIRE(m.get(1,0) == 4); REQUIRE(m.get(1,1) == 5); REQUIRE(m.get(1,2) == 6);

		// Example 3 (Vectorwise Matrix)
		//
		// m = [1 2 3]	=	|1 2 3|
		//     				|4 5 6|
		m = "[1 2 3]";
		REQUIRE(m.rows() == 1); REQUIRE(m.cols() == 3);
		REQUIRE(m.get(0,0) == 1); REQUIRE(m.get(0,1) == 2); REQUIRE(m.get(0,2) == 3);
	}
	SECTION("Test boundary conditions."){
		// Example 1 (Empty matrix)
		// m = []	=	||
		m = "[]";
		REQUIRE(m.rows() == 0); REQUIRE(m.cols() == 0);
		REQUIRE_THROWS( m.get(0,0) );

		// Example 2 (Input-string with numbers and letters)
		// m = [1 a 4; 2 4 5]	=	exception
		REQUIRE_THROWS( m = "[1 a 4; 2 4 5]");

		// Example 3 (Input-string with letter only)
		// m = [a b c; d e f]	=	||
		m = "[a b c; d e f]";
		REQUIRE(m.rows() == 0); REQUIRE(m.cols() == 0);
		REQUIRE_THROWS( m.get(0,0) );

		// Example 4 (Input-string with rows with different amount of elements)
		// m = [1 2 3;4 5]	=	exception
		REQUIRE_THROWS( m = "[1 2 3;4 5]" );
	}
}


TEST_CASE( " Test 'mat::operator()(size_t i, size_t j)' " ){
	mat m;
	SECTION("Test normal conditions"){
		// Example 1
		//		|1 2 3|
		//   m =|4 5 6| ,
		//
		m = "[1 2 3;4 5 6]";
		REQUIRE(m(0,0) == 1); REQUIRE(m(0,1) == 2); REQUIRE(m(0,2) == 3);
		REQUIRE(m(1,0) == 4); REQUIRE(m(1,1) == 5); REQUIRE(m(1,2) == 6);
	}
	SECTION("Test boundary conditions."){
		// Example 1 (access NULL VECTOR)
		REQUIRE_THROWS(m(0,0));

		// Example 2 (out of range index)
		m = "[1 2 3;4 5 6]";
		REQUIRE_THROWS(m(2,1)); // i > rows
		REQUIRE_THROWS(m(1,3)); // j > cols
		REQUIRE_THROWS(m(-1,1));// i < 0
		REQUIRE_THROWS(m(1,-1));// j < 0
	}
}

TEST_CASE( " Test 'mat::operator()(size_t r1, size_t r2, size_t c1, size_t c2)' " ){
	mat m, p;
	SECTION("Test normal conditions"){
		// Example 1
		//		|1 2 3|
		//   m =|4 5 6| , p = m.get(0, 0, 1, 2) = [2 3]
		//
		m = "[1 2 3;4 5 6]";
		p = m(0, 0, 1, 2);
		REQUIRE( p.rows() == 1 ); REQUIRE( p.cols() == 2 );
		REQUIRE(p(0,0) == 2); REQUIRE(p(0,1) == 3);
	}
	SECTION("Test boundary conditions."){
		// Example 1 (access NULL VECTOR)
		REQUIRE_THROWS( p = m(0, 0, 1, 2) );
		m = "[1 2 3;4 5 6]";
		REQUIRE_THROWS( p = m(2,1,1,1) ); // r2 < r1
		REQUIRE_THROWS( p = m(0,1,2,1) ); // c2 < c1
		// no further testing required
		// it's already tested in get(size_t r1, size_t r2, size_t c1, size_t c2)
	}
}

TEST_CASE( " Test 'mat::operator+(const mat& m)' " ){
	mat m,p,b;
	SECTION("Test normal conditions"){
		// Example 1
		//		|1 2 3|		  |-4 -2 1|				 |-3 0 4 |
		//   m =|4 5 6| , p = | 2  8 7|, b = m + p = |6 13 13|
		//
		m = "[1 2 3;4 5 6]";
		p = "[-4 -2 1;2 8 7]";
		b = m + p;
		REQUIRE( b.rows() == 2 ); REQUIRE( b.cols() == 3 );
		REQUIRE(b(0,0) == -3); REQUIRE(b(0,1) == 0); REQUIRE(b(0,2) == 4);
		REQUIRE(b(1,0) == 6); REQUIRE(b(1,1) == 13); REQUIRE(b(1,2) == 13);
	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE_THROWS( b = m + p );
		// Example 2: dimension mismatch
		m = "[1 2 3;4 5 6]";
		p = "[-4 -2;2 8]";
		REQUIRE_THROWS( b = m + p );
	}
}

TEST_CASE( " Test 'mat::operator+(float t)' " ){
	mat m,b; float t;
	SECTION("Test normal conditions"){
		// Example 1
		//		|1 2 3|		  				   |11 12 13|
		//   m =|4 5 6| , t = 10, b = m + t =  |14 15 16|
		//
		m = "[1 2 3;4 5 6]";
		t = 10;
		b = m + t;
		REQUIRE( b.rows() == 2 ); REQUIRE( b.cols() == 3 );
		REQUIRE(b(0,0) == 11); REQUIRE(b(0,1) == 12); REQUIRE(b(0,2) == 13);
		REQUIRE(b(1,0) == 14); REQUIRE(b(1,1) == 15); REQUIRE(b(1,2) == 16);
	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE_THROWS( b = m + t );
	}
}

TEST_CASE( " Test 'mat::operator-(const mat& m)' " ){
	mat m,p,b;
	SECTION("Test normal conditions"){
		// Example 1
		//		|1 2 3|		  |-4 -2 1|				 |5 4  2|
		//   m =|4 5 6| , p = | 2  8 7|, b = m - p = |2 -3 -1|
		//
		m = "[1 2 3;4 5 6]";
		p = "[-4 -2 1;2 8 7]";
		b = m - p;
		REQUIRE( b.rows() == 2 ); REQUIRE( b.cols() == 3 );
		REQUIRE(b(0,0) == 5); REQUIRE(b(0,1) == 4); REQUIRE(b(0,2) == 2);
		REQUIRE(b(1,0) == 2); REQUIRE(b(1,1) == -3); REQUIRE(b(1,2) == -1);
	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE_THROWS( b = m - p );
		// Example 2: dimension mismatch
		m = "[1 2 3;4 5 6]";
		p = "[-4 -2;2 8]";
		REQUIRE_THROWS( b = m - p );
	}
}

TEST_CASE( " Test 'mat::operator-(float t)' " ){
	mat m,b; float t;
	SECTION("Test normal conditions"){
		// Example 1
		//		|1 2 3|		  				   |-9 -8 -7|
		//   m =|4 5 6| , t = 10, b = m - t =  |-6 -5 -4|
		//
		m = "[1 2 3;4 5 6]";
		t = 10;
		b = m - t;
		REQUIRE( b.rows() == 2 ); REQUIRE( b.cols() == 3 );
		REQUIRE(b(0,0) == -9); REQUIRE(b(0,1) == -8); REQUIRE(b(0,2) == -7);
		REQUIRE(b(1,0) == -6); REQUIRE(b(1,1) == -5); REQUIRE(b(1,2) == -4);
	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE_THROWS( b = m - t );
	}
}

TEST_CASE( " Test 'mat::operator*(float t)' " ){
	mat m,b; float t;
	SECTION("Test normal conditions"){
		// Example 1
		//		|1 2 3|		  				   |10 20 30|
		//   m =|4 5 6| , t = 10, b = m * t =  |40 50 60|
		//
		m = "[1 2 3;4 5 6]";
		t = 10;
		b = m * t;
		REQUIRE( b.rows() == 2 ); REQUIRE( b.cols() == 3 );
		REQUIRE(b(0,0) == 10); REQUIRE(b(0,1) == 20); REQUIRE(b(0,2) == 30);
		REQUIRE(b(1,0) == 40); REQUIRE(b(1,1) == 50); REQUIRE(b(1,2) == 60);
	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE_THROWS( b = m * t );
	}
}

TEST_CASE( " Test 'mat::operator*(const mat& m)' " ){
	mat m,p,b;
	SECTION("Test normal conditions"){
		// Example 1
		//		|1 2 3|		  |-4  1|			   |3 38|
		//   m =|4 5 6| , p = | 2  8|, b = m * p = |0 86|
		//					  | 1  7|
		m = "[1 2 3;4 5 6]";
		p = "[-4 1;2 8;1 7]";
		b = m * p;
		REQUIRE( b.rows() == m.rows() ); REQUIRE( b.cols() == p.cols() );
		REQUIRE(b(0,0) == 3); REQUIRE(b(0,1) == 38);
		REQUIRE(b(1,0) == 0); REQUIRE(b(1,1) == 86);

		// Example 2
		//		|1|		  			  			 |2 8 -3|
		//   m =|2| , p = |2  8 -3|, b = m * p = |4 16 -6|
		//		|3|			  					 |1 24 -9|

		m = "[1;2;3]";
		p = "[2 8 -3]";
		b = m * p;
		REQUIRE( b.rows() == m.rows() ); REQUIRE( b.cols() == p.cols() );
		REQUIRE(b(0,0) == 2); REQUIRE(b(0,1) == 8); REQUIRE(b(0,2) == -3);
		REQUIRE(b(1,0) == 4); REQUIRE(b(1,1) == 16); REQUIRE(b(1,2) == -6);
		REQUIRE(b(2,0) == 6); REQUIRE(b(2,1) == 24); REQUIRE(b(2,2) == -9);
	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE_THROWS( b = m * p );
		// Example 2: dimension mismatch
		m = "[1 2 3;4 5 6]";
		p = "[-4 -2;2 8]";
		REQUIRE_THROWS( b = m * p );
	}
}

TEST_CASE( " Test 'mat::operator*(const vec& v)' " ){
	mat m; vec a, b;
	SECTION("Test normal conditions"){
		// Example 1
		//		|1 2 3|		  |-1|				 |8 |
		//   m =|4 5 6| , a = |0 | , b = m * a = |14|
		//					  |3 |
		m = "[1 2 3;4 5 6]";
		a = "[-1 0 3]";
		b = m * a;
		REQUIRE( b.size() == m.rows() );
		REQUIRE(b(0) == 8);
		REQUIRE(b(1) == 14);
	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE_THROWS( b = m * a );
	}
}

TEST_CASE( " Test 'mat::operator/(float t)' " ){
	mat m,b; float t;
	SECTION("Test normal conditions"){
		// Example 1
		//		|1 2 3|		  				   |0.1 0.2 0.3|
		//   m =|4 5 6| , t = 10, b = m / t =  |0.4 0.5 0.6|
		//
		m = "[1 2 3;4 5 6]";
		t = 10;
		b = m / t;
		REQUIRE( b.rows() == m.rows() ); REQUIRE( b.cols() == m.cols() );
		REQUIRE(b(0,0) == Approx(0.1)); REQUIRE(b(0,1) == Approx(0.2)); REQUIRE(b(0,2) == Approx(0.3));
		REQUIRE(b(1,0) == Approx(0.4)); REQUIRE(b(1,1) == Approx(0.5)); REQUIRE(b(1,2) == Approx(0.6));
	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE_THROWS( b = m / t );
	}
}

// *************************** TEST FRIEND FUNCTIONS *************************
// ***************************************************************************

TEST_CASE( " Test 'strassen(const mat& m)' " ){
	mat a, b, c, c_test;
	SECTION("Test normal conditions"){
		// Example 1: Small sized square matrices
		//		|1 2 3|		  |1 -2 3|
		//   a =|4 5 6| , b = |4 -5 2|,
		//		|7 8 9|		  |-9 -2 8|
		//
		//						  |-18 -18 31|
		//   c = strassen(a, b) = |-30 -45 70|
		//						  |-42 -72 109|

		a = "[1 2 3;4 5 6;7 8 9]";
		b = "[1 -2 3;4 -5 2;-9 -2 8]";
		c = strassen(a, b);

		REQUIRE( c.rows() == a.rows() ); REQUIRE( c.cols() == a.cols() );
		REQUIRE( c(0, 0) == -18 ); REQUIRE( c(0, 1) == -18 ); REQUIRE( c(0, 2) == 31 );
		REQUIRE( c(1, 0) == -30 ); REQUIRE( c(1, 1) == -45 ); REQUIRE( c(1, 2) == 70 );
		REQUIRE( c(2, 0) == -42 ); REQUIRE( c(2, 1) == -72 ); REQUIRE( c(2, 2) == 109 );

		// Example 2: Square matrices with larger size
		// One can test the speed for 1000 size and above by calling the
		// function test_matrix_strassen_multiplication_performance ( size ).

		a = eye(40)*0.5;
		b = eye(40)*8;
		c = strassen(a, b);

		REQUIRE( c.rows() == a.rows() ); REQUIRE( c.cols() == a.cols() );

		size_t i,j, size = c.rows();
		for(i = size; i--;){
			for(j = size; j--;){
				if(i == j){
					REQUIRE( c(i,j) == Approx(4.0) );
				}else{
					REQUIRE( c(i,j) == 0 );
				}
			}
		}
	}
	SECTION("Test boundary conditions."){
		// Example 1: NON-SQUARE MATRIX
		a = ones(40, 30);
		b = ones(30, 60);
		REQUIRE_THROWS( strassen(a, b) );
	}
}

TEST_CASE( " Test 'transpose(const mat& a)' " ){
	mat m, m_t; vec r0, r1, r2;
	SECTION("Test normal conditions."){
		//Example 1 (Square matrix):
		//		|1 7 14 |						 |1 3 -1  |
		//  m = |3 12 99| , m_t = m.tranpose() = |7 12 8  |
		//		|-1 8 56|						 |14 99 56|

		m.set_size(3,3);
		r0 = "[1 7 14]"; r1 = "[3 12 99]"; r2 = "[-1 8 56]";
		m.set_row(0, r0); m.set_row(1, r1); m.set_row(2, r2);
		m_t = transpose(m);
		REQUIRE( m_t.get(0,0) == 1 ); REQUIRE( m_t.get(0,1) == 3 ); REQUIRE( m_t.get(0,2) == -1 );
		REQUIRE( m_t.get(1,0) == 7 ); REQUIRE( m_t.get(1,1) == 12 ); REQUIRE( m_t.get(1,2) == 8 );
		REQUIRE( m_t.get(2,0) == 14 ); REQUIRE( m_t.get(2,1) == 99 ); REQUIRE( m_t.get(2,2) == 56 );
		REQUIRE( m_t.rows() == m.cols() ); REQUIRE( m_t.cols() == m.rows() );

		//Example 2 (Non-Square matrix):
		//		|1 7 14 |						 |1 3  |
		//  m = |3 12 99| , m_t = m.tranpose() = |7 12 |
		//								         |14 99|

		m.set_size(2,3);
		r0 = "[1 7 14]"; r1 = "[3 12 99]";
		m.set_row(0, r0); m.set_row(1, r1);
		m_t = transpose(m);
		REQUIRE( m_t.get(0,0) == 1 ); REQUIRE( m_t.get(0,1) == 3 );
		REQUIRE( m_t.get(1,0) == 7 ); REQUIRE( m_t.get(1,1) == 12 );
		REQUIRE( m_t.get(2,0) == 14 ); REQUIRE( m_t.get(2,1) == 99 );
		REQUIRE( m_t.rows() == m.cols() ); REQUIRE( m_t.cols() == m.rows() );

		//Example 3 (Non-Square matrix):
		//		|1|
		//  m = |3| , m_t = m.tranpose() = |1 3 5|
		//		|5|

		m.set_size(3,1);
		r0 = "[1 3 5]";
		m.set_col(0, r0);
		m_t = transpose(m);
		REQUIRE( m_t.get(0,0) == 1 ); REQUIRE( m_t.get(0,1) == 3 ); REQUIRE( m_t.get(0,2) == 5 );
		REQUIRE( m_t.rows() == m.cols() ); REQUIRE( m_t.cols() == m.rows() );
	}
	SECTION("Test boundary conditions."){
		m.set_size(0,0);
		m_t = transpose(m);
		REQUIRE(m_t.size() == 0);
	}
}

TEST_CASE( " Test 'inv(const mat& a)' " ){
	mat m, m_inv, confirm;
	SECTION("Test normal conditions."){
		//Example 1 (3x3 matrix):
		//		|7  2  1|					  |-2 8  -5|					  |1 0 0|
		//  m = |0  3 -1| , m_inv = inv(m)  = |3 -11  7|, m*m_inv = m_inv*m = |0 1 0|
		//		|-3 4 -2|					  |9 -34 21|					  |0 0 1|

		m = "[7 2 1;0 3 -1;-3 4 -2]";
		m_inv = inv(m);
		confirm = m*m_inv;
		size_t i, j, size = m.rows();
		for(i = size; i--;){
			for(j = size; j--;){
				if(i ==j){
					REQUIRE( confirm(i, j) == Approx(1.0) );
				}else{
					REQUIRE( confirm(i, j) == Approx(0.0) );
				}
			}
		}

	}
	//Example 2 (5x5 matrix):
	//		|1  2  0 -8  1|
	//  m = |-2 3  4  0 -7|
	//		|0  0  9  8  0|
	//		|0  2 17 32 -4|
	//		|44 0 -5  0 -6|
	//
	m = "[1 2 0 -8 1;-2 3 4 0 -7;0 0 9 8 0;0 2 17 32 -4;44 0 -5 0 -6]";
	m_inv = inv(m);
	confirm = m*m_inv;
	size_t i, j, size = m.rows();
	for(i = size; i--;){
		for(j = size; j--;){
			if(i ==j){
				REQUIRE( confirm(i, j) == Approx(1.0) );
			}else{
				REQUIRE( confirm(i, j) == Approx(0.0) );
			}
		}
	}

	SECTION("Test boundary conditions."){
		m = "[7 2 1;0 3 -1]"; // non-sqaure
		REQUIRE_THROWS( inv(m) );
		//Example with singular matrix
		//		|1  0  0|
		//  m = |-2 0  0|
		//		|4  6  1|
		m = "[1 0 0;-2 0 0;4 6 1]";
		m_inv = inv(m); // a log should be saved with a warning
		size_t i, j;
		for(i = m.rows(); i--;){
			for(j = m.cols(); j--;){
				REQUIRE( isnan(m_inv(i,j)) == 1 );
			}
		}
	}
}

TEST_CASE( " Test 'pinv(const mat& a)' " ){
	mat m, m_pinv, confirm;
	SECTION("Test normal conditions."){
		//Example 1 (2x3 matrix and right inverse):
		//		|1 1 1 1|					   |2    -0.25|
		//  m = |5 7 7 9| , m_inv = pinv(m)  = |0.25     0| (Right Inverse)
		//							           |0.25     0|
		//									   |-1.5  0.25|
		m = "[1 1 1 1;5 7 7 9]";
		m_pinv = pinv(m);
		confirm = m*m_pinv;
		size_t i, j, size = confirm.rows();
		for(i = size; i--;){
			for(j = size; j--;){
				if(i ==j){
					REQUIRE( confirm(i, j) == Approx(1.0) );
				}else{
					REQUIRE( confirm(i, j) == Approx(0.0) );
				}
			}
		}

		//Example 2 (3x2 matrix and left inverse):
		//		|-3 -4|
		//  m = | 4  6| , m_inv = pinv(m) = (1/9)*|-11 -10 16| (Left Inverse)
		//		| 1	 1|				              | 7   8 -11|
		//
		m = "[-3 -4;4 6;1 1]";
		m_pinv = pinv(m);
		confirm = m_pinv*m;
		size = confirm.rows();
		for(i = size; i--;){
			for(j = size; j--;){
				if(i ==j){
					REQUIRE( confirm(i, j) == Approx(1.0) );
				}else{
					REQUIRE( confirm(i, j) == Approx(0.0) );
				}
			}
		}

	}

	SECTION("Test boundary conditions."){
		m = "[1 0 1 0;1 0 1 0]"; // neither full row rank, nor full column rank
		m_pinv = pinv(m);
	}
}


// *************************** TEST INLINE FUNCTIONS *************************
// ***************************************************************************

TEST_CASE( " Test 'mat::mat2vec(const mat& m)' " ){
	mat m; vec v;
	SECTION("Test normal conditions"){
		// Example 1
		//		|1 2 3|
		//   m =|4 5 6| , v = mat2vec(m) = [1 2 3 4 5 6]
		//
		m = "[1 2 3;4 5 6]";
		v = mat2vec(m);
		REQUIRE( v.size() == m.size());
		REQUIRE(v(0) == 1); REQUIRE(v(1) == 2); REQUIRE(v(2) == 3);
		REQUIRE(v(3) == 4); REQUIRE(v(4) == 5); REQUIRE(v(5) == 6);
	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE( mat2vec(m).size() == 0);
	}
}


TEST_CASE( " Test 'mat::max(const mat& m)' " ){
	mat m;
	SECTION("Test normal conditions"){
		// Example 1
		//		|1 -2 3|
		//   m =|4 77 6| , t = max(m) = 77
		//
		m = "[1 -2 3;4 77 6]";
		REQUIRE( max(m) == 77);
	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE_THROWS( max(m) );
	}
}

TEST_CASE( " Test 'mat::min(const mat& m)' " ){
	mat m;
	SECTION("Test normal conditions"){
		// Example 1
		//		|1 -2 3|
		//   m =|4 77 6| , t = min(m) = -2
		//
		m = "[1 -2 3;4 77 6]";
		REQUIRE( min(m) == -2);
	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE_THROWS( min(m) );
	}
}

TEST_CASE( " Test 'mat::abs(const mat& m)' " ){
	mat m, t;
	SECTION("Test normal conditions"){
		// Example 1
		//		|1 -2 -3|				 |1  2 3|
		//   m =|4 77 -6| , t = abs(m) = |4 77 6|
		//
		m = "[1 -2 -3;4 77 -6]";
		t = abs(m);
		REQUIRE( t(0,0) == 1); REQUIRE( t(0,1) == 2); REQUIRE( t(0,2) == 3);
		REQUIRE( t(1,0) == 4); REQUIRE( t(1,1) == 77); REQUIRE( t(1,2) == 6);

	}
}

TEST_CASE( " Test algebra::find_non_zero(const mat& m) function" ){
	mat m1, m2;
	SECTION(" Test normal conditions. "){
		m1 = "[3 4 0;0 1 0;-9 0 5]";
		m2 = find_non_zero(m1);
		size_t i, j, rows = m1.rows(), cols = m1.cols();
		for(i = rows; i--;){
			for(j = cols; j--;){
				if(m1(i,j) != 0){
					REQUIRE(m2(i,j) == 1);
				}else{
					REQUIRE(m2(i,j) == 0);
				}
			}
		}
	}
}

TEST_CASE( " Test algebra::find_zero(const mat& m) function" ){
	mat m1, m2;
	SECTION(" Test normal conditions. "){
		m1 = "[3 4 0;0 1 0;-9 0 5]";
		m2 = find_zero(m1);
		size_t i, j, rows = m1.rows(), cols = m1.cols();
		for(i = rows; i--;){
			for(j = cols; j--;){
				if(m1(i,j) == 0){
					REQUIRE(m2(i,j) == 1);
				}else{
					REQUIRE(m2(i,j) == 0);
				}
			}
		}
	}
}

TEST_CASE( " Test 'mat<int>::rand_symmetric(size_t n)' " ){
	imat m;
	SECTION("Test normal conditions"){
		size_t i, j, n = 4;
		m = rand_symmetric_i(n);
		for(i = n; i--;){
			for(j = n; j--;){
				if(i != j){
					REQUIRE( m(i,j) == m(j,i) );
				}
			}
		}
	}
	SECTION("Test boundary conditions."){
		REQUIRE_THROWS( m = rand_symmetric_i(MAX_ACCEPTABLE_VECTOR_SIZE + 1) );
	}
}

TEST_CASE( " Test 'mat<double>::rand_symmetric(size_t n)' " ){
	mat m;
	SECTION("Test normal conditions"){
		size_t i, j, n = 4;
		m = rand_symmetric(n);
		for(i = n; i--;){
			for(j = n; j--;){
				if(i != j){
					REQUIRE( m(i,j) == m(j,i) );
				}
			}
		}
	}
	SECTION("Test boundary conditions."){
		REQUIRE_THROWS( m = rand_symmetric(MAX_ACCEPTABLE_VECTOR_SIZE + 1) );
	}
}


TEST_CASE( " Test 'mat::zeros(size_t n, size_t m)' " ){
	mat m;
	SECTION("Test normal conditions"){

		m = zeros(3,2);
		REQUIRE( m.rows() == 3 ); REQUIRE( m.cols() == 2 );
		for(size_t i = 0; i < m.rows(); i++){
			for(size_t j = 0; j < m.cols(); j++){
				REQUIRE( m(i,j) == 0 );
			}
		}

	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE_THROWS( m = zeros(-3,1) );
		REQUIRE_THROWS( m = zeros(3,-1) );
		REQUIRE_THROWS( m = zeros(MAX_ACCEPTABLE_VECTOR_SIZE, MAX_ACCEPTABLE_VECTOR_SIZE + 1) );
	}
}

TEST_CASE( " Test 'mat::ones(size_t n, size_t m)' " ){
	mat m;
	SECTION("Test normal conditions"){

		m = ones(3,2);
		REQUIRE( m.rows() == 3 ); REQUIRE( m.cols() == 2 );
		for(size_t i = 0; i < m.rows(); i++){
			for(size_t j = 0; j < m.cols(); j++){
				REQUIRE( m(i,j) == 1 );
			}
		}

	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE_THROWS( m = ones(-3,1) );
		REQUIRE_THROWS( m = ones(3,-1) );
		REQUIRE_THROWS( m = ones(MAX_ACCEPTABLE_VECTOR_SIZE, MAX_ACCEPTABLE_VECTOR_SIZE + 1) );
	}
}

TEST_CASE( " Test 'mat::eye(size_t k)' " ){
	mat m;
	SECTION("Test normal conditions"){

		m = eye(3);
		REQUIRE( m.rows() == 3 ); REQUIRE( m.cols() == 3 );
		for(size_t i = 0; i < m.rows(); i++){
			REQUIRE( m(i,i) == 1 );
		}

	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE_THROWS( m = eye(-3) );
		REQUIRE_THROWS( m = eye(MAX_ACCEPTABLE_VECTOR_SIZE + 1) );
	}
}

TEST_CASE( " Test 'mat::diag(const mat& m)' " ){
	mat m; vec v;
	SECTION("Test normal conditions"){
		// Example 1:
		//     |1 2 3|
		// m = |4 5 6|, v = diag(m) = [1 5 9]
		//     |7 8 9|
		m = "[1 2 3;4 5 6;7 8 9]";
		v = diag(m);
		REQUIRE( v.size() == m.rows() ); REQUIRE( v.size() == m.cols() );
		for(size_t i = 0; i < v.size(); i++){
			REQUIRE( v(i) == m(i, i) );
		}

	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE( diag(m).size() == 0 );
	}
}

TEST_CASE( " Test 'mat::diag(const vec& v)' " ){
	mat m; vec v;
	SECTION("Test normal conditions"){
		// Example 1:
		//     						  |1 0 0|
		// v = [1 5 9], m = diag(v) = |0 5 0|
		//                            |0 0 9|
		v = "[1 5 9]";
		m = diag(v);
		REQUIRE( m.rows() == v.size() ); REQUIRE( m.cols() == v.size() );
		for(size_t i = 0; i < m.rows(); i++){
			for(size_t j = 0; j < m.cols(); j++){
				if(i == j){
					REQUIRE( m(i, j) == v(i) );
				}else{
					REQUIRE( m(i, j) == 0 );
				}
			}
		}

	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE( diag(v).size() == 0 );
	}
}

TEST_CASE( " Test 'mat::concat_hor(const mat& m1, const mat& m2)' " ){
	mat a,b,c;
	SECTION("Test normal conditions"){
		// Example 1:
		//     |2 3 4|		|9 4|					     |2 3 4 9 4|
		// a = |1 5 9|, b = |0 5|, c = concat_hor(a,b) = |1 5 9 0 5|
		//
		a = "[2 3 4;1 5 9]";
		b = "[9 4;0 5]";
		c = concat_hor(a, b);
		REQUIRE( c.rows() == a.rows() ); REQUIRE( c.cols() == a.cols() + b.cols() );
		REQUIRE(c(0,0) == 2); REQUIRE(c(0,1) == 3); REQUIRE(c(0,2) == 4); REQUIRE(c(0,3) == 9); REQUIRE(c(0,4) == 4);
		REQUIRE(c(1,0) == 1); REQUIRE(c(1,1) == 5); REQUIRE(c(1,2) == 9); REQUIRE(c(1,3) == 0); REQUIRE(c(1,4) == 5);
	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE( concat_hor(a, b).size() == 0 );
	}
}

TEST_CASE( " Test 'mat::concat_ver(const mat& m1, const mat& m2)' " ){
	mat a,b,c;
	SECTION("Test normal conditions"){
		// Example 1:
		//     |2 3|		|9 4|					     |2 3|
		// a = |1 5|,   b = |0 5|, c = concat_ver(a,b) = |1 5|
		//												 |9 4|
		//												 |0 5|
		a = "[2 3;1 5]";
		b = "[9 4;0 5]";
		c = concat_ver(a, b);
		REQUIRE( c.rows() == a.rows() + b.rows() ); REQUIRE( c.cols() == a.cols()  );
		REQUIRE(c(0,0) == 2); REQUIRE(c(0,1) == 3);
		REQUIRE(c(1,0) == 1); REQUIRE(c(1,1) == 5);
		REQUIRE(c(2,0) == 9); REQUIRE(c(2,1) == 4);
		REQUIRE(c(3,0) == 0); REQUIRE(c(3,1) == 5);

	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE( concat_hor(a, b).size() == 0 );
	}
}

TEST_CASE( " Test 'outer_product(const vec& v1, const vec& v2)' function" ){
	vec a, b; mat c;
	SECTION(" Test normal conditions. "){
		//Example 1: Vectors of equal size ==> sqaure matrix
		//														 |-2 4 12|
		// a = [2 -1 4], b = [-1 2 6]: c = outer_product(a, b) = |1 -2 -6|
		//														 |-4 8 24|
		a = "[2 -1 4]"; b = "[-1 2 6]";
		c = outer_product(a, b);
		REQUIRE( c.rows() == a.size() ); REQUIRE( c.cols() == b.size() );
		REQUIRE(c(0,0) == -2); REQUIRE(c(0,1) == 4); REQUIRE(c(0,2) == 12);
		REQUIRE(c(1,0) == 1); REQUIRE(c(1,1) == -2); REQUIRE(c(1,2) == -6);
		REQUIRE(c(2,0) == -4); REQUIRE(c(2,1) == 8); REQUIRE(c(2,2) == 24);

		//Example 2: Vectors of different size ==> non-square matrix
		//														 |-2 4|
		// a = [2 -1 4], b = [-1 2]: c = outer_product(a, b) =   |1 -2|
		//														 |-4 8|
		a = "[2 -1 4]"; b = "[-1 2]";
		c = outer_product(a, b);
		REQUIRE( c.rows() == a.size() ); REQUIRE( c.cols() == b.size() );
		REQUIRE(c(0,0) == -2); REQUIRE(c(0,1) == 4);
		REQUIRE(c(1,0) == 1); REQUIRE(c(1,1) == -2);
		REQUIRE(c(2,0) == -4); REQUIRE(c(2,1) == 8);

	}
	SECTION(" Test boundary conditions. "){
		REQUIRE_THROWS( outer_product(a, b) );
	}
}

TEST_CASE( " Test 'mat::determinant(const mat& m)' " ){
	mat m;
	SECTION("Test normal conditions"){
		// Example 1: 4x4 matrix (found online)

		//     |3   2  -1   4|
		//     |2   1   5   7|
		// m = |0   5   2  -6|,  determinant(m) = -418
		//     |-1  2   1   0|

		m = "[3 2 -1 4;2 1 5 7;0 5 2 -6;-1 2 1 0]";
		REQUIRE( determinant(m) == Approx(-418) );
	}
	SECTION("Test boundary conditions."){
		// Example 1: NULL MATRIX
		REQUIRE_THROWS( determinant(m) );
	}
}

TEST_CASE( " Test 'mat::magic_square(int n)' " ){
	imat m;
	SECTION("Test normal conditions"){
		// Example 1: 5x5 magic square
		/*|  9    3  22  16  15 |
          |  2   21  20  14   8 |
          |  25  19  13   7   1 |
          |  18  12   6   5  24 |
          |  11  10   4  23  17 |*/
		int size = 5;
		m = magic_square(size);
		int sum_row_0 = sum(m.get_row(0));
		int i;
		// Verify that the sum of each row and each column
		// as well as the sum of the diagonal elements is the same
		for(i = size; i--;){
			REQUIRE( sum(m.get_row(i)) == sum_row_0 );
			REQUIRE( sum(m.get_col(i)) == sum_row_0 );
		}
		REQUIRE( sum(diag(m)) == sum_row_0 );
	}
	SECTION("Test boundary conditions."){
		REQUIRE_THROWS( magic_square(-1) );
		REQUIRE_THROWS( magic_square(4) );
		REQUIRE_THROWS( magic_square(MAX_ACCEPTABLE_VECTOR_SIZE + 1) );
		log_error(" ");
		log_error("===========================================================================================================");
		log_error("================================= UNIT-TEST FOR MATRIX CLASS FINISHED ==================================");
		log_error("===========================================================================================================");
	}
}

} /* namespace algebra */
