/*==============================================================================
 * Name         : vec_test.cpp implements a unit-test for
 *                the 'vec' class of the LinearAlgebra library.
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
===================================================================================*/


// The catch.hpp file was taken from: https://raw.githubusercontent.com/philsquared/Catch/master/single_include/catch.hpp

// catch.hpp tutorial: https://github.com/philsquared/Catch/blob/master/docs/tutorial.md

// Just for your information the REQUIRE_THROWS(expr) functions checks if an exception is thrown
// when one is expected. The UNIT TEST is run for each SECTION separately. It is guaranteed that
// all the thrown exceptions are caught by the current UNIT TEST and therefore there are no memory leaks
// since all the destructors of the instantiated objects are called properly. Take a look at this
// piece of code from the 'catch.hpp' file:
// #define REQUIRE_THROWS( expr ) INTERNAL_CATCH_THROWS( "REQUIRE_THROWS", Catch::ResultDisposition::Normal, "", expr )

#include "../include/catch.hpp"
#include "../../../include/base.h"

namespace algebra {

// *********************** TEST MEMBER FUNCTIONS *************************
// ***********************************************************************

TEST_CASE( " Test vec constructor." ){

	SECTION(" Test boundary conditions."){
		log_error(" ");
		log_error("===========================================================================================================");
		log_error("================================= NEW UNIT-TEST FOR VECTOR CLASS STARTED ==================================");
		log_error("===========================================================================================================");
		log_error(" ");

		INFO("Unit test failed in algebra::vec(size_t n)");  // Only appears on a FAIL
		REQUIRE_THROWS(vec(MAX_ACCEPTABLE_VECTOR_SIZE + 1));
	}
	SECTION(" Test user-defined constructor."){
		vec v(2);
		INFO("Unit test failed in algebra::vec(size_t n)");  // Only appears on a FAIL
		REQUIRE(v.get(0) == 0);
		REQUIRE(v.get(1) == 0);
		REQUIRE(v.size() == 2);
	}
	SECTION(" Test default constructor."){
		vec v;
		INFO("Unit test failed in algebra::vec()");  // Only appears on a FAIL
		REQUIRE(v.size() == 0);
	}
}

TEST_CASE( " Test 'vec::set(size_t i, float k)' function." ){
	vec v(2);
	SECTION(" Test for normal index"){
		v.set(0, 3.51001);
		REQUIRE(v.get(0) == Approx(3.51));
	}
	SECTION(" Test for index out of bounds."){
		REQUIRE_THROWS(v.set(2, 3.4));
		REQUIRE_THROWS(v.set(-1, -5.7));
	}
}

TEST_CASE( " Test 'vec::get(size_t r) const' function." ){
	vec v(2);
	v.set(0, -2.54); v.set(1, 3.51001);
	SECTION(" Test for normal index"){
		REQUIRE(v.get(0) == Approx(-2.54));
		REQUIRE(v.get(1) == Approx(3.51));
	}
	SECTION(" Test for index out of bounds."){
		// The .at(n) throws an std::out_of_range exception if n is an invalid index
		REQUIRE_THROWS(v.get(2));
		REQUIRE_THROWS(v.get(-1));
	}
}

TEST_CASE( " Test 'vec::get(size_t i, size_t j) const' function." ){
	// v1 = [1 4 5 7]; ==> v2 = v1.get(1,2) = [4 5]
	vec v1(4);
	v1.set(0, 1); v1.set(1, 4); v1.set(2, 5); v1.set(3, 7);
	SECTION(" Test for normal index"){
		vec v2 = v1.get(1,2);
		REQUIRE(v2.get(0) == Approx(4));
		REQUIRE(v2.get(1) == Approx(5));
	}
	SECTION(" Test for normal index for complex vector"){
		cvec v3(2); v3(0) = 6. + 2i; v3(1) = 3. - 1i;
		REQUIRE(v3.get(0).real() == 6); REQUIRE(v3.get(0).imag() == 2);
		REQUIRE(v3.get(1).real() == 3); REQUIRE(v3.get(1).imag() == -1);
	}
	SECTION(" Test for boundary conditions."){
		REQUIRE_THROWS(v1.get(1, -3));     // case j < 0
		REQUIRE_THROWS(v1.get(2, 4));      // case j > v1.size()
		REQUIRE_THROWS(v1.get(100, 200));  // case i > v1.size() && j > v1.size()
		REQUIRE_THROWS(v1.get(2, 0));      // case j < i
	}
}

TEST_CASE( " Test 'vec::resize(size_t size)' ." ) {
	vec v(4);
	SECTION(" Decrease size of vector."){
		v.set_size(3);
		REQUIRE(v.size() == 3);
	}
	SECTION(" Increase size of vector."){
		v.set_size(5);
		REQUIRE(v.size() == 5);
	}
	SECTION(" Test for boundary conditions."){
		REQUIRE_THROWS(v.set_size(MAX_ACCEPTABLE_VECTOR_SIZE + 1));
		REQUIRE_THROWS(v.set_size(-1));
	}
}

TEST_CASE( " Test 'vec::set_subvector(size_t start, const vec& v)' ." ){
	vec v(4), v2(2);
	v.set(0,1); v.set(1,2); v.set(2,3); v.set(3,4);
	v2.set(0,44); v2.set(1,55);
	SECTION(" Test normal conditions. "){
		v.set_subvector(1,v2);
		REQUIRE(v.get(1) == 44);
		REQUIRE(v.get(2) == 55);
	}
	SECTION(" Test boundary conditions. "){
		REQUIRE_THROWS(v.set_subvector(3,v2));
		REQUIRE_THROWS(v.set_subvector(4,v2));
		REQUIRE_THROWS(v.set_subvector(-1,v2));
	}
}

TEST_CASE( " Test 'vec::zeros()' " ){
	SECTION("Test for double vector"){
		vec v(2);
		v.set(0,44); v.set(1,55);
		v.zeros();
		REQUIRE(v.get(0) == 0);
		REQUIRE(v.get(1) == 0);
	}
	SECTION("Test for complex vector"){
		cvec v(2);
		v.set(0, 6. + 2i); v.set(1, 3. - 0.5i);
		v.zeros();
		REQUIRE(v.get(0).real() == 0); REQUIRE(v.get(0).imag() == 0);
		REQUIRE(v.get(1).real() == 0); REQUIRE(v.get(1).imag() == 0);
	}
}

TEST_CASE( " Test 'vec::clear()' " ){
	SECTION("Test for double vector"){
		vec v(2);
		v.set(0,44); v.set(1,55);
		v.clear();
		REQUIRE(v.get(0) == 0);
		REQUIRE(v.get(1) == 0);
	}
	SECTION("Test for complex vector"){
		cvec v(2);
		v.set(0, 6. + 2i); v.set(1, 3. - 0.5i);
		v.clear();
		REQUIRE(v.get(0).real() == 0); REQUIRE(v.get(0).imag() == 0);
		REQUIRE(v.get(1).real() == 0); REQUIRE(v.get(1).imag() == 0);
	}
}

TEST_CASE( " Test vec::ones() function" ){
	SECTION("Test for double vector"){
		vec v(2);
		v.set(0,44); v.set(1,55);
		v.ones();
		REQUIRE(v.get(0) == 1);
		REQUIRE(v.get(1) == 1);
	}
	SECTION("Test for complex vector"){
		cvec v(2);
		v.set(0, 6. + 2i); v.set(1, 3. - 0.5i);
		v.ones();
		REQUIRE(v.get(0).real() == 1); REQUIRE(v.get(0).imag() == 0);
		REQUIRE(v.get(1).real() == 1); REQUIRE(v.get(1).imag() == 0);
	}
}

TEST_CASE( " Test vec::add(const vec& v1) function" ){
	vec a(2);
	SECTION(" Test normal conditions. "){
		// c = a + b = [2 3] + [1 2.1] = [3 5.1]
		vec b(2);
		a.set(0,2); a.set(1,3);
		b.set(0,1); b.set(1,2.1);
		vec c = a.add(b);
		REQUIRE(c.size() == 2);
		REQUIRE(c.get(0) == 3);
		REQUIRE(c.get(1) == Approx(5.1));
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// v3 = v1 + v2 = [6+2i  3-0.5i] + [3-1i  -3-2.5i] = [9+1i  0-3i]
		cvec v1(2), v2(2), v3;
		v1(0) = 6. + 2i; v1(1) = 3. -0.5i;
		v2(0) = 3. - 1i; v2(1) = -3. - 2.5i;
		v3 = v1.add(v2);
		REQUIRE(v3.size() == 2);
		REQUIRE(v3.get(0).real() == 9); REQUIRE(v3.get(0).imag() == 1);
		REQUIRE(v3.get(1).real() == 0); REQUIRE(v3.get(1).imag() == -3);
	}
	SECTION(" Test boundary conditions. "){
		vec b(3);
		a.set(0,2); a.set(1,3);
		REQUIRE_THROWS(a.add(b));
		vec c(1);
		REQUIRE_THROWS(a.add(c));
	}
}

TEST_CASE( " Test vec::sub(const vec& v1) function" ){
	vec a(2);
	SECTION(" Test normal conditions. "){
		// c = a + b = [2 3] - [5 1.2] = [-3 1.8]
		vec b(2);
		a.set(0, 2); a.set(1, 3);
		b.set(0, 5); b.set(1, 1.2);
		vec c = a.sub(b);
		REQUIRE(c.size() == 2);
		REQUIRE(c.get(0) == -3);
		REQUIRE(c.get(1) == Approx(1.8));
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// v3 = v1 - v2 = [6+2i  3-0.5i] - [3-1i  -3-2.5i] = [3+3i  6+2i]
		cvec v1(2), v2(2), v3;
		v1(0) = 6. + 2i; v1(1) = 3. -0.5i;
		v2(0) = 3. - 1i; v2(1) = -3. - 2.5i;
		v3 = v1.sub(v2);
		REQUIRE(v3.size() == 2);
		REQUIRE(v3.get(0).real() == 3); REQUIRE(v3.get(0).imag() == 3);
		REQUIRE(v3.get(1).real() == 6); REQUIRE(v3.get(1).imag() == 2);
	}
	SECTION(" Test boundary conditions. "){
		vec b(3);
		a.set(0,2); a.set(1,3);
		REQUIRE_THROWS(a.sub(b));
		vec c(1);
		REQUIRE_THROWS(a.sub(c));
	}
}

TEST_CASE( " Test vec::dot(const vec& v1) function" ){
	vec a(3), b(3);
	SECTION(" Test normal conditions. "){
		// c = a * b = [1 2 3] * [1 -2 5] = 1*1 + 2*(-2) + 3*5 = 12
		a.set(0,1); a.set(1,2); a.set(2,3);
		b.set(0,1); b.set(1,-2); b.set(2,5);
		REQUIRE(a.dot(b) == 12);
		REQUIRE(b.dot(a) == 12);
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// v3 = v1 * v2 = [6+2i  3-0.5i] * [3-1i  -3-2.5i] = 9.75-6i
		cvec v1(2), v2(2), v3;
		v1(0) = 6. + 2i; v1(1) = 3. -0.5i;
		v2(0) = 3. - 1i; v2(1) = -3. - 2.5i;
		REQUIRE(v1.dot(v2).real() == 9.75);
		REQUIRE(v1.dot(v2).imag() == -6);
	}
	SECTION(" Test boundary conditions. "){
		vec c(4);
		REQUIRE_THROWS(a.dot(c));
		vec d(1);
		REQUIRE_THROWS(a.dot(d));
	}
}

TEST_CASE( " Test vec::cross(const vec& v1) function" ){
	SECTION(" Test normal conditions. "){
		vec a(3), b(3), c1, c2;
		// c1 = a x b = [3 -3 1] x [4 9 2] = [-15 -2 39]
		// c2 = b x a = [4 9 2] x [3 -3 1] = [15 2 -39] = -c1
		a.set(0, 3); a.set(1, -3); a.set(2, 1);
		b.set(0, 4); b.set(1, 9); b.set(2, 2);
		c1 = a.cross(b);
		c2 = b.cross(a);
		REQUIRE(c1.get(0) == -15);
		REQUIRE(c1.get(1) == -2);
		REQUIRE(c1.get(2) == 39);
		REQUIRE(c2.get(0) == -c1.get(0));
		REQUIRE(c2.get(1) == -c1.get(1));
		REQUIRE(c2.get(2) == -c1.get(2));
	}
	SECTION(" Test boundary conditions. "){
		vec a(4); vec b(4);
		REQUIRE_THROWS(a.cross(b));
		vec d(1); vec e(3);
		REQUIRE_THROWS(d.cross(e));
	}
}

TEST_CASE( " Test vec::sort() function" ){
	vec a(6);
	SECTION(" Test normal conditions. "){
		// a        = [2 5.6 -12.1 77 32 5.6]
		// a.sort() = [-12.1 2 5.6 32 77 5.6]
		a.set(0,2); a.set(1,5.6); a.set(2,-12.1); a.set(3, 77); a.set(4, 32); a.set(5, 5.6);
		a.sort();
		REQUIRE(a.get(0) == Approx(-12.1));
		REQUIRE(a.get(1) == Approx(2));
		REQUIRE(a.get(2) == Approx(5.6));
		REQUIRE(a.get(3) == Approx(5.6));
		REQUIRE(a.get(4) == Approx(32));
		REQUIRE(a.get(5) == Approx(77));
	}
}

TEST_CASE( " Test vec::swap(size_t i, size_t j) function" ){
	vec a, tmp;
	SECTION(" Test normal conditions. "){
		a = rand(6);
		tmp = a;
		a.swap(1, 3);

		REQUIRE(tmp(0) == a(0)); REQUIRE(tmp(1) == a(3)); REQUIRE(tmp(2) == a(2));
		REQUIRE(tmp(3) == a(1)); REQUIRE(tmp(4) == a(4)); REQUIRE(tmp(5) == a(5));
	}
	SECTION(" Test normal conditions for complex vectors. "){
		cvec b, c;
		b = rand_c(3);
		c = b;
		b.swap(1, 2);

		REQUIRE(c(0).real() == b(0).real()); REQUIRE(c(1).real() == b(2).real()); REQUIRE(c(2).real() == b(1).real());
		REQUIRE(c(0).imag() == b(0).imag()); REQUIRE(c(1).imag() == b(2).imag()); REQUIRE(c(2).imag() == b(1).imag());
	}

	SECTION(" Test normal conditions. "){
		a = rand(6);
		REQUIRE_THROWS( a.swap(1, 6) );
		REQUIRE_THROWS( a.swap(-1, 3) );
	}
}


// *********************** TEST OVERLOADED OPERATORS *************************
// ***************************************************************************

TEST_CASE( " Test vec::overload=(const char* a) function" ){
	vec a;
	SECTION(" Test normal conditions. "){
		// a = [2 5.6 -12.1 77 32]
		a = "[2 5.6 -12.1 77 32]";
		REQUIRE(a.get(0) == Approx(2));
		REQUIRE(a.get(1) == Approx(5.6));
		REQUIRE(a.get(2) == Approx(-12.1));
		REQUIRE(a.get(3) == Approx(77));
		REQUIRE(a.get(4) == Approx(32));

		a = "[1 2 3]";
		REQUIRE(a.size() == 3);
		a = "[1 2]";
		REQUIRE(a.size() == 2);
		a = "[1 2 3 4 5]";
		REQUIRE(a.size() == 5);
		a = "[]";
		REQUIRE(a.size() == 0);
	}
	SECTION(" Test boundary conditions. "){
		// a = []
		a = "[]";
		REQUIRE(a.size() == 0);
		a = "[ ]";
		REQUIRE(a.size() == 0);
		a = "";
		REQUIRE(a.size() == 0);
		a = " ";
		REQUIRE(a.size() == 0);
		// It should throw exception if a letter occurs
		REQUIRE_THROWS(a = "[2 3.9 -1k]");
	}
}

TEST_CASE( " Test vec::overload+(const vec& a) function" ){
	vec a ,b;
	SECTION(" Test normal conditions. "){
		// a = [3 -5 9], b = [4 3 8] ==> c = a + b = b + a = [7 -2 17]
		a = "[3 -5 9]", b = "[4 3 8]";
		vec c = a + b;
		vec d = b + a; // Test commutative property.
		REQUIRE(c.get(0) == 7);
		REQUIRE(c.get(1) == -2);
		REQUIRE(c.get(2) == 17);
		REQUIRE(c.get(0) == d.get(0));
		REQUIRE(c.get(1) == d.get(1));
		REQUIRE(c.get(2) == d.get(2));
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// v3 = v1 + v2 = [6+2i  3-0.5i] + [3-1i  -3-2.5i] = [9+1i  0-3i]
		cvec v1(2), v2(2), v3;
		v1(0) = 6. + 2i; v1(1) = 3. -0.5i;
		v2(0) = 3. - 1i; v2(1) = -3. - 2.5i;
		v3 = v1 + v2;
		REQUIRE(v3.size() == 2);
		REQUIRE(v3.get(0).real() == 9); REQUIRE(v3.get(0).imag() == 1);
		REQUIRE(v3.get(1).real() == 0); REQUIRE(v3.get(1).imag() == -3);
	}
	SECTION(" Test boundary conditions. "){
		REQUIRE((a + b).size() == 0);
		// a = [3 -5 9], b = [4 3]
		a = "[3 -5 9]", b = "[4 3]";
		REQUIRE_THROWS(a + b);
	}
}

TEST_CASE( " Test vec::overload+(T t) function" ){
	vec a;
	float p;
	SECTION(" Test normal conditions. "){
		// a = [3 -5 9], p = 3.5 ==> c = a + p = [6.5 -1.5 12.5]
		a = "[3 -5 9]", p = 3.5;
		vec c = a + p;
		REQUIRE(c.get(0) == Approx(6.5));
		REQUIRE(c.get(1) == Approx(-1.5));
		REQUIRE(c.get(2) == Approx(12.5));
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// v2 = v1 + t = [6+2i  3-0.5i] + 1+2i = [7+4i  4+1.5i]
		cvec v1(2), v2;
		std::complex<double> t = 1. + 2i;
		v1(0) = 6. + 2i; v1(1) = 3. -0.5i;
		v2 = v1 + t;
		REQUIRE(v2.size() == 2);
		REQUIRE(v2.get(0).real() == 7); REQUIRE(v2.get(0).imag() == 4);
		REQUIRE(v2.get(1).real() == 4); REQUIRE(v2.get(1).imag() == 1.5);
	}
	SECTION(" Test boundary conditions. "){
		p = 3.5;
		REQUIRE_THROWS(a + p);
	}
}

TEST_CASE( " Test vec::overload-(const vec& a) function" ){
	vec a ,b;
	SECTION(" Test normal conditions. "){
		// a = [3 -5 9], b = [4 3 8] ==> c = a - b = [-1 -8 1]
		a = "[3 -5 9]", b = "[4 3 8]";
		vec c = a - b;
		REQUIRE(c.get(0) == -1);
		REQUIRE(c.get(1) == -8);
		REQUIRE(c.get(2) == 1);
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// v3 = v1 - v2 = [6+2i  3-0.5i] - [3-1i  -3-2.5i] = [3+3i  6+2i]
		cvec v1(2), v2(2), v3;
		v1(0) = 6. + 2i; v1(1) = 3. -0.5i;
		v2(0) = 3. - 1i; v2(1) = -3. - 2.5i;
		v3 = v1 - v2;
		REQUIRE(v3.size() == 2);
		REQUIRE(v3.get(0).real() == 3); REQUIRE(v3.get(0).imag() == 3);
		REQUIRE(v3.get(1).real() == 6); REQUIRE(v3.get(1).imag() == 2);
	}
	SECTION(" Test boundary conditions. "){
		REQUIRE((a + b).size() == 0);
		// a = [3 -5 9], b = [4 3]
		a = "[3 -5 9]", b = "[4 3]";
		REQUIRE_THROWS(a - b);
	}
}

TEST_CASE( " Test vec::overload-(T t) function" ){
	vec a;
	float p;
	SECTION(" Test normal conditions. "){
		// a = [3 -5 9], p = 3.5 ==> c = a - p = [-0.5 -8.5 5.5]
		a = "[3 -5 9]", p = 3.5;
		vec c = a - p;
		REQUIRE(c.get(0) == Approx(-0.5));
		REQUIRE(c.get(1) == Approx(-8.5));
		REQUIRE(c.get(2) == Approx(5.5));
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// v2 = v1 - t = [6+2i  3-0.5i] - 1+2i = [5+0i  2-2.5i]
		cvec v1(2), v2;
		std::complex<double> t = 1. + 2i;
		v1(0) = 6. + 2i; v1(1) = 3. -0.5i;
		v2 = v1 - t;
		REQUIRE(v2.size() == 2);
		REQUIRE(v2.get(0).real() == 5); REQUIRE(v2.get(0).imag() == 0);
		REQUIRE(v2.get(1).real() == 2); REQUIRE(v2.get(1).imag() == -2.5);
	}
	SECTION(" Test boundary conditions. "){
		p = 3.5;
		REQUIRE_THROWS(a - p);
	}
}


TEST_CASE( " Test vec::overload*(const vec& a) function" ){
	vec a ,b;
	SECTION(" Test normal conditions. "){
		// c = a + b = [1 2 3] * [1 -2 5] = 1*1 + 2*(-2) + 3*5 = 12
		a = "[1 2 3]", b = "[1 -2 5]";
		float c = a*b;
		float d = b*a;
		REQUIRE(c == 12);
		REQUIRE(d == c); // Commutative property
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// v3 = v1 * v2 = [6+2i  3-0.5i] * [3-1i  -3-2.5i] = 9.75-6i
		cvec v1(2), v2(2);
		std::complex<double> tmp;
		v1(0) = 6. + 2i; v1(1) = 3. -0.5i;
		v2(0) = 3. - 1i; v2(1) = -3. - 2.5i;
		tmp = v1 * v2;
		REQUIRE(tmp.real() == 9.75);
		REQUIRE(tmp.imag() == -6);
	}
	SECTION(" Test boundary conditions. "){
		REQUIRE((a + b).size() == 0);
		// a = [3 -5 9], b = [4 3]
		a = "[3 -5 9]", b = "[4 3]";
		REQUIRE_THROWS(a*b);
	}
}

TEST_CASE( " Test vec::overload*(T t) function" ){
	vec a,c;
	float p;
	SECTION(" Test normal conditions. "){
		// Example 1:
		// a = [3 -5 9], p = 1 ==> c = a*p = a = [3 -5 9]
		a = "[3 -5 9]", p = 1;
		c = a*p;
		REQUIRE(c.get(0) == a.get(0));
		REQUIRE(c.get(1) == a.get(1));
		REQUIRE(c.get(2) == a.get(2));

		// Example 2:
		// a = [3 -5 9 7], p = 0.5 ==> c = a*p = [1.5 -2.5 4.5 3.5]
		a = "[3 -5 9 7]";
		p = 0.5;
		c = a*p;
		REQUIRE(c.get(0) == Approx(1.5));
		REQUIRE(c.get(1) == Approx(-2.5));
		REQUIRE(c.get(2) == Approx(4.5));
		REQUIRE(c.get(3) == Approx(3.5));
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// v2 = v1 * t = [6+2i  3-0.5i] * 1+2i = [2+16i  4+5.5i]
		cvec v1(2), v2;
		std::complex<double> t = 1. + 2i;
		v1(0) = 6. + 2i; v1(1) = 3. -0.5i;
		v2 = v1 * t;
		REQUIRE(v2.size() == 2);
		REQUIRE(v2.get(0).real() == 2); REQUIRE(v2.get(0).imag() == 14);
		REQUIRE(v2.get(1).real() == 4); REQUIRE(v2.get(1).imag() == 5.5);
	}
	SECTION(" Test vec::boundary conditions. "){
		// Try to multiply with null vector.
		p = 3.5;
		REQUIRE( (a*p).size() == 0 );
	}
}

TEST_CASE( " Test vec::overload/(float t) function" ){
	vec a;
	float p;
	SECTION(" Test normal conditions. "){
		// a = [3 -5 9], p = 1 ==> c = a/p = a = [3 -5 9]
		// a = [3 -5 9], p = 0.5 ==> c = a/p = [6 -10 18]
		a = "[3 -5 9]", p = 1;
		vec c = a/p;
		REQUIRE(c.get(0) == a.get(0));
		REQUIRE(c.get(1) == a.get(1));
		REQUIRE(c.get(2) == a.get(2));

		p = 0.5;
		c = a/p;
		REQUIRE(c.get(0) == Approx(6));
		REQUIRE(c.get(1) == Approx(-10));
		REQUIRE(c.get(2) == Approx(18));
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// v2 = v1 / t = [6+2i  3-0.5i] / 1+2i
		cvec v1(2), v2, test(2);
		std::complex<double> t = 1. + 2i;
		v1(0) = 6. + 2i; v1(1) = 3. -0.5i;

		// Polar coordinates according to Euler's formula
		// if z = a + ib ==> z = |z|( cos(phi)+isin(phi) )
		// where |z| = abs(z) = sqrt(a^2+b^2)
		// and phi = arg(z) = arctan(b/a)

		double magn_0, magn_1, magn_t, arg_0, arg_1, arg_t;
		magn_0 = std::abs(v1(0)); arg_0 = atan2(v1(0).imag(), v1(0).real());
		magn_1 = std::abs(v1(1)); arg_1 = atan2(v1(1).imag(), v1(1).real());
		magn_t = std::abs(t); arg_t = atan2(t.imag(), t.real());
		test(0).real( (magn_0/magn_t)*(cos(arg_0-arg_t)) );
		test(0).imag( (magn_0/magn_t)*(sin(arg_0-arg_t)) );
		test(1).real( (magn_1/magn_t)*(cos(arg_1-arg_t)) );
		test(1).imag( (magn_1/magn_t)*(sin(arg_1-arg_t)) );

		v2 = v1 / t;
		REQUIRE(v2.size() == 2);
		REQUIRE(v2(0).real() == Approx( test(0).real()) ); REQUIRE(v2(0).imag() == Approx( test(0).imag()) );
		REQUIRE(v2(1).real() == Approx( test(1).real()) ); REQUIRE(v2(1).imag() == Approx( test(1).imag()) );
	}
	SECTION(" Test boundary conditions. "){
		p = 3.5;
		// Try to access null vector.
		REQUIRE( (a/p).size() == 0);
		a = "[1 2 3]";
		p = 0;
		REQUIRE_THROWS(a/p);
	}
}


TEST_CASE( " Test vec::overload()(const size_t) function" ){
	vec a;
	float p;
	SECTION(" Test normal conditions. "){
		// a = [3 -5 9]
		a = "[3 -5 9]";
		REQUIRE(a(0) == a.get(0));
		REQUIRE(a(1) == a.get(1));
		REQUIRE(a(2) == a.get(2));
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// b = [3+i  -1+2i]
		cvec b(2);
		b(0).real(3); b(0).imag(1);
		b(1).real(-1); b(0).imag(2);

		REQUIRE(b(0).real() == b.get(0).real()); REQUIRE(b(0).imag() == b.get(0).imag());
		REQUIRE(b(1).real() == b.get(1).real()); REQUIRE(b(1).imag() == b.get(1).imag());
	}
	SECTION(" Test boundary conditions. "){
		// Try to access null vector.
		REQUIRE_THROWS(a(1));
		a = "[1 2 3]";
		REQUIRE_THROWS(a(3));
	}
}

TEST_CASE( " Test vec::overload[](const size_t) function" ){
	vec a;
	float p;
	SECTION(" Test normal conditions. "){
		// a = [3 -5 9]
		a = "[3 -5 9]";
		REQUIRE(a[0] == a.get(0));
		REQUIRE(a[1] == a.get(1));
		REQUIRE(a[2] == a.get(2));
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// b = [3+i  -1+2i]
		cvec b(2);
		b(0).real(3); b(0).imag(1);
		b(1).real(-1); b(0).imag(2);

		REQUIRE(b[0].real() == b.get(0).real()); REQUIRE(b[0].imag() == b.get(0).imag());
		REQUIRE(b[1].real() == b.get(1).real()); REQUIRE(b[1].imag() == b.get(1).imag());
	}
	SECTION(" Test boundary conditions. "){
		// Try to access null vector.
		REQUIRE_THROWS(a[1]);
		a = "[1 2 3]";
		REQUIRE_THROWS(a[3]);
	}
}

// *************************** TEST INLINE FUNCTIONS *************************
// ***************************************************************************

TEST_CASE( " Test algebra::zeros(size_t n) function" ){
	SECTION(" Test normal conditions. "){
		vec a = zeros(2);
		REQUIRE(a.size() == 2);
		REQUIRE(a[0] == 0);
		REQUIRE(a[1] == 0);
	}
	SECTION(" Test boundary conditions. "){
		REQUIRE_THROWS(zeros(MAX_ACCEPTABLE_VECTOR_SIZE + 1));
		REQUIRE_THROWS(zeros(-1));
	}
}

TEST_CASE( " Test algebra::ones(size_t n) function" ){
	SECTION(" Test normal conditions. "){
		vec a = ones(2);
		REQUIRE(a.size() == 2);
		REQUIRE(a[0] == 1);
		REQUIRE(a[1] == 1);
	}
	SECTION(" Test boundary conditions. "){
		REQUIRE_THROWS(ones(MAX_ACCEPTABLE_VECTOR_SIZE + 1));
		REQUIRE_THROWS(ones(-1));
	}
}

TEST_CASE( " Test algebra::find_non_zero(const vec& v) function" ){
	vec v1, v2;
	SECTION(" Test normal conditions. "){
		v1 = "[3 4 0 0 1 0]";
		v2 = find_non_zero(v1);
		size_t i, size = v1.size();
		for(i = size; i--;){
			if(v1(i) != 0){
				REQUIRE(v2(i) == 1);
			}else{
				REQUIRE(v2(i) == 0);
			}
		}
	}
}

TEST_CASE( " Test algebra::find_zero(const vec& v) function" ){
	vec v1, v2;
	SECTION(" Test normal conditions. "){
		v1 = "[3 4 0 0 1 0]";
		v2 = find_zero(v1);
		size_t i, size = v1.size();
		for(i = size; i--;){
			if(v1(i) == 0){
				REQUIRE(v2(i) == 1);
			}else{
				REQUIRE(v2(i) == 0);
			}
		}
	}
}

TEST_CASE( " Test algebra::dot(const vec& v1, const vec& v2) function" ){
	vec a, b;
	SECTION(" Test normal conditions. "){
		// c = a + b = [1 2 3] * [1 -2 5] = 1*1 + 2*(-2) + 3*5 = 12
		a = "[1 2 3]"; b = "[1 -2 5]";
		float c = dot(a,b);
		float d = dot(b,a);
		REQUIRE(c == 12);
		REQUIRE(d == c);
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// v3 = v1 * v2 = [6+2i  3-0.5i] * [3-1i  -3-2.5i] = 9.75-6i
		cvec v1(2), v2(2), v3;
		v1(0) = 6. + 2i; v1(1) = 3. -0.5i;
		v2(0) = 3. - 1i; v2(1) = -3. - 2.5i;
		REQUIRE(dot(v1, v2).real() == 9.75);
		REQUIRE(dot(v1, v2).imag() == -6);
	}
	SECTION(" Test boundary conditions. "){
		a = "[1 2 3]"; b = "[1 -2]";
		REQUIRE_THROWS(dot(a,b));
		a.set_size(0); b = "[1 -2]";
		REQUIRE(dot(a,b) == 0); // dot product with at least one null vector is 0.
	}
}


TEST_CASE( " Test algebra::mean(const Vec<T>& v1) function" ){
	vec a;
	SECTION(" Test normal conditions. "){
		// a = [2 1 4 3] ==> mean(a) = 2.5
		a = "[2 1 4 3]";
		float c = mean(a);
		REQUIRE(c == Approx(2.5));
	}
	SECTION(" Test normal conditions for complex vectors. "){
		cvec b = rand_c(4);
		std::complex<double> c = mean(b);
		double real = 0 , imaginary = 0;
		for(size_t i = 0; i < b.size(); i++)
		{
			real += b(i).real();
			imaginary += b(i).imag();
		}
		real /= b.size();
		imaginary /= b.size();

		REQUIRE(c.real() == Approx(real));
		REQUIRE(c.imag() == Approx(imaginary));
	}
	SECTION(" Test boundary conditions. "){
		REQUIRE_THROWS(mean(a));
	}
}

TEST_CASE( " Test algebra::min(const vec& v1) function" ){
	vec a;
	SECTION(" Test normal conditions. "){
		// a = [2 -1 4 -3] ==> min(a) = -3
		a = "[2 -1 4 -3]";
		float c = min(a);
		REQUIRE(c == Approx(-3));
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// b = [3-2i  -3+2i  4+i]

		// b(0) = 3-2i = |b(0)|*exp(i*phi) = sqrt(13)*exp(i*(-0.588))
		// b(1) =-3+2i = |b(1)|*exp(i*phi) = sqrt(13)*exp(i*(-0.588+pi))
		//                                 = sqrt(13)*exp(i*2.553)
		// b(2) = 4+i  = |b(2)|*exp(i*phi) = sqrt(15)*exp(i*0.2449)

		// So, it should be min(b) = b(0)

		cvec b(3);
		b(0).real(3); b(0).imag(-2);
		b(1).real(-3); b(1).imag(2);
		b(2).real(4); b(2).imag(1);

		// Remember that for a complex number z = x + yi

		//             | arctan(y/x)      , x > 0
		// atan2(y,x)= | arctan(y/x) + pi , x < 0, y >= 0
		//             | arctan(y/x) - pi , x < 0, y < 0

		// see: https://en.wikipedia.org/wiki/Atan2

		REQUIRE( min(b) == b(0) );
	}
	SECTION(" Test boundary conditions. "){
		REQUIRE_THROWS(min(a));
	}
}

TEST_CASE( " Test algebra::min(const vec& v, size_t &index) function" ){
	vec a;
	SECTION(" Test normal conditions. "){
		// a = [2 -1 4 -3] ==> min(a, ind) = -3 and ind = 3
		a = "[2 -1 4 -3]";
		float minimum_value; size_t minimum_index;
		minimum_value = min(a, minimum_index);
		REQUIRE(minimum_value == Approx(-3));
		REQUIRE(minimum_index == 3);
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// b = [3-2i  -3+2i  4+i]

		// b(0) = 3-2i = |b(0)|*exp(i*phi) = sqrt(13)*exp(i*(-0.588))
		// b(1) =-3+2i = |b(1)|*exp(i*phi) = sqrt(13)*exp(i*(-0.588+pi))
		//                                 = sqrt(13)*exp(i*2.553)
		// b(2) = 4+i  = |b(2)|*exp(i*phi) = sqrt(15)*exp(i*0.2449)

		// So, it should be min(b, ind) = b(0) and ind = 0

		cvec b(3);
		b(0).real(3); b(0).imag(-2);
		b(1).real(-3); b(1).imag(2);
		b(2).real(4); b(2).imag(1);

		// Remember that for a complex number z = x + yi

		//             | arctan(y/x)      , x > 0
		// atan2(y,x)= | arctan(y/x) + pi , x < 0, y >= 0
		//             | arctan(y/x) - pi , x < 0, y < 0

		// see: https://en.wikipedia.org/wiki/Atan2

		std::complex<double> minimum_value; size_t minimum_index;
		minimum_value = min(b, minimum_index);

		REQUIRE( minimum_value == b(0) );
		REQUIRE( minimum_index == 0 );
	}
	SECTION(" Test boundary conditions. "){
		size_t minimum_index;
		REQUIRE_THROWS(min(a, minimum_index));
	}
}


TEST_CASE( " Test algebra::max(const vec& v1) function" ){
	vec a;
	SECTION(" Test normal conditions. "){
		// a = [2 -1 4 -3] ==> max(a) = 4
		a = "[2 -1 4 -3]";
		float c = max(a);
		REQUIRE(c == Approx(4));
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// b = [3-2i  -3+2i  4+i]

		// b(0) = 3-2i = |b(0)|*exp(i*phi) = sqrt(13)*exp(i*(-0.588))
		// b(1) =-3+2i = |b(1)|*exp(i*phi) = sqrt(13)*exp(i*(-0.588+pi))
		//                                 = sqrt(13)*exp(i*2.553)
		// b(2) = 4+i  = |b(2)|*exp(i*phi) = sqrt(15)*exp(i*0.2449)

		// So, it should be max(b) = b(2)

		cvec b(3);
		b(0).real(3); b(0).imag(-2);
		b(1).real(-3); b(1).imag(2);
		b(2).real(4); b(2).imag(1);

		// Remember that for a complex number z = x + yi

		//             | arctan(y/x)      , x > 0
		// atan2(y,x)= | arctan(y/x) + pi , x < 0, y >= 0
		//             | arctan(y/x) - pi , x < 0, y < 0

		// see: https://en.wikipedia.org/wiki/Atan2

		REQUIRE( max(b) == b(2) );
	}
	SECTION(" Test boundary conditions. "){
		REQUIRE_THROWS(max(a));
	}
}

TEST_CASE( " Test algebra::max(const vec& v, size_t &index) function" ){
	vec a;
	SECTION(" Test normal conditions. "){
		// a = [2 -1 4 -3] ==> max(a, ind) = 4 and ind = 2
		a = "[2 -1 4 -3]";
		float maximum_value; size_t maximum_index;
		maximum_value = max(a, maximum_index);
		REQUIRE(maximum_value == Approx(4));
		REQUIRE(maximum_index == 2);
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// b = [3-2i  -3+2i  4+i]

		// b(0) = 3-2i = |b(0)|*exp(i*phi) = sqrt(13)*exp(i*(-0.588))
		// b(1) =-3+2i = |b(1)|*exp(i*phi) = sqrt(13)*exp(i*(-0.588+pi))
		//                                 = sqrt(13)*exp(i*2.553)
		// b(2) = 4+i  = |b(2)|*exp(i*phi) = sqrt(15)*exp(i*0.2449)

		// So, it should be max(b, ind) = b(2) and ind = 2

		cvec b(3);
		b(0).real(3); b(0).imag(-2);
		b(1).real(-3); b(1).imag(2);
		b(2).real(4); b(2).imag(1);

		// Remember that for a complex number z = x + yi

		//             | arctan(y/x)      , x > 0
		// atan2(y,x)= | arctan(y/x) + pi , x < 0, y >= 0
		//             | arctan(y/x) - pi , x < 0, y < 0

		// see: https://en.wikipedia.org/wiki/Atan2

		std::complex<double> maximum_value; size_t maximum_index;
		maximum_value = max(b, maximum_index);

		REQUIRE( maximum_value == b(2) );
		REQUIRE( maximum_index == 2 );
	}
	SECTION(" Test boundary conditions. "){
		size_t maximum_index;
		REQUIRE_THROWS(max(a, maximum_index));
	}
}

TEST_CASE( " Test algebra::cross(const vec& v1, const vec& v2) function" ){
	SECTION(" Test normal conditions. "){
		vec a(3), b(3), c1, c2;
		// c1 = a x b = [3 -3 1] x [4 9 2] = [-15 -2 39]
		// c2 = b x a = [4 9 2] x [3 -3 1] = [15 2 -39] = -c1
		a = "[3 -3 1]"; b = "[4 9 2]";
		c1 = cross(a, b); c2 = cross(b, a);
		REQUIRE(c1(0) == -15);
		REQUIRE(c1(1) == -2);
		REQUIRE(c1(2) == 39);
		REQUIRE(c2(0) == -c1(0));
		REQUIRE(c2(1) == -c1(1));
		REQUIRE(c2(2) == -c1(2));
	}
	SECTION(" Test boundary conditions. "){
		vec a(4); vec b(4);
		REQUIRE_THROWS(cross(a, b));
		vec d(1); vec e(3);
		REQUIRE_THROWS(cross(d, e));
	}
}

TEST_CASE( " Test algebra::concat(const Vec<T>& v1, T t) function" ){
	vec a; double t;
	SECTION(" Test normal conditions. "){
		// a = [2 -1 4 -3], t = 99: a = concat(a, t) = [2 -1 4 -3 99]
		a = "[2 -1 4 -3]"; t = 99;
		a = concat(a, t);
		REQUIRE(a.size() == 5);
		REQUIRE(a(0) == 2);
		REQUIRE(a(1) == -1);
		REQUIRE(a(2) == 4);
		REQUIRE(a(3) == -3);
		REQUIRE(a(4) == 99);
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// b = [-2+5i  3-4i], k = 1+88i, b = concat(b, k) = [-2+5i  3-4i  1+88i]
		cvec b(2); std::complex<double> k = 1. + 88i;
		b(0).real(-2); b(0).imag(5);
		b(1).real(3); b(1).imag(-4);
		b = concat(b, k);
		REQUIRE( b.size() == 3 );
		REQUIRE( b(0).real() == -2 ); REQUIRE( b(0).imag() == 5 );
		REQUIRE( b(1).real() == 3 ); REQUIRE( b(1).imag() == -4 );
		REQUIRE( b(2).real() == 1 ); REQUIRE( b(2).imag() == 88 );

	}
	SECTION(" Test boundary conditions. "){
		t = 3;
		a = concat(a,t); // From null vector to vector of size 1
		REQUIRE(a.size() == 1);
		REQUIRE(a(0) == t);
		a.set_size(MAX_ACCEPTABLE_VECTOR_SIZE);
		REQUIRE_THROWS( a = concat(a,t) );
	}
}

TEST_CASE( " Test algebra::concat(float t, const vec& v) function" ){
	vec a; double t;
	SECTION(" Test normal conditions. "){
		// a = [2 -1 4 -3], t = 99: a = concat(t, a) = [99 2 -1 4 -3]
		a = "[2 -1 4 -3]"; t = 99;
		a = concat(t, a);
		REQUIRE(a.size() == 5);
		REQUIRE(a(0) == 99);
		REQUIRE(a(1) == 2);
		REQUIRE(a(2) == -1);
		REQUIRE(a(3) == 4);
		REQUIRE(a(4) == -3);
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// b = [-2+5i  3-4i], k = 1+88i, b = concat(k, b) = [1+88i  -2+5i  3-4i]
		cvec b(2); std::complex<double> k = 1. + 88i;
		b(0).real(-2); b(0).imag(5);
		b(1).real(3); b(1).imag(-4);
		b = concat(k, b);
		REQUIRE( b.size() == 3 );
		REQUIRE( b(0).real() == 1 ); REQUIRE( b(0).imag() == 88 );
		REQUIRE( b(1).real() == -2 ); REQUIRE( b(1).imag() == 5 );
		REQUIRE( b(2).real() == 3 ); REQUIRE( b(2).imag() == -4 );

	}
	SECTION(" Test boundary conditions. "){
		a = concat(t, a);
		REQUIRE(a.size() == 1);
		REQUIRE(a(0) == t);
		a.set_size(MAX_ACCEPTABLE_VECTOR_SIZE);
		REQUIRE_THROWS( a = concat(t, a) );
	}
}

TEST_CASE( " Test algebra::concat(const vec& v1, const vec& v2) function" ){
	vec a, b, c;
	SECTION(" Test normal conditions. "){
		// a = [2 -1 4 -3], b = [44 55]: c = concat(a, b) = [2 -1 4 -3 44 55]
		a = "[2 -1 4 -3]"; b = "[44 55]";
		c = concat(a, b);
		REQUIRE(c.size() == a.size() + b.size());
		REQUIRE(c(0) == 2);
		REQUIRE(c(1) == -1);
		REQUIRE(c(2) == 4);
		REQUIRE(c(3) == -3);
		REQUIRE(c(4) == 44);
		REQUIRE(c(5) == 55);
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// b = [-2+5i  3-4i], k = [1+88i  9-12i], b = concat(b, k) = [-2+5i  3-4i  1+88i  9-12i]
		cvec b(2), k(2);
		b(0).real(-2); b(0).imag(5);
		b(1).real(3); b(1).imag(-4);
		k(0).real(1); k(0).imag(88);
		k(1).real(9); k(1).imag(-12);

		b = concat(b, k);

		REQUIRE( b.size() == 4 );
		REQUIRE( b(0).real() == -2 ); REQUIRE( b(0).imag() == 5 );
		REQUIRE( b(1).real() == 3 ); REQUIRE( b(1).imag() == -4 );
		REQUIRE( b(2).real() == 1 ); REQUIRE( b(2).imag() == 88 );
		REQUIRE( b(3).real() == 9 ); REQUIRE( b(3).imag() == -12 );
	}
	SECTION(" Test boundary conditions. "){
		c = concat(a,b);
		REQUIRE(a.size() == 0);
		a.set_size(floor(MAX_ACCEPTABLE_VECTOR_SIZE/2.0) + 2);
		b.set_size(floor(MAX_ACCEPTABLE_VECTOR_SIZE/2.0) + 2);
		REQUIRE_THROWS( c = concat(a, b) );
	}
}

TEST_CASE( " Test algebra::linspace(float from, float to, size_t step) function" ){
	vec a;
	SECTION(" Test normal conditions. "){
		// a = linspace(-3, 10, 2) = [-3 -1 1 3 5 7]
		a = linspace(-3.0, 8.0, 2);
		REQUIRE(a.size() == 6);
		REQUIRE(a(0) == -3);
		REQUIRE(a(1) == -1);
		REQUIRE(a(2) == 1);
		REQUIRE(a(3) == 3);
		REQUIRE(a(4) == 5);
		REQUIRE(a(5) == 7);
	}
	SECTION(" Test boundary conditions. "){
		REQUIRE_THROWS(linspace(3, 2, 1));
		REQUIRE_THROWS(linspace(-5, -10, 1));
		REQUIRE_THROWS(linspace(1, MAX_ACCEPTABLE_VECTOR_SIZE + 10, 1));
	}
}

TEST_CASE( " Test algebra::elem_mult(const vec& v1, const vec& v2) function" ){
	vec a, b, c;
	SECTION(" Test normal conditions. "){
		// a = [2 -1 4 -3], b = [0.5 -2 0.75 -3]: c = elem_mult(a, b) = [1 2 3 9]
		a = "[2 -1 4 -3]"; b = "[0.5 -2 0.75 -3]";
		c = elem_mult(a, b);
		REQUIRE(c.size() == a.size());
		REQUIRE(c(0) == 1);
		REQUIRE(c(1) == 2);
		REQUIRE(c(2) == 3);
		REQUIRE(c(3) == 9);
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// b = [2+i  3-i], c = [2-i 3+i], d = elem_mult(b, c) = [5 10]
		cvec b(2), c(2), d;
		b(0).real(2); b(0).imag(1); b(1).real(3); b(1).imag(-1);
		c(0).real(2); c(0).imag(-1); c(1).real(3); c(1).imag(1);
		d = elem_mult(b, c);
		REQUIRE( d.size() == 2 );
		REQUIRE( d(0).real() == 5 ); REQUIRE( d(0).imag() == 0 );
		REQUIRE( d(1).real() == 10 ); REQUIRE( d(1).imag() == 0 );
	}
	SECTION(" Test boundary conditions. "){
		a = "[2 -1 4 -3]"; b = "[4 5]";
		REQUIRE_THROWS( c = elem_mult(a, b) );
	}
}

TEST_CASE( " Test algebra::sum(const vec& v1) function" ){
	vec a;
	SECTION(" Test normal conditions. "){
		// sum of NULL VECTOR should be 0
		REQUIRE( sum(a) == 0);
		// a = [2 -1 4 -3], k = sum(a) = 2
		a = "[2 -1 4 -3]";
		REQUIRE( sum(a) == 2);
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// b = [2+i  3-i] ==> sum(b) = 5
		cvec b(2);
		b(0).real(2); b(0).imag(1); b(1).real(3); b(1).imag(-1);
		REQUIRE( sum(b).real() == 5 );
		REQUIRE( sum(b).imag() == 0 );
	}
}

TEST_CASE( " Test algebra::cumsum(const vec& v1) function" ){
	vec a, b;
	SECTION(" Test normal conditions. "){
		// a = [2 -1 4 -3], b = cumsum(a) = [2 1 5 2]
		a = "[2 -1 4 -3]";
		b = cumsum(a);
		REQUIRE( b(0) == 2);
		REQUIRE( b(1) == 1);
		REQUIRE( b(2) == 5);
		REQUIRE( b(3) == 2);
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// v = [1 3i 2+4i -6+8i 9-2i] ==> cumsum(v) = [1  1+3i  3+7i  -3+15i  6+13i]
		cvec v(5), v2;
		v(0).real(1); v(0).imag(0);
		v(1).real(0); v(1).imag(3);
		v(2).real(2); v(2).imag(4);
		v(3).real(-6); v(3).imag(8);
		v(4).real(9); v(4).imag(-2);
		v2 = cumsum(v);
		REQUIRE( v2(0).real() == 1 ); REQUIRE( v2(0).imag() == 0 );
		REQUIRE( v2(1).real() == 1 ); REQUIRE( v2(1).imag() == 3 );
		REQUIRE( v2(2).real() == 3 ); REQUIRE( v2(2).imag() == 7 );
		REQUIRE( v2(3).real() == -3 ); REQUIRE( v2(3).imag() == 15 );
		REQUIRE( v2(4).real() == 6 ); REQUIRE( v2(4).imag() == 13 );
	}
	SECTION(" Test boundary conditions. "){
		b = cumsum(a);
		REQUIRE( b.size() == 0 );
		a.set_size(1);
		b = cumsum(a);
		REQUIRE( b.size() == 1 );
	}
}

TEST_CASE( " Test algebra::norm(const vec& v) function" ){
	vec a;
	SECTION(" Test normal conditions. "){
		// a = [2 -1 4 -3], norm(a) = 5.4777225
		a = "[2 -1 4 -3]";
		REQUIRE( norm(a) == Approx(5.477).epsilon(0.001));
	}
	SECTION(" Test normal conditions for complex vectors. "){
		// b = [2+i  3-2i  -1-i], norm(b) = sqrt(20)
		cvec b(3);
		b(0).real(2); b(0).imag(1);
		b(1).real(3); b(1).imag(-2);
		b(2).real(-1); b(2).imag(-1);
		REQUIRE( norm(b) == sqrt(20) );
	}
	SECTION(" Test boundary conditions. "){
		REQUIRE( norm(a) == 0 );
	}
}

TEST_CASE(" Test algebra::conj(const cvec& v) function" ){
	// a = [1+2i  -4-5i], conj(a) = [1-2i  -4+5i]
	cvec a(2), a2;
	a(0).real(1); a(0).imag(2);
	a(1).real(-4); a(1).imag(-5);
	a2 = conj(a);
	REQUIRE( a2(0).real() == a(0).real() ); REQUIRE( a2(0).imag() == -a(0).imag() );
	REQUIRE( a2(1).real() == a(1).real() ); REQUIRE( a2(1).imag() == -a(1).imag() );
}

TEST_CASE( " Test algebra::abs(const vec& v) function" ){
	vec a, b;
	SECTION(" Test normal conditions. "){
		// a = [2 -1 4 -3], b = abs(a) = [2 1 4 3]
		a = "[2 -1 4 -3]";
		b = abs(a);
		REQUIRE( b(0) == 2);
		REQUIRE( b(1) == 1);
		REQUIRE( b(2) == 4);
		REQUIRE( b(3) == 3);
	}
	SECTION(" Test boundary conditions. "){
		b = abs(a);
		REQUIRE( b.size() == 0 );
		log_error(" ");
		log_error("===========================================================================================================");
		log_error("================================= UNIT-TEST FOR VECTOR CLASS FINISHED ==================================");
		log_error("===========================================================================================================");
	}
}

} /* namespace algebra */
