/*============================================================================
 * Name         : vec.h implements basic vector operations
 *                and can be used for linear algebra projects.
 * Version      : 1.0.0, 11 Sep 2017
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
================================================================================*/


#ifndef VEC_H_
#define VEC_H_

#include <stdio.h>      // printf, scanf, NULL
#include <stdlib.h>     // malloc, free, rand
#include <stdint.h>     // unsigned integers uint8_t, size_t
#include <vector>
#include <limits>       // numeric limits
#include <cmath>
#include <math.h>       // pow
#include <stdexcept>    // for exception, runtime_error, out_of_range

#include <stdarg.h>     // va_list, va_start, va_arg, va_end */
#include <iostream>     // for cout
#include <sstream>      // for peek
#include <algorithm>    // std::min, std::sort()

#include "utilities/mylog.h"

#include <typeinfo>
#include <memory>       // for smart pointer: unique_ptr
#include <random>       // C++11 feature for random number generation

#include <complex>      // C++14 feature for complex numbers

using namespace std::complex_literals;

#define M_PI 3.14159265358979323846

// At the matrix inversion function we have: size_t P[size + 1];
#define SIZE_T_MAX std::numeric_limits<size_t>::max()-1

// Don't increase that unless you have lots of RAM available.
#define MAX_ACCEPTABLE_VECTOR_SIZE 16000

#define MAX(a) ( std::numeric_limits<a>::max() )
#define NaN(a) ( std::numeric_limits<a>::quiet_NaN() )
#define Inf(a) ( std::numeric_limits<a>::infinity() )

#define EPSILON 1e-10  // This will be used for comparison with zero in find_zero() and find_non_zero()

// String to be used for logging errors.
#define FILE_LINE_ERROR std::string(__FILE__) + std::string(":") + std::to_string(__LINE__) + std::string(": ")

namespace algebra {

// Declaration of Vec
template<class T> class Vec;
// Forward Declaration of Mat
template<class T> class Mat;

template <class T>
Mat<T> lup_invert(Mat<T>&, const Vec<T>&);

template <class T>
class Vec {
public:

	explicit Vec();
	Vec(size_t);
	~Vec();

	size_t size() const noexcept;
	size_t max_size() const noexcept;
	size_t size_in_memory() const noexcept; // For debugging purposes
	size_t capacity() const noexcept;
	void set(size_t, T);
	T get(size_t) const;
	Vec<T> get(size_t, size_t) const;
	void set_size(size_t);
	void set_subvector(size_t, const Vec<T>&);
	void zeros();
	void clear();
	void ones();

	Vec<T> add(const Vec<T>&);
	Vec<T> sub(const Vec<T>&);
	T dot(const Vec<T>&);
	Vec<T> cross(const Vec<T>&);
	void sort();
	void swap(size_t, size_t);


	/********** OVERLOAD OPERATORS ***********/

	void operator=(const char*);

	Vec<T> operator+(const Vec<T>&);
	Vec<T> operator+(T);

	Vec<T> operator-(const Vec<T>&);
	Vec<T> operator-(T);

	T operator*(const Vec<T>&);
	Vec<T> operator*(T);
	Vec<T> operator/(T);

	T& operator()(size_t k);
	T& operator[](size_t k);

	// Friend function for lu-decomposition. It will be defined
	// in class Mat.
	friend Mat<T> lup_invert<T>(Mat<T>&, const Vec<T>&);

protected:

private:
	// Since we are not using pointers the destructor will
	// destroy the allocated memory for the object.
	std::vector<T> data_;
	size_t length_ = NaN(size_t);
};

// DEFAULT CONSTRUCTOR
template <class T>
Vec<T>::Vec()
{
	// members should be initialized in the order they were declared;
	data_.resize(0);
	length_ = data_.size();
}

template <class T>
Vec<T>::Vec(size_t n)
{
	if ( n > MAX_ACCEPTABLE_VECTOR_SIZE )
	{
		std::string msg = FILE_LINE_ERROR + " exception in vec(size_t n): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		data_.resize(n);
		length_ = data_.size();
		for (size_t i = 0; i < n; i++)
		{
			data_[i] = T(0);
		}
	}
}

// Destructor
template <class T>
Vec<T>::~Vec() {}


// ##################################################################################################
// ################################# DEFINITIONS OF vec, ivec, cvec #################################

typedef Vec<double> vec;
typedef Vec<int> ivec;
typedef Vec<std::complex<double>> cvec;


// ##################################################################################################
// ############################ FRIENDS AND MEMBER FUNCTIONS DECLARATION ############################


// It computes the number of elements of the vector.
template <class T>
size_t Vec<T>::size() const noexcept { return data_.size(); }

// It computes the memory the vector occupies in Bytes
template <class T>
size_t Vec<T>::size_in_memory() const noexcept{ return sizeof(T)*data_.capacity(); }

// It computes the maximum number of elements vector can contain
template <class T>
size_t Vec<T>::max_size() const noexcept{ return data_.max_size(); }

// It computes the size of the vector's allocated storage capacity
template <class T>
size_t Vec<T>::capacity() const noexcept{ return data_.capacity(); }

// It assigns the i^th element of the vector the value k
template <class T>
void Vec<T>::set(size_t i, T k)
{
	if (i >= length_)
	{
		std::string msg = FILE_LINE_ERROR + " exception in set(size_t k, double v)";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
	else
	{
		data_[i] = k;
	}
}

// It returns the value of the r^th element of the vector
// at(r) from vector class will throw exception if r is
// larger than the size of vector
template <class T>
T Vec<T>::get(size_t r) const
{
	if (length_ == 0)
	{
		std::string msg = FILE_LINE_ERROR + " exception in get(size_t r):: tried to access NULL VECTOR";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
	else
	{
		return data_.at(r);
	}
}

// It returns a vector containing the elements
// from i to j of the current vector.
template <class T>
Vec<T> Vec<T>::get(size_t i, size_t j) const
{
	if (i >= length_ || j >= length_)
	{
		std::string msg = FILE_LINE_ERROR + " exception in get(size_t i, size_t j)";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
	else if(j < i)
	{
		std::string msg = FILE_LINE_ERROR + " exception in get(size_t i, size_t j) ==> j < i";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		Vec<T> result(j - i + 1);
		for(size_t m = i; m < j+1; m++)
		{
			result.data_[m - i] = data_.at(m);
		}
		return result;
	}
}

// It sets the size of the vector
template <class T>
void Vec<T>::set_size(size_t new_size)
{
	if ( new_size > MAX_ACCEPTABLE_VECTOR_SIZE )
	{
		std::string msg = FILE_LINE_ERROR + " exception in set_size(size_t new_size): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else if(length_ == new_size){
		return;
	}
	else
	{
		data_.resize(new_size);
		length_ = data_.size();
	}
}

// Substitutes elements of vector from start-element with
// vector v. They should fit.
template <class T>
void Vec<T>::set_subvector(size_t start, const Vec<T>& v)
{
	if (start >= length_ || length_ - start < v.size())
	{
		std::string msg = FILE_LINE_ERROR + " exception in set_subvector(size_t start, const Vec<T>& v)";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
	else
	{
		size_t i, size = v.size();
		for(i = size; i--;)
		{
			data_[i+start] = v.data_[i];
		}
	}
}

// It sets all elements of vector to 0.
template <class T>
void Vec<T>::zeros()
{
	size_t i = 0;
	for(i = length_; i--;){
		data_[i] = T(0);
	}
}

// It sets all elements of vector to 0.
template <class T>
void Vec<T>::clear(){ zeros(); }

// It sets all elements of vector to 1.
template <class T>
void Vec<T>::ones(){
	size_t i = 0;
	for (i = length_; i--;)
	{
		data_[i] = T(1);
	}
}

// It adds two vectors.
template <class T>
Vec<T> Vec<T>::add(const Vec<T>& v1)
{
	if (this->length_ != v1.length_)
	{
		std::string msg = FILE_LINE_ERROR + " dimension mismatch in add(const vec& v1)";
		log_error(msg.c_str());
		throw std::length_error(msg);
	}
	else
	{
		Vec<T> result(v1.length_);
		size_t i, size = result.size();
		for (i = size; i--;)
		{
			result.data_[i] = data_[i] + v1.data_[i];
		}
		return result;
	}
}

//It subtracts two vectors.
template <class T>
Vec<T> Vec<T>::sub(const Vec<T>& v1)
{
	if (this->length_ != v1.length_)
	{
		std::string msg = FILE_LINE_ERROR + " dimension mismatch in sub(const vec& v1)";
		log_error(msg.c_str());
		throw std::length_error(msg);
	}
	else
	{
		Vec<T> result(v1.length_);
		size_t i, size = result.size();
		for (i = size; i--;)
		{
			result.data_[i] = data_[i] - v1.data_[i];
		}
		return result;
	}
}

// It computes the inner product of the vectors.
template <class T>
T Vec<T>::dot(const Vec<T>& v1)
{
	if (this->length_ != v1.length_)
	{
		std::string msg = FILE_LINE_ERROR + " dimension mismatch in dot(const vec& v1)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		T result = 0;
		size_t i = 0, length = v1.length_;
		for (i = length; i--;)
		{
			result += data_[i] * v1.data_[i];
		}
		return result;
	}
}

// It computes the outer product of two vectors. The size of the
// vectors has to be 3 or 7. I only consider the trivial case of 3.
// For more details see: https://en.wikipedia.org/wiki/Cross_product
// Cross product is not defined for complex vectors.
template <class T>
Vec<T> Vec<T>::cross(const Vec<T>& v1)
{
	if ( this->length_ == 3 && v1.length_ == 3 )
	{
		Vec<T> result(3);
		result.data_[0] = data_[1] * v1.data_[2] - v1.data_[1] * data_[2];
		result.data_[1] = data_[2] * v1.data_[0] - v1.data_[2] * data_[0];
		result.data_[2] = data_[0] * v1.data_[1] - v1.data_[0] * data_[1];
		return result;
	}
	else
	{
		std::string msg = FILE_LINE_ERROR + " Vectors should be of size 3 in cross(const vec& v1)";
		log_error(msg.c_str());
		throw std::length_error(msg);
	}
}


// It computes a vector with sorted elements.
template <class T>
void Vec<T>::sort(){ std::sort(data_.begin(), data_.end()); }

// It swaps elements i and j of the vector
template <class T>
void Vec<T>::swap(size_t i, size_t j)
{
	if (i < length_ && j < length_ && i >= 0 && j >= 0)
	{
		T tmp = data_[i];
		data_[i] = data_[j];
		data_[j] = tmp;
	}
	else
	{
		std::string msg = FILE_LINE_ERROR + " exception in vec::swap(size_t i, size_t j): out of range indices.";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
}

/********** OVERLOAD OPERATORS ***********/

// It pass the values of the string in a vector.
// The input must be of the form "1 2 3" or "[1 2 3]".
template <class T>
void Vec<T>::operator=(const char* a)
{
	// Clear vector from any previous values
	if ( (*this).size() != 0 )
	{
		data_.resize(0);
		length_ = 0;
	}

	//Convert input to string
	std::string str = std::string(a);

	// Remove brackets if there are any.
	if (str.size() >= 1)
	{
		std::string start = str.substr(0,1);
		std::string end = str.substr(str.size()-1, str.size());
		if ( start.compare("[") == 0 && end.compare("]") == 0)
		{
			str.erase(0,1);
			str.erase(str.size()-1, str.size());
		}
	}

	// Check if str has only numbers. If not, an exception should be thrown.
	std::string acceptable_characters = " 0123456789.+-";
	if (str.find_first_not_of(acceptable_characters) != std::string::npos)
	{
		std::string msg = FILE_LINE_ERROR + " exception in operator =(const char* a): vec can contain only numbers.";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}

	// Check if it's an empty string (if there is no number in the string, consider it an empty vector)
	std::string numbers = "0123456789";
	if ( str.find_first_of(numbers) == std::string::npos )
	{
		data_.resize(0);
		length_ = 0;
		return;
	}

	std::string delimiter = " "; // DEFAULT DELIMITER
	size_t pos = 0;
	size_t k = 0;
	std::string tmp_str = str;

	// Determine vector size
	while ((pos = tmp_str.find(delimiter)) != std::string::npos)
	{
		tmp_str.erase(0, pos + delimiter.length());
		k++;
	}

	// Initialize vector
	data_.resize(k+1);
	length_ = data_.size();

	k = 0;
	pos = 0;
	std::string token;

	// Pass each number into the vector
	while ((pos = str.find(delimiter)) != std::string::npos)
	{
		token = str.substr(0, pos);
		data_[k] = (T) strtod((token).c_str(),0);
		str.erase(0, pos + delimiter.length());
		k++;
	}
	data_[k] = (T) strtod((str).c_str(),0);
}

// It adds vectors v to the current vector (Commutative property holds)
template <class T>
Vec<T> Vec<T>::operator+(const Vec<T>& v)
{
	if (this->length_ != v.length_)
	{
		std::string msg = FILE_LINE_ERROR + " Dimension mismatch for operator+(const vec& v1)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		Vec<T> result = *this;
		size_t i = 0, size = result.size();
		for (i = size; i--;)
		{
			result.data_[i] += v.data_[i];
		}
		return result;
	}
}

// It adds the value 't' to each and every element of the vector.
template <class T>
Vec<T> Vec<T>::operator+(T t)
{
	if ( (*this).size() == 0 )
	{
		std::string msg = FILE_LINE_ERROR + " exception in operator+(double t): NULL VECTOR";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
	else
	{
		Vec<T> result = *this;
		size_t i = 0, size = result.size();
		for (i = size; i--;)
		{
			result.data_[i] += t;
		}
		return result;
	}
}

// It subtracts vector v from the current vector
template <class T>
Vec<T> Vec<T>::operator-(const Vec<T>& v)
{
	if (this->length_ != v.length_)
	{
		std::string msg = FILE_LINE_ERROR + " Dimension mismatch for operator-(const vec& v1)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		Vec<T> result = *this;
		size_t i = 0, size = result.size();
		for (i = size; i--;)
		{
			result.data_[i] -= v.data_[i];
		}
		return result;
	}
}

// It subtracts the value 't' from each and every element of the vector.
template <class T>
Vec<T> Vec<T>::operator-(T t)
{
	if ( (*this).size() == 0 )
	{
		std::string msg = FILE_LINE_ERROR + " exception in operator-(double t): NULL VECTOR";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
	else
	{
		Vec<T> result = *this;
		size_t i = 0, size = result.size();
		for (i = size; i--;)
		{
			result.data_[i] -= t;
		}
		return result;
	}
}

// It returns the inner product of the vectors
template <class T>
T Vec<T>::operator*(const Vec<T>& v)
{
	if (length_ != v.length_)
	{
		std::string msg = FILE_LINE_ERROR + " Dimension mismatch for operator*(const vec& v1)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);

	}
	else
	{
		return this->dot(v);
	}
}

// It multiplies each elements of the vector with x
template <class T>
Vec<T> Vec<T>::operator*(T x)
{
	Vec<T> result = *this;
	size_t i = 0, size = result.size();
	for (i = size; i--;)
	{
		result.data_[i] *= x;
	}
	return result;
}

// It devides each element of the vector with t.
// I throw exception if t = 0 (undefined behavior in some cases)
template <class T>
Vec<T> Vec<T>::operator/(T t)
{
	if (t == T(0))
	{
		std::string msg = FILE_LINE_ERROR + " 'std::invalid_argument' thrown in operator/(T t): DIVISION BY ZERO ";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		Vec<T> result = *this;
		size_t i = 0, size = result.size();
		for (i = size; i--;)
		{
			result.data_[i] /= t;
		}
		return result;
	}
}

// It returns the k^th element of the vector.
template <class T>
T& Vec<T>::operator()(size_t k)
{
	if ( (*this).size() == 0 )
	{
		std::string msg = FILE_LINE_ERROR + " exception in operator()(const size_t k): tried to access NULL VECTOR ";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
	else if (k >= length_)
	{
		std::string msg = FILE_LINE_ERROR + " exception in operator()(const size_t k): index > vector size ";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		return data_.at(k);;
	}
}

// It returns the k^th element of the vector.
template <class T>
T& Vec<T>::operator[](size_t k)
{
	if ( (*this).size() == 0 )
	{
		std::string msg = FILE_LINE_ERROR + " exception in operator[](const size_t k): tried to access NULL VECTOR ";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
	else if (k >= length_)
	{
		std::string msg = FILE_LINE_ERROR + " exception in operator[](const size_t k): index > vector size ";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		return data_.at(k);;
	}
}

// ##################################################################################################
// ############################ MISCELLANEOUS OPERATIONS AND FUNCTIONS ##############################


// It returns a 'double'-vector of size n with all elements equal to 0.
inline vec zeros(size_t n)
{
	if ( n > MAX_ACCEPTABLE_VECTOR_SIZE )
	{
		std::string msg = FILE_LINE_ERROR + " exception in vec::zeros(size_t n): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		vec a(n);
		size_t i;
		for (i = n; i--;)
		{
			a(i) = 0.0;
		}
		return a;
	}
}

// It returns an 'integer'-vector of size n with all elements equal to 0.
inline ivec zeros_i(size_t n)
{
	if ( n > MAX_ACCEPTABLE_VECTOR_SIZE )
	{
		std::string msg = FILE_LINE_ERROR + " exception in vec::zeros(size_t n): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		ivec a(n);
		size_t i;
		for (i = n; i--;)
		{
			a(i) = 0;
		}
		return a;
	}
}

// It returns a 'double'-vector of size n with all elements equal to 1.
inline vec ones(size_t n)
{
	if ( n > MAX_ACCEPTABLE_VECTOR_SIZE )
	{
		std::string msg = FILE_LINE_ERROR + " exception in ones(size_t n): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		vec a(n);
		size_t i;
		for (i = n; i--;)
		{
			a(i) = 1.0;
		}
		return a;
	}
}

// It returns a 'integer'-vector of size n with all elements equal to 1.
inline ivec ones_i(size_t n)
{
	if ( n > MAX_ACCEPTABLE_VECTOR_SIZE )
	{
		std::string msg = FILE_LINE_ERROR + " exception in ones(size_t n): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		ivec a(n);
		size_t i;
		for (i = n; i--;)
		{
			a(i) = 1;
		}
		return a;
	}
}

// It computes a mask with ones where the input vector
// has non-zero elements and zero elsewhere.
template <class T>
inline Vec<T> find_non_zero(const Vec<T>& v)
{
	size_t i, size = v.size();
	Vec<T> result(size);
	for (i = size; i--;)
	{
		if ( fabs(v.get(i) - T(0)) < EPSILON )
		{
			result.set(i, T(0));
		}
		else
		{
			result.set(i, T(1));
		}
	}
	return result;
}

// It computes a mask with ones where the input vector
// has zero elements and zero elsewhere.
template <class T>
inline Vec<T> find_zero(const Vec<T>& v)
{
	size_t i, size = v.size();
	Vec<T> result(size);
	for (i = size; i--;)
	{
		if ( fabs(v.get(i) - T(0)) < EPSILON )
		{
			result.set(i, T(1));
		}
		else
		{
			result.set(i, T(0));
		}
	}
	return result;
}

// It calculates a 'double'-vector with elements taking
// random values within the range [-10, 10].
inline vec rand(size_t n)
{
	if ( n > MAX_ACCEPTABLE_VECTOR_SIZE )
	{
		std::string msg = FILE_LINE_ERROR + " exception in rand(size_t n): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		vec a(n);
		size_t i = 0;

		// C++11 feature: It will be used to obtain a seed for the random number engine
		std::random_device rd;
		// C++11 feature: Standard mersenne_twister_engine seeded with rd()
		std::mt19937 gen(rd());

		// Define range
		double min = -10, max = 10;

		std::uniform_real_distribution<double> dis(0, 2*max);
		for (i = n; i--;)
		{
			// Generate random number within the range [-10 10]
			a(i) = min + (double) dis(gen);
		}
		return a;
	}
}

// It calculates an 'integer'-vector with elements taking
// random values within the range [-10, 10].
inline ivec rand_i(size_t n)
{
	if ( n > MAX_ACCEPTABLE_VECTOR_SIZE )
	{
		std::string msg = FILE_LINE_ERROR + " exception in rand_i(size_t n): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		ivec a(n);
		size_t i = 0;

		// C++11 feature: It will be used to obtain a seed for the random number engine
		std::random_device rd;
		// C++11 feature: Standard mersenne_twister_engine seeded with rd()
		std::mt19937 gen(rd());

		// Define range
		int min = -10, max = 10;

		std::uniform_int_distribution<int> dis(0, 2*max);
		for (i = n; i--;)
		{
			// Generate random number within the range [-10 10]
			a(i) = min + (int) dis(gen);
		}
		return a;
	}
}

// It computes the inner product of v1, v2.
template <class T>
inline T dot(const Vec<T>& v1, const Vec<T>& v2)
{
	if ( v1.size() != v1.size() )
	{
		std::string msg = FILE_LINE_ERROR + " dimension mismatch in dot(const vec& v1, const vec& v2)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		T result = 0;
		size_t i = 0, size = v1.size();
		for (i = size; i--;)
		{
			result += v1.get(i)*v2.get(i);
		}
		return result;
	}
}

// It computes the mean value of the vector
template <class T>
inline double mean(const Vec<T>& v)
{
	if (v.size() == 0)
	{
		std::string msg = FILE_LINE_ERROR + " NULL VECTOR in mean(const vec& v)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		double result = 0;
		size_t i = 0, size = v.size();
		for (i = size; i--;)
		{
			result += v.get(i);
		}
		result /= size;
		return result;
	}
}

// It computes the min value of a vector
template <class T>
inline T min(const Vec<T>& v)
{
	if (v.size() == 0)
	{
		std::string msg = FILE_LINE_ERROR + " NULL VECTOR in min(const vec& v)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		T min = MAX(T);
		size_t i = 0, size = v.size();
		for (i = size; i--;)
		{
			if (v.get(i) < min )
			{
				min = v.get(i);
			}
		}
		return min;
	}
}

// It computes the min value of a vector together with its index
template <class T>
inline T min(const Vec<T>& v, size_t &index)
{
	if (v.size() == 0)
	{
		std::string msg = FILE_LINE_ERROR + " NULL VECTOR in min(const vec& v, size_t &index)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		index = 0;
		T min = MAX(T);;
		// Calculate the index of the minimum value
		size_t i = 0, size = v.size();
		for (i = size; i--;)
		{
			if (v.get(i) < min )
			{
				min = v.get(i);
				index = i;
			}
		}
		return v.get(index);
	}
}

// It computes the max value of a vector
template <class T>
inline T max(const Vec<T>& v)
{
	if (v.size() == 0)
	{
		std::string msg = FILE_LINE_ERROR + " NULL VECTOR in max(const vec& v)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		T max = 0;
		size_t i = 0, size = v.size();
		for (i = size; i--;){
			if(v.get(i) > max ){
				max = v.get(i);
			}
		}
		return max;
	}
}

// It computes the max value of a vector together with its index
template <class T>
inline T max(const Vec<T>& v, size_t &index)
{
	if (v.size() == 0)
	{
		std::string msg = FILE_LINE_ERROR + " NULL VECTOR in max(const vec& v, size_t &index)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		index = 0;
		T max = 0;
		// Calculate the index of the maximum value
		size_t i = 0, size = v.size();
		for (i = size; i--;)
		{
			if (v.get(i) > max )
			{
				max = v.get(i);
				index = i;
			}
		}
		return max;
	}
}

// It computes the cross product of v1, v2. Vectors must be of size 3
template <class T>
inline Vec<T> cross(const Vec<T>& v1, const Vec<T>& v2)
{
	if ( v1.size() == 3 && v2.size() == 3 )
	{
		vec result(3);
		result.set(0, v1.get(1)*v2.get(2) - v2.get(1)*v1.get(2));
		result.set(1, v1.get(2)*v2.get(0) - v2.get(2)*v1.get(0));
		result.set(2, v1.get(0)*v2.get(1) - v2.get(0)*v1.get(1));
		return result;
	}
	else
	{
		std::string msg = FILE_LINE_ERROR + " Vectors should be of size 3 in cross(const vec& v1, const vec& v2)";
		log_error(msg.c_str());
		throw std::length_error(msg);
	}
}

// It concatenates the current vector with t.
template <class T>
inline Vec<T> concat(const Vec<T>& v, T t)
{
	if ( v.size() == 0 )
	{
		Vec<T> result(1);
		result.set(0, t);
		return result;
	}
	else if ( v.size() == MAX_ACCEPTABLE_VECTOR_SIZE )
	{
		std::string msg = FILE_LINE_ERROR + " Inputs define out of range vector in concat(const vec& v, double t)";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
	else
	{
		Vec<T> result(v.size() + 1);
		size_t i = 0, size = result.size();
		for (i = size - 1; i--;)
		{
			result(i) = v.get(i);
		}
		result(size - 1) = t;
		return result;
	}
}

// It concatenates t with the current vector.
template <class T>
inline Vec<T> concat(T t, const Vec<T>& v)
{
	if ( v.size() == 0 )
	{
		Vec<T> result(1);
		result.set(0, t);
		return result;
	}
	else if ( v.size() == MAX_ACCEPTABLE_VECTOR_SIZE )
	{
		std::string msg = FILE_LINE_ERROR + " Inputs define out of range vector in concat(double t, const vec& v)";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
	else
	{
		Vec<T> result(v.size() + 1);
		size_t i = 0, size = result.size();
		result(0) = t;
		for (i = size - 1; i--;)
		{
			result(i+1) = v.get(i);
		}
		return result;
	}
}

// It concatenates vectors v1 and v2.
template <class T>
inline Vec<T> concat(const Vec<T>& v1, const Vec<T>& v2)
{
	if (v1.size() == 0 && v2.size() == 0)
	{
		// Concatenation of null vectors should yield back a null vector
		return Vec<T>(0);
	}
	else if (v1.size() + v2.size() > MAX_ACCEPTABLE_VECTOR_SIZE)
	{
		std::string msg = FILE_LINE_ERROR + " Inputs define out of range vector in concat(const vec& v1, const vec& v2)";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
	else
	{
		Vec<T> result(v1.size() + v2.size());
		size_t i = 0;
		size_t size = result.size(), size_v1 = v1.size();
		for (i = size; i--;)
		{
			if (i < size_v1)
			{
				result(i) = v1.get(i);
			}
			else
			{
				result(i) = v2.get(i - size_v1);
			}
		}
		return result;
	}
}

// It returns a vector with elements equally spaced by step.
// E.g.: linspace(0, 5 , 1) = [0 1 2 3 4 5]
template <class T>
inline Vec<T> linspace(T from, T to, size_t step)
{
	size_t size = (size_t) (std::floor((to -from)/step) + 1);
	if ( (to - from) < 0 || step < 0 )
	{
		std::string msg = FILE_LINE_ERROR + " Invalid arguments in linspace(double from, double to, size_t step)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else if (size > MAX_ACCEPTABLE_VECTOR_SIZE)
	{
		std::string msg = FILE_LINE_ERROR + " Inputs define out of range vector in linspace(double from, double to, size_t step)";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
	else
	{
		Vec<T> result(size);
		size_t i = 0, size = result.size();
		for (i = size; i--;)
		{
			result.set(i, from + step*i);
		}
		return result;
	}
}

// It computes the element-wise vector multiplication
template <class T>
inline Vec<T> elem_mult(const Vec<T>& v1, const Vec<T>& v2)
{
	if (v1.size() != v2.size())
	{
		std::string msg = FILE_LINE_ERROR + " Dimension mismatch in elem_mult(const vec& v1, const vec& v2)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		Vec<T> result(v1.size());
		size_t i = 0, size = result.size();
		for (i = size; i--;)
		{
			result.set(i, v1.get(i)*v2.get(i));
		}
		return result;
	}
}

// It computes the sum of all elements of vector v
template <class T>
inline T sum(const Vec<T>& v)
{
	T result = 0;
	size_t i = 0, size = v.size();
	for (i = size; i--;)
	{
		result += v.get(i);
	}
	return result;
}

// It computes the cumulative sum of v1
template <class T>
inline Vec<T> cumsum(const Vec<T>& v1)
{
	if (v1.size() <= 1)
	{
		return Vec<T>(v1.size());
	}
	else
	{
		Vec<T> result(v1.size());
		result.set(0, v1.get(0));
		size_t i = 0, size = result.size();
		for (i = 1; i < size; i++)
		{
			result(i) = result(i-1) + v1.get(i);
		}
		return result;
	}
}

// It computes the 2-norm for Euclidean space
template <class T>
inline double norm(const Vec<T>& v)
{
	double result = 0;
	size_t i = 0, size = v.size();
	for (i = size; i--;)
	{
		result += v.get(i)*v.get(i);
	}
	return std::sqrt(result);
}

// It returns a vector whose i^th element equals
// to the absolute value of i^th element of vector v
template <class T>
inline Vec<T> abs(const Vec<T>& v)
{
	Vec<T> result(v.size());
	size_t i = 0, size = v.size();
	for (i = size; i--;)
	{
		result(i) = std::abs(v.get(i));
	}
	return result;
}

// It prints the elements of the vector
template <class T>
inline void print(const Vec<T>& v)
{
	if (v.size() == 0)
	{
		std::cout << "[ ]" << std::endl;
		return;
	}
	else
	{
		std::cout << "[ ";
		size_t i;
		for (i = 0; i < v.size(); i++)
		{
			std::cout << v.get(i) << " ";
		}
		std::cout << "]" << std::endl;
	}
}

// ##################################################################################################
// ########################### COMPLEX NUMBER OPERATIONS AND FUNCTIONS ##############################

// It returns an 'complex'-vector of size n with all elements equal to 0.
inline cvec zeros_c(size_t n)
{
	if ( n > MAX_ACCEPTABLE_VECTOR_SIZE )
	{
		std::string msg = FILE_LINE_ERROR + " exception in vec::zeros(size_t n): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		cvec a(n);
		size_t i;
		for (i = n; i--;)
		{
			a(i).real(0);
			a(i).imag(0);
		}
		return a;
	}
}

// It calculates an 'complex'-vector with elements taking
// random values within the range [-10, 10].
inline cvec rand_c(size_t n)
{
	if ( n > MAX_ACCEPTABLE_VECTOR_SIZE )
	{
		std::string msg = FILE_LINE_ERROR + " exception in rand_c(size_t n): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		cvec a(n);
		size_t i = 0;

		// C++11 feature: It will be used to obtain a seed for the random number engine
		std::random_device rd;
		// C++11 feature: Standard mersenne_twister_engine seeded with rd()
		std::mt19937 gen(rd());

		// Define range
		double min = -10, max = 10;

		std::uniform_real_distribution<double> dis(0, 2*max);
		for (i = n; i--;)
		{
			// Generate random number within the range [-10 10]
			a(i).real(min + (double) dis(gen));
			a(i).imag(min + (double) dis(gen));
		}
		return a;
	}
}

// It computes the mean value of the complex vector
inline std::complex<double> mean(const cvec& v)
{
	if (v.size() == 0)
	{
		std::string msg = FILE_LINE_ERROR + " NULL VECTOR in mean(const cvec& v)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		std::complex<double> result = 0. + 0i;
		size_t i = 0, size = v.size();
		for (i = size; i--;)
		{
			result += v.get(i);
		}
		result /= size;
		return result;
	}
}

// It computes the min value of a complex vector
// Based on MATLAB's implementation I compare
// the modulus of each element and in case they
// are equal I pick the one with the minimum phase
inline std::complex<double> min(const cvec& v)
{
	if (v.size() == 0)
	{
		std::string msg = FILE_LINE_ERROR + " NULL VECTOR in min(const vec& v)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		std::complex<double> min;
		min.real( MAX(double) );
		min.imag( MAX(double) );
		size_t i = 0, size = v.size();
		double modulus = 0;
		double phase = 0;
		for (i = size; i--;)
		{
			modulus = std::abs(v.get(i));
			phase = atan2(v.get(i).imag(), v.get(i).real());
			if (modulus < std::abs(min) )
			{
				min = v.get(i);
			}
			if (std::abs(modulus - std::abs(min)) < EPSILON && phase < atan2(min.imag(), min.real()))
			{
				min = v.get(i);
			}
		}
		return min;
	}
}

// It computes the min value of a complex vector
// together with its index.
inline std::complex<double> min(const cvec& v, size_t &index)
{
	if (v.size() == 0)
	{
		std::string msg = FILE_LINE_ERROR + " NULL VECTOR in min(const vec& v)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		std::complex<double> min;
		min.real( MAX(double) );
		min.imag( MAX(double) );
		size_t i = 0, size = v.size();
		double modulus = 0;
		double phase = 0;
		index = 0;
		for (i = size; i--;)
		{
			modulus = std::abs(v.get(i));
			phase = atan2(v.get(i).imag(), v.get(i).real());
			if (modulus < std::abs(min) )
			{
				min = v.get(i);
				index = i;
			}
			if (std::abs(modulus - std::abs(min)) < EPSILON && phase < atan2(min.imag(), min.real()))
			{
				min = v.get(i);
				index = i;
			}
		}
		return min;
	}
}

// It computes the max value of a complex vector
// Based on MATLAB's implementation I compare
// the modulus of each element and in case they
// are equal I pick the one with the minimum phase
inline std::complex<double> max(const cvec& v)
{
	if (v.size() == 0)
	{
		std::string msg = FILE_LINE_ERROR + " NULL VECTOR in min(const vec& v)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		std::complex<double> max;
		max.real( 0 );
		max.imag( 0 );
		size_t i = 0, size = v.size();
		double modulus = 0;
		double phase = 0;
		for (i = size; i--;)
		{
			modulus = std::abs(v.get(i));
			phase = atan2(v.get(i).imag(), v.get(i).real());
			if (modulus > std::abs(max) )
			{
				max = v.get(i);
			}
			if (std::abs(modulus - std::abs(max)) < EPSILON && phase > atan2(max.imag(), max.real()))
			{
				max = v.get(i);
			}
		}
		return max;
	}
}

// It computes the max value of a complex vector
// together with its index.
inline std::complex<double> max(const cvec& v, size_t &index)
{
	if (v.size() == 0)
	{
		std::string msg = FILE_LINE_ERROR + " NULL VECTOR in min(const vec& v)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	else
	{
		std::complex<double> max;
		max.real( 0 );
		max.imag( 0 );
		size_t i = 0, size = v.size();
		double modulus = 0;
		double phase = 0;
		index = 0;
		for (i = size; i--;)
		{
			modulus = std::abs(v.get(i));
			phase = atan2(v.get(i).imag(), v.get(i).real());
			if (modulus > std::abs(max) )
			{
				max = v.get(i);
				index = i;
			}
			if (std::abs(modulus - std::abs(max)) < EPSILON && phase > atan2(max.imag(), max.real()))
			{
				max = v.get(i);
				index = i;
			}
		}
		return max;
	}
}


// It computes the conjugate of a complex vector
inline cvec conj(const cvec& v)
{
	size_t i, size = v.size();
	cvec result = v;
	for (i = size; i--;)
	{
		result(i).imag(result(i).imag()*-1);
	}
	return result;
}

// It computes the 2-norm for complex space
inline double norm(const cvec& v)
{
	cvec tmp = v;
	return sqrt( (tmp*conj(tmp)).real() );
}





inline void print(const cvec& v)
{
	if (v.size() == 0)
	{
		std::cout << "[ ]" << std::endl;
		return;
	}
	else
	{
		std::cout << "[ ";
		size_t i;
		for (i = 0; i < v.size()-1; i++)
		{
			if(v.get(i).imag() >= 0)
			{
				std::cout << v.get(i).real() << "+" << std::abs(v.get(i).imag()) << "i" << "  ";
			}
			else
			{
				std::cout << v.get(i).real() << "-" << std::abs(v.get(i).imag()) << "i" << "  ";
			}

		}
		if(v.get(v.size()-1).imag() >= 0)
		{
			std::cout << v.get(v.size()-1).real() << "+" << std::abs(v.get(v.size()-1).imag()) << "i" << " ";
		}
		else
		{
			std::cout << v.get(v.size()-1).real() << "-" << std::abs(v.get(v.size()-1).imag()) << "i" << " ";
		}
		std::cout << "]" << std::endl;
	}
}

// ======================== UNIT TEST IS WORKING FINE UP TO THIS POINT ========================

} /* namespace algebra */

#endif /* VEC_H_ */
