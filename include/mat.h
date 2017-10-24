/*===========================================================================
 * Name         : mat.h implements basic matrix operations and
 *                can be used for linear algebra projects.
 * Version      : 1.0.0, 16 Sep 2017
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
===============================================================================*/

#ifndef MAT_H_
#define MAT_H_


#include "vec.h"

// The absolute value of the determinant should
//be above that threshold for a matrix to be invertible.
#define SINGULARITY_THRESHOLD 1e-9

namespace algebra {

// Declaration of vec
template<class T> class Vec;
// Declaration of mat
template<class T> class Mat;

// Declaration of friend functions
template <class T>
Mat<T> transpose(const Mat<T>&);
template <class T>
Vec<T> lup_decompose(Mat<T>&, bool& is_singular);
template <class T>
Mat<T> lup_invert(Mat<T>&, const Vec<T>&);
template <class T>
Mat<T> inv(const Mat<T>&);
template <class T>
Mat<T> pinv(Mat<T>&); // pseudoinverse
template <class T>
Mat<T> strassen_algorithm(const Mat<T>& , const Mat<T>&, size_t leafsize );
template <class T>
Mat<T> strassen(const Mat<T>&, const Mat<T>& );

template <class T>
class Mat {
public:
	typedef T type;

	explicit Mat();
	Mat(size_t, size_t);
	~Mat();

	T get(size_t, size_t) const;

	size_t size() const noexcept;
	size_t rows() const noexcept;
	size_t cols() const noexcept;
	size_t size_in_memory() const noexcept; // For debugging purposes
	void set_size(size_t, size_t);
	void set(size_t, size_t, T);
	void set_row(size_t, const Vec<T>&);
	void set_col(size_t, const Vec<T>&);
	void set_rows(size_t, const Mat<T>&);
	void set_cols(size_t, const Mat<T>&);

	void set_submatrix(size_t, size_t, const Mat<T>&);

	Vec<T> get_col(size_t) const;
	Mat<T> get_cols(size_t, size_t) const;
	Vec<T> get_row(size_t) const;
	Mat<T> get_rows(size_t, size_t) const;
	Mat<T> get(size_t, size_t, size_t, size_t) const;
	Mat<T> get(const Vec<T>&, const Vec<T>&) const;

	void zeros();
	void clear();
	void ones();

	void swap_rows(size_t, size_t);
	void swap_cols(size_t, size_t);


	/********** OVERLOAD OPERATORS ***********/
	void operator=(const char* a);

	T& operator()(size_t i, size_t j);
	Mat<T> operator()(size_t r1, size_t r2, size_t c1, size_t c2);

	// overload +
	// add input matrix to current matrix
	Mat<T> operator+(const Mat<T>& );
	Mat<T> operator+(T);

	// overload -
	// subtract input matrix from current matrix
	Mat<T> operator-(const Mat<T>&);
	Mat<T> operator-(T);

	// overload *
	// multiply current matrix with scalar
	Mat<T> operator*(T);

	Mat<T> operator*(const Mat<T>&);

	Vec<T> operator*(const Vec<T>&);

	Mat<T> operator/(T t);

	void print(int precision);
	void print();

	// Declaration of friend functions
	friend Mat<T> transpose<>(const Mat<T>&);
	friend Vec<T> lup_decompose<>(Mat<T>&, bool& is_singular);
	friend Mat<T> lup_invert<>(Mat<T>&, const Vec<T>&);
	friend Mat<T> inv<>(const Mat<T>&);
	friend Mat<T> pinv<>(Mat<T>&); // pseudoinverse
	friend Mat<T> strassen_algorithm<T>(const Mat<T>&, const Mat<T>&, size_t leafsize );
	friend Mat<T> strassen<>(const Mat<T>&, const Mat<T>& );

protected:

private:
	std::vector< std::vector<T> > data_;
	size_t rows_ = NaN_NUMERIC_LIMIT(size_t);
	size_t cols_ = NaN_NUMERIC_LIMIT(size_t);
};


// ##################################################################################################
// #################################### DEFINITIONS OF vec, ivec ####################################

typedef Mat<double> mat;
typedef Mat<int> imat;


// ##################################################################################################
// ############################ FRIENDS AND MEMBER FUNCTIONS DECLARATION ############################



// DEFAULT CONSTUCTOR
template <class T>
Mat<T>::Mat() {
	data_.resize(0, std::vector<T> (0));
	rows_ = 0;
	cols_ = 0;
}

template <class T>
Mat<T>::Mat(size_t r, size_t c) {
	if(( (r*c) > MAX_ACCEPTABLE_VECTOR_SIZE*MAX_ACCEPTABLE_VECTOR_SIZE) ){
		std::string msg = FILE_LINE_ERROR + "exception in mat(size_t r, size_t c): (r*c) should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{
		data_.resize(r, std::vector<T> (c));
		rows_ = r;
		cols_ = c;
		// Initialization
		size_t i,j;
		for(i = r; i--;){
			for(j = c; j--;){
				data_[i][j] = 0;
			}
		}
	}
}

template <class T>
Mat<T>::~Mat() {
	// TODO Auto-generated destructor stub
}

template <class T>
T Mat<T>::get(size_t r1, size_t c1) const{
	if(rows_ == 0 || cols_ == 0){
		std::string msg = FILE_LINE_ERROR + "exception in  mat::get(size_t r1, size_t c1): NULL MATRIX";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else if(r1 >= rows_ || c1 >= cols_){
		std::string msg = FILE_LINE_ERROR + "exception in  mat::get(size_t r1, size_t c1): index exceeds size of matrix";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		return data_[r1][c1];
	}
}

template <class T>
size_t Mat<T>::size() const noexcept{
	return rows_*cols_;
}

template <class T>
size_t Mat<T>::rows() const noexcept{
	return rows_;
}

template <class T>
size_t Mat<T>::cols() const noexcept{
	return cols_;
}

template <class T>
size_t Mat<T>::size_in_memory() const noexcept{
	return (*this).size()*sizeof(double);
}

template <class T>
void Mat<T>::set_size(size_t r, size_t c){
	if( r < 0 || c < 0 || (r*c) >  MAX_ACCEPTABLE_VECTOR_SIZE*MAX_ACCEPTABLE_VECTOR_SIZE ){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::set_size(size_t r1, size_t c1): (r*c) should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else if( r == 0 || c == 0){
		data_.resize(0, std::vector<double> (0));
		rows_ = 0;
		cols_ = 0;
	}else{
		data_.resize(r, std::vector<double> (c));
		rows_ = r;
		cols_ = c;
		size_t i,j;
		for(i = r; i--;){
			for(j = c; j--;){
				data_[i][j] = 0;
			}
		}
	}
}

template <class T>
void Mat<T>::set(size_t r, size_t c, T value){
	if((r < 0) || (c < 0) || (r >= rows_) || (c >= cols_)){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::set(size_t r1, size_t c1): Indices out of bounds";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		data_[r][c] = value;
	}
}

template <class T>
void Mat<T>::set_row(size_t r, const Vec<T>& v1){
	if( (r < 0) || (r >= rows_) || (v1.size() != cols_) ){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::set_row(size_t r, const vec& v1): Erroneous index or dimension mismatch";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		size_t j;
		for(j = cols_; j--;){
			data_[r][j] = v1.get(j);
		}
	}
}

template <class T>
void Mat<T>::set_col(size_t c, const Vec<T>& v1){
	if((c < 0) || (c >= cols_) || (v1.size() != rows_)){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::set_col(size_t r, const vec& v1): Erroneous index or dimension mismatch";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		size_t i;
		for(i = rows_; i--;){
			data_[i][c] = v1.get(i);
		}
	}
}

template <class T>
void Mat<T>::set_rows(size_t r0, const Mat<T>& m1){
	if((r0 < 0) || (r0 >= rows_) || (m1.rows() > rows_ - r0) || (m1.cols() > cols_)){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::set_rows(size_t r0, const mat& m1): Erroneous index or dimension mismatch";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		size_t i, j, rows = m1.rows(), cols = m1.cols();
		for(i = rows; i--;){
			for(j = cols; j--;){
				data_[i+r0][j] = m1.data_[i][j];
			}
		}
	}
}

template <class T>
void Mat<T>::set_cols(size_t c0, const Mat<T>& m1){
	if((c0 < 0) || (c0 >= cols_) || (m1.cols() > cols_ - c0) || (m1.rows() > rows_)){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::set_cols(size_t r0, const mat& m1): Erroneous index or dimension mismatch";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		size_t i, j, rows = m1.rows(), cols = m1.cols();
		for( i = rows; i--;){
			for(j = cols; j--;){
				data_[i][j+c0] = m1.data_[i][j];
			}
		}
	}
}

template <class T>
void Mat<T>::set_submatrix(size_t r0, size_t c0, const Mat<T>& m){
	if(r0 < 0 || r0 > rows_ || (rows_ - r0) < m.rows() || c0 < 0 || c0 > cols_ || (cols_ - c0) < m.cols()){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::set_submatrix(size_t r0, size_t c0, const mat& m): Erroneous index or dimension mismatch";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		size_t i, j, rows = m.rows(), cols = m.cols();
		for(i = rows; i--;){
			for(j = cols; j--;){
				data_[i+c0][j+c0] = m.data_[i][j];
			}
		}
	}
}

template <class T>
Vec<T> Mat<T>::get_col(size_t c) const{
	if(c < 0 || c >= cols_){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::get_col(size_t c1) const: Index exceeds matrix dimensions";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		Vec<T> result(rows_);
		size_t i;
		for(i = rows_; i--;){
			result.set(i,  data_[i][c]);
		}
		return result;
	}
}

template <class T>
Mat<T> Mat<T>::get_cols(size_t c1, size_t c2) const{
	if(c1 < 0 || c1 >= cols_ || c2 < 0 || c2 >= cols_ || c2 < c1){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::get_cols(size_t c1, size_t c2) const: Index exceeds matrix dimensions";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		Mat<T> result(rows_, c2 - c1 + 1);
		size_t i, j, rows = result.rows(), cols = result.cols();
		for(i = rows; i--;){
			for(j = cols; j--;){
				result.data_[i][j] = data_[i][j + c1];
			}
		}
		return result;
	}
}

template <class T>
Vec<T> Mat<T>::get_row(size_t r) const{
	if(r < 0 || r >= rows_){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::get_row(size_t r) const: Index exceeds matrix dimensions";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		Vec<T> result(cols_);
		size_t i;
		for(i = cols_; i--;){
			result.set(i, data_[r][i]);
		}
		return result;
	}
}

template <class T>
Mat<T> Mat<T>::get_rows(size_t r1, size_t r2) const{
	if(r1 < 0 || r1 >= rows_ || r2 < 0 || r2 >= rows_ || r2 < r1){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::get_rows(size_t r1, size_t r2) const: Index exceeds matrix dimensions";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		Mat<T> result(r2 - r1 + 1, cols_);
		size_t i, j, rows = result.rows(), cols = result.cols();
		for(i = rows; i--;){
			for(j = cols; j--;){
				result.data_[i][j] = data_[i + r1][j];
			}
		}
		return result;
	}
}

template <class T>
Mat<T> Mat<T>::get(size_t r1, size_t r2, size_t c1, size_t c2) const{
	if((*this).size() == 0){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator()(size_t r1, size_t r2, size_t c1, size_t c2): tried to access NULL MATRIX";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else if(r1 < 0 || r2 < 0 || r1 >= rows_ || r2 >= rows_ || r1 > r2 || c1 < 0 || c2 < 0 || c1 >= cols_ || c2 >= cols_ || c1 > c2){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::get_rows(size_t r1, size_t r2, size_t c1, size_t c2) const: Index exceeds matrix dimensions";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		Mat<T> result(r2 - r1 + 1, c2 - c1 + 1);
		size_t i, j, rows = result.rows(), cols = result.cols();
		for(i = rows; i--;){
			for(j = cols; j--;){
				result.data_[i][j] = data_[i + r1][j + c1];
			}
		}
		return result;
	}
}

// vec r is the vector with the row-indices and c is the vector with col-indices
template <class T>
Mat<T> Mat<T>::get(const Vec<T>& r, const Vec<T>& c) const{

	//Test if there are no elements in the mask vector
	if(r.size() == 0 || c.size() == 0){
		return mat(0,0);
	}

	//Test if size of index vectors exceed matrix dimensions
	if(r.size() > rows_ || c.size() > cols_){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::get(const vec& r, const vec& c): Index exceeds matrix dimensions";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}

	//Test if the elements of the  row-indices exceed bounds, and throw exception if they do.
	if(r.get(0) < 0 || r.get(0) >= rows_ ){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::get(const vec& r, const vec& c): Row-index exceeds matrix dimensions";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
	for(size_t i = 1; i < r.size(); i++){
		if( r.get(i) < 0 || r.get(i) >= rows_ || r.get(i) < r.get(i-1)){
			std::string msg = FILE_LINE_ERROR + " exception in  mat::get(const vec& r, const vec& c): Row-index exceeds matrix dimensions";
			log_error(msg.c_str());
			throw std::out_of_range(msg);
		}
	}

	//Test if the elements of the  col-indices exceed bounds, and throw exception if they do.
	if(c.get(0) < 0 || c.get(0) >= cols_ ){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::get(const vec& r, const vec& c): Col-index exceeds matrix dimensions";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
	for(size_t i = 1; i < c.size(); i++){
		if( c.get(i) < 0 || c.get(i) >= rows_ || c.get(i) < c.get(i-1)){
			std::string msg = FILE_LINE_ERROR + " exception in  mat::get(const vec& r, const vec& c1): Col-index exceeds matrix dimensions";
			log_error(msg.c_str());
			throw std::out_of_range(msg);
		}
	}

	Mat<T> result(r.size(), c.size());
	size_t i, j, rows = result.rows(), cols = result.cols();
	for(i = rows; i--;) {
		for(j = cols; j--;) {
			result.data_[i][j] = data_[r.get(i)][c.get(j)];
		}
	}
	return result;
}

template <class T>
void Mat<T>::zeros(){
	size_t i,j;
	for(i = rows_; i--;){
		for(j = cols_; j--;){
			data_[i][j] = T(0);
		}
	}
}

template <class T>
void Mat<T>::clear(){
	zeros();
}

template <class T>
void Mat<T>::ones(){
	size_t i,j;
	for(i = rows_; i--;){
		for(j = cols_; j--;){
			data_[i][j] = T(1);
		}
	}
}

// swap rows i and j
template <class T>
void Mat<T>::swap_rows(size_t i, size_t j) {
	if(i < rows_ && j < rows_ && i >= 0 && j >= 0){
		Vec<T> tmp = (*this).get_row(i);
		(*this).set_row(i, (*this).get_row(j));
		(*this).set_row(j, tmp);
	}else{
		std::string msg = FILE_LINE_ERROR + " exception in mat::swap_rows(size_t i, size_t j): out of range indices.";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
}

// swap columns i and j
template <class T>
void Mat<T>::swap_cols(size_t i, size_t j) {
	if(i < cols_ && j < cols_ && i >= 0 && j >= 0){
		Vec<T> tmp = (*this).get_col(i);
		(*this).set_col(i, (*this).get_col(j));
		(*this).set_col(j, tmp);
	}else{
		std::string msg = FILE_LINE_ERROR + " exception in mat::swap_cols(size_t i, size_t j): out of range indices.";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}
}

// ========= Overload basic operators ===========
template <class T>
void Mat<T>::operator=(const char* a){

	// Clear matrix from any previous values
	if( (*this).size() != 0 ){
		data_.resize(0, std::vector<double> (0));
		rows_ = 0;
		cols_ = 0;
	}

	//Convert input to string
	std::string str = std::string(a);

	// Remove brackets if there are any.
	if(str.size() >= 1){
		std::string start = str.substr(0,1);
		std::string end = str.substr(str.size()-1, str.size());
		if( start.compare("[") == 0 && end.compare("]") == 0){
			str.erase(0,1);
			str.erase(str.size()-1, str.size());
		}
	}

	// If there is no number in the string, consider it a null matrix
	std::string numbers = "0123456789";
	if( str.find_first_of(numbers) == std::string::npos ){
		return;
	}

	std::string delimiter = ";"; // DEFAULT DELIMITER FOR ROWS' SEPARATION
	size_t pos = 0;
	std::vector<std::string> row_vector;

	// Store the rows of the matrix in a string vector
	while ((pos = str.find(delimiter)) != std::string::npos) {
		row_vector.push_back(str.substr(0, pos));
		str.erase(0, pos + delimiter.length());
	}
	row_vector.push_back(str.substr(0, pos));

	// Find row dimension (it's equal to the number of the detected rows stored in the row_vector)
	size_t mat_row = row_vector.size();
	// Find column dimension
	Vec<T> first_row; first_row = row_vector[0].c_str();
	size_t mat_col = first_row.size();

	// Initialize matrix
	data_.resize(mat_row, std::vector<double> (mat_col));
	rows_ = mat_row;
	cols_ = mat_col;
	size_t i, j;
	for(i = rows_; i--;){
		for(j = cols_; j--;){
			data_[i][j] = 0;
		}
	}

	// Check that all the detected row vectors are of equal size,
	// throw exception otherwise.
	Vec<T> tmp;
	for(i = rows_; i--;){
		tmp = row_vector[i].c_str();
		if(tmp.size() != first_row.size()){
			std::string msg = FILE_LINE_ERROR + " exception in  mat::operator=(const char* a): rows must be of same length";
			log_error(msg.c_str());
			throw std::out_of_range(msg);
		}else{
			for(j = cols_; j--;){
				data_[i][j] = tmp[j];
			}
		}
	}
}

template <class T>
T& Mat<T>::operator()(size_t i, size_t j){
	if((*this).size() == 0){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator()(size_t i, size_t j): tried to access NULL MATRIX";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else if(i < 0 || i >= rows_ || j < 0 || j >= cols_){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator()(size_t i, size_t j): index out of range";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		return data_[i][j];
	}
}

template <class T>
Mat<T> Mat<T>::operator()(size_t r1, size_t r2, size_t c1, size_t c2){
	return (*this).get(r1, r2, c1, c2);
}

template <class T>
Mat<T> Mat<T>::operator+(const Mat<T>& m){
	if((*this).size() == 0 || m.size() == 0){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator+(const mat& m): tried to add NULL MATRIX";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else if(rows_ != m.rows() || cols_ != m.cols()){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator+(const mat& m): dimension mismatch";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		Mat<T> result = *this;
		size_t rows = result.rows(), cols = result.cols(), i = 0, j = 0;
		for(i = rows; i--;){
			for(j = cols; j--;){
				result.data_[i][j] += m.data_[i][j];
			}
		}
		return result;
	}
}

template <class T>
Mat<T> Mat<T>::operator+(T t){
	if((*this).size() == 0){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator+(double t): tried to add NULL MATRIX";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		Mat<T> result = *this;
		size_t rows = result.rows(), cols = result.cols(), i = 0, j = 0;
		for(i = rows; i--;){
			for(j = cols; j--;){
				result.data_[i][j] += t;
			}
		}
		return result;
	}
}

template <class T>
Mat<T> Mat<T>::operator-(const Mat<T>& m){
	if((*this).size() == 0 || m.size() == 0){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator-(const mat& m): tried to subtract NULL MATRIX";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else if(rows_ != m.rows() || cols_ != m.cols()){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator-(const mat& m): dimension mismatch";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		Mat<T> result = *this;
		size_t rows = result.rows(), cols = result.cols(), i = 0, j = 0;
		for(i = rows; i--;){
			for(j = cols; j--;){
				result.data_[i][j] -= m.data_[i][j];
			}
		}
		return result;
	}
}

template <class T>
Mat<T> Mat<T>::operator-(T t){
	if((*this).size() == 0){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator+(double t): tried to subtract NULL MATRIX";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		Mat<T> result = *this;
		size_t rows = result.rows(), cols = result.cols(), i = 0, j = 0;
		for(i = rows; i--;){
			for(j = cols; j--;){
				result.data_[i][j] -= t;
			}
		}
		return result;
	}
}


template <class T>
Mat<T> Mat<T>::operator*(T t){
	if((*this).size() == 0){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator*(double t): tried to multiply NULL MATRIX";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		Mat result = *this;
		size_t rows = result.rows(), cols = result.cols(), i = 0, j = 0;
		for(i = rows; i--;){
			for(j = cols; j--;){
				result.data_[i][j] *= t;
			}
		}
		return result;
	}
}

template <class T>
Mat<T> Mat<T>::operator*(const Mat<T>& m){
	if((*this).size() == 0 || m.size() == 0){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator*(const mat& m): tried to multiply NULL MATRIX";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else if(cols_ != m.rows()){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator*(const mat& m): dimension mismatch";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		mat result(rows_, m.cols());
		size_t rows = result.rows(), cols = result.cols(), i = 0, j = 0, k = 0;
		size_t common_dimension = cols_;
		double tmp = 0;
		for(i = rows; i--;){
			for(k = common_dimension; k--;){
				tmp = data_[i][k];
				for(j = cols; j--;){
					result.data_[i][j] += tmp * m.data_[k][j];
				}
			}
		}
		return result;
	}
}

template <class T>
Vec<T> Mat<T>::operator*(const Vec<T>& v){
	//convert vec into an v.size()x1 matrix depending
	// a similar multiplication should exist in vec class
	//
	if((*this).size() == 0 || v.size() == 0){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator*(const vec& v): tried to multiply NULL MATRIX or NULL VECTOR";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else if(cols_ != v.size()){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator*(const vec& v): dimension mismatch";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		Vec<T> result(rows_);
		size_t i, j;
		for(i = rows_; i--;){
			for(j = cols_; j--;){
				result.set(i, result.get(i) + data_[i][j]*v.get(j));
			}
		}
		return result;
	}
}

template <class T>
Mat<T> Mat<T>::operator/(T t){
	if((*this).size() == 0){
		std::string msg = FILE_LINE_ERROR + " exception in  mat::operator*(double t): tried to divide NULL MATRIX";
		log_error(msg.c_str());
		throw std::out_of_range(msg);
	}else{
		Mat result = *this;
		size_t rows = result.rows(), cols = result.cols(), i = 0, j = 0;
		for(i = rows; i--;){
			for(j = cols; j--;){
				result.data_[i][j] /= t;
			}
		}
		return result;
	}
}

template <class T>
void Mat<T>::print(int precision){

	if((*this).size() == 0){
		printf("| | \n");
		return;
	}
	if(typeid(data_[0][0]) == typeid(double)){
		for(size_t i = 0; i < rows_; i++){
			printf("| ");
			for(size_t j = 0; j < cols_; j++){
				printf("%.*f ", precision, data_[i][j]);
			}
			printf("|\n");
		}
	}else if(typeid(data_[0][0]) == typeid(int)){
		for(size_t i = 0; i < rows_; i++){
			printf("| ");
			for(size_t j = 0; j < cols_; j++){
				printf("%d ", data_[i][j]);
			}
			printf("|\n");
		}
	}
}

template <class T>
void Mat<T>::print(){

	if((*this).size() == 0){
		printf("| | \n");
		return;
	}
	if(typeid(data_[0][0]) == typeid(double)){
		for(size_t i = 0; i < rows_; i++){
			printf("| ");
			for(size_t j = 0; j < cols_; j++){
				printf("%.3f ", data_[i][j]);
			}
			printf("|\n");
		}
	}else if(typeid(data_[0][0]) == typeid(int)){
		for(size_t i = 0; i < rows_; i++){
			printf("| ");
			for(size_t j = 0; j < cols_; j++){
				printf("%d ", data_[i][j]);
			}
			printf("|\n");
		}
	}
}

// ***************** DEFINITION OF FRIEND FUNCTIONS ********************************************

// Implementation of the algorithm described in https://en.wikipedia.org/wiki/Strassen_algorithm
// This algorithm is beneficial for large square matrices only.
// ATTENTION: the Strassen multiplication algorithm requires more memory than the traditional one.

template <class T>
Mat<T> strassen_algorithm(const Mat<T> &a, const Mat<T> &b, size_t leafsize ) {

	size_t size = a.rows();

	if (size <= leafsize) {
		Mat<T> c(a.rows(), b.cols());
		c = a;
		return c * b;
	}else{

		size_t new_size = size/2;
		Vec<T> inner(new_size);
		Mat<T>
		a11(new_size, new_size), a12(new_size, new_size), a21(new_size, new_size), a22(new_size, new_size),
		b11(new_size, new_size), b12(new_size, new_size), b21(new_size, new_size), b22(new_size, new_size),
		c11(new_size, new_size), c12(new_size, new_size), c21(new_size, new_size), c22(new_size, new_size),
		m1(new_size, new_size), m2(new_size, new_size), m3(new_size, new_size), m4(new_size, new_size),
		m5(new_size, new_size), m6(new_size, new_size), m7(new_size, new_size),
		tmp_a(new_size, new_size), tmp_b(new_size, new_size);

		size_t i, j;

		//dividing the matrices in 4 sub-matrices:
		for (i = new_size; i--;) {
			for (j = new_size; j--;) {
				a11.data_[i][j] = a.data_[i][j];
				a12.data_[i][j] = a.data_[i][j + new_size];
				a21.data_[i][j] = a.data_[i + new_size][j];
				a22.data_[i][j] = a.data_[i + new_size][j + new_size];

				b11.data_[i][j] = b.data_[i][j];
				b12.data_[i][j] = b.data_[i][j + new_size];
				b21.data_[i][j] = b.data_[i + new_size][j];
				b22.data_[i][j] = b.data_[i + new_size][j + new_size];
			}
		}

		// Calculating m1 to m7:

		tmp_a = a11 + a22; // a11 + a22
		tmp_b = b11 + b22; // b11 + b22
		m1 = strassen_algorithm( tmp_a, tmp_b, leafsize ); // m1 = (a11+a22) * (b11+b22)

		tmp_a = a21 + a22; // a21 + a22
		m2 = strassen_algorithm( tmp_a, b11, leafsize ); // m2 = (a21+a22) * (b11)

		tmp_b = b12 - b22; // b12 - b22
		m3 = strassen_algorithm( a11, tmp_b, leafsize ); // m3 = (a11) * (b12 - b22)

		tmp_b = b21 - b11; // b21 - b11
		m4 = strassen_algorithm( a22, tmp_b, leafsize ); // m4 = (a22) * (b21 - b11)

		tmp_a = a11 + a12; // a11 + a12
		m5 = strassen_algorithm( tmp_a, b22, leafsize ); // m5 = (a11+a12) * (b22)

		tmp_a = a21 - a11; // a21 - a11
		tmp_b = b11 + b12; // b11 + b12
		m6 = strassen_algorithm( tmp_a, tmp_b, leafsize ); // m6 = (a21-a11) * (b11+b12)

		tmp_a = a12 - a22; // a12 - a22
		tmp_b = b21 + b22; // b21 + b22
		m7 = strassen_algorithm( tmp_a, tmp_b, leafsize ); // m7 = (a12-a22) * (b21+b22)

		// Calculate c21, c21, c11 and c22:

		c11 = m1 + m4 - m5 + m7;
		c12 = m3 + m5;
		c21 = m2 + m4;
		c22 = m1 - m2 + m3 + m6;

		// Shape the C matrix from its calculated 4 submatrices:
		Mat<T> c(a.rows(), b.cols());
		for (i = new_size; i--;) {
			for (j = new_size; j--;) {
				c.data_[i][j] = c11.data_[i][j];
				c.data_[i][j + new_size] = c12.data_[i][j];
				c.data_[i + new_size][j] = c21.data_[i][j];
				c.data_[i + new_size][j + new_size] = c22.data_[i][j];
			}
		}
		return c;
	}
}

template <class T>
Mat<T> strassen(const Mat<T> &a, const Mat<T> &b) {

	if(a.rows() != a.cols() || a.rows() != b.cols() || b.rows() != b.cols()){
		std::string msg = FILE_LINE_ERROR + " exception in strassen(const mat &a, const mat &b): NON-SQUARE MATRICES";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{

		size_t n = a.rows();

		// If the matrices A, B are not of type 2^n × 2^n,
		// we fill the missing rows and columns with zeros.
		// Thus, we create new matrices with size equal to
		// the next closest power of two.
		size_t m = (size_t) std::pow(2, (ceil(std::log2(n))));
		Mat<T> a_new(m, m), b_new(m, m), c_new(m, m);

		// Calculate leafsize. How deep the strassen algorithm will recurse
		size_t leafsize = ceil((m/32.0));;

		size_t i,j;
		// copy the elements of the small matrices
		// into the larger ones.
		for(i = n; i--;) {
			for (j = n; j--;) {
				a_new.data_[i][j] = a.data_[i][j];
				b_new.data_[i][j] = b.data_[i][j];
			}
		}

		c_new = strassen_algorithm(a_new, b_new, leafsize);
		Mat<T> c(a.rows(), a.cols());
		// pass the elements of the large matrix
		// into the smaller matrix.
		for(i = n; i--;) {
			for (j = n; j--;) {
				c.data_[i][j] = c_new.data_[i][j];
			}
		}
		return c;
	}
}

template <class T>
Mat<T> transpose(const Mat<T>& m){
	Mat<T> result(m.cols(), m.rows());
	size_t i, j, rows = result.rows(), cols = result.cols();
	for(i = rows; i--;){
		for(j = cols; j--;){
			result.data_[i][j] = m.data_[j][i];
		}
	}
	return result;
}

// LU-Decomposition taken from: https://en.wikipedia.org/wiki/LU_decomposition
template <class T>
Vec<T> lup_decompose(Mat<T>& a, bool& is_singular) {

	size_t n = a.rows();
	Vec<T> pivot(n + 1); //Unit permutation vector
	size_t i, j, k, imax;
	double maxA, absA;
	Vec<T> row_pivot(n);

	for (i = 0; i <= n; i++)
		pivot[i] = i; // pivot[n] initialized with n

	for (i = 0; i < n; i++) {
		maxA = 0.0;
		imax = i;

		for (k = i; k < n; k++)
			if ((absA = fabs(a.data_[k][i])) > maxA) {
				maxA = absA;
				imax = k;
			}

		if (maxA < SINGULARITY_THRESHOLD){
			is_singular = true;
		}

		if (imax != i) {
			//pivoting P
			j = pivot[i];
			pivot[i] = pivot[imax];
			pivot[imax] = j;

			//pivoting rows of A
			row_pivot = a.get_row(i);
			a.set_row(i, a.get_row(imax));
			a.set_row(imax,  row_pivot);

			//counting pivots starting from N (for determinant)
			pivot[n]++;
		}

		for (j = i + 1; j < n; j++) {
			a.data_[j][i] /= a.data_[i][i];
			for (k = i + 1; k < n; k++){
				a.data_[j][k] -= a.data_[j][i] * a.data_[i][k];
			}
		}
	}
	return pivot;

	//decomposition done
}

template <class T>
Mat<T> lup_invert(Mat<T>& a, const Vec<T>& pivot) {

	size_t N = a.rows();
	Mat<T> a_inv(N, N);
	for (size_t j = 0; j < N; j++) {
		for (size_t i = 0; i < N; i++) {
			if (pivot.data_[i] == j)
				a_inv.data_[i][j] = 1.0;
			else
				a_inv.data_[i][j] = 0.0;

			for (size_t k = 0; k < i; k++)
				a_inv.data_[i][j] -= a.data_[i][k] * a_inv.data_[k][j];
		}
		// this was initially the first loop
		// for (int i =  N - 1; i >= 0; i--)
		// But with unsigned integers like size_t there is
		// no point expecting i to get less than 0 since it
		// will wrap around and jump to the size_t_MAX
		for (size_t i =  N - 1; i < SIZE_T_MAX; i--) {
			for (size_t k = i + 1; k < N; k++)
				a_inv.data_[i][j] -= a.data_[i][k] * a_inv.data_[k][j];

			a_inv.data_[i][j] /= a.data_[i][i];
		}
	}
	return a_inv;
}

template <class T>
Mat<T> inv(const Mat<T>& a){

	if(a.rows() == 1 && a.cols() == 1){
		Mat<T> a_inv(1,1);
		a_inv(0,0) = 1.0/a.data_[0][0];
		return a_inv;
	}else if(a.rows() == a.cols()){
		// copy 'a' in a temporary matrix 'tmp'
		Mat<T> tmp = a;
		bool is_singular = false;
		Vec<T> pivot = lup_decompose(tmp, is_singular);
		if(!is_singular){
			Mat<T> a_inv = lup_invert(tmp, pivot);
			return a_inv;
		}else{
			std::string msg = FILE_LINE_ERROR + "warning in  mat::inv(const mat& m): SINGULAR MATRIX.";
			warning(msg.c_str());
			return abs(a)*NaN_NUMERIC_LIMIT(T);
		}
	}else{
		std::string msg = FILE_LINE_ERROR + " exception in  mat::inv(const mat& m): NON-SQUARE MATRIX: use pinv(const mat& m)";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
}

// PseudoInverse: calculates the inverse of a non-square matrix
// Based on Moore–Penrose pseudoinverse
//
// Some info for intuition:
// Thin matrix case A(m,n) with m > n:
// a) rank(A) <= n
// b) The homogeneous system of equations Ax = 0 has a unique (trivial) solution
//    if and only if rank(A) = n.
// c) the system of equations Ax = b is inconsistent for all b in R^m
// d) the system of equations Ax = b has at most one solution for every b in R^m
//    if and only if rank(A) = n.
//
// Something similar holds for fat matrix case A(m,n) with m < n
template <class T>
Mat<T> pinv(Mat<T>& a){

	Mat<T> a_inv;

	if( a.rows() == a.cols() ){
		a_inv = inv(a); // Normal inverse
	}else if( std::abs(determinant(transpose(a)*a)) > std::abs(determinant(a*transpose(a))) ){
		a_inv = inv(transpose(a)*a)*transpose(a); // Left inverse (full column rank)
	}else if( std::abs(determinant(transpose(a)*a)) < std::abs(determinant(a*transpose(a))) ){
		a_inv = transpose(a)*inv(a*transpose(a)); // Right inverse (full row rank)
	}else{
		std::string msg = FILE_LINE_ERROR + "warning in 'mat::pinv(const mat& a)': ILL-DEFINED MATRIX !!!";
		warning(msg.c_str());
		Mat<T> tmp(a.cols(), a.rows());
		a_inv = abs(tmp)*NaN_NUMERIC_LIMIT(T);
	}

	return a_inv;
}




// ##################################################################################################
// ############################ MISCELLANEOUS OPERATIONS AND FUNCTIONS ##############################

/* zeros(), ones(), eye(), concat_hor(), concat_ver(), outer_prod(), diag(), inv(), eig(), det() */

// it returns the vectorized version of the matrix
template <class T>
inline Vec<T> mat2vec(const Mat<T>& m){
	if(m.size() == 0){
		return Vec<T>(0);
	}else{
		Vec<T> result;
		size_t i, size = m.rows();
		for(i = 0; i < size; i++){
			result = concat(result, m.get_row(i));
		}
		return result;
	}
}

template <class T>
inline T max(const Mat<T>& m){
	return max(mat2vec(m));
}

template <class T>
inline T min(const Mat<T>& m){
	return min(mat2vec(m));
}

template <class T>
inline Mat<T> abs(const Mat<T>& m){
	Mat<T> result(m.rows(), m.cols());
	size_t i, j, rows = m.rows(), cols = m.cols();
	for(i = rows; i--; ){
		for(j = cols; j--;){
			result.set(i, j, std::abs(m.get(i,j)));
		}
	}
	return result;
}

// This will return a mask with ones where the input matrix
// has non-zero elements and zero elsewhere.
template <class T>
inline Mat<T> find_non_zero(const Mat<T>& m){
	size_t i, j, rows = m.rows(), cols = m.cols();
	Mat<T> result(rows, cols);
	for(i = rows; i--;){
		for(j = cols; j--;){
			if(m.get(i,j) != 0){
				result.set(i, j, 1.0);
			}else{
				result.set(i, j, 0);
			}
		}
	}
	return result;
}

// This will return a mask with ones where the input
// matrix has zero elements and zero elsewhere.
template <class T>
inline Mat<T> find_zero(const Mat<T>& m){
	size_t i, j, rows = m.rows(), cols = m.cols();
	Mat<T> result(rows, cols);
	for(i = rows; i--;){
		for(j = cols; j--;){
			if(m.get(i,j) == 0){
				result.set(i, j, 1.0);
			}else{
				result.set(i, j, 0);
			}
		}
	}
	return result;
}

inline mat rand(size_t m, size_t n){
	if( (n*m) > MAX_ACCEPTABLE_VECTOR_SIZE*MAX_ACCEPTABLE_VECTOR_SIZE ){
		std::string msg = FILE_LINE_ERROR + " exception in mat::rand_double(size_t n, size_t m): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{
		mat a(m, n);
		size_t i = 0, j = 0;

		// C++11 feature: It will be used to obtain a seed for the random number engine
		std::random_device rd;
		// C++11 feature: Standard mersenne_twister_engine seeded with rd()
		std::mt19937 gen(rd());

		// Define range
		double min = -10, max = 10;
		std::uniform_real_distribution<double> dis(0, 2*max);
		for(i = m; i--;){
			for(j = n; j--;){
				// Generate random double number within the range [-10 10]
				a(i, j) = min + (double) dis(gen);
			}
		}
		return a;
	}
}

inline imat rand_i(size_t m, size_t n){
	if( (n*m) > MAX_ACCEPTABLE_VECTOR_SIZE*MAX_ACCEPTABLE_VECTOR_SIZE ){
		std::string msg = FILE_LINE_ERROR + " exception in mat::rand_double(size_t n, size_t m): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{
		imat a(m, n);
		size_t i = 0, j = 0;

		// C++11 feature: It will be used to obtain a seed for the random number engine
		std::random_device rd;
		// C++11 feature: Standard mersenne_twister_engine seeded with rd()
		std::mt19937 gen(rd());

		// Define range
		int min = -10, max = 10;

		std::uniform_int_distribution<int> dis(0, 2*max);
		for(i = m; i--;){
			for(j = n; j--;){
				// Generate random double number within the range [-10 10]
				a(i, j) = min + (int) dis(gen);
			}
		}
		return a;
	}
}

inline mat rand_symmetric(size_t n){

	if( n > MAX_ACCEPTABLE_VECTOR_SIZE ){
		std::string msg = FILE_LINE_ERROR + " exception in rand_int_symmetric(size_t n): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{

		mat result(n, n);
		size_t i = 0, j = 0;

		// C++11 feature: It will be used to obtain a seed for the random number engine
		std::random_device rd;
		// C++11 feature: Standard mersenne_twister_engine seeded with rd()
		std::mt19937 gen(rd());

		// Define range
		double min = -10, max = 10;
		std::uniform_real_distribution<double> dis(0, 2*max);

		for(i = n; i--;){
			for(j = n; j--;){
				// Generate random double number within the range [-10 10]
				result(i, j) = min + (double) dis(gen);
			}
		}

		for(i = n; i--;){
			for(j = n; j--;){
				if(i < j){
					result(i,j) = result(j,i);
				}
			}
		}
		return result;
	}
}

inline imat rand_symmetric_i(size_t n){

	if( n > MAX_ACCEPTABLE_VECTOR_SIZE ){
		std::string msg = FILE_LINE_ERROR + " exception in rand_int_symmetric(size_t n): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{

		imat result(n, n);
		size_t i = 0, j = 0;

		// C++11 feature: It will be used to obtain a seed for the random number engine
		std::random_device rd;
		// C++11 feature: Standard mersenne_twister_engine seeded with rd()
		std::mt19937 gen(rd());

		// Define range
		int min = -10, max = 10;
		std::uniform_int_distribution<int> dis(0, 2*max);

		for(i = n; i--;){
			for(j = n; j--;){
				// Generate random double number within the range [-10 10]
				result(i, j) = min + (int) dis(gen);
			}
		}

		for(i = n; i--;){
			for(j = n; j--;){
				if(i < j){
					result(i,j) = result(j,i);
				}
			}
		}
		return result;
	}
}

inline mat zeros(size_t n, size_t m){
	if( (n*m) > (MAX_ACCEPTABLE_VECTOR_SIZE*MAX_ACCEPTABLE_VECTOR_SIZE) ){
		std::string msg = FILE_LINE_ERROR + " exception in mat::zeros(size_t n, size_t m): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{
		mat a(n,m);
		size_t i, j;
		for(i = n; i--;){
			for(j = m; j--;){
				a.set(i, j, 0.0);
			}
		}
		return a;
	}
}

inline imat zeros_i(size_t n, size_t m){
	if( (n*m) > (MAX_ACCEPTABLE_VECTOR_SIZE*MAX_ACCEPTABLE_VECTOR_SIZE) ){
		std::string msg = FILE_LINE_ERROR + " exception in mat::zeros(size_t n, size_t m): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{
		imat a(n,m);
		size_t i, j;
		for(i = n; i--;){
			for(j = m; j--;){
				a.set(i, j, 0);
			}
		}
		return a;
	}
}

inline mat ones(size_t n, size_t m){
	if( (n*m) > MAX_ACCEPTABLE_VECTOR_SIZE*MAX_ACCEPTABLE_VECTOR_SIZE ){
		std::string msg = FILE_LINE_ERROR + " exception in mat::ones(size_t n, size_t m): n*m should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{
		mat a(n,m);
		size_t i, j;
		for(i = n; i--;){
			for(j = m; j--;){
				a.set(i, j, 1.0);
			}
		}
		return a;
	}
}

inline imat ones_i(size_t n, size_t m){
	if( (n*m) > MAX_ACCEPTABLE_VECTOR_SIZE*MAX_ACCEPTABLE_VECTOR_SIZE ){
		std::string msg = FILE_LINE_ERROR + " exception in mat::ones(size_t n, size_t m): n*m should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{
		imat a(n,m);
		size_t i, j;
		for(i = n; i--;){
			for(j = m; j--;){
				a.set(i, j, 1);
			}
		}
		return a;
	}
}

inline mat eye(size_t k){
	if( k > MAX_ACCEPTABLE_VECTOR_SIZE){
		std::string msg = FILE_LINE_ERROR + " exception in mat::eye(size_t k): k should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{
		mat result(k,k);
		size_t i;
		for(i = k; i--;){
			result(i,i) = 1.0;
		}
		return result;
	}
}

inline imat eye_i(size_t k){
	if( k > MAX_ACCEPTABLE_VECTOR_SIZE){
		std::string msg = FILE_LINE_ERROR + " exception in mat::eye(size_t k): k should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{
		imat result(k,k);
		size_t i;
		for(i = k; i--;){
			result(i,i) = 1;
		}
		return result;
	}
}

// It returns a vector containing the diagonal elements of a matrix
template <class T>
inline Vec<T> diag(const Mat<T>& m){
	if(m.rows() != m.cols()){
		std::string msg = FILE_LINE_ERROR + " exception in mat::diag(const mat& m)): diagonal is defined only for square matrices";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
	Vec<T> result(m.rows());
	size_t i, size = result.size();
	for(i = size; i--;){
		result[i] = m.get(i,i);
	}
	return result;
}

// It returns a matrix with diagonal elements the elements of the input-vector and zero all the rest
template <class T>
inline Mat<T> diag(const Vec<T>& v){
	Mat<T> result(v.size(), v.size());
	size_t i, size = result.rows();
	for(i = size; i--;){
		result(i,i) = v.get(i);
	}
	return result;
}

template <class T>
inline Mat<T> concat_hor(const Mat<T>& m1, const Mat<T>& m2){
	if(m1.rows() != m2.rows()){
		std::string msg = FILE_LINE_ERROR + " exception in mat::concat_hor(const mat& m1, const mat& m2): dimension mismatch";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{
		Mat<T> result(m1.rows(), m1.cols() + m2.cols());
		size_t i, j, rows = result.rows(), cols = result.cols(), m1_cols = m1.cols();
		for(i = rows; i--;){
			for(j = 0; j < m1_cols; j++){
				result(i,j) = m1.get(i,j);
			}
			for(j = m1.cols(); j < cols; j++){
				result(i,j) = m2.get(i,j - m1.cols());
			}
		}
		return result;
	}
}

template <class T>
inline Mat<T> concat_ver(const Mat<T>& m1, const Mat<T>& m2){
	if(m1.cols() != m2.cols()){
		std::string msg = FILE_LINE_ERROR + " exception in mat::concat_ver(const mat& m1, const mat& m2): dimension mismatch";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{
		Mat<T> result(m1.rows() + m2.rows(), m1.cols());
		size_t i, j, rows = result.rows(), cols = result.cols(), m1_rows = m1.rows();
		for(j = cols; j--;){
			for(i = 0; i < m1_rows; i++){
				result(i,j) = m1.get(i,j);
			}
			for(i = m1_rows; i < rows; i++){
				result(i,j) = m2.get(i - m1_rows,j);
			}
		}
		return result;
	}
}

template <class T>
inline Mat<T> outer_product(const Vec<T>& v1, const Vec<T>& v2){
	if(v1.size() == 0 || v2.size() == 0){
		std::string msg = FILE_LINE_ERROR + " exception in mat::outer_product(const vec& v1, const vec& v2): NULL VECTOR";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{
		Mat<T> result(v1.size(), v2.size());
		size_t i, j, rows = result.rows(), cols = result.cols();

		for(i = rows; i--;){
			for(j = cols; j--;){
				result.set(i, j, v1.get(i)* v2.get(j));
			}
		}
		return result;
	}
}

template <class T>
inline T determinant(const Mat<T>& m){

	if( m.size() == 0 ){
		std::string msg = FILE_LINE_ERROR + " exception in mat::determinant(const mat& m): Not defined for NULL MATRIX";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{
		T det = 1;
		if(m.rows() == m.cols()){
			size_t size = m.rows();
			Mat<T> tmp(m.rows(), m.cols());
			for (size_t i = 0; i < m.rows(); i++){
				for (size_t j = 0; j < m.cols(); j++){
					tmp(i,j) = m.get(i,j);
				}
			}
			double k = 1, con, kc;
			size_t p;
			for(p = 0; p < size - 1; p++){
				if(tmp(p,p) != 0)//prepares the pivots by making them 1
				{
					k *= (1 /tmp.get(p,p));
					kc = (1 / tmp.get(p,p));
					for(size_t j = p; j < size; j++){
						tmp(j,p) *= kc;
					}
				}
				for(size_t c = p + 1; c < size; c++)//makes the triangle matrix
				{                             //with zeros above the diagonal
					con = -1 * tmp.get(p,c);
					for(size_t i = 0; i < size; i++){
						tmp(i,c) += (tmp.get(i,p) * con);
					}
				}
			}
			if(tmp(p,p) != 0)// makes the elemnt n,n 1 to end the pivots
			{
				k *= (1 / tmp.get(p,p));
				for(size_t j = p; j < size; j++)
				{
					tmp(p,j) /= tmp.get(p,p);
				}
			}
			for (size_t z = 0; z < size; z++)
			{
				det *= tmp.get(z,z);
			}
			det /= k;  //divides by the connstant acumulated so far
			//where det|A| = det|B|/k where B is the triangular matrix
		}
		return det;
	}
}


/*
   In any magic square, the first number i.e. 1 is stored at position (n/2, n-1).
   Let this position be (i,j). The next number is stored at position (i-1, j+1)
   where we can consider each row & column as circular array i.e. they wrap around.

   Three conditions hold:

    1. The position of next number is calculated by decrementing row number of previous
       number by 1, and incrementing the column number of previous number by 1. At any time,
       if the calculated row position becomes -1, it will wrap around to n-1. Similarly, if
       the calculated column position becomes n, it will wrap around to 0.

    2. If the magic square already contains a number at the calculated position, calculated
       column position will be decremented by 2, and calculated row position will be incremented by 1.

    3. If the calculated row position is -1 & calculated column position is n, the new position
       would be: (0, n-2).
 *
 */

template <class T>
inline Mat<T> magic_square(int n){

	if(n < 0 || n > MAX_ACCEPTABLE_VECTOR_SIZE){
		std::string msg = FILE_LINE_ERROR + "exception in magic_square(int n): n should lie in [0," + std::to_string(MAX_ACCEPTABLE_VECTOR_SIZE) +"]";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else if((n & 1) == 0){
		std::string msg = FILE_LINE_ERROR + "exception in magic_square(int n): n should be an odd number";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}else{
		Mat<T> result(n,n);

		// Initialize position for 1
		int i = n/2, j = n - 1, num;

		// One by one put all values in magic square
		for (num = 1; num <= n*n; ){
			if (i == -1 && j == n){  //3rd condition
				j = n-2;
				i = 0;
			}else{
				//1st condition helper if next number
				// goes to out of square's right side
				if (j == n){
					j = 0;
				}
				//1st condition helper if next number
				// is goes to out of square's upper side
				if (i < 0){
					i=n-1;
				}
			}

			if (result(i, j) != 0){ //2nd condition
				j -= 2;
				i++;
				continue;
			}else{
				result(i, j) = num++; //set number
			}
			j++; i--; // 1st condition
		}
		return result;
	}
}


} /* namespace algebra */

#endif /* MAT_H_ */
