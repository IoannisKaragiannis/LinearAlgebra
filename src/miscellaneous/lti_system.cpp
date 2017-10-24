/*================================================================================
 * Name         : lti_system.cpp implements a discrete linear time-invariant system.
 * Version      : 1.0.0, 15 Oct 2017
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
==================================================================================*/

#include "../../include/miscellaneous/lti_system.h"

namespace algebra {

lti_system::lti_system() {
	x = zeros(0,0);
	F = zeros(0,0);
	B = zeros(0,0);
	Q = zeros(0,0);

	z = zeros(0,0);
	H = zeros(0,0);
	R = zeros(0,0);
	dt = 0;
	initial_conditions = false;
}

lti_system::~lti_system() {
	// TODO Auto-generated destructor stub
}


mat lti_system::get_state_transition_matrix(){
	return F;
}

mat lti_system::get_control_matrix(){
	return B;
}


mat lti_system::get_process_noise_variance(){
	return Q;
}

mat lti_system::get_observation_matrix(){
	return H;
}

mat lti_system::get_observation_noise_variance(){
	return R;
}

double lti_system::get_sampling_period(){
	return dt;
}


// The order of how you set your lti system is very important.
// If you're not sure, write down the equations describing your
// system and you will figure out what the right dimensions should be.
void lti_system::set_system(const mat& m1, const mat& m2, const mat& m3,
		const mat& m4, const mat& m5, double sampling_period){
	F = m1; B = m2; Q = m3;
	H = m4; R = m5;


	if(sampling_period > 0 && std::abs(sampling_period) != Inf){
		dt = sampling_period;
	}else{
		std::string msg = FILE_LINE_ERROR + "exception in lti_system::set_system(...): frequency has to be positive";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}

	check_dimension_mismatch();

	u.set_size(B.cols(), 1);
	z.set_size(H.rows(), 1);

}

void lti_system::check_dimension_mismatch(){

	// =====  SYSTEM DYNAMIC MODEL ========

	// Check that transition matrix F is square
	if( F.rows() != F.cols() ){
		std::string msg = FILE_LINE_ERROR + "exception in lti_system::set_system(...): F has to be square";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}

	// Check that input matrix has correct columns
	if( B.rows() != F.rows() ){
		std::string msg = FILE_LINE_ERROR + "exception in lti_system::set_system(...): B.rows() != F.rows()";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}

	// Check that Q is square and has dimensions equal to F
	if( Q.rows() != Q.cols() || Q.rows() != F.rows() ){
		std::string msg = FILE_LINE_ERROR + "exception in lti_system::set_system(...): Q.rows() != Q.cols() || Q.rows() ! F.rows()";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}

	// =====  OBSERVATION MODEL ========

	// check that the observation matrix has the right dimensions
	// Usually in most systems we can only observe a number of states
	// from the state vector; not all of them.
	if(H.cols() != F.rows() || H.rows() > F.rows()){
		std::string msg = FILE_LINE_ERROR + "exception in lti_system::set_system(...): H.cols() != F.rows() || H.rows() > F.rows()";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}

	// Check observation noise matrix. It has to be square
	if(R.rows() != R.cols() || R.rows() != H.rows()){
		std::string msg = FILE_LINE_ERROR + "exception in lti_system::set_system(...): R.rows() != R.cols() || R.rows() != H.rows()";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
}


void lti_system::run_model(const vec& x0, const vec& input){

	if(!initial_conditions){
		set_initial_conditions(x0);
	}

	if( input.size() == u.rows()){
		u.set_col(0, input);
		x = F*x + B*u;
		z = H*x;
	}else{
		std::string msg = FILE_LINE_ERROR + "exception in lti_system::run_deterministic_model(const vec& input): Erroneous input dimension ";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
}

void lti_system::set_initial_conditions(const vec& x0){
	//Set initial conditions of the state
	if( x0.size() == F.rows() ){
		x.set_size(x0.size(), 1);
		x.set_col(0, x0);
		initial_conditions = true;
	}else{
		std::string msg = FILE_LINE_ERROR + "exception in lti_system::set_system(...): erroneous dimensions of initial conditions x0";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
}


vec lti_system::get_state(){
	return mat2vec(x);
}

vec lti_system::get_output(){
	return mat2vec(z);
}


} /* namespace algebra */
