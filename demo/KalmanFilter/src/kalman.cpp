/*====================================================================================================
 * Name         : kalman.cpp implements an online Kalman filter.
 * Version      : 1.0.0, 15 Oct 2017
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
=====================================================================================================*/

#include "../include/kalman.h"

namespace algebra {

kalman::kalman() {
	x_hat.set_size(0,0);
	u.set_size(0,0);
	F.set_size(0,0);
	B.set_size(0,0);
	Q.set_size(0,0);
	P.set_size(0, 0);

	z.set_size(0,0);
	H.set_size(0,0);
	R.set_size(0,0);
	I = eye(0);
	K.set_size(0,0);
	initialize = false;
	initial_conditions = false;

}

kalman::~kalman() {
	// TODO Auto-generated destructor stub
}

void kalman::set_initial_conditions(const vec& x0, const mat& P0){
	x_hat.set_size(x0.size(), 1);
	x_hat.set_col(0, x0);
	P.set_size(x0.size(), x0.size());
	P = P0;
	initial_conditions = true;
}

void kalman::update( lti_system &sys, const vec& input, const vec& measurement ){

	// Initialize filter
	if(!initialize){
		initialize_filter(sys);
	}

	// Set default initial conditions x_hat[0] = 0, P_hat[0] = 0
	// Without them the algorithm cannot start
	if(!initial_conditions){
		size_t size = sys.get_state_transition_matrix().rows();
		// Set default initial conditions to zero
		vec x0 = zeros(size);
		mat P0 = eye(size);
		set_initial_conditions(x0, P0);
	}

	// update input(u) and measurement(z) vectors
	update_input_and_measurement(input, measurement);

	// =======  TIME UPDATE ============

	// Update the state estimate (x[k|k-1] = F*x[k-1|k-1] + B*u[k-1]).
	x_hat = F * x_hat + B * u;

	// Time update covariance (P[k|k-1] = F*P[k-1|k-1]*F' + Q)
	P = F * P * transpose(F) + Q;\

	// ========== MEASUREMENT UPDATE =======

	// Kalman gain calculation K = (P*H')*inv(H*P*H' + R)
	K = P*transpose(H)*inv(H * P * transpose(H) + R);

	// Calculate the filter estimate x(k|k) based on innovation
	x_hat = x_hat + K*(z - H * x_hat);

	// Update covariance matrix P(k|k)
	P = (I - K * H) * P;

	// For numerical stability (the covariance matrix has to be symmetric)
	P = (P + transpose(P)) * 0.5;
}

vec kalman::get_estimate(){
	if(x_hat.cols() == 1){
		return mat2vec(x_hat);
	}else{
		std::string msg = FILE_LINE_ERROR + "exception in kalman::get_estimate(): x_hat must be a vector => x_hat.cols() = 1 ";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
}

vec kalman::get_cov_error(){
	return diag(P);
}

void kalman::update_input_and_measurement(const vec& input, const vec& measurement){
	if( input.size() == u.rows()){
		u.set_col(0, input);
	}else{
		std::string msg = FILE_LINE_ERROR + "exception in kalman::update(...): Erroneous input dimension ";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}

	if(measurement.size() == z.rows()){
		z.set_col(0, measurement);
	}else{
		std::string msg = FILE_LINE_ERROR + "exception in kalman::update(...): Erroneous measurement dimension ";
		log_error(msg.c_str());
		throw std::invalid_argument(msg);
	}
}


void kalman::initialize_filter(lti_system& sys){

	F = sys.get_state_transition_matrix();
	B = sys.get_control_matrix();
	Q = sys.get_process_noise_variance();
	u.set_size(B.cols(), 1);

	H = sys.get_observation_matrix();
	R = sys.get_observation_noise_variance();
	z.set_size(H.rows(), 1);

	K.set_size(B.rows(), B.cols());

	I = eye(F.rows());

	initialize = true;

}

} /* namespace algebra */
