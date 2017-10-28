/*===========================================================================
 * Name         : KalmanFilter.cpp implements an online Kalman filter
 *                for the simple case of a free falling ball.
 * Version      : 1.0.0, 15 Oct 2017
 *
 * Copyright (c) 2017 Ioannis Karagiannis
 * All rights reserved

 * This file is part of the LinearAlgebra library.

 * KalmanFilter is free software: you can redistribute it and/or modify
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
 *               ioanniskaragiannis1987@gmail.com
===============================================================================*/


#include <iostream>
#include "../include/kalman.h"


int main(){
	clear_file(LOG_ERROR_FILE);
	clear_file(LOG_FILE);
	clear_file(WARNING_FILE);
	try{

		// ========================== FREE FALLING BALL EXAMPLE =============================
		/*| State vector:                                                                   |
		 *| x[k] = [pos[k] vel[k]]                                                          |
		 *|                                                                                 |
		 *|          |1   dt|          |dt^2/2|                                             |
		 *| x[k+1] = |0    1|*x[k] +   |  dt  |*(-g) + w[k] = F*x[k] + B*u[k] + w[k]        |
		 *|                                                                                 |
		 *| It is worth-mentioning that in this model there is no noise since the input is  |
		 *| the gravitational force. However, we could introduce some noise which is related|
		 *| to unmodeled dynamics (e.g.: friction). For dt = 1[s] (f = 1[Hz]) we get:       |
		 *|                                                                                 |
		 *|          |1    1|          |0.5|                                                |
		 *| x[k+1] = |0    1|*x[k] +   | 1 |*(-g) + w[k] = F*x[k] + B*u[k] + w[k]           |
		 *|                                                                                 |
		 *| Assume that we measure directly the height of the ball.                         |
		 *| Our measurements come with some noise.                                          |
		 *|                                                                                 |
		 *| z[k] = [1 0]*x[k]+v[k]                                                          |
		 *===================================================================================
		 */

		// Declare system's matrices
		algebra::mat F, B, Q, H, R;

		// I will use the parameters from the example given at page 15 of
		// http://biorobotics.ri.cmu.edu/papers/sbp_papers/integrated3/kleeman_kalman_basics.pdf

		float freq = 1; // [Hz]
		float dt = 1.0/freq; //sampling period

		// Gravitational force
		float g = 1;

		// State transition matrix
		F = algebra::eye(2);
		F(0,1) = dt;

		// Standard deviation of process noise
		float q_std = 0;

		// Process noise variance
		algebra::mat q(1,1); q(0,0) = q_std*q_std;

		// Control input matrix
		B.set_size(2,1); B(0,0) = dt*dt/2.0; B(1,0) = dt;

		// Auxiliary matrix to shape proper noise variance matrix.
		algebra::mat G(2,1); G(0,0) = dt*dt/2.0; G(1,0) = dt;

		// Process noise variance matrix
		Q = G*q*transpose(G);

		// Observation matrix
		H = "[1 0]";

		// Observation noise variance
		R.set_size(1,1); R(0,0) = 1;

		// True initial conditions
		float pos_0 = 100;
		float vel_0 = 0;
		algebra::vec x0(2); x0(0) = pos_0; x0(1) = vel_0;

		// Create an LTI system
		algebra::lti_system sys;
		sys.set_system(F, B, Q, H, R, dt);

		// Set size of position and velocity vectors
		size_t N = 6;

		// Declare vectors for true position and velocity
		algebra::vec pos_true(N), vel_true(N);

		// Declare position-measurement vector
		algebra::vec pos_meas(N);

		// Set true initial conditions
		pos_true(0) = pos_0;
		vel_true(0) = vel_0;

		// Set input to be u = -g
		algebra::vec u(1); u(0) = -g;

		for(size_t i = 1; i < N; i++){
			sys.run_model(x0, u);
			pos_true(i) = (sys.get_output())(0);
			vel_true(i) = (sys.get_state())(1);
		}

		// Normally velocity is not observable in this system
		// so, we shouldn't be able to access it. But, for the
		// sake of this example I assume I can get it.

		// Position measurements taken from page 24 of
		// http://biorobotics.ri.cmu.edu/papers/sbp_papers/integrated3/kleeman_kalman_basics.pdf
		pos_meas = "[0 100 97.9 94.4 92.7 87.3]";

		// There is no measurement at time 0
		pos_meas(0) = NaN(double);

		// Introduce Kalman filter object
		algebra::kalman kalman;
		algebra::vec pos_hat(N), vel_hat(N), est_err_pos(N), est_err_vel(N);

		// Initial guess to feed the Kalman filter
		algebra::vec x_hat0; x_hat0 = "[95 1]";

		// Uncertainty of the initial guess
		algebra::mat P0 = algebra::eye(2);
		P0(0,0) = 10; P0(1,1) = 1;

		kalman.set_initial_conditions(x_hat0, P0);

		pos_hat(0) = x_hat0(0);
		vel_hat(0) = x_hat0(1);
		est_err_pos(0) = P0(0,0);
		est_err_vel(0) = P0(1,1);

		// Introduce measurement vector to feed the Kalman filter
		algebra::vec y(1);

		for(size_t i = 1; i < N; i++){
			y(0) = pos_meas(i);
			kalman.update(sys, u, y);
			pos_hat(i) = kalman.get_estimate()(0);
			vel_hat(i) = kalman.get_estimate()(1);
			est_err_pos(i) = kalman.get_cov_error()[0];
			est_err_vel(i) = kalman.get_cov_error()[1];
		}

		printf("pos_true = "); print(pos_true);
		printf("pos_meas = "); print(pos_meas);
		printf("pos_hat  = "); print(pos_hat);
		printf("\n");
		printf("vel_true = "); print(vel_true);
		printf("vel_hat  = "); print(vel_hat);
		printf("\n");
		printf("est_err_pos = "); print(est_err_pos);
		printf("est_err_vel = "); print(est_err_vel);

	}catch(const std::exception& e){
		std::cerr << "EXCEPTION CAUGHT: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}
