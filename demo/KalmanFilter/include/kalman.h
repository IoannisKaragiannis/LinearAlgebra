/*==========================================================================
 * Name         : kalman.h implements an online Kalman filter.
 * Version      : 1.0.0, 15 Oct 2017
 *
 * Copyright (c) 2017 Ioannis Karagiannis
 *
 * All rights reserved
 *
 * This file is part of the LinearAlgebra library.
 *
 * LinearAlgebra is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * any later version.
 *
 * You are free to use this library under the terms of the GNU General
 * Public License, but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LinearAlgebra.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact info: https://www.linkedin.com/in/ioannis-karagiannis-7129394a/
 * 				ioanniskaragiannis1987@gmail.com
=============================================================================*/

#ifndef KALMAN_H_
#define KALMAN_H_

#include "lti_system.h"

namespace algebra {

class kalman {
public:
	kalman();
	~kalman();

	/* =====================================================================
	 * ======================  SET PROPER x[0], P[0]  ======================
	 * =====================================================================
	 * | If your are certain of the initial conditions of x_hat[0], then   |
	 * | set P[0]=0. If not, then use some large values for P[0] so that   |
	 * | the filter will prefer the information from the first measurements|
	 * | over the information already in the model.	In some sense, you     |
	 * | allow the filter to escape from your erroneous guess.             |
	 * =====================================================================
	 */
	void set_initial_conditions(const vec& x0, const mat& P0);

	// Method to update the state estimate, the covariance
	void update(lti_system &sys, const vec& input, const vec& measurement );

	vec get_estimate();

	vec get_cov_error();

protected:
	// Initialization of the Kalman filter parameters
	void initialize_filter(lti_system& sys);

	// Initialize the matrices describing the model
	int state_matrix(double dt);

	// This is only for cases where the time-step is not fixed.
	void update_state_matrix(double dt);

	// Test input and measurement vectors
	void update_input_and_measurement(const vec& input, const vec& measurement);

private:

	// I try to stick with the most common notations from the signal processing field

	mat x_hat;              // State vector
	mat u;                  // Input vector

	mat F;                  // State transition matrix.
	mat B;                  // Control matrix
	mat Q;                  // Process noise variance
	mat P;                  // Error covariance matrix

	mat H;                  // Observation matrix
	mat R;                  // Observation noise variance
	mat z;                  // measurement vector

	mat I;                  // Identity matrix
	mat K;                  // kalman gain


	/* ====================================================================
	 * ======================  ESTIMATION OF Q AND R  =====================
	 * ====================================================================
	 * | The tuning of Q and R is supposed to be done when you defined    |
	 * | your system. For more information on how to tune the parameters  |
	 * | look at the 'lti_system.h' file. No self-tuning is included in   |
	 * | this class. Tuning of Q and R is always going to be a trade-off. |
	 * |                                                                  |
	 * | (a)  Q/R << 1 : Slow tracking but filter robust in noise         |
	 * | (b)  Q/R >> 1 : Fast tracking but noise-sensitive filter         |
	 * ====================================================================
	 */

	bool initialize;
	bool initial_conditions;

};


} /* namespace algebra */

#endif /* KALMAN_H_ */
