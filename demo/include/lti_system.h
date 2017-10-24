/*===============================================================================
 * Name         : lti_system.h implements a discrete linear time-invariant system.
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
=================================================================================*/

/* For the sake of consistency with the most common notation,
 *  I used capital letters for most of the member variables of this class. */

#ifndef LTI_SYSTEM_H_
#define LTI_SYSTEM_H_

#include </usr/local/LinearAlgebra/include/base.h>

namespace algebra {

class lti_system {
public:
	lti_system();
	~lti_system();

	mat get_state_transition_matrix();          // F
	mat get_control_matrix();                   // B
	mat get_process_noise_variance();           // Q

	mat get_observation_matrix();               // H
	mat get_observation_noise_variance();       // R
	double get_sampling_period();               // dt

	void set_system(const mat& m1, const mat& m2, const mat& m3,
			const mat& m4, const mat& m5, double sampling_period);

	void run_model(const vec& x0, const vec& input);

	vec get_state();
	vec get_output();

protected:

	void set_initial_conditions(const vec& x0);
	void check_dimension_mismatch();

private:
	/* ================= GENERIC STATE SPACE MODEL =========================
	 *
	 * ======================================================
	 * ============  SYSTEM DYNAMIC MODEL  ==================
	 * ======================================================
	 * | x[k+1] = f(x[k],u[k]) + w[k], where w~N(0, Q)      |
	 * |                                                    |
	 * | For Linear/Linearized Time Invariant systems we get|
	 * |                                                    |
	 * | x[k+1] = F*x[k] + B*u[k] + w[k], where w~N(0, Q)   |
	 * |                                                    |
	 * | Remember: F = exp(A*dt). This is the discretized   |
	 * |           version of a continuous system           |
	 * ======================================================
	 *   */

	mat x;          // state vector of the system
	mat u;          // input of the system

	mat F;          // State transition matrix.
	mat B;          // Control matrix
	mat Q;          // Process noise variance


	/* ====================================================================
	 * ======================  OBSERVATION MODEL  =========================
	 * ====================================================================
	 * | z[k] = h(x[k]) + v[k] ,  where v~N(0, R)                         |
	 * |                                                                  |
	 * |Observation model can be nonlinear, e.g.: z[k] = cos(x[k]) + v[k] |
	 * |For Linear/Linearized time invariant observation models we get:   |
	 * |                                                                  |
	 * |z[k] = H*x[k] + v[k],  where v~N(0, R)                            |
	 * |                                                                  |
	 * |it is assumed that E[w[i]*v[j]'] = 0 for all i,j =>               |
	 * |=> The two noises are statistically independent.                  |
	 * ====================================================================
	 *  */

	mat z;          // Output of the system

	mat H;          // Observation matrix
	mat R;          // Observation noise variance

	/* ====================================================================
	 * ======================  ESTIMATION OF Q AND R  =====================
	 * ====================================================================
	 * | The covariance matrices Q,R are only going to be used in the     |
	 * | Kalman filter. It's very difficult to get a good estimate of     |
	 * | these matrices. Extensive research has been done in this field   |
	 * | to estimate these covariances from data. One practical approach  |
	 * | to do this is the autocovariance least-squares (ALS) technique   |
	 * | that uses the time-lagged autocovariances of routine operating   |
	 * | data to estimate the covariances. The GNU Octave and Matlab code |
	 * | used to calculate the noise covariance matrices using the ALS    |
	 * | technique is available online under the GNU General Public       |
	 * | License. So, the tuning of these matrices entirely depends on    |
	 * | you.                                                             |
	 * |                                                                  |
	 * | See: http://jbrwww.che.wisc.edu/software/als/                    |
	 * ====================================================================
	 *  */

	// Sampling period:
	double dt;       // [Hz]

	bool initial_conditions;

	std::random_device rd;

};

} /* namespace algebra */

#endif /* LTI_SYSTEM_H_ */
