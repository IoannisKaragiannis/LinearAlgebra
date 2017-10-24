# LinearAlgebraLibrary

This is a linear algebra library. It consists of two major classes, _**vec**_ standing for vector
and _**mat**_ standing for matrix together with a set of accompanying functions and operations.
This library might be useful for projects within the field of signal/image processing,
machine learning, and in general any kind of academic/scientific project.

Being a static library it will increase the overall size of your project's binary, but it means that
you won't need to carry along a copy of the library that is being used. As the code is connected
at compile time there are not any additional run-time loading costs. The code is simply there.
Personally, I prefer static libraries to ensure that the binary does not have many external
dependencies that may be difficult to meet, such as specific versions of the C++ standard library
or specific versions of the Boost C++ library.


## GETTING STARTED

I have only been working on Eclipse running on Ubuntu, thus I will only provide you with instructions
assuming your OS and IDE is Linux and Eclipse respectively. You should find it out yourself how to adopt
these instructions for different OS or IDE.

### PREREQUISITES


 

### ENABLE C++11

**_Step 1_:**
```
On your Eclipse environment, configure syntax parser:

Window -> Preferences -> C/C++ -> Build -> Settings -> Discovery -> CDT GCC Build-in Compiler Settings
In the text box entitled "Command to get compiler specs" append -std=c++11
The following screenshot displays the above mentioned settings
```

![Screenshot](/images/LinearAlgebraLibrary/c++11_in_eclipse.png)

**_Step 2_:**
```
Right click on your project:

Properties -> C/C++ Build -> Settings -> Tool Settings -> GCC C++ Compiler -> Dialect
Put -std=c++11 into text box entitled 'other dialect flags' or select ISO C++11 from the 
language standard drop down as shown below:
```

![Screenshot](/images/LinearAlgebraLibrary/c++11_in_project.png "con comment")

**_Step 3_:**
```

Right click on your project
Properties -> C/C++ General -> Path and symbols (Tab) -> GNU C++ and press add
Remember to add symbol '__cplusplus' with value '201103L' as shown below:

```
![Screenshot](/images/LinearAlgebraLibrary/cplusplus_in_project.png)


### COMPILE AT LEAST WITH g++-4.8

**_Step 1_:**

This library was compiled with a g++-4.8 compiler (which includes C++11), so make sure
you have at least this version on your linux-machine. 
```
$ sudo apt-get install python-software-properties
$ sudo add-apt-repository ppa:ubuntu-toolchain-r/test
$ sudo apt-get update
$ sudo apt-get install gcc-4.8
$ sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 50
```

**_Step 2_:**

```
In your Eclipse project force the compiler and the linker to be g++-4.8 as shown below
```
![Screenshot](/images/LinearAlgebraLibrary/g++-48_compiler.png)
![Screenshot](/images/LinearAlgebraLibrary/g++-4.8_linker.png)


### DOWNLOAD AND INSTALL *LinearAlgebra* LIB ON YOUR MACHINE

**_Step 1_:**
```
Download the current repository from the 'Clone or Download' button as as shown below
```

![Screenshot](/images/LinearAlgebraLibrary/download_repo.png)


**_Step 2_:**

If the unzip command isn't already installed on your system, then run:
```
$ sudo apt-get install unzip
```

Unzip the file. Make sure you are root so you have all transferring/deleting permissions).

```
$ sudo su
$ cd /home/user_name/Downloads
$ sudo unzip LinearAlgebraLibrary-master.zip 
$ mv /home/user_name/Downloads/LinearAlgebraLibrary-master/LinearAlgebra/ /usr/local/lib/
```

Now the *LinearAlgebra*  containing the source code, as well as 
the header files should lay under the */usr/local/lib* directory.

**_Step 3_:**

In order to build the library for either Debug or Release configuration type
the following on your terminal:

```
$ cd /usr/local/lib/LinearAlgebra/Debug
$ make clean
$ make all
```
or
```
$ cd /usr/local/lib/LinearAlgebra/Release
$ make clean
$ make all
```
respectively.

### LINK YOUR PROJECT WITH *LinearAlgebra* LIBRARY.

**_Step 1_:**

```
Right click on your project:
Properties -> C/C++ Build -> Settings -> GCC C++ Linker -> Libraries
and add the proper library and library path as shown below.

Depending on your configuration (Debug/Release), you should link
your project accordingly to the proper library that was built as
described above.

```
![Screenshot](/images/LinearAlgebraLibrary/Linker_libraries.png)

**_Step 2_:**

```
Right click on your project:
Properties -> C/C++ Build -> Settings -> GCC C++ Compiler -> Includes
and add the proper path for the header files of the current library
as shown below:

```
![Screenshot](/images/LinearAlgebraLibrary/compiler_includes.png)



## DEMO

```c++
#include <iostream>
#include <stdio.h>
#include </usr/local/lib/LinearAlgebra/base.h>

int main(){
  clear_file(LOG_ERROR_FILE);
  clear_file(LOG_FILE);
  clear_file(WARNING_FILE);
  try{
	//Vectors' and matricies' declaration:
	algebra::vec v1, v2, v3;
	algebra::mat m1, m2, m3;
	
	//Use a string of values to define a vector:
	v1 = "[0.5 -2 4.7 -0.9 8.6 11.3]";
	v2 = "[-1 4.5 -5.9 44 -2.3 6.7]";
		
	//Add vectors v1 and v2
	v3 = v1 + v2;

	//Calculate the mean value of v1
	float mean_v1 = algebra::mean(v1);

	//Calculate the dot product of v1, v2
	float dot_v1v2 = algebra::dot(v1, v2);
		
	//Calculate the norm of v3
	float norm_v3 = algebra::norm(v3);

	// print out the results
	printf("v1 = "); v1.print();
	printf("v2 = "); v2.print();
	printf("v3 = "); v3.print();
	printf("mean(v1) = %.3f \n", mean_v1);
	printf("dot(v1,v2) = %.3f \n", dot_v1v2);
	printf("norm(v3) = %.3f \n", norm_v3);

	//Use a string to define a matrix
	m1 = "[7 2 1;0 3 -1;-3 4 -2]";

	//Calculate the inverse of matrix m1:
	m2 = inv(m1);

	// Confirm that m1*m2 equals to the identity matrix
	m3 = m1 * m2;

	// print out the results
	printf("\n");
	printf("m1 = \n"); m1.print();
	printf("\n");
	printf("inv(m1) = \n"); m2.print();
	printf("\n");
	printf("m1*inv(m1) = \n"); m3.print();

   }catch(const std::exception& e){
	std::cerr << "EXCEPTION CAUGHT: " << e.what() << std::endl;
	return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}


```
If you run this program and see the output below, you have successfully installed and linked
the *LinearAlgebra* library to your project. Congratulations! Enjoy!

```
v1 = [ 0.500 -2.000 4.700 -0.900 8.600 11.300] 
v2 = [ -1.000 4.500 -5.900 44.000 -2.300 6.700] 
v3 = [ -0.500 2.500 -1.200 43.100 6.300 18.000] 
mean(v1) = 3.700 
dot(v1,v2) = -20.900 
norm(v3) = 47.215 

m1 = 
|  7.000 2.000 1.000 |
|  0.000 3.000 -1.000 |
|  -3.000 4.000 -2.000 |

inv(m1) = 
|  -2.000 8.000 -5.000 |
|  3.000 -11.000 7.000 |
|  9.000 -34.000 21.000 |

m1*inv(m1) = 
|  1.000 0.000 0.000 |
|  0.000 1.000 0.000 |
|  0.000 -0.000 1.000 |

```
If for any reason an exception is thrown you can navigate to the */tmp/LinearAlgebra/* directory
where the *error_log.txt* file will guide you to debug your code. 

## CROSS COMPILE FOR ARM PROCESSORS

Install the ARM cross compiler toolchain on your Linux Ubuntu PC:
```
$ sudo apt-get install libc6-armel-cross libc6-dev-armel-cross
$ sudo apt-get install binutils-arm-linux-gnueabi
$ sudo apt-get install libncurses5-dev
```
If you are using an Arietta, Aria or FOX G20 board:
```
$ sudo apt-get install gcc-arm-linux-gnueabi
$ sudo apt-get install g++-arm-linux-gnueabi
```
If you are using an Acqua or RoadRunner board:
```
$ sudo apt-get install gcc-arm-linux-gnueabihf
$ sudo apt-get install g++-arm-linux-gnueabihf
```

If you have to cross-compile your project using this library, you have to make sure that you link the project to the correct cross-compiled version of the library as well. I have already cross compiled the library for *AM3358 ARM Cortex-A8* processor (BeagleBone). That is why I used gcc-arm-linux-gnueabihf.  It's in the DebugARM or ReleaseARM folder. But if for some reason you need a newer gcc/g++ version or another processor, you have to cross-compile it yourself. 

Cross Compile Settings for your project.

Right click on your project:
Properties -> C/C++ Build -> Settings -> Manage Configurations
Press New.. and name it as you wish. For simplicity I just call them DebugARM/ReleaseARM.

Now, in the *GCC C++ Compiler* and *GCC C++ Linker* command insert **arm-linux-gnueabihf-g++-4.8** or a newer version if you prefer. Respectively insert **arm-linux-gnueabihf-gcc-4.8** in the *GCC C Compiler*. We don't need to specify the path for them as they are stored under */usr/bin/* during installation. This path is by default included in Eclipse. If that's not your case, make sure you inlude the respective path in front of the compiler's name. Find where your cross-compiler is stored by typing the following on your terminal:

```
$ find /usr/ -name "arm-linux-gnueabihf-g++-4.8"
```

In the *GCC C++ Compiler* add the following path in the includes: **/usr/arm-linux-gnueabihf/include**. In the *GCC C++ Linker* add this in the library search path: **/usr/arm-linux-gnueabihf/lib**. Last but not least in *GCC Assembler* change **as** to **/usr/arm-linux-gnueabihf/bin/as**

YOU ARE READY TO CROSS COMPILE YOUR PROJECT!!!



## RUNNING TESTS


### UNIT TESTS

**_Step 1_:**
```
Take the project called UnitTests and load it on Eclipse.
```
**_Step 2_:**
```
Link the project with LinearAlgebra library as described above.
```
**_Step 3_:**
```
Compile and Run the tests. Done!
```
An error log file will be created under the /tmp/LinearAlgebra directory. This file contains all the exceptions thrown 
while running the unit-test. It reveals all the extreme cases I have taken into consideration. If you can think of a
counter-example that would not pass the test feel free to contact me.


### SPEED TESTS

The following script reveals how simple it is to run the speed tests.

```c++
#include <iostream>
#include <stdio.h>
#include </usr/local/lib/LinearAlgebra/base.h>

int main(){
  clear_file(LOG_ERROR_FILE);
  clear_file(LOG_FILE);
  clear_file(WARNING_FILE);
  try{
	test_speed_of_basic_operations(3);
	test_speed_of_basic_operations(4000);

   }catch(const std::exception& e){
	std::cerr << "EXCEPTION CAUGHT: " << e.what() << std::endl;
	return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
```

If you run the above mentioned script you should get something similar to this. It depends on your machine.
MACHINE DETAILS: Intel® Core™ i3 CPU M 330 @ 2.13GHz × 4 (64-bit), 8GB RAM.

```
addition of two squared (3x3) matrices lasted 0.006 [ms]
transposed of a squared (3x3) matrix lasted 0.005 [ms]
normal-multiplication of two squared (3x3) matrices lasted 0.008 [ms]
strassen-multiplication of two squared (3x3) matrices lasted 0.861 [ms]
inversion of a squared (3x3) matrix lasted 0.018 [ms]
addition of two squared (4000x4000) matrices lasted 433.425 [ms]
transposed of a squared (4000x4000) matrix lasted 845.269 [ms]
normal-multiplication of two squared (4000x4000) matrices lasted 1338710.784 [ms]
strassen-multiplication of two squared (4000x4000) matrices lasted 856570.816 [ms]
inversion of a squared (4000x4000) matrix lasted 3010656.512 [ms]
```
In conclusion, one should think twice before multiplying or inverting large matrices both time-wise and memory-wise speaking. Consider the memory these matrices occupy. A 4000x4000 matrix consisting of floats occupies 64MB. For a single multiplication we would need three times this memory. However, one can notice how useful Strassen multiplication becomes for large matrices. In the above example for two 4000x4000 matrices the normal multiplication lasted 22.3 minutes when it only took 14.3 for Strassen.

## LOG FILE

Assuming you are using Linux, an error-log file will be stored under the
/tmp/LinearAlgebra/ directory. This will be helpful to detect bugs in your
code related to *LinearAlgebra* lib.


## EXTENSIVE EXAMPLE USING *lti_system* and *kalman* classes.

The following example was taken from [Kalman example](http://biorobotics.ri.cmu.edu/papers/sbp_papers/integrated3/kleeman_kalman_basics.pdf). You can compare the results of the following code with the table at page 24 of the aforementioned example. In order to run the following code make sure you have linked it with the LinearAlgebra library as described above.


```c++
#include <iostream>
#include <stdio.h>
#include </usr/local/lib/LinearAlgebra/base.h>

int main(){
  clear_file(LOG_ERROR_FILE);
  clear_file(LOG_FILE);
  clear_file(WARNING_FILE);
  try{
      	//test_speed_of_basic_operations(3);

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
	pos_meas(0) = NaN;
	
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

	printf("pos_true = "); pos_true.print(2);
	printf("pos_meas = "); pos_meas.print(2);
	printf("pos_hat  = "); pos_hat.print(2);
	printf("\n");
	printf("vel_true = "); vel_true.print(2);
	printf("vel_hat  = "); vel_hat.print(2);
	printf("\n");
	printf("est_err_pos = "); est_err_pos.print(2);
	printf("est_err_vel = "); est_err_vel.print(2);

   }catch(const std::exception& e){
	std::cerr << "EXCEPTION CAUGHT: " << e.what() << std::endl;
	return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
```
If you run this code it should yield back.
```
pos_true = [ 100.00 99.50 98.00 95.50 92.00 87.50] 
pos_meas = [ nan 100.00 97.90 94.40 92.70 87.30] 
pos_hat  = [ 95.00 99.62 98.43 95.21 92.35 87.68] 

vel_true = [ 0.00 -1.00 -2.00 -3.00 -4.00 -5.00] 
vel_hat  = [ 1.00 0.38 -1.16 -2.90 -3.69 -4.84] 

est_err_pos = [ 10.00 0.92 0.67 0.66 0.61 0.55] 
est_err_vel = [ 1.00 0.92 0.58 0.30 0.15 0.08] 
```

So, now you must have familiriazed yourself with the LinearAlgebra library. Start your exciting project using this library and don't forget to say a good word for me.

P.S.: There are several linear algebra libraries out there way more complete and probably better performing than mine (e.g.: itpp, armadillo, eigen, etc.). But this is just a project I wanted to work on and share it with you. Cheers.

## AUTHORS

* **Ioannis Karagiannis** 

## LICENSE

This project is licensed under the GNU General Public License - see the [LICENSE.md](https://github.com/IoannisKaragiannis/LinearAlgebraLibrary/blob/master/LICENSE) file for details
