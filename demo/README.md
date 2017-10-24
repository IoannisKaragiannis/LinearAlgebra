# DEMO

This is a simple demo to illustrate how the LinearAlgebra works. The two main features of this demo is the *lti_system* class and the *kalman* class. In this demo we will implement a Kalman filter to track the position of a free falling ball given that we measure its position with a noisy instrument. The example was based upon [Kalman](http://biorobotics.ri.cmu.edu/papers/sbp_papers/integrated3/kleeman_kalman_basics.pdf). All functions of this demo are selfexplanatory, so feel free to read the respective header and source files.

## GETTING STARTED

These were the two steps to make my demo collaborate with the *LinearAlgebra* library.

**_Step 1_:**
```
Enabled C++11 and introduced the g++-4.8 compiler as explained in the LinearAlgebra README.md. 
```
**_Step 2_:**
```
Linked demo with the base.h file.
```

I did that only once in the *lti_system.h* file:

```c++
#include </usr/local/LinearAlgebra/include/base.h>
```
 
## RUN THIS DEMO

In order for this demo to be ran you need to store the *LinearAlgebra* library under the */usr/local/LinearAlgebra* directory as desribed in the main README.md. From your terminal type the following:

```
/usr/local/LinearAlgebra/demo/Debug/./KalmanFilter
```

You should see:
```
pos_true = [ 100.00 99.50 98.00 95.50 92.00 87.50 ] 
pos_meas = [ nan 100.00 97.90 94.40 92.70 87.30 ] 
pos_hat  = [ 95.00 99.62 98.43 95.21 92.35 87.68 ] 

vel_true = [ 0.00 -1.00 -2.00 -3.00 -4.00 -5.00 ] 
vel_hat  = [ 1.00 0.38 -1.16 -2.90 -3.69 -4.84 ] 

est_err_pos = [ 10.00 0.92 0.67 0.66 0.61 0.55 ] 
est_err_vel = [ 1.00 0.92 0.58 0.30 0.15 0.08 ] 
```

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

Cross Compile Settings for your project.

Right click on your project:
Properties -> C/C++ Build -> Settings -> Manage Configurations
Press New.. and name it as you wish. For simplicity I usually call them DebugARM/ReleaseARM.

Now, in the *GCC C++ Compiler* and *GCC C++ Linker* command insert **arm-linux-gnueabihf-g++-4.8** or a newer version if you prefer. Respectively insert **arm-linux-gnueabihf-gcc-4.8** in the *GCC C Compiler*. We don't need to specify the path for them as they are stored under */usr/bin/* during installation. This path is by default included in Eclipse. If that's not your case, make sure you inlude the respective path in front of the compiler's name. Find where your cross-compiler is stored by typing the following on your terminal:

```
$ find /usr/ -name "arm-linux-gnueabihf-g++-4.8"
```

In the *GCC C++ Compiler* add the following path in the includes: **/usr/arm-linux-gnueabihf/include**. In the *GCC C++ Linker* add this in the library search path: **/usr/arm-linux-gnueabihf/lib**. Last but not least in *GCC Assembler* change **as** to **/usr/arm-linux-gnueabihf/bin/as**

YOU ARE READY TO CROSS COMPILE YOUR PROJECT!!!


## AUTHORS

* **Ioannis Karagiannis** 

## LICENSE

This project is licensed under the GNU General Public License - see the [LICENSE.md](https://github.com/IoannisKaragiannis/LinearAlgebra/blob/master/LICENSE) file for details
