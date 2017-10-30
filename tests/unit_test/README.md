# UNIT TEST

This is a unit test based upon *catch.hpp* file. For more details see [CATCH](https://github.com/philsquared/Catch).

## PREREQUISITES

Make sure you have enabled C++14 as described in the main README file.
 
## RUN THE UNIT TEST

**_Step 1_:**

Download and store the *LinearAlgebra* library  under any directory you prefer as desribed in the main README file.

**_Step 2_:**

Navigate to the unit test folder and type the following on your terminal

```
cd /your_path/LinearAlgebra/tests/unit_test
make clean
make all
make run
```

You should see the following:

```
********* RUN UNIT TEST *************
*************************************
 
===============================================================================
All tests passed (2975 assertions in 102 test cases)


```
Be patient. The unit test, unlike any of your projects, is testing all the functions defined in the *vec.h* and *mat.h* files. So, it takes some time to create the executable of the unit test. It's roughly 1.1 [MB]; way larger than the respective executable of the demo which was around 136 [KB]. That's becauase the demo code only used a few functions of the *LinearAlgebra* library, not all of them. 

An error log file will be created under the /tmp/LinearAlgebra directory. This file contains all the exceptions thrown 
while running the unit-test. It reveals all the extreme cases I have taken into consideration. If you can think of a
counter-example that would fail the test feel free to contact me.


## CROSS-COMPILE UNIT TEST

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

Now, you can simply edit the makefile and modify line 11 accordingly. In order to be able to compile with arm-linux-gnueabihf-g++-4.9 you need to work on Ubuntu 16.06. If you are working on Ubuntu 14.04, like I do for instance, you will only be able to install arm-linux-gnueabihf-g++-4.8 which will complain if you try for example to exploit specific C++14 features. I prefer compiling my project directly on the machine, but that's just how I do it. Cross compilers might restrict you some times.


## AUTHORS

* **Ioannis Karagiannis** 

## LICENSE

This project is licensed under the GNU General Public License - see the [LICENSE.md](https://github.com/IoannisKaragiannis/LinearAlgebra/blob/master/LICENSE) file for details
