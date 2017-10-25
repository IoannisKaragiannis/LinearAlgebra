# UNIT TEST

This is a unit test based upon *catch.hpp* file. For more details see [CATCH](https://github.com/philsquared/Catch).

## PREREQUISITES

Make sure you have enabled C++11 as described in the main README file.
 
## RUN THE UNIT TEST

**_Step 1_:**

Download and store the *LinearAlgebra* library  under any directory you prefer as desribed in the main README file.

**_Step 2_:**

Navigate to the unit test folder and type the following on your terminal

```
cd /your_path/LinearAlgebra/tests/unit_test
make clean
make all
```

You should see the following:

```
********* RUN UNIT TEST *************
*************************************
 
===============================================================================
All tests passed (2479 assertions in 98 test cases)

```
Be patient. The unit test, unlike any of your projects, is testing all the functions defined in the *vec.h* and *mat.h* files. So, it takes some time to create the executable of the unit test. It's roughly 1.1 [MB]; way larger than the respective executable of the demo which was around 136 [KB]. That's becauase the demo code only used a few functions of the *LinearAlgebra* library. Not all of them. 


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

Now, you can simply edit the makefile and modify line 11 accordingly.


## AUTHORS

* **Ioannis Karagiannis** 

## LICENSE

This project is licensed under the GNU General Public License - see the [LICENSE.md](https://github.com/IoannisKaragiannis/LinearAlgebra/blob/master/LICENSE) file for details
