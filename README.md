# LinearAlgebra

This is an easy-to-install linear algebra library. It consists of two major classes, _**vec**_ standing for vector
and _**mat**_ standing for matrix together with a set of accompanying functions and operations.
This library might be useful for projects within the field of signal/image processing,
machine learning, and in general any kind of academic/scientific project.

This is a header-only library which means you don't have to worry about different platforms where the library might be used.
This header-only library simplifies the build process. You don't need to build the library, and you don't need to specify the compiled library during the link step of the build. If you do have a compiled library, you will probably want to build multiple versions of it: One compiled with debugging enabled, another with optimization enabled, and possibly yet another stripped of symbols. And maybe even more for a multi-platform system.

## GETTING STARTED

I have only been working on Eclipse running on Ubuntu, thus I will only provide you with instructions
assuming your OS and IDE is Linux and Eclipse respectively. You should find it out yourself how to adopt
these instructions for different OS or IDE.

### PREREQUISITES

At least g++-4.8 and C++11 enabled.
 

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
$ sudo unzip LinearAlgebra-master.zip 
$ cp -r LinearAlgebra-master/. /usr/local/LinearAlgebra/
```
Clean your *Downloads* folder

```
$ rm -r LinearAlgebra-master
$ rm -r LinearAlgebra-master.zip
```

Now the *LinearAlgebra*  containing all the necessary header files should lay under the */usr/local/LinearAlgebra* directory. Could it be more simple? You can always store the LinearAlgebra anywhere you prefer. Even in your project's directory. It's totally up to you.

### LINK YOUR PROJECT WITH *LinearAlgebra* LIBRARY.


Remember that this is a header-only library based on templates. Therefore, there is no need to compile anything at this step. Go to your main project and include this line 

```c++
#include </usr/local/LinearAlgebra/include/base.h>
```
That was it. Enjoy!

P.S.: If you stored the *LinearAlgebra* in a directory other than */usr/local/* you should make the approriate change like this:

```c++
#include </your_path/LinearAlgebra/include/base.h>
```

## TUTORIAL EXAMPLE

```c++
#include <iostream>
#include </usr/local/LinearAlgebra/include/base.h>

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
#include </usr/local/LinearAlgebra/include/base.h>

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

# DEMO

Look at the demo folder.


So, now you must have familiriazed yourself with the LinearAlgebra library. Start your exciting project using this library and don't forget to say a good word for me.

P.S.: There are several linear algebra libraries out there way more complete and probably better performing than mine (e.g.: itpp, armadillo, eigen, etc.). But this is just a project I wanted to work on and share it with you. Cheers.

## AUTHORS

* **Ioannis Karagiannis** 

## LICENSE

This project is licensed under the GNU General Public License - see the [LICENSE.md](https://github.com/IoannisKaragiannis/LinearAlgebra/blob/master/LICENSE) file for details
