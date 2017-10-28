################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../tests/unit_test/src/main.cpp \
../tests/unit_test/src/mat_test.cpp \
../tests/unit_test/src/vec_test.cpp 

OBJS += \
./tests/unit_test/src/main.o \
./tests/unit_test/src/mat_test.o \
./tests/unit_test/src/vec_test.o 

CPP_DEPS += \
./tests/unit_test/src/main.d \
./tests/unit_test/src/mat_test.d \
./tests/unit_test/src/vec_test.d 


# Each subdirectory must supply rules for building sources it contributes
tests/unit_test/src/%.o: ../tests/unit_test/src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++-4.9 -std=c++14 -D__cplusplus=201103L -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


