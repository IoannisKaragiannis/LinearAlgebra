################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/KalmanFilter.cpp \
../src/kalman.cpp \
../src/lti_system.cpp 

OBJS += \
./src/KalmanFilter.o \
./src/kalman.o \
./src/lti_system.o 

CPP_DEPS += \
./src/KalmanFilter.d \
./src/kalman.d \
./src/lti_system.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++-4.8 -std=c++0x -D__cplusplus=201103L -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


