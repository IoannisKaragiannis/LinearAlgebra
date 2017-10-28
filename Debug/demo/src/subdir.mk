################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../demo/src/KalmanFilter.cpp \
../demo/src/kalman.cpp \
../demo/src/lti_system.cpp 

OBJS += \
./demo/src/KalmanFilter.o \
./demo/src/kalman.o \
./demo/src/lti_system.o 

CPP_DEPS += \
./demo/src/KalmanFilter.d \
./demo/src/kalman.d \
./demo/src/lti_system.d 


# Each subdirectory must supply rules for building sources it contributes
demo/src/%.o: ../demo/src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++-4.9 -std=c++14 -D__cplusplus=201103L -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


