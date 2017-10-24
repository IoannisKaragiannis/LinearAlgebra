
#configuration settings


# If you want to cross-compile your library for ARM processor
# uncomment/comment accordingly the following 4 lines.
CXX = g++-4.8
AR = ar
#CXX = arm-linux-gnueabihf-g++-4.8
#AR = arm-linux-gnueabihf-ar


CXXFLAGS = -std=c++0x -D__cplusplus=201103L -O3 -Wall -fmessage-length=0

TARGET = KalmanFilter
BUILD_FOLDER = build
TARGET_DEST = $(BUILD_FOLDER)/bin
TARGET_FULL_PATH = $(TARGET_DEST)/$(TARGET)
OBJ_DEST = build/objects

# Command to remove files after calling "make clean"
RM = sudo rm -rf

# Create directory. -p checkes whether the directory already exists.
# -m 777 assigns the folder root permissions
MKDIR_P = sudo mkdir -p -m 777



# Every subdirectory with source files must be described here


O_SRCS := 
CPP_SRCS := 
C_UPPER_SRCS := 
C_SRCS := 
S_UPPER_SRCS := 
OBJ_SRCS := 
ASM_SRCS := 
CXX_SRCS := 
C++_SRCS := 
CC_SRCS := 
OBJS := 
C++_DEPS := 
C_DEPS := 
CC_DEPS := 
CPP_DEPS := 
EXECUTABLES := 
CXX_DEPS := 
C_UPPER_DEPS := 


# Add inputs and outputs from these tool invocations to the build variables 

CPP_SRCS = $(wildcard *.cpp)

OBJS = ./$(OBJ_DEST)/kalman.o ./$(OBJ_DEST)/lti_system.o  ./$(OBJ_DEST)/KalmanFilter.o

CPP_DEPS = ./$(OBJ_DEST)/kalman.d ./$(OBJ_DEST)/lti_system.d ./$(OBJ_DEST)/KalmanFilter.d


# Each subdirectory must supply rules for building sources it contributes
$(OBJ_DEST)/%.o: src/%.cpp
	@$(MKDIR_P) $(OBJ_DEST)/
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -c -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


USER_OBJS =

LIBS =


# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: $(TARGET)

# Tool invocations
$(TARGET): $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	@$(MKDIR_P) $(TARGET_DEST)
	$(CXX) -o $(TARGET_FULL_PATH) $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(EXECUTABLES)$(CPP_DEPS)$(CXX_DEPS)$(C_UPPER_DEPS) $(TARGET_FULL_PATH) 
	-$(RM) $(BUILD_FOLDER)
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
