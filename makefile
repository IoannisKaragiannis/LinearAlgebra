
#configuration settings


# If you want to cross-compile your library for ARM processor
# uncomment/comment accordingly the following 4 lines.
CXX = g++-4.8
AR = ar
#CXX = arm-linux-gnueabihf-g++-4.8
#AR = arm-linux-gnueabihf-ar


CXXFLAGS = -std=c++0x -D__cplusplus=201103L -O3 -Wall -fmessage-length=0

# Command to remove files after calling "make clean"
RM = sudo rm -rf

# Create directory. -p checkes whether the directory already exists.
# -m 777 assigns the folder root permissions
MKDIR_P = sudo mkdir -p -m 777

LIB_NAME = "libLinearAlgebra.a"
LIB_DIR = lib
LIB_FULL_PATH = $(LIB_DIR)/$(LIB_NAME)


# All of the sources participating in the build are defined here

# Every subdirectory with source files must be described here

BIN_DIR = bin

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
ARCHIVES := 
CPP_DEPS := 
CXX_DEPS := 
C_UPPER_DEPS := 


# Add inputs and outputs from these tool invocations to the build variables 

CPP_SRCS = $(wildcard *.cpp)

OBJS = ./$(BIN_DIR)/kalman.o ./$(BIN_DIR)/lti_system.o 

CPP_DEPS = ./$(BIN_DIR)/kalman.d ./$(BIN_DIR)/lti_system.d 


# Each subdirectory must supply rules for building sources it contributes
$(BIN_DIR)/%.o: src/miscellaneous/%.cpp
	@$(MKDIR_P) $(BIN_DIR)/
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -c -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


USER_OBJS =

LIBS =


# Add inputs and outputs from these tool invocations to the build variables 

# All Target
all: libLinearAlgebra.a

# Tool invocations
libLinearAlgebra.a: $(OBJS) $(USER_OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC Archiver'
	@$(MKDIR_P) lib/
	$(AR) -r $(LIB_FULL_PATH) $(OBJS) $(USER_OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

# Other Targets
clean:
	-$(RM) $(OBJS)$(C++_DEPS)$(C_DEPS)$(CC_DEPS)$(ARCHIVES)$(CPP_DEPS)$(CXX_DEPS)$(C_UPPER_DEPS) $(LIB_FULL_PATH) 
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
