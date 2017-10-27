# ------------------------------------------------
# Generic Makefile
#
# Author: Ioannis Karagiannis
# email : ioanniskaragiannis1987@gmail.com
# Date  : 2017-10-25
# ------------------------------------------------

# Choose only between the two platforms: [LINUX, ARM]
# Change the PLATFORM variable to ARM for cross-compilation
PLATFORM = LINUX

# In case you need to use different cross-compiler for your
# platform, then extend the code below with an extra else if()
ifeq ($(PLATFORM), LINUX)
CXX = g++-4.9
CXX_LINKER = g++-4.9
# Define the flags for your compiler
CXXFLAGS = -std=c++14 -O3 -Wall -fmessage-length=0
else
ifeq ($(PLATFORM), ARM)
CXX = arm-linux-gnueabihf-g++-4.8
CXX_LINKER = arm-linux-gnueabihf-g++-4.9
# Define the flags for your compiler
CXXFLAGS = -std=c++14 -O3 -Wall -fmessage-length=0
endif
endif

# Define the build folder where the objects and 
# the binary file will be stored.
BUILD_FOLDER = ./build

TARGET = LinearAlgebraUnitTest
TARGET_DIR = $(BUILD_FOLDER)/bin
TARGET_FULL_PATH = $(TARGET_DIR)/$(TARGET)
OBJ_DIR = $(BUILD_FOLDER)/objects

# include the source folder of your project (with all the .cpp files)
SRC_DIR = src

# Command to remove files after calling "make clean"
RM = sudo rm -rf

# Create directory. -p checkes whether the directory already exists.
# -m 777 assigns the folder root permissions
MKDIR_P = sudo mkdir -p -m 777


# This line collects all the source files (.cpp) in the current directory
# and saves them in the CPP_SRCS variable.
CPP_SRCS = $(wildcard $(SRC_DIR)/*.cpp)

# This line transforms the content of the CPP_SRCS variable, changing all
# file suffixes from .cpp to .o, thus constructing the object file list we need.
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(CPP_SRCS))

# This line includes one dependency file (.d) for each source file (.cpp).
# Remember, each .o file is both a target and a dependency.
CPP_DEPS = $(patsubst $(OBJ_DIR)/%.o,$(OBJ_DIR)/%.d,$(OBJS))

-include $(CPP_DEPS)   # include all dep files in the makefile

# Each subdirectory must supply rules for building sources it contributes
./$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	@$(MKDIR_P) $(OBJ_DIR)/
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -c -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '
	
LIBS =

# All Target
all: $(TARGET)

# Tool invocations
$(TARGET): $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	@$(MKDIR_P) $(TARGET_DIR)
	$(CXX_LINKER) -o $(TARGET_FULL_PATH) $(OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '


# Other Targets
clean:
	-$(RM) $(OBJS) $(CPP_DEPS) $(TARGET_FULL_PATH) 
	-$(RM) $(BUILD_FOLDER)
	-@echo ' '
	
	
# Run executable	
run:
	@echo '********* RUN UNIT TEST *************'
	@echo '*************************************'
	@echo ' '
	@build/bin/./LinearAlgebraUnitTest
	@echo ' '

.PHONY: all clean dependents
.SECONDARY:
