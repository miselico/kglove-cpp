################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/graph/LabeledGraph.cpp 

OBJS += \
./src/graph/LabeledGraph.o 

CPP_DEPS += \
./src/graph/LabeledGraph.d 


# Each subdirectory must supply rules for building sources it contributes
src/graph/%.o: ../src/graph/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++1y -I/home/hongyu/Data/libs/boost_1_74_0/ -O0 -g3 -pedantic -Wall -Wextra -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


