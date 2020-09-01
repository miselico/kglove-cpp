################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/BCA.cpp \
../src/GraphWeigher.cpp \
../src/KGloVe.cpp \
../src/Main.cpp \
../src/nTripleParser.cpp 

OBJS += \
./src/BCA.o \
./src/GraphWeigher.o \
./src/KGloVe.o \
./src/Main.o \
./src/nTripleParser.o 

CPP_DEPS += \
./src/BCA.d \
./src/GraphWeigher.d \
./src/KGloVe.d \
./src/Main.d \
./src/nTripleParser.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++1y -I/home/hongyu/Data/libs/boost_1_74_0/ -O0 -g3 -pedantic -Wall -Wextra -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


