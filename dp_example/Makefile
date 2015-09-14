# the compiler: gcc for C program, define as g++ for C++
CC=g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
#  -std=gnu++11 adds compiler and library support for the ISO C++ 2011
CFLAGS=-O3

# the build target executable:
TARGET = emsf

all: $(TARGET)

$(TARGET): $(TARGET).cpp $(TARGET).h
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cpp

blackjack: blackjack.cpp emsf.cpp
	$(CC) $(CFLAGS) -o blackjack blackjack.cpp emsf.cpp mdp.cpp pisf.cpp

clean:
	$(RM) $(TARGET) blackjack
