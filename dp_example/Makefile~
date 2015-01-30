# the compiler: gcc for C program, define as g++ for C++
CC = g++

# compiler flags:
#  -g    adds debugging information to the executable file
#  -Wall turns on most, but not all, compiler warnings
#  -std=gnu++11 adds compiler and library support for the ISO C++ 2011
CFLAGS  = -g -Wall

# the build target executable:
TARGET = emsf

all: $(TARGET)

$(TARGET): $(TARGET).cc
	$(CC) $(CFLAGS) -o $(TARGET) $(TARGET).cc

clean:
	$(RM) $(TARGET)
