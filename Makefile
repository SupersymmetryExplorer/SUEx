# Compiler and flags
CXX = g++
CXXFLAGS = -O2 -Wall -std=c++17

# Executable name
TARGET = susy.exe

# Source and object files
SRCS = susy.c functions.cpp variables.cpp readwrite.cpp electroweak.cpp loop.cpp complex.cpp radiative.cpp initcond.cpp numericx.cpp higgs.cpp
OBJS = $(SRCS:.cpp=.o)

# Default target
all: $(TARGET)

# Link the final executable
$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS)

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up build files
clean:
	rm -f $(OBJS) $(TARGET)

# Optional: run your program
run: $(TARGET)
	./$(TARGET)
