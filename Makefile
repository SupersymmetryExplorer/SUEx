# Compiler and flags
CXX = g++
CXXFLAGS = -O2 -Wall -std=c++17
DEPFLAGS = -MMD -MP

# Executable name
TARGET = susy.exe

# Source and object files
SRCS = susy.c functions.cpp variables.cpp readwrite.cpp electroweak.cpp loop.cpp complex.cpp radiative.cpp initcond.cpp numericx.cpp higgs.cpp
OBJS = $(patsubst %.c,%.o,$(SRCS:.cpp=.o))
DEPS = $(OBJS:.o=.d)

# Quote the generated filenames for native Windows cleanup.
empty :=
space := $(empty) $(empty)
comma := ,
CLEAN_FILES = '$(subst $(space),'$(comma)',$(OBJS) $(DEPS) $(TARGET))'

# Default target
all: $(TARGET)

# Link the final executable
$(TARGET): $(OBJS)
	$(CXX) -o $@ $(OBJS)

# Rebuild objects when compiler flags or build rules change.
$(OBJS): Makefile

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c $< -o $@

%.o: %.c
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -x c++ -c $< -o $@

# Clean up build files
clean:
ifeq ($(OS),Windows_NT)
	powershell -NoProfile -NonInteractive -Command "trap { exit 1 }; Get-Item -LiteralPath $(CLEAN_FILES) -Force -ErrorAction SilentlyContinue | Remove-Item -Force -ErrorAction Stop; exit 0"
else
	rm -f $(OBJS) $(DEPS) $(TARGET)
endif

# Optional: run your program
run: $(TARGET)
	./$(TARGET)

-include $(DEPS)

.PHONY: all clean run
