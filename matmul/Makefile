# Build directory
BUILD_DIR := build
# Compiler
CXX := g++

# Compiler flags
CXXFLAGS := -O3 -Wall -Wextra -pedantic -std=c++17 -pthread

# Source files
SRCS := $(wildcard src/*.cxx)
# Header files
HDRS := $(wildcard include/*.hpp)
# Create build directory if it doesn't exist
$(shell mkdir -p $(BUILD_DIR))
# Object files
OBJS := $(SRCS)
# Main object file
MAIN_OBJ := $(BUILD_DIR)/main.o

# Executable
TARGET := $(BUILD_DIR)/MatMulBenchmark

# Default target
all: $(TARGET)

# Linking
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

# Compilation
$(BUILD_DIR)/%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean
clean:
	rm -rf $(BUILD_DIR)