# Compiler and compiler flags
CC := g++
CFLAGS := -std=c++17 -Wall -Wextra

# Directories
SRCDIR := ./src
INCDIR := ./include
BUILDDIR := ./build
BINDIR := ./bin

# Names of the executable files
TARGET := main
TARGET_DEBUG := debug

# Names of the source files
SOURCES := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS := $(SOURCES:$(SRCDIR)/%.cpp=$(BUILDDIR)/%.o)

# Include directories
INCLUDES := -I$(INCDIR)

# Build version
$(BINDIR)/$(TARGET): $(OBJECTS)
	$(CC) $(CFLAGS) -O3 -fopenmp $(INCLUDES) $^ -o $@

# Debug version
$(BINDIR)/$(TARGET_DEBUG): $(OBJECTS)
	$(CC) $(CFLAGS) -O0 $(INCLUDES) $^ -o $@ -g

# Object files
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Clean
clean:
	rm -f $(BUILDDIR)/*.o $(BINDIR)/$(TARGET) $(BINDIR)/$(TARGET_DEBUG)

# Debug
# debug:
# 	$(BINDIR)/$(TARGET_DEBUG): $(OBJECTS)
# 		$(CC) $(CFLAGS) $(INCLUDES) $^ -o $@ -g -O0

.PHONY: clean
