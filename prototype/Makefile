CC = g++
CFLAGS = -std=c++17 -Wall -Wextra

BIN_DIR = ./bin
SRC_DIR = ./src
INCLUDE_DIR = ./include

SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(BIN_DIR)/%.o,$(SRCS))

TARGET = proto

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) $^ -o $(BIN_DIR)/$@

$(BIN_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(BIN_DIR)/$(TARGET)

.PHONY: clean
