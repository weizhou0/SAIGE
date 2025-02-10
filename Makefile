CXX = g++
CXXFLAGS = -std=c++11 -O2 -Wall -fopenmp -I/usr/local/include -I/usr/local/opt/libomp/include -I/usr/local/opt/armadillo/include
LDFLAGS = -L/usr/local/opt/libomp/lib -L/usr/local/opt/armadillo/lib -larmadillo -fopenmp

SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin

SRCS = $(wildcard $(SRC_DIR)/**/*.cpp $(SRC_DIR)/*.cpp)
OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))

TARGET = $(BIN_DIR)/saige

all: $(TARGET)

$(TARGET): $(OBJS)
    @mkdir -p $(BIN_DIR)
    $(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
    @mkdir -p $(OBJ_DIR)
    $(CXX) $(CXXFLAGS) -c -o $@ $<

clean:
    rm -rf $(OBJ_DIR) $(BIN_DIR)

.PHONY: all clean
