CXX      := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -O2 -I./include

SRC_DIR  := src
OBJ_DIR  := build
BIN      := colebrook_solver

SRCS := $(SRC_DIR)/main.cpp      \
        $(SRC_DIR)/colebrook.cpp \
        $(SRC_DIR)/methods.cpp   \
        $(SRC_DIR)/utils.cpp

OBJS := $(SRCS:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

.PHONY: all clean run

all: $(OBJ_DIR) $(BIN)

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@
	@echo "  [OK] Собран: $(BIN)"

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@
	@echo "  [CC] $<"

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# printf гарантирует \n на любом sh/bash/dash (в отличие от echo)
run: all
	@printf "100000\n0.0002\n0.05\n" | ./$(BIN)

clean:
	rm -rf $(OBJ_DIR) $(BIN)
	@echo "  [clean] Удалены объектные файлы и бинарник."