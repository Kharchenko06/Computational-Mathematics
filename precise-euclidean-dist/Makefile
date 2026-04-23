#  Makefile
#  High-Precision Euclidean Distance (IEEE 754 Double)
#
#  Команды:
#    make          — собрать все
#    make run      — собрать и запустить все тесты
#    make demo     — запустить CLI с примером (3-4-5 треугольник)
#    make clean    — удалить все артефакты сборки

CXX      := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -Wpedantic -O2 \
             -fno-unsafe-math-optimizations \
             -ffp-contract=on \
             -march=native

INCLUDES := -Iinclude
LIBS     := -lgmp

#	Директории
SRC_DIR   := src
TEST_DIR  := tests
BUILD_DIR := build

#	Исходники библиотеки
LIB_SRCS := $(SRC_DIR)/ulp_utils.cpp       \
             $(SRC_DIR)/multiply.cpp        \
             $(SRC_DIR)/compensated_sum.cpp \
             $(SRC_DIR)/sqrt_refined.cpp    \
             $(SRC_DIR)/euclidean.cpp       \
             $(SRC_DIR)/gmp_reference.cpp

LIB_OBJS := $(patsubst $(SRC_DIR)/%.cpp, $(BUILD_DIR)/%.o, $(LIB_SRCS))

#	Исполняемые файлы
MAIN_BIN  := $(BUILD_DIR)/euclidean_distance

TEST_SRCS := $(TEST_DIR)/test_multiply.cpp        \
             $(TEST_DIR)/test_compensated_sum.cpp  \
             $(TEST_DIR)/test_sqrt_refined.cpp     \
             $(TEST_DIR)/test_permutation.cpp      \
             $(TEST_DIR)/test_euclidean_full.cpp

TEST_BINS := $(patsubst $(TEST_DIR)/%.cpp, $(BUILD_DIR)/%, $(TEST_SRCS))

#	Все цели
ALL_BINS  := $(MAIN_BIN) $(TEST_BINS)

.PHONY: all run demo clean

all: $(BUILD_DIR) $(ALL_BINS)

#	Создать build/ если не существует
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

#	Компиляция объектных файлов библиотеки
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

#	Линковка main
$(MAIN_BIN): main.cpp $(LIB_OBJS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@ $(LIBS)

#	Линковка тестов
$(BUILD_DIR)/test_%: $(TEST_DIR)/test_%.cpp $(LIB_OBJS) | $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $^ -o $@ $(LIBS)

#	Запуск всех тестов
run: all
	@echo ""
	@echo "══════════════════════════════════════════"
	@echo "  test_multiply"
	@echo "══════════════════════════════════════════"
	@$(BUILD_DIR)/test_multiply

	@echo ""
	@echo "══════════════════════════════════════════"
	@echo "  test_compensated_sum"
	@echo "══════════════════════════════════════════"
	@$(BUILD_DIR)/test_compensated_sum

	@echo ""
	@echo "══════════════════════════════════════════"
	@echo "  test_sqrt_refined"
	@echo "══════════════════════════════════════════"
	@$(BUILD_DIR)/test_sqrt_refined

	@echo ""
	@echo "══════════════════════════════════════════"
	@echo "  test_permutation"
	@echo "══════════════════════════════════════════"
	@$(BUILD_DIR)/test_permutation

	@echo ""
	@echo "══════════════════════════════════════════"
	@echo "  test_euclidean_full"
	@echo "══════════════════════════════════════════"
	@$(BUILD_DIR)/test_euclidean_full

#	Демо: CLI с примером
demo: $(MAIN_BIN)
	@echo ""
	@echo ">>> 3-4-5 triangle: P1=(0,0,0) P2=(3,4,0), expected dist=5"
	@$(MAIN_BIN) 0 0 0 3 4 0
	@echo ""
	@echo ">>> Cancellation: P1=(1e15,1e15,1e15) P2=(1e15+1, ...), expected dist=sqrt(3)"
	@$(MAIN_BIN) 1000000000000000 1000000000000000 1000000000000000 \
	             1000000000000001 1000000000000001 1000000000000001

#	Очистка
clean:
	rm -rf $(BUILD_DIR)