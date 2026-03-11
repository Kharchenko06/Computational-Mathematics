CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Iinclude
TARGET = lorenz_chaos
SRCS = src/main.cpp src/lorenz.cpp \
       src/rk4.cpp src/rk38.cpp src/rk23.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ -lstdc++fs

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

run: all
	./$(TARGET)

plot: run
	python3 python/plot_all.py

clean:
	rm -f $(OBJS) $(TARGET)
	-rm -rf data/ 2>/dev/null || true
	mkdir -p data 
	touch data/.gitkeep

.PHONY: all run plot clean