cmake_minimum_required(VERSION 3.14)
project(lorenz_chaos CXX)
set(CMAKE_CXX_STANDART 17)

if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O2 -march=native")
set(CMAKE_CXX_FLAGS_DEBUG   "-g -O0 -Wall -Wextra")

include_directories(include)

add_executable(lorenz_chaos
	src/main.cpp
	src/lorenz.cpp
	src/rk4.cpp
	src/rk38.cpp
	src/rk23.cpp
)

if(UNIX AND NOT APPLE)
	target_link_libraries(lorenz_chaos std++fc)
endif()