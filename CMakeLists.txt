cmake_minimum_required(VERSION 3.12.1)
project(Regenerator)

include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR}/spdlog)
set(CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

file(GLOB_RECURSE REGEN_SOURCES "src/*.cpp" "src/*.h" "spdlog/spdlog.h")
	
add_executable(Regenerator ${REGEN_SOURCES})