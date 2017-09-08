cmake_minimum_required(VERSION 3.6)
project(kc_kmeans)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp kmeans.h kmeans.cpp)
add_executable(kmeans ${SOURCE_FILES})