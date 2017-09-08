cmake_minimum_required(VERSION 3.6)
project(kc_kmeans)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp kc_kmeans.h kc_kmeans.cpp)
add_executable(kc_kmeans ${SOURCE_FILES})