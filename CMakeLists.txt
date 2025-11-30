cmake_minimum_required(VERSION 3.10)
project(LBM_D2Q9)

set(CMAKE_CXX_STANDARD 17)

add_executable(lbm_simulation src/main.cpp src/LBM.cpp)
