cmake_minimum_required(VERSION 3.16)
project(DP_Sliding_window)

set(CMAKE_CXX_STANDARD 14)

add_executable(DP_Sliding_window main.cpp Mechanism/new_alg.cpp Mechanism/new_alg.h Mechanism/md5.h Mechanism/md5.cpp Mechanism/PCC.cpp Mechanism/PCC.h Mechanism/DPSW.cpp Mechanism/DPSW.h Mechanism/Upadhyay.cpp Mechanism/Upadhyay.h Mechanism/PCC_counter.h Mechanism/PCC_counter.cpp)
