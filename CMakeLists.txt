cmake_minimum_required(VERSION 2.0)

project(Netflix)

FIND_PACKAGE( OpenMP REQUIRED ) 

if(OPENMP_FOUND)
 set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
 set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_executable(Netflix Constant.cpp DotProductGaussianAddNG.cpp GammaNode.cpp GaussianNode.cpp InferenceEngineNG.cpp MassGaussianNode.cpp MD5.cpp NetflixDataNG.cpp PreProc.cpp)
