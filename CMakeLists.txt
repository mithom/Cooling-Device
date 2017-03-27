cmake_minimum_required(VERSION 3.6)
project(Cooling_Device)

set(CMAKE_CXX_STANDARD 14)

# add libmesh with all correct flags
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lz -Wl,-rpath,/usr/local/lib -L/usr/local/lib -lmesh_opt -I/usr/local/include -pthread -O2 -felide-constructors -funroll-loops -fstrict-aliasing -Wdisabled-optimization -Wno-variadic-macros -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/home/thomas/CoinIpopt/build/lib -g")

set(SOURCE_FILES Main.cpp)

add_executable(Cooling_Device ${SOURCE_FILES})

find_package(BLAS)
find_package(PythonLibs REQUIRED)
find_package(LAPACK)

find_package(Eigen3 REQUIRED)
include_directories(Eigen3_INCLUDE_DIR)
include_directories(~/CoinIpopt/build/include/coin)

target_link_libraries(Cooling_Device blas)
target_link_libraries(Cooling_Device lapack)
target_link_libraries(Cooling_Device lapacke)
target_link_libraries(Cooling_Device ipopt)
target_link_libraries(Cooling_Device m)
target_link_libraries(Cooling_Device dl)

include_directories(${PYTHON_INCLUDE_DIR})
target_link_libraries(Cooling_Device ${PYTHON_LIBRARIES})