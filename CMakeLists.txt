cmake_minimum_required(VERSION 3.6)
project(HubbardModel2D)

SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_CXX_STANDARD 11)
find_package(OpenMP)

set(boost_req chrono timer system random)
find_package (Boost 1.54.0 COMPONENTS ${boost_req} REQUIRED)
message(STATUS "Found Boost includes: ${Boost_INCLUDE_DIRS}" )
message(STATUS "Found Boost libs: ${Boost_LIBRARIES}" )
include_directories(${Boost_INCLUDE_DIRS})

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  message(STATUS "using GNU compiler, need external lapack")
  if($ENV{CRAYPE_VERSION} MATCHES "." )
    message("On a cray systen, will use libsci")
  else()
  find_package(LAPACK REQUIRED)
endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
  message(STATUS "using llvm compiler, need external lapack")
  find_package(LAPACK REQUIRED)
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  message(STATUS "using intel compiler and parallel mkl")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mkl=parallel -std=c++11")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3 -pthread -DEIGEN_NO_DEBUG") 

set(SOURCE_FILES
        Hubbard2D.cpp
        Hubbard2D.hpp
        )
#find_package (Eigen3 3.3 REQUIRED)
if(NOT IETL_DIR)
  set(IETL_DIR ~/software/ietl)
endif()
if(NOT Spectra_DIR)
  set(Spectra_DIR ~/software/spectra/include)
endif()
include_directories(${EIGEN3_INCLUDE_DIR} ${IETL_DIR} ${Spectra_DIR})

add_library(hubbard ${SOURCE_FILES})
add_executable(HubbardModel2D main.cpp )
target_link_libraries(HubbardModel2D hubbard ${Boost_LIBRARIES} ${LAPACK_LIBRARIES})

add_executable(HubbardModel2D_mult main_mult.cpp)
target_link_libraries(HubbardModel2D_mult hubbard ${Boost_LIBRARIES} ${LAPACK_LIBRARIES})
