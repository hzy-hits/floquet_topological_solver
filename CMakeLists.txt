cmake_minimum_required(VERSION 3.10.0)
project(Mat_Solver VERSION 0.1.0 LANGUAGES CXX)

include(CTest)
enable_testing()

# Find OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
endif()

# Find Eigen
find_package(Eigen3 3.3 REQUIRED)

# Include directories
include_directories(include/cpp)

# Gather source files
file(GLOB SRC_FILES "src/*.cpp")

# Add executable
add_executable(Mat_Solver ${SRC_FILES})

# Link OpenMP
if(OpenMP_CXX_FOUND)
    target_link_libraries(Mat_Solver PUBLIC OpenMP::OpenMP_CXX)
endif()

# Link Eigen
target_link_libraries(Mat_Solver PUBLIC Eigen3::Eigen)

# Set default build type to Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Set Debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-g")

# Set Release flags with optimizations
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
