cmake_minimum_required(VERSION 3.10.0)
project(Mat_Solver VERSION 0.1.0 LANGUAGES C CXX)

include(CTest)
enable_testing()

# 找到OpenMP
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found")
endif()
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
# 添加include目录
include_directories(include/cpp)

# 将src目录下的所有源文件赋值给变量SRC_FILES
file(GLOB SRC_FILES "src/*.cpp")

# 添加可执行文件，使用变量SRC_FILES
add_executable(Mat_Solver ${SRC_FILES})

# 链接OpenMP
if(OpenMP_CXX_FOUND)
    target_link_libraries(Mat_Solver PUBLIC OpenMP::OpenMP_CXX)
endif()

# 设置构建类型的默认值
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# 设置Debug模式的编译标志
set(CMAKE_CXX_FLAGS_DEBUG "-g")

# 设置Release模式的编译标志（包括O2优化）
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
