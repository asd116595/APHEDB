cmake_minimum_required(VERSION 3.16.3)

project(fhe_test CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Debug)
option(BUILD_STATIC "Set to ON to include static versions of the library" OFF)

find_package(OpenFHE CONFIG REQUIRED)
if(OpenFHE_FOUND)
    message(STATUS "FOUND PACKAGE OpenFHE")
    message(STATUS "OpenFHE Version: ${BASE_OPENFHE_VERSION}")
    message(STATUS "OpenFHE installed as shared libraries: ${OpenFHE_SHARED}")
    message(STATUS "OpenFHE include files location: ${OpenFHE_INCLUDE}")
    message(STATUS "OpenFHE lib files location: ${OpenFHE_LIBDIR}")
    message(STATUS "OpenFHE Native Backend size: ${OpenFHE_NATIVE_SIZE}")
else()
    message(FATAL_ERROR "PACKAGE OpenFHE NOT FOUND")
endif()

set(OPENFHE_SOURCE_DIR ${CMAKE_SOURCE_DIR}/thirdparty/openfhe-development)

set(CMAKE_CXX_FLAGS ${OpenFHE_CXX_FLAGS})

include_directories(${OPENMP_INCLUDES})
include_directories(${OpenFHE_INCLUDE})
include_directories(${OpenFHE_INCLUDE}/third-party/include)
include_directories(${OpenFHE_INCLUDE}/core)
include_directories(${OpenFHE_INCLUDE}/pke)
include_directories(${OpenFHE_INCLUDE}/binfhe)

# Include project headers
include_directories(${CMAKE_SOURCE_DIR}/include)

if(UNIX AND NOT APPLE)
    add_link_options(-Wl,--no-as-needed)
endif()

link_directories(${OpenFHE_LIBDIR})
link_directories(${OPENMP_LIBRARIES})

if(BUILD_STATIC)
    set(CMAKE_EXE_LINKER_FLAGS "${OpenFHE_EXE_LINKER_FLAGS} -static")
    link_libraries(${OpenFHE_STATIC_LIBRARIES})
else()
    set(CMAKE_EXE_LINKER_FLAGS ${OpenFHE_EXE_LINKER_FLAGS})
    link_libraries(${OpenFHE_SHARED_LIBRARIES})
endif()

# Source files
set(SOURCES
    src/equal_bootstrapping.cpp
    src/search.cpp
    src/comparsion.cpp
    src/group_equal.cpp
    src/aggregation.cpp
    src/range.cpp
    src/homlike.cpp
    src/multi_precision.cpp
)

# Create a library from source files
add_library(fhe_lib ${SOURCES})

# Test executables
add_executable(test_equal_bootstrapping test/test_equal_bootstrapping.cpp)
add_executable(test_search test/test_search.cpp)
add_executable(test_comparsion test/test_comparsion.cpp)
add_executable(test_group_equal test/test_group_equal.cpp)
add_executable(test_aggregation test/test_aggregation.cpp)
add_executable(test_range test/test_range.cpp)
add_executable(test_homlike test/test_homlike.cpp)
add_executable(test_multi_precision test/test_multi_precision.cpp)


# Link test executables with the library
target_link_libraries(test_equal_bootstrapping fhe_lib)
target_link_libraries(test_search fhe_lib)
target_link_libraries(test_comparsion fhe_lib)
target_link_libraries(test_group_equal fhe_lib)
target_link_libraries(test_aggregation fhe_lib)
target_link_libraries(test_range fhe_lib)
target_link_libraries(test_homlike fhe_lib)
target_link_libraries(test_multi_precision fhe_lib)
