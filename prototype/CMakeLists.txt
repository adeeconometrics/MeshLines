cmake_minimum_required(VERSION 3.16.3)
project(MeshLines)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)

# Set output directories
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Set compiler flags for debug and release versions
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Wall -Wextra -pedantic")

# Add source and include directories
file(GLOB_RECURSE SRC_FILES src/*.cpp)

add_executable(${PROJECT_NAME} ${SRC_FILES})
target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Fetch GTest library
find_package(GTest REQUIRED CONFIG)

include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)

if(WIN32)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)
endif()

# Add test directory
enable_testing()

add_executable(mesh_test ${CMAKE_CURRENT_SOURCE_DIR}/test/unittest.cpp)
target_link_libraries(mesh_test PRIVATE GTest::gtest_main)
add_test(mesh_test ${CMAKE_CURRENT_SOURCE_DIR}/test/unittest.cpp)

include(GoogleTest)
gtest_discover_tests(mesh_test)

# Specify the packaging information
set(CPACK_PACKAGE_NAME "my_project")
set(CPACK_PACKAGE_VERSION "1.0.0")
set(CPACK_GENERATOR "ZIP")
include(CPack)
