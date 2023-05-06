cmake_minimum_required(VERSION 3.16)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -pedantic")

project("MeshLines"
        VERSION 1.0
        DESCRIPTION "Linear Algebra Routines"
        LANGUAGES CXX)

file(GLOB_RECURSE SRC_FILES src/*.cpp)

add_executable(${PROJECT_NAME} ${SRC_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# trigger automate doxygenation
# find_package(Doxygen
#             REQUIRED dot
#             OPTIONAL_COMPONENTS mscgen dia)