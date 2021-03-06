PROJECT(dp_parking)

cmake_minimum_required(VERSION 2.8.0)
# cmake_minimum_required(VERSION 3.3)

file(GLOB SRCS 
	main.cpp 
	"vendor/lodepng/lodepng.cpp" )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

SET(BASEPATH "${CMAKE_SOURCE_DIR}") 
add_definitions ( -DBASEPATH="${BASEPATH}" )
add_subdirectory(${BASEPATH}/vendor/gflags)

# Use the basepath as a search location for includes
include_directories(${BASEPATH})

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# set(CMAKE_CXX_FLAGS "-march=native -fverbose-asm -ffast-math -fno-math-errno -g -save-temps ${CMAKE_CXX_FLAGS}")
# set(CMAKE_CXX_FLAGS "-march=native -ffast-math -fno-math-errno -g -rdynamic -save-temps ${CMAKE_CXX_FLAGS}")
# set(CMAKE_CXX_FLAGS "-fverbose-asm -ffast-math -fno-math-errno -save-temps ${CMAKE_CXX_FLAGS}")

# set(CMAKE_CXX_FLAGS "-g -save-temps ${CMAKE_CXX_FLAGS}")
# set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-march=native -ffast-math -fno-math-errno ${CMAKE_CXX_FLAGS} -O3")

add_executable(dp_parking ${SRCS})
target_link_libraries(dp_parking gflags pthread)