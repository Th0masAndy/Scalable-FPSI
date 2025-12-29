cmake_minimum_required(VERSION 3.10)
project(FPSI)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_TYPE Release)

set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/thirdparty/out/install")
include_directories(include)

find_package(volePSI REQUIRED)

file(GLOB SOURCES src/*.cpp)

add_executable(fpsi ${SOURCES})

target_link_libraries(fpsi PRIVATE oc::cryptoTools oc::libOTe visa::volePSI visa::volePSI)
target_link_libraries(fpsi PRIVATE ${CMAKE_SOURCE_DIR}/thirdparty/out/install/lib/libblake3.a)

if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|i.86")
    target_compile_options(fpsi PRIVATE -mpclmul -maes -msse2 -msse4.1)
endif()

add_definitions(-DNDEBUG)