cmake_minimum_required(VERSION 3.13)
project (rANS LANGUAGES CXX DESCRIPTION "Streaming rANS")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if (APPLE)
	if (CMAKE_HOST_SYSTEM_PROCESSOR MATCHES "^arm64")
		if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
			message(STATUS "add -mcpu=apple-m1")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mcpu=apple-m1 -mtune=native")
		endif()
	endif()
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif()

if (CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS} -O0 -g")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS} -O3")
	set(CMAKE_CXX_FLAGS_RelWithDebInfo "${CMAKE_CXX_FLAGS} -O2 -g")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(st-rans Streaming_rANS.cpp BitStream.cpp)
add_executable(bs-test bs_test.cpp BitStream.cpp)
#add_executable(rans rANS.cpp)
