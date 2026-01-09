cmake_minimum_required( VERSION 3.6.2)

#Project name:
set( PROJECT_NAME matrixSD )

project( ${PROJECT_NAME})

#Type of build - 2 options to choose: Debug is easier for writting and debuging
#Release better for optimalization, speed, etc.
set( CMAKE_BUILD_TYPE Debug )
#set( CMAKE_BUILD_TYPE Release )

#Complier (ver. c++20 should be stable)
set ( CMAKE_CXX_FLAGS "-std=c++20 -Wall")
set ( CMAKE_CXX_FLAGS_DEBUG "-std=c++20 -Wall")

#Include files:
include_directories( include )

#Source files:
file( GLOB SOURCES "./src/*.cpp" "./include/*.hpp" )
add_executable( ${PROJECT_NAME} ${SOURCES} )

message( "CMAKE_BUIL_TYPE is ${CMAKE_BUILD_TYPE}")