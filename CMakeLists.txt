# Definition of the minimum required cmake Version
cmake_minimum_required(VERSION 3.0)

# Definition of the Project
project(FLUIDCHEN_EM VERSION 1.0)

# Define all configuration options
option(gpp9 "compile with gpp9 filesystem" ON)

# Definition of the C++ Standard 
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(gpp9 True)

# flags for the compiler in debug mode.
set( CMAKE_CXX_FLAGS_DEBUG " -pg -g" )

# flags for the compiler in release mode
set( CMAKE_CXX_FLAGS_RELEASE " -O3 -s -DNDEBUG" )

# flags for the compiler in precice mode
set( CMAKE_CXX_FLAGS_PRECICE " -O3 -s -DPRECICE" )

if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE RELEASE )
endif()

# Creating the executable of our project and the required dependencies
# the executable is called fluidchen_em
file(GLOB files src/*.cpp)
add_executable(fluidchen_em ${files})

target_compile_options(fluidchen_em PUBLIC "-Wno-trigraphs")
target_compile_definitions(fluidchen_em PUBLIC -Dsolution_liddriven)
target_compile_definitions(fluidchen_em PUBLIC -Dsolution_energy)
target_compile_definitions(fluidchen_em PUBLIC -Dsolution_parallelization)


# find_package(MPI)
find_package(MPI REQUIRED)
# if( CMAKE_BUILD_TYPE EQUAL PRECICE)
# find_package(precice REQUIRED CONFIG)
# message (STATUS "FOUND PREICE")
# endif()
find_package(precice REQUIRED CONFIG)
# VTK Library
find_package(VTK REQUIRED)
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
include(${VTK_USE_FILE})

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "9")
    message("g++ Version is lower than Version 9")
    if (NOT APPLE)  
      target_link_libraries(fluidchen_em PUBLIC stdc++fs)
    endif()
    else()
    message("g++ Version is 9 or higher")
    target_compile_definitions(fluidchen_em PUBLIC gpp9)
    target_compile_definitions(fluidchen_em PUBLIC -DGCC_VERSION_9_OR_HIGHER)
  endif()
endif()

# Add include directory
target_include_directories(fluidchen_em PUBLIC include)

if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX /usr/local)
endif()

target_link_libraries(fluidchen_em PRIVATE MPI::MPI_CXX)
target_link_libraries(fluidchen_em PRIVATE ${VTK_LIBRARIES})
target_link_libraries(fluidchen_em PRIVATE precice::precice)

install(TARGETS fluidchen_em DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

