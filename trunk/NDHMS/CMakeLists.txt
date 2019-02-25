cmake_minimum_required (VERSION 2.8)
project (National_Water_Model)

#set versions number
set (National_Water_Model_VERSION_MAJOR 2)
set (National_Water_model_VERSION_MINOR 2)

# set cmake to work with fortran
enable_language (Fortran)

#try to find the installed MPI library
find_package(MPI REQUIRED)

#set compile flags based on compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    message ( "Using gfortran")
    set (CMAKE_Fortran_FLAGS "-O2 -w -c -ffree-form -ffree-line-length-none -fconvert=big-endian -frecord-marker=4" -I ${MPI_INCLUDE_PATH})
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
    message ( "Using ifort")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif ( Fortran_COMPILER_NAME MATCHES "gfortran.*" )

set(CMAKE_BINARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mods)

include_directories(AFTER ${PROJECT_BINARY_DIR}/mods)

add_subdirectory("Routing/Overland")
add_subdirectory("Routing/Subsurface")
add_subdirectory("Routing/Groundwater")
add_subdirectory("Data_Rec")
add_subdirectory("MPP")
add_subdirectory("Routing")

