cmake_minimum_required (VERSION 2.8)

#set(CMAKE_C_COMPILER "mpicc")
#set(CMAKE_CXX_COMPILER "mpic++")
#set(CMAKE_Fortran_COMPILER "mpif90")

project (National_Water_Model)

#set versions number
set (National_Water_Model_VERSION_MAJOR 2)
set (National_Water_model_VERSION_MINOR 2)

# set cmake to work with fortran
enable_language (Fortran)

#try to find the installed MPI library
find_package(MPI REQUIRED)
#message ("-- MPI Include directories: " ${MPI_INCLUDE_PATH} )
#message ("-- MPI C COMPILER : " ${MPI_C_COMPILER} )
#message ("-- MPI CXX COMPILER : " ${MPI_CXX_COMPILER} )
#message ("-- MPI Fortran COMPILER : " ${MPI_Fortran_COMPILER} )
#message ("-- MPI COMPILE FLAGS : " ${MPI_COMPILE_FLAGS} )
#message ("-- MPI LINK FLAGS : " ${MPI_LINK_FLAGS} )
#message ("-- MPI Fortran LINK FLAGS : " ${MPI_Fortran_LINK_FLAGS} )
#message ("-- MPI LIBRARY : " ${MPI_LIBRARY} )
#message ("-- MPI EXTRA LIBRARY : " ${MPI_EXTRA_LIBRARY} )
#message ("-- MPI LIBRARIES : " ${MPI_LIBRARIES} )

# netcdf does not have a built in package locator so add custom module directory
# that contains FindNetCDF.cmake
set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake-modules )

#try to find the installed NETCDF library
set (NETCDF_F90 "YES")
set (NETCDF_F77 "YES")
find_package(NetCDF REQUIRED)
message("-- NetCDF Include Dir: ${NETCDF_INCLUDES}")

#set user controled enviorment variables
set (NWM_LSM $ENV{NWM_LSM} CACHE STRING "Name of the Land Surface Model to Use")

# set enviorment variables if they have not been set
if ("${NWM_LSM}" STREQUAL "")
	set (NWM_LSM "NoahMP")
	message ("-- Setting LSM to: NoahMP")
endif("${NWM_LSM}" STREQUAL "")

#set compile flags based on compiler id

if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU.*")
    # set compile flags for gfortran
    message ( "-- Using gfortran")
    set (FPP_FLAGS "-DMPP_LAND -DHYDRO_D -DWRF_HYDRO -DSPATIAL_SOIL -DWRFIO_NCD_LARGE_FILE_SUPPORT -DWRF_HYDRO_NUDGING")
    set (CMAKE_Fortran_FLAGS " ${FPP_FLAGS} -O2 -w -c -ffree-form -ffree-line-length-none -fconvert=big-endian -frecord-marker=4")
elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel.*")
    # set compile flags for ifort
    message ( "Using ifort")
else (CMAKE_Fortran_COMPILER_ID MATCHES "GNU.*")
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif ( CMAKE_Fortran_COMPILER_ID MATCHES "GNU.*" )

#set output directories for libraries binaries and fortran .mod files
set(CMAKE_BINARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/mods)

#add common include directores need in the build
include_directories(AFTER ${PROJECT_BINARY_DIR}/mods)
include_directories(AFTER ${MPI_INCLUDE_PATH})
include_directories(AFTER ${NETCDF_INCLUDES})
include_directories(AFTER ${PROJECT_SOURCE_DIR}/Data_Rec)

# build the various sup-projects
add_subdirectory("utils")
add_subdirectory("MPP")
add_subdirectory("Debug_Utilities")
add_subdirectory("Routing/Overland")
add_subdirectory("Routing/Subsurface")
add_subdirectory("Routing/Groundwater")
add_subdirectory("Data_Rec")
add_subdirectory("nudging")
add_subdirectory("Routing")
add_subdirectory("HYDRO_drv")

if (NWM_LSM MATCHES "NoahMP")
	message ("-- Building NoahMP LSM")
	add_subdirectory("Land_models/NoahMP")	
elseif(NWM_LSM MATCHES "Noah")
	message ("-- Building Noah LSM")
else()
	message ("Unknown Land Surface Model:" ${NWM_LSM} )
endif (NWM_LSM MATCHES "NoahMP")

