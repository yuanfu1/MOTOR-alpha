######## This file sets up some common environment for developers:
#-------------------------------------------------------------------------------
# Use Fortran compiler:
ENABLE_LANGUAGE( Fortran )
ENABLE_LANGUAGE( C )
ENABLE_LANGUAGE( CXX )

#-------------------------------------------------------------------------------
# Enable tests:
ENABLE_testing()

#-------------------------------------------------------------------------------
# Build type:
SET ( CMAKE_BUILD_TYPE "Debug" )
SET ( CMAKE_C_COMPILER $ENV{CMAKE_C_COMPILER} )
SET (CMAKE_CXX_STANDARD 11)

#-------------------------------------------------------------------------------
# Common library directory where all lib files are stored:
# Use of an environmental variable of "COMMON_LIBS":
SET (LIBRARY_OUTPUT_PATH $ENV{COMMON_LIBS})
SET (EXECUTABLE_OUTPUT_PATH $ENV{BIN_DIR})

IF (IS_ABSOLUTE $ENV{COMMON_LIBS})
  MESSAGE(STATUS "Lib target path is set!")
ELSE ($ENV{COMMON_LIBS} LESS 1)
  MESSAGE(FATAL_ERROR "Env COMMON_LIBS is not set, please set and rerun!")
ENDIF (IS_ABSOLUTE $ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Common mods directory where Fortran mod files are stored:
SET(CMAKE_Fortran_MODULE_DIRECTORY $ENV{COMMON_MODS})
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include)
# install(FILES ${HEADERS} DESTINATION include/mylib)

IF (IS_ABSOLUTE $ENV{COMMON_MODS})
  MESSAGE(STATUS "Mod target path is set!")
ELSE ($ENV{COMMON_MODS} LESS 1)
  MESSAGE(FATAL_ERROR "Env COMMON_MODS is not set, please set and rerun!")
ENDIF (IS_ABSOLUTE $ENV{COMMON_MODS})

#-------------------------------------------------------------------------------
## Get Fortran compiler name:
MESSAGE(STATUS "Compiler name: " ${CMAKE_Fortran_COMPILER})
MESSAGE(STATUS "Version is: " ${CMAKE_Fortran_COMPILER_VERSION})

# Currently PGI only: expand later:
IF (CMAKE_Fortran_COMPILER MATCHES "pgf")
  SET (CMAKE_Fortran_FLAGS_RELEASE "-D_REAL8_ -fopenmp -fopenmp")
  SET (CMAKE_Fortran_FLAGS_DEBUG "-D_REAL8_ -fopenmp -fopenmp")
ELSEIF (CMAKE_Fortran_COMPILER MATCHES "gfortran")
  IF (CMAKE_Fortran_COMPILER_VERSION LESS 10.0)
    SET (CMAKE_Fortran_FLAGS_RELEASE " -D_REAL8_ -g -fbounds-check -fopenmp")
    SET(CMAKE_Fortran_FLAGS_DEBUG " -D_REAL8_ -g -fbounds-check -DDEBUG -fopenmp")
  ELSE ()
    SET(CMAKE_Fortran_FLAGS_RELEASE " -D_REAL8_ -g -fbounds-check -fopenmp -fallow-argument-mismatch")
    SET(CMAKE_Fortran_FLAGS_DEBUG " -D_REAL8_ -g -fbounds-check -DDEBUG -fopenmp -fallow-argument-mismatch")
  ENDIF ()
  SET (CMAKE_Fortran_FLAGS "-DUSE_MPI -fpic -O3 -g -cpp -ffree-line-length-none -D_REAL8_ -fopenmp -fdec-math")
ELSEIF (CMAKE_Fortran_COMPILER MATCHES "ifort")
  SET (CMAKE_Fortran_FLAGS "-DUSE_MPI -fpic -O3 -g -cpp -D_REAL8_ -qopenmp")
ENDIF ()

# Double precision
# set(CMAKE_Fortran_Compiler,"ifort -m64")
# -ffree-form 
#-------------------------------------------------------------------------------

# set(CMAKE_LINKER ${CMAKE_Fortran_COMPILER})
# set(CMAKE_C_LINK_EXECUTABLE ${CMAKE_Fortran_COMPILER})
# set(CMAKE_CXX_LINK_EXECUTABLE ${CMAKE_Fortran_COMPILER})

# Add loac module path
SET (CMAKE_MODULE_PATH $ENV{CMAKETMP_DIR};${CMAKE_MODULE_PATH})

# Find packages 
# set(MPI_C_COMPILER "/opt/homebrew/Cellar/mpich/4.2.2/bin/mpicc")
# set(MPI_CXX_COMPILER "/opt/homebrew/Cellar/mpich/4.2.2/bin/mpicxx")
# set(MPI_Fortran_COMPILER "/opt/homebrew/Cellar/mpich/4.2.2/bin/mpifort")
find_package(MPI REQUIRED)
message(STATUS "MPI_Fortran_INCLUDE_PATH: ${MPI_Fortran_INCLUDE_PATH}")
INCLUDE_DIRECTORIES(SYSTEM ${MPI_Fortran_INCLUDE_PATH})

find_package(netcdf REQUIRED)
find_package(hdf5 REQUIRED)
find_package(Slint REQUIRED)

MESSAGE(STATUS "NETCDF_Fortran_INCLUDE_DIR is: " ${NETCDF_Fortran_INCLUDE_DIR})
MESSAGE(STATUS "HDF5HL_Fortran_LIBRARY is: " ${HDF5HL_Fortran_LIBRARY})
MESSAGE(STATUS "NETCDF_LIBRARY_Shared is: " ${NETCDF_LIBRARY_Shared})
MESSAGE(STATUS "NETCDF_LIBRARY_C_Shared is: " ${NETCDF_LIBRARY_C_Shared})

MESSAGE(STATUS "HDF5_HL_LIBRARY is: " ${HDF5_HL_LIBRARY})
MESSAGE(STATUS "HDF5_Fortran_LIBRARY is: " ${HDF5_Fortran_LIBRARY})
MESSAGE(STATUS "HDF5_LIBRARY is: " ${HDF5_LIBRARY})
MESSAGE(STATUS "NETCDF_LIBRARY is: " ${NETCDF_LIBRARY})
MESSAGE(STATUS "NETCDF_LIBRARY_C is: " ${NETCDF_LIBRARY_C})
MESSAGE(STATUS "CMAKE_MODULE_PATH is: " ${CMAKE_MODULE_PATH})
MESSAGE(STATUS "MPI_Fortran_INCLUDE_PATH is: " ${MPI_Fortran_INCLUDE_PATH})
MESSAGE(STATUS "MPI_Fortran_COMPILER_INCLUDE_DIRS is: " ${MPI_Fortran_COMPILER_INCLUDE_DIRS})
MESSAGE(STATUS "ZLIB_LIBRARY is: " ${ZLIB_LIBRARY})

MESSAGE(STATUS "SLINT_Fortran_INCLUDE_DIR is: " ${SLINT_Fortran_INCLUDE_DIR})
MESSAGE(STATUS "SLINT_LIBRARY is: " ${SLINT_LIBRARY})


# The subsequent is make for enabling the openmp for clang on macos
if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
  set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
  set(OpenMP_C_LIB_NAMES "omp")
  set(OpenMP_omp_LIBRARY omp)
  LINK_DIRECTORIES(/opt/homebrew/opt/libomp/lib)
endif()

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
  set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
  set(OpenMP_CXX_LIB_NAMES "omp")
  set(OpenMP_omp_LIBRARY omp)
  LINK_DIRECTORIES(/opt/homebrew/opt/libomp/lib)
endif()

find_package(OpenMP REQUIRED)

IF ( ${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
  find_package(rttov REQUIRED)
ELSE() 
  find_package(rttov)
ENDIF()

# Download and compile third party modules
# find_package(Git QUIET)
# if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# # Update submodules as needed
#     option(GIT_SUBMODULE "Check submodules during build" ON)
#     if(GIT_SUBMODULE)
#         message(STATUS "Submodule update")
#         execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --remote --recursive
#                         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#                         RESULT_VARIABLE GIT_SUBMOD_RESULT)
#         MESSAGE(STATUS ${GIT_EXECUTABLE})
#         if(NOT GIT_SUBMOD_RESULT EQUAL "0")
#             message(STATUS "git submodule update --init --remote --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
#         endif()
#     endif()
# endif()

# if(NOT EXISTS "${PROJECT_SOURCE_DIR}/yaml-cpp/CMakeLists.txt")
#     message(FATAL_ERROR "The submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
# endif()

# add_subdirectory(yaml-cpp EXCLUDE_FROM_ALL)

######## END the common setup:
