# This script is built for finding the netcdf library automatically.

find_path(NETCDF_Fortran_INCLUDE_DIR netcdf.mod $ENV{INCLUDE_NETCDF_FORTRAN})

find_library(NETCDF_LIBRARY NAMES libnetcdff.a PATHS $ENV{LIBRARY_NETCDF_FORTRAN})
find_library(NETCDF_LIBRARY_C NAMES libnetcdf.a PATHS $ENV{LIBRARY_NETCDF_FORTRAN})
find_library(NETCDF_LIBRARY_Shared NAMES libnetcdf.dylib PATHS $ENV{LIBRARY_NETCDF_FORTRAN})
find_library(NETCDF_LIBRARY_C_Shared NAMES libnetcdff.dylib PATHS $ENV{LIBRARY_NETCDF_FORTRAN})
find_library(NETCDF_LIBRARY_Shared NAMES libnetcdf.so PATHS $ENV{LIBRARY_NETCDF_FORTRAN})
find_library(NETCDF_LIBRARY_C_Shared NAMES libnetcdff.so PATHS $ENV{LIBRARY_NETCDF_FORTRAN})

if (NETCDF_Fortran_INCLUDE_DIR AND NETCDF_LIBRARY AND NETCDF_LIBRARY_C
    AND NETCDF_LIBRARY_Shared AND NETCDF_LIBRARY_C_Shared)
    # Set the variables for the netcdf library
    set(NETCDF_FOUND TRUE)
endif(NETCDF_Fortran_INCLUDE_DIR AND NETCDF_LIBRARY AND NETCDF_LIBRARY_C
    AND NETCDF_LIBRARY_Shared AND NETCDF_LIBRARY_C_Shared)
