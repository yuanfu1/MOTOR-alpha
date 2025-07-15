cmake_minimum_required(VERSION 3.5)

PROJECT(ObsSatellite)

INCLUDE_DIRECTORIES($ENV{COMMON_MODS})
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

# Add Slint libraries.
INCLUDE_DIRECTORIES(${SLINT_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
# Source codes to this build:
include (srcFiles.txt)

# Add rttov libraries.
INCLUDE_DIRECTORIES(${RTTOV_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${RTTOV_MOD_DIR})

# Add HDF5 libraries.
INCLUDE_DIRECTORIES(${HDF5_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(ObsSatellite STATIC ${ObsSatellite_SRC})
TARGET_LINK_LIBRARIES(ObsSatellite ${NETCDF_LIBRARY} Geometry MPDD YAMLRead Field State ObsSet 
Mock State2NC Utility RTTOV ObsRadar ObsBase
${RTTOV_OTHER_LIBRARY} ${RTTOV_COEF_IO_LIBRARY} ${RTTOV_PARALLEL_LIBRARY} 
${RTTOV_HDF_LIBRARY} ${RTTOV_LIBRARY} ${RTTOV_EMIS_ATLAS_LIBRARY} 
${HDF5HL_Fortran_LIBRARY} ${HDF5_HL_LIBRARY} ${HDF5_Fortran_LIBRARY} ${HDF5_LIBRARY})

