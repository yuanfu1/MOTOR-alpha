cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
PROJECT(LBFGSB)


# Add netcdf libraries.
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})
LINK_DIRECTORIES(${NETCDF_LIBRARY})

#-------------------------------------------------------------------------------
# Source codes to this build:
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(LBFGSB STATIC ${LBFGSB_SRC})

# # 2. Build test executables:
# SET (EXE1 "drver1.exe")
# # SET (EXE2 "test_rotate.exe")

# # 3. Add code to the executables:
# ADD_EXECUTABLE(${EXE1} driver1.f90)
# TARGET_LINK_LIBRARIES(${EXE1} LBFGSB)
