cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(GrapesIO)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Source codes to this build:

# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

# #-------------------------------------------------------------------------------
# # Generate lib file:
ADD_LIBRARY(GrapesIO STATIC ${GrapesIO_SRC} )
TARGET_LINK_LIBRARIES(GrapesIO Geometry Utility NMLRead ${MPI_Fortran_LIBRARIES})
 
# 2. Build test executables:
#SET (EXE1 "Test_GrapesIO.exe")
#SET (EXE2 "Test_qcqrIO.exe")
#SET (EXE3 "Test_GrapesModelvarIO.exe")
#SET (EXE4 "Test_GrapesModelvarIO_500m.exe")

# 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE1} Test_GrapesIO.F90)
# TARGET_LINK_LIBRARIES(${EXE1} Utility GrapesIO ${MPI_Fortran_LIBRARIES} IOGrapes)
# ADD_EXECUTABLE(${EXE2} Test_qcqrIO.F90)
# TARGET_LINK_LIBRARIES(${EXE2} Utility GrapesIO ${MPI_Fortran_LIBRARIES} IOGrapes)
# # 2. Build test executables:
# SET (EXE1 "Test_GrapesIO.exe")

# SET (EXE3 "Test_GrapesPostvarIO.exe")

#3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_GrapesIO.F90)
# TARGET_LINK_LIBRARIES(${EXE1} Utility GrapesIO ${MPI_Fortran_LIBRARIES} IOGrapes)

# ADD_EXECUTABLE(${EXE3} Test_GrapesModelvarIO.F90)
# TARGET_LINK_LIBRARIES(${EXE3} Utility GrapesIO ${MPI_Fortran_LIBRARIES} IOGrapes)

#ADD_EXECUTABLE(${EXE4} Test_GrapesModelvarIO_500m.F90)
# TARGET_LINK_LIBRARIES(${EXE4} Utility GrapesIO ${MPI_Fortran_LIBRARIES} IOGrapes)

# ADD_EXECUTABLE(${EXE3} Test_GrapesPostvarIO.F90)
# TARGET_LINK_LIBRARIES(${EXE3} Utility GrapesIO ${MPI_Fortran_LIBRARIES} IOGrapes)
