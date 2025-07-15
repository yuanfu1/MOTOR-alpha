cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(IOGrapes)


#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})
# INCLUDE_DIRECTORIES("/opt/hpc/software/mathlib/netcdf/4.4.1/gcc/include/") 
# INCLUDE_DIRECTORIES(${SLINT_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})
# LINK_DIRECTORIES("/opt/hpc/software/mathlib/netcdf/4.4.1/gcc/lib")

# Add netcdf libraries.
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})
LINK_DIRECTORIES(${NETCDF_LIBRARY})

#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

# #-------------------------------------------------------------------------------
# # Generate lib file:
ADD_LIBRARY(IOGrapes STATIC ${IOGrapes_SRC})
TARGET_LINK_LIBRARIES(IOGrapes PUBLIC 
GrapesIO Utility MPDD State IOModel YAMLRead InterpHP ModelCoupler M2CV ${MPI_Fortran_LIBRARIES} YAMLRead MGOpts)

# # 2. Build test executables:
# SET (EXE2 "Test_IOGrapes.exe")
# SET (EXE3 "Test_IOGrapesModelvar.exe")
# SET (EXE4 "Test_IOGrapesModelvar_Ens.exe")
# SET (EXE5 "Test_IOGrapesPostvar.exe")
#SET (EXE6 "Test_IOGrapesModelvar_500m_to_nc.exe")

# # 3. Add codes to the executables:
# ADD_EXECUTABLE(${EXE2} Test_IOGrapes.F90)
# TARGET_LINK_LIBRARIES(${EXE2} Utility GrapesIO MPDD Geometry IOGrapes State State2NC ${MPI_Fortran_LIBRARIES} MGOpts)
# ADD_EXECUTABLE(${EXE3} Test_IOGrapesModelvar.F90)
# TARGET_LINK_LIBRARIES(${EXE3} Utility GrapesIO MPDD Geometry IOGrapes State YAMLRead State2NC ${MPI_Fortran_LIBRARIES} MGOpts)
# ADD_EXECUTABLE(${EXE4} Test_IOGrapesModelvar_Ens.F90)
# TARGET_LINK_LIBRARIES(${EXE4} Utility GrapesIO MPDD Geometry IOGrapes State YAMLRead State2NC ${MPI_Fortran_LIBRARIES} MGOpts)
# ADD_EXECUTABLE(${EXE5} Test_IOGrapesPostvar.F90)
# TARGET_LINK_LIBRARIES(${EXE5} Utility GrapesIO MPDD Geometry IOGrapes State YAMLRead State2NC ${MPI_Fortran_LIBRARIES} MGOpts)
#ADD_EXECUTABLE(${EXE6} Test_IOGrapesModelvar_500m_to_nc.F90)
# TARGET_LINK_LIBRARIES(${EXE6}  Utility GrapesIO MPDD Geometry IOGrapes State YAMLRead State2NC PostProc  ${MPI_Fortran_LIBRARIES} MGOpts)
