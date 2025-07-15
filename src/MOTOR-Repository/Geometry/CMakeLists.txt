cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(geometry)


# include_directories(SYSTEM ${MPI_Fortran_COMPILER_INCLUDE_DIRS})

#-------------------------------------------------------------------------------
#  To include another path to -I in cmake:
#  For references, please see Yuanfu's "Notes on cmake and ctest" word docs.
#-------------------------------------------------------------------------------
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

# Add Slint libraries.
INCLUDE_DIRECTORIES(${SLINT_Fortran_INCLUDE_DIR})

#-------------------------------------------------------------------------------
# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
# ADD_LIBRARY(my_project_name STATIC ${my_project_name_SRC})
ADD_LIBRARY(Geometry STATIC ${Geometry_SRC})
TARGET_LINK_LIBRARIES(Geometry Utility MGGen MPDD YAMLRead ${MPI_Fortran_LIBRARIES} ${SLINT_LIBRARY})

# 2. Build test executables:
#SET (EXE1 "Test_Geometry.exe")
#SET (EXE2 "Test_Halo.exe")

# 3. Add codes to the executables:
#ADD_EXECUTABLE(${EXE1} Test_Geometry.F90)
#ADD_EXECUTABLE(${EXE2} Test_ImmersedBoundary.F90)
TARGET_LINK_LIBRARIES(${EXE1} Geometry Utility MGGen MPDD ${MPI_Fortran_LIBRARIES})
#TARGET_LINK_LIBRARIES(${EXE2} Utility MGGen MPDD ${MPI_Fortran_LIBRARIES} Geometry)

#ADD_EXECUTABLE(${EXE2} Test_Halo.F90)
TARGET_LINK_LIBRARIES(${EXE2} Geometry Utility  MGGen MPDD ${MPI_Fortran_LIBRARIES} State2NC IOGrapes)

# 4. Add tests: (MIGRATED TO UnitTests/MOTOR-Repository/CMakeLists.txt)
# ADD_TEST(NAME Test_Geometry
#         COMMAND mpirun -np 4 --oversubscribe Test_Geometry.exe
#         WORKING_DIRECTORY $ENV{BIN_DIR})

# 5. Set test criteria: (MIGRATED TO UnitTests/MOTOR-Repository/CMakeLists.txt)
# SET_TESTS_PROPERTIES(Test_Geometry PROPERTIES PASS_REGULAR_EXPRESSION "Test passed!")
