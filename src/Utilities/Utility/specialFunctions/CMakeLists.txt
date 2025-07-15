cmake_minimum_required(VERSION 3.1)

# Change to your project:
# build my project:
PROJECT(specialFunctions)

INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

# Add netcdf libraries.
INCLUDE_DIRECTORIES(${NETCDF_Fortran_INCLUDE_DIR})
LINK_DIRECTORIES(${NETCDF_LIBRARY})

# Add netcdf libraries.
# INCLUDE_DIRECTORIES(${SLINT_Fortran_INCLUDE_DIR})
# LINK_DIRECTORIES(${SLINT_LIBRARY})
#-------------------------------------------------------------------------------
#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

#-------------------------------------------------------------------------------
# Source codes to this build:
include (srcFiles.txt)

#-------------------------------------------------------------------------------
# Generate lib file:
ADD_LIBRARY(specialFunctions STATIC ${specialFunctions_SRC})
TARGET_LINK_LIBRARIES(specialFunctions PUBLIC lapack blas ${NETCDF_LIBRARY} ${NETCDF_LIBRARY_Shared} ${NETCDF_LIBRARY_C_Shared})

# 2. Build test executables:
#SET (EXE1 "Test_Specials.exe")

# 3. Add code to the executables:
#ADD_EXECUTABLE(${EXE1} Test_Specials.F90)
#  LINK_DIRECTORIES(dir1 dir2 ...)
#-------------------------------------------------------------------------------
#  To link another lib like -l in cmake:
#-------------------------------------------------------------------------------
#  TARGET_LINK_LIBRARIES(target <list of libnames>)
# TARGET_LINK_LIBRARIES(${EXE1} specialFunctions Utility Geometry zgrid Utility ${NETCDF_LIBRARY})

# # 4. Add tests:
# ADD_TEST(NAME specials
#         COMMAND Test_Specials.exe 
#         WORKING_DIRECTORY $ENV{BIN_DIR})

# # 5. Set test criteria:
# SET_TESTS_PROPERTIES(specials PROPERTIES FAIL_REGULAR_EXPRESSION "Test failed")
# SET_TESTS_PROPERTIES(specials PROPERTIES PASS_REGULAR_EXPRESSION "Test passed")
