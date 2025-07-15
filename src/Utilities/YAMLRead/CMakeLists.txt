cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
PROJECT(YAMLRead)

# #-------------------------------------------------------------------------------
# # Source codes to this build:
# include (srcFiles.txt)
LINK_DIRECTORIES($ENV{COMMON_LIBS})
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

# #-------------------------------------------------------------------------------
# # Generate lib file:
ADD_LIBRARY(YAMLRead STATIC YAMLRead.F90)
TARGET_LINK_LIBRARIES(YAMLRead PUBLIC Utility fortran-yaml-cpp OpenMP::OpenMP_CXX yaml-cpp stdc++)

# 2. Build test executables:
#SET (EXE1 "Test_YAMLRead.exe")

# 3. Add code to the executables:
#ADD_EXECUTABLE(${EXE1} Test_YAMLRead.F90)
# TARGET_LINK_LIBRARIES(${EXE1} YAMLRead Utility fortran-yaml-cpp)
