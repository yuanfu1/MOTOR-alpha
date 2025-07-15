# For developer

There are some basic rule in developing the MOTOR-DA.

## Compilation with CMake

MOTOR-DA is compiled with CMake currently. These is a shared included CMake list files in `/ctemplates/CMakeCommon.CMake`, this file will be included by CMakeLists.txt in every sub modules. So the compilation options, env variables can be written in this files.

When adding a new module, a simple way is to copy the template in the `/src/!pendingDevelop/yourModule` and amend the file list, exe names to make it your modules. 

```CMake
cmake_minimum_required(VERSION 3.5)

# Change to your project:
# Note: this CMake builds the OO design of Multiscale Data Assimilation System:
PROJECT(yourModule)

# Include the common template:
IF (IS_ABSOLUTE $ENV{CTEMPLATES})
  MESSAGE(STATUS "CTEMPLATES path is set!")
ELSE ($ENV{CTEMPLATES} LESS 1)
  MESSAGE(FATAL_ERROR "Env CTEMPLATES is not set, please set and rerun!")
ENDIF (IS_ABSOLUTE $ENV{CTEMPLATES})
include($ENV{CTEMPLATES}/CMakeCommon.CMake)

#  To include another path to -I in CMake:
INCLUDE_DIRECTORIES($ENV{COMMON_MODS})

#  For libraries saved in a place, import them:
LINK_DIRECTORIES($ENV{COMMON_LIBS})

# Source codes to this build:
# SET (my_project_name_SRC f90-files)
# Or include a textfile for list of filenames if the file list is too long
# in srcFile.txt:
# SET (yourModule_SRC	
#      module1.F90 
#      module2.F90
#     )
include (srcFiles.txt)

# 1. Generate lib file:
ADD_LIBRARY(yourModule STATIC ${yourModule_SRC})
TARGET_LINK_LIBRARIES(yourModule mpdd geometry utility ${MPI_Fortran_LIBRARIES})

# 2. Build test executables:
SET (EXE1 "test_PS.exe")

# 3. Add codes to the executables:
ADD_EXECUTABLE(${EXE1} test_ym.F90)
TARGET_LINK_LIBRARIES(${EXE1} yourModule mpdd geometry utility ${MPI_Fortran_LIBRARIES})
```

Finally, add your new module to `CMakeLists.txt` in the MOTOR-DA root dir:

```CMake
# MOTOR-DA module list, 
# Each module can be compile independently, comment if you do not wish to compile
# any other modules.

# states
add_subdirectory(./src/states/state)
add_subdirectory(./src/states/State_Xm)

# IO
add_subdirectory(./src/IO/grapesIO)

# Tempplate:
add_subdirectory(yourModule)     # Add your module here.
```
## Git Configs

Use 

```
git config core.ignorecase false
```

 Set the git, in case the repository is file name case sensitive.

## Sphinx and Doxygen Docs

### Writing doxygen and sphinx

The `Doxygen` docs will be generated automatically. The `Sphinx` docs can be written in format of markdown, and put into the path of `etc/docs`. To finally add `.md` file in the compilation, the file name and path should be added to the script `etc/docs/modules.rst`.
### compilation
The docs in MOTOR-DA is compiled using the compileDocs.sh in the root directory. Once the compilation is success, the html docs can be found in the `build/sphinx_html` and `build/doxygen_html`.
## Consistency of type precision

A fixed computing precision which is independent of the compiler is important, thus it is necessary to use a pre-defined type among the whole repository. 
One should use the follow codes to define the INTEGER and REAL numbers in the program:

```fortran
USE kinds_m, ONLY: i_kind, r_kind, r_Double, i_llong

    INTEGER(i_kind) :: A
    REAL(r_kind) :: B
```
## Configuration files

All the configurations are written and invoked in one namelist file (.nl) in the STATIC_DIR. There are two default config files in this directory: test.nl and template.nl: test.nl is only used for unit test configuration; and template.nl a guided template for users to create new applications.
```
     |-- static                      % Namelist files 
        |-- test.nl 
        |-- template.nl
```
In the formal applications, the config files should be an input argument from the command line, e. g.:
```
./3DVar.exe grapesJob_01.nl
```

## Naming conventions

Use suffix `_m` to any module names.
Use suffic `_t` to any types.

Use CamelCase for any types.

Use snake_case for any member varible or functions.

Use g_ for global variables.

Use m_ for member variables.

## ALLOCATE and DEALLOCATE

Always finish the DEALLOCATE Sentence when ALLOCATE arrays:
```fortran
    IF(ALLOCATED(memory)) DEALLOCATE(memory) ! Check whether the memory is allocated.
```
## Always initialize your variables

USE:
```fortran
REAL(r_kind) :: vLevel = 0
```

RATHER THAN:
```fortran
REAL(r_kind) :: vLevel = 0
```
## Asynchronous

## OOP Principles

+ Open-Closed Principle,OCP
+ Liskov Substitution Principle,LSP
+ Dependence Inversion Principle,DIP
+ Composite/Aggregate Reuse Principle，CARP
+ Least Knowledge Principle
+ Interface Segregation Principle,ISP
+ Single Response Principle，SRP
  

[Link](https://blog.csdn.net/weixin_34235135/article/details/93624883)
