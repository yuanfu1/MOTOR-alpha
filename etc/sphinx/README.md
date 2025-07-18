# Meteorological Object-oriented Tools and Operators Repository - Data Assimilation (MOTOR-DA)

## Module list of MOTOR-DA so far

| module | version | developer | package dependencies| other dependencies| support llmesh | support itmesh|
| :-----:| :----:  | :----: | :----: | :----: | :----: | :----: |
| utility | 0.0 | Yuanfu Xie | | | Yes | Yes |
| mgGen | 0.0 | Yuanfu Xie | utility | | Yes | Yes |
| mpdd | 0.0 | Zilong Qin | utility | MPI | Yes | No |
| geometry | 0.1 | Zilong Qin | utility, mpdd, mggen| MPI| Yes | No |
| poissonSolver_Kp | 0.1 | Zilong Qin | utility, mpdd, geometry, mggen| MPI| Yes | No |

## Basic rules to developer

### Version control and repository

MOTOR-DA use the repositories on coding.net temporarily:

`https://e.coding.net/opensimi/MOTOR-DA/MOTOR-DA.git` [[Click to open]](https://e.coding.net/opensimi/MOTOR-DA/MOTOR-DA.git)

### Hierarchy of directories                     
```
root:/
     |-- build                       % Build directory   
     |-- templates                   % Common used cmake file
     |-- etc                         % Docs and so on.
         |-- docs                    % Doxygen files
         |-- historicalCode
     |-- grid                        % Grid files
     |-- libs                        % Static libraries
     |-- mods                        % .mod files generated by .F90   
     |-- src                         % Source files
     |-- static                      % Namelist files 
```
### How to build this project

Under the root path:
1. `source pathEnv.sh`
Add the static variables to the system env.

2. `cd ./build`
Direct the current path to the build directory.

3. `cmake ..`
Configure the project, and generate the makefile with CMake tools.

4. `make`
Build the libs and programs with makefile.

### Consistency of type precision

A fixed computing precision independent of compiler is importrant, thus it is necessory to use a pre-defined type among the whole repository. 
One should use the follow codes to define the INTEGER and REAL numbers in the program:

```
USE kinds_m, ONLY: i_kind, r_kind, r_Double, i_llong

    INTEGER(i_kind) :: A
    REAL(r_kind) :: B
```

### Configuration files

All the configurations are written and invoked in one namelist file (.nl) in the STATIC_DIR. There are two default config files in this directory: test.nl and template.nl: test.nl is only used for unit test configuration; and template.nl a guided template for users to create new applications.
```
     |-- static                      % Namelist files 
        |-- test.nl 
        |-- template.nl
```
In the formal applications, the config files should be a input arguement from the command line, e. g.:
```
./3DVar.exe grapesJob_01.nl
```

## Introduction for each module
### mgGen

mgGen is a grid generator written by Dr. Yuanfu Xie, which provides the unconstructed mesh struct of both regional lat-lon grid and global icosaheral-triangle grid. 

### geometry

geometry is a tool module which initialzes a base class providing a series of operators for basic data manipulation on the model/anlysis space. The tool is optimized for 2D unconstructed mesh grid. And the key points below are supported:

+ Unconstructed grid indices in horizontal direction;
+ Multigrid is natively supported;
+ MPI is supported, and high-level parallel operators are enclosed;
+ The 2D grid field has one boundary cell surranding.

### mpdd (Multiprocess for domain decomposition)

mpdd is a tool module contains the low-level manipulations required by geometry.

### poissonSolver_Kp

poissonSolver_Kp is a multigrid solver for solving stream and potential functions with vorticity and divergence and Dirichlet boundary condition.


