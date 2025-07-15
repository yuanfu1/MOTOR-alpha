cmake_minimum_required(VERSION 3.5)
project(gzm_upd)

# Include directories
include_directories($ENV{COMMON_MODS})
include_directories(${CMAKE_SOURCE_DIR}/src/MOTOR-PS/PreSys)
include_directories(${CMAKE_BINARY_DIR}/src/MOTOR-PS/gzm_upd)  # Ensure the module files are found
include_directories(${CMAKE_SOURCE_DIR}/src/MOTOR-PS/poissonsolver_kp)


# Link directories
link_directories($ENV{COMMON_LIBS})


# Source files
set(gzm_upd_SRC
    ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/gzm_tlm.F90
    ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/gzm_adj.F90
    ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/Objective_Function.F90
    ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/Matrix_Utils.F90
    ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/CalVerDer_TL.F90
    ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/CalVerDer_adj.F90
)

# ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/psMatrix_adj.F90
# ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/psSolutn_adj.F90
# ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/namelist_ps_adj.F90
# ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/poissonSolver_adjoint.F90
# Generate library
add_library(gzm_upd STATIC ${gzm_upd_SRC})
target_link_libraries(gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid )

# Build executables
set(EXE1 "verify_divergen_adjoint.exe")
set(EXE2 "verify_jacobian_adjoint.exe")
set(EXE3 "verify_laplacia_adjoint.exe")
# set(EXE4 "Test_CalVerDer_AD_without_TL.exe")
# set (EXE5 "verify_adj_calverder_generic.exe")
#set (EXE6 "unit_test_gzm_tlm_ad_grid.exe")
#set (EXE7 "test_gzm_tl_ad_consistency.exe")
#set (EXE8 "test_adjoint_operator.exe")
#set (EXE9 "taylor_expansion_test.exe")
# set (EXE10 "test_adjoint_without_tl.exe")
# set (EXE11 "test_grid_point_adjoint_without_tl.exe")
set (EXE12 "verify_adj_calverder_firstorder.exe")
set (EXE13 "verify_adj_calverder_secondorder.exe")
set (EXE14 "verify_adj_calverder_firstorder_half.exe")
set (EXE15 "verify_grad_lat_adjoint.exe")
set (EXE16 "verify_grad_lon_adjoint.exe")
# set (EXE17 "test_PS_Adjoint.exe")
#set (EXE15 "verify_adj_precals.exe")


add_executable(${EXE1} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/verify_divergen_adjoint.F90)
target_link_libraries(${EXE1} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

 add_executable(${EXE2} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/verify_jacobian_adjoint.F90)
 target_link_libraries(${EXE2} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

add_executable(${EXE3} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/verify_laplacia_adjoint.F90) 
target_link_libraries(${EXE3} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

# add_executable(${EXE4} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/Test_CalVerDer_AD_without_TL.F90)
# target_link_libraries(${EXE4} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

# add_executable(${EXE5} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/verify_adj_calverder_generic.F90)
# target_link_libraries(${EXE5} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

#add_executable(${EXE6} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/unit_test_gzm_tlm_ad_grid.F90)
#target_link_libraries(${EXE6} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

#add_executable(${EXE7} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/test_gzm_tl_ad_consistency.F90)
#target_link_libraries(${EXE7} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

#add_executable(${EXE8} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/test_adjoint_operator.F90)
#target_link_libraries(${EXE8} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

#add_executable(${EXE9} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/taylor_expansion_test.F90)
#target_link_libraries(${EXE9} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

# add_executable(${EXE10} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/test_adjoint_without_tl.F90)
# target_link_libraries(${EXE10} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

# add_executable(${EXE11} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/test_grid_point_adjoint_without_tl.F90)
# target_link_libraries(${EXE11} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)
# add_executable(${EXE8} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/taylor_expansion_test.F90)
# target_link_libraries(${EXE7} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)
# Ensure the module files are created in the binary directory

 add_executable(${EXE12} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/verify_adj_calverder_firstorder.F90)
 target_link_libraries(${EXE12} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

add_executable(${EXE13} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/verify_adj_calverder_secondorder.F90)
target_link_libraries(${EXE13} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

 add_executable(${EXE14} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/verify_adj_calverder_firstorder_half.F90)
 target_link_libraries(${EXE14} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

add_executable(${EXE15} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/verify_grad_lat_adjoint.F90)
target_link_libraries(${EXE15} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

add_executable(${EXE16} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/verify_grad_lon_adjoint.F90)
target_link_libraries(${EXE16} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

# add_executable(${EXE17} ${CMAKE_SOURCE_DIR}/src/MOTOR-PS/gzm_upd/test_PS_Adjoint.F90)
# target_link_libraries(${EXE17} gzm_upd MPDD Geometry Utility YAMLRead ${MPI_Fortran_LIBRARIES} PreSys zgrid)

set_target_properties(gzm_upd PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/src/MOTOR-PS/gzm_upd)
