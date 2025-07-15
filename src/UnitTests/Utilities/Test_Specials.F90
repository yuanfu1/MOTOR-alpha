!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.TestMixedSolver.F90
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring 
!                     Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie on 2025-03-18
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This test program tests the harmonic functions, currently supporting spherical and plane harmonic.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
PROGRAM Test_Specials
  USE harmonicBase_m
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int

  USE mpddGlob_m, ONLY: mpddGlob_t
  USE geometry_m, ONLY: geometry_t
  USE gzm_m, ONLY: gzm_t
  USE State_m, ONLY: State_t
  USE State2NC_m
  USE YAMLRead_m
  USE parameters_m, ONLY: EarthRadius, pi

  IMPLICIT NONE
  
  TYPE(sphericalHarmonics_t), ALLOCATABLE :: sphericalObj
  TYPE(planeHarmonics_t), ALLOCATABLE :: planeObj
  CHARACTER(LEN=1024) :: yamlFile, ncOutputFile
  CHARACTER(LEN=20) :: task
  INTEGER(c_int) :: i, l, m, istatus
  REAL(c_double) :: theta, phi, x, y
  REAL(c_double) :: sph_result, plane_result

  ! MPDD and Geometry:
  TYPE(mpddGlob_t) :: mpddGlob
  TYPE(geometry_t) :: geo
  TYPE(gzm_t) :: gzm
  TYPE(State_t) :: state

  ! Get the yaml file:
  CALL getarg(1, yamlFile)
  IF (TRIM(yamlFile) .EQ. '') THEN
    PRINT*,'Usage: <this executable> <yaml file>, please and find and supply your yamlfile and rerun'
    STOP
  ELSE
    PRINT*,'yaml: ',TRIM(yamlFile)
  END IF

  ! Initialize test parameters for simple test suggested by Copilot:
  l = 2
  m = 1
  theta = 0.5D0  ! Polar angle in radians
  phi = 1.0D0    ! Azimuthal angle in radians
  x = 1.0D0      ! x-coordinate for plane harmonic
  y = 1.0D0      ! y-coordinate for plane harmonic

  CALL mpddGlob%initialize()
  CALL geo%initialize(yamlFile,mpddGlob)

  ! Initialize a gzm model for using its Laplacian operator:
  ASSOCIATE(sg => geo%mg%sg(geo%mg%mg_finest), mgEnd => geo%mg%mg_finest)
    gzm = gzm_t(sg)

    ! Test the functions:
    istatus = yaml_get_var(TRIM(yamlFile), 'IO', 'output_dir', ncOutputFile)
    istatus = yaml_get_var(TRIM(yamlFile), 'RunMode', 'Task', task)
    CALL state%initialize(yamlFile,sg)

    ! Create and test spherical harmonic object
    ALLOCATE(sphericalObj)
    ALLOCATE(planeObj)

    ! Obtain the spherical harmonic object values at grid cell center:
    DO i=1,sg%num_cell
        theta = sg%cell_cntr(1,i)
        phi = sg%cell_cntr(2,i)
        state%fields(1)%DATA(1,i,1) = sphericalObj%harmonicFunction(l, m, theta, phi)*EarthRadius
    END DO
    CALL gzm%Laplacia(state%fields(1)%DATA(:,:,1), state%fields(2)%DATA(:,:,1))
    WRITE(*,1) MAXVAL(ABS(state%fields(2)%DATA(1,:,1))),sg%glevel
1   FORMAT('Max Laplacian of Spherical Harmonic: ',E14.6,' at G:',I1)
    CALL Output_NC_State_AV(state, TRIM(ncOutputFile), TRIM(task)//"harmonics", .TRUE., .TRUE.)
  END ASSOCIATE

  ! Create and test plane harmonic object
!   plane_result = planeObj%harmonicFunction(l, m, x, y)
!   PRINT*, 'Plane Harmonic Result (l=2, x=1.0, y=1.0): ', plane_result

!   ! Test sine function (existing test)
!   PRINT*, 'Testing sine function:'
!   PRINT*, 'sin(0.5) = ', DSIN(0.5D0)

  ! Clean up
  DEALLOCATE(sphericalObj, planeObj)
END PROGRAM Test_Specials