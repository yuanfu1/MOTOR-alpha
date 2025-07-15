!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-PS.PlotAState.F90
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring 
!                     Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie on 2025-02-07
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This test program tests a mixed solver of Psi and Chi with U and V boundary conditions.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
PROGRAM PlotAState
  USE kinds_m, ONLY: i_kind, r_kind
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE geometry_m, ONLY: geometry_t
  USE State_m, ONLY: State_t
  USE State2NC_m
  USE YAMLRead_m

  TYPE(mpddGlob_t) :: mpddGlob
  TYPE(geometry_t) :: geo
  CHARACTER(LEN=1024) :: yamlFile, ncOutputFile
  CHARACTER(LEN=20) :: task
  INTEGER(i_kind) :: glevel,vlevel
  REAL(r_kind), ALLOCATABLE :: data(:,:)
  TYPE(State_t) :: X

  ! Get the configuration file
  CALL getarg(1, yamlFile)

  ! Set up MPDD and Geometry:
  CALL mpddGlob%initialize()
  CALL geo%initialize(yamlFile,mpddGlob)

  ! Plot the solutions:
  istatus = yaml_get_var(TRIM(yamlFile), 'IO', 'output_dir', ncOutputFile)
  istatus = yaml_get_var(TRIM(yamlFile), 'RunMode', 'Task', task)

  ! Read the data file:
  glevel = geo%mg%mg_finest; vlevel = geo%mg%sg(glevel)%vLevel
  CALL X%initialize(yamlFile,geo%mg%sg(glevel))
  READ(10) X%fields(1)%DATA(:,:,1)
  READ(10) X%fields(2)%DATA(:,:,1)
  
  CALL Output_NC_State_AV(X, TRIM(ncOutputFile), TRIM(task)//"_PlotAState", .TRUE., .TRUE.)

END PROGRAM PlotAState