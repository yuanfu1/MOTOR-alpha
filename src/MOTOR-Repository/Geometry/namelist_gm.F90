!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2020/12/18, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This is a module reading in all needed namelist files.
!! @author Zilong Qin, Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
MODULE namelist_gm_m

  USE kinds_m, ONLY: i_kind, r_kind
  !USE NMLRead_m, ONLY: namelist_read
  USE YAMLRead_m
  USE AdvanceTime_m, ONLY: Time_GMT_to_Unix

CONTAINS

  SUBROUTINE namelist_gm(configFile, proc_layout, t_steps_mg, start_time_uTime_sec, end_time_uTime_sec, &
                         mg_coarsest, mg_finest, dimCell_global, vLevel)
    ! Model namelist variables:
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER(i_kind)                       :: fu, rc
    INTEGER(i_kind), INTENT(OUT) :: dimCell_global(2)        !< Dimension of grid cells in Lat-lon grid
    INTEGER(i_kind), INTENT(in) :: &
      mg_finest, &
      mg_coarsest
    INTEGER(i_kind), INTENT(OUT) :: &
      proc_layout(mg_coarsest:mg_finest, 2), t_steps_mg(mg_coarsest:mg_finest)
    REAL(r_kind), INTENT(OUT) :: start_time_uTime_sec, end_time_uTime_sec

    INTEGER(i_kind), ALLOCATABLE :: mpi_layout_g14(:), mpi_layout_g58(:), mpi_layout_g9L(:), &
                                    num_grid(:), domain_latlon(:), start_time(:), end_time(:)
    INTEGER(i_kind) :: time_steps_g14, time_steps_g58, time_steps_g9L
    ! Yuanfu Xie adds a new time step variable for reading default time step setting from a yaml
    INTEGER(i_kind), ALLOCATABLE :: time_steps_all(:)
    INTEGER(i_kind), ALLOCATABLE :: mpi_layout_all_dim1(:), mpi_layout_all_dim2(:)

    INTEGER(i_kind) :: mpi_layout(20, 2), t_steps_mg_file(20)
    INTEGER(i_kind) :: uTime_sec !start_time(6), end_time(6), uTime_sec

    ! INTEGER(i_kind) :: num_grid(2)         ! Numbers of grid points in latlon directions
    ! REAL(r_kind) :: domain_latlon(2, 2)         ! The domain range in degrees from the origin
    INTEGER(i_kind), INTENT(INOUT) :: vLevel
    INTEGER(i_kind) :: i, istatus

    !NAMELIST /geometry/ mpi_layout_g14, &
    !  mpi_layout_g58, &
    !  mpi_layout_g9L, &
    !  time_steps_g14, &
    !  time_steps_g58, &
    !  time_steps_g9L

    !NAMELIST /latlon_grid/ num_grid, domain_latlon
    !NAMELIST /analysis_para/ start_time, end_time

    ! Check whether file exists.
    INQUIRE (file=TRIM(configFile), iostat=rc)

    ! Read the namelist file
    IF (rc .EQ. 0) THEN
      !OPEN (10, FILE=configFile)
      !READ (10, NML=geometry)
      !CLOSE (10)
      istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mpi_layout_g14', mpi_layout_g14)
      istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mpi_layout_g58', mpi_layout_g58)
      istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mpi_layout_g9L', mpi_layout_g9L)

      istatus = yaml_get_var(TRIM(configFile), 'geometry', 'time_steps_g14', time_steps_g14)
      istatus = yaml_get_var(TRIM(configFile), 'geometry', 'time_steps_g58', time_steps_g58)
      istatus = yaml_get_var(TRIM(configFile), 'geometry', 'time_steps_g9L', time_steps_g9L)

      ! Yuanfu Xie adds a new time step variable for reading default time step setting from a yaml
      IF (ALLOCATED(time_steps_all)) THEN
        WRITE (*, 2) UBOUND(time_steps_all, 1), (time_steps_all(i), i=1, UBOUND(time_steps_all, 1))
2       FORMAT('Time step setting for ', I2, ' levels starting from 1: ', 20I4)
      END IF

      !OPEN (10, FILE=configFile)
      !READ (10, NML=latlon_grid)
      !CLOSE (10)
      istatus = yaml_get_var(TRIM(configFile), 'latlon_grid', 'num_grid', num_grid)
      istatus = yaml_get_var(TRIM(configFile), 'latlon_grid', 'domain_latlon', domain_latlon)

      !OPEN (10, FILE=configFile)
      !READ (10, NML=analysis_para)
      !CLOSE (10)
      istatus = yaml_get_var(TRIM(configFile), 'analysis_para', 'start_time', start_time)
      istatus = yaml_get_var(TRIM(configFile), 'analysis_para', 'end_time', end_time)

    ELSE
      ERROR STOP "Namelist file of *geometry* does not exist!"
    END IF

    FORALL (i=1:4) mpi_layout(i, :) = mpi_layout_g14
    FORALL (i=5:8) mpi_layout(i, :) = mpi_layout_g58
    FORALL (i=9:20) mpi_layout(i, :) = mpi_layout_g9L

    FORALL (i=1:4) t_steps_mg_file(i) = time_steps_g14
    FORALL (i=5:8) t_steps_mg_file(i) = time_steps_g58
    FORALL (i=9:20) t_steps_mg_file(i) = time_steps_g9L

    ! Yuanfu Xie adds a new time step variable for reading default time step setting from a yaml
    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'time_steps_all', time_steps_all)
    IF (istatus == 0) THEN
      IF (UBOUND(time_steps_all, 1) .GT. 20) WRITE (*, 3) UBOUND(time_steps_all, 1)
3     FORMAT('+------------------------------------------------------------+', /, &
             '| Warning the time steps are set for time frames over 20: ', I3, ' |', /, &
             '+------------------------------------------------------------+')
      IF (UBOUND(time_steps_all, 1) .GT. 0) t_steps_mg_file(1:UBOUND(time_steps_all, 1)) = time_steps_all(:)
    END IF

    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mpi_layout_all_dim1', mpi_layout_all_dim1)
    IF (istatus == 0) THEN
      IF (UBOUND(mpi_layout_all_dim1, 1) .GT. 20) WRITE (*, 4) UBOUND(mpi_layout_all_dim1, 1)
4     FORMAT('+------------------------------------------------------------+', /, &
             '| Warning the mpi layout are set for time frames over 20: ', I3, ' |', /, &
             '+------------------------------------------------------------+')
      IF (UBOUND(mpi_layout_all_dim1, 1) .GT. 0) mpi_layout(1:UBOUND(mpi_layout_all_dim1, 1), 1) = mpi_layout_all_dim1(:)
    END IF

    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mpi_layout_all_dim2', mpi_layout_all_dim2)
    IF (istatus == 0) THEN
      IF (UBOUND(mpi_layout_all_dim2, 1) .GT. 20) WRITE (*, 5) UBOUND(mpi_layout_all_dim2, 1)
5     FORMAT('+------------------------------------------------------------+', /, &
             '| Warning the mpi layout are set for time frames over 20: ', I3, ' |', /, &
             '+------------------------------------------------------------+')
      IF (UBOUND(mpi_layout_all_dim2, 1) .GT. 0) mpi_layout(1:UBOUND(mpi_layout_all_dim2, 1), 2) = mpi_layout_all_dim2(:)
    END IF

    proc_layout(mg_coarsest:mg_finest, :) = mpi_layout(mg_coarsest:mg_finest, :)
    t_steps_mg(mg_coarsest:mg_finest) = t_steps_mg_file(mg_coarsest:mg_finest)

    dimCell_global = num_grid + 1

    !CALL namelist_read(configFile, 'vLevel', vLevel)
    istatus = yaml_get_var(TRIM(configFile), 'modelState', 'vLevel', vLevel)

    ! Yuanfu Xie changed output style and added more info on the output:
    WRITE (*, 1) start_time(1:5), end_time(1:5)
1   FORMAT('namelist_mg - Analysis window:', /, &
           '    start time: ', I4, '-', I2.2, '-', I2.2, ' @', I2.2, ':', I2.2, /, &
           '    end   time: ', I4, '-', I2.2, '-', I2.2, ' @', I2.2, ':', I2.2)

    CALL Time_GMT_to_Unix(start_time, uTime_sec)
    start_time_uTime_sec = uTime_sec
    CALL Time_GMT_to_Unix(end_time, uTime_sec)
    end_time_uTime_sec = uTime_sec
  END SUBROUTINE namelist_gm

END MODULE namelist_gm_m
