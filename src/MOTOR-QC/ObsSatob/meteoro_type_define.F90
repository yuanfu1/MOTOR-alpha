MODULE Meteoro_type_define
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Description:
! Definitions of derived data types

! History:
! Version Date Comment
! ------- ---- -------
! HISTORY           : Origionally from RCNMP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !use obs_size_constants
  USE Meteoro_constants

  IMPLICIT NONE

  INTEGER :: yyyy0, mon0, dd0, hh0           ! assimilation or QC time  == QC_date
!-----------------------------------------------------------------------------
! [1.0] Observation structure definition:
!-----------------------------------------------------------------------------
  TYPE type_level_info
    INTEGER            :: level_id      !  id for level
    REAL                 :: time_d ! displacement  information for time
    REAL                 :: lat_d  ! displacement  information for latitude
    REAL                 :: lon_d  ! displacement  information for longitude
    REAL                 :: wind_shear
    REAL                 :: humid_interp
    REAL                 :: d_p
    REAL                 :: T_interp
    REAL                 :: T_lapse      ! temperature lapse rate
    REAL                 :: T_invers     ! inverse temperature
    REAL                 :: dthick  ! thick_o- thick_c
  END TYPE type_level_info
!--------------------------------------------------------------------------------
!
  TYPE type_obs_field
    REAL    :: VALUE       ! Observation
    REAL    :: corret_value       ! corret_value
    REAL    :: corret_mothod      ! corret_mothod
    REAL    :: error       ! Observational error
    INTEGER :: flag        ! Observational quality identifier
    REAL    :: o_b          ! o_b
    REAL    :: o_a         ! o_a
    REAL    :: bias      ! the bias correction
  END TYPE type_obs_field
!--------------------------------------------------------------------------------
  TYPE type_header_info
    CHARACTER(LEN=30)  :: name          ! Station name
    CHARACTER(LEN=9)   :: call_name     ! call of ships or plane
    INTEGER               :: itime(5)      ! yyyy-mm-dd-hh-mm
    REAL                  :: stn_id        ! 5 digit station identifer
    REAL                  :: obs_type_id   !
    REAL                  :: lat           ! Latitude in degree
    REAL                  :: lon           ! Longitude in degree
    INTEGER               :: position      ! the lat lon check flag
    TYPE(type_obs_field) :: elv           ! Elevation in m
    TYPE(type_obs_field) :: ps            ! Surface pressure
    INTEGER               :: levels        !(for temp it means number of levels)
    INTEGER               :: instru_code   ! instrument
    INTEGER               :: radia_correct ! the radiation correction method for temp
    REAL                  :: baseline_bias ! the first thick bias
    REAL                  :: QI            ! the qi index for satob
  END TYPE type_header_info

!---------------------------------------------------------------

  TYPE type_level_data
    TYPE(type_obs_field) :: p ! Pressure in Pa
    TYPE(type_obs_field) :: u ! Wind x-component in m/s
    TYPE(type_obs_field) :: v ! Wind y-component in m/s
    TYPE(type_obs_field) :: dd ! Wind direction
    TYPE(type_obs_field) :: ff ! Wind speed
    TYPE(type_obs_field) :: rh ! relative humidity in %
    TYPE(type_obs_field) :: phi ! Height in gm
    TYPE(type_obs_field) :: T ! Temperature in K
    TYPE(type_obs_field) :: Td !
    TYPE(type_obs_field) :: Vtan ! Tangent wind for bda
  END TYPE type_level_data
!-----------------------------------------------------------------
  TYPE type_ob

    INTEGER :: total_obs, &
               num_cosmic, num_synop, num_ships, &
               num_airep, num_satob, num_modis
    INTEGER :: num_temp, num_surx, num_pilot

    TYPE(type_temp), POINTER    :: temp(:)
    TYPE(type_cosmic), POINTER :: cosmic(:)
    TYPE(type_synop), POINTER  :: synop(:)
    TYPE(type_airep), POINTER  :: airep(:)
    TYPE(type_ships), POINTER  :: ships(:)
    TYPE(type_satob), POINTER  :: satob(:)
    TYPE(type_buoy), POINTER  :: buoy(:)
    TYPE(type_pilot), POINTER  :: pilot(:)
  END TYPE type_ob
!-------------------------------------------------------------------------
  TYPE type_temp
    TYPE(type_level_info)  :: level_info(max_mut_levels)
    TYPE(type_header_info) :: header_info
    TYPE(type_level_data)  :: level_data(max_mut_levels)
  END TYPE type_temp
!------------------------------------------------------------------------
!  follow for frist guess check

  TYPE type_y_temp
    TYPE(type_level_data) :: level_data(max_mut_levels)
  END TYPE type_y_temp

!------------------
  TYPE type_pilot
    TYPE(type_level_info)  :: level_info(max_mut_levels)
    TYPE(type_header_info) :: header_info
    TYPE(type_level_data)  :: level_data(max_mut_levels)
    INTEGER                 :: corect_method
    REAL                    :: sea_T
    INTEGER                 :: trace_method
    INTEGER                 :: p_levs, h_levs
  END TYPE type_pilot
!------------------------------------------------------------------------
!  follow for frist guess check

  TYPE type_y_pilot
    TYPE(type_level_data) :: level_data(max_mut_levels)

  END TYPE type_y_pilot
!-----------------
  TYPE type_y_synop
    TYPE(type_level_data) :: level_data
    TYPE(type_obs_field)   :: Pmsl
  END TYPE type_y_synop

  TYPE type_cosmic
    TYPE(type_header_info) :: header_info
    TYPE(type_level_data)  :: level_data(max_mut_levels)
  END TYPE type_cosmic
!---------------------------------------------------
  TYPE type_synop
    TYPE(type_header_info) :: header_info
    TYPE(type_level_data)  :: level_data
    TYPE(type_obs_field)   :: Pmsl
  END TYPE type_synop
!---------------------------------------------------
  TYPE type_surx
    TYPE(type_header_info) :: header_info
    TYPE(type_level_data)  :: level_data
  END TYPE type_surx

!----------------------------------------------------
  TYPE type_airep
    TYPE(type_header_info) :: header_info
    TYPE(type_level_data)  :: level_data
    INTEGER :: flight_state
  END TYPE type_airep

  TYPE type_y_airep
    TYPE(type_level_data)  :: level_data
  END TYPE type_y_airep
!----------------------------------------------------
  TYPE type_ships
    TYPE(type_header_info) :: header_info
    TYPE(type_level_data)  :: level_data

  END TYPE type_ships

  TYPE type_y_ships
    TYPE(type_level_data) :: level_data
  END TYPE type_y_ships

!-------------------------------------------------------
  TYPE type_buoy
    TYPE(type_header_info) :: header_info
    TYPE(type_level_data)  :: level_data
  END TYPE type_buoy
!--------------------------------------------------------

  TYPE type_satob
    TYPE(type_header_info) :: header_info
    TYPE(type_level_data)  :: level_data
  END TYPE type_satob

  TYPE type_y_satob
    TYPE(type_level_data)   :: level_data
  END TYPE type_y_satob

!-------------------------------------------------------------
!---------------------------------------------------------------
!=============================================================================
!---------------------------------------------------------------
![2.0] Background structure definition:
!---------------------------------------------------------------
!2.1 Background field
  TYPE type_xb
    INTEGER :: mdate(5) ! date of background field
    INTEGER :: miy ! y: 1st dimension of background grid.
    INTEGER :: mjx ! x: 2nd dimension of background grid.
    INTEGER :: mkp ! p: 3rd dimension of background grid.

    REAL :: dy ! grid mesh in y direction
    REAL, POINTER :: dx(:) ! grid mesh in x direction
    REAL, POINTER :: glat(:) !latitude
    REAL, POINTER :: glon(:) !longitude
    REAL, POINTER :: cori(:) !Corilis parameter f(y)

    REAL, POINTER :: u(:, :, :) ! u component of wind
    REAL, POINTER :: v(:, :, :) ! v component of wind
    REAL, POINTER :: phi(:, :, :)! geopotential height g*z
    REAL, POINTER :: t(:, :, :) ! temperature
    REAL, POINTER :: rh(:, :, :) ! relative humidity
    REAL, POINTER :: q(:, :, :) ! specific humidity
    REAL, POINTER  :: qs(:, :, :) ! saturated specific humidity

    REAL, POINTER :: ps(:, :) ! surface pressure
    REAL, POINTER :: ts(:, :) ! surface temperature
    REAL, POINTER :: q2m(:, :) ! 2m specific humidity
    REAL, POINTER :: t2m(:, :) ! 2m temperature
    REAL, POINTER :: rh2m(:, :) ! 2m relative humidity
    REAL, POINTER :: u10m(:, :)
    REAL, POINTER :: v10m(:, :)
  END TYPE type_xb
!-------------------------
!   2.2 Background field
!----------------------
  TYPE type_xberr
    REAL, POINTER  :: u(:, :, :)  ! u component of wind
    REAL, POINTER  :: v(:, :, :)  ! v component of wind
    REAL, POINTER  :: phi(:, :, :)! geopotential height g*z
    REAL, POINTER  :: t(:, :, :)  ! temperature
    REAL, POINTER  :: rh(:, :, :) ! relative humidity
  END TYPE type_xberr

END MODULE Meteoro_type_define

