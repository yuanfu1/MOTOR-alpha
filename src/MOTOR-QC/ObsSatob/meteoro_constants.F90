!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Description:
! Define constants
! 1. fundamental constants are assigned directly,
! 2. Users-defined constants are specified declaration in this module
! and are assigned in namelist.
!
! History:
! Version Date Comment
! ------- ---- -------
! HISTORY           : Origionally from RCNMP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
MODULE Meteoro_constants

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! [1.0] Physical parameter constants
!-----------------------------------------------------------------------------

! Fundamental constants

  REAL, PARAMETER   :: pi = 3.1415926535897932346
  REAL, PARAMETER   :: pi180 = PI / 180.0
  REAL, PARAMETER   :: gravity_s = 9.80665
  REAL, PARAMETER   :: earth_rad = 6378.15E+3
  REAL, PARAMETER   :: omg = 7.292E-5
  REAL, PARAMETER   :: t_kelvin = 273.15
  REAL, PARAMETER   :: rd = 2.8705E+2
  REAL, PARAMETER   :: rv = 4.6150E+2
  REAL, PARAMETER   :: xlv = 2.50E+6
  REAL, PARAMETER   :: cp = 1004.64

! Saturation vapour pressure constants (Rogers & Yau, 1989) for calculating
! relative humidty from specific humidity

  REAL, PARAMETER   :: es_alpha = 6.112      ! in hPa
  REAL, PARAMETER   :: es_beta = 17.67
  REAL, PARAMETER   :: es_gamma = 243.5

! for all obs types
  INTEGER, PARAMETER :: DMA_flag = 1
  INTEGER, PARAMETER :: element_lib_flag = 2
  INTEGER, PARAMETER :: extreme_flag = 3
  INTEGER, PARAMETER :: inter_consist_flag = 4

  INTEGER, PARAMETER :: background_flag = 8
  INTEGER, PARAMETER :: thin_flag = 9
  INTEGER, PARAMETER :: blk_flag = 9
! for temp
  INTEGER, PARAMETER :: lapse_rat_flag = 5
  INTEGER, PARAMETER :: wind_shear_flag = 5
  INTEGER, PARAMETER :: vertical_inter_flag = 5
  INTEGER, PARAMETER :: invers_flag = 6
  INTEGER, PARAMETER :: hydro_flag = 7
! synop and airep
  INTEGER, PARAMETER :: time_consist_flag = 5
  INTEGER, PARAMETER :: time_persist_flag = 6

  INTEGER, PARAMETER :: num_flag = 9
  INTEGER, PARAMETER :: flag_missing = 999999999 ! flag intial value
!-----------------------------------------------------------------------------
  INTEGER, PARAMETER :: NlevP = 29
  INTEGER, PARAMETER :: NlevP_satob = 15
  INTEGER, PARAMETER :: NlevU_satob = 15
  INTEGER, PARAMETER :: NlevV_satob = 15

  REAL, SAVE :: P_err_satob(NlevP_satob) = (/0.7994, 0.9104, 0.7975, 0.7301, 0.7766, 0.7968, &
                                             0.7997, 0.9342, 0.8498, 0.7796, 0.6595, 0.5840, &
                                             0.4766, 0.3247, 0.2938/)

  REAL, SAVE :: U_err_satob(NlevU_satob) = (/2.0, 2.0, 2.0, 3.5, 4.3, 5.0, &
                                             5.0, 5.0, 5.0, 5.0, 5.0, 5.0, &
                                             5.0, 5.0, 5.7/)

  REAL, SAVE :: V_err_satob(NlevV_satob) = (/2.0, 2.0, 2.0, 3.5, 4.3, 5.0, &
                                             5.0, 5.0, 5.0, 5.0, 5.0, 5.0, &
                                             5.0, 5.0, 5.7/)

  REAL, SAVE :: p_mandat(NlevP) = (/1000.0, 962.5, 925.0, 887.5, 850.0, 800.0, 750.0, 700.0, &
                                    650.0, 600.0, 550.0, 500.0, 450.0, 400.0, 350.0, 300.0, 275.0, 250.0, 225.0, &
                                    200.0, 175.0, 150.0, 125.0, 100.0, 70.0, 50.0, 30.0, 20.0, 10.0/)

! Define  levels and pressure value, extreme value are assigned directly
  INTEGER, PARAMETER   :: max_mut_levels = 1000  ! Maximum levels for single observation
  INTEGER, PARAMETER :: maxstl = 21

  REAL :: stdp(maxstl) = (/1000., 925., 850., 700., 500., 400., 300., 250., &
                           200., 150., 100., 70., 50., 30., 20., 10., 7., 5., 3., 2., 1./)
  REAL :: stdh(maxstl) = (/300., 624., 1500., 3000., 5500., 7000., 9000., 10000., &
                           12000., 14000., 16500., 18500., 20000., 22000., 26000., &
                           30000., 33000., 36000., 39000., 42000., 48000./)

! Missing values and the index number of the quality control
  INTEGER, PARAMETER :: iobsmissing = -888888
  REAL, PARAMETER :: robsmissing = -888888.0
  REAL, PARAMETER :: inobsmissing = 999999.0
  REAL, PARAMETER ::  elvmissing = 9999.0
  REAL, PARAMETER :: threshold = 0.000001

  INTEGER            :: max_satob = 1200000   ! ... of satellite wind obs.

END MODULE Meteoro_Constants
