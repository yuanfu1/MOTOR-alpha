! Description:
!> @file
!!   Defines routines to load geomagnetic field look-up table, and to compute
!!   the geomagnetic field and its angles relative to the wave propagation
!!   direction.
!
!> @brief
!!   Defines routines to load geomagnetic field look-up table, and to compute
!!   the geomagnetic field and its angles relative to the wave propagation
!!   direction.
!!
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
MODULE rttov_zutility

#include "throw.h"
  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : realtol

  IMPLICIT NONE

#include "rttov_errorreport.interface"

  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC :: load_bfield_lut   ! function to load Geomagnetic field LUT
  PUBLIC :: compute_bfield    ! subroutine to calculate the Geomagnetic field
                              ! using a LUT and two cosine angles of the field
                              ! relative to the wave propagation direction k.
  PUBLIC :: compute_kb_angles ! subroutine to compute two cosine angles of the
                              ! geomagnetic field relative to the wave
                              ! propagation direction k.

  INTERFACE compute_bfield
    MODULE PROCEDURE Compute_bfield_F1
    MODULE PROCEDURE compute_bfield_F2
    MODULE PROCEDURE compute_bfield_F3
  END INTERFACE compute_bfield

  ! Array for Earth's magnetic field
  INTEGER, PARAMETER   :: n_lat = 91, n_lon = 181
  INTEGER, SAVE        :: BField(3, n_lat, n_lon)
  REAL(jprb), PARAMETER :: DEGREES_TO_RADIANS  = 3.141592653589793238462643_JPRB/180.0_JPRB
CONTAINS

  !--------------------------------------------------------------------------------
  !
  ! NAME:
  !       load_bfield_lut
  !
  ! PURPOSE:
  !       Function to the geomagnetic field LUT
  !
  ! CALLING SEQUENCE:
  !   errorstatus = load_bfield_lut(filename_LUT)
  !
  ! INPUT ARGUMENT:
  !
  !   filename_LUT:       file name for the file containing the LUT of
  !                       the Earth magnetic field.
  !                       UNITS:      N/A
  !                       TYPE:       character string
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !----------------------------------------------------------------------------

  !> Load geomagnetic field LUT, returns errorstatus
  !! @param[in]     filename_LUT    LUT file name
  FUNCTION load_bfield_lut(filename_LUT) RESULT(err)
    CHARACTER(*), INTENT(IN)  :: filename_LUT
    INTEGER(jpim)             :: err

    CHARACTER(*), PARAMETER :: NameOfRoutine = 'rttov_load_bfield_lut '
    INTEGER(jpim) :: file_id
    LOGICAL       :: existence, is_open
    INTEGER(jpim) :: i, j, iskip

    TRY

    ! Find a free logical unit number for file access
    file_id = 9
    ! Start open loop for Lun Search
    lun_search: DO
      file_id = file_id + 1
      INQUIRE(UNIT = file_id, EXIST = existence)
      IF (.NOT. existence) THEN
        file_id = -1
        EXIT lun_search
      ENDIF
      INQUIRE(UNIT = file_id, OPENED = is_open)
      IF (.NOT. is_open) EXIT lun_search
    ENDDO lun_search

    OPEN(file_id, FILE = filename_LUT, STATUS = 'OLD', IOSTAT = err)
    THROWM(err .NE. 0,"Error opening "//TRIM(filename_LUT))

    READ(file_id, *, IOSTAT = err) iskip, iskip, &
                      ((BField(:, i, j), j = 1, n_lon), i = 1, n_lat)
    THROWM(err .NE. 0,"Error reading "//TRIM(filename_LUT))
    CLOSE(UNIT = file_id)

    CATCH
  END FUNCTION load_bfield_lut

  !--------------------------------------------------------------------------------
  !
  ! NAME:
  !       compute_bfield
  !
  ! PURPOSE:
  !       Subroutine to calculate the Geomagnetic field using a LUT and
  !       two cosine angles of the field relative to the wave propagation
  !       direction k.
  !
  ! CALLING SEQUENCE:
  !
  !        CALL compute_bfield(latitude, longitude, & ! Inputs
  !                            Bx, By, Bz, Be)        ! Outputs
  !
  !   OR
  !
  !        CALL compute_bfield(latitude,          &   ! input
  !                            longitude,         &   ! input
  !                            sensor_zenang,     &   ! input
  !                            sensor_aziang,     &   ! input
  !                            Be,                &   ! output
  !                            cos_bkang,         &   ! output
  !                            cos_baziang)           ! output
  !
  !   OR
  !
  !        CALL compute_bfield(latitude,               &   ! input
  !                            longitude,              &   ! input
  !                            sensor_zenang,          &   ! input
  !                            sensor_relative_aziang, &   ! input
  !                            Julian_day,             &   ! input
  !                            utc_time,               &   ! input
  !                            Be,                     &   ! output
  !                            cos_bkang,              &   ! output
  !                            cos_baziang)                ! output
  !
  !
  ! INPUT ARGUMENTS:
  !
  !
  !      latitude :       Latitude, -90 - 90 (90 North Pole; -90 - South Pole)
  !                       UNITS:      Degree
  !                       TYPE:       real(jprb)
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !
  !      longitude:       longitude, 0 - 360 East
  !                       UNITS:      Degree
  !                       TYPE:       real(jprb)
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !
  !   sensor_zenang:      sensor zenith angle
  !                       UNITS:      Degree
  !                       TYPE:       real(jprb)
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !
  !   sensor_aziang:      sensor azimuth angle (starts from North,
  !                                             positive clockwise)
  !                       UNITS:      Degree
  !                       TYPE:       real(jprb)
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !
  ! sensor_relative_aziang: sensor azimuth angle relative to the sun azimuth
  !                         angle
  !    Sun_azimuth_angle - from north, East positive
  !    sensor_relative_aziang = sun_azimuth_angle - Sensor_aziang
  !                       UNITS:      degree
  !                       TYPE:       real(jprb)
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !
  !        Julian_day:    Julian_Day 1=Jan 1, 365=Dec31 (366 leap year)
  !                       UNITS:      day
  !                       TYPE:       integer(JPIM)
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !
  !          utc_time:    Universal_Time 0.00-23.999,(GMT,Z time)
  !                       UNITS:      hour
  !                       TYPE:       real(jprb)
  !                       DIMENSION:  Scale
  !                       ATTRIBUTES: INTENT(IN)
  !
  ! OUTPUT ARGUMENTS:
  !       Bx:             Magetic field East component
  !                       UNITS:      Gauss
  !                       TYPE:       real(jprb)
  !                       DIMENSION:  Scalar
  !                       ATTRIBUTES: INTENT(OUT)
  !
  !       By:             Magetic field North component
  !                       UNITS:      Gauss
  !                       TYPE:       real(jprb)
  !                       DIMENSION:  Scalar
  !                       ATTRIBUTES: INTENT(OUT)
  !
  !       Bz:             Magetic field zenith component (positive upward)
  !                       UNITS:      Gauss
  !                       TYPE:       real(jprb)
  !                       DIMENSION:  Scalar
  !                       ATTRIBUTES: INTENT(OUT)
  !
  !       Be:             Magetic field strength (sqrt(BxBx + ByBy + BzBz))
  !                       UNITS:      Gauss
  !                       TYPE:       real(jprb)
  !                       DIMENSION:  Scalar
  !                       ATTRIBUTES: INTENT(OUT)
  !
  !     cos_bkang:        cosine of the angle between the magnetic field Be
  !                       vector and the wave propagation direction k.
  !                       UNITS:      N/A
  !                       TYPE:       real(jprb)
  !                       DIMENSION:  Scalar
  !                       ATTRIBUTES: INTENT(OUT)
  !
  !   cos_baziang:        cosine of the azimuth angle of the Be vector in the
  !                       (v, h, k) coordinates system, where v, h and k comprise
  !                       a right-hand orthogonal system, similar to the (x, y, z)
  !                       Catesian coordinates. The h vector is normal to the
  !                       plane containing the k and z vectors, where k points
  !                       to the wave propagation direction and z points
  !                       to the zenith. h = (z cross k)/|z cross k|. The
  !                       azimuth angle is the angle on the (v, h) plane
  !                       from the positive v axis to the projected line of the
  !                       Be vector on this plane, positive counterclockwise.
  !                       UNITS:      N/A
  !                       TYPE:       real(jprb)
  !                       DIMENSION:  Scalar
  !                       ATTRIBUTES: INTENT(OUT)
  !
  !--------------------------------------------------------------------------------

  !> Calculate geomagnetic field from LUT based on latitude and longitude
  !! @param[in]     latitude    latitude (-90 - 90, degrees)
  !! @param[in]     longitude   longitude (0 - 360, degrees)
  !! @param[out]    Bx          East-component of magnetic field (Gauss)
  !! @param[out]    By          North-component of magnetic field (Gauss)
  !! @param[out]    Bz          Zenith-component of magnetic field, positive upward (Gauss)
  !! @param[out]    Be          Magnetic field strength (Gauss)
  SUBROUTINE compute_bfield_F1(latitude, longitude, & ! Inputs
                               Bx, By, Bz, Be)        ! Outputs
    REAL(jprb), INTENT(IN)  :: latitude, longitude
    REAL(jprb), INTENT(OUT) :: Bx, By, Bz, Be

    REAL(jprb)                :: lat, lon
    REAL(jprb)                :: x2, w1_lat, w1_lon
    INTEGER                   :: idx_lat, idx_lon
    REAL(jprb), PARAMETER     :: dlat = 2.0_jprb, dlon = 2.0_jprb

    lat = 90.0_jprb - latitude  ! lat: 0 - 180 (0 North Pole; 180 South Pole)
    lon = longitude             ! lon: 0 - 360 East
    IF(lon < 0.0_jprb)lon = lon + 360.0_jprb

    idx_lat = INT(lat / dlat) + 1
    IF (idx_lat >= n_lat) idx_lat = n_lat - 1
    idx_lon = INT(lon / dlat) + 1
    IF (idx_lon >= n_lon) idx_lon = n_lon - 1

    x2 = REAL(idx_lat, jprb) * dlat
    w1_lat = (x2 - lat) / dlat

    x2 = REAL(idx_lon, jprb) * dlat
    w1_lon = (x2 - lon) / dlat

    Bx = BField_Component(1, Bfield, w1_lat, w1_lon, idx_lat, idx_lon)
    By = BField_Component(2, Bfield, w1_lat, w1_lon, idx_lat, idx_lon)
    Bz = BField_Component(3, Bfield, w1_lat, w1_lon, idx_lat, idx_lon)

    Be = SQRT(Bx*Bx + By*By + Bz*Bz)

  CONTAINS

    FUNCTION BField_Component(comp, Bfield, w1_lat, w1_lon, &
                              idx_lat, idx_lon) RESULT(B)

     INTEGER,    INTENT(IN)  :: comp
     INTEGER,    INTENT(IN)  :: BField(:,:,:)
     REAL(jprb), INTENT(IN)  :: w1_lat, w1_lon
     INTEGER,    INTENT(IN)  :: idx_lat, idx_lon

     REAL(jprb) :: B

     REAL(jprb), PARAMETER :: scale = 0.001_jprb
     REAL(jprb)            :: w2_lat, w2_lon

     w2_lat = 1.0_jprb - w1_lat
     w2_lon = 1.0_jprb - w1_lon
     B = (w1_lon*(w1_lat*REAL(BField(comp,idx_lat,   idx_lon),  jprb) + &
                  w2_lat*REAL(BField(comp,idx_lat+1, idx_lon),  jprb))  &
         + w2_lon*(w1_lat*REAL(BField(comp,idx_lat,   idx_lon+1),jprb) + &
                   w2_lat*REAL(BField(comp,idx_lat+1, idx_lon+1),jprb))) * scale

     END FUNCTION BField_Component

  END SUBROUTINE compute_bfield_F1

  !> Calculate geomagnetic field from LUT based on latitude and longitude and
  !! sensor geometry. Cos_baziang is the cosine of the azimuth angle of the
  !! Be vector in the (v, h, k) coordinates system, where v, h and k comprise
  !! a right-hand orthogonal system, similar to the (x, y, z) Catesian
  !! coordinates. The h vector is normal to the plane containing the k and z
  !! vectors, where k points to the wave propagation direction and z points
  !! to the zenith. h = (z cross k)/|z cross k|. The azimuth angle is the angle
  !! on the (v, h) plane from the positive v axis to the projected line of the
  !! Be vector on this plane, positive counterclockwise.
  !! @param[in]     latitude        latitude (-90 - 90, degrees)
  !! @param[in]     longitude       longitude (0 - 360, degrees)
  !! @param[in]     sensor_zenang   sensor zenith angle (degrees)
  !! @param[in]     sensor_aziang   sensor azimuth angle (degrees)
  !! @param[out]    Be              magnetic field strength (Gauss)
  !! @param[out]    cos_bkang       cosine of the angle between the magnetic field
  !!                                  Be vector and the wave propagation direction k
  !! @param[out]    cos_baziang     cosine of the azimuth angle of the  Be vector in
  !!                                  the (v, h, k) coordinates system (see above)
  SUBROUTINE compute_bfield_F2(latitude, longitude, sensor_zenang, sensor_aziang, &
                               Be, cos_bkang, cos_baziang)

    REAL(jprb), INTENT(IN)  :: latitude
    REAL(jprb), INTENT(IN)  :: longitude
    REAL(jprb), INTENT(IN)  :: sensor_zenang
    REAL(jprb), INTENT(IN)  :: sensor_aziang
    REAL(jprb), INTENT(OUT) :: Be
    REAL(jprb), INTENT(OUT) :: cos_bkang
    REAL(jprb), INTENT(OUT) :: cos_baziang

    REAL(jprb) :: Bx, By, Bz

    ! get Earth magnetic filed from LUT
    CALL compute_bfield_F1(latitude, longitude,  & ! inputs
                           Bx, By, Bz, Be)         ! outputs

    ! compute the cosines of the angle between the magnetic field Be and
    ! propagation direction and Be's azimuth angle
    CALL Compute_kb_Angles(Bx, By, Bz, sensor_zenang, sensor_aziang,  &  ! Input
                           cos_bkang, cos_baziang)                       ! Output

  END SUBROUTINE compute_bfield_F2

  !> Calculate geomagnetic field from LUT based on latitude and longitude, sensor
  !! geometry and time. Cos_baziang is the cosine of the azimuth angle of the
  !! Be vector in the (v, h, k) coordinates system, where v, h and k comprise
  !! a right-hand orthogonal system, similar to the (x, y, z) Catesian
  !! coordinates. The h vector is normal to the plane containing the k and z
  !! vectors, where k points to the wave propagation direction and z points
  !! to the zenith. h = (z cross k)/|z cross k|. The azimuth angle is the angle
  !! on the (v, h) plane from the positive v axis to the projected line of the
  !! Be vector on this plane, positive counterclockwise.
  !! @param[in]     latitude                 latitude (-90 - 90, degrees)
  !! @param[in]     longitude                longitude (0 - 360, degrees)
  !! @param[in]     sensor_zenang            sensor zenith angle (degrees)
  !! @param[in]     sensor_relative_aziang   sun azimuth angle minus sensor azimuth angle (degrees)
  !! @param[in]     julian_day               Julian day 1=Jan 1, 365=Dec31 (366 leap year)
  !! @param[in]     utc_time                 universal time 0.00-23.999 (GMT/Z time)
  !! @param[out]    Be                       magnetic field strength (Gauss)
  !! @param[out]    cos_bkang                cosine of the angle between the magnetic field
  !!                                           Be vector and the wave propagation direction k
  !! @param[out]    cos_baziang              cosine of the azimuth angle of the  Be vector in
  !!                                           the (v, h, k) coordinates system (see above)
  SUBROUTINE compute_bfield_F3(latitude, longitude, sensor_zenang, &
                               sensor_relative_aziang, julian_day, utc_time, &
                               Be, cos_bkang, cos_baziang)

    REAL(jprb),    INTENT(IN)  :: latitude
    REAL(jprb),    INTENT(IN)  :: longitude
    REAL(jprb),    INTENT(IN)  :: sensor_zenang
    REAL(jprb),    INTENT(IN)  :: sensor_relative_aziang
    INTEGER(jpim), INTENT(IN)  :: julian_day
    REAL(jprb),    INTENT(IN)  :: utc_time
    REAL(jprb),    INTENT(OUT) :: Be
    REAL(jprb),    INTENT(OUT) :: cos_bkang
    REAL(jprb),    INTENT(OUT) :: cos_baziang

    REAL(jprb)    :: Bx, By, Bz, Solar_ZA, Solar_AZ, lat, lon, sensor_aziang

    ! compute the cosines of the angle between the magnetic field Be and
    ! propagation direction and Be's azimuth angle
    CALL compute_bfield_F1(latitude, longitude,  & ! inputs
                            Bx, By, Bz, Be)        ! outputs

    lat = latitude
    lon = longitude
    IF (lon > 180.0_jprb) lon = lon - 360.0_jprb
    ! Compute Solar azimuth angle
    CALL Solar_ZA_AZ(lat, lon, REAL(julian_day, jprb), utc_time, &
                     Solar_ZA, Solar_AZ)
    ! Compute satellite azimuth angle (zero north, east +90)
    ! NB relative azimuth is measured as (sun_aziang - sensor_aziang)
    sensor_aziang = solar_AZ - sensor_relative_aziang

    ! compute for cos_bkangle, cos_baziangle
    CALL Compute_kb_Angles(Bx, By, Bz, sensor_zenang, sensor_aziang,  &  ! Input
                           cos_bkang, cos_baziang)                       ! Output

  END SUBROUTINE compute_bfield_F3


!--------------------------------------------------------------------------------
!
! NAME:
!       Compute_kb_Angles
!
! PURPOSE:
!       Subroutine to calculate the cosine of the angle between the Geomagnetic
!       field (Be) and wave propagation direction (k) and the cosine of the
!       azimuth angle of the Be vector in the (v, h, k) coordinates system (see
!       more detailed description below)
!
!
! CALLING SEQUENCE:
!      CALL Compute_kb_Angles(Bx, By, Bz,   &                       ! Input
!                             sensor_zenang, sensor_aziang, &       ! Input
!                             cos_bk_Angle, cos_baziang)            ! Output
! INPUT ARGUMENTS:
!
!       Bx:             Magetic field East component
!                       UNITS:      Gauss
!                       TYPE:       Real(jprb)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       By:             Magetic field North component
!                       UNITS:      Gauss
!                       TYPE:       Real(jprb)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Bz:             Magetic field zenith component (positive upward)
!                       UNITS:      Gauss
!                       TYPE:       Real(jprb)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!   sensor_zenang :     sensor zenith angle
!                       UNITS:      Degree
!                       TYPE:       Real(jprb)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!   sensor_aziang :     sensor zenith angle defined as the
!                       angle from North, positive clockwise.
!                       UNITS:      Degree
!                       TYPE:       Real(jprb)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!     cos_bkang:        cosine of the angle between the magnetic field Be
!                       vector and the wave propagation direction k.
!                       UNITS:      N/A
!                       TYPE:       real(jprb)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(OUT)
!
!   cos_baziang:        cosine of the azimuth angle of the Be vector in the
!                       (v, h, k) coordinates system, where v, h and k comprise
!                       a right-hand orthogonal system, similar to the (x, y, z)
!                       Catesian coordinates. The h vector is normal to the
!                       plane containing the k and z vectors, where k points
!                       to the wave propagation direction and z points
!                       to the zenith. h = (z cross k)/|z cross k|. The
!                       azimuth angle is the angle on the (v, h) plane
!                       from the positive v axis to the projected line of the
!                       Be vector on this plane, positive counterclockwise.
!                       UNITS:      N/A
!                       TYPE:       Real(jprb)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(OUT)
!
!--------------------------------------------------------------------------------

  !> Calculate the cosine of the angle between the geomagnetic field (Be) and
  !! wave propagation direction (k) and the cosine of the azimuth angle of the
  !! Be vector in the (v, h, k) coordinates system, where v, h and k comprise
  !! a right-hand orthogonal system, similar to the (x, y, z) Catesian
  !! coordinates. The h vector is normal to the plane containing the k and z
  !! vectors, where k points to the wave propagation direction and z points
  !! to the zenith. h = (z cross k)/|z cross k|. The azimuth angle is the angle
  !! on the (v, h) plane from the positive v axis to the projected line of the
  !! Be vector on this plane, positive counterclockwise.
  !! @param[in]     Bx              East-component of magnetic field (Gauss)
  !! @param[in]     By              North-component of magnetic field (Gauss)
  !! @param[in]     Bz              Zenith-component of magnetic field, positive upward (Gauss)
  !! @param[in]     sensor_zenang   sensor zenith angle (degrees)
  !! @param[in]     sensor_aziang   sensor azimuth angle (degrees)
  !! @param[out]    cos_bkang       cosine of the angle between the magnetic field
  !!                                  Be vector and the wave propagation direction k
  !! @param[out]    cos_baziang     cosine of the azimuth angle of the  Be vector in
  !!                                  the (v, h, k) coordinates system (see above)
  SUBROUTINE Compute_kb_Angles(Bx, By, Bz,   &                       ! Input
                               sensor_zenang, sensor_aziang, &       ! Input
                               cos_bkang, cos_baziang)   ! Output
    REAL(jprb), INTENT(IN)  :: Bx, By, Bz
    REAL(jprb), INTENT(IN)  :: sensor_zenang, sensor_aziang
    REAL(jprb), INTENT(OUT) :: cos_bkang, cos_baziang

    ! Local
    REAL(jprb) :: B, B_v, B_h, B_p, kx, ky, kz, azi, &
                  SIN_SenZA, COS_SenAZ, SIN_SenAZ, COS_SenZA

    ! Input azimuth is zero north, east +90. The code below expects zero east, north +90
    azi = 90._jprb - sensor_aziang

    SIN_SenZA = SIN(sensor_zenang*DEGREES_TO_RADIANS)
    COS_SenZA = COS(sensor_zenang*DEGREES_TO_RADIANS)
    SIN_SenAZ = SIN(azi*DEGREES_TO_RADIANS)
    COS_SenAZ = COS(azi*DEGREES_TO_RADIANS)

    ! compute k directional vector from satellite's zenith and azimuth angles
    kx = SIN_SenZA * COS_SenAZ
    ky = SIN_SenZA * SIN_SenAZ
    kz = COS_SenZA

    ! compute consine of the angle between the magnetic field B and k
    B = SQRT(bx*bx + by*by + bz*bz)
    cos_bkang = (kx*bx + ky*by + kz*bz) / B

    ! Project the B vector on the V and H plane: B_v - B component on the V
    ! axis; B_h - B component on the H axis.
    B_v = bx*kx*kz + by*ky*kz - bz*(ky*ky + kx*kx)
    B_h = -bx*ky + by*kx

    ! compute the cosine of the azimuth angle
    B_p = SQRT(B_v*B_v + B_h*B_h)
    IF (ABS(B_p) > realtol) THEN
      cos_baziang = B_v / B_p
    ELSE
      cos_baziang = 0.0   ! not defined (take an arbitrary value)
    ENDIF

  END SUBROUTINE Compute_kb_Angles

  !> Calculate solar zenith and azimuth angles for a given place and time
  !! @param[in]     latitude        latitude (-90 - 90, degrees)
  !! @param[in]     longitude       longitude (0 - 360, degrees)
  !! @param[in]     julian_day      Julian day 1=Jan 1, 365=Dec31 (366 leap year)
  !! @param[in]     universal_time  universal time 0.00-23.999 (GMT/Z time)
  !! @param[out]    za              solar zenith angle
  !! @param[out]    az              solar azimuth angle
  SUBROUTINE Solar_ZA_Az(latitude,        &
                         longitude,       &
                         julian_day,      &
                         universal_time,  &
                         ZA,              &
                         Az)
!
!************************************************************************
!*                                                                      *
!*  Module Name:  Solar_Az_Za                                           *
!*                                                                      *
!*  Language:  Fortran    Library:                                      *
!*  Version.Rev:  1.0  22 Feb 91  Programmer: Kleespies                 *
!*      1.1  28 Feb 91  Kleespies                                       *
!*           Put equation of time into hour angle.                      *
!*
!*      updated to f95, 4, September 2007, Y. Han                       *
!*                                                                      *
!*  Calling Seq:    Call Solar_Az_ZA(                                   *
!*     &        latitude,                                               *
!*     &        longitude,                                              *
!*     &        julian_day,                                             *
!*     &        universal_time,                                         *
!*     &        ZA,                                                     *
!*     &        Az)                                                     *
!*                                                                      *
!*                                                                      *
!*  Description:  Computes solar azimuth and zenith angle for a         *
!*    given place and time.  Solar azimuth is the angle                 *
!*    positive east of north form a point to the center of              *
!*    the solar disc. Solar zenith angle is the angle                   *
!*    in degrees from zenith to the center of the solar disk.           *
!*    The algorithms are taken from "Introduction to Solar              *
!*    Radiation" by Muhamad Iqbal, Academic Press, 1983.                *
!*    Note that lat=0,lon=0, 21Mar 12z will not give sun                *
!*    overhead, because that is not the way things work.                *
!*                                                                      *
!*  Input Args:  R*4  Latitude, +NH, -SH degrees                        *
!*      R*4  Longitude,-W, +E                                           *
!*      R*4  Julian_Day 1=Jan 1, 365=Dec31 (366 leap yr)                *
!*      R*4  Universal_Time 0.00-23.99,(GMT,Z time)                     *
!*                                                                      *
!*  Output Args:  R*4  ZA  Solar Zenith Angle                           *
!*      R*4  AZ  Solar Azmuth Angle                                     *
!*                                                                      *
!*  Common Blks:  none                                                  *
!*  Include:  none                                                      *
!*  Externals:  none                                                    *
!*  Data Files:  none                                                   *
!*                                                                      *
!*  Restrictions:  Accurate to within .1 deg.                           *
!*      No checking made to the validity of input                       *
!*      parameters.                                                     *
!*      Since solar zenith angle is a conic angle,                      *
!*      it has no sign.                                                 *
!*      No correction made for refraction.                              *
!*                                                                      *
!*  Error Codes:  none                                                  *
!*                                                                      *
!*  Error Messages:                                                     *
!*                                                                      *
!************************************************************************
!
        IMPLICIT NONE

        REAL(jprb), INTENT(IN)  :: latitude,longitude,julian_day,universal_time
        REAL(jprb), INTENT(OUT) :: ZA, Az

        REAL(jprb) :: local_sun_time,solar_elevation,equation_of_time
        REAL(jprb) :: cosza,cosaz
        REAL(jprb) :: hour_angle,day_angle
        REAL(jprb) :: solar_declination
        REAL(jprb) :: rlatitude
        REAL(jprb), PARAMETER :: one_eighty_over_pi  = 1.0/DEGREES_TO_RADIANS
        REAL(jprb), PARAMETER :: threesixty_over_24 = 15.0
        REAL(jprb), PARAMETER ::  threesixty_over_365 = 0.98630137
        REAL(jprb), PARAMETER ::  min_declination = -23.433
        REAL(jprb), PARAMETER ::  day_offset = 10.0  !  original equation had this nine

!*  Compute day angle

        day_angle = threesixty_over_365 * (julian_day - 1.) * DEGREES_TO_RADIANS

        rlatitude  = latitude * DEGREES_TO_RADIANS

!*  Compute equation of Time

        Equation_of_Time =  &
               ( 0.000075   &
               + 0.001868 * COS(Day_Angle)  &
               - 0.032077 * SIN(Day_Angle)  &
               - 0.014615 * COS(2. * Day_Angle) &
               - 0.040890 * SIN(2. * Day_Angle)) * 229.18 / 60 ! in hours

!*  Compute local sun time
        local_sun_time  = universal_time  &
                       + Equation_of_Time &
                       + longitude / threesixty_over_24

!*  Compute solar declination

        solar_declination = &
         (       0.006918  &
               - 0.399912 * COS(day_angle)      &
               + 0.070257 * SIN(day_angle)      &
               - 0.006758 * COS(2. * day_angle)  &
               + 0.000907 * SIN(2. * day_angle)  &
               - 0.002697 * COS(3. * day_angle)  &
               + 0.001480 * SIN(3. * day_angle))

!*  Compute hour angle
        hour_angle = threesixty_over_24 * MOD(local_sun_time + 12.0_jprb, 24.0_jprb)

!*  Compute solar zenith angle
        cosza = SIN(rlatitude) * SIN(solar_declination)  &
              + COS(rlatitude) * COS(solar_declination) * COS(hour_angle*DEGREES_TO_RADIANS)

        ZA = ACOS(cosza) * one_eighty_over_pi

!*  Compute solar azimuth angle
        solar_elevation = 90.0_jprb - ZA

        IF (ABS(ZA) < realtol) Then
          Az = 180.0    ! handle arbitrary case
        ELSE
          cosaz = (SIN(solar_elevation*DEGREES_TO_RADIANS) * SIN(rlatitude) - &
                   SIN(solar_declination)) / &
                  (COS(solar_elevation*DEGREES_TO_RADIANS) * COS(rlatitude))

          IF (cosaz .LT. -1.0_jprb) cosaz = -1.0
          IF (cosaz .GT.  1.0_jprb) cosaz =  1.0

          Az = ACOS(cosaz) * one_eighty_over_pi

!    The above formula produces azimuth positive east, zero south.
!    We want positive east, zero north.

!
          IF (Az .GE. 0.0_jprb) THEN
              Az = 180.0 - Az
          ELSE
              Az = -180.0 + Az
          ENDIF

          IF (hour_angle .LT. 180.0) Az = - Az
          IF (Az .LT. 0) Az = 360.0 + Az

        ENDIF

   END SUBROUTINE Solar_ZA_Az

END MODULE rttov_zutility
