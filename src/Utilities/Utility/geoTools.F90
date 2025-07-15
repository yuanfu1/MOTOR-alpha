!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/4/12, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the manipulating functions in geo-coordinate calculations.
!! @author Zilong Qin
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
MODULE geoTools_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE parameters_m

  TYPE GeoBox_t
    REAL(r_kind) :: maxLon ! In radius
    REAL(r_kind) :: minLon
    REAL(r_kind) :: maxLat
    REAL(r_kind) :: minLat
  CONTAINS
    PROCEDURE :: inBox
    PROCEDURE, NOPASS :: inBoxElem
  END TYPE GeoBox_t

  INTERFACE GeoBox_t
    PROCEDURE :: constructor
  END INTERFACE GeoBox_t

CONTAINS

  FUNCTION inBox(this, ll)
    CLASS(GeoBox_t) :: this
    REAL(r_kind) :: ll(2)     !< Lat/Lon
    LOGICAL :: inBox

    inBox = .TRUE.
    IF (ll(1) > this%maxLat) THEN
      inBox = .FALSE.
      RETURN
    END IF

    IF (ll(1) < this%minLat) THEN
      inBox = .FALSE.
      RETURN
    END IF

    IF (ll(2) > this%maxLon) THEN
      inBox = .FALSE.
      RETURN
    END IF

    IF (ll(2) < this%minLon) THEN
      inBox = .FALSE.
      RETURN
    END IF
  END FUNCTION inBox

  PURE ELEMENTAL FUNCTION inBoxElem(lat, lon, maxLat, minLat, maxLon, minLon)
    REAL(r_kind), INTENT(IN) :: lat, lon, maxLat, minLat, maxLon, minLon     !< Lat/Lon
    LOGICAL :: inBoxElem

    inBoxElem = .TRUE.
    IF (lat > maxLat) THEN
      inBoxElem = .FALSE.
      RETURN
    END IF

    IF (lat < minLat) THEN
      inBoxElem = .FALSE.
      RETURN
    END IF

    IF (lon > maxLon) THEN
      inBoxElem = .FALSE.
      RETURN
    END IF

    IF (lon < minLon) THEN
      inBoxElem = .FALSE.
      RETURN
    END IF
  END FUNCTION

  FUNCTION constructor(maxLat, minLat, maxLon, minLon) RESULT(this)
    TYPE(GeoBox_t) :: this
    REAL(r_kind), INTENT(IN) :: maxLat, minLat, maxLon, minLon
    this%maxLat = maxLat
    this%minLat = minLat
    this%maxLon = maxLon
    this%minLon = minLon
  END FUNCTION

!> @brief
!! Calculate the distance of two points on the Earth sphere.
  FUNCTION distance(locSrc, locDst)
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: locSrc(2), &        !> in lat/lon radius
                                locDst(2)           !> in lat/lon radius
    REAL(r_kind) :: distance

    distance = EarthRadius * DACOS(DSIN(locSrc(1)) * DSIN(locDst(1)) + &
                                   DCOS(locSrc(1)) * DCOS(locDst(1)) * DCOS(locSrc(2) - locDst(2)))
  END FUNCTION distance

!> @brief   BY PJM
!! Calculate the distance of two points on the Earth sphere via Haversine formula.
  FUNCTION distance_hav(lat1, lon1, lat2, lon2)
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: lat1, lon1, &       !> in radius
                                lat2, lon2          !> in radius
    REAL(r_kind) :: dlat, dlon, a
    REAL(r_kind) :: distance_hav

    ! lat1 = locSrc(1)
    ! lon1 = locSrc(2)
    ! lat2 = locDst(1)
    ! lon2 = locDst(2)
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    a = (DSIN(dlat / 2))**2 + DCOS(lat1) * DCOS(lat2) * (DSIN(dlon / 2))**2
    distance_hav = 2 * EarthRadius * DASIN(DSQRT(a))
  END FUNCTION distance_hav

!> @brief   BY PJM
!! Calculate the distance of two points in 3D.
  FUNCTION distance_3d(lat1, lon1, hgt1, lat2, lon2, hgt2)
    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: lat1, lon1, &       !> in radius
                                lat2, lon2, &       !> in radius
                                hgt1, hgt2          !> height, in meters
    REAL(r_kind) :: distance_2d
    REAL(r_kind) :: distance_3d

    distance_2d = distance_hav(lat1, lon1, lat2, lon2)
    distance_3d = DSQRT(distance_2d**2 + (hgt2 - hgt1)**2)

  END FUNCTION distance_3d

!> @brief   BY PJM
!! Calculate the locolization coeff between two point via Gaspari_Cohn scheme.
  ! FUNCTION loccoeffs_GC(lat1,lon1,hgt1, lat2,lon2,hgt2, lozRadius)
  FUNCTION loccoeffs_GC(dist, lozRadius)
    IMPLICIT NONE
    ! REAL(r_kind), INTENT(IN) :: lat1, lon1, &       !> in radius
    !                             lat2, lon2, &       !> in radius
    !                             hgt1, hgt2, &       !> height, in meters, &
    REAL(r_kind), INTENT(IN) :: dist, &       !> distance between two points, in meters
                                lozRadius           !> radius of locolization, in meters
    REAL(r_kind) :: a, r1
    REAL(r_kind) :: loccoeffs_GC

    ! dist = distance_3d(lat1,lon1,hgt1, lat2,lon2,hgt2)
    a = lozRadius * DSQRT(10.0D0 / 3.0D0)
    r1 = dist / a
    IF (dist .GE. 0.0D0 .AND. dist .LE. a) THEN
      loccoeffs_GC = -1.0D0 / 4.0D0 * (r1**5.0D0) + 1.0D0 / 2.0D0 * (r1**4.0D0) + 5.0D0 / 8.0D0 * (r1**3.0D0) - 5.0D0 / 3.0D0 * (r1**2.0D0) + 1.0D0
    ELSE IF (dist .GT. a .AND. dist .LE. a * 2.0D0) THEN
      loccoeffs_GC = 1.0D0 / 12.0D0 * (r1**5.0D0) - 1.0D0 / 2.0D0 * (r1**4.0D0) + 5.0D0 / 8.0D0 * (r1**3.0D0) + 5.0D0 / 3.0D0 * (r1**2.0D0) &
                     - 5.0D0 * r1 + 4.0D0 - 2.0D0 / 3.0D0 * (r1**(-1.0D0))
    ELSE
      loccoeffs_GC = 0.0D0
    END IF

  END FUNCTION loccoeffs_GC

END MODULE geoTools_m
