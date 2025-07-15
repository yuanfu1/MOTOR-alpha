! Description:
!> @file
!!   Calculates information related to the viewing and solar geometry.
!
!> @brief
!!   Calculates information related to the viewing and solar geometry.
!!
!! @details
!!   Angle-related calculations are done here. There are no
!!   associated TL/AD/K calculations for these.
!!
!!   The input profiles(:)%zenangle is the geometrical (no refraction)
!!   zenith angle to satellite at surface
!!
!!   The calculated angles(:)%viewang is the geometrical (no refraction)
!!   nadir angle to surface at satellite: viewang is calculated from
!!   zenangle by application of the sine rule.
!!
!! @param[in]     plane_parallel  flag for strict plane parallel geometry
!! @param[in]     profiles        profiles structure
!! @param[in]     coef            rttov_coef structure
!! @param[out]    angles          geometry structure
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
SUBROUTINE rttov_setgeometry( &
              plane_parallel, &
              profiles,       &
              coef,           &
              angles)

  USE rttov_types, ONLY :   &
         rttov_profile,     &
         rttov_coef,        &
         rttov_geometry
  USE parkind1, ONLY : jplm
!INTF_OFF
  USE rttov_const, ONLY : deg2rad
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE

  LOGICAL(jplm),           INTENT(IN)  :: plane_parallel
  TYPE(rttov_profile),     INTENT(IN)  :: profiles(:)
  TYPE(rttov_coef),        INTENT(IN)  :: coef
  TYPE(rttov_geometry),    INTENT(OUT) :: angles(SIZE(profiles))
!INTF_END

  INTEGER(KIND=jpim) :: i, nprofiles
  REAL   (KIND=jprb) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------

!Notes on notation:
! zen  => zenith angle
!   (definition: angle at surface between view path to satellite and zenith)
! view => view angle
!   (definition: angle at the satellite between view path and nadir)
! _sq = square of given value
! _sqrt = square root of given value
! _minus1 = given value - 1
! trigonometric function abbreviations have their usual meanings
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETGEOMETRY', 0_jpim, ZHOOK_HANDLE)
  nprofiles = SIZE(profiles)
  DO i = 1, nprofiles
    angles(i)%sinzen     = SIN(profiles(i)%zenangle * deg2rad)
    angles(i)%sinzen_sq  = angles(i)%sinzen * angles(i)%sinzen
    angles(i)%coszen     = COS(profiles(i)%zenangle * deg2rad)
    angles(i)%coszen_sq  = angles(i)%coszen * angles(i)%coszen
    angles(i)%seczen     = 1.0_jprb / ABS(angles(i)%coszen)
    angles(i)%seczen_sq  = angles(i)%seczen * angles(i)%seczen
    IF (plane_parallel) THEN
      angles(i)%sinview  = angles(i)%sinzen
    ELSE
      angles(i)%sinview  = angles(i)%sinzen / coef%ratoe
    ENDIF
    angles(i)%sinview_sq = angles(i)%sinview * angles(i)%sinview
    angles(i)%cosview_sq = 1.0_jprb - angles(i)%sinview_sq
    angles(i)%normzen    = profiles(i)%zenangle / 60.0_jprb   ! normalized zenith angle for ISEM

    angles(i)%coszen_sun = COS(profiles(i)%sunzenangle * deg2rad)
    angles(i)%sinzen_sun = SIN(profiles(i)%sunzenangle * deg2rad)

    angles(i)%sinlat     = SIN(profiles(i)%latitude * deg2rad)
    angles(i)%coslat     = COS(profiles(i)%latitude * deg2rad)
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETGEOMETRY', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setgeometry
