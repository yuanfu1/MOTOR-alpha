! Description:
!> @file
!!   K of non-LTE bias correction.
!
!> @brief
!!   K of non-LTE bias correction.
!!
!! @param[in]     coef            RTTOV optical depth coefficient structure
!! @param[in]     profiles        input atmospheric profiles and surface variables
!! @param[in,out] profiles_k      output gradient wrt atmospheric profile and surface variables
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in]     rad_k           radiance increments
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
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_nlte_bias_correction_k(coef, profiles, profiles_k, chanprof, rad_k)

  USE rttov_types, ONLY : rttov_coef, rttov_profile, rttov_chanprof, rttov_radiance
!INTF_OFF
  USE rttov_const, ONLY : deg2rad
  USE rttov_types, ONLY : rttov_nlte_coef
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_coef),     INTENT(IN)    :: coef
  TYPE(rttov_profile),  INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),  INTENT(INOUT) :: profiles_k(:)
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_radiance), INTENT(IN)    :: rad_k
!INTF_END

  TYPE(rttov_nlte_coef), POINTER :: nlte_coef  ! nlte_coef structure

  REAL(jprb)    :: p1, p2, p3, p4, p
  REAL(jprb)    :: t_k(3,SIZE(chanprof))
  REAL(jprb)    :: x_k(coef%nlte_coef%ncoef,SIZE(chanprof)), path
  INTEGER(jpim) :: bounds(2,2,SIZE(profiles))
  INTEGER(jpim) :: i, j, k
  INTEGER(jpim) :: iprof, prof, chan, npred
  LOGICAL(jplm) :: do_nlte(SIZE(profiles))

! Initialise variables and define macros
!========================================

#define lb1 bounds(1,1,iprof)
#define ub1 bounds(2,1,iprof)
#define lb2 bounds(1,2,iprof)
#define ub2 bounds(2,2,iprof)

  nlte_coef => coef%nlte_coef
  do_nlte(:) = .TRUE.

  p1 = 0.005_jprb; p2 = 0.2244_jprb; p3 = 0.3454_jprb ; p4 = 51.5278_jprb
  bounds = -99_jpim

!===============================================================================!
! Static code (same for FWD/TL/AD/K) - calc interpolation coefs per profile (z) !
!===============================================================================!
  DO iprof = 1, SIZE(profiles)

    IF (profiles(iprof)%sunzenangle > (90._jprb - 1.E-6_jprb) .OR. profiles(iprof)%sunzenangle < 0._jprb) THEN
      do_nlte(iprof) = .FALSE.
      CYCLE
    ENDIF

    ! find pressure bounds for finding average temperature
    p = 0.; j = 1
    DO WHILE (p < p3)
      p = profiles(iprof)%p(j)

      IF ( p < p1 ) THEN
        j = j + 1
        CYCLE
      ELSE IF (p >= p1 .AND. p < p2) THEN ! contribution to t1_avg
        IF (lb1 < 0) lb1 = j
        ub1 = j
      ELSE IF (p >= p3 .AND. p <= p4) THEN! contribution to t2_avg
        IF (lb2 < 0) lb2 = j
        ub2 = j
      ENDIF
      j = j + 1
    ENDDO
  ENDDO
! END static code


!=====================================================================!
!      Calculate radiance bias correction for required channels       !
!=====================================================================!
  x_k = 0._jprb
  DO i = 1, SIZE(chanprof)
    prof = chanprof(i)%prof

    IF (do_nlte(prof)) THEN
      chan = chanprof(i)%chan

      IF (chan >= nlte_coef%start_chan .AND. chan <= nlte_coef%end_chan) THEN

        ! corresponding array element for chan
        k = chan - nlte_coef%start_chan + 1

        ! Only consider the total radiance increment as the input clear and
        ! cloudy increments should be zero
        DO npred = 1, coef%nlte_coef%ncoef
          x_k(npred,i) = x_k(npred,i) + nlte_coef%coef(npred,1,1,k) * rad_k%total(i)
        ENDDO
      ENDIF
    ENDIF
  ENDDO

!================================!
! Calc t1 and t2 avg per profile !
!================================!
  t_k = 0._jprb
  DO i = 1, SIZE(chanprof)
    iprof = chanprof(i)%prof
    IF (do_nlte(iprof)) THEN
      path     = 1._jprb / COS(profiles(iprof)%zenangle * deg2rad)
      t_k(2,i) = t_k(2,i) + path * x_k(9,i)
      t_k(1,i) = t_k(1,i) + path * x_k(8,i)
      t_k(2,i) = t_k(2,i) + COS(profiles(iprof)%sunzenangle * deg2rad) * x_k(7,i)

      t_k(1,i) = t_k(1,i) + COS(profiles(iprof)%sunzenangle * deg2rad) * x_k(5,i)

      profiles_k(i)%t(lb1:ub1) = profiles_k(i)%t(lb1:ub1) + t_k(1,i) / (ub1 - lb1 + 1)
      profiles_k(i)%t(lb2:ub2) = profiles_k(i)%t(lb2:ub2) + t_k(2,i) / (ub2 - lb2 + 1)
    ENDIF
  ENDDO

END SUBROUTINE rttov_nlte_bias_correction_k
