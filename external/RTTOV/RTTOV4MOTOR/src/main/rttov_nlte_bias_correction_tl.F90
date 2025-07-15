! Description:
!> @file
!!   TL of non-LTE bias correction.
!
!> @brief
!!   TL of non-LTE bias correction.
!!
!! @param[in]     opts            options to configure the simulations
!! @param[in]     coef            RTTOV optical depth coefficient structure
!! @param[in]     profiles        input atmospheric profiles and surface variables
!! @param[in]     profiles_tl     input atmospheric profile and surface variable perturbations
!! @param[in]     chanprof        specifies channels and profiles to simulate
!! @param[in,out] rad_tl          radiance perturbations
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
SUBROUTINE rttov_nlte_bias_correction_tl(opts, coef, profiles, profiles_tl, chanprof, rad_tl)

  USE rttov_types, ONLY : rttov_options, rttov_coef, rttov_profile, rttov_chanprof, rttov_radiance
!INTF_OFF
  USE rttov_const, ONLY : deg2rad
  USE rttov_types, ONLY : rttov_nlte_coef
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_options),  INTENT(IN)    :: opts
  TYPE(rttov_coef),     INTENT(IN)    :: coef
  TYPE(rttov_profile),  INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),  INTENT(IN)    :: profiles_tl(:)
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_radiance), INTENT(INOUT) :: rad_tl
!INTF_END

  TYPE(rttov_nlte_coef), POINTER :: nlte_coef  ! nlte_coef structure

  REAL(jprb)    :: p1, p2, p3, p4, p
  REAL(jprb)    :: t_tl(3,SIZE(profiles))
  REAL(jprb)    :: x_tl(coef%nlte_coef%ncoef,SIZE(profiles)), path, znlte_tl
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

!================================!
! Calc t1 and t2 avg per profile !
!================================!
  DO iprof = 1, SIZE(profiles)
    IF (do_nlte(iprof)) THEN
      t_tl(1,iprof) = SUM(profiles_tl(iprof)%t(lb1:ub1)) / (ub1 - lb1 + 1)
      t_tl(2,iprof) = SUM(profiles_tl(iprof)%t(lb2:ub2)) / (ub2 - lb2 + 1)
      t_tl(3,iprof) = 0._jprb
      path          = 1._jprb / COS(profiles(iprof)%zenangle * deg2rad)
      x_tl(1,iprof) = 0._jprb
      x_tl(2,iprof) = 0._jprb
      x_tl(3,iprof) = 0._jprb
      x_tl(4,iprof) = 0._jprb
      x_tl(5,iprof) = COS(profiles(iprof)%sunzenangle * deg2rad) * t_tl(1,iprof)
      x_tl(6,iprof) = 0._jprb
      x_tl(7,iprof) = COS(profiles(iprof)%sunzenangle * deg2rad) * t_tl(2,iprof)
      x_tl(8,iprof) = path * t_tl(1,iprof)
      x_tl(9,iprof) = path * t_tl(2,iprof)
    ENDIF
  ENDDO

!=====================================================================!
!      Calculate radiance bias correction for required channels       !
!=====================================================================!
  DO i = 1, SIZE(chanprof)
    prof = chanprof(i)%prof

    IF (do_nlte(prof)) THEN
      chan = chanprof(i)%chan

      IF (chan >= nlte_coef%start_chan .AND. chan <= nlte_coef%end_chan) THEN

        ! corresponding array element for chan
        k = chan - nlte_coef%start_chan + 1

        znlte_tl = 0._jprb
        DO npred = 1, coef%nlte_coef%ncoef
          znlte_tl = znlte_tl + nlte_coef%coef(npred,1,1,k) * x_tl(npred,prof)
        ENDDO

        rad_tl%clear(i)  = rad_tl%clear(i)  + znlte_tl
        rad_tl%total(i)  = rad_tl%total(i)  + znlte_tl
        IF (.NOT. opts%rt_ir%pc%addpc) rad_tl%cloudy(i) = rad_tl%cloudy(i) + znlte_tl
      ENDIF
    ENDIF
  ENDDO

END SUBROUTINE rttov_nlte_bias_correction_tl
