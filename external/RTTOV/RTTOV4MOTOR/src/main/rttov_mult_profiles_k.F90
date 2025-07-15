! Description:
!> @file
!!   Compute final PC-RTTOV Jacobians
!
!> @brief
!!   Compute final PC-RTTOV Jacobians
!!
!! @details
!!   The PC-RTTOV K model computes first the Jacobians of the predictor channel
!!   radiances wrt each profile variable (in profiles_k) and then the Jacobians
!!   of the PC scores or reconstructed radiances/BTs wrt the predictor channels
!!   (in total_k_pc). This subroutine multiplies Jacobians together to obtain
!!   the final output Jacobians of PC scores or reconstructed radiances/BTs wrt
!!   each profile variable (profiles_k_rec).
!!
!! @param[in,out] profiles_k_rec    Output Jacobians
!! @param[in]     profiles_k        RTTOV optical depth coefficient structure
!! @param[in]     total_k_pc        input atmospheric profiles and surface variables
!! @param[in]     opts              options to configure the simulations
!! @param[in]     coef_pccomp       PC-RTTOV coefficients structure
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
SUBROUTINE rttov_mult_profiles_k(profiles_k_rec, profiles_k, total_k_pc, opts, coef_pccomp)

  USE parkind1, ONLY : jprb
  USE rttov_types, ONLY : rttov_profile, rttov_options, rttov_coef_pccomp
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
  TYPE(rttov_profile),     INTENT(INOUT) :: profiles_k_rec(:)
  TYPE(rttov_profile),     INTENT(IN)    :: profiles_k(:)
  REAL(KIND=jprb),         INTENT(IN)    :: total_k_pc(:,:,:)
  TYPE(rttov_options),     INTENT(IN)    :: opts
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
!INTF_END
  REAL(kind=jprb) :: prof_k_array(SIZE(total_k_pc, 1) * SIZE(total_k_pc, 3), &
                                  SIZE(profiles_k(1)%p), 8)
  REAL(kind=jprb) :: prof_k_array_surf(SIZE(total_k_pc, 1) * SIZE(total_k_pc, 3), 10)
  REAL(kind=jprb), ALLOCATABLE :: prof_k_array_aer(:,:,:)
  REAL(kind=jprb), ALLOCATABLE :: prof_k_array_cfrac(:,:), prof_k_array_cld(:,:,:)
  REAL(kind=jprb), ALLOCATABLE :: prof_k_array_clwde(:,:), prof_k_array_icede(:,:)
  INTEGER(KIND=jpim) :: nchannels_rec_p
  INTEGER(KIND=jpim) :: nchannels_p
  INTEGER(KIND=jpim) :: nchannels
  INTEGER(KIND=jpim) :: nlevels
  INTEGER(KIND=jpim) :: ichan_in1
  INTEGER(KIND=jpim) :: ichan1
  INTEGER(KIND=jpim) :: ichan_in
  INTEGER(KIND=jpim) :: ichan
  INTEGER(KIND=jpim) :: iprof
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: ilev, lay, iae, icld

  nchannels_p     = SIZE(total_k_pc, 1)
  nchannels_rec_p = SIZE(total_k_pc, 2)
  nprofiles       = SIZE(total_k_pc, 3)
  nlevels         = SIZE(profiles_k(1)%p)
  nchannels       = nchannels_p * nprofiles

! DAR: Re-pack predictor channel Jacobians to improve data locality
  DO iprof = 1, nprofiles
    DO ichan = 1, nchannels_p
      ichan1 = ichan + (iprof - 1) * nchannels_p
      prof_k_array(ichan1, :, 1) = profiles_k(ichan1)%t
      prof_k_array(ichan1, :, 2) = profiles_k(ichan1)%q
    ENDDO
  ENDDO

  IF (ASSOCIATED(profiles_k_rec(1)%o3)) THEN
    DO iprof = 1, nprofiles
      DO ichan = 1, nchannels_p
        ichan1 = ichan + (iprof - 1) * nchannels_p
        prof_k_array(ichan1, :, 3) = profiles_k(ichan1)%o3
      ENDDO
    ENDDO
  ENDIF

  IF (ASSOCIATED(profiles_k_rec(1)%co2)) THEN
    DO iprof = 1, nprofiles
      DO ichan = 1, nchannels_p
        ichan1 = ichan + (iprof - 1) * nchannels_p
        prof_k_array(ichan1, :, 4) = profiles_k(ichan1)%co2
      ENDDO
    ENDDO
  ENDIF

  IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN

    IF (ASSOCIATED(profiles_k_rec(1)%n2o)) THEN
      DO iprof = 1, nprofiles
        DO ichan = 1, nchannels_p
          ichan1 = ichan + (iprof - 1) * nchannels_p
          prof_k_array(ichan1, :, 5) = profiles_k(ichan1)%n2o
        ENDDO
      ENDDO
    ENDIF

    IF (ASSOCIATED(profiles_k_rec(1)%co)) THEN
      DO iprof = 1, nprofiles
        DO ichan = 1, nchannels_p
          ichan1 = ichan + (iprof - 1) * nchannels_p
          prof_k_array(ichan1, :, 6) = profiles_k(ichan1)%co
        ENDDO
      ENDDO
    ENDIF

    IF (ASSOCIATED(profiles_k_rec(1)%ch4)) THEN
      DO iprof = 1, nprofiles
        DO ichan = 1, nchannels_p
          ichan1 = ichan + (iprof - 1) * nchannels_p
          prof_k_array(ichan1, :, 7) = profiles_k(ichan1)%ch4
        ENDDO
      ENDDO
    ENDIF

  ENDIF

  IF (opts%interpolation%lgradp) THEN
    DO iprof = 1, nprofiles
      DO ichan = 1, nchannels_p
        ichan1 = ichan + (iprof - 1) * nchannels_p
        prof_k_array(ichan1, :, 8) = profiles_k(ichan1)%p
      ENDDO
    ENDDO
  ENDIF

  DO iprof = 1, nprofiles
     DO ichan = 1, nchannels_p
        ichan1 = ichan + (iprof - 1) * nchannels_p
        prof_k_array_surf(ichan1, 1) = profiles_k(ichan1)%ctp
        prof_k_array_surf(ichan1, 2) = profiles_k(ichan1)%cfraction
        prof_k_array_surf(ichan1, 3) = profiles_k(ichan1)%s2m%p
        prof_k_array_surf(ichan1, 4) = profiles_k(ichan1)%s2m%t
        prof_k_array_surf(ichan1, 5) = profiles_k(ichan1)%s2m%q
        prof_k_array_surf(ichan1, 6) = profiles_k(ichan1)%s2m%o
        prof_k_array_surf(ichan1, 7) = profiles_k(ichan1)%s2m%u
        prof_k_array_surf(ichan1, 8) = profiles_k(ichan1)%s2m%v
        prof_k_array_surf(ichan1, 9) = profiles_k(ichan1)%s2m%wfetc
        prof_k_array_surf(ichan1, 10) = profiles_k(ichan1)%skin%t
     ENDDO
  ENDDO

  IF (ASSOCIATED(profiles_k_rec(1)%aerosols)) THEN
    ALLOCATE(prof_k_array_aer(nchannels, &
                              SIZE(profiles_k(1)%aerosols,dim=1), &
                              SIZE(profiles_k(1)%aerosols,dim=2)))
    DO iprof = 1, nprofiles
      DO ichan = 1, nchannels_p
        ichan1 = ichan + (iprof - 1) * nchannels_p
        prof_k_array_aer(ichan1, :, :) = profiles_k(ichan1)%aerosols
      ENDDO
    ENDDO
  ENDIF

  IF (ASSOCIATED(profiles_k_rec(1)%cloud)) THEN
    ALLOCATE(prof_k_array_cld(nchannels, &
                              SIZE(profiles_k(1)%cloud, dim=1), &
                              SIZE(profiles_k(1)%cloud, dim=2)))
    DO iprof = 1, nprofiles
      DO ichan = 1, nchannels_p
        ichan1 = ichan + (iprof - 1) * nchannels_p
        prof_k_array_cld(ichan1, :, :) = profiles_k(ichan1)%cloud
      ENDDO
    ENDDO
  ENDIF

  IF (ASSOCIATED(profiles_k_rec(1)%cfrac)) THEN
    ALLOCATE(prof_k_array_cfrac(nchannels, SIZE(profiles_k(1)%cfrac)))

    DO iprof = 1, nprofiles
      DO ichan = 1, nchannels_p
        ichan1 = ichan + (iprof - 1) * nchannels_p
        prof_k_array_cfrac(ichan1, :) = profiles_k(ichan1)%cfrac
      ENDDO
    ENDDO
  ENDIF

  IF (ASSOCIATED(profiles_k_rec(1)%clwde)) THEN
    ALLOCATE(prof_k_array_clwde(nchannels, SIZE(profiles_k(1)%clwde)))

    DO iprof = 1, nprofiles
      DO ichan = 1, nchannels_p
        ichan1 = ichan + (iprof - 1) * nchannels_p
        prof_k_array_clwde(ichan1, :) = profiles_k(ichan1)%clwde
      ENDDO
    ENDDO
  ENDIF

  IF (ASSOCIATED(profiles_k_rec(1)%icede)) THEN
    ALLOCATE(prof_k_array_icede(nchannels, SIZE(profiles_k(1)%icede)))

    DO iprof = 1, nprofiles
      DO ichan = 1, nchannels_p
        ichan1 = ichan + (iprof - 1) * nchannels_p
        prof_k_array_icede(ichan1, :) = profiles_k(ichan1)%icede
      ENDDO
    ENDDO
  ENDIF

  DO iprof = 1, nprofiles
    DO ichan_in = 1, nchannels_rec_p
      ichan_in1 = ichan_in + (iprof - 1) * nchannels_rec_p

      ! DAR: Only calculate reconstructed p Jacobian if necessary
      ! Massively improved performance by looping over levels and using
      ! contiguous data to reconstruct Jacobians
      IF (opts%interpolation%lgradp) THEN
        DO ilev = 1, nlevels
          profiles_k_rec(ichan_in1)%p(ilev) = &
            DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
              prof_k_array((iprof - 1)*nchannels_p+1 : iprof*nchannels_p,ilev,8))
        ENDDO
      ENDIF

      DO ilev = 1, nlevels
        profiles_k_rec(ichan_in1)%t(ilev) = &
          DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
            prof_k_array((iprof - 1)*nchannels_p+1 : iprof*nchannels_p,ilev,1))
      ENDDO

      DO ilev = 1, nlevels
        profiles_k_rec(ichan_in1)%q(ilev) = &
          DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
            prof_k_array((iprof - 1)*nchannels_p+1 : iprof*nchannels_p,ilev,2))
      ENDDO

      IF (ASSOCIATED(profiles_k_rec(ichan_in1)%o3)) THEN
        DO ilev = 1, nlevels
          profiles_k_rec(ichan_in1)%o3(ilev) = &
            DOT_PRODUCT(total_k_pc(1:nchannels_p,ichan_in, iprof), &
              prof_k_array((iprof - 1)*nchannels_p+1 : iprof*nchannels_p,ilev,3))
        ENDDO
      ENDIF

      IF (ASSOCIATED(profiles_k_rec(ichan_in1)%co2)) THEN
        DO ilev = 1, nlevels
          profiles_k_rec(ichan_in1)%co2(ilev) = &
            DOT_PRODUCT(total_k_pc(1:nchannels_p,ichan_in, iprof), &
              prof_k_array((iprof - 1)*nchannels_p+1 : iprof*nchannels_p,ilev,4))
        ENDDO
      ENDIF

      IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN

        IF (ASSOCIATED(profiles_k_rec(ichan_in1)%n2o)) THEN
          DO ilev = 1, nlevels
            profiles_k_rec(ichan_in1)%n2o(ilev) = &
              DOT_PRODUCT(total_k_pc(1:nchannels_p,ichan_in, iprof), &
                prof_k_array((iprof - 1)*nchannels_p+1 : iprof*nchannels_p,ilev,5))
          ENDDO
        ENDIF
        IF (ASSOCIATED(profiles_k_rec(ichan_in1)%co)) THEN
          DO ilev = 1, nlevels
            profiles_k_rec(ichan_in1)%co(ilev) = &
              DOT_PRODUCT(total_k_pc(1:nchannels_p,ichan_in, iprof), &
                prof_k_array((iprof - 1)*nchannels_p+1 : iprof*nchannels_p,ilev,6))
          ENDDO
        ENDIF
        IF (ASSOCIATED(profiles_k_rec(ichan_in1)%ch4)) THEN
          DO ilev = 1, nlevels
            profiles_k_rec(ichan_in1)%ch4(ilev) = &
              DOT_PRODUCT(total_k_pc(1:nchannels_p,ichan_in, iprof), &
                prof_k_array((iprof - 1)*nchannels_p+1 : iprof*nchannels_p,ilev,7))
          ENDDO
        ENDIF

      ELSE

        IF (ASSOCIATED(profiles_k_rec(ichan_in1)%n2o)) &
          profiles_k_rec(ichan_in1)%n2o = 0._jprb
        IF (ASSOCIATED(profiles_k_rec(ichan_in1)%co)) &
          profiles_k_rec(ichan_in1)%co = 0._jprb
        IF (ASSOCIATED(profiles_k_rec(ichan_in1)%ch4)) &
          profiles_k_rec(ichan_in1)%ch4 = 0._jprb

      ENDIF

      IF (ASSOCIATED(profiles_k_rec(ichan_in1)%so2)) &
        profiles_k_rec(ichan_in1)%so2 = 0._jprb
      IF (ASSOCIATED(profiles_k_rec(ichan_in1)%clw)) &
        profiles_k_rec(ichan_in1)%clw = 0._jprb

      profiles_k_rec(ichan_in1)%ctp         = &
        DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
        prof_k_array_surf((iprof - 1) * nchannels_p + 1 :iprof * nchannels_p,1))
      profiles_k_rec(ichan_in1)%cfraction   = &
        DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
        prof_k_array_surf((iprof - 1) * nchannels_p + 1 :iprof * nchannels_p,2))
      profiles_k_rec(ichan_in1)%s2m%p       = &
        DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
        prof_k_array_surf((iprof - 1) * nchannels_p + 1 :iprof * nchannels_p,3))
      profiles_k_rec(ichan_in1)%s2m%t       = &
        DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
        prof_k_array_surf((iprof - 1) * nchannels_p + 1 :iprof * nchannels_p,4))
      profiles_k_rec(ichan_in1)%s2m%q       = &
        DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
        prof_k_array_surf((iprof - 1) * nchannels_p + 1 :iprof * nchannels_p,5))
      profiles_k_rec(ichan_in1)%s2m%o       = &
        DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
        prof_k_array_surf((iprof - 1) * nchannels_p + 1 :iprof * nchannels_p,6))
      profiles_k_rec(ichan_in1)%s2m%u       = &
        DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
        prof_k_array_surf((iprof - 1) * nchannels_p + 1 :iprof * nchannels_p,7))
      profiles_k_rec(ichan_in1)%s2m%v       = &
        DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
        prof_k_array_surf((iprof - 1) * nchannels_p + 1 :iprof * nchannels_p,8))
      profiles_k_rec(ichan_in1)%s2m%wfetc   = &
        DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
        prof_k_array_surf((iprof - 1) * nchannels_p + 1 :iprof * nchannels_p,9))
      profiles_k_rec(ichan_in1)%skin%t      = &
        DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
        prof_k_array_surf((iprof - 1) * nchannels_p + 1 :iprof * nchannels_p,10))

      IF (ASSOCIATED(profiles_k_rec(ichan_in1)%aerosols)) THEN
        DO lay = 1, SIZE(profiles_k_rec(ichan_in1)%aerosols, dim=2)
          DO iae = 1, SIZE(profiles_k_rec(ichan_in1)%aerosols, dim=1)
            profiles_k_rec(ichan_in1)%aerosols(iae, lay) = &
              DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
              prof_k_array_aer((iprof - 1)*nchannels_p+1:iprof*nchannels_p,iae,lay))
          ENDDO
        ENDDO
      ENDIF

      IF (ASSOCIATED(profiles_k_rec(ichan_in1)%cloud)) THEN
        DO lay = 1, SIZE(profiles_k_rec(ichan_in1)%cloud, dim=2)
          DO icld = 1, SIZE(profiles_k_rec(ichan_in1)%cloud, dim=1)
            profiles_k_rec(ichan_in1)%cloud(icld, lay) = &
              DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
              prof_k_array_cld((iprof - 1)*nchannels_p+1:iprof*nchannels_p,icld,lay))
          ENDDO
        ENDDO
      ENDIF

      IF (ASSOCIATED(profiles_k_rec(ichan_in1)%cfrac)) THEN
        DO lay = 1, SIZE(profiles_k_rec(ichan_in1)%cfrac)
          profiles_k_rec(ichan_in1)%cfrac(lay) = &
            DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
            prof_k_array_cfrac((iprof-1)*nchannels_p+1:iprof*nchannels_p,lay))
        ENDDO
      ENDIF

      IF (ASSOCIATED(profiles_k_rec(ichan_in1)%clwde)) THEN
        DO lay = 1, SIZE(profiles_k_rec(ichan_in1)%clwde)
          profiles_k_rec(ichan_in1)%clwde(lay) = &
            DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
            prof_k_array_clwde((iprof-1)*nchannels_p+1:iprof*nchannels_p,lay))
        ENDDO
      ENDIF

      IF (ASSOCIATED(profiles_k_rec(ichan_in1)%icede)) THEN
        DO lay = 1, SIZE(profiles_k_rec(ichan_in1)%icede)
          profiles_k_rec(ichan_in1)%icede(lay) = &
            DOT_PRODUCT(total_k_pc(1:nchannels_p, ichan_in, iprof), &
            prof_k_array_icede((iprof-1)*nchannels_p+1:iprof*nchannels_p,lay))
        ENDDO
      ENDIF
    ENDDO
  ENDDO

  IF (ALLOCATED(prof_k_array_aer)) DEALLOCATE(prof_k_array_aer)
  IF (ALLOCATED(prof_k_array_cld)) DEALLOCATE(prof_k_array_cld)
  IF (ALLOCATED(prof_k_array_cfrac)) DEALLOCATE(prof_k_array_cfrac)
  IF (ALLOCATED(prof_k_array_clwde)) DEALLOCATE(prof_k_array_clwde)
  IF (ALLOCATED(prof_k_array_icede)) DEALLOCATE(prof_k_array_icede)
END SUBROUTINE rttov_mult_profiles_k
