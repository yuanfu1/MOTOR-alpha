! Description:
!> @file
!!   Allocate/deallocate profiles structure.
!
!> @brief
!!   Allocate/deallocate profiles structure.
!!
!! @details
!!   The profiles structure contains the input profiles for the direct
!!   model, the input profile perturbations for the TL model, and the
!!   output gradients and Jacobians for the AD and K models. This
!!   subroutine allocates all array members of this structure for
!!   each profile.
!!
!!   The profiles argument should be declared as an array of size
!!   nprofiles. For the K model profiles_k should have size nchanprof
!!   (i.e. the total number of channels being simulated).
!!
!!   For PC-RTTOV the profiles_k_pc array should be of size
!!   (npcscores * nprofiles) and the profiles_k_rec array should be
!!   of size (nchannels_rec * nprofiles).
!!
!!   The default value of gas_units is set here to kg/kg over moist air.
!!   The default value of mmr_cldaer is set here to TRUE (kg/kg).
!!
!! @param[out]    err            status on exit
!! @param[in]     nprofiles      number of profiles being simulated
!! @param[in,out] profiles       input atmospheric profiles and surface variables
!! @param[in]     nlevels        number of levels in input profiles
!! @param[in]     opts           options to configure the simulations
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     coefs          coefficients structure for instrument to simulate
!! @param[in]     init           set .TRUE. to initialise newly allocated structures, optional
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
SUBROUTINE rttov_alloc_prof( &
              err,       &
              nprofiles, &
              profiles,  &
              nlevels,   &
              opts,      &
              asw,       &
              coefs,     &
              init)
!INTF_OFF
#include "throw.h"
  USE rttov_const, ONLY : gas_unit_specconc, ncldtyp
!INTF_ON

  USE rttov_types, ONLY : rttov_options, rttov_profile, rttov_coefs
  USE parkind1, ONLY : jpim, jplm

  IMPLICIT NONE

  INTEGER(KIND=jpim),  INTENT(OUT)            :: err
  INTEGER(KIND=jpim),  INTENT(IN)             :: nprofiles
  TYPE(rttov_profile), INTENT(INOUT)          :: profiles(nprofiles)
  INTEGER(KIND=jpim),  INTENT(IN)             :: nlevels
  TYPE(rttov_options), INTENT(IN)             :: opts
  INTEGER(KIND=jpim),  INTENT(IN)             :: asw            ! 1=allocate, 0=deallocate
  TYPE(rttov_coefs),   INTENT(IN),   OPTIONAL :: coefs
  LOGICAL(KIND=jplm),  INTENT(IN),   OPTIONAL :: init
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_init_prof.interface"

  INTEGER(KIND=jpim) :: j
  INTEGER(KIND=jpim) :: nlayers
  INTEGER(KIND=jpim) :: naertyp
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------

  TRY
  nlayers = nlevels - 1
  init1   = .FALSE.
  IF (PRESENT(init)) init1 = init

  IF (asw == 1) THEN
    IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
      IF (.NOT. PRESENT(coefs)) THEN
        err = errorstatus_fatal
        THROWM(err .NE. 0, "Dummy argument coefs is required")
      ENDIF
      naertyp = MAX(1_jpim, coefs%coef_scatt%optp_aer%ntypes)
    ENDIF

    DO j = 1, nprofiles
      NULLIFY (profiles(j)%p)
      NULLIFY (profiles(j)%t)
      NULLIFY (profiles(j)%q)
      NULLIFY (profiles(j)%o3)
      NULLIFY (profiles(j)%co2)
      NULLIFY (profiles(j)%n2o)
      NULLIFY (profiles(j)%co)
      NULLIFY (profiles(j)%ch4)
      NULLIFY (profiles(j)%so2)
      NULLIFY (profiles(j)%clw)
      NULLIFY (profiles(j)%aerosols)
      NULLIFY (profiles(j)%cloud)
      NULLIFY (profiles(j)%cfrac)
      NULLIFY (profiles(j)%clwde)
      NULLIFY (profiles(j)%icede)
    ENDDO

    DO j = 1, nprofiles
      profiles(j)%gas_units = gas_unit_specconc ! kg/kg over moist air
      profiles(j)%mmr_cldaer = .TRUE.
      profiles(j)%nlevels = nlevels
      profiles(j)%nlayers = nlevels - 1
      ALLOCATE (profiles(j)%p(nlevels), STAT = err)
      THROWM(err .NE. 0, "allocation of profiles%p")
      ALLOCATE (profiles(j)%t(nlevels), STAT = err)
      THROWM(err .NE. 0, "allocation of profiles%t")
      ALLOCATE (profiles(j)%q(nlevels), STAT = err)
      THROWM(err .NE. 0, "allocation of profiles%q")
      IF (opts%rt_all%ozone_data) THEN
        ALLOCATE (profiles(j)%o3(nlevels), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%o3 ")
      ENDIF
      IF (opts%rt_all%co2_data) THEN
        ALLOCATE (profiles(j)%co2(nlevels), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%co2")
      ENDIF
      IF (opts%rt_all%n2o_data) THEN
        ALLOCATE (profiles(j)%n2o(nlevels), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%n2o")
      ENDIF
      IF (opts%rt_all%co_data) THEN
        ALLOCATE (profiles(j)%co(nlevels), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%co")
      ENDIF
      IF (opts%rt_all%ch4_data) THEN
        ALLOCATE (profiles(j)%ch4(nlevels), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%ch4")
      ENDIF
      IF (opts%rt_all%so2_data) THEN
        ALLOCATE (profiles(j)%so2(nlevels), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%so2")
      ENDIF
      IF (opts%rt_mw%clw_data) THEN
        ALLOCATE (profiles(j)%clw(nlevels), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%clw")
      ENDIF

      IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
        ALLOCATE (profiles(j)%aerosols(naertyp, nlayers), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%aerosols")
      ENDIF
      IF (opts%rt_ir%addclouds) THEN
        ALLOCATE (profiles(j)%cfrac(nlayers), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%cfrac")
        IF (.NOT. opts%rt_ir%user_cld_opt_param) THEN
          ALLOCATE (profiles(j)%cloud(ncldtyp, nlayers), STAT = err)
          THROWM(err .NE. 0, "allocation of profiles%cloud")
          ALLOCATE (profiles(j)%clwde(nlayers), STAT = err)
          THROWM(err .NE. 0, "allocation of profiles%clwde")
          ALLOCATE (profiles(j)%icede(nlayers), STAT = err)
          THROWM(err .NE. 0, "allocation of profiles%icede")
        ENDIF
      ENDIF
    ENDDO
    IF (init1) CALL rttov_init_prof(profiles)
  ENDIF ! asw == 1

  IF (asw == 0) THEN
    DO j = 1, nprofiles
      DEALLOCATE (profiles(j)%p, STAT = err)
      THROWM(err .NE. 0, "deallocation of profiles%p")
      DEALLOCATE (profiles(j)%t, STAT = err)
      THROWM(err .NE. 0, "deallocation of profiles%t")
      DEALLOCATE (profiles(j)%q, STAT = err)
      THROWM(err .NE. 0, "deallocation of profiles%q")
      IF (ASSOCIATED(profiles(j)%o3)) THEN
        DEALLOCATE (profiles(j)%o3, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%o3")
      ENDIF
      IF (ASSOCIATED(profiles(j)%co2)) THEN
        DEALLOCATE (profiles(j)%co2, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%co2")
      ENDIF
      IF (ASSOCIATED(profiles(j)%n2o)) THEN
        DEALLOCATE (profiles(j)%n2o, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%n2o")
      ENDIF
      IF (ASSOCIATED(profiles(j)%co)) THEN
        DEALLOCATE (profiles(j)%co, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%co")
      ENDIF
      IF (ASSOCIATED(profiles(j)%ch4)) THEN
        DEALLOCATE (profiles(j)%ch4, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%ch4")
      ENDIF
      IF (ASSOCIATED(profiles(j)%so2)) THEN
        DEALLOCATE (profiles(j)%so2, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%so2")
      ENDIF
      IF (ASSOCIATED(profiles(j)%clw)) THEN
        DEALLOCATE (profiles(j)%clw, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%clw")
      ENDIF
      IF (ASSOCIATED(profiles(j)%aerosols)) THEN
        DEALLOCATE (profiles(j)%aerosols, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%aerosols")
      ENDIF
      IF (ASSOCIATED(profiles(j)%cloud)) THEN
        DEALLOCATE (profiles(j)%cloud, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%cloud")
      ENDIF
      IF (ASSOCIATED(profiles(j)%cfrac)) THEN
        DEALLOCATE (profiles(j)%cfrac, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%cfrac")
      ENDIF
      IF (ASSOCIATED(profiles(j)%clwde)) THEN
        DEALLOCATE (profiles(j)%clwde, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%clwde")
      ENDIF
      IF (ASSOCIATED(profiles(j)%icede)) THEN
        DEALLOCATE (profiles(j)%icede, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%icede")
      ENDIF
      NULLIFY (profiles(j)%p)
      NULLIFY (profiles(j)%t)
      NULLIFY (profiles(j)%q)
      NULLIFY (profiles(j)%o3)
      NULLIFY (profiles(j)%co2)
      NULLIFY (profiles(j)%n2o)
      NULLIFY (profiles(j)%co)
      NULLIFY (profiles(j)%ch4)
      NULLIFY (profiles(j)%so2)
      NULLIFY (profiles(j)%clw)
      NULLIFY (profiles(j)%aerosols)
      NULLIFY (profiles(j)%cloud)
      NULLIFY (profiles(j)%cfrac)
      NULLIFY (profiles(j)%clwde)
      NULLIFY (profiles(j)%icede)
    ENDDO
  ENDIF ! asw == 0
  CATCH
END SUBROUTINE rttov_alloc_prof
