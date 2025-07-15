! Description:
!> @file
!!   Allocate/deallocate internal sunglint structure.
!
!> @brief
!!   Allocate/deallocate internal sunglint structure.
!!
!! @details
!!   The sunglint structure contains information related to the
!!   internal sea surface BRDF model. Not all members of the
!!   need to be allocated for the TL/AD/K structures hence the
!!   optional direct argument.
!!
!! @param[out]    err            status on exit
!! @param[in,out] sunglint       sunglint structure
!! @param[in]     opts           options to configure the simulations
!! @param[in]     nprofiles      number of profiles
!! @param[in]     nomega         size of the wave spectrum array
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     init           set .TRUE. to initialise newly allocated structures, optional
!! @param[in]     direct         flag indicating whether this is the direct model structure, optional
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
SUBROUTINE rttov_alloc_sunglint( &
              err,       &
              sunglint,  &
              opts,      &
              nprofiles, &
              nomega,    &
              asw,       &
              init,      &
              direct)
!INTF_OFF
#include "throw.h"
  USE rttov_const, ONLY : nk_elf
!INTF_ON
  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_sunglint, rttov_options
  IMPLICIT NONE
  INTEGER(KIND=jpim),   INTENT(INOUT)        :: err
  TYPE(rttov_sunglint), INTENT(INOUT)        :: sunglint
  TYPE(rttov_options),  INTENT(IN)           :: opts
  INTEGER(KIND=jpim),   INTENT(IN)           :: nprofiles
  INTEGER(KIND=jpim),   INTENT(IN)           :: nomega
  INTEGER(KIND=jpim),   INTENT(IN)           :: asw
  LOGICAL(KIND=jplm),   INTENT(IN), OPTIONAL :: init
  LOGICAL(KIND=jplm),   INTENT(IN), OPTIONAL :: direct
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_sunglint.interface"

  LOGICAL(KIND=jplm) :: init1, direct1

  TRY

  init1 = .FALSE.
  IF (PRESENT(init)) init1 = init
  direct1 = .FALSE.
  IF (PRESENT(direct)) direct1 = direct
  IF (asw == 1) THEN
    CALL nullify_struct()

    ALLOCATE(sunglint%s(nprofiles), STAT = err)
    THROWM(err.NE.0,"Allocation of sunglint failed")
    IF (direct1) THEN
      IF (opts%rt_ir%solar_sea_brdf_model == 1) THEN
        ALLOCATE(sunglint%beta(nomega, nprofiles), &
                 sunglint%psi(nomega, nprofiles), STAT = err)
        THROWM(err.NE.0,"Allocation of sunglint failed")
      ELSE
        ALLOCATE(sunglint%c(nk_elf, nprofiles),         &
                 sunglint%lpm(nk_elf, nprofiles),       &
                 sunglint%gamma_exp(nk_elf, nprofiles), &
                 sunglint%jp(nk_elf, nprofiles),        &
                 sunglint%fpexp(nk_elf, nprofiles),     &
                 sunglint%fm(nk_elf, nprofiles),        &
                 sunglint%dk(nk_elf, nprofiles),        &
                 sunglint%sk2(0:nk_elf, nprofiles), STAT = err)
        THROWM(err.NE.0,"Allocation of sunglint failed")
      ENDIF
    ENDIF
    IF (init1) CALL rttov_init_sunglint(sunglint)
  ENDIF
  IF (asw == 0) THEN
    IF (ASSOCIATED(sunglint%s)) DEALLOCATE(sunglint%s, STAT = err)
    THROWM(err.NE.0,"Deallocation of sunglint%s failed")
    IF (ASSOCIATED(sunglint%beta)) DEALLOCATE(sunglint%beta, STAT = err)
    THROWM(err.NE.0,"Deallocation of sunglint%beta failed")
    IF (ASSOCIATED(sunglint%psi)) DEALLOCATE(sunglint%psi, STAT = err)
    THROWM(err.NE.0,"Deallocation of sunglint%psi failed")
    IF (ASSOCIATED(sunglint%c)) DEALLOCATE(sunglint%c, STAT = err)
    THROWM(err.NE.0,"Deallocation of sunglint%c failed")
    IF (ASSOCIATED(sunglint%lpm)) DEALLOCATE(sunglint%lpm, STAT = err)
    THROWM(err.NE.0,"Deallocation of sunglint%lpm failed")
    IF (ASSOCIATED(sunglint%gamma_exp)) DEALLOCATE(sunglint%gamma_exp, STAT = err)
    THROWM(err.NE.0,"Deallocation of sunglint%gamma_exp failed")
    IF (ASSOCIATED(sunglint%jp)) DEALLOCATE(sunglint%jp, STAT = err)
    THROWM(err.NE.0,"Deallocation of sunglint%jp failed")
    IF (ASSOCIATED(sunglint%fpexp)) DEALLOCATE(sunglint%fpexp, STAT = err)
    THROWM(err.NE.0,"Deallocation of sunglint%fpexp failed")
    IF (ASSOCIATED(sunglint%fm)) DEALLOCATE(sunglint%fm, STAT = err)
    THROWM(err.NE.0,"Deallocation of sunglint%fm failed")
    IF (ASSOCIATED(sunglint%dk)) DEALLOCATE(sunglint%dk, STAT = err)
    THROWM(err.NE.0,"Deallocation of sunglint%dk failed")
    IF (ASSOCIATED(sunglint%sk2)) DEALLOCATE(sunglint%sk2, STAT = err)
    THROWM(err.NE.0,"Deallocation of sunglint%sk2 failed")

    CALL nullify_struct()
  ENDIF

  CATCH

CONTAINS
  SUBROUTINE nullify_struct()
    NULLIFY(sunglint%s,         &
            sunglint%beta,      &
            sunglint%psi,       &
            sunglint%c,         &
            sunglint%lpm,       &
            sunglint%gamma_exp, &
            sunglint%jp,        &
            sunglint%fpexp,     &
            sunglint%fm,        &
            sunglint%dk,        &
            sunglint%sk2)
  END SUBROUTINE nullify_struct

END SUBROUTINE 
