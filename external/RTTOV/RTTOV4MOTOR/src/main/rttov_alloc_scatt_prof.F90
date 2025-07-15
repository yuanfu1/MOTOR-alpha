! Description:
!> @file
!!   Allocate/deallocate RTTOV-SCATT cloud profiles structure.
!
!> @brief
!!   Allocate/deallocate RTTOV-SCATT cloud profiles structure.
!!
!! @details
!!   The cloud profiles structure contains the cloud liquid and
!!   ice water and hydrometeor profiles for the RTTOV-SCATT direct
!!   model, the input cloud profile perturbations for the TL model,
!!   and the output gradients and Jacobians for the AD and K models.
!!   This subroutine allocates all array members of this structure
!!   for each cloud profile.
!!
!!   The profiles argument should be declared as an array of size
!!   nprof. For the K model cld_profiles_k should have size nchanprof
!!   (i.e. the total number of channels being simulated).
!!
!!
!! @param[out]    err              status on exit
!! @param[in]     nprof            number of profiles being simulated
!! @param[in,out] cld_profiles     input cloud profiles
!! @param[in]     nlev             number of levels in cloud profiles (must match the RTTOV profiles structure)
!! @param[in]     nhydro           number of hydrometeors in cloud profiles
!! @param[in]     nhydro_frac      number of hydrometeor fractions in cloud profiles (should be 1 or nhydro)
!! @param[in]     asw              1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     init             set .TRUE. to initialise newly allocated structures, optional
!! @param[in]     flux_conversion  input units: 0 (default) => kg/kg, 1,2 => kg/m2/s, optional for rain, snow
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
SUBROUTINE rttov_alloc_scatt_prof(err, nprof, cld_profiles, nlev, nhydro, nhydro_frac, asw, init, flux_conversion)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY: jpim, jplm
  USE rttov_types, ONLY : rttov_profile_cloud
!INTF_OFF
  USE parkind1, ONLY: jprb
  USE YOMHOOK,  ONLY: LHOOK , DR_HOOK
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),             INTENT(OUT)   :: err
  INTEGER(jpim),             INTENT(IN)    :: nprof
  TYPE(rttov_profile_cloud), INTENT(INOUT) :: cld_profiles (nprof)
  INTEGER(jpim),             INTENT(IN)    :: nlev
  INTEGER(jpim),             INTENT(IN)    :: nhydro
  INTEGER(jpim),             INTENT(IN)    :: nhydro_frac
  INTEGER(jpim),             INTENT(IN)    :: asw
  LOGICAL(jplm), OPTIONAL,   INTENT(IN)    :: init
  INTEGER(jpim), OPTIONAL,   INTENT(IN)    :: flux_conversion(nhydro)
!INTF_END

#include "rttov_init_scatt_prof.interface"
#include "rttov_errorreport.interface"

  INTEGER(jpim) :: iprof, flux_local(nhydro)
  LOGICAL(jplm) :: init1
  REAL(jprb)    :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------

  TRY
  IF (lhook) CALL dr_hook('RTTOV_ALLOC_SCATT_PROF',0_jpim,zhook_handle)

  init1 = .FALSE.
  IF (PRESENT(init)) init1 = init

  flux_local = 0_jpim
  IF (PRESENT(flux_conversion)) flux_local = flux_conversion

  IF (asw == 1) THEN
    DO iprof = 1, nprof

      NULLIFY( cld_profiles(iprof)%flux_conversion )
      ALLOCATE( cld_profiles(iprof)%flux_conversion (nhydro), STAT = err )
      THROWM(err.NE.0, "error allocating cld_profiles")

      cld_profiles(iprof)%nlevels      = nlev
      cld_profiles(iprof)%nhydro       = nhydro
      cld_profiles(iprof)%nhydro_frac  = nhydro_frac
      cld_profiles(iprof)%flux_conversion = flux_local
      cld_profiles(iprof)%cfrac        = 0.0_jprb

      NULLIFY( cld_profiles(iprof)%ph )
      NULLIFY( cld_profiles(iprof)%hydro )
      NULLIFY( cld_profiles(iprof)%hydro_frac )

      ALLOCATE( cld_profiles(iprof)%ph (nlev+1), STAT = err )
      THROWM(err.NE.0, "error allocating cld_profiles")
      ALLOCATE( cld_profiles(iprof)%hydro (nlev,nhydro), STAT = err )
      THROWM(err.NE.0, "error allocating cld_profiles")
      ALLOCATE( cld_profiles(iprof)%hydro_frac (nlev,nhydro_frac), STAT = err )
      THROWM(err.NE.0, "error allocating cld_profiles")

    ENDDO
    IF (init1) CALL rttov_init_scatt_prof(cld_profiles)
  ELSE
    DO iprof = 1, nprof

      DEALLOCATE( cld_profiles(iprof)%flux_conversion, STAT = err )
      THROWM(err.NE.0, "error deallocating cld_profiles")
      DEALLOCATE( cld_profiles(iprof)%ph, STAT = err )
      THROWM(err.NE.0, "error deallocating cld_profiles")
      DEALLOCATE( cld_profiles(iprof)%hydro, STAT = err )
      THROWM(err.NE.0, "error deallocating cld_profiles")
      DEALLOCATE( cld_profiles(iprof)%hydro_frac, STAT = err )
      THROWM(err.NE.0, "error deallocating cld_profiles")

      NULLIFY( cld_profiles(iprof)%flux_conversion )
      NULLIFY( cld_profiles(iprof)%ph )
      NULLIFY( cld_profiles(iprof)%hydro )
      NULLIFY( cld_profiles(iprof)%hydro_frac )

    ENDDO
  ENDIF

  IF (lhook) CALL dr_hook('RTTOV_ALLOC_SCATT_PROF',1_jpim,zhook_handle)
  CATCH
  IF (lhook) CALL dr_hook('RTTOV_ALLOC_SCATT_PROF',1_jpim,zhook_handle)
END SUBROUTINE rttov_alloc_scatt_prof

