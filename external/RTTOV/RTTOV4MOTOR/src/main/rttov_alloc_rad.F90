! Description:
!> @file
!!   Allocate/deallocate radiance structure.
!
!> @brief
!!   Allocate/deallocate radiance structure.
!!
!! @details
!!   The radiance structure contains the output radiances, BTs and
!!   reflectances for the direct model, the output radiance
!!   perturbations for the TL model and the input gradients and
!!   perturbations for the AD and K models.
!!
!!   The radiance structure also contains a geometric_height array member which
!!   contains the altitude of each level calculated within the RTTOV direct
!!   model. For consistency with other RTTOV outputs these heights are written
!!   out per-channel, but the values are the same for all channels associated
!!   with a given input profile.
!!
!! @param[out]    err            status on exit
!! @param[in]     nchanprof      total number of channels being simulated (SIZE(chanprof))
!! @param[in,out] radiance       radiances and corresponding BTs and BRFs
!! @param[in]     nlevels        number of nlevels in input profiles
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
!! @param[in,out] radiance2      secondary radiances calculated by direct model only, optional
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
SUBROUTINE rttov_alloc_rad( &
              err,       &
              nchanprof, &
              radiance,  &
              nlevels,   &
              asw,       &
              radiance2, &
              init)
!INTF_OFF
#include "throw.h"
!INTF_ON

  USE rttov_types, ONLY : rttov_radiance, rttov_radiance2
  USE parkind1, ONLY : jpim, jplm

  IMPLICIT NONE

  INTEGER(KIND=jpim),    INTENT(OUT)             :: err
  INTEGER(KIND=jpim),    INTENT(IN)              :: nchanprof
  TYPE(rttov_radiance),  INTENT(INOUT)           :: radiance
  INTEGER(KIND=jpim),    INTENT(IN)              :: nlevels
  INTEGER(KIND=jpim),    INTENT(IN)              :: asw
  TYPE(rttov_radiance2), INTENT(INOUT), OPTIONAL :: radiance2
  LOGICAL(KIND=jplm),    INTENT(IN),    OPTIONAL :: init
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_init_rad.interface"

  INTEGER(KIND=jpim) :: nlayers
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------

  TRY
  nlayers = nlevels - 1
  init1   = .FALSE.
  IF (PRESENT(init)) init1 = init

  IF (asw == 1) THEN
    ALLOCATE(radiance%clear(nchanprof),             &
             radiance%total(nchanprof),             &
             radiance%bt_clear(nchanprof),          &
             radiance%bt(nchanprof),                &
             radiance%refl_clear(nchanprof),        &
             radiance%refl(nchanprof),              &
             radiance%cloudy(nchanprof),            &
             radiance%overcast(nlayers, nchanprof), &
             radiance%quality(nchanprof),           &
             radiance%geometric_height(nlevels, nchanprof), STAT = err)
    THROWM(err.NE.0, "allocation of radiance")

    IF (PRESENT(radiance2)) THEN
      ALLOCATE(radiance2%upclear(nchanprof),       &
               radiance2%dnclear(nchanprof),       &
               radiance2%refldnclear(nchanprof),   &
               radiance2%up(nlayers, nchanprof),   &
               radiance2%down(nlayers, nchanprof), &
               radiance2%surf(nlayers, nchanprof), STAT = err)
      THROWM(err.NE.0, "allocation of radiance2")
    ENDIF

    IF (init1) CALL rttov_init_rad(radiance, radiance2)
  ENDIF

  IF (asw == 0) THEN
    DEALLOCATE(radiance%clear,      &
               radiance%total,      &
               radiance%bt_clear,   &
               radiance%bt,         &
               radiance%refl_clear, &
               radiance%refl,       &
               radiance%cloudy,     &
               radiance%overcast,   &
               radiance%quality,    &
               radiance%geometric_height, STAT = err)
    THROWM(err.NE.0, "deallocation of radiance")

    IF (PRESENT(radiance2)) THEN
      DEALLOCATE(radiance2%upclear,     &
                 radiance2%dnclear,     &
                 radiance2%refldnclear, &
                 radiance2%up,          &
                 radiance2%down,        &
                 radiance2%surf, STAT = err)
      THROWM(err.NE.0, "deallocation of radiance2")
    ENDIF

    NULLIFY(radiance%clear,      &
            radiance%total,      &
            radiance%bt_clear,   &
            radiance%bt,         &
            radiance%refl_clear, &
            radiance%refl,       &
            radiance%cloudy,     &
            radiance%overcast,   &
            radiance%quality,    &
            radiance%geometric_height)

    IF (PRESENT(radiance2)) THEN
      NULLIFY(radiance2%upclear,     &
              radiance2%dnclear,     &
              radiance2%refldnclear, &
              radiance2%up,          &
              radiance2%down,        &
              radiance2%surf)
    ENDIF
  ENDIF
  CATCH
END SUBROUTINE rttov_alloc_rad
