! Description:
!> @file
!!   Allocate/deallocate internal profiles structure.
!
!> @brief
!!   Allocate/deallocate internal profiles structure.
!!
!! @details
!!   The internal profiles structure contains the gas and, when relevant, 
!!   cloud and aerosol profiles in the units used internally by RTTOV.
!!   Units are:
!!     gases    - ppmv over dry air
!!     clouds   - g/m3
!!     aerosols - number density cm-3
!!
!!   Only gas, cloud and aerosol components are required. This is only
!!   intended for use internally by RTTOV.
!!
!! @param[out]    err            status on exit
!! @param[in,out] profiles_int   internal profile structure
!! @param[in]     nlevels        number of levels in input profiles
!! @param[in]     opts           options to configure the simulations
!! @param[in]     coefs          coefficients structure for instrument to simulate
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
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
SUBROUTINE rttov_alloc_prof_internal( &
              err,          &
              profiles_int, &
              nlevels,      &
              opts,         &
              coefs,        &
              asw)
!INTF_OFF
#include "throw.h"
  USE rttov_const, ONLY : gas_unit_ppmvdry, ncldtyp
!INTF_ON

  USE rttov_types, ONLY : rttov_options, rttov_coefs, rttov_profile
  USE parkind1, ONLY : jpim

  IMPLICIT NONE

  INTEGER(KIND=jpim),  INTENT(OUT)            :: err
  TYPE(rttov_profile), INTENT(INOUT)          :: profiles_int(:)
  INTEGER(KIND=jpim),  INTENT(IN)             :: nlevels
  TYPE(rttov_options), INTENT(IN)             :: opts
  TYPE(rttov_coefs),   INTENT(IN)             :: coefs
  INTEGER(KIND=jpim),  INTENT(IN)             :: asw            ! 1=allocate, 0=deallocate
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(KIND=jpim) :: j
  INTEGER(KIND=jpim) :: nprofiles, nlayers
  INTEGER(KIND=jpim) :: naertyp
!- End of header --------------------------------------------------------

  TRY
  nprofiles = SIZE(profiles_int)
  nlayers = nlevels - 1

  IF (asw == 1) THEN
    IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) &
          naertyp = MAX(1_jpim, coefs%coef_scatt%optp_aer%ntypes)

    DO j = 1, nprofiles
      NULLIFY (profiles_int(j)%p)
      NULLIFY (profiles_int(j)%t)
      NULLIFY (profiles_int(j)%q)
      NULLIFY (profiles_int(j)%o3)
      NULLIFY (profiles_int(j)%co2)
      NULLIFY (profiles_int(j)%n2o)
      NULLIFY (profiles_int(j)%co)
      NULLIFY (profiles_int(j)%ch4)
      NULLIFY (profiles_int(j)%so2)
      NULLIFY (profiles_int(j)%clw)
      NULLIFY (profiles_int(j)%aerosols)
      NULLIFY (profiles_int(j)%cloud)
      NULLIFY (profiles_int(j)%cfrac)
      NULLIFY (profiles_int(j)%clwde)
      NULLIFY (profiles_int(j)%icede)
    ENDDO

    DO j = 1, nprofiles
      profiles_int(j)%gas_units = gas_unit_ppmvdry
      profiles_int(j)%mmr_cldaer = .FALSE.
      profiles_int(j)%nlevels = nlevels
      profiles_int(j)%nlayers = nlayers

      ALLOCATE (profiles_int(j)%q(nlevels), STAT = err)
      THROWM(err .NE. 0, "allocation of profiles%q")
      IF (opts%rt_all%ozone_data) THEN
        ALLOCATE (profiles_int(j)%o3(nlevels), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%o3 ")
      ENDIF
      IF (opts%rt_all%co2_data) THEN
        ALLOCATE (profiles_int(j)%co2(nlevels), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%co2")
      ENDIF
      IF (opts%rt_all%n2o_data) THEN
        ALLOCATE (profiles_int(j)%n2o(nlevels), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%n2o")
      ENDIF
      IF (opts%rt_all%co_data) THEN
        ALLOCATE (profiles_int(j)%co(nlevels), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%co")
      ENDIF
      IF (opts%rt_all%ch4_data) THEN
        ALLOCATE (profiles_int(j)%ch4(nlevels), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%ch4")
      ENDIF
      IF (opts%rt_all%so2_data) THEN
        ALLOCATE (profiles_int(j)%so2(nlevels), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%so2")
      ENDIF

      IF (opts%rt_ir%addaerosl .AND. .NOT. opts%rt_ir%user_aer_opt_param) THEN
        ALLOCATE (profiles_int(j)%aerosols(naertyp, nlayers), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%aerosols")
      ENDIF
      IF (opts%rt_ir%addclouds .AND. .NOT. opts%rt_ir%user_cld_opt_param) THEN
        ALLOCATE (profiles_int(j)%cloud(ncldtyp, nlayers), STAT = err)
        THROWM(err .NE. 0, "allocation of profiles%cloud")
      ENDIF
    ENDDO

  ENDIF ! asw == 1

  IF (asw == 0) THEN
    DO j = 1, nprofiles
      DEALLOCATE (profiles_int(j)%q, STAT = err)
      THROWM(err .NE. 0, "deallocation of profiles%q")
      IF (ASSOCIATED(profiles_int(j)%o3)) THEN
        DEALLOCATE (profiles_int(j)%o3, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%o3")
      ENDIF
      IF (ASSOCIATED(profiles_int(j)%co2)) THEN
        DEALLOCATE (profiles_int(j)%co2, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%co2")
      ENDIF
      IF (ASSOCIATED(profiles_int(j)%n2o)) THEN
        DEALLOCATE (profiles_int(j)%n2o, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%n2o")
      ENDIF
      IF (ASSOCIATED(profiles_int(j)%co)) THEN
        DEALLOCATE (profiles_int(j)%co, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%co")
      ENDIF
      IF (ASSOCIATED(profiles_int(j)%ch4)) THEN
        DEALLOCATE (profiles_int(j)%ch4, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%ch4")
      ENDIF
      IF (ASSOCIATED(profiles_int(j)%so2)) THEN
        DEALLOCATE (profiles_int(j)%so2, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%so2")
      ENDIF
      IF (ASSOCIATED(profiles_int(j)%aerosols)) THEN
        DEALLOCATE (profiles_int(j)%aerosols, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%aerosols")
      ENDIF
      IF (ASSOCIATED(profiles_int(j)%cloud)) THEN
        DEALLOCATE (profiles_int(j)%cloud, STAT = err)
        THROWM(err .NE. 0, "deallocation of profiles%cloud")
      ENDIF
      NULLIFY (profiles_int(j)%p)
      NULLIFY (profiles_int(j)%t)
      NULLIFY (profiles_int(j)%q)
      NULLIFY (profiles_int(j)%o3)
      NULLIFY (profiles_int(j)%co2)
      NULLIFY (profiles_int(j)%n2o)
      NULLIFY (profiles_int(j)%co)
      NULLIFY (profiles_int(j)%ch4)
      NULLIFY (profiles_int(j)%so2)
      NULLIFY (profiles_int(j)%clw)
      NULLIFY (profiles_int(j)%aerosols)
      NULLIFY (profiles_int(j)%cloud)
      NULLIFY (profiles_int(j)%cfrac)
      NULLIFY (profiles_int(j)%clwde)
      NULLIFY (profiles_int(j)%icede)
    ENDDO
  ENDIF ! asw == 0
  CATCH
END SUBROUTINE rttov_alloc_prof_internal
