! Description:
!> @file
!!   Allocate/deallocate internal predictors structure.
!
!> @brief
!!   Allocate/deallocate internal predictors structure.
!!
!! @details
!!   This structure holds the values of the predictors.
!!
!! @param[out]    err                    status on exit
!! @param[in]     nprofiles              number of input profiles
!! @param[in,out] predictors             predictors structure to (de)allocate
!! @param[in]     coef                   optical depth coefficient structure
!! @param[in]     asw                    1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     addsolar               flag to indicate if solar radiation is turned on
!! @param[in]     init                   set .TRUE. to initialise newly allocated structures, optional
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
SUBROUTINE rttov_alloc_predictor( &
              err,         &
              nprofiles,   &
              predictors,  &
              coef,        &
              asw,         &
              addsolar,    &
              init)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef, rttov_predictors
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE

  INTEGER(KIND=jpim),     INTENT(OUT)            :: err
  INTEGER(KIND=jpim),     INTENT(IN)             :: nprofiles
  TYPE(rttov_predictors), INTENT(INOUT)          :: predictors
  TYPE(rttov_coef),       INTENT(IN)             :: coef
  INTEGER(KIND=jpim),     INTENT(IN)             :: asw         ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm),     INTENT(IN)             :: addsolar
  LOGICAL(KIND=jplm),     INTENT(IN),   OPTIONAL :: init
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_init_predictor.interface"

  INTEGER(KIND=jpim) :: nlayers, iprof
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------
  TRY
  init1 = .FALSE.
  IF (PRESENT(init)) init1 = init
  nlayers = coef%nlayers

  IF (asw == 1) THEN
    predictors%nlevels = coef%nlevels
    predictors%nmixed  = coef%nmixed
    predictors%nwater  = coef%nwater
    predictors%nozone  = coef%nozone
    predictors%nwvcont = coef%nwvcont
    predictors%nco2    = coef%nco2
    predictors%nn2o    = coef%nn2o
    predictors%nco     = coef%nco
    predictors%nch4    = coef%nch4
    predictors%nso2    = coef%nso2
    predictors%npmc    = coef%pmc_nvar

    ALLOCATE( predictors%path1(nprofiles))
    If (addsolar) ALLOCATE( predictors%path2(nprofiles))
    CALL nullify_struct()

    DO iprof = 1, nprofiles
      ALLOCATE (predictors%path1(iprof)%mixedgas(coef%nmixed, nlayers), STAT = err)
      THROWM(err.NE.0, "allocation of predictors%path1(iprof)%mixedgas")
      ALLOCATE (predictors%path1(iprof)%watervapour(coef%nwater, nlayers), STAT = err)
      THROWM(err.NE.0, "allocation of predictors%path1(iprof)%watervapour")
      IF (coef%nozone > 0) THEN
        ALLOCATE (predictors%path1(iprof)%ozone(coef%nozone, nlayers), STAT = err)
        THROWM(err.NE.0, "allocation of predictors%path1(iprof)%ozone")
      ENDIF
      IF (coef%nwvcont > 0) THEN
        ALLOCATE (predictors%path1(iprof)%wvcont(coef%nwvcont, nlayers), STAT = err)
        THROWM(err.NE.0, "allocation of predictors%path1(iprof)%wvcont")
      ENDIF
      IF (coef%nco2 > 0) THEN
        ALLOCATE (predictors%path1(iprof)%co2(coef%nco2, nlayers), STAT = err)
        THROWM(err.NE.0, "allocation of predictors%path1(iprof)%co2")
      ENDIF
      IF (coef%nn2o > 0) THEN
        ALLOCATE (predictors%path1(iprof)%n2o(coef%nn2o, nlayers), STAT = err)
        THROWM(err.NE.0, "allocation of predictors%path1(iprof)%n2o")
      ENDIF
      IF (coef%nco > 0) THEN
        ALLOCATE (predictors%path1(iprof)%co(coef%nco + 1, nlayers), STAT = err)
        THROWM(err.NE.0, "allocation of predictors%path1(iprof)%co")
      ENDIF
      IF (coef%nch4 > 0) THEN
        ALLOCATE (predictors%path1(iprof)%ch4(coef%nch4, nlayers), STAT = err)
        THROWM(err.NE.0, "allocation of predictors%path1(iprof)%ch4")
      ENDIF
      IF (coef%nso2 > 0) THEN
        ALLOCATE (predictors%path1(iprof)%so2(coef%nso2, nlayers), STAT = err)
        THROWM(err.NE.0, "allocation of predictors%path1(iprof)%so2")
      ENDIF
      IF (coef%pmc_shift) THEN
        ALLOCATE (predictors%path1(iprof)%pmc(coef%pmc_nvar, nlayers+1, coef%fmv_chn), STAT = err)
        THROWM(err.NE.0, "allocation of predictors%path1(iprof)%pmc")
      ENDIF

      IF (addsolar) THEN
        ALLOCATE (predictors%path2(iprof)%mixedgas(coef%nmixed, nlayers), STAT = err)
        THROWM(err.NE.0, "allocation of predictors%path2(iprof)%mixedgas")
        ALLOCATE (predictors%path2(iprof)%watervapour(coef%nwater, nlayers), STAT = err)
        THROWM(err.NE.0, "allocation of predictors%path2(iprof)%watervapour")
        IF (coef%nozone > 0) THEN
          ALLOCATE (predictors%path2(iprof)%ozone(coef%nozone, nlayers), STAT = err)
          THROWM(err.NE.0, "allocation of predictors%path2(iprof)%ozone")
        ENDIF
        IF (coef%nwvcont > 0) THEN
          ALLOCATE (predictors%path2(iprof)%wvcont(coef%nwvcont, nlayers), STAT = err)
          THROWM(err.NE.0, "allocation of predictors%path2(iprof)%wvcont")
        ENDIF
        IF (coef%nco2 > 0) THEN
          ALLOCATE (predictors%path2(iprof)%co2(coef%nco2, nlayers), STAT = err)
          THROWM(err.NE.0, "allocation of predictors%path2(iprof)%co2")
        ENDIF
        IF (coef%nn2o > 0) THEN
          ALLOCATE (predictors%path2(iprof)%n2o(coef%nn2o, nlayers), STAT = err)
          THROWM(err.NE.0, "allocation of predictors%path2(iprof)%n2o")
        ENDIF
        IF (coef%nco > 0) THEN
          ALLOCATE (predictors%path2(iprof)%co(coef%nco + 1, nlayers), STAT = err)
          THROWM(err.NE.0, "allocation of predictors%path2(iprof)%co")
        ENDIF
        IF (coef%nch4 > 0) THEN
          ALLOCATE (predictors%path2(iprof)%ch4(coef%nch4, nlayers), STAT = err)
          THROWM(err.NE.0, "allocation of predictors%path2(iprof)%ch4")
        ENDIF
        IF (coef%nso2 > 0) THEN
          ALLOCATE (predictors%path2(iprof)%so2(coef%nso2, nlayers), STAT = err)
          THROWM(err.NE.0, "allocation of predictors%path2(iprof)%so2")
        ENDIF
      ENDIF
    ENDDO

    IF (init1) CALL rttov_init_predictor(addsolar, predictors)
  ENDIF

  IF (asw == 0) THEN
    DO iprof = 1, nprofiles
      DEALLOCATE (predictors%path1(iprof)%mixedgas, STAT = err)
      THROWM(err.NE.0, "deallocation of predictors%path1(iprof)%mixedgas")
      DEALLOCATE (predictors%path1(iprof)%watervapour, STAT = err)
      THROWM(err.NE.0, "deallocation of predictors%path1(iprof)%watervapour")
      IF (coef%nozone > 0) THEN
        DEALLOCATE (predictors%path1(iprof)%ozone, STAT = err)
        THROWM(err.NE.0, "deallocation of predictors%path1(iprof)%ozone")
      ENDIF
      IF (coef%nwvcont > 0) THEN
        DEALLOCATE (predictors%path1(iprof)%wvcont, STAT = err)
        THROWM(err.NE.0, "deallocation of predictors%path1(iprof)%wvcont")
      ENDIF
      IF (coef%nco2 > 0) THEN
        DEALLOCATE (predictors%path1(iprof)%co2, STAT = err)
        THROWM(err.NE.0, "deallocation of predictors%path1(iprof)%co2")
      ENDIF
      IF (coef%nn2o > 0) THEN
        DEALLOCATE (predictors%path1(iprof)%n2o, STAT = err)
        THROWM(err.NE.0, "deallocation of predictors%path1(iprof)%n2o")
      ENDIF
      IF (coef%nco > 0) THEN
        DEALLOCATE (predictors%path1(iprof)%co, STAT = err)
        THROWM(err.NE.0, "deallocation of predictors%path1(iprof)%co")
      ENDIF
      IF (coef%nch4 > 0) THEN
        DEALLOCATE (predictors%path1(iprof)%ch4, STAT = err)
        THROWM(err.NE.0, "deallocation of predictors%path1(iprof)%ch4")
      ENDIF
      IF (coef%nso2 > 0) THEN
        DEALLOCATE (predictors%path1(iprof)%so2, STAT = err)
        THROWM(err.NE.0, "deallocation of predictors%path1(iprof)%so2")
      ENDIF
      IF (coef%pmc_shift) THEN
        DEALLOCATE (predictors%path1(iprof)%pmc, STAT = err)
        THROWM(err.NE.0, "deallocation of predictors%path1(iprof)%pmc")
      ENDIF

      IF (addsolar) THEN
        DEALLOCATE (predictors%path2(iprof)%mixedgas, STAT = err)
        THROWM(err.NE.0, "deallocation of predictors%path2(iprof)%mixedgas")
        DEALLOCATE (predictors%path2(iprof)%watervapour, STAT = err)
        THROWM(err.NE.0, "deallocation of predictors%path2(iprof)%watervapour")
        IF (coef%nozone > 0) THEN
          DEALLOCATE (predictors%path2(iprof)%ozone, STAT = err)
          THROWM(err.NE.0, "deallocation of predictors%path2(iprof)%ozone")
        ENDIF
        IF (coef%nwvcont > 0) THEN
          DEALLOCATE (predictors%path2(iprof)%wvcont, STAT = err)
          THROWM(err.NE.0, "deallocation of predictors%path2(iprof)%wvcont")
        ENDIF
        IF (coef%nco2 > 0) THEN
          DEALLOCATE (predictors%path2(iprof)%co2, STAT = err)
          THROWM(err.NE.0, "deallocation of predictors%path2(iprof)%co2")
        ENDIF
        IF (coef%nn2o > 0) THEN
          DEALLOCATE (predictors%path2(iprof)%n2o, STAT = err)
          THROWM(err.NE.0, "deallocation of predictors%path2(iprof)%n2o")
        ENDIF
        IF (coef%nco > 0) THEN
          DEALLOCATE (predictors%path2(iprof)%co, STAT = err)
          THROWM(err.NE.0, "deallocation of predictors%path2(iprof)%co")
        ENDIF
        IF (coef%nch4 > 0) THEN
          DEALLOCATE (predictors%path2(iprof)%ch4, STAT = err)
          THROWM(err.NE.0, "deallocation of predictors%path2(iprof)%ch4")
        ENDIF
        IF (coef%nso2 > 0) THEN
          DEALLOCATE (predictors%path2(iprof)%so2, STAT = err)
          THROWM(err.NE.0, "deallocation of predictors%path2(iprof)%so2")
        ENDIF
      ENDIF
    ENDDO

    CALL nullify_struct()
    DEALLOCATE(predictors%path1)
    IF (addsolar) DEALLOCATE(predictors%path2)

  ENDIF
  CATCH

CONTAINS

  SUBROUTINE nullify_struct()
    DO iprof = 1, nprofiles
      NULLIFY (predictors%path1(iprof)%mixedgas)
      NULLIFY (predictors%path1(iprof)%watervapour)
      NULLIFY (predictors%path1(iprof)%ozone)
      NULLIFY (predictors%path1(iprof)%wvcont)
      NULLIFY (predictors%path1(iprof)%co2)
      NULLIFY (predictors%path1(iprof)%n2o)
      NULLIFY (predictors%path1(iprof)%co)
      NULLIFY (predictors%path1(iprof)%ch4)
      NULLIFY (predictors%path1(iprof)%so2)
      NULLIFY (predictors%path1(iprof)%pmc)
      IF(addsolar) THEN
        NULLIFY (predictors%path2(iprof)%mixedgas)
        NULLIFY (predictors%path2(iprof)%watervapour)
        NULLIFY (predictors%path2(iprof)%ozone)
        NULLIFY (predictors%path2(iprof)%wvcont)
        NULLIFY (predictors%path2(iprof)%co2)
        NULLIFY (predictors%path2(iprof)%n2o)
        NULLIFY (predictors%path2(iprof)%co)
        NULLIFY (predictors%path2(iprof)%ch4)
        NULLIFY (predictors%path2(iprof)%so2)
      ENDIF
    ENDDO
  END SUBROUTINE nullify_struct
END SUBROUTINE rttov_alloc_predictor
