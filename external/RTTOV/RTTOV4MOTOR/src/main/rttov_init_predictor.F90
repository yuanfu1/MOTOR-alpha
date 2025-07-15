! Description:
!> @file
!!   Initialise internal predictors structure.
!
!> @brief
!!   Initialise internal predictors structure.
!!
!! @param[in]     addsolar            flag to indicate if solar radiation is turned on
!! @param[in,out] predictors          predictors structure to initialise
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
SUBROUTINE rttov_init_predictor(addsolar, predictors)

  USE rttov_types, ONLY : rttov_predictors
  USE parkind1,    ONLY : jplm
!INTF_OFF
  USE parkind1, ONLY : jprb, jpim
!INTF_ON
  IMPLICIT NONE
  LOGICAL(KIND=jplm),     INTENT(IN)    :: addsolar
  TYPE(rttov_predictors), INTENT(INOUT) :: predictors
!INTF_END
  INTEGER(KIND=jpim) :: iprof

  IF (ASSOCIATED(predictors%path1)) THEN
    DO iprof = LBOUND(predictors%path1, DIM=1), UBOUND(predictors%path1, DIM=1)
      IF (ASSOCIATED(predictors%path1(iprof)%mixedgas)   ) predictors%path1(iprof)%mixedgas    = 0._jprb
      IF (ASSOCIATED(predictors%path1(iprof)%watervapour)) predictors%path1(iprof)%watervapour = 0._jprb
      IF (ASSOCIATED(predictors%path1(iprof)%ozone)      ) predictors%path1(iprof)%ozone       = 0._jprb
      IF (ASSOCIATED(predictors%path1(iprof)%wvcont)     ) predictors%path1(iprof)%wvcont      = 0._jprb
      IF (ASSOCIATED(predictors%path1(iprof)%co2)        ) predictors%path1(iprof)%co2         = 0._jprb
      IF (ASSOCIATED(predictors%path1(iprof)%n2o)        ) predictors%path1(iprof)%n2o         = 0._jprb
      IF (ASSOCIATED(predictors%path1(iprof)%co)         ) predictors%path1(iprof)%co          = 0._jprb
      IF (ASSOCIATED(predictors%path1(iprof)%ch4)        ) predictors%path1(iprof)%ch4         = 0._jprb
      IF (ASSOCIATED(predictors%path1(iprof)%so2)        ) predictors%path1(iprof)%so2         = 0._jprb
      IF (ASSOCIATED(predictors%path1(iprof)%pmc)        ) predictors%path1(iprof)%pmc         = 0._jprb
    ENDDO
  ENDIF

  IF (addsolar) THEN
    IF (ASSOCIATED(predictors%path2)) THEN
      DO iprof = LBOUND(predictors%path2, DIM=1), UBOUND(predictors%path2, DIM=1)
        IF (ASSOCIATED(predictors%path2(iprof)%mixedgas)   ) predictors%path2(iprof)%mixedgas    = 0._jprb
        IF (ASSOCIATED(predictors%path2(iprof)%watervapour)) predictors%path2(iprof)%watervapour = 0._jprb
        IF (ASSOCIATED(predictors%path2(iprof)%ozone)      ) predictors%path2(iprof)%ozone       = 0._jprb
        IF (ASSOCIATED(predictors%path2(iprof)%wvcont)     ) predictors%path2(iprof)%wvcont      = 0._jprb
        IF (ASSOCIATED(predictors%path2(iprof)%co2)        ) predictors%path2(iprof)%co2         = 0._jprb
        IF (ASSOCIATED(predictors%path2(iprof)%n2o)        ) predictors%path2(iprof)%n2o         = 0._jprb
        IF (ASSOCIATED(predictors%path2(iprof)%co)         ) predictors%path2(iprof)%co          = 0._jprb
        IF (ASSOCIATED(predictors%path2(iprof)%ch4)        ) predictors%path2(iprof)%ch4         = 0._jprb
        IF (ASSOCIATED(predictors%path2(iprof)%so2)        ) predictors%path2(iprof)%so2         = 0._jprb
      ENDDO
    ENDIF
  ENDIF
END SUBROUTINE rttov_init_predictor

