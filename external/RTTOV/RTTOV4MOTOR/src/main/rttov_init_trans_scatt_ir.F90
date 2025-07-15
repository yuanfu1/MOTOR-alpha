! Description:
!> @file
!!   Initialise internal transmission_scatt_ir structure.
!
!> @brief
!!   Initialise internal transmission_scatt_ir structure.
!!
!! @param[in,out] trans_scatt_ir  transmission_scatt_ir structure to initialise
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
SUBROUTINE rttov_init_trans_scatt_ir(trans_scatt_ir)

  USE rttov_types, ONLY : rttov_transmission_scatt_ir
!INTF_OFF
  USE rttov_types, ONLY : rttov_scatt_ir_aercld
  USE parkind1, ONLY : jprb, jpim
!INTF_ON
  IMPLICIT NONE

  TYPE(rttov_transmission_scatt_ir), INTENT(INOUT) :: trans_scatt_ir
!INTF_END
  INTEGER(jpim) :: i

  IF (ASSOCIATED(trans_scatt_ir%phasefn)) THEN
    DO i = 1, SIZE(trans_scatt_ir%phasefn)
      IF (ASSOCIATED(trans_scatt_ir%phasefn(i)%legcoef)) trans_scatt_ir%phasefn(i)%legcoef = 0._jprb
    ENDDO
  ENDIF
  IF (ASSOCIATED(trans_scatt_ir%opdpacl))     trans_scatt_ir%opdpacl     = 0._jprb
  IF (ASSOCIATED(trans_scatt_ir%opdpac))      trans_scatt_ir%opdpac      = 0._jprb
  IF (ASSOCIATED(trans_scatt_ir%opdpaclsun))  trans_scatt_ir%opdpaclsun  = 0._jprb
  IF (ASSOCIATED(trans_scatt_ir%opdpacsun))   trans_scatt_ir%opdpacsun   = 0._jprb
  IF (ASSOCIATED(trans_scatt_ir%opdpabs))     trans_scatt_ir%opdpabs     = 0._jprb
  IF (ASSOCIATED(trans_scatt_ir%opdpsca))     trans_scatt_ir%opdpsca     = 0._jprb
  IF (ASSOCIATED(trans_scatt_ir%ssa_solar))   trans_scatt_ir%ssa_solar   = 0._jprb
  IF (ASSOCIATED(trans_scatt_ir%ssa_thermal)) trans_scatt_ir%ssa_thermal = 0._jprb
  IF (ASSOCIATED(trans_scatt_ir%opdpext))     trans_scatt_ir%opdpext     = 0._jprb
  IF (ASSOCIATED(trans_scatt_ir%phup))        trans_scatt_ir%phup        = 0._jprb
  IF (ASSOCIATED(trans_scatt_ir%phdo))        trans_scatt_ir%phdo        = 0._jprb
  IF (ASSOCIATED(trans_scatt_ir%ray_sca))     trans_scatt_ir%ray_sca     = 0._jprb

  IF (ASSOCIATED(trans_scatt_ir%aer)) CALL init_scatt_ir_aercld(trans_scatt_ir%aer)
  IF (ASSOCIATED(trans_scatt_ir%cld)) CALL init_scatt_ir_aercld(trans_scatt_ir%cld)

CONTAINS

  SUBROUTINE init_scatt_ir_aercld(aercld)
    TYPE(rttov_scatt_ir_aercld), INTENT(INOUT) :: aercld

    IF (ASSOCIATED(aercld%opdpabs))    aercld%opdpabs    = 0._jprb
    IF (ASSOCIATED(aercld%opdpsca))    aercld%opdpsca    = 0._jprb
    IF (ASSOCIATED(aercld%opdpscabpr)) aercld%opdpscabpr = 0._jprb
    IF (ASSOCIATED(aercld%opdp))       aercld%opdp       = 0._jprb
    IF (ASSOCIATED(aercld%opdpsun))    aercld%opdpsun    = 0._jprb
    IF (ASSOCIATED(aercld%partsca))    aercld%partsca    = 0._jprb
    IF (ASSOCIATED(aercld%partbpr))    aercld%partbpr    = 0._jprb
    IF (ASSOCIATED(aercld%sca))        aercld%sca        = 0._jprb
    IF (ASSOCIATED(aercld%phintup))    aercld%phintup    = 0._jprb
    IF (ASSOCIATED(aercld%phintdo))    aercld%phintdo    = 0._jprb
    IF (ASSOCIATED(aercld%phtotup))    aercld%phtotup    = 0._jprb
    IF (ASSOCIATED(aercld%phtotdo))    aercld%phtotdo    = 0._jprb

  END SUBROUTINE init_scatt_ir_aercld

END SUBROUTINE
