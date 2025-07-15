! Description:
!> @file
!!   Nullify/zero an HTFRTC coefficients structure.
!
!> @brief
!!   Nullify/zero an HTFRTC coefficients structure.
!!
!!
!! @param[in,out]  coef_htfrtc     the HTFRTC coefficients structure to nullify/zero
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
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_nullify_coef_htfrtc(coef_htfrtc)
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON
  USE rttov_types, ONLY : rttov_coef_htfrtc
  IMPLICIT NONE
  TYPE(rttov_coef_htfrtc), INTENT(INOUT) :: coef_htfrtc
!INTF_END

  coef_htfrtc%n_f           = 0_jpim
  coef_htfrtc%n_gas_l       = 0_jpim
  coef_htfrtc%n_p           = 0_jpim
  coef_htfrtc%n_val_l       = 0_jpim
  coef_htfrtc%n_b           = 0_jpim
  coef_htfrtc%n_lt          = 0_jpim
  coef_htfrtc%n_cont        = 0_jpim
  coef_htfrtc%n_ssemp       = 0_jpim
  coef_htfrtc%n_iremis      = 0_jpim
  coef_htfrtc%n_pc          = 0_jpim
  coef_htfrtc%n_pc_oc       = 0_jpim
  coef_htfrtc%n_ch          = 0_jpim
  coef_htfrtc%n_mftlb       = 0_jpim

  NULLIFY (coef_htfrtc%freq)
  NULLIFY (coef_htfrtc%gasid_l)
  NULLIFY (coef_htfrtc%p)
  NULLIFY (coef_htfrtc%val_b)
  NULLIFY (coef_htfrtc%val_lt)
  NULLIFY (coef_htfrtc%coef_l)
  NULLIFY (coef_htfrtc%coef_ct)
  NULLIFY (coef_htfrtc%coef_ctt)
  NULLIFY (coef_htfrtc%coef_b)
  NULLIFY (coef_htfrtc%coef_lt)
  NULLIFY (coef_htfrtc%coef_ssemp)
  NULLIFY (coef_htfrtc%coef_iremis)
  NULLIFY (coef_htfrtc%coef_pdt)
  NULLIFY (coef_htfrtc%val_mean)
  NULLIFY (coef_htfrtc%val_norm)
  NULLIFY (coef_htfrtc%sensor_freq)
  NULLIFY (coef_htfrtc%ch_mean)
  NULLIFY (coef_htfrtc%pc)
  NULLIFY (coef_htfrtc%mixed_ref_frac)
  NULLIFY (coef_htfrtc%mftlb)
  NULLIFY (coef_htfrtc%addf)
  NULLIFY (coef_htfrtc%addch)
END SUBROUTINE
