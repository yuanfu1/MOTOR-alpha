! Description:
!> @file
!!   Nullify/zero the VIS/IR aerosol/cloud optical properties structure.
!
!> @brief
!!   Nullify/zero the VIS/IR aerosol/cloud optical properties structure.
!!
!!
!! @param[in,out]  coef_scatt  the cloud/aerosol optical properties structure to nullify/zero
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
SUBROUTINE rttov_nullify_coef_scatt(coef_scatt)
!INTF_OFF
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  USE rttov_types, ONLY : rttov_coef_scatt
  IMPLICIT NONE
  TYPE(rttov_coef_scatt), INTENT(INOUT) :: coef_scatt
!INTF_END

  CALL nullify_optp(coef_scatt%optp_aer)
  CALL nullify_optp(coef_scatt%optp_wcl_opac)
  CALL nullify_optp(coef_scatt%optp_wcl_deff)
  CALL nullify_optp(coef_scatt%optp_icl_baum)
  CALL nullify_optp_baran(coef_scatt%optp_icl_baran2014)
  CALL nullify_optp_baran(coef_scatt%optp_icl_baran2018)

CONTAINS

  SUBROUTINE nullify_phfn_int(phfn_int)
    USE rttov_types, ONLY : rttov_phasefn_int

    TYPE(rttov_phasefn_int), INTENT(INOUT) :: phfn_int

    phfn_int%zminphadiff = 0._jprb
    NULLIFY(phfn_int%iphangle)
    NULLIFY(phfn_int%cosphangle)

  END SUBROUTINE nullify_phfn_int

  SUBROUTINE nullify_optp(optp)
    USE rttov_types, ONLY : rttov_optp

    TYPE(rttov_optp), INTENT(INOUT) :: optp

    optp%version   = 0_jpim
    optp%id        = 0_jpim
    optp%nchan     = 0_jpim
    optp%nchan_pha = 0_jpim
    optp%maxnmom   = 0_jpim
    optp%nphangle  = 0_jpim
    optp%ntypes    = 0_jpim

    NULLIFY(optp%chan_pha)
    NULLIFY(optp%chan_pha_index)
    NULLIFY(optp%phangle)
    NULLIFY(optp%data)
    CALL nullify_phfn_int(optp%phfn_int)

  END SUBROUTINE nullify_optp

  SUBROUTINE nullify_optp_baran(optp)
    USE rttov_types, ONLY : rttov_optp_baran

    TYPE(rttov_optp_baran), INTENT(INOUT) :: optp

    NULLIFY(optp%iwn)
    NULLIFY(optp%jwn)
    NULLIFY(optp%dx_dwn)
    NULLIFY(optp%q)
    NULLIFY(optp%w)
    CALL nullify_phfn_int(optp%phfn_int)

  END SUBROUTINE nullify_optp_baran

END SUBROUTINE rttov_nullify_coef_scatt
