! Description:
!> @file
!!   Initialise a VIS/IR cloud optical properties structure.
!!   This should usually be called via rttov_init_coefs.
!
!> @brief
!!   Initialise a VIS/IR cloud optical properties structure.
!!   This should usually be called via rttov_init_coefs.
!!
!! @details
!!   This subroutine precalculates some data related to the Baran
!!   ice scheme for VIS/IR cloud scattering simulations.
!!
!! @param[out]     err          status on exit
!! @param[in]      coef         the optical depth coefficient structure
!! @param[in,out]  coef_scatt   the cloud/aerosol optical properties structure to initialise
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
SUBROUTINE rttov_init_coef_scatt(err, coef, coef_scatt)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef, rttov_coef_scatt
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE mod_rttov_baran2014_icldata, ONLY : baran2014_wvn, n_baran2014_wn
  USE mod_rttov_baran2018_icldata, ONLY : baran2018_wvn, n_baran2018_wn
  USE parkind1, ONLY : jprb, jplm
  USE rttov_const, ONLY : phangle_hires, baran_ngauss
  USE rttov_scattering_mod, ONLY : gauss_quad
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),          INTENT(OUT)   :: err
  TYPE(rttov_coef),       INTENT(IN)    :: coef
  TYPE(rttov_coef_scatt), INTENT(INOUT) :: coef_scatt
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_alloc_phfn_int.interface"

!- End of header --------------------------------------------------------
  TRY

  ! For the phase functions use the Baran 2018 arrays for Baran 2014
  CALL init_optp_baran(err, coef_scatt%optp_icl_baran2014, baran2014_wvn, n_baran2014_wn, .FALSE._jplm)
  THROWM(err.NE.0, 'Error initialising Baran 2014 data')

  CALL init_optp_baran(err, coef_scatt%optp_icl_baran2018, baran2018_wvn, n_baran2018_wn, .TRUE._jplm)
  THROWM(err.NE.0, 'Error initialising Baran 2018 data')

  CATCH

CONTAINS

  SUBROUTINE init_optp_baran(err, optp_icl_baran, baran_wvn, n_baran_wn, init_phfn_int)

    USE rttov_types, ONLY : rttov_optp_baran

    INTEGER(jpim), INTENT(OUT) :: err
    TYPE(rttov_optp_baran), INTENT(INOUT) :: optp_icl_baran
    REAL(jprb),             INTENT(IN)    :: baran_wvn(:)
    INTEGER(jpim),          INTENT(IN)    :: n_baran_wn
    LOGICAL(jplm),          INTENT(IN)    :: init_phfn_int

    INTEGER(jpim) :: ichn, iwn, jwn
    REAL(jprb)    :: dx_dwn

    TRY

    ALLOCATE (optp_icl_baran%iwn(coef%fmv_chn), STAT=err)
    THROWM(err.NE.0, "allocation of optp_icl_baran%iwn" )

    ALLOCATE (optp_icl_baran%jwn(coef%fmv_chn), STAT=err)
    THROWM(err.NE.0, "allocation of optp_icl_baran%jwn" )

    ALLOCATE (optp_icl_baran%dx_dwn(coef%fmv_chn), STAT=err)
    THROWM(err.NE.0, "allocation of optp_icl_baran%dx_dwn" )

    DO ichn = 1, coef%fmv_chn

      IF (baran_wvn(1) .GE. coef%ff_cwn(ichn)) THEN
        iwn = 1_jpim
        jwn = 1_jpim
        dx_dwn  = 0.0_jprb
      ELSEIF (baran_wvn(n_baran_wn) .LE. coef%ff_cwn(ichn)) THEN
        iwn = n_baran_wn
        jwn = n_baran_wn
        dx_dwn  = 0.0_jprb
      ELSE
        iwn = 1_jpim
        DO WHILE (baran_wvn(iwn) .LE. coef%ff_cwn(ichn))
         iwn = iwn + 1_jpim
        ENDDO
        iwn = iwn - 1_jpim
        jwn = iwn + 1_jpim
        dx_dwn  = (coef%ff_cwn(ichn) - baran_wvn(iwn)) / (baran_wvn(jwn)  - baran_wvn(iwn))
      ENDIF

      optp_icl_baran%iwn(ichn) = iwn
      optp_icl_baran%jwn(ichn) = jwn
      optp_icl_baran%dx_dwn(ichn) = dx_dwn

      ! Interpolations will be done like this:
      ! value = value(iwn) + ( value(jwn) - value(iwn) ) * dx_dwn

    ENDDO

    IF (init_phfn_int) THEN
      CALL rttov_alloc_phfn_int(err, phangle_hires, optp_icl_baran%phfn_int, 1_jpim)
      THROW(err.NE.0)

      ALLOCATE(optp_icl_baran%q(baran_ngauss), optp_icl_baran%w(baran_ngauss), STAT=err)
      THROWM(err.NE.0, "allocation of Baran Gaussian quadrature arrays")

      CALL gauss_quad(-1._jprb, 1._jprb, optp_icl_baran%q, optp_icl_baran%w)
    ENDIF

    CATCH
  END SUBROUTINE init_optp_baran
END SUBROUTINE rttov_init_coef_scatt
