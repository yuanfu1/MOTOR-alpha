! Description:
!> @file
!!   Extract data for given channel list from a PC coefficients structure.
!
!> @brief
!!   Extract data for given channel list from a PC coefficients structure.
!!
!! @details
!!   This is used by HDF5 I/O code to read in a subset of channels from a
!!   coefficient file. The first coef argument contains the coefficients
!!   from the file. The second argument is an uninitialised structure
!!   which contains the extracted coefficients on exit.
!!
!!   NB With PC-RTTOV coefficients the channel list MUST be a superset of
!!   predictor channels for the band/predictor set being used.
!!
!!   If channels_rec is not supplied coefficients are extracted to enable
!!   reconstructed radiances for all channels.
!!
!!   It is possible to extract a subset channels_rec for reconstructed
!!   radiances while reading all channels for the PC predictor calculations.
!!
!! @param[out]     err            status on exit
!! @param[in]      coef_pccomp1   input coefficients read from file
!! @param[in,out]  coef_pccomp2   output coefficients, uninitialised on entry
!! @param[in]      channels       list of channels to extract, optional
!! @param[in]      channels_rec   list of channels to extract for reconstructed radiances, optional
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
SUBROUTINE rttov_channel_extract_pccoef(err, coef_pccomp1, coef_pccomp2, channels, channels_rec)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY : rttov_coef_pccomp
!INTF_OFF
  USE parkind1, ONLY : jplm
!INTF_ON
  IMPLICIT NONE

  INTEGER(jpim),           INTENT(OUT)   :: err
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp1
  TYPE(rttov_coef_pccomp), INTENT(INOUT) :: coef_pccomp2
  INTEGER(jpim), OPTIONAL, INTENT(IN)    :: channels(:)
  INTEGER(jpim), OPTIONAL, INTENT(IN)    :: channels_rec(:)
!INTF_END

#include "rttov_errorreport.interface"

  LOGICAL(jplm) :: all_channels, all_channels_rec
  INTEGER(jpim) :: m, n
! ----------------------------------------------------------------------------

TRY

  all_channels = .NOT. PRESENT(channels)
  all_channels_rec = .NOT. PRESENT(channels_rec)

  ! Scalar variables

  ! A few are different...

  IF (all_channels_rec) THEN
    coef_pccomp2%fmv_pc_nchn = coef_pccomp1%fmv_pc_nchn
  ELSE
    coef_pccomp2%fmv_pc_nchn = SIZE(channels_rec)
  ENDIF

  IF (all_channels) THEN
    coef_pccomp2%fmv_pc_nchn_noise = coef_pccomp1%fmv_pc_nchn_noise
    coef_pccomp2%fmv_pc_nche = coef_pccomp1%fmv_pc_nche
  ELSE
    coef_pccomp2%fmv_pc_nchn_noise = SIZE(channels)
    coef_pccomp2%fmv_pc_nche = SIZE(channels)
  ENDIF

  ! ...but the rest are the same

  coef_pccomp2%fmv_pc_comp_pc    = coef_pccomp1%fmv_pc_comp_pc
  coef_pccomp2%fmv_pc_cld        = coef_pccomp1%fmv_pc_cld
  coef_pccomp2%fmv_pc_aer        = coef_pccomp1%fmv_pc_aer
  coef_pccomp2%fmv_pc_naer_types = coef_pccomp1%fmv_pc_naer_types
  coef_pccomp2%fmv_pc_nlte       = coef_pccomp1%fmv_pc_nlte
  coef_pccomp2%fmv_pc_msets      = coef_pccomp1%fmv_pc_msets
  coef_pccomp2%fmv_pc_bands      = coef_pccomp1%fmv_pc_bands
  coef_pccomp2%fmv_pc_mnum       = coef_pccomp1%fmv_pc_mnum
  coef_pccomp2%fmv_pc_mchn       = coef_pccomp1%fmv_pc_mchn
  coef_pccomp2%fmv_pc_gas        = coef_pccomp1%fmv_pc_gas
  coef_pccomp2%fmv_pc_gas_lim    = coef_pccomp1%fmv_pc_gas_lim
  coef_pccomp2%fmv_pc_nlev       = coef_pccomp1%fmv_pc_nlev

  coef_pccomp2%lim_pc_prfl_pmin  = coef_pccomp1%lim_pc_prfl_pmin
  coef_pccomp2%lim_pc_prfl_pmax  = coef_pccomp1%lim_pc_prfl_pmax
  coef_pccomp2%lim_pc_prfl_tsmin = coef_pccomp1%lim_pc_prfl_tsmin
  coef_pccomp2%lim_pc_prfl_tsmax = coef_pccomp1%lim_pc_prfl_tsmax
  coef_pccomp2%lim_pc_prfl_skmin = coef_pccomp1%lim_pc_prfl_skmin
  coef_pccomp2%lim_pc_prfl_skmax = coef_pccomp1%lim_pc_prfl_skmax
  coef_pccomp2%lim_pc_prfl_wsmin = coef_pccomp1%lim_pc_prfl_wsmin
  coef_pccomp2%lim_pc_prfl_wsmax = coef_pccomp1%lim_pc_prfl_wsmax


  ! PC predictors

  ALLOCATE(coef_pccomp2%fmv_pc_sets(coef_pccomp2%fmv_pc_bands), stat=err)
  THROWM(err.NE.0, 'allocation of fmv_pc_sets')
  coef_pccomp2%fmv_pc_sets = coef_pccomp1%fmv_pc_sets

  ALLOCATE(coef_pccomp2%pcreg(coef_pccomp2%fmv_pc_bands,coef_pccomp2%fmv_pc_msets), stat=err)
  THROWM(err.NE.0, 'allocation of pcreg')

  DO m = 1, coef_pccomp2%fmv_pc_bands
    DO n = 1, coef_pccomp2%fmv_pc_sets(m)

      coef_pccomp2%pcreg(m,n)%fmv_pc_npred = coef_pccomp1%pcreg(m,n)%fmv_pc_npred

      ALLOCATE(coef_pccomp2%pcreg(m,n)%predictindex(coef_pccomp2%pcreg(m,n)%fmv_pc_npred), stat=err)
      THROWM(err.NE.0, 'allocation of pcreg(m,n)%predictindex')
      coef_pccomp2%pcreg(m,n)%predictindex = coef_pccomp1%pcreg(m,n)%predictindex

    ENDDO
  ENDDO


  ! Eigenvectors

  ALLOCATE(coef_pccomp2%eigen(coef_pccomp2%fmv_pc_bands), stat=err)
  THROWM(err.NE.0, 'allocation of eigen')

  DO m = 1, coef_pccomp2%fmv_pc_bands

    ALLOCATE(coef_pccomp2%eigen(m)%eigenvectors(coef_pccomp2%fmv_pc_nchn,coef_pccomp2%fmv_pc_mnum), stat=err)
    THROWM(err.NE.0, 'allocation of eigen(m)%eigenvectors')

    IF (all_channels_rec) THEN
      coef_pccomp2%eigen(m)%eigenvectors = coef_pccomp1%eigen(m)%eigenvectors
    ELSE
      coef_pccomp2%eigen(m)%eigenvectors = coef_pccomp1%eigen(m)%eigenvectors(channels_rec,:)
    ENDIF

    NULLIFY(coef_pccomp2%eigen(m)%eigenvectors_t)

  ENDDO


  ! PC coefficients

  DO m = 1, coef_pccomp2%fmv_pc_bands
    DO n = 1, coef_pccomp2%fmv_pc_sets(m)

      ALLOCATE(coef_pccomp2%pcreg(m,n)%coefficients(coef_pccomp2%pcreg(m,n)%fmv_pc_npred, &
                                                    coef_pccomp2%fmv_pc_mnum), stat=err)
      THROWM(err.NE.0, 'allocation of pcreg(m,n)%coefficients')
      coef_pccomp2%pcreg(m,n)%coefficients = coef_pccomp1%pcreg(m,n)%coefficients

      NULLIFY(coef_pccomp2%pcreg(m,n)%coefficients_t)

    ENDDO
  ENDDO


  ! Emissivity coefficients

  ALLOCATE(coef_pccomp2%emiss_chn(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err.NE.0, 'allocation of emiss_chn')
  IF (all_channels) THEN
    coef_pccomp2%emiss_chn = coef_pccomp1%emiss_chn
  ELSE
    coef_pccomp2%emiss_chn = coef_pccomp1%emiss_chn(channels)
  ENDIF

  ALLOCATE(coef_pccomp2%emiss_c1(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err.NE.0, 'allocation of emiss_c1')
  IF (all_channels) THEN
    coef_pccomp2%emiss_c1 = coef_pccomp1%emiss_c1
  ELSE
    coef_pccomp2%emiss_c1 = coef_pccomp1%emiss_c1(channels)
  ENDIF

  ALLOCATE(coef_pccomp2%emiss_c2(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err.NE.0, 'allocation of emiss_c2')
  IF (all_channels) THEN
    coef_pccomp2%emiss_c2 = coef_pccomp1%emiss_c2
  ELSE
    coef_pccomp2%emiss_c2 = coef_pccomp1%emiss_c2(channels)
  ENDIF

  ALLOCATE(coef_pccomp2%emiss_c3(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err.NE.0, 'allocation of emiss_c3')
  IF (all_channels) THEN
    coef_pccomp2%emiss_c3 = coef_pccomp1%emiss_c3
  ELSE
    coef_pccomp2%emiss_c3 = coef_pccomp1%emiss_c3(channels)
  ENDIF

  ALLOCATE(coef_pccomp2%emiss_c4(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err.NE.0, 'allocation of emiss_c4')
  IF (all_channels) THEN
    coef_pccomp2%emiss_c4 = coef_pccomp1%emiss_c4
  ELSE
    coef_pccomp2%emiss_c4 = coef_pccomp1%emiss_c4(channels)
  ENDIF

  ALLOCATE(coef_pccomp2%emiss_c5(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err.NE.0, 'allocation of emiss_c5')
  IF (all_channels) THEN
    coef_pccomp2%emiss_c5 = coef_pccomp1%emiss_c5
  ELSE
    coef_pccomp2%emiss_c5 = coef_pccomp1%emiss_c5(channels)
  ENDIF

  ALLOCATE(coef_pccomp2%emiss_c6(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err.NE.0, 'allocation of emiss_c6')
  IF (all_channels) THEN
    coef_pccomp2%emiss_c6 = coef_pccomp1%emiss_c6
  ELSE
    coef_pccomp2%emiss_c6 = coef_pccomp1%emiss_c6(channels)
  ENDIF

  ALLOCATE(coef_pccomp2%emiss_c7(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err.NE.0, 'allocation of emiss_c7')
  IF (all_channels) THEN
    coef_pccomp2%emiss_c7 = coef_pccomp1%emiss_c7
  ELSE
    coef_pccomp2%emiss_c7 = coef_pccomp1%emiss_c7(channels)
  ENDIF

  ALLOCATE(coef_pccomp2%emiss_c8(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err.NE.0, 'allocation of emiss_c8')
  IF (all_channels) THEN
    coef_pccomp2%emiss_c8 = coef_pccomp1%emiss_c8
  ELSE
    coef_pccomp2%emiss_c8 = coef_pccomp1%emiss_c8(channels)
  ENDIF

  ALLOCATE(coef_pccomp2%emiss_c9(coef_pccomp2%fmv_pc_nche), stat=err)
  THROWM(err.NE.0, 'allocation of emiss_c9')
  IF (all_channels) THEN
    coef_pccomp2%emiss_c9 = coef_pccomp1%emiss_c9
  ELSE
    coef_pccomp2%emiss_c9 = coef_pccomp1%emiss_c9(channels)
  ENDIF


  ! PC reference profiles and limits

  ALLOCATE(coef_pccomp2%ref_pc_prfl_p(coef_pccomp2%fmv_pc_nlev), stat=err)
  THROWM(err.NE.0, 'allocation of ref_pc_prfl_p')
  coef_pccomp2%ref_pc_prfl_p = coef_pccomp1%ref_pc_prfl_p

  IF (coef_pccomp1%fmv_pc_comp_pc < 5) THEN
    ALLOCATE(coef_pccomp2%ref_pc_prfl_mr(coef_pccomp2%fmv_pc_nlev,coef_pccomp2%fmv_pc_gas), stat=err)
    THROWM(err.NE.0, 'allocation of ref_pc_prfl_mr')
    coef_pccomp2%ref_pc_prfl_mr = coef_pccomp1%ref_pc_prfl_mr
  ENDIF

  ALLOCATE(coef_pccomp2%lim_pc_prfl_tmin(coef_pccomp2%fmv_pc_nlev), stat=err)
  THROWM(err.NE.0, 'allocation of lim_pc_prfl_tmin')
  coef_pccomp2%lim_pc_prfl_tmin = coef_pccomp1%lim_pc_prfl_tmin

  ALLOCATE(coef_pccomp2%lim_pc_prfl_tmax(coef_pccomp2%fmv_pc_nlev), stat=err)
  THROWM(err.NE.0, 'allocation of lim_pc_prfl_tmax')
  coef_pccomp2%lim_pc_prfl_tmax = coef_pccomp1%lim_pc_prfl_tmax

  ALLOCATE(coef_pccomp2%lim_pc_prfl_qmin(coef_pccomp2%fmv_pc_nlev), stat=err)
  THROWM(err.NE.0, 'allocation of lim_pc_prfl_qmin')
  coef_pccomp2%lim_pc_prfl_qmin = coef_pccomp1%lim_pc_prfl_qmin

  ALLOCATE(coef_pccomp2%lim_pc_prfl_qmax(coef_pccomp2%fmv_pc_nlev), stat=err)
  THROWM(err.NE.0, 'allocation of lim_pc_prfl_qmax')
  coef_pccomp2%lim_pc_prfl_qmax = coef_pccomp1%lim_pc_prfl_qmax

  ALLOCATE(coef_pccomp2%lim_pc_prfl_ozmin(coef_pccomp2%fmv_pc_nlev), stat=err)
  THROWM(err.NE.0, 'allocation of lim_pc_prfl_ozmin')
  coef_pccomp2%lim_pc_prfl_ozmin = coef_pccomp1%lim_pc_prfl_ozmin

  ALLOCATE(coef_pccomp2%lim_pc_prfl_ozmax(coef_pccomp2%fmv_pc_nlev), stat=err)
  THROWM(err.NE.0, 'allocation of lim_pc_prfl_ozmax')
  coef_pccomp2%lim_pc_prfl_ozmax = coef_pccomp1%lim_pc_prfl_ozmax

  IF (coef_pccomp1%fmv_pc_comp_pc >= 5) THEN
    ALLOCATE(coef_pccomp2%lim_pc_prfl_gasmin(coef_pccomp2%fmv_pc_nlev,coef_pccomp2%fmv_pc_gas_lim), &
             coef_pccomp2%lim_pc_prfl_gasmax(coef_pccomp2%fmv_pc_nlev,coef_pccomp2%fmv_pc_gas_lim), stat=err)
    THROWM(err.NE.0, 'allocation of lim_pc_prfl_gasmin/max')
    coef_pccomp2%lim_pc_prfl_gasmin = coef_pccomp1%lim_pc_prfl_gasmin
    coef_pccomp2%lim_pc_prfl_gasmax = coef_pccomp1%lim_pc_prfl_gasmax
  ENDIF

  IF (coef_pccomp1%fmv_pc_aer /= 0) THEN
    ALLOCATE(coef_pccomp2%lim_pc_prfl_aermin(coef_pccomp1%fmv_pc_naer_types,coef_pccomp2%fmv_pc_nlev-1), &
             coef_pccomp2%lim_pc_prfl_aermax(coef_pccomp1%fmv_pc_naer_types,coef_pccomp2%fmv_pc_nlev-1), stat=err)
    THROWM(err.NE.0, 'allocation of lim_pc_prfl_aermin/max')
    coef_pccomp2%lim_pc_prfl_aermin = coef_pccomp1%lim_pc_prfl_aermin
    coef_pccomp2%lim_pc_prfl_aermax = coef_pccomp1%lim_pc_prfl_aermax
  ENDIF


  ! Instrument noise

  ALLOCATE(coef_pccomp2%noise(coef_pccomp2%fmv_pc_nchn_noise), stat=err)
  THROWM(err.NE.0, 'allocation of noise')
  IF (all_channels) THEN
    coef_pccomp2%noise = coef_pccomp1%noise
  ELSE
    coef_pccomp2%noise = coef_pccomp1%noise(channels)
  ENDIF

  ALLOCATE(coef_pccomp2%noise_in(coef_pccomp2%fmv_pc_nchn), stat=err)
  THROWM(err.NE.0, 'allocation of noise_in')

  ALLOCATE(coef_pccomp2%ff_ori_chn_in(coef_pccomp2%fmv_pc_nchn), stat=err)
  THROWM(err.NE.0, 'allocation of ff_ori_chn_in')

  ALLOCATE(coef_pccomp2%ff_cwn_in(coef_pccomp2%fmv_pc_nchn), stat=err)
  THROWM(err.NE.0, 'allocation of ff_cwn_in')

  ALLOCATE(coef_pccomp2%ff_bco_in(coef_pccomp2%fmv_pc_nchn), stat=err)
  THROWM(err.NE.0, 'allocation of ff_bco_in')

  ALLOCATE(coef_pccomp2%ff_bcs_in(coef_pccomp2%fmv_pc_nchn), stat=err)
  THROWM(err.NE.0, 'allocation of ff_bcs_in')

  IF (all_channels_rec) THEN
    coef_pccomp2%noise_in      = coef_pccomp1%noise_in
    coef_pccomp2%ff_ori_chn_in = coef_pccomp1%ff_ori_chn_in
    coef_pccomp2%ff_cwn_in     = coef_pccomp1%ff_cwn_in
    coef_pccomp2%ff_bco_in     = coef_pccomp1%ff_bco_in
    coef_pccomp2%ff_bcs_in     = coef_pccomp1%ff_bcs_in
  ELSE
    coef_pccomp2%noise_in      = coef_pccomp1%noise_in(channels_rec)
    coef_pccomp2%ff_ori_chn_in = coef_pccomp1%ff_ori_chn_in(channels_rec)
    coef_pccomp2%ff_cwn_in     = coef_pccomp1%ff_cwn_in(channels_rec)
    coef_pccomp2%ff_bco_in     = coef_pccomp1%ff_bco_in(channels_rec)
    coef_pccomp2%ff_bcs_in     = coef_pccomp1%ff_bcs_in(channels_rec)
  ENDIF

CATCH
END SUBROUTINE rttov_channel_extract_pccoef
