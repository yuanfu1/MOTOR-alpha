! Description:
!> @file
!!   Write an ASCII PC coefficient file.
!
!> @brief
!!   Write an ASCII PC coefficient file.
!!
!! @details
!!   The file unit must be open when this subroutine is called.
!!
!! @param[out]    err             status on exit
!! @param[in]     coef_pccomp     PC-RTTOV coefficient structure
!! @param[in]     file_id         logical unit for output pccoef file
!! @param[in]     verbose         flag to switch verbose output on/off (default TRUE), optional
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
SUBROUTINE rttov_write_ascii_pccoef(err, coef_pccomp, file_id, verbose)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef_pccomp
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE

  INTEGER(KIND=jpim)     , INTENT(OUT)          :: err
  TYPE(rttov_coef_pccomp), INTENT(IN)           :: coef_pccomp
  INTEGER(KIND=jpim)     , INTENT(IN)           :: file_id
  LOGICAL(KIND=jplm)     , INTENT(IN), OPTIONAL :: verbose
!INTF_END
#include "rttov_errorreport.interface"

  CHARACTER(LEN=32)   :: section
  CHARACTER(LEN = 80) :: errMessage
  INTEGER(KIND=jpim)  :: m, n, i, j
  LOGICAL(KIND=jplm)  :: lverbose
!- End of header --------------------------------------------------------

  TRY
  IF (PRESENT(verbose)) THEN
    lverbose = verbose
  ELSE
    lverbose = .TRUE._jplm
  END IF

  !ASCII file
  IF (lverbose) THEN
    WRITE (errMessage, '( "write coefficient to file_id ", i2, " in ASCII format")')file_id
    INFO(errMessage)
  END IF

  section = 'PRINCOMP_PREDICTORS'
  WRITE (file_id, "(a)", iostat = err) TRIM(section)
  THROW(err.NE.0)

  WRITE (file_id,  '(i5,1x,a)', iostat=err)coef_pccomp%fmv_pc_comp_pc, '!Principal components coefficient file version number'
  THROWM(err.NE.0, 'io status while writing section '//section)

  WRITE (file_id,  '(i5,1x,a)', iostat=err)coef_pccomp%fmv_pc_cld, '!Compatibility with cloudy computations'
  THROWM(err.NE.0, 'io status while writing section '//section)

  IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN
    WRITE (file_id,  '(i5,1x,a)', iostat=err)coef_pccomp%fmv_pc_aer, '!Compatibility with aerosol computations'
    THROWM(err.NE.0, 'io status while writing section '//section)

    WRITE (file_id,  '(i5,1x,a)', iostat=err)coef_pccomp%fmv_pc_nlte, '!Compatibility with NLTE computations'
    THROWM(err.NE.0, 'io status while writing section '//section)
  ENDIF

  WRITE (file_id,  '(i5,1x,a)', iostat=err)coef_pccomp%fmv_pc_bands, '!Number of bands'
  THROWM(err.NE.0, 'io status while writing section '//section)

  WRITE (file_id,  '(i5,1x,a)', iostat=err)coef_pccomp%fmv_pc_msets, '!Maximum number of regression sets'
  THROWM(err.NE.0, 'io status while writing section '//section)

  WRITE (file_id,  '(i5,1x,a)', iostat=err)coef_pccomp%fmv_pc_mnum, ' !Maximum number of eigenvectors'
  THROWM(err.NE.0, 'io status while writing section '//section)

  WRITE (file_id,  '(i5,1x,a)', iostat=err)coef_pccomp%fmv_pc_mchn, '!Maximum number of channels'
  THROWM(err.NE.0, 'io status while writing section '//section)

! loop on bands
  DO m = 1, coef_pccomp%fmv_pc_bands
    WRITE (file_id,  '(i5,1x,a,i5)', iostat=err)coef_pccomp%fmv_pc_sets(m), '!Number of regression sets in band', m
    THROWM(err.NE.0, 'io status while writing section '//section)

! loop on predictor sets
    DO n = 1, coef_pccomp%fmv_pc_sets(m)

      WRITE (file_id,  '(i5,1x,a,i5)', iostat=err)coef_pccomp%pcreg(m,n)%fmv_pc_npred, '!Number of predictors in set', n
      THROWM(err.NE.0, 'io status while writing section '//section)

      WRITE (file_id, '(a,i3)', iostat = err) '!Channel indices for set', n
      THROWM(err.NE.0, 'io status while writing section '//section)

      WRITE (file_id,  '(10i6)', iostat=err) &
          (coef_pccomp%pcreg(m,n)%predictindex(i), i = 1, coef_pccomp%pcreg(m,n)%fmv_pc_npred)
      THROWM(err.NE.0, 'io status while writing section '//section)
    ENDDO
  ENDDO

  section = 'PRINCOMP_EIGENVECTORS'
  WRITE (file_id, "(a)", iostat = err) TRIM(section)
  THROW(err.NE.0)

! Take care of the user list of channels
! file_channels store the number of channels in the file
! coef % fmv_chn is the number of channels that the user requests

  DO m = 1, coef_pccomp%fmv_pc_bands
    DO n = 1, coef_pccomp%fmv_pc_mnum
      WRITE (file_id,  '(5e23.14)', iostat=err) &
          (coef_pccomp%eigen(m)%eigenvectors(i, n), i = 1, coef_pccomp%fmv_pc_nchn)
      THROWM(err.NE.0, 'io status while writing section '//section)
   ENDDO
  ENDDO

  section = 'PRINCOMP_COEFFICIENTS'
  WRITE (file_id, "(a)", iostat = err) TRIM(section)
  THROW(err.NE.0)

  DO m = 1, coef_pccomp%fmv_pc_bands
    DO n = 1, coef_pccomp%fmv_pc_sets(m)
      WRITE (file_id, '(a,1x,i5)', iostat = err) '! Predictor set ', n
      THROW(err.NE.0)
      DO j = 1, coef_pccomp%fmv_pc_mnum
        WRITE (file_id, '(5e23.14)', iostat=err) &
            (coef_pccomp%pcreg(m,n)%coefficients(i, j), i = 1, coef_pccomp%pcreg(m,n)%fmv_pc_npred)
        THROWM(err.NE.0, 'io status while writing section '//section)
      ENDDO
    ENDDO
  ENDDO

  section = 'EMISSIVITY_COEFFICIENTS'
  WRITE (file_id, "(a)", iostat = err) TRIM(section)
  THROW(err.NE.0)

  WRITE (file_id, "(a)", iostat = err) '! Emissivity coefficients for RTIASI routine'
  THROW(err.NE.0)

  WRITE (file_id, '(i5,1x,a)', iostat=err)coef_pccomp%fmv_pc_nche, '!Number of channels for which coefficients are available'
  THROWM(err.NE.0, 'io status while writing section '//section)

! Take care of the user list of channels
! file_channels store the number of channels in the file
! coef % fmv_chn is the number of channels that the user requests

  DO i = 1, coef_pccomp%fmv_pc_nche
    WRITE (file_id,  '(i5,5e17.9,/,5x,4e17.9)', iostat=err) &
        coef_pccomp%emiss_chn(i), coef_pccomp%emiss_c1(i), coef_pccomp%emiss_c2(i),                         &
        coef_pccomp%emiss_c3(i), coef_pccomp%emiss_c4(i), coef_pccomp%emiss_c5(i), coef_pccomp%emiss_c6(i), &
        coef_pccomp%emiss_c7(i), coef_pccomp%emiss_c8(i), coef_pccomp%emiss_c9(i)
    THROWM(err.NE.0, 'io status while writing section '//section)
  ENDDO

  IF (coef_pccomp%fmv_pc_comp_pc < 5) THEN
    section = 'PC_REFERENCE_PROFILE'
    WRITE (file_id, "(a)", iostat = err) TRIM(section)
    THROW(err.NE.0)

    WRITE (file_id,  '(i5,1x,a)', iostat=err)coef_pccomp%fmv_pc_gas, &
        '! Number of gases for which reference profiles must be used (CO2, N2O, CO, CH4)'
    THROWM(err.NE.0, 'io status while writing section '//section)

    WRITE (file_id,  '(i5,1x,a)', iostat=err)coef_pccomp%fmv_pc_nlev, '! Number of levels'
    THROWM(err.NE.0, 'io status while writing section '//section)

    DO n = 1, coef_pccomp%fmv_pc_gas
      WRITE (file_id,  '("! ",i5)', iostat=err) n
      THROWM(err.NE.0, 'io status while writing section '//section)
      DO i = 1, coef_pccomp%fmv_pc_nlev
        WRITE (file_id,  '(2g26.16)', iostat=err)coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%ref_pc_prfl_mr(i,n)
        THROWM(err.NE.0, 'io status while writing section '//section)
      ENDDO
    ENDDO
  ENDIF

  section = 'PC_PROFILE_LIMITS'
  WRITE (file_id, "(a)", iostat = err) TRIM(section)
  THROW(err.NE.0)

  IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN
    WRITE (file_id,  '(i5,1x,a)', iostat=err)coef_pccomp%fmv_pc_gas_lim, &
        '! Number of additional gases (beyond water vapour and ozone) for which limit profiles are stored (CO2, N2O, CO, CH4)'
    THROWM(err.NE.0, 'io status while writing section '//section)

    WRITE (file_id,  '(i5,1x,a)', iostat=err)coef_pccomp%fmv_pc_nlev, '! Number of levels'
    THROWM(err.NE.0, 'io status while writing section '//section)
  ENDIF

  WRITE (file_id, "(a)", iostat = err) '! Ref.pressure (hPa)'
  THROW(err.NE.0)
  WRITE (file_id, "(a)", iostat = err) '! Temperature Max and Min [K]'
  THROW(err.NE.0)

  DO i = 1, coef_pccomp%fmv_pc_nlev
    WRITE (file_id,  '(3g26.16)', iostat=err)     &
        coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_tmin(i), coef_pccomp%lim_pc_prfl_tmax(i)
    THROWM(err.NE.0, 'io status while writing section '//section)
  ENDDO

  WRITE (file_id, "(a)", iostat = err) '! Ref.pressure (hPa)'
  THROW(err.NE.0)
  WRITE (file_id, "(a)", iostat = err) '! H2O Max and Min [ppmv]'
  THROW(err.NE.0)

  DO i = 1, coef_pccomp%fmv_pc_nlev
    WRITE (file_id, '(3g26.16)', iostat=err)     &
        coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_qmin(i), coef_pccomp%lim_pc_prfl_qmax(i)
    THROWM(err.NE.0, 'io status while writing section '//section)
  ENDDO

  WRITE (file_id, "(a)", iostat = err) '! Ref.pressure (hPa)'
  THROW(err.NE.0)
  WRITE (file_id, "(a)", iostat = err) '! O3 Max and Min [ppmv]'
  THROW(err.NE.0)

  DO i = 1, coef_pccomp%fmv_pc_nlev
    WRITE (file_id,  '(3g26.16)', iostat=err)     &
        coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_ozmin(i), coef_pccomp%lim_pc_prfl_ozmax(i)
    THROWM(err.NE.0, 'io status while writing section '//section)
  ENDDO

  IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN
    DO n = 1, coef_pccomp%fmv_pc_gas_lim
      WRITE (file_id,  '("! ",i5)', iostat=err) n
      THROWM(err.NE.0, 'io status while writing section '//section)
      DO i = 1, coef_pccomp%fmv_pc_nlev
        WRITE (file_id,  '(3g26.16)', iostat=err)coef_pccomp%ref_pc_prfl_p(i), &
                                                 coef_pccomp%lim_pc_prfl_gasmin(i,n), &
                                                 coef_pccomp%lim_pc_prfl_gasmax(i,n)
        THROWM(err.NE.0, 'io status while writing section '//section)
      ENDDO
    ENDDO
  ENDIF

  WRITE (file_id,  '(a)', iostat = err) '! -------------------------------------------------------'
  THROW(err.NE.0)
  WRITE (file_id,  '(a)', iostat = err) '!  PC_PROFILE_LIMITS Surface pressure'
  THROW(err.NE.0)
  WRITE (file_id,  '(a)', iostat = err) '!  Max and Min [hPa]'
  THROW(err.NE.0)

  WRITE (file_id,  '(2g26.16)', iostat=err)coef_pccomp%lim_pc_prfl_pmin, coef_pccomp%lim_pc_prfl_pmax
  THROWM(err.NE.0, 'io status while writing section '//section)

  WRITE (file_id,  '(a)', iostat = err) '! -------------------------------------------------------'
  THROW(err.NE.0)
  WRITE (file_id,  '(a)', iostat = err) '!  PC_PROFILE_LIMITS Surface temperature'
  THROW(err.NE.0)
  WRITE (file_id,  '(a)', iostat = err) '!  Max and Min [K]'
  THROW(err.NE.0)

  WRITE (file_id,  '(2g26.16)', iostat=err)coef_pccomp%lim_pc_prfl_tsmin, coef_pccomp%lim_pc_prfl_tsmax
  THROWM(err.NE.0, 'io status while writing section '//section)

  WRITE (file_id,  '(a)', iostat = err) '! -------------------------------------------------------'
  THROW(err.NE.0)
  WRITE (file_id,  '(a)', iostat = err) '!  PC_PROFILE_LIMITS Skin temperature'
  THROW(err.NE.0)
  WRITE (file_id,  '(a)', iostat = err) '!  Max and Min [K]'
  THROW(err.NE.0)

  WRITE (file_id,  '(2g26.16)', iostat=err)coef_pccomp%lim_pc_prfl_skmin, coef_pccomp%lim_pc_prfl_skmax
  THROWM(err.NE.0, 'io status while writing section '//section)

  WRITE (file_id,  '(a)', iostat = err) '! -------------------------------------------------------'
  THROW(err.NE.0)
  WRITE (file_id,  '(a)', iostat = err) '!  PC_PROFILE_LIMITS 10m wind speed'
  THROW(err.NE.0)
  WRITE (file_id,  '(a)', iostat = err) '!  Max and Min [m/s]'
  THROW(err.NE.0)

  WRITE (file_id,  '(2g26.16)', iostat=err)coef_pccomp%lim_pc_prfl_wsmin, coef_pccomp%lim_pc_prfl_wsmax
  THROWM(err.NE.0, 'io status while writing section '//section)

  IF (coef_pccomp%fmv_pc_aer /= 0) THEN
    WRITE (file_id,  '(a)', iostat = err) '! -------------------------------------------------------'
    THROW(err.NE.0)
    WRITE (file_id,  '(a)', iostat = err) '!  PC_PROFILE_LIMITS Aerosols'
    THROW(err.NE.0)
    WRITE (file_id,  '(i5,1x,a)', iostat=err)coef_pccomp%fmv_pc_naer_types, '!Number of aerosol types'
    THROWM(err.NE.0, 'io status while writing section '//section)

    WRITE (file_id,  '(a)', iostat = err) '!  Aerosol Min [cm^-3]'
    THROW(err.NE.0)
    DO i = 1, coef_pccomp%fmv_pc_nlev-1
      WRITE (file_id, '(10f14.6)', iostat=err)coef_pccomp%lim_pc_prfl_aermin(:,i)
      THROWM(err.NE.0, 'io status while writing section '//section)
    ENDDO
    WRITE (file_id,  '(a)', iostat = err) '!  Aerosol Max [cm^-3]'
    THROW(err.NE.0)
    DO i = 1, coef_pccomp%fmv_pc_nlev-1
      WRITE (file_id, '(10f14.6)', iostat=err)coef_pccomp%lim_pc_prfl_aermax(:,i)
      THROWM(err.NE.0, 'io status while writing section '//section)
    ENDDO
  ENDIF
  WRITE (file_id,  '(a)', iostat = err) '! -------------------------------------------------------'
  THROW(err.NE.0)

  section = 'INSTRUMENT_NOISE'
  WRITE (file_id, "(a)", iostat = err) TRIM(section)
  THROW(err.NE.0)

  DO n = 1, coef_pccomp%fmv_pc_nchn
    WRITE (file_id,  '(i8,4g26.16)', iostat=err)coef_pccomp%ff_ori_chn_in(n), coef_pccomp%ff_cwn_in(n), &
        coef_pccomp%ff_bco_in(n), coef_pccomp%ff_bcs_in(n), coef_pccomp%noise_in(n)
    THROWM(err.NE.0, 'io status while writing section '//section)
  ENDDO

  WRITE(file_id, '(a)', iostat = err) 'END'
  THROW(err.NE.0)

  IF (lverbose) INFO("end of write coefficient")

  CATCH
END SUBROUTINE
