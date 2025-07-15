! Description:
!> @file
!!   Write a binary PC coefficient file.
!
!> @brief
!!   Write a binary PC coefficient file.
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
SUBROUTINE rttov_write_binary_pccoef(err, coef_pccomp, file_id, verbose)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef_pccomp
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE rttov_const, ONLY : rttov_magic_string, rttov_magic_number
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim)     , INTENT(OUT)          :: err
  TYPE(rttov_coef_pccomp), INTENT(IN)           :: coef_pccomp
  INTEGER(KIND=jpim)     , INTENT(IN)           :: file_id
  LOGICAL(KIND=jplm)     , INTENT(IN), OPTIONAL :: verbose
!INTF_END
#include "rttov_errorreport.interface"
  INTEGER(KIND=jpim)  :: n, m, i, j
  LOGICAL(KIND=jplm)  :: lverbose
  CHARACTER(LEN = 80) :: errMessage
!- End of header --------------------------------------------------------

  TRY
  IF (PRESENT(verbose)) THEN
    lverbose = verbose
  ELSE
    lverbose = .TRUE._jplm
  END IF

! Binary file
  IF (lverbose) THEN
    WRITE (errMessage, '( "write coefficient to file_id ", i2, " in binary format")')file_id
    INFO(errMessage)
  END IF
! Write a string that could be displayed
! Write a real number to be able to check single/double precision
  WRITE (file_id, iostat=err)rttov_magic_string, rttov_magic_number
  THROW(err.NE.0)

! PRINCOMP_PREDICTORS
  WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_comp_pc
  THROW(err.NE.0)

  WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_cld
  THROW(err.NE.0)

  IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN
    WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_aer
    THROW(err.NE.0)

    WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_nlte
    THROW(err.NE.0)
  ENDIF

  WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_bands
  THROW(err.NE.0)

  WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_msets
  THROW(err.NE.0)

  WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_mnum
  THROW(err.NE.0)

  WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_mchn
  THROW(err.NE.0)

  DO m = 1, coef_pccomp%fmv_pc_bands
    WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_sets(m)
    THROW(err.NE.0)

    DO n = 1, coef_pccomp%fmv_pc_sets(m)
      WRITE (file_id, iostat=err)coef_pccomp%pcreg(m,n)%fmv_pc_npred
      THROW(err.NE.0)

      WRITE (file_id, iostat=err) &
         (coef_pccomp%pcreg(m,n)%predictindex(i), i = 1, coef_pccomp%pcreg(m,n)%fmv_pc_npred)
     THROW(err.NE.0)
    ENDDO
  ENDDO

!PRINCOMP_EIGENVECTORS

  DO m = 1, coef_pccomp%fmv_pc_bands
    DO n = 1, coef_pccomp%fmv_pc_mnum
      WRITE (file_id, iostat=err)(coef_pccomp%eigen(m)%eigenvectors(i, n), i = 1, coef_pccomp%fmv_pc_nchn)
      THROW(err.NE.0)
    ENDDO
  ENDDO

!PRINCOMP_COEFFICIENTS
  DO m = 1, coef_pccomp%fmv_pc_bands
    DO n = 1, coef_pccomp%fmv_pc_sets(m)
      DO j = 1, coef_pccomp%fmv_pc_mnum
        WRITE (file_id, iostat=err)     &
            (coef_pccomp%pcreg(m,n)%coefficients(i, j), i = 1, coef_pccomp%pcreg(m,n)%fmv_pc_npred)
        THROW(err.NE.0)
      ENDDO
    ENDDO
  ENDDO

!EMISSIVITY_COEFFICIENTS

  WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_nche
  THROW(err.NE.0)

  DO i = 1, coef_pccomp%fmv_pc_nche
    WRITE (file_id, iostat=err)   &
        coef_pccomp%emiss_chn(i), &
        coef_pccomp%emiss_c1(i),  &
        coef_pccomp%emiss_c2(i),  &
        coef_pccomp%emiss_c3(i),  &
        coef_pccomp%emiss_c4(i),  &
        coef_pccomp%emiss_c5(i),  &
        coef_pccomp%emiss_c6(i),  &
        coef_pccomp%emiss_c7(i),  &
        coef_pccomp%emiss_c8(i),  &
        coef_pccomp%emiss_c9(i)
    THROW(err.NE.0)
  ENDDO


  IF (coef_pccomp%fmv_pc_comp_pc < 5) THEN
!PC_REFERENCE_PROFILE

    WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_gas
    THROW(err.NE.0)

    WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_nlev
    THROW(err.NE.0)

    DO n = 1, coef_pccomp%fmv_pc_gas
      DO i = 1, coef_pccomp%fmv_pc_nlev
        WRITE (file_id, iostat=err)coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%ref_pc_prfl_mr(i, n)
        THROW(err.NE.0)
      ENDDO
    ENDDO
  ENDIF

!PC_PROFILE_LIMITS

  IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN
    WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_gas_lim
    THROW(err.NE.0)

    WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_nlev
    THROW(err.NE.0)
  ENDIF

  DO i = 1, coef_pccomp%fmv_pc_nlev
    WRITE (file_id, iostat=err)          &
        coef_pccomp%ref_pc_prfl_p(i),    &
        coef_pccomp%lim_pc_prfl_tmin(i), &
        coef_pccomp%lim_pc_prfl_tmax(i)
    THROW(err.NE.0)
  ENDDO

  DO i = 1, coef_pccomp%fmv_pc_nlev
    WRITE (file_id, iostat=err)         &
        coef_pccomp%ref_pc_prfl_p(i),   &
        coef_pccomp%lim_pc_prfl_qmin(i),&
        coef_pccomp%lim_pc_prfl_qmax(i)
    THROW(err.NE.0)
  ENDDO

  DO i = 1, coef_pccomp%fmv_pc_nlev
    WRITE (file_id, iostat=err)           &
        coef_pccomp%ref_pc_prfl_p(i),     &
        coef_pccomp%lim_pc_prfl_ozmin(i), &
        coef_pccomp%lim_pc_prfl_ozmax(i)
    THROW(err.NE.0)
  ENDDO

  IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN
    WRITE (file_id, iostat=err)coef_pccomp%lim_pc_prfl_gasmin
    THROW(err.NE.0)

    WRITE (file_id, iostat=err)coef_pccomp%lim_pc_prfl_gasmax
    THROW(err.NE.0)
  ENDIF

  WRITE (file_id, iostat=err)coef_pccomp%lim_pc_prfl_pmin, coef_pccomp%lim_pc_prfl_pmax
  THROW(err.NE.0)

  WRITE (file_id, iostat=err)coef_pccomp%lim_pc_prfl_tsmin, coef_pccomp%lim_pc_prfl_tsmax
  THROW(err.NE.0)

  WRITE (file_id, iostat=err)coef_pccomp%lim_pc_prfl_skmin, coef_pccomp%lim_pc_prfl_skmax
  THROW(err.NE.0)

  WRITE (file_id, iostat=err)coef_pccomp%lim_pc_prfl_wsmin, coef_pccomp%lim_pc_prfl_wsmax
  THROW(err.NE.0)

  IF (coef_pccomp%fmv_pc_aer /= 0) THEN
    WRITE (file_id, iostat=err)coef_pccomp%fmv_pc_naer_types
    THROW(err.NE.0)

    WRITE (file_id, iostat=err)coef_pccomp%lim_pc_prfl_aermin
    THROW(err.NE.0)

    WRITE (file_id, iostat=err)coef_pccomp%lim_pc_prfl_aermax
    THROW(err.NE.0)
  ENDIF

!INSTRUMENT_NOISE
  DO n = 1, coef_pccomp%fmv_pc_nchn
    WRITE (file_id, iostat=err)       &
        coef_pccomp%ff_ori_chn_in(n), &
        coef_pccomp%ff_cwn_in(n),     &
        coef_pccomp%ff_bco_in(n),     &
        coef_pccomp%ff_bcs_in(n),     &
        coef_pccomp%noise_in(n)
    THROW(err.NE.0)
  ENDDO

  IF (lverbose) INFO("end of write coefficient")
  CATCH
END SUBROUTINE 
