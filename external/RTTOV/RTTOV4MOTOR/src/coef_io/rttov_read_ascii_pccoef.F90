! Description:
!> @file
!!   Read an ASCII PC-RTTOV coefficient file, optionally extracting a subset
!!   of channels.
!
!> @brief
!!   Read an ASCII PC-RTTOV coefficient file, optionally extracting a subset
!!   of channels.
!!
!! @details
!!   The file unit must be open when this subroutine is called.
!!
!!   Note that after reading a subset of channels RTTOV will identify them by
!!   indexes 1...SIZE(channels), not by the original channel numbers.
!!
!!   Any channel subset extracted MUST be a superset of the set of predictor
!!   channels being used in the PC-RTTOV simulations.
!!
!! @param[out]    err             status on exit
!! @param[in]     opts            RTTOV options structure
!! @param[in]     coef            RTTOV optical depth coefficient structure
!! @param[in,out] coef_pccomp     PC-RTTOV coefficient structure
!! @param[in]     file_id         logical unit for input pccoef file
!! @param[in]     channels        list of PC-RTTOV predictor channels to read, optional
!! @param[in]     channels_rec    list of reconstructed radiance channels to read, optional
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
SUBROUTINE rttov_read_ascii_pccoef( &
              err,           &
              opts,          &
              coef,          &
              coef_pccomp,   &
              file_id,       &
              channels,      &
              channels_rec)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef, rttov_coef_pccomp, rttov_options
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_hi, lensection, &
                          pccoef_version_compatible_min, &
                          pccoef_version_compatible_max
  USE parkind1, ONLY : jprb, jplm
!INTF_ON
  IMPLICIT NONE

  INTEGER(KIND=jpim),      INTENT(OUT)          :: err
  TYPE(rttov_options),     INTENT(IN)           :: opts
  TYPE(rttov_coef),        INTENT(IN)           :: coef
  TYPE(rttov_coef_pccomp), INTENT(INOUT)        :: coef_pccomp
  INTEGER(KIND=jpim),      INTENT(IN)           :: file_id
  INTEGER(KIND=jpim),      INTENT(IN), OPTIONAL :: channels(:)
  INTEGER(KIND=jpim),      INTENT(IN), OPTIONAL :: channels_rec(:)
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_skipcommentline.interface"
#include "rttov_findnextsection.interface"
#include "rttov_nullify_coef_pccomp.interface"

  INTEGER(KIND=jpim) :: file_channels
  INTEGER(KIND=jpim) :: file_channels_rec
  LOGICAL(KIND=jplm) :: all_channels
  LOGICAL(KIND=jplm) :: all_channels_rec
  INTEGER(KIND=jpim) :: io_status
  INTEGER(KIND=jpim) :: i, j, m, n

  REAL   (KIND=jprb), POINTER :: values0      (:)
  REAL   (KIND=jprb), POINTER :: values1      (:)
  REAL   (KIND=jprb), POINTER :: values2      (:)
  REAL   (KIND=jprb), POINTER :: values3      (:)
  REAL   (KIND=jprb), POINTER :: values4      (:)
  REAL   (KIND=jprb), POINTER :: values5      (:)
  REAL   (KIND=jprb), POINTER :: values6      (:)
  REAL   (KIND=jprb), POINTER :: values7      (:)
  REAL   (KIND=jprb), POINTER :: values8      (:)
  INTEGER(KIND=jpim), POINTER :: ivalues0     (:)
  REAL   (KIND=jprb), POINTER :: eigenarray   (:,:)
  REAL   (KIND=jprb), POINTER :: noisearray   (:)
  REAL   (KIND=jprb), POINTER :: pcbcoarray   (:)
  REAL   (KIND=jprb), POINTER :: pcbcsarray   (:)
  REAL   (KIND=jprb), POINTER :: pccwnarray   (:)
  INTEGER(KIND=jpim), POINTER :: pcchnarray   (:)
  CHARACTER(LEN=lensection)   :: section
!- End of header --------------------------------------------------------
  TRY

  IF (coef%id_sensor /= sensor_id_hi) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0, "PC-RTTOV only compatible with hi-res sounders")
  ENDIF

  all_channels = .NOT. PRESENT(channels)
  all_channels_rec = .NOT. PRESENT(channels_rec)

  CALL rttov_nullify_coef_pccomp(coef_pccomp)

  readfile : DO
    CALL rttov_findnextsection(file_id, io_status, section)
    IF (io_status < 0) EXIT!end-of-file

    SELECT CASE (TRIM(section))
    CASE ('PRINCOMP_PREDICTORS')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_comp_pc
      THROWM(err.NE.0, 'io status while reading section '//section)

      IF (coef%id_comp_pc /= coef_pccomp%fmv_pc_comp_pc) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0, "Version of PC coef file is incompatible with RTTOV rtcoef file")
      ENDIF

      IF (coef_pccomp%fmv_pc_comp_pc < pccoef_version_compatible_min .OR. &
          coef_pccomp%fmv_pc_comp_pc > pccoef_version_compatible_max) THEN
        err = errorstatus_fatal
        THROWM(err.NE.0, "Version of PC coef file is incompatible with RTTOV library")
      ENDIF

      READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_cld
      THROWM(err.NE.0, 'io status while reading section '//section)

      IF (opts%rt_ir%pc%addpc .AND. opts%rt_ir%addclouds) THEN
        IF (coef_pccomp%fmv_pc_cld == 0)THEN
          err = errorstatus_fatal
          THROWM(err.NE.0, "PC coef file is incompatible with cloudy computations")
        ENDIF
      ENDIF

      IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN
        READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_aer
        THROWM(err.NE.0, 'io status while reading section '//section)

        READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_nlte
        THROWM(err.NE.0, 'io status while reading section '//section)
      ELSE
        coef_pccomp%fmv_pc_aer = 0
        coef_pccomp%fmv_pc_nlte = 0  ! See rttov_init_coef_pccomp
      ENDIF
      IF (opts%rt_ir%pc%addpc .AND. opts%rt_ir%addaerosl) THEN
        IF (coef_pccomp%fmv_pc_aer == 0) THEN
          err = errorstatus_fatal
          THROWM(err.NE.0, "PC coef file is incompatible with aerosol computations")
        ENDIF
      ENDIF
      IF (opts%rt_ir%pc%addpc .AND. opts%rt_ir%do_nlte_correction) THEN
        IF (coef_pccomp%fmv_pc_nlte == 0) THEN
          err = errorstatus_fatal
          THROWM(err.NE.0, "PC coef file is incompatible with NLTE computations")
        ENDIF
      ENDIF

      READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_bands
      THROWM(err.NE.0, 'io status while reading section '//section)

      ALLOCATE (coef_pccomp%fmv_pc_sets(coef_pccomp%fmv_pc_bands), STAT = err)
      THROWM(err.NE.0, "allocation of coef_pccomp%fmv_pc_sets array")

      ALLOCATE (coef_pccomp%eigen(coef_pccomp%fmv_pc_bands), STAT = err)
      THROWM(err.NE.0, "allocation of coef_pccomp%eigen")

      READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_msets
      THROWM(err.NE.0, 'io status while reading section '//section)

      ALLOCATE (coef_pccomp%pcreg(coef_pccomp%fmv_pc_bands,coef_pccomp%fmv_pc_msets), STAT = err)
      THROWM(err.NE.0, "allocation of coef_pccomp%pcreg arrays")

      READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_mnum
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_mchn
      THROWM(err.NE.0, 'io status while reading section '//section)

! loop on bands 
      DO m = 1, coef_pccomp%fmv_pc_bands
        READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_sets(m)
        THROWM(err.NE.0, 'io status while reading section '//section)

! loop on predictor sets
        DO n = 1,  coef_pccomp%fmv_pc_sets(m)
          READ (file_id,  * , iostat=err)coef_pccomp%pcreg(m,n)%fmv_pc_npred
          THROWM(err.NE.0, 'io status while reading section '//section)

          CALL rttov_skipcommentline(file_id, err)
          THROWM(err.NE.0, 'io status while reading section '//section)

          ALLOCATE (coef_pccomp%pcreg(m,n)%predictindex(coef_pccomp%pcreg(m,n)%fmv_pc_npred), STAT = err)
          THROWM(err.NE.0, "allocation of predictindex arrays")

          READ (file_id,  * , iostat=err) &
              (coef_pccomp%pcreg(m,n)%predictindex(i), i = 1, coef_pccomp%pcreg(m,n)%fmv_pc_npred)
          THROWM(err.NE.0, "io status while reading predictindex arrays")
        ENDDO
      ENDDO


    CASE ('PRINCOMP_EIGENVECTORS')
      coef_pccomp%fmv_pc_nchn = coef_pccomp%fmv_pc_mchn

! Take care of the user list of channels
! file_channels store the number of channels in the file
! coef%fmv_chn is the number of channels that the user requests
      file_channels_rec = coef_pccomp%fmv_pc_nchn

      IF (.NOT. all_channels_rec) THEN
        coef_pccomp%fmv_pc_nchn = SIZE(channels_rec)
      ENDIF

      IF (.NOT. all_channels) THEN
        coef_pccomp%fmv_pc_nchn_noise = SIZE(channels)
      ELSE
        coef_pccomp%fmv_pc_nchn_noise = coef_pccomp%fmv_pc_mchn
      ENDIF

      DO m = 1,  coef_pccomp%fmv_pc_bands
        CALL rttov_skipcommentline(file_id, err)
        THROWM(err.NE.0, 'io status while reading section '//section)

        NULLIFY(coef_pccomp%eigen(m)%eigenvectors_t)
        ALLOCATE (coef_pccomp%eigen(m)%eigenvectors(coef_pccomp%fmv_pc_nchn, coef_pccomp%fmv_pc_mnum), STAT = err)
        THROWM(err.NE.0, "allocation of eigenvectors")

        IF (all_channels_rec) THEN
          DO n = 1, coef_pccomp%fmv_pc_mnum
            READ (file_id,  * , iostat=err)(coef_pccomp%eigen(m)%eigenvectors(i, n), i = 1, coef_pccomp%fmv_pc_nchn)
            THROWM(err.NE.0, 'io status while reading section '//section)
          ENDDO
        ELSE
          ALLOCATE (eigenarray(file_channels_rec, coef_pccomp%fmv_pc_mnum), STAT = err)
          THROWM(err.NE.0, "allocation of eigenarray")

          DO n = 1, coef_pccomp%fmv_pc_mnum
            READ (file_id,  * , iostat=err)(eigenarray(i, n), i = 1, file_channels_rec)
            THROWM(err.NE.0, 'io status while reading section '//section)
          ENDDO

          coef_pccomp%eigen(m)%eigenvectors(:,:) = eigenarray(channels_rec(:), :)
          DEALLOCATE (eigenarray, STAT = err)
          THROWM(err.NE.0, "deallocation of eigenarray")
        ENDIF
      ENDDO


    CASE ('PRINCOMP_COEFFICIENTS')

      DO m=1,coef_pccomp%fmv_pc_bands
        DO n = 1, coef_pccomp%fmv_pc_sets(m)
          NULLIFY(coef_pccomp%pcreg(m,n)%coefficients_t)
          ALLOCATE(coef_pccomp%pcreg(m,n)%coefficients(coef_pccomp%pcreg(m,n)%fmv_pc_npred, coef_pccomp%fmv_pc_mnum), &
                   STAT = err)
          THROWM(err.NE.0, "allocation of pcreg(n)%coefficients")

          CALL rttov_skipcommentline(file_id, err)
          THROWM(err.NE.0, 'io status while reading section '//section)

          DO j = 1, coef_pccomp%fmv_pc_mnum
            READ (file_id,  * , iostat=err)     &
                (coef_pccomp%pcreg(m,n)%coefficients(i, j), i = 1, coef_pccomp%pcreg(m,n)%fmv_pc_npred)
            THROWM(err.NE.0, 'io status while reading section '//section)
          ENDDO
        ENDDO
      ENDDO


    CASE ('EMISSIVITY_COEFFICIENTS')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_nche
      THROWM(err.NE.0, 'io status while reading section '//section)

! Take care of the user list of channels
! file_channels store the number of channels in the file
! coef%fmv_chn is the number of channels that the user requests
      file_channels = coef_pccomp%fmv_pc_nche

      IF (.NOT. all_channels) THEN
        coef_pccomp%fmv_pc_nche = SIZE(channels)
      ENDIF

      ALLOCATE (coef_pccomp%emiss_chn(coef_pccomp%fmv_pc_nche), &
                coef_pccomp%emiss_c1(coef_pccomp%fmv_pc_nche),  &
                coef_pccomp%emiss_c2(coef_pccomp%fmv_pc_nche),  &
                coef_pccomp%emiss_c3(coef_pccomp%fmv_pc_nche),  &
                coef_pccomp%emiss_c4(coef_pccomp%fmv_pc_nche),  &
                coef_pccomp%emiss_c5(coef_pccomp%fmv_pc_nche),  &
                coef_pccomp%emiss_c6(coef_pccomp%fmv_pc_nche),  &
                coef_pccomp%emiss_c7(coef_pccomp%fmv_pc_nche),  &
                coef_pccomp%emiss_c8(coef_pccomp%fmv_pc_nche),  &
                coef_pccomp%emiss_c9(coef_pccomp%fmv_pc_nche), STAT = err)
      THROWM(err.NE.0, "allocation of coef_pccomp%emiss_*")

      IF (all_channels) THEN
        DO i = 1, coef_pccomp%fmv_pc_nche
          READ (file_id,  * , iostat=err)coef_pccomp%emiss_chn(i), coef_pccomp%emiss_c1(i), coef_pccomp%emiss_c2(i), &
              coef_pccomp%emiss_c3(i), coef_pccomp%emiss_c4(i), coef_pccomp%emiss_c5(i), coef_pccomp%emiss_c6(i),    &
              coef_pccomp%emiss_c7(i), coef_pccomp%emiss_c8(i), coef_pccomp%emiss_c9(i)
          THROWM(err.NE.0, 'io status while reading section '//section)
        ENDDO
      ELSE
      ALLOCATE (ivalues0(file_channels), &
                values0(file_channels),  &
                values1(file_channels),  &
                values2(file_channels),  &
                values3(file_channels),  &
                values4(file_channels),  &
                values5(file_channels),  &
                values6(file_channels),  &
                values7(file_channels),  &
                values8(file_channels), STAT = err)
        THROWM(err.NE.0, "allocation of values*")

        DO i = 1, file_channels
          READ (file_id,  * , iostat=err)ivalues0(i), values0(i), values1(i), values2(i), values3(i), values4(i), &
              values5(i), values6(i), values7(i), values8(i)
          THROWM(err.ne.0,"io status while reading section "//section)
        ENDDO

        coef_pccomp%emiss_chn(:) = ivalues0(channels(:))
        coef_pccomp%emiss_c1(:)  = values0(channels(:))
        coef_pccomp%emiss_c2(:)  = values1(channels(:))
        coef_pccomp%emiss_c3(:)  = values2(channels(:))
        coef_pccomp%emiss_c4(:)  = values3(channels(:))
        coef_pccomp%emiss_c5(:)  = values4(channels(:))
        coef_pccomp%emiss_c6(:)  = values5(channels(:))
        coef_pccomp%emiss_c7(:)  = values6(channels(:))
        coef_pccomp%emiss_c8(:)  = values7(channels(:))
        coef_pccomp%emiss_c9(:)  = values8(channels(:))
        DEALLOCATE (ivalues0, values0, values1, &
                    values2,  values3, values4, &
                    values5,  values6, values7, &
                    values8, STAT = err)
        THROWM(err.NE.0, "deallocation of values*")
      ENDIF


    CASE ('PC_REFERENCE_PROFILE')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_gas
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_nlev
      THROWM(err.NE.0, 'io status while reading section '//section)

      ALLOCATE (coef_pccomp%ref_pc_prfl_p(coef_pccomp%fmv_pc_nlev), STAT = err)
      THROWM(err.NE.0, "allocation of pccomp%ref_pc_prfl_p")

      ALLOCATE (coef_pccomp%ref_pc_prfl_mr(coef_pccomp%fmv_pc_nlev, coef_pccomp%fmv_pc_gas), STAT = err)
      THROWM(err.NE.0, "allocation of pccomp%ref_pc_prfl_mr")

      DO n = 1, coef_pccomp%fmv_pc_gas
        CALL rttov_skipcommentline(file_id, err)
        THROWM(err.NE.0, 'io status while reading section '//section)

        DO i = 1, coef_pccomp%fmv_pc_nlev
          READ (file_id,  * , iostat=err)coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%ref_pc_prfl_mr(i,n)
          THROWM(err.NE.0, 'io status while reading section '//section)
        ENDDO
      ENDDO


    CASE ('PC_PROFILE_LIMITS')
      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN
        CALL rttov_skipcommentline(file_id, err)
        THROWM(err.NE.0, 'io status while reading section '//section)

        READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_gas_lim
        THROWM(err.NE.0, 'io status while reading section '//section)

        READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_nlev
        THROWM(err.NE.0, 'io status while reading section '//section)

        ALLOCATE (coef_pccomp%lim_pc_prfl_gasmax(coef_pccomp%fmv_pc_nlev, coef_pccomp%fmv_pc_gas_lim), &
                  coef_pccomp%lim_pc_prfl_gasmin(coef_pccomp%fmv_pc_nlev, coef_pccomp%fmv_pc_gas_lim), STAT = err)
        THROWM(err.NE.0, "allocation of PC regression limit arrays")

        ALLOCATE (coef_pccomp%ref_pc_prfl_p(coef_pccomp%fmv_pc_nlev), STAT = err)
        THROWM(err.NE.0, "allocation of pccomp%ref_pc_prfl_p")
      ENDIF

      ALLOCATE (coef_pccomp%lim_pc_prfl_tmax(coef_pccomp%fmv_pc_nlev),  &
                coef_pccomp%lim_pc_prfl_tmin(coef_pccomp%fmv_pc_nlev),  &
                coef_pccomp%lim_pc_prfl_qmax(coef_pccomp%fmv_pc_nlev),  &
                coef_pccomp%lim_pc_prfl_qmin(coef_pccomp%fmv_pc_nlev),  &
                coef_pccomp%lim_pc_prfl_ozmax(coef_pccomp%fmv_pc_nlev), &
                coef_pccomp%lim_pc_prfl_ozmin(coef_pccomp%fmv_pc_nlev), STAT = err)
      THROWM(err.NE.0, "allocation of PC regression limit arrays")

      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      DO i = 1, coef_pccomp%fmv_pc_nlev
        READ (file_id,  * , iostat=err)     &
            coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_tmin(i), coef_pccomp%lim_pc_prfl_tmax(i)
        THROWM(err.NE.0, 'io status while reading section '//section)
      ENDDO

      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      DO i = 1, coef_pccomp%fmv_pc_nlev
        READ (file_id,  * , iostat=err)     &
            coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_qmin(i), coef_pccomp%lim_pc_prfl_qmax(i)
        THROWM(err.NE.0, 'io status while reading section '//section)
      ENDDO

      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      DO i = 1, coef_pccomp%fmv_pc_nlev
        READ (file_id,  * , iostat=err)     &
            coef_pccomp%ref_pc_prfl_p(i), coef_pccomp%lim_pc_prfl_ozmin(i), coef_pccomp%lim_pc_prfl_ozmax(i)
        THROWM(err.NE.0, 'io status while reading section '//section)
      ENDDO

      IF (coef_pccomp%fmv_pc_comp_pc >= 5) THEN
        DO n = 1, coef_pccomp%fmv_pc_gas_lim
          CALL rttov_skipcommentline(file_id, err)
          THROWM(err.NE.0, 'io status while reading section '//section)

          DO i = 1, coef_pccomp%fmv_pc_nlev
            READ (file_id,  * , iostat=err)coef_pccomp%ref_pc_prfl_p(i), &
                                           coef_pccomp%lim_pc_prfl_gasmin(i,n), &
                                           coef_pccomp%lim_pc_prfl_gasmax(i,n)
            THROWM(err.NE.0, 'io status while reading section '//section)
          ENDDO
        ENDDO
      ENDIF

      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)coef_pccomp%lim_pc_prfl_pmin, coef_pccomp%lim_pc_prfl_pmax
      THROWM(err.NE.0, 'io status while reading section '//section)

      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)coef_pccomp%lim_pc_prfl_tsmin, coef_pccomp%lim_pc_prfl_tsmax
      THROWM(err.NE.0, 'io status while reading section '//section)

      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)coef_pccomp%lim_pc_prfl_skmin, coef_pccomp%lim_pc_prfl_skmax
      THROWM(err.NE.0, 'io status while reading section '//section)

      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      READ (file_id,  * , iostat=err)coef_pccomp%lim_pc_prfl_wsmin, coef_pccomp%lim_pc_prfl_wsmax
      THROWM(err.NE.0, 'io status while reading section '//section)

      IF (coef_pccomp%fmv_pc_aer /= 0) THEN
        CALL rttov_skipcommentline(file_id, err)
        THROWM(err.NE.0, 'io status while reading section '//section)

        READ (file_id,  * , iostat=err)coef_pccomp%fmv_pc_naer_types
        THROWM(err.NE.0, 'io status while reading section '//section)

        ALLOCATE (coef_pccomp%lim_pc_prfl_aermax(coef_pccomp%fmv_pc_naer_types,coef_pccomp%fmv_pc_nlev-1), &
                  coef_pccomp%lim_pc_prfl_aermin(coef_pccomp%fmv_pc_naer_types,coef_pccomp%fmv_pc_nlev-1), STAT = err)

        CALL rttov_skipcommentline(file_id, err)
        THROWM(err.NE.0, 'io status while reading section '//section)

        DO i = 1, coef_pccomp%fmv_pc_nlev-1
          READ (file_id,  * , iostat=err)coef_pccomp%lim_pc_prfl_aermin(:,i)
          THROWM(err.NE.0, 'io status while reading section '//section)
        ENDDO

        CALL rttov_skipcommentline(file_id, err)
        THROWM(err.NE.0, 'io status while reading section '//section)

        DO i = 1, coef_pccomp%fmv_pc_nlev-1
          READ (file_id,  * , iostat=err)coef_pccomp%lim_pc_prfl_aermax(:,i)
          THROWM(err.NE.0, 'io status while reading section '//section)
        ENDDO
      ENDIF

    CASE ('INSTRUMENT_NOISE')
      ALLOCATE (coef_pccomp%noise_in(coef_pccomp%fmv_pc_nchn), STAT = err)
      THROWM(err.NE.0, "allocation of pccomp%noise_in")

      ALLOCATE (coef_pccomp%ff_ori_chn_in(coef_pccomp%fmv_pc_nchn), STAT = err)
      THROWM(err.NE.0, "allocation of pccomp%ff_ori_chn_in")

      ALLOCATE (coef_pccomp%ff_cwn_in(coef_pccomp%fmv_pc_nchn), STAT = err)
      THROWM(err.NE.0, "allocation of pccomp%ff_cwn_in")

      ALLOCATE (coef_pccomp%ff_bco_in(coef_pccomp%fmv_pc_nchn), STAT = err)
      THROWM(err.NE.0, "allocation of pccomp%ff_bco_in")

      ALLOCATE (coef_pccomp%ff_bcs_in(coef_pccomp%fmv_pc_nchn), STAT = err)
      THROWM(err.NE.0, "allocation of pccomp%ff_bcs_in")

      ALLOCATE (coef_pccomp%noise(coef_pccomp%fmv_pc_nchn_noise), STAT = err)
      THROWM(err.NE.0, "allocation of pccomp%noise")

      CALL rttov_skipcommentline(file_id, err)
      THROWM(err.NE.0, 'io status while reading section '//section)

      IF (all_channels_rec) THEN
        DO n = 1, coef_pccomp%fmv_pc_nchn
          READ (file_id,  * , iostat=err)coef_pccomp%ff_ori_chn_in(n), coef_pccomp%ff_cwn_in(n),      &
              coef_pccomp%ff_bco_in(n), coef_pccomp%ff_bcs_in(n), coef_pccomp%noise_in(n)
          THROWM(err.NE.0, 'io status while reading section '//section)
        ENDDO

        IF (all_channels) THEN
          coef_pccomp%noise = coef_pccomp%noise_in
        ELSE
          coef_pccomp%noise = coef_pccomp%noise_in(coef%ff_ori_chn(:))
        ENDIF
      ELSE
        ALLOCATE (noisearray(file_channels_rec), STAT = err)
        THROWM(err.NE.0, "allocation of noisearray")

        ALLOCATE (pcchnarray(file_channels_rec), STAT = err)
        THROWM(err.NE.0, "allocation of pcchnarray")

        ALLOCATE (pccwnarray(file_channels_rec), STAT = err)
        THROWM(err.NE.0, "allocation of pccwnarray")

        ALLOCATE (pcbcoarray(file_channels_rec), STAT = err)
        THROWM(err.NE.0, "allocation of pcbcoarray")

        ALLOCATE (pcbcsarray(file_channels_rec), STAT = err)
        THROWM(err.NE.0, "allocation of pcbcsarray")

        DO n = 1, file_channels_rec
          READ (file_id,  * , iostat=err)pcchnarray(n), pccwnarray(n), pcbcoarray(n), pcbcsarray(n), noisearray(n)
          THROWM(err.NE.0, 'io status while reading section '//section)
        ENDDO

        coef_pccomp%ff_cwn_in(:)     = pccwnarray(channels_rec(:))
        coef_pccomp%ff_ori_chn_in(:) = pcchnarray(channels_rec(:))
        coef_pccomp%ff_bco_in(:)     = pcbcoarray(channels_rec(:))
        coef_pccomp%ff_bcs_in(:)     = pcbcsarray(channels_rec(:))
        coef_pccomp%noise_in(:)      = noisearray(channels_rec(:))
        coef_pccomp%noise(:)         = noisearray(coef%ff_ori_chn(:))

        DEALLOCATE (noisearray, pccwnarray, pcchnarray, pcbcoarray, pcbcsarray, STAT = err)
        THROWM(err.NE.0, "deallocation of pcbcsarray")
      ENDIF

    CASE ('END')
      RETURN
    CASE DEFAULT
      CYCLE readfile
    END SELECT

  ENDDO readfile

  CATCH
END SUBROUTINE rttov_read_ascii_pccoef
