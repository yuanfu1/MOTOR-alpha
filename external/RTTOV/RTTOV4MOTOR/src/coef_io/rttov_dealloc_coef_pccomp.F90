! Description:
!> @file
!!   Deallocate a PC coefficients structure.
!!   This should usually be called via rttov_dealloc_coefs.
!
!> @brief
!!   Deallocate a PC coefficients structure.
!!   This should usually be called via rttov_dealloc_coefs.
!!
!! @param[out]     err           status on exit
!! @param[in,out]  coef_pccomp   the PC coefficient structure to deallocate
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
SUBROUTINE rttov_dealloc_coef_pccomp (err, coef_pccomp)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef_pccomp
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  INTEGER(KIND=jpim),      INTENT(OUT)   :: err
  TYPE(rttov_coef_pccomp), INTENT(INOUT) :: coef_pccomp
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_nullify_coef_pccomp.interface"
  INTEGER(KIND=jpim) :: m, n
!- End of header --------------------------------------------------------
  TRY

    IF (ASSOCIATED(coef_pccomp%emiss_chn)) &
        DEALLOCATE(coef_pccomp%emiss_chn, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%emiss_c1)) &
        DEALLOCATE(coef_pccomp%emiss_c1, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%emiss_c2)) &
        DEALLOCATE(coef_pccomp%emiss_c2, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%emiss_c3)) &
        DEALLOCATE(coef_pccomp%emiss_c3, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%emiss_c4)) &
        DEALLOCATE(coef_pccomp%emiss_c4, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%emiss_c5)) &
        DEALLOCATE(coef_pccomp%emiss_c5, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%emiss_c6)) &
        DEALLOCATE(coef_pccomp%emiss_c6, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%emiss_c7)) &
        DEALLOCATE(coef_pccomp%emiss_c7, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%emiss_c8)) &
        DEALLOCATE(coef_pccomp%emiss_c8, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%emiss_c9)) &
        DEALLOCATE(coef_pccomp%emiss_c9, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%ff_cwn_in)) &
        DEALLOCATE(coef_pccomp%ff_cwn_in, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%ff_bco_in)) &
        DEALLOCATE(coef_pccomp%ff_bco_in, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%ff_bcs_in)) &
        DEALLOCATE(coef_pccomp%ff_bcs_in, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%planck1_in)) &
        DEALLOCATE(coef_pccomp%planck1_in, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%planck2_in)) &
        DEALLOCATE(coef_pccomp%planck2_in, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%ff_ori_chn_in)) DEALLOCATE (coef_pccomp%ff_ori_chn_in, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%noise_in)) DEALLOCATE (coef_pccomp%noise_in, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%noise)) DEALLOCATE (coef_pccomp%noise, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%noise_r)) DEALLOCATE (coef_pccomp%noise_r, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%ref_pc_prfl_mr)) DEALLOCATE (coef_pccomp%ref_pc_prfl_mr, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%ref_pc_prfl_p)) DEALLOCATE (coef_pccomp%ref_pc_prfl_p, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%lim_pc_prfl_tmin)) DEALLOCATE (coef_pccomp%lim_pc_prfl_tmin, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%lim_pc_prfl_tmax)) DEALLOCATE (coef_pccomp%lim_pc_prfl_tmax, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%lim_pc_prfl_qmin)) DEALLOCATE (coef_pccomp%lim_pc_prfl_qmin, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%lim_pc_prfl_qmax)) DEALLOCATE (coef_pccomp%lim_pc_prfl_qmax, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%lim_pc_prfl_ozmin)) DEALLOCATE (coef_pccomp%lim_pc_prfl_ozmin, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%lim_pc_prfl_ozmax)) DEALLOCATE (coef_pccomp%lim_pc_prfl_ozmax, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%lim_pc_prfl_gasmin)) DEALLOCATE (coef_pccomp%lim_pc_prfl_gasmin, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%lim_pc_prfl_gasmax)) DEALLOCATE (coef_pccomp%lim_pc_prfl_gasmax, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%lim_pc_prfl_aermin)) DEALLOCATE (coef_pccomp%lim_pc_prfl_aermin, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%lim_pc_prfl_aermax)) DEALLOCATE (coef_pccomp%lim_pc_prfl_aermax, STAT = err)
    THROW(err.NE.0)


    IF (ASSOCIATED(coef_pccomp%co2_pc_ref)) DEALLOCATE (coef_pccomp%co2_pc_ref, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%co_pc_ref)) DEALLOCATE (coef_pccomp%co_pc_ref, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%ch4_pc_ref)) DEALLOCATE (coef_pccomp%ch4_pc_ref, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%n2o_pc_ref)) DEALLOCATE (coef_pccomp%n2o_pc_ref, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%co2_pc_min)) DEALLOCATE (coef_pccomp%co2_pc_min, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%co_pc_min)) DEALLOCATE (coef_pccomp%co_pc_min, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%ch4_pc_min)) DEALLOCATE (coef_pccomp%ch4_pc_min, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%n2o_pc_min)) DEALLOCATE (coef_pccomp%n2o_pc_min, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%co2_pc_max)) DEALLOCATE (coef_pccomp%co2_pc_max, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%co_pc_max)) DEALLOCATE (coef_pccomp%co_pc_max, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%ch4_pc_max)) DEALLOCATE (coef_pccomp%ch4_pc_max, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%n2o_pc_max)) DEALLOCATE (coef_pccomp%n2o_pc_max, STAT = err)
    THROW(err.NE.0)

    IF (ASSOCIATED(coef_pccomp%eigen)) THEN
      DO m = 1, coef_pccomp%fmv_pc_bands
        IF (ASSOCIATED(coef_pccomp%eigen(m)%eigenvectors)) &
          DEALLOCATE (coef_pccomp%eigen(m)%eigenvectors, STAT = err)
        THROW(err.NE.0)

        IF (ASSOCIATED(coef_pccomp%eigen(m)%eigenvectors_t)) &
          DEALLOCATE (coef_pccomp%eigen(m)%eigenvectors_t, STAT = err)
        THROW(err.NE.0)
      ENDDO

      DEALLOCATE(coef_pccomp%eigen, STAT = err)
      THROW(err.NE.0)
    ENDIF

    IF (ASSOCIATED(coef_pccomp%pcreg)) THEN
      DO m = 1, coef_pccomp%fmv_pc_bands
        DO n = 1, coef_pccomp%fmv_pc_sets(m)
          IF (ASSOCIATED(coef_pccomp%pcreg(m,n)%coefficients)) &
            DEALLOCATE (coef_pccomp%pcreg(m,n)%coefficients, STAT = err)
          THROW(err.NE.0)

          IF (ASSOCIATED(coef_pccomp%pcreg(m,n)%coefficients_t)) &
            DEALLOCATE (coef_pccomp%pcreg(m,n)%coefficients_t, STAT = err)
          THROW(err.NE.0)

          IF (ASSOCIATED(coef_pccomp%pcreg(m,n)%predictindex)) &
            DEALLOCATE (coef_pccomp%pcreg(m,n)%predictindex, STAT = err)
          THROW(err.NE.0)

        ENDDO
      ENDDO

      DEALLOCATE (coef_pccomp%pcreg, STAT = err)
      THROW(err.NE.0)
    ENDIF

    IF (ASSOCIATED(coef_pccomp%fmv_pc_sets)) &
      DEALLOCATE(coef_pccomp%fmv_pc_sets, STAT = err)
    THROW(err.NE.0)


  CALL rttov_nullify_coef_pccomp(coef_pccomp)

  CATCH
END SUBROUTINE
