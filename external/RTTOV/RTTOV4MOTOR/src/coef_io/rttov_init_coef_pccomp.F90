! Description:
!> @file
!!   Initialise a PC coefficients structure.
!!   This should usually be called via rttov_init_coefs.
!
!> @brief
!!   Initialise a PC coefficients structure.
!!   This should usually be called via rttov_init_coefs.
!!
!! @param[out]     err            status on exit
!! @param[in]      coef           the optical depth coefficient structure
!! @param[in,out]  coef_pccomp    the PC coefficients structure to initialise
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_init_coef_pccomp(err, coef, coef_pccomp)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef, rttov_coef_pccomp
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim),      INTENT(OUT)   :: err
  TYPE(rttov_coef),        INTENT(INOUT) :: coef
  TYPE(rttov_coef_pccomp), INTENT(INOUT) :: coef_pccomp
!INTF_END
  INTEGER(KIND=jpim)  :: i, j

#include "rttov_errorreport.interface"

!- End of header --------------------------------------------------------
  TRY

  ALLOCATE (coef_pccomp%planck1_in(coef_pccomp%fmv_pc_nchn), &
            coef_pccomp%planck2_in(coef_pccomp%fmv_pc_nchn), STAT = err)
  THROWM(err.NE.0, "allocation of Planck1/2_in" )

  coef_pccomp%planck1_in(:) = coef%fc_planck_c1 * coef_pccomp%ff_cwn_in(:) ** 3
  coef_pccomp%planck2_in(:) = coef%fc_planck_c2 * coef_pccomp%ff_cwn_in(:)

  IF (coef_pccomp%fmv_pc_comp_pc < 5) THEN

    IF (coef_pccomp%fmv_pc_comp_pc == 4) THEN ! Original v12.1 PC NLTE coef file
      coef_pccomp%fmv_pc_nlte = 1
    ELSE
      coef_pccomp%fmv_pc_nlte = 0
    ENDIF

    ALLOCATE (coef_pccomp%co2_pc_ref(coef%nlevels), &
              coef_pccomp%n2o_pc_ref(coef%nlevels), &
              coef_pccomp%co_pc_ref(coef%nlevels),  &
              coef_pccomp%ch4_pc_ref(coef%nlevels), STAT = err)
    THROWM(err.NE.0, "allocation of PC gas ref profiles" )

    IF (coef%nco2 > 0) coef_pccomp%co2_pc_ref(:) = coef_pccomp%ref_pc_prfl_mr(:,1)
    IF (coef%nn2o > 0) coef_pccomp%n2o_pc_ref(:) = coef_pccomp%ref_pc_prfl_mr(:,2)
    IF (coef%nco > 0)  coef_pccomp%co_pc_ref(:)  = coef_pccomp%ref_pc_prfl_mr(:,3)
    IF (coef%nch4 > 0) coef_pccomp%ch4_pc_ref(:) = coef_pccomp%ref_pc_prfl_mr(:,4)
  ELSE
    ALLOCATE (coef_pccomp%co2_pc_min(coef%nlevels), &
              coef_pccomp%n2o_pc_min(coef%nlevels), &
              coef_pccomp%co_pc_min(coef%nlevels),  &
              coef_pccomp%ch4_pc_min(coef%nlevels), &
              coef_pccomp%co2_pc_max(coef%nlevels), &
              coef_pccomp%n2o_pc_max(coef%nlevels), &
              coef_pccomp%co_pc_max(coef%nlevels),  &
              coef_pccomp%ch4_pc_max(coef%nlevels), STAT = err)
    THROWM(err.NE.0, "allocation of PC gas min/max profiles" )

    IF (coef%nco2 > 0) THEN
      coef_pccomp%co2_pc_min(:) = coef_pccomp%lim_pc_prfl_gasmin(:,1)
      coef_pccomp%co2_pc_max(:) = coef_pccomp%lim_pc_prfl_gasmax(:,1)
    ENDIF
    IF (coef%nn2o > 0) THEN
      coef_pccomp%n2o_pc_min(:) = coef_pccomp%lim_pc_prfl_gasmin(:,2)
      coef_pccomp%n2o_pc_max(:) = coef_pccomp%lim_pc_prfl_gasmax(:,2)
    ENDIF
    IF (coef%nco > 0) THEN
      coef_pccomp%co_pc_min(:)  = coef_pccomp%lim_pc_prfl_gasmin(:,3)
      coef_pccomp%co_pc_max(:)  = coef_pccomp%lim_pc_prfl_gasmax(:,3)
    ENDIF
    IF (coef%nch4 > 0) THEN
      coef_pccomp%ch4_pc_min(:) = coef_pccomp%lim_pc_prfl_gasmin(:,4)
      coef_pccomp%ch4_pc_max(:) = coef_pccomp%lim_pc_prfl_gasmax(:,4)
    ENDIF
  ENDIF

! DAR: Added for RTTOV 11.2 - Allocate memory for transposed arrays for PC calcs
! We're doing it this way in order to maintain compatibility with the current
! coef files. The additional memory requirements for IASI are approx 40MB.
  DO j = 1, SIZE(coef_pccomp%pcreg(:,:), DIM=2)
    DO i = 1, SIZE(coef_pccomp%pcreg(:,:), DIM=1)
      ALLOCATE(coef_pccomp%pcreg(i,j)%coefficients_t( &
        SIZE(coef_pccomp%pcreg(i,j)%coefficients(1,:)), &
        SIZE(coef_pccomp%pcreg(i,j)%coefficients(:,1))))

      coef_pccomp%pcreg(i,j)%coefficients_t = &
        TRANSPOSE(coef_pccomp%pcreg(i,j)%coefficients)
    ENDDO
  ENDDO

  DO i = 1, SIZE(coef_pccomp%eigen(:))
    ALLOCATE(coef_pccomp%eigen(i)%eigenvectors_t( &
      SIZE(coef_pccomp%eigen(i)%eigenvectors(1,:)), &
      SIZE(coef_pccomp%eigen(i)%eigenvectors(:,1))))

    coef_pccomp%eigen(i)%eigenvectors_t = &
      TRANSPOSE(coef_pccomp%eigen(i)%eigenvectors)
  ENDDO

  ALLOCATE(coef_pccomp%noise_r(SIZE(coef_pccomp%noise)))
  coef_pccomp%noise_r = 1._jprb / coef_pccomp%noise

  CATCH
END SUBROUTINE 
