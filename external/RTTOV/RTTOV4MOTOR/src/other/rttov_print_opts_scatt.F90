! Description:
!> @file
!!   Print out the contents of the rttov_options_scatt structure containing
!!   the configuration options for RTTOV-SCATT.
!
!> @brief
!!   Print out the contents of the rttov_options_scatt structure containing
!!   the configuration options for RTTOV-SCATT.
!!
!! @details
!!   If not supplied the output is written to the error_unit
!!   as set by rttov_errorhandling or the default if unset.
!!
!!   The optional text argument is printed at the top of the
!!   output.
!!
!! @param[in]   opts_scatt  options to configure the RTTOV-SCATT simulations
!! @param[in]   lu          logical unit for output, optional
!! @param[in]   text        additional text to print, optional
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
SUBROUTINE rttov_print_opts_scatt(opts_scatt, lu, text)

  USE rttov_types, ONLY : rttov_options_scatt
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_global, ONLY : error_unit
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_options_scatt), INTENT(IN)           :: opts_scatt  ! options
  INTEGER(KIND=jpim),        INTENT(IN), OPTIONAL :: lu          ! logical unit for print
  CHARACTER(LEN=*),          INTENT(IN), OPTIONAL :: text        ! additional text to print
!INTF_END

  INTEGER(KIND=jpim)  :: iu       ! logical unit for print
  CHARACTER(LEN=20)   :: tmp_text ! temporary string for formatting output

  iu = error_unit
  IF (PRESENT(lu)) iu = lu

  IF (PRESENT(text)) THEN
    WRITE(iu,'(/,a,a)') "RTTOV-SCATT options structure: ", TRIM(text)
  ELSE
    WRITE(iu,'(/,a)') "RTTOV-SCATT options structure"
  END IF

  WRITE(iu,'(a)')           "General config options"
  WRITE(iu,'(2x,a,l1)')     "apply_reg_limits      ", opts_scatt%config%apply_reg_limits
  WRITE(iu,'(2x,a,l1)')     "verbose               ", opts_scatt%config%verbose
  WRITE(iu,'(2x,a,l1)')     "do_checkinput         ", opts_scatt%config%do_checkinput
  WRITE(iu,'(2x,a,l1)')     "fix_hgpl              ", opts_scatt%config%fix_hgpl
  WRITE(iu,'(a)',advance='yes')

  WRITE(iu,'(a)')           "RT options"
  WRITE(iu,'(2x,a,l1)')     "lusercfrac            ", opts_scatt%lusercfrac
  WRITE(tmp_text,'(f10.4)') opts_scatt%cc_threshold
  WRITE(iu,'(2x,a,a)')      "cc_threshold          ", TRIM(ADJUSTL(tmp_text))
  WRITE(tmp_text,'(f10.4)') opts_scatt%ice_polarisation
  WRITE(iu,'(2x,a,a)')      "ice_polarisation      ", TRIM(ADJUSTL(tmp_text))
  WRITE(iu,'(2x,a,l1)')     "hydro_cfrac_tlad      ", opts_scatt%hydro_cfrac_tlad
  WRITE(iu,'(2x,a,l1)')     "zero_hydro_tlad       ", opts_scatt%zero_hydro_tlad
  WRITE(iu,'(2x,a,l1)')     "ozone_data            ", opts_scatt%ozone_data
  WRITE(iu,'(2x,a,l1)')     "use_t2m_opdep         ", opts_scatt%use_t2m_opdep
  WRITE(iu,'(2x,a,l1)')     "use_q2m               ", opts_scatt%use_q2m
  WRITE(iu,'(2x,a,l1)')     "addrefrac             ", opts_scatt%addrefrac
  WRITE(iu,'(2x,a,l1)')     "rad_down_lin_tau      ", opts_scatt%rad_down_lin_tau
  WRITE(iu,'(2x,a,l1)')     "dtau_test             ", opts_scatt%dtau_test
  WRITE(tmp_text,'(i3)') opts_scatt%fastem_version
  WRITE(iu,'(2x,a,a)')      "fastem_version        ", TRIM(ADJUSTL(tmp_text))
  WRITE(iu,'(2x,a,l1)')     "supply_foam_fraction  ", opts_scatt%supply_foam_fraction
  WRITE(tmp_text,'(i3)') opts_scatt%interp_mode
  WRITE(iu,'(2x,a,a)')      "interp_mode           ", TRIM(ADJUSTL(tmp_text))
  WRITE(iu,'(2x,a,l1)')     "reg_limit_extrap      ", opts_scatt%reg_limit_extrap
  WRITE(iu,'(2x,a,l1)')     "lgradp                ", opts_scatt%lgradp
  WRITE(iu,'(a)',advance='yes')

END SUBROUTINE rttov_print_opts_scatt
