! Description:
!> @file
!!   Print out the contents of an RTTOV-SCATT cloud profile structure.
!
!> @brief
!!   Print out the contents of an RTTOV-SCATT cloud profile structure.
!!
!! @details
!!   If not supplied the output is written to the error_unit
!!   as set by rttov_errorhandling or the default if unset.
!!
!!   The optional text argument is printed at the top of the
!!   output.
!!
!! @param[in]   cld_profile   RTTOV-SCATT cloud profile structure
!! @param[in]   lu            logical unit for output, optional
!! @param[in]   text          additional text to print, optional
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
SUBROUTINE rttov_print_cld_profile(cld_profile, lu, text)

  USE rttov_types, ONLY : rttov_profile_cloud
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_global, ONLY : error_unit
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_profile_cloud), INTENT(IN)           :: cld_profile ! cloud profile
  INTEGER(KIND=jpim),        INTENT(IN), OPTIONAL :: lu          ! logical unit for print
  CHARACTER(LEN=*),          INTENT(IN), OPTIONAL :: text        ! text for print
!INTF_END

  INTEGER(KIND=jpim)  :: iu  ! logical unit for print
  INTEGER(KIND=jpim)  :: l   ! level
  INTEGER(KIND=jpim)  :: i_type

  iu = error_unit
  IF (PRESENT(lu)) iu = lu

  IF (PRESENT(text)) THEN
    WRITE(iu,'(/,a,a)') "RTTOV-SCATT cloud profile structure: ", TRIM(text)
  ELSE
    WRITE(iu,'(/,a)') "RTTOV-SCATT cloud profile structure"
  ENDIF
  WRITE(iu,'(2x,a,i4)') "number of levels ", cld_profile%nlevels
  WRITE(iu,'(2x,a,e14.6)') "user average cloud fraction (0 - 1)    ", cld_profile%cfrac
  WRITE(iu,'(2x,a,30i4)') "hydro units (0=default=kg/kg; 1,2=kg/m2/s)", cld_profile%flux_conversion

  WRITE(iu,'(a5,1x,a21)',advance='no') &
    "level", "Pressure  top  bottom"
  DO i_type = 1, cld_profile%nhydro
    WRITE(iu,'(a15)', advance='no') "    hydro      "
  ENDDO
  DO i_type = 1, cld_profile%nhydro_frac
    WRITE(iu,'(a15)', advance='no') "    frac       "
  ENDDO
  WRITE(iu,'(a)', advance='yes')

  DO l = 1, cld_profile%nlevels
    WRITE(iu,'(1x,i4,2(1x,f10.4),1x,f9.6)', advance='no') &
      l, cld_profile%ph(l), cld_profile%ph(l+1)
    DO i_type = 1, cld_profile%nhydro
      WRITE(iu,'(1x,e14.6)', advance='no') cld_profile%hydro(l,i_type)
    ENDDO
    DO i_type = 1, cld_profile%nhydro_frac
      WRITE(iu,'(1x,e14.6)', advance='no') cld_profile%hydro_frac(l,i_type)
    ENDDO
    WRITE(iu,'(a)', advance='yes')
  ENDDO
END SUBROUTINE rttov_print_cld_profile
