! Description:
!> @file
!!   Print out the contents of a profile structure.
!
!> @brief
!!   Print out the contents of a profile structure.
!!
!! @details
!!   If not supplied the output is written to the error_unit
!!   as set by rttov_errorhandling or the default if unset.
!!
!!   The optional text argument is printed at the top of the
!!   output.
!!
!! @param[in]   profile   RTTOV profile structure
!! @param[in]   lu        logical unit for output, optional
!! @param[in]   text      additional text to print, optional
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
SUBROUTINE rttov_print_profile(profile, lu, text)

  USE rttov_types, ONLY : rttov_profile
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE rttov_global, ONLY : error_unit
  USE rttov_const, ONLY : gas_unit_specconc, gas_unit_ppmv
!INTF_ON

  IMPLICIT NONE

  TYPE(rttov_profile), INTENT(IN)           :: profile ! profile
  INTEGER(KIND=jpim),  INTENT(IN), OPTIONAL :: lu      ! logical unit for print
  CHARACTER(LEN=*),    INTENT(IN), OPTIONAL :: text    ! text for print
!INTF_END

  INTEGER(KIND=jpim)  :: iu  ! logical unit for print
  INTEGER(KIND=jpim)  :: l   ! level
  CHARACTER(LEN=7)    :: gas_units

  iu = error_unit
  IF (PRESENT(lu)) iu = lu

  IF (profile%gas_units == gas_unit_specconc) THEN
    gas_units = '(kg/kg)'
  ELSE
    gas_units = '(ppmv) '
  ENDIF

  IF (PRESENT(text)) THEN
    WRITE(iu,'(/,a,a)') "RTTOV profile structure: ", TRIM(text)
  ELSE
    WRITE(iu,'(/,a)') "RTTOV profile structure"
  ENDIF
  WRITE(iu,'(2x,a,a)')                      "id        ", TRIM(profile%id)
  WRITE(iu,'(2x,a,i4.4,"/",i2.2,"/",i2.2)') "date      ", profile%date
  WRITE(iu,'(2x,a,i2.2,":",i2.2,":",i2.2)') "time      ", profile%time
  WRITE(iu,'(2x,a,i4)') "number of levels ", profile%nlevels
  WRITE(iu,'(2x,a,i4)') "number of layers ", profile%nlayers
  WRITE(iu,'(a)') "Viewing conditions "
  WRITE(iu,'(2x,a,f10.3)') "zenangle     ", profile%zenangle
  WRITE(iu,'(2x,a,f10.3)') "azangle      ", profile%azangle
  WRITE(iu,'(2x,a,f10.3)') "sunzenangle  ", profile%sunzenangle
  WRITE(iu,'(2x,a,f10.3)') "sunazangle   ", profile%sunazangle
  WRITE(iu,'(2x,a,f10.3)') "elevation    ", profile%elevation
  WRITE(iu,'(2x,a,f10.3)') "latitude     ", profile%latitude
  WRITE(iu,'(2x,a,f10.3)') "longitude    ", profile%longitude
  WRITE(iu,'(2x,a,e14.6)') "Earth magnetic field   ", profile%be
  WRITE(iu,'(2x,a,e14.6)') "cos of Earth magnetic field and wave propagation ", profile%cosbk

  WRITE(iu,'(a)') "Black body cloud"
  WRITE(iu,'(2x,a,e14.6)') "cloud top pressure  (hPa) ", profile%ctp
  WRITE(iu,'(2x,a,e14.6)') "cloud fraction (0 - 1)    ", profile%cfraction

  WRITE(iu,'(a)') "Skin parameters"
  WRITE(iu,'(2x,a,i4)')    "surface type 0=land, 1=sea, 2=sea-ice    ", profile%skin%surftype
  WRITE(iu,'(2x,a,i4)')    "water type  0=fresh water, 1=ocean water ", profile%skin%watertype
  WRITE(iu,'(2x,a,e14.6)') "radiative skin temperature (K)           ", profile%skin%t
  WRITE(iu,'(2x,a,e14.6)') "snow_fraction (0 - 1)                    ", profile%skin%snow_fraction
  WRITE(iu,'(2x,a,e14.6)') "soil_moisture (m^3/m^3)                  ", profile%skin%soil_moisture
  WRITE(iu,'(2x,a,e14.6)') "salinity (%o)                            ", profile%skin%salinity
  WRITE(iu,'(2x,a,e14.6)') "foam fraction (0 - 1)                    ", profile%skin%foam_fraction
  WRITE(iu,'(2x,a)')       "land/sea-ice surface parameters for fastem"
  WRITE(iu,'(10x,10f10.4)')  profile%skin%fastem(:)

  WRITE(iu,'(a)') "2m air parameters"
  WRITE(iu,'(2x,a,e14.6)') "surface pressure (hPa) ", profile%s2m%p
  WRITE(iu,'(2x,a,e14.6)') "temperature (K)        ", profile%s2m%t
  WRITE(iu,'(2x,a,e14.6)') "water vapour "//gas_units//"   ", profile%s2m%q
  WRITE(iu,'(2x,a,e14.6)') "ozone "//gas_units//"          ", profile%s2m%o
  WRITE(iu,'(2x,a,e14.6)') "U wind component (m/s) ", profile%s2m%u
  WRITE(iu,'(2x,a,e14.6)') "V wind component (m/s) ", profile%s2m%v
  WRITE(iu,'(2x,a,e14.6)') "Wind fetch (m)         ", profile%s2m%wfetc

  IF (profile%gas_units == gas_unit_ppmv) THEN
    WRITE(iu, '(/,a)') "Gas units: ppmv wrt moist air"
  ELSE IF (profile%gas_units == gas_unit_specconc) THEN
    WRITE(iu, '(/,a)') "Gas units: kg/kg wrt moist air"
  ELSE
    WRITE(iu, '(/,a)') "Gas units: ppmv wrt dry air"
  ENDIF
  WRITE(iu,'(a5,a9,1x,a14,1x,a14)',advance='no') "level", "Pressure", "Temp", "WV "//gas_units
  IF (ASSOCIATED(profile%o3))  WRITE(iu,'(1x,a14)', advance='no') "O3 "//gas_units
  IF (ASSOCIATED(profile%co2)) WRITE(iu,'(1x,a14)', advance='no') "CO2 "//gas_units
  IF (ASSOCIATED(profile%n2o)) WRITE(iu,'(1x,a14)', advance='no') "N2O "//gas_units
  IF (ASSOCIATED(profile%co))  WRITE(iu,'(1x,a14)', advance='no') "CO "//gas_units
  IF (ASSOCIATED(profile%ch4)) WRITE(iu,'(1x,a14)', advance='no') "CH4 "//gas_units
  IF (ASSOCIATED(profile%so2)) WRITE(iu,'(1x,a14)', advance='no') "SO2 "//gas_units
  IF (ASSOCIATED(profile%clw)) WRITE(iu,'(1x,a14)', advance='no') "CLW (kg/kg)"
  WRITE(iu,'(a)', advance='yes')

  DO l = 1, profile%nlevels
    WRITE(iu,'(i4,f12.4,e14.6,1x,e14.6)', advance='no') l, profile%p(l), profile%t(l), profile%q(l)
    IF (ASSOCIATED(profile%o3)) &
    WRITE(iu,'(1x,e14.6)', advance='no') profile%o3(l)
    IF (ASSOCIATED(profile%co2)) &
    WRITE(iu,'(1x,e14.6)', advance='no') profile%co2(l)
    IF (ASSOCIATED(profile%n2o)) &
    WRITE(iu,'(1x,e14.6)', advance='no') profile%n2o(l)
    IF (ASSOCIATED(profile%co)) &
    WRITE(iu,'(1x,e14.6)', advance='no') profile%co(l)
    IF (ASSOCIATED(profile%ch4)) &
    WRITE(iu,'(1x,e14.6)', advance='no') profile%ch4(l)
    IF (ASSOCIATED(profile%so2)) &
    WRITE(iu,'(1x,e14.6)', advance='no') profile%so2(l)
    IF (ASSOCIATED(profile%clw)) &
    WRITE(iu,'(1x,e14.6)', advance='no') profile%clw(l)
    WRITE(iu,'(a)', advance='yes')
  ENDDO

  IF (ASSOCIATED(profile%aerosols)) THEN
    WRITE(iu,'(/,a)') "Aerosols"
    IF (profile%mmr_cldaer) THEN
      WRITE(iu,'(a5,1x,a21,2x,a)') "layer","Pressure  top  bottom", &
                                   "Aerosols mass mixing ratio in units of kg/kg"
    ELSE
      WRITE(iu,'(a5,1x,a21,2x,a)') "layer","Pressure  top  bottom", &
                                   "Aerosols mean number density in units of cm-3"
    ENDIF
    DO l = 1, profile%nlayers
      WRITE(iu,'(1x,i4,1x,f10.4,1x,f10.4,2x,20e14.6)') &
           & l, profile%p(l), profile%p(l+1), profile%aerosols(:,l)
    ENDDO
  ENDIF

  IF (ASSOCIATED(profile%cloud)) THEN
    WRITE(iu,'(/,a)') "Clouds"
    WRITE(iu,'(2x,a,i4)') "liquid cloud scheme       ", profile%clw_scheme
    WRITE(iu,'(2x,a,i4)') "clw to eff diameter param ", profile%clwde_param
    WRITE(iu,'(2x,a,i4)') "ice cloud scheme          ", profile%ice_scheme
    WRITE(iu,'(2x,a,i4)') "iwc to eff diameter param ", profile%icede_param
    IF (profile%mmr_cldaer) THEN
      WRITE(iu,'(a5,1x,a21,a,73x,3a14)') "layer", "Pressure  top  bottom", "  Cloud (kg/kg)", &
                                         " cfrac   ", "  clwde (um)", "  icede (um)"
    ELSE
      WRITE(iu,'(a5,1x,a21,a,73x,3a14)') "layer", "Pressure  top  bottom", "  Cloud (g/m3)", &
                                         " cfrac   ", "  clwde (um)", "  icede (um)"
    ENDIF
    DO l = 1, profile%nlayers
      WRITE(iu,'(1x,i4,1x,f10.4,1x,f10.4,2x,6e14.6,3f14.6)') l, profile%p(l), profile%p(l+1), profile%cloud(:,l), &
                                                             profile%cfrac(l), profile%clwde(l), profile%icede(l)
    ENDDO
  ENDIF
END SUBROUTINE rttov_print_profile
