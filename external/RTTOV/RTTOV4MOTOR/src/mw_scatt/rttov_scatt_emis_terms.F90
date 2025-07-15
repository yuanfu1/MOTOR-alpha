!
Subroutine rttov_scatt_emis_terms( &
        & ccthres,           &! in
        & chanprof,          &! in
        & coef_rttov,        &! in
        & scatt_aux,         &! in
        & emissivity,        &! in
        & transmission,      &! in
        & radiance,          &! in
        & sfc_terms,         &! in
        & emis_retrieval_terms) ! out

  ! Description:
  !
  ! Extract surface and TOA terms from the radiative transfer equation, for the purposes of 
  ! all-sky dynamic emissivity retrieval (see Baordo and Geer, 2016, DOI:10.1002/qj.2873)
  ! Caller must allocate all emis_retrieval_terms%arrays with size(chanprof)
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
  !    Copyright 2013, EUMETSAT, All Rights Reserved.
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !   1.0    02/2018   Initial version     Alan Geer
  
  use rttov_types, only: rttov_coefs, rttov_radiance2, rttov_chanprof, rttov_emissivity 
  use rttov_types, only: rttov_profile_scatt_aux, rttov_transmission, eddington_sfc_type, &
    & rttov_scatt_emis_retrieval_type
  use parkind1, only: jprb
!INTF_OFF
  use rttov_const, only : sensor_id_po
  use parkind1, only: jpim, jplm
  use yomhook, only: lhook , dr_hook
!INTF_ON
  implicit none

  !* Subroutine arguments:
  Real(Kind=jprb),          Intent (in)  :: ccthres
  Type(rttov_chanprof),     Intent (in)  :: chanprof(:) ! Channel and profile indices
  Type (rttov_coefs),       Intent (in)  :: coef_rttov               ! RTTOV Coefficients
  Type (rttov_profile_scatt_aux), Intent (in) :: scatt_aux  
  Type (rttov_emissivity),  Intent (in)  :: emissivity(size(chanprof))
  Type (rttov_transmission), Intent (in) :: transmission  
  Type (rttov_radiance2),    Intent (in) :: radiance                 ! Radiances
  Type (eddington_sfc_type),Intent (in)  :: sfc_terms ! Upward and downward radiance source terms, Total transmittances
  Type (rttov_scatt_emis_retrieval_type), Intent (inout)  :: emis_retrieval_terms

!INTF_END

  !* Local variables:
  Integer (Kind=jpim) :: nchannels
  Integer (Kind=jpim) :: iprof, ichan, ichan_act         
  Logical(Kind=jplm)  :: lpolarimetric(size(chanprof))
  real(kind=jprb) :: zhook_handle
  
  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_SCATT_EMIS_TERMS',0_jpim,zhook_handle)
 
  nchannels = size(chanprof)

  ! Clear sky terms
  emis_retrieval_terms%tau_clr(:)  = transmission%tau_total(:)
  
  emis_retrieval_terms%up_clr(:)   = radiance%upclear(:) - &
    & scatt_aux%bsfc(:) * transmission%tau_total(:) * emissivity(:)%emis_out ! AJGDB scatt_aux%emis_cld used (wrongly?) at ECMWF

  emis_retrieval_terms%down_clr(:) = radiance%dnclear(:)
 
  ! Surface black body radiance (possibly, as with other radiances, incorporating band correction)
  emis_retrieval_terms%bsfc(:) = scatt_aux%bsfc(:)
  
  do ichan = 1, nchannels
    ichan_act = chanprof(ichan)%chan
    iprof = chanprof(ichan)%prof 
    
    ! Identify polarimetric channels 
    lpolarimetric(ichan) = ( (coef_rttov % coef % id_sensor == sensor_id_po) &
                   & .and.   (coef_rttov % coef % fastem_polar(ichan_act) + 1_jpim >= 6_jpim) )

    if (scatt_aux%cfrac(iprof) > ccthres .and. .not. lpolarimetric(ichan)) then    
      emis_retrieval_terms%cfrac   (ichan) = scatt_aux%cfrac(iprof)
      emis_retrieval_terms%tau_cld (ichan) = sfc_terms%tau(ichan)
      emis_retrieval_terms%up_cld  (ichan) = sfc_terms%up(ichan) - &
        & scatt_aux%bsfc(ichan) * sfc_terms%tau(ichan) * scatt_aux%ems_cld(ichan)
      ! AJGDB cosmic contribution already considered within rttov_eddington
      emis_retrieval_terms%down_cld(ichan) = sfc_terms%down(ichan)
    else
      ! Revert to clear-sky R/T
      emis_retrieval_terms%cfrac   (ichan) = 0.0_JPRB
      emis_retrieval_terms%tau_cld (ichan) = 0.0_JPRB
      emis_retrieval_terms%up_cld  (ichan) = 0.0_JPRB
      emis_retrieval_terms%down_cld(ichan) = 0.0_JPRB
    endif   
  enddo
  
  if (lhook) call dr_hook('RTTOV_SCATT_EMIS_TERMS',1_jpim,zhook_handle)

end subroutine rttov_scatt_emis_terms
