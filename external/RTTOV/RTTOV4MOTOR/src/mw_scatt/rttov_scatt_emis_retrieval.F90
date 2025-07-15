! Description:
!> @file
!!   Computes all-sky retrieved emissivities
!
!> @brief
!!   Computes all-sky retrieved emissivities
!!
!! @details
!!   See Baordo and Geer, 2016, DOI:10.1002/qj.2873
!!
!!   If a retrieval fails the corresponding output emissivity is negative.
!!
!! @param[in]     chanprof      channels and profiles simulated by RTTOV-SCATT
!! @param[in]     coefs         RTTOV coefficients structure
!! @param[in]     emis_terms    output radiances and corresponding BTs
!! @param[in]     obs_tb        observed BTs corresponding to simulated BTs
!! @param[out]    land_emis     output retrieved emissivities
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
subroutine rttov_scatt_emis_retrieval(chanprof, coefs, emis_terms, obs_tb, land_emis)

  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !   1.0    03/2013   Initial version     Fabrizio Baordo
  !   2.0    06/2013   RTTOV-11            Alan Geer
  !          01/2018   Radiances, 12.2, move out of RTTOV-SCATT  Alan Geer
  
  use parkind1, only : jpim, jprb
  use rttov_types, only : rttov_scatt_emis_retrieval_type, rttov_chanprof, rttov_coefs
!INTF_OFF    
  use yomhook, only: lhook , dr_hook
  use rttov_math_mod, only: planck
!INTF_ON
  implicit none

  !* Subroutine arguments:
  type (rttov_chanprof), intent(in) :: chanprof(:)
  type (rttov_coefs),    intent(in) :: coefs
  type (rttov_scatt_emis_retrieval_type),target,intent(in) :: emis_terms ! Down & up radiance source terms, Total transmittance
  real (kind=jprb), intent (in)  :: obs_tb (:)     ! Observed TB 
  real (kind=jprb), intent (out) :: land_emis (:)  ! Retrieved emissivity

!INTF_END

  ! local variables:
  integer (kind=jpim) :: nchannels

  integer (kind=jpim) :: ichan, jchan      
  
  real(kind=jprb) :: radobs    
  real(kind=jprb) :: tot_t_up, tot_t_down, tot_tau, denom
  real(kind=jprb), pointer :: tau_cld, up_cld, down_cld, tau_clr, up_clr, down_clr, c, bsfc

  real(kind=jprb) :: zhook_handle
  
  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_SCATT_EMIS_RETRIEVAL',0_jpim,zhook_handle)
 
  nchannels = size(chanprof)

  do ichan = 1, nchannels
    jchan = chanprof(ichan)%chan
    
    if ((obs_tb(ichan) > 1.0_JPRB) .and. (obs_tb(ichan) < 1000.0_JPRB)) then
    
      ! Convert observed TB to radiance
      call planck(coefs%coef%planck1(jchan), coefs%coef%planck2(jchan), obs_tb(ichan), radobs)
                  
      ! All-sky emissivity retrieval (reduces to clear-sky retreival when cfrac = 0.0)
      c         => emis_terms%cfrac(ichan)
      bsfc      => emis_terms%bsfc(ichan)
      tau_cld   => emis_terms%tau_cld(ichan)
      up_cld    => emis_terms%up_cld(ichan)
      down_cld  => emis_terms%down_cld(ichan)
      tau_clr   => emis_terms%tau_clr(ichan)
      up_clr    => emis_terms%up_clr(ichan)
      down_clr  => emis_terms%down_clr(ichan)

      ! Combine clear + cloudy
      tot_t_up   = c*up_cld           + (1.0_JPRB - c)*up_clr
      tot_t_down = c*down_cld*tau_cld + (1.0_JPRB - c)*down_clr*tau_clr
      tot_tau    = c*tau_cld          + (1.0_JPRB - c)*tau_clr

      denom = (tot_tau*bsfc - tot_t_down)
      if( abs(denom) > 1E-20_JPRB) then    
        ! Emissivity retrieval
        land_emis(ichan) = (radobs - tot_t_up - tot_t_down) / denom       
      else
        ! Protect against divison by zero, and return an obviously bad emissivity value
        land_emis(ichan) = -1.0_JPRB
      endif
      
    else
      land_emis(ichan) = -1.0_JPRB
    endif
  enddo

  if (lhook) call dr_hook('RTTOV_SCATT_EMIS_RETRIEVAL',1_jpim,zhook_handle)

end subroutine rttov_scatt_emis_retrieval
