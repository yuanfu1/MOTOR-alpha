! Description:
!> @file
!!   Set up channel/profile and frequency indices.
!
!> @brief
!!   Set up channel/profile and frequency indices.
!!
!! @details
!!   This subroutine populates the chanprof structure and the frequencies
!!   array which are inputs to RTTOV-SCATT.
!!
!!   If the lchannel_subset array is not specified then nchannels must equal
!!   nprofiles x n_chan. If the lchannel_subset array is supplied then
!!   nchannels must be equal to the number of .TRUE. elements of
!!   lchannel_subset. These .TRUE. elements correspond to the channels to be
!!   simulated for each profile.
!!
!!
!! @param[out]     errorstatus      RTTOV error status
!! @param[in]      nprofiles        Number of profiles
!! @param[in]      n_chan           Number of channels in coefficient structure
!! @param[in]      coef_rttov       RTTOV coefficients structure
!! @param[in]      coef_scatt       RTTOV-SCATT coefficients structure
!! @param[in]      nchannels        Total number of channels to be simulated
!! @param[out]     chanprof         List of channels/profiles to simulate
!! @param[out]     frequencies      Frequency number for each channel being simulated
!! @param[in]      lchannel_subset  Logical flags to select a subset of channels, optional
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
subroutine rttov_scatt_setupindex (&
  & errorstatus,   & ! error status
  & nprofiles,     & ! number of profiles
  & n_chan,        & ! number of channels 
  & coef_rttov,    & ! coef structure read in from rttov coef file
  & coef_scatt,    & ! RTTOV_SCATT Coefficients
  & nchannels,     & ! number of calculated channels
  & chanprof,      & ! channels and profile numbers
  & frequencies,   & ! array, frequency number for each "channel"
  & lchannel_subset) ! OPTIONAL array of logical flags to indicate a subset of channels

! P. Bauer, P. Lopez, E. Moreau, D. Salmond   ECMWF    May 2004

! Modifications:
! 10/10/2006 Chris O'Dell : Made to work for any sensor (removed SSM/I hardcoding).
!                         : Changed the way lsprofiles2 is assigned to be slightly 
!                         : more general. Changed order of passed variables for 
!                         : style consistency.
! 11/11/2007 Alan Geer    : RTTOV9 version, cleaned.
! 06/06/2008 Alan Geer    : Added optional channel subsetting 
! 26/02/2009 W. Bell      : Changes for Windsat.
! 16/06/2020 Alan Geer    : Allow polarised scattering properties

!* KIND     
use parkind1,    only : jpim, jplm
use rttov_types, only : rttov_coefs, rttov_scatt_coef, rttov_chanprof
!INTF_OFF
use parkind1,    only : jprb
use rttov_const, only : errorstatus_success, errorstatus_fatal
!INTF_ON
implicit none

integer (kind=jpim), intent (out) :: errorstatus
integer (kind=jpim), intent ( in) :: nprofiles
integer (kind=jpim), intent ( in) :: n_chan 
type  (rttov_coefs), intent ( in) :: coef_rttov
type (rttov_scatt_coef), intent (in) :: coef_scatt ! RTTOV_SCATT Coefficients
integer (kind=jpim), intent ( in) :: nchannels
logical (kind=jplm), optional, intent ( in) :: lchannel_subset(nprofiles, n_chan)

 
integer  (kind=jpim), intent (out), dimension (nchannels) :: frequencies
type(rttov_chanprof), Intent (out), dimension (nchannels) :: chanprof ! Channel and profile indices

!INTF_END

#include "rttov_errorreport.interface"

integer (kind=jpim) :: i_prof, i_chan, j_chan, i_freq, polid
integer (kind=jpim) :: old_polid
real (kind=jprb)    :: cwn, old_cwn
logical (kind=jplm) :: luse(nprofiles, n_chan)
logical (kind=jplm) :: lpolarised_scattering

!- End of header ------------------------------------------------------

errorstatus = errorstatus_success

luse(:,:) = .true.
if (present(lchannel_subset)) luse = lchannel_subset

!* Set index arrays
j_chan = 0  ! counter to store calculated channels
old_cwn = 0.0_jprb
old_polid = -1

lpolarised_scattering = any(coef_scatt%mpol /= -1)

do i_prof = 1, nprofiles
  i_freq = 0
  do i_chan = 1, n_chan
   
    polid = coef_rttov % coef % fastem_polar (i_chan) + 1 ! polarisation ID of this channel
    cwn = coef_rttov % coef % ff_cwn(i_chan)

    if (lpolarised_scattering) then

      ! Scattering coefficients are tabulated per-frequency and per-polarisation

      i_freq = i_freq + 1
      if ( (coef_scatt%mpol(i_freq)+1) /= polid) then
        call rttov_errorreport (errorstatus_fatal, 'Incorrect channel polarisations in hydrotable', 'rttov_scatt ')
        return
      endif

    else

      ! Scattering coefficients are tabulated per-frequency and are independent of polarisation

      ! Below is a test to see if this channel represents a new frequency.
      ! The following 3 conditions must all hold in order to be same frequency as
      ! previous channel:
      !   1) Same central wave number as previous channel
      !      (NB check the absolute diff because some SSMI/S channel pairs have
      !          slightly different central wavenumbers: 1E-3 cm-1 == 0.03 GHz)
      !   2) Polarization ID eq 4,5,6 or 7 (single V or H polarisation only)
      !   3) Polarization ID different from last channel
      if ( (abs(cwn - old_cwn) > 1.E-3_jprb) .or. &
        & ((polid /= 4) .and. (polid /= 5) .and. (polid /= 6) .and. (polid /= 7)) &
        & .or. (polid == old_polid) ) then

        i_freq = i_freq + 1

      endif
    endif

    if( luse(i_prof,i_chan) ) then

      ! Profile and frequency number for each calculated channel
      j_chan = j_chan + 1
      chanprof   (j_chan)%chan  = i_chan
      frequencies(j_chan)       = i_freq
      chanprof   (j_chan)%prof  = i_prof

    endif

    old_polid = polid
    old_cwn = cwn

  end do
end do 

end subroutine rttov_scatt_setupindex
