!
Subroutine rttov_mieproc_tl (&        
     & nlevels,           &! in
     & nchannels,         &! in
     & nprofiles,         &! in
     & frequencies,       &! in
     & lprofiles,         &! in
     & lreflectivity,     &! in
     & profiles,          &! in
     & cld_profiles,      &! in
     & cld_profiles_tl,   &! in
     & opts_scatt,        &! in
     & coef_scatt,        &! in
     & coef_rttov,        &! in
     & chanprof,          &! in
     & scatt_aux,         &! inout
     & scatt_aux_tl)       ! inout 
  !
  ! Description:
  ! Calculates scattering parameters from hydrotables
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
  !    Copyright 2002, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  ! - Bauer, P., 2002: Microwave radiative transfer modeling in clouds and precipitation.
  !     Part I: Model description.
  !     NWP SAF Report No. NWPSAF-EC-TR-005, 27 pp.
  ! - Chevallier, F., and P. Bauer, 2003:
  !     Model rain and clouds over oceans: comparison with SSM/I observations.
  !     Mon. Wea. Rev., 131, 1240-1255.
  ! - Moreau, E., P. Bauer and F. Chevallier, 2002: Microwave radiative transfer modeling in clouds and precipitation.
  !     Part II: Model evaluation.
  !     NWP SAF Report No. NWPSAF-EC-TR-006, 27 pp.
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0       09/2002   Initial version     (P. Bauer, E. Moreau)
  !  1.1       05/2003   RTTOV7.3 compatible (F. Chevallier)
  !  1.2       11/2004   Clean-up            (P. Bauer)
  !  1.3       10/2006   Introduce interpolation to zero for scatt coeffs 
  !                      for small LWC (U.O'Keeffe). Cleaner version (A. Geer)
  !  1.4       11/2007   RTTOV9 version      (A. Geer)
  !  1.5       06/2008   Performance enhancements (D. Salmond)
  !  1.6       06/2008   Fix minor bug in small LWC interpolation (A. Geer)
  !  1.7       07/2008   Clear sky speed-ups (A. Geer)
  !  1.8       03/2009   Prevent extrapolation beyond table bounds (A. Geer)
  !  1.9       03/2010   Optimisation (A. Geer)
  !  1.10      06/2020   Simple polarised scattering (V. Barlakas and A. Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  !
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:

  Use rttov_types, Only :         &
       & rttov_profile_scatt_aux, &
       & rttov_profile,           &
       & rttov_profile_cloud,     &
       & rttov_scatt_coef,        &
       & rttov_options_scatt,     &
       & rttov_coef,              &
       & rttov_chanprof

  Use parkind1, Only : jpim, jplm
!INTF_OFF
  Use parkind1, Only : jprb
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
  use rttov_const, only : min_reflectivity
!INTF_ON
  Implicit None

!* Subroutine arguments:
  Integer (Kind=jpim), Intent (in) :: nlevels                 ! Number of  levels
  Integer (Kind=jpim), Intent (in) :: nchannels               ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent (in) :: nprofiles               ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: frequencies (nchannels) ! Frequency indices
  Integer (Kind=jpim), Intent (in) :: lprofiles   (nchannels) ! Profile indices
  Logical (kind=jplm), Intent (in) :: lreflectivity            ! radar simulator active?

  Type (rttov_profile),           Intent (in)    :: profiles        (nprofiles) ! profiles
  Type (rttov_profile_cloud),     Intent (in)    :: cld_profiles    (nprofiles) ! Cloud profiles
  Type (rttov_profile_cloud),     Intent (in)    :: cld_profiles_tl (nprofiles) ! Cloud profiles
  Type (rttov_options_scatt),     Intent (in)    :: opts_scatt              ! RTTOV_SCATT options
  Type (rttov_scatt_coef),        Intent (in)    :: coef_scatt              ! RTTOV_SCATT Coefficients

  Type (rttov_coef),              Intent (in)    :: coef_rttov           ! RTTOV Coefficients
  Type (rttov_chanprof),          Intent (in)    :: chanprof (nchannels) ! Channel and profile indices

  Type (rttov_profile_scatt_aux), Intent (inout) :: scatt_aux               ! Auxiliary profile variables
  Type (rttov_profile_scatt_aux), Intent (inout) :: scatt_aux_tl            ! Auxiliary profile variables

!INTF_END

!* Local variables:
  Integer (Kind=jpim) :: iwc(coef_scatt % nhydro), itemp(coef_scatt % nhydro), itype, ichan, ifreq, iprof, ilayer
  Integer (Kind=jpim) :: ifreq_last, iprof_last
  Real    (Kind=jprb) :: wc(coef_scatt % nhydro), temp, kp, ap, gp, zp, s_k, s_a, s_g, s_z, zln10
  Real    (Kind=jprb) :: kpp,app,gpp,zpp,kpp_tl,app_tl,gpp_tl,zpp_tl
  Real    (Kind=jprb) :: wc_tl(coef_scatt % nhydro), kp_tl, ap_tl, gp_tl, zp_tl
  Real    (Kind=jprb) :: cont(coef_scatt % nhydro), cont_tl(coef_scatt % nhydro), cont_min
  Real    (Kind=jprb) :: cfrac(coef_scatt % nhydro), cfrac_tl(coef_scatt % nhydro)
  Real    (Kind=jprb) :: pol_alpha
  Logical (Kind=jplm) :: ladd_polarisation
  Integer (Kind=jpim) :: ichanid, ichanid_last

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------
  
  if (lhook) call dr_hook('RTTOV_MIEPROC_TL',0_jpim,zhook_handle)

  zln10 = log (10.0_JPRB)
  cont_min = 10.0_JPRB ** ( (1.0_JPRB + coef_scatt % offset_water) / coef_scatt % scale_water )

  ladd_polarisation = opts_scatt % ice_polarisation > 0.0_jprb
  if (ladd_polarisation) then
    pol_alpha = (opts_scatt % ice_polarisation - 1.0_jprb)/(opts_scatt % ice_polarisation + 1.0_jprb)
  endif

  ! Loops over channels, levels, hydrometeor types
  nlayer_loop: do ilayer = 1, nlevels

    ifreq_last = -1
    iprof_last = -1
    ichanid_last = -1

    nchan_loop: do ichan = 1, nchannels

      iprof = lprofiles   (ichan)     
      ifreq = frequencies (ichan)  
      ichanid = chanprof(ichan) % chan

      if( scatt_aux % cfrac(iprof) > opts_scatt % cc_threshold ) then   

        if(iprof /= iprof_last ) then

          cont(:)     = scatt_aux % hydro (iprof,ilayer,:)
          cont_tl(:)  = scatt_aux_tl % hydro (iprof,ilayer,:)

          if (cld_profiles(iprof) % nhydro_frac == coef_scatt % nhydro) then
            cfrac(:)    = cld_profiles   (iprof) % hydro_frac(ilayer,:)
            cfrac_tl(:) = cld_profiles_tl(iprof) % hydro_frac(ilayer,:)
          else
            cfrac(:)    = cld_profiles   (iprof) % hydro_frac(ilayer,1)
            cfrac_tl(:) = cld_profiles_tl(iprof) % hydro_frac(ilayer,1)
          endif

          ntype_loop1: do itype = 1, coef_scatt % nhydro

            !* nearest index for hydrotable: LWC/IWC
            wc(itype)    = 0.0_JPRB
            wc_tl(itype) = 0.0_JPRB
            if (cont(itype) >= cont_min) then
              wc(itype)    = coef_scatt % scale_water * log10 (cont(itype)) - coef_scatt % offset_water
              wc_tl(itype) = coef_scatt % scale_water * cont_tl(itype) / (zln10 * cont(itype))   
            else if (cont(itype) < cont_min .and. cont(itype) >= 0.0_JPRB) then
              if (opts_scatt%zero_hydro_tlad .or. cont(itype) > 0.0_JPRB) then
                wc(itype) = cont(itype) / cont_min
                wc_tl(itype) = cont_tl(itype) / cont_min
              endif
            endif

            iwc(itype) = floor (wc(itype))
            if (iwc(itype) > coef_scatt % mwc - 1) then
              ! Prevent extrapolation 
              iwc(itype) = coef_scatt % mwc - 1
              wc(itype)  = coef_scatt % mwc
            endif

            !* nearest index for hydrotable: T (w/o melting layer)
            temp = profiles (iprof) % t (ilayer) - coef_scatt % offset_temp(itype)

            itemp(itype) = nint (temp)
            if (itemp(itype) <                      1) itemp(itype) = 1
            if (itemp(itype) > coef_scatt % mtemp - 1) itemp(itype) = coef_scatt % mtemp - 1

          enddo ntype_loop1

        endif

        if(ifreq /= ifreq_last .or. iprof /= iprof_last  .or. (ladd_polarisation .and. (ichanid /= ichanid_last))) then

          kpp=0.0_JPRB
          app=0.0_JPRB
          gpp=0.0_JPRB
          zpp=0.0_JPRB
          kpp_tl=0.0_JPRB
          app_tl=0.0_JPRB
          gpp_tl=0.0_JPRB
          zpp_tl=0.0_JPRB

          ntype_loop: do itype = 1, coef_scatt % nhydro 

            if (iwc(itype) >= 1) then

              s_k   = coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype)+1) &
                  & - coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype))
              s_a   = coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype)+1) &
                  & - coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype))
              s_g   = coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype)+1) &
                  & - coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype))
              if (associated(coef_scatt % zef)) then
                s_z = coef_scatt % zef (ifreq,itype,itemp(itype),iwc(itype)+1) &
                  & - coef_scatt % zef (ifreq,itype,itemp(itype),iwc(itype))
              endif

              kp    = coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype)) + s_k * (wc(itype) - iwc(itype))
              kp_tl = s_k * wc_tl(itype)

              ap    = coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype)) + s_a * (wc(itype) - iwc(itype))
              ap_tl = s_a * wc_tl(itype)

              gp    = coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype)) + s_g * (wc(itype) - iwc(itype))
              gp_tl = s_g * wc_tl(itype)

              if (associated(coef_scatt % zef)) then
                zp    = coef_scatt % zef (ifreq,itype,itemp(itype),iwc(itype)) + s_z * (wc(itype) - iwc(itype))
                zp_tl = s_z * wc_tl(itype)
              endif
                                   
            else

              ! For small water contents, interpolate linearly to zero 
              s_k   = coef_scatt % ext (ifreq,itype,itemp(itype),1) 
              s_a   = coef_scatt % ssa (ifreq,itype,itemp(itype),1) 
              s_g   = coef_scatt % asp (ifreq,itype,itemp(itype),1) 
              if (associated(coef_scatt % zef)) then
                s_z = coef_scatt % zef (ifreq,itype,itemp(itype),1)
              endif

              if (opts_scatt%zero_hydro_tlad) then
                kp  = s_k * wc(itype)
                kp_tl = s_k * wc_tl(itype)
              else
                kp  = max( s_k * wc(itype), 1E-10_JPRB )
                if (kp > 1E-10_JPRB) then
                  kp_tl = s_k * wc_tl(itype)
                else
                  kp_tl = 0.0_JPRB
                endif
              endif

              if (opts_scatt%zero_hydro_tlad .or. wc(itype) > 1E-10_JPRB) then
                ap    = s_a * wc(itype)
                ap_tl = s_a * wc_tl(itype)
                gp    = s_g * wc(itype) 
                gp_tl = s_g * wc_tl(itype)   
                if (associated(coef_scatt % zef)) then
                  zp    = s_z * wc(itype)
                  zp_tl = s_z * wc_tl(itype)
                endif
              else
                ap    = 0.0_JPRB 
                ap_tl = 0.0_JPRB
                gp    = 0.0_JPRB 
                gp_tl = 0.0_JPRB
                zp    = 0.0_JPRB
                zp_tl = 0.0_JPRB
              endif

            endif

            ! Simple way to get polarised scattering
            if (ladd_polarisation .and. coef_scatt%is_frozen(itype)) then
              if (coef_rttov%fastem_polar(ichanid) == 4) then ! h
                kp    = kp    * (1.0_jprb + pol_alpha)
                kp_tl = kp_tl * (1.0_jprb + pol_alpha)
              elseif (coef_rttov%fastem_polar(ichanid) == 3) then ! v
                kp    = kp    * (1.0_jprb - pol_alpha)
                kp_tl = kp_tl * (1.0_jprb - pol_alpha)
              endif
            endif

            if (lreflectivity) then
              ! Apply relevant cloud fraction (~random cloud overlap in attenuation computation)
              kpp   =kpp+kp*cfrac(itype)
              kpp_tl=kpp_tl+kp_tl*cfrac(itype)+kp*cfrac_tl(itype)
              zpp   =zpp+zp*cfrac(itype)
              zpp_tl=zpp_tl+zp_tl*cfrac(itype)+zp*cfrac_tl(itype)
              app=0.0_JPRB ! AJGDB pending move of the radar cloud fraction to a better location
              gpp=0.0_JPRB
              app_tl=0.0_JPRB
              gpp_tl=0.0_JPRB
            else
              kpp   =kpp+kp
              kpp_tl=kpp_tl+kp_tl
              app   =app+kp*ap
              gpp   =gpp+kp*ap*gp
              app_tl=app_tl+kp_tl*ap+kp*ap_tl
              gpp_tl=gpp_tl+kp_tl*ap*gp+kp*ap_tl*gp+kp*ap*gp_tl
            endif

          enddo ntype_loop
        endif

        ifreq_last=ifreq
        iprof_last=iprof
        ichanid_last=ichanid

        scatt_aux    % ext (ichan,ilayer) = scatt_aux    % ext (ichan,ilayer) + kpp
        scatt_aux_tl % ext (ichan,ilayer) = scatt_aux_tl % ext (ichan,ilayer) + kpp_tl

        scatt_aux    % ssa (ichan,ilayer) = scatt_aux    % ssa (ichan,ilayer) + app
        scatt_aux_tl % ssa (ichan,ilayer) = scatt_aux_tl % ssa (ichan,ilayer) + app_tl

        scatt_aux    % asm (ichan,ilayer) = scatt_aux    % asm (ichan,ilayer) + gpp
        scatt_aux_tl % asm (ichan,ilayer) = scatt_aux_tl % asm (ichan,ilayer) + gpp_tl

        if (zpp > 0.0_JPRB) then
          scatt_aux    % zef (ichan,ilayer) = 10.0_JPRB*log10(zpp)
          scatt_aux_tl % zef (ichan,ilayer) = 10.0_JPRB/log(10.0_JPRB)/zpp*zpp_tl
        else
          scatt_aux    % zef (ichan,ilayer) = min_reflectivity
          scatt_aux_tl % zef (ichan,ilayer) = 0.0_JPRB
        endif

      endif
    enddo nchan_loop
  enddo nlayer_loop

  do ilayer = 1, nlevels
  do ichan = 1, nchannels
    if (scatt_aux    % asm (ichan,ilayer) > 0.0_JPRB) then
    scatt_aux_tl % asm (ichan,ilayer) = scatt_aux_tl % asm (ichan,ilayer) / scatt_aux % ssa (ichan,ilayer) &
                        & - scatt_aux_tl % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer) &
                        & / scatt_aux    % ssa (ichan,ilayer) / scatt_aux % ssa (ichan,ilayer) 
    scatt_aux    % asm (ichan,ilayer) = scatt_aux    % asm (ichan,ilayer) / scatt_aux % ssa (ichan,ilayer)
     endif
    if (scatt_aux    % ssa (ichan,ilayer) > 0.0_JPRB) then
    scatt_aux_tl % ssa (ichan,ilayer) = (scatt_aux_tl % ssa (ichan,ilayer) * scatt_aux    % ext (ichan,ilayer) &
                        & -  scatt_aux    % ssa (ichan,ilayer) * scatt_aux_tl % ext (ichan,ilayer)) &
                        & / (scatt_aux    % ext (ichan,ilayer) * scatt_aux    % ext (ichan,ilayer)) 
    scatt_aux    % ssa (ichan,ilayer) =  scatt_aux    % ssa (ichan,ilayer) / scatt_aux    % ext (ichan,ilayer)
     endif
    if (scatt_aux % ext (ichan,ilayer) >= 20.0_JPRB) then
         scatt_aux_tl % ext (ichan,ilayer) =  0.0_JPRB
         scatt_aux    % ext (ichan,ilayer) = 20.0_JPRB
     endif
  enddo
  enddo

  if (lhook) call dr_hook('RTTOV_MIEPROC_TL',1_jpim,zhook_handle)

End subroutine rttov_mieproc_tl
