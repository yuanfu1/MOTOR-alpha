!
Subroutine rttov_mieproc_ad (&        
     & nlevels,           &! in
     & nchannels,         &! in
     & nprofiles,         &! in
     & nprofilesad,       &! in
     & frequencies,       &! in
     & lprofiles,         &! in
     & lreflectivity,     &! in
     & profiles,          &! in
     & cld_profiles,      &! in
     & cld_profiles_ad,   &! inout
     & opts_scatt,        &! in
     & coef_scatt,        &! in
     & coef_rttov,        &! in
     & chanprof,          &! in  
     & scatt_aux,         &! inout
     & scatt_aux_ad)       ! inout 
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
  ! - Bauer, P., 2002: Microwave radiative transfer modeling in clouds and 
  !     precipitation. Part I: Model description.
  !     NWP SAF Report No. NWPSAF-EC-TR-005, 27 pp.
  ! - Chevallier, F., and P. Bauer, 2003:
  !     Model rain and clouds over oceans: comparison with SSM/I observations.
  !     Mon. Wea. Rev., 131, 1240-1255.
  ! - Moreau, E., P. Bauer and F. Chevallier, 2002: Microwave radiative transfer 
  !     modeling in clouds and precipitation. Part II: Model evaluation.
  !     NWP SAF Report No. NWPSAF-EC-TR-006, 27 pp.
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0    09/2002      Initial version         (E. Moreau)
  !  1.1    05/2003      RTTOV7.3 compatible     (F. Chevallier)
  !  1.3    03/2004      Polarimetry code added  (R. Saunders)
  !  1.4    11/2004      Clean-up                (P. Bauer)
  !  1.5    10/2006      Introduce interpolation to zero for scatt coeffs 
  !                      for small LWC (U.O'Keeffe). Cleaner version (A. Geer)
  !  1.6    11/2007      RTTOV9 versions         (A. Geer)
  !  1.7    06/2008      Fix 2 rarely occurring adjoint bugs (A. Geer)
  !  1.8    06/2008      Performance enhancements (D. Salmond)
  !  1.9    06/2008      Fix minor bug in small LWC interpolation (A. Geer)
  !  1.10   07/2008      Clear sky speed-ups     (A. Geer)
  !  1.11   03/2009      Prevent extrapolation beyond table bounds (A. Geer)
  !  1.12   03/2010      Optimisation (A. Geer)
  !  1.13   06/2020      Simple polarised scattering (V. Barlakas and A. Geer)
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
  use rttov_const, only: adk_adjoint, adk_k, min_reflectivity
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None

!* Subroutine arguments:
  Integer (Kind=jpim), Intent (in) :: nlevels                 ! Number oflevels
  Integer (Kind=jpim), Intent (in) :: nchannels               ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent (in) :: nprofiles               ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nprofilesad             ! Number of profiles in adjoint variables
  Integer (Kind=jpim), Intent (in) :: frequencies (nchannels) ! Frequency indices
  Integer (Kind=jpim), Intent (in) :: lprofiles   (nchannels) ! Profile indices
  Logical (kind=jplm), Intent (in) :: lreflectivity            ! radar simulator active?

  Type (rttov_profile),           Intent (in)    :: profiles        (nprofiles)
  Type (rttov_profile_cloud),     Intent (in)    :: cld_profiles    (nprofiles)   ! Cloud profiles
  Type (rttov_profile_cloud),     Intent (inout) :: cld_profiles_ad (nprofilesad) ! Cloud profiles
  Type (rttov_options_scatt),     Intent (in)    :: opts_scatt               ! RTTOV_SCATT options
  Type (rttov_scatt_coef),        Intent (in)    :: coef_scatt               ! RTTOV_SCATT Coefficients

  Type (rttov_coef),              Intent (in)    :: coef_rttov           ! RTTOV Coefficients
  Type (rttov_chanprof),          Intent (in)    :: chanprof (nchannels) ! Channel and profile indices

  Type (rttov_profile_scatt_aux), Intent (inout) :: scatt_aux                ! Auxiliary profile variables
  Type (rttov_profile_scatt_aux), Intent (inout) :: scatt_aux_ad             ! Auxiliary profile variables

!INTF_END

!* Local variables:
  Integer (Kind=jpim) :: iwc(coef_scatt % nhydro), itemp(coef_scatt % nhydro), itype, ichan, iprof, ifreq, ilayer
  Integer (Kind=jpim) :: ifreq_last, iprof_last
  Integer (Kind=jpim) :: adk, iprofad
  Real    (Kind=jprb) :: wc(coef_scatt % nhydro), temp, kp, ap, gp, zp, s_k, s_a, s_g, s_z, zln10
  Real    (Kind=jprb) :: kpp, app, gpp
  Real    (Kind=jprb) :: zpp, zpp_store(nchannels,nlevels), zpp_ad
  Real    (Kind=jprb) :: wc_ad, kp_ad, ap_ad, gp_ad, zp_ad
  Real    (Kind=jprb) :: cont(coef_scatt % nhydro), cont_ad(coef_scatt % nhydro), cont_min
  Real    (Kind=jprb) :: cfrac(coef_scatt % nhydro)
  Real    (Kind=jprb), pointer :: cfrac_ad
  Real    (Kind=jprb) :: pol_alpha
  Logical (Kind=jplm) :: ladd_polarisation
  Integer (Kind=jpim) :: ichanid, ichanid_last

  Real    (Kind=jprb) :: ext, ssa, asm!, zef

  REAL (KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------      

  if (lhook) call dr_hook('RTTOV_MIEPROC_AD',0_jpim,zhook_handle)

  if (nprofilesad == nprofiles) then 
    adk = adk_adjoint   ! Adjoint mode
  else if (nprofilesad == nchannels) then
    adk = adk_k         ! K mode
  endif 

  zln10 = log (10.0_JPRB)
  cont_min = 10.0_JPRB ** ( (1.0_JPRB + coef_scatt % offset_water) / coef_scatt % scale_water )

  ladd_polarisation = opts_scatt % ice_polarisation > 0.0_jprb
  if (ladd_polarisation) then
    pol_alpha = (opts_scatt % ice_polarisation - 1.0_jprb)/(opts_scatt % ice_polarisation + 1.0_jprb)
  endif

  !* Loops over channels, levels, hydrometeor types
  nlayer_loop1: do ilayer = 1, nlevels

    ifreq_last = -1
    iprof_last = -1
    ichanid_last = -1

    nchan_loop1: do ichan = 1, nchannels

      iprof = lprofiles   (ichan)
      ifreq = frequencies (ichan)  
      ichanid = chanprof(ichan) % chan

      if( scatt_aux % cfrac(iprof) > opts_scatt % cc_threshold ) then

        if(iprof /= iprof_last) then

          cont(:)  = scatt_aux % hydro (iprof,ilayer,:)

          if (cld_profiles(iprof) % nhydro_frac == coef_scatt % nhydro) then
            cfrac(:) = cld_profiles(iprof) % hydro_frac(ilayer,:)
          else
            cfrac(:) = cld_profiles(iprof) % hydro_frac(ilayer,1)
          endif

          ntype_loop1: do itype = 1, coef_scatt % nhydro

            !* Nearest index for hydrotable: LWC/IWC
            wc(itype) = 0.0_JPRB
            if (cont(itype) >= cont_min) then
              wc(itype) = coef_scatt % scale_water * log10 (cont(itype)) - coef_scatt % offset_water
            else if (cont(itype) < cont_min .and. cont(itype) >= 0.0_JPRB) then
              if (opts_scatt%zero_hydro_tlad .or. cont(itype) > 0.0_JPRB) then
                wc(itype) = cont(itype) / cont_min
              endif
            endif

            iwc(itype) = floor (wc(itype))
            if (iwc(itype) > coef_scatt % mwc - 1) then
              ! Prevent extrapolation 
              iwc(itype) = coef_scatt % mwc - 1
              wc(itype)  = coef_scatt % mwc
            endif

            !* Nearest index for hydrotable: T (w/o melting layer)
            temp = profiles (iprof) % t (ilayer) - coef_scatt % offset_temp(itype)

            itemp(itype) = nint (temp)
            if (itemp(itype) <                      1) itemp(itype) = 1
            if (itemp(itype) > coef_scatt % mtemp - 1) itemp(itype) = coef_scatt % mtemp - 1

          enddo ntype_loop1
        endif

        if(ifreq /= ifreq_last .or. iprof /= iprof_last .or. (ladd_polarisation .and. (ichanid /= ichanid_last))) then

          kpp=0.0_JPRB
          app=0.0_JPRB
          gpp=0.0_JPRB
          zpp=0.0_JPRB

          ntype_loop2: do itype = 1, coef_scatt % nhydro

            if (iwc(itype) >= 1) then

              s_k = coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype)+1) &
                & - coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype))
              s_a = coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype)+1) &
                & - coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype))
              s_g = coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype)+1) &
                & - coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype))
              if (associated(coef_scatt % zef)) then
                s_z = coef_scatt % zef (ifreq,itype,itemp(itype),iwc(itype)+1) &
                  & - coef_scatt % zef (ifreq,itype,itemp(itype),iwc(itype))
              endif

              kp  = coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype)) + s_k * (wc(itype) - iwc(itype))
              ap  = coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype)) + s_a * (wc(itype) - iwc(itype))
              gp  = coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype)) + s_g * (wc(itype) - iwc(itype))
              if (associated(coef_scatt % zef)) then
                zp = coef_scatt % zef (ifreq,itype,itemp(itype),iwc(itype)) + s_z * (wc(itype) - iwc(itype))
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
              else
                kp  = max( s_k * wc(itype), 1E-10_JPRB )
              endif
              if (opts_scatt%zero_hydro_tlad .or. wc(itype) > 1E-10_JPRB) then
                ap = s_a * wc(itype)
                gp = s_g * wc(itype)
                if (associated(coef_scatt % zef)) then
                  zp = s_z * wc(itype)
                endif
              else
                ap = 0.0_JPRB
                gp = 0.0_JPRB
                zp = 0.0_JPRB
              endif

            endif

            ! Simple way to get polarised scattering
            if (ladd_polarisation .and. coef_scatt%is_frozen(itype)) then
              if (coef_rttov%fastem_polar(ichanid) == 4) then ! h
                kp = kp * (1.0_jprb + pol_alpha)
              elseif (coef_rttov%fastem_polar(ichanid) == 3) then ! v
                kp = kp * (1.0_jprb - pol_alpha)
              endif
            endif

            if (lreflectivity) then
              ! Apply relevant cloud fraction (~random cloud overlap in attenuation computation)
              kpp=kpp+kp*cfrac(itype)
              zpp=zpp+zp*cfrac(itype)
              app=0.0_JPRB ! AJGDB pending move of the radar cloud fraction to a better location
              gpp=0.0_JPRB
            else
              kpp=kpp+kp
              app=app+kp * ap
              gpp=gpp+kp * ap * gp
            endif

          enddo ntype_loop2
        endif

        ifreq_last=ifreq
        iprof_last=iprof
        ichanid_last=ichanid

        scatt_aux % ext (ichan,ilayer) = scatt_aux % ext (ichan,ilayer) + kpp
        scatt_aux % ssa (ichan,ilayer) = scatt_aux % ssa (ichan,ilayer) + app
        scatt_aux % asm (ichan,ilayer) = scatt_aux % asm (ichan,ilayer) + gpp

        zpp_store(ichan,ilayer)=zpp
        if (zpp > 0.0_JPRB) then
          scatt_aux % zef (ichan,ilayer) = 10.0_JPRB*log10(zpp)
        else
          scatt_aux % zef (ichan,ilayer) = min_reflectivity
        endif

      endif
    enddo nchan_loop1
  enddo nlayer_loop1

  do ilayer = 1, nlevels
   do ichan = 1, nchannels
     ext  = scatt_aux % ext (ichan,ilayer)
     ssa  = scatt_aux % ssa (ichan,ilayer)
     asm  = scatt_aux % asm (ichan,ilayer)
     ! zef  = scatt_aux % zef (ichan,ilayer)

!* ADJOINT PART
     if (ext  >= 20.0_JPRB) then
        scatt_aux % ext (ichan,ilayer)  = 20.0_JPRB
        scatt_aux_ad % ext (ichan,ilayer) = 0.0_JPRB
     endif

     if (ssa  > 0.0_JPRB ) then
        scatt_aux % ssa (ichan,ilayer) =  ssa /  ext 
        scatt_aux_ad % ext (ichan,ilayer) = scatt_aux_ad % ext (ichan,ilayer) &
          & - scatt_aux_ad % ssa (ichan,ilayer) * ssa   &
          & / (ext * ext ) 
        scatt_aux_ad % ssa (ichan,ilayer) = scatt_aux_ad % ssa (ichan,ilayer) / ext 
     endif

     if (asm  > 0.0_JPRB ) then
        scatt_aux % asm (ichan,ilayer) = asm /  ssa 
        scatt_aux_ad % ssa (ichan,ilayer) = scatt_aux_ad % ssa (ichan,ilayer) &
          & - scatt_aux_ad % asm (ichan,ilayer) * asm        &
          & / (ssa * ssa )
        scatt_aux_ad     % asm (ichan,ilayer) = scatt_aux_ad  % asm (ichan,ilayer) / ssa 
     endif
     
   enddo
  enddo
  
!* Loops over channels, levels, hydrometeor types
  nlayer_loop2: do ilayer = nlevels, 1, -1

    ifreq_last = -1
    iprof_last = -1
    ichanid_last = -1

    nchan_loop2: do ichan = 1, nchannels

      iprof = lprofiles   (ichan)    
      if (adk == adk_adjoint) then
        iprofad = iprof  
      else if (adk == adk_k) then
        iprofad = ichan  
      endif
      ifreq = frequencies (ichan)  
      ichanid = chanprof(ichan) % chan

      if (scatt_aux % cfrac (iprof) > opts_scatt % cc_threshold) then

        if (zpp_store(ichan,ilayer) > 0.0_JPRB) then
          zpp_ad = 10.0_JPRB/log(10.0_JPRB)/zpp_store(ichan,ilayer) * scatt_aux_ad % zef (ichan,ilayer)
        else
          zpp_ad = 0.0_JPRB
        endif

        if(iprof /= iprof_last) then

          cont(:)  = scatt_aux % hydro (iprof,ilayer,:)

          if (cld_profiles(iprof) % nhydro_frac == coef_scatt % nhydro) then
            cfrac(:) = cld_profiles (iprof) % hydro_frac(ilayer,:)
          else
            cfrac(:) = cld_profiles (iprof) % hydro_frac(ilayer,1)
          endif

          ntype_loop3: do itype = coef_scatt % nhydro, 1, -1

            !* Nearest index for hydrotable: LWC/IWC
            wc(itype) = 0.0_JPRB
            if (cont(itype) >= cont_min) then
              wc(itype) = coef_scatt % scale_water * log10 (cont(itype)) - coef_scatt % offset_water
            else if (cont(itype) < cont_min .and. cont(itype) >= 0.0_JPRB) then
              if (opts_scatt%zero_hydro_tlad .or. cont(itype) > 0.0_JPRB) then
                wc(itype) = cont(itype) / cont_min
              endif
            endif

            iwc(itype) = floor (wc(itype))
            if (iwc(itype) > coef_scatt % mwc - 1) then
              ! Prevent extrapolation 
              iwc(itype) = coef_scatt % mwc - 1
              wc(itype)  = coef_scatt % mwc
            endif

            !* Nearest index for hydrotable: T (w/o melting layer)
            temp = profiles (iprof) % t (ilayer) - coef_scatt % offset_temp(itype)

            itemp(itype) = nint (temp)
            if (itemp(itype) <                      1) itemp(itype) = 1
            if (itemp(itype) > coef_scatt % mtemp - 1) itemp(itype) = coef_scatt % mtemp - 1

          enddo ntype_loop3

        endif

        ntype_loop4: do itype = coef_scatt % nhydro, 1, -1

          if (iwc(itype) >= 1) then
            s_k = coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype)+1) &
              & - coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype))
            s_a = coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype)+1) &
              & - coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype))
            s_g = coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype)+1) &
              & - coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype))
            if (associated(coef_scatt % zef)) then
              s_z = coef_scatt % zef (ifreq,itype,itemp(itype),iwc(itype)+1) &
                & - coef_scatt % zef (ifreq,itype,itemp(itype),iwc(itype))
            endif

            kp  = coef_scatt % ext (ifreq,itype,itemp(itype),iwc(itype)) + s_k * (wc(itype) - iwc(itype))
            ap  = coef_scatt % ssa (ifreq,itype,itemp(itype),iwc(itype)) + s_a * (wc(itype) - iwc(itype))
            gp  = coef_scatt % asp (ifreq,itype,itemp(itype),iwc(itype)) + s_g * (wc(itype) - iwc(itype))
            if (associated(coef_scatt % zef)) then
              zp  = coef_scatt % zef (ifreq,itype,itemp(itype),iwc(itype)) + s_z * (wc(itype) - iwc(itype))
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
            else
              kp  = max( s_k * wc(itype), 1E-10_JPRB )
            endif
            if (opts_scatt%zero_hydro_tlad .or. wc(itype) > 1E-10_JPRB) then
              ap = s_a * wc(itype)
              gp = s_g * wc(itype)
              if (associated(coef_scatt % zef)) then
                zp = s_z * wc(itype)
              endif
            else
              ap = 0.0_JPRB
              gp = 0.0_JPRB
              zp = 0.0_JPRB
            endif

          endif

          ! Simple way to get polarised scattering
          if (ladd_polarisation .and. coef_scatt%is_frozen(itype)) then
            if (coef_rttov%fastem_polar(ichanid) == 4) then ! h
              kp = kp * (1.0_jprb + pol_alpha)
            elseif (coef_rttov%fastem_polar(ichanid) == 3) then ! v
              kp = kp * (1.0_jprb - pol_alpha)
            endif
          endif

          kp_ad   = 0.0_JPRB
          ap_ad   = 0.0_JPRB
          gp_ad   = 0.0_JPRB
          zp_ad   = 0.0_JPRB

          wc_ad   = 0.0_JPRB

          if (lreflectivity) then

            ! Apply relevant cloud fraction (~random cloud overlap in attenuation computation)
            if (cld_profiles(iprof) % nhydro_frac == coef_scatt % nhydro) then
              cfrac_ad=>cld_profiles_ad (iprofad) % hydro_frac(ilayer,itype)
            else
              cfrac_ad=>cld_profiles_ad (iprofad) % hydro_frac(ilayer,1)
            endif

            zp_ad    = zp_ad    + zpp_ad*cfrac(itype)
            cfrac_ad = cfrac_ad + zpp_ad*zp
            kp_ad    = kp_ad    + scatt_aux_ad % ext (ichan,ilayer)*cfrac(itype)
            cfrac_ad = cfrac_ad + scatt_aux_ad % ext (ichan,ilayer)*kp

          else
            kp_ad = kp_ad + ap * gp * scatt_aux_ad % asm (ichan,ilayer)
            ap_ad = ap_ad + kp * gp * scatt_aux_ad % asm (ichan,ilayer)
            gp_ad = gp_ad + kp * ap * scatt_aux_ad % asm (ichan,ilayer)

            kp_ad = kp_ad +      ap * scatt_aux_ad % ssa (ichan,ilayer)
            ap_ad = ap_ad +      kp * scatt_aux_ad % ssa (ichan,ilayer)

            kp_ad = kp_ad + scatt_aux_ad % ext (ichan,ilayer)
          endif

          ! Simple way to get polarised scattering
          if (ladd_polarisation .and. coef_scatt%is_frozen(itype)) then
            if (coef_rttov%fastem_polar(ichanid) == 4) then ! h
              kp_ad = kp_ad * (1.0_jprb + pol_alpha)
            elseif (coef_rttov%fastem_polar(ichanid) == 3) then ! v
              kp_ad = kp_ad * (1.0_jprb - pol_alpha)
            endif
          endif

          if (iwc(itype) >= 1) then

            wc_ad = wc_ad + s_g * gp_ad
            gp_ad = 0.0_JPRB
            wc_ad = wc_ad + s_a * ap_ad
            ap_ad = 0.0_JPRB
            wc_ad = wc_ad + s_k * kp_ad
            kp_ad = 0.0_JPRB
            if (associated(coef_scatt % zef)) then
              wc_ad = wc_ad + s_z * zp_ad
              zp_ad = 0.0_JPRB
            endif

          else

            if (opts_scatt%zero_hydro_tlad .or. wc(itype) > 1E-10_JPRB) then
              wc_ad = wc_ad + s_g * gp_ad
              gp_ad = 0.0_JPRB
              wc_ad = wc_ad + s_a * ap_ad
              ap_ad = 0.0_JPRB
              if (associated(coef_scatt % zef)) then
                wc_ad = wc_ad + s_z * zp_ad
                zp_ad = 0.0_JPRB
              endif
            endif
            if (opts_scatt%zero_hydro_tlad) then
              wc_ad = wc_ad + s_k * kp_ad
            else
              if( kp > 1E-10_JPRB) wc_ad = wc_ad + s_k * kp_ad
            endif
            kp_ad = 0.0_JPRB

          endif

          cont_ad(itype) = 0.0_JPRB
          if (cont(itype) >= cont_min) then
            cont_ad(itype) = cont_ad(itype) + coef_scatt % scale_water * wc_ad / (zln10 * cont(itype))
          else if (cont(itype) < cont_min .and. cont(itype) >= 0.0_JPRB) then
            if (opts_scatt%zero_hydro_tlad .or. cont(itype) > 0.0_JPRB) then
              cont_ad(itype) = cont_ad(itype) + wc_ad / cont_min
            endif
!           else
!             cont_ad(itype) = cont_ad(itype) + 0.0_JPRB
          endif

          wc_ad = 0.0_JPRB

        enddo ntype_loop4

        where ((scatt_aux % hydro (iprof,ilayer,:) > 0.0_JPRB) .or. opts_scatt%zero_hydro_tlad) &
         & scatt_aux_ad % hydro (iprofad,ilayer,:) = scatt_aux_ad % hydro (iprofad,ilayer,:) + cont_ad(:)

        ifreq_last=ifreq
        iprof_last=iprof
        ichanid_last=ichanid

      endif 
    enddo nchan_loop2
  enddo nlayer_loop2

  if (lhook) call dr_hook('RTTOV_MIEPROC_AD',1_jpim,zhook_handle)

End subroutine rttov_mieproc_ad
