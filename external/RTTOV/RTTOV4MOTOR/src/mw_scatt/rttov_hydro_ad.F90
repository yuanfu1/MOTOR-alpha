!
Subroutine rttov_hydro_ad (&
      & ccthres,           &! in
      & nlevels,           &! in
      & nprofiles,         &! in
      & nprofilesad,       &! in
      & nchannels,         &! in
      & usercfrac,         &! in
      & lreflectivity,     &! in
      & chanprof,          &! in
      & presf,             &! in
      & presfad,           &! inout
      & profiles,          &! in  
      & profiles_ad,       &! inout  
      & cld_profiles,      &! in 
      & cld_profiles_ad,   &! inout 
      & opts_scatt,        &! in
      & coef_scatt,        &! in
      & scatt_aux,         &! inout 
      & scatt_aux_ad)       ! inout 

  !
  ! Description:
  ! Initialises the hydrometeor profile in the cloudy sub-column.
  !
  ! Copyright:
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 7 December 2016, between
  !    EUMETSAT and the Met Office, UK, by one or more painiscattrtners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, DWD and MeteoFrance.
  !
  !    Copyright 2012, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  ! 
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !   1.0    09/2012   Initial version     (A Geer)
  !   2.0    10/2018   Flexible hydrometeors (A Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  ! 
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:
  Use rttov_types, Only :         &
       & rttov_scatt_coef,        &
       & rttov_options_scatt,     &
       & rttov_profile_scatt_aux, &
       & rttov_profile,           &
       & rttov_profile_cloud,     &
       & rttov_chanprof 

  Use parkind1, Only : jpim, jplm, jprb

!INTF_OFF
  Use rttov_const, Only:      &
       & rgc,                 &
       & mair,                &
       & rho_rain,            &
       & rho_snow,            &
       & conv_rain,           &
       & conv_sp,             &
       & adk_adjoint,         &
       & adk_k,               &
       & min_cfrac_radar

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None
  
!* Subroutine arguments:
  Real    (Kind=jprb), Intent (in) :: ccthres
  Integer (Kind=jpim), Intent (in) :: nlevels   ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nprofilesad  ! Number of profiles in adjoint
  Integer (Kind=jpim), Intent (in) :: nchannels  ! Number of channels
  Logical (Kind=jplm), Intent (in) :: usercfrac               ! User has supplied cloud fraction
  Logical (kind=jplm), Intent (in) :: lreflectivity            ! radar simulator active?
  Type(rttov_chanprof), Intent(in) :: chanprof    (nchannels) ! Channel and profile indices
  Real (Kind=jprb),                 Intent (in)    :: presf(nprofiles,nlevels) ! Pressure levels [hPa]
  Real (Kind=jprb),                 Intent (inout) :: presfad(nprofilesad,nlevels) ! Pressure levels [hPa]
  Type (rttov_profile),             Intent (in)    :: profiles (nprofiles)     ! Atmospheric profiles
  Type (rttov_profile),             Intent (inout) :: profiles_ad (nprofilesad)     ! Atmospheric profiles
  Type (rttov_options_scatt),       Intent (in)    :: opts_scatt               ! RTTOV_SCATT options
  Type (rttov_scatt_coef),          Intent (in)    :: coef_scatt               ! RTTOV_SCATT Coefficients
  Type (rttov_profile_cloud),       Intent (in)    :: cld_profiles (nprofiles) ! Cloud profiles
  Type (rttov_profile_cloud),       Intent (inout) :: cld_profiles_ad (nprofilesad) ! Cloud profiles AD
  Type (rttov_profile_scatt_aux),   Intent (inout) :: scatt_aux                ! Auxiliary profile variables
  Type (rttov_profile_scatt_aux),   Intent (inout) :: scatt_aux_ad             ! Auxiliary profile variables
!INTF_END

!* Local variables
  Integer (Kind=jpim) :: ilayer, iprof, itype
  Real    (Kind=jprb) :: rho, conv(2)
  Integer (Kind=jpim) :: iprofad, adk
  Real    (Kind=jprb) :: de2mr, mmr_to_density(nprofiles,nlevels), mmr_to_density_ad
  Real    (Kind=jprb) :: cfrac(nprofiles)
  Real    (Kind=jprb) :: hydro_frac_weighted_sum(nprofiles)                  ! Hydrometeor fraction times amount, summed
  Real    (Kind=jprb) :: hydro_frac_weighted(nlevels,coef_scatt % nhydro)    ! Hydrometeor fraction times amount
  Real    (Kind=jprb) :: hydro_weights(nlevels,coef_scatt % nhydro,nprofiles)! Hydrometeor amounts [g m^-2]
  Real    (Kind=jprb) :: hydro_column(nprofiles)                             ! Total hydrometeor amounts [g m^-2]
  Real    (Kind=jprb) :: hydro_frac_weighted_ad(nlevels,coef_scatt % nhydro) ! Hydrometeor fraction times amount
  Real    (Kind=jprb) :: hydro_weights_ad(nlevels,coef_scatt % nhydro)       ! Hydrometeor amounts [g m^-2]
  Real    (Kind=jprb) :: hydro_column_ad                                     ! Total hydrometeor amounts [g m^-2]
  Real    (Kind=jprb), dimension (nprofiles,nlevels,coef_scatt%nhydro) :: hydro_scale, hydro_precf
  Real    (Kind=jprb), pointer :: cfrac_column(:), cfrac_column_ad(:)

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_HYDRO_AD',0_jpim,zhook_handle)

  if (nprofilesad == nprofiles) then 
    adk = adk_adjoint   ! Adjoint mode
  else if (nprofilesad == nchannels) then
    adk = adk_k         ! K mode
  endif 

  de2mr =  1.0E+02_JPRB * mair / rgc

  !* Initialise cloud and rain properties of the cloudy/rainy column
  do ilayer=1,nlevels
    do iprof = 1, nprofiles  
      scatt_aux % hydro (iprof,ilayer,:) = cld_profiles (iprof) % hydro (ilayer,:)
    enddo
  enddo

  !* Change units
  do ilayer = 1, nlevels
    do iprof = 1, nprofiles

      mmr_to_density(iprof,ilayer) = presf (iprof,ilayer) * de2mr / profiles (iprof) % t (ilayer)

      do itype = 1, coef_scatt % nhydro
        if( cld_profiles (iprof) % flux_conversion (itype) == 0 ) then

          !* Save for adjoint
          hydro_scale (iprof,ilayer,itype) = scatt_aux % hydro (iprof,ilayer,itype)

          !* Condensate from g/g to g/m^3
          scatt_aux % hydro (iprof,ilayer,itype) = scatt_aux % hydro (iprof,ilayer,itype) * mmr_to_density(iprof,ilayer)

        else

          !* Rates from kg/m^2/s to g/m^3
          if( cld_profiles (iprof) % flux_conversion (itype) == 1) then
            rho  = rho_rain
            conv = conv_rain
          else
            rho  = rho_snow
            conv = conv_sp
          endif

          scatt_aux % hydro (iprof,ilayer,itype) = scatt_aux    % hydro (iprof,ilayer,itype) / rho

          !* Save for adjoint
          hydro_scale (iprof,ilayer,itype) = scatt_aux % hydro (iprof,ilayer,itype)

          if (scatt_aux % hydro (iprof,ilayer,itype) > 0.0_JPRB) then
            scatt_aux % hydro (iprof,ilayer,itype) = (scatt_aux % hydro (iprof,ilayer,itype) * conv(1))**(conv(2))
          endif

        endif
      enddo
    enddo
  enddo

  !* Save for adjoint
  hydro_precf = scatt_aux % hydro

  scatt_aux % cfrac (:)   = 0.0_JPRB

  do iprof = 1, nprofiles

    if (lreflectivity) then

      ! Cfrac has no meaning for radar simulator, but set it to 1 to pass checks elsewhere
      scatt_aux % cfrac (iprof) = 1.0_JPRB

      ! Each hydrometeor is weighted by its own subgrid fraction
      do itype = 1, coef_scatt % nhydro

        if (cld_profiles (iprof) % nhydro_frac == 1) then
          cfrac_column => cld_profiles (iprof) % hydro_frac (:,1)
        else
          cfrac_column => cld_profiles (iprof) % hydro_frac (:,itype)
        endif

        where(cfrac_column > min_cfrac_radar)
          scatt_aux % hydro (iprof,:,itype) = scatt_aux % hydro (iprof,:,itype) / cfrac_column
        elsewhere
          scatt_aux % hydro (iprof,:,itype) = 0.0_JPRB
        endwhere

      enddo

    else

      if( usercfrac ) then

        !* User-supplied cloud fraction
        scatt_aux % cfrac (iprof) = cld_profiles (iprof) % cfrac

      else

        !* Calculate a hydrometeor-weighted average cloudy sky fraction

        !* Partial column of hydrometeors in g/m^2
        do itype = 1, coef_scatt % nhydro
          hydro_weights (:,itype,iprof) = scatt_aux % hydro (iprof,:,itype) * scatt_aux % dz (iprof,:)
        enddo

        hydro_column (iprof) = sum(hydro_weights(:,:,iprof))

        !* Weighted mean cloud fraction
        if (hydro_column (iprof) > 1e-10_JPRB) then

          if (cld_profiles (iprof) % nhydro_frac == 1) then
            ! One cloud fraction applies to all hydrometeors
            do itype = 1, coef_scatt % nhydro
              hydro_frac_weighted (:,itype) = &
                cld_profiles (iprof) % hydro_frac (:,1) * hydro_weights(:,itype,iprof)
            enddo
          else
            ! One hydrometeor fraction per hydrometeor
            hydro_frac_weighted (:,:) = cld_profiles (iprof) % hydro_frac (:,:) * hydro_weights (:,:,iprof)
          endif

          !* Save for adjoint
          hydro_frac_weighted_sum(iprof) = sum(hydro_frac_weighted)

          scatt_aux % cfrac (iprof) = sum(hydro_frac_weighted) / hydro_column (iprof)
          cfrac(iprof) = scatt_aux % cfrac (iprof)
          if ( scatt_aux % cfrac (iprof) < 0.0_JPRB ) scatt_aux % cfrac (iprof) = 0.0_JPRB
          if ( scatt_aux % cfrac (iprof) > 1.0_JPRB ) scatt_aux % cfrac (iprof) = 1.0_JPRB

        else
          scatt_aux % cfrac (iprof) = 0.0_JPRB
        endif

      endif

      !* Partition all cloud and rain into the cloudy column
      if (scatt_aux % cfrac (iprof) > ccthres) then
        scatt_aux % hydro (iprof,:,:) = scatt_aux % hydro (iprof,:,:) / scatt_aux % cfrac (iprof)
      else
        scatt_aux % hydro (iprof,:,:) = 0.0_JPRB
      endif

    endif
  enddo

  ! ADJOINT PART

  ! Calculate a hydrometeor-weighted average cloudy sky fraction 
  Do iprofad = 1, nprofilesad

    if (adk == adk_adjoint) then
      iprof = iprofad
    else if (adk == adk_k) then
      iprof = chanprof(iprofad) % prof
    endif

    if (lreflectivity) then

      ! Each hydrometeor is weighted by its own subgrid fraction
      do itype = 1, coef_scatt % nhydro

        if (cld_profiles (iprof) % nhydro_frac == 1) then
          cfrac_column    => cld_profiles    (iprof)   % hydro_frac (:,1)
          cfrac_column_ad => cld_profiles_ad (iprofad) % hydro_frac (:,1)
        else
          cfrac_column    => cld_profiles    (iprof)   % hydro_frac (:,itype)
          cfrac_column_ad => cld_profiles_ad (iprofad) % hydro_frac (:,itype)
        endif

        where(cfrac_column > min_cfrac_radar)

          cfrac_column_ad = cfrac_column_ad - scatt_aux_ad % hydro (iprofad,:,itype) &
                        & * hydro_precf (iprof,:,itype) / (cfrac_column**2)

          scatt_aux_ad % hydro (iprofad,:,itype) = scatt_aux_ad % hydro (iprofad,:,itype) / cfrac_column

        elsewhere
          scatt_aux_ad % hydro (iprofad,:,itype) = 0.0_JPRB
        endwhere

      enddo

      scatt_aux_ad % cfrac (iprofad) = 0.0_JPRB

    else

      !* Partition all cloud and rain into the cloudy column
      if (scatt_aux % cfrac (iprof) > ccthres) Then

        scatt_aux_ad % cfrac (iprofad) = scatt_aux_ad % cfrac (iprofad)        &
          & - sum ( scatt_aux_ad % hydro (iprofad,:,:) * hydro_precf (iprof,:,:) ) &
          & / ( scatt_aux % cfrac (iprof) ** 2 )

        scatt_aux_ad % hydro (iprofad,:,:) = scatt_aux_ad % hydro (iprofad,:,:) / scatt_aux % cfrac (iprof)

      else
        scatt_aux_ad % hydro (iprofad,:,:) = 0.0_JPRB
      endif

      if( usercfrac ) then

        !* User-supplied cloud fraction
        cld_profiles_ad (iprofad) % cfrac = scatt_aux_ad % cfrac (iprofad)
        scatt_aux_ad % cfrac (iprofad) = 0.0_JPRB

      else

        !* Weighted mean cloud fraction
        if (hydro_column(iprof) > 1e-10_JPRB) then

          if ( cfrac(iprof) < 0.0_JPRB ) scatt_aux_ad % cfrac (iprofad) = 0.0_JPRB
          if ( cfrac(iprof) > 1.0_JPRB ) scatt_aux_ad % cfrac (iprofad) = 0.0_JPRB

          hydro_frac_weighted_ad = scatt_aux_ad % cfrac (iprofad) / hydro_column (iprof)

          hydro_column_ad = -1.0_JPRB * scatt_aux_ad % cfrac (iprofad) * &
            hydro_frac_weighted_sum(iprof) / ( hydro_column(iprof) ** 2 )

          if (cld_profiles (iprof) % nhydro_frac == 1) then

            ! One cloud fraction applies to all hydrometeors
            do itype = 1, coef_scatt % nhydro
              cld_profiles_ad (iprofad) % hydro_frac (:,1) = cld_profiles_ad (iprofad) % hydro_frac (:,1) + &
                hydro_weights (:,itype,iprof) * hydro_frac_weighted_ad (:,itype)

              hydro_weights_ad(:,itype) = &
                cld_profiles (iprof) % hydro_frac (:,1) * hydro_frac_weighted_ad (:,itype)
            enddo

          else

            ! One hydrometeor fraction per hydrometeor
            cld_profiles_ad (iprofad) % hydro_frac (:,:) = cld_profiles_ad (iprofad) % hydro_frac (:,:) + &
              hydro_weights (:,:,iprof) * hydro_frac_weighted_ad (:,:)

            hydro_weights_ad(:,:) = cld_profiles (iprof) % hydro_frac (:,:) * hydro_frac_weighted_ad (:,:)

          endif

        else
          hydro_weights_ad(:,:) = 0.0_JPRB
          hydro_column_ad       = 0.0_JPRB
        endif

        scatt_aux_ad % cfrac (iprofad) = 0.0_JPRB

        !* Partial column of hydrometeors in g/m^2

        hydro_weights_ad(:,:) = hydro_weights_ad(:,:) + hydro_column_ad
        hydro_column_ad = 0.0_JPRB

        do itype = 1, coef_scatt % nhydro
          if (opts_scatt%hydro_cfrac_tlad) then
            scatt_aux_ad % hydro (iprofad,:,itype) = scatt_aux_ad % hydro (iprofad,:,itype) + &
              hydro_weights_ad(:,itype) * scatt_aux % dz (iprof,:)
          endif

          scatt_aux_ad % dz (iprofad,:) = scatt_aux_ad % dz (iprofad,:) + &
            hydro_weights_ad(:,itype) * hydro_precf (iprof,:,itype)
        enddo

        hydro_weights_ad(:,:) = 0.0_JPRB

      endif

    endif
  enddo

  do ilayer = 1,nlevels
    do iprofad = 1, nprofilesad

      if (adk == adk_adjoint) then
        iprof = iprofad
      else if (adk == adk_k) then
        iprof = chanprof(iprofad) % prof
      endif

      mmr_to_density_ad = 0.0_JPRB

      do itype = 1, coef_scatt % nhydro

        if( cld_profiles (iprof) % flux_conversion (itype) == 0 ) then

          !* Condensate from g/g to g/m^3
          mmr_to_density_ad = mmr_to_density_ad + &
            & hydro_scale (iprof,ilayer,itype) * scatt_aux_ad % hydro (iprofad,ilayer,itype)
          scatt_aux_ad % hydro (iprofad,ilayer,itype) = &
            & scatt_aux_ad % hydro (iprofad,ilayer,itype) * mmr_to_density(iprof,ilayer)

        else

          !* Rates from kg/m^2/s to g/m^3
          if( cld_profiles (iprof) % flux_conversion (itype) == 1) then
            rho  = rho_rain
            conv = conv_rain
          else
            rho  = rho_snow
            conv = conv_sp
          endif

          if (hydro_scale (iprof,ilayer,itype) > 0.0_JPRB) then
            scatt_aux_ad % hydro (iprofad,ilayer,itype) = scatt_aux_ad % hydro (iprofad,ilayer,itype) &
              & * (conv(2)) * (hydro_scale (iprof,ilayer,itype)**(conv(2) - 1.0_JPRB)) &
              & * (conv(1))**(conv(2))
          else
            scatt_aux_ad % hydro (iprofad,ilayer,itype) = 0.0_JPRB
          endif

          scatt_aux_ad % hydro (iprofad,ilayer,itype) = scatt_aux_ad % hydro (iprofad,ilayer,itype) / rho

        endif

      enddo

      presfad (iprofad,ilayer) = presfad (iprofad,ilayer) &
                              & + de2mr / profiles (iprof) % t (ilayer) * mmr_to_density_ad
      profiles_ad (iprofad) % t (ilayer) = profiles_ad (iprofad) % t (ilayer) &
                                        & - presf (iprof,ilayer) * de2mr &
                                        & / (profiles (iprof) % t (ilayer) * profiles (iprof) % t (ilayer)) &
                                        & * mmr_to_density_ad

      mmr_to_density_ad = 0.0_JPRB

    enddo
  enddo

  !* Initialise cloud and rain properties of the cloudy/rainy column
  do ilayer=1,nlevels
    do iprofad = 1, nprofilesad 

      if (adk == adk_adjoint) then
        iprof = iprofad  
      else if (adk == adk_k) then
        iprof = chanprof(iprofad) % prof  
      endif

      cld_profiles_ad (iprofad) % hydro (ilayer,:) = cld_profiles_ad (iprofad) % hydro (ilayer,:) &
                                            & + scatt_aux_ad % hydro (iprofad,ilayer,:)

      scatt_aux_ad % hydro (iprofad,ilayer,:) = 0.0_JPRB

    enddo
  enddo

  if (lhook) call dr_hook('RTTOV_HYDRO_AD',1_jpim,zhook_handle)

end subroutine rttov_hydro_ad
