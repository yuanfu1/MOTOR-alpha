!
Subroutine rttov_hydro_tl (&
      & ccthres,           &! in
      & nlevels,           &! in
      & nprofiles,         &! in
      & usercfrac,         &! in
      & lreflectivity,     &! in
      & presf,             &! in
      & presftl,           &! in
      & profiles,          &! in  
      & profiles_tl,       &! in  
      & cld_profiles,      &! in 
      & cld_profiles_tl,   &! in 
      & opts_scatt,        &! in
      & coef_scatt,        &! in
      & scatt_aux,         &! inout 
      & scatt_aux_tl)       ! inout 

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
       & rttov_profile_cloud

  Use parkind1, Only : jpim, jplm, jprb

!INTF_OFF
  Use rttov_const, Only:      &
       & rgc,                 &
       & mair,                &
       & rho_rain,            &
       & rho_snow,            &
       & conv_rain,           &
       & conv_sp,             &
       & min_cfrac_radar

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
!INTF_ON
  Implicit None
  
!* Subroutine arguments:
  Real    (Kind=jprb), Intent (in) :: ccthres
  Integer (Kind=jpim), Intent (in) :: nlevels   ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles ! Number of profiles
  Logical (Kind=jplm), Intent (in) :: usercfrac               ! User has supplied cloud fraction
  Logical (kind=jplm),              Intent (in)    :: lreflectivity            ! radar simulator active?
  Real (Kind=jprb),                 Intent (in)    :: presf(nprofiles,nlevels) ! Pressure levels [hPa]
  Real (Kind=jprb),                 Intent (in)    :: presftl(nprofiles,nlevels) ! Pressure levels [hPa]
  Type (rttov_profile),             Intent (in)    :: profiles (nprofiles)     ! Atmospheric profiles
  Type (rttov_profile),             Intent (in)    :: profiles_tl (nprofiles)  ! Atmospheric profiles
  Type (rttov_options_scatt),       Intent (in)    :: opts_scatt               ! RTTOV_SCATT options
  Type (rttov_scatt_coef),          Intent (in)    :: coef_scatt               ! RTTOV_SCATT Coefficients
  Type (rttov_profile_cloud),       Intent (in)    :: cld_profiles (nprofiles) ! Cloud profiles
  Type (rttov_profile_cloud),       Intent (in)    :: cld_profiles_tl (nprofiles) ! Cloud profiles TL
  Type (rttov_profile_scatt_aux),   Intent (inout) :: scatt_aux                ! Auxiliary profile variables
  Type (rttov_profile_scatt_aux),   Intent (inout) :: scatt_aux_tl             ! Auxiliary profile variables
!INTF_END

!* Local variables
  Integer (Kind=jpim) :: ilayer, iprof, itype
  Real    (Kind=jprb) :: rho, conv(2)
  Real    (Kind=jprb) :: de2mr, mmr_to_density, mmr_to_density_tl
  Real    (Kind=jprb) :: hydro_frac_weighted(nlevels,coef_scatt % nhydro)    ! Hydrometeor fraction times amount
  Real    (Kind=jprb) :: hydro_weights(nlevels,coef_scatt % nhydro)          ! Hydrometeor amounts [g m^-2]
  Real    (Kind=jprb) :: hydro_column                                        ! Total hydrometeor amounts [g m^-2]
  Real    (Kind=jprb) :: hydro_frac_weighted_tl(nlevels,coef_scatt % nhydro) ! Hydrometeor fraction times amount
  Real    (Kind=jprb) :: hydro_weights_tl(nlevels,coef_scatt % nhydro)       ! Hydrometeor amounts [g m^-2]
  Real    (Kind=jprb) :: hydro_column_tl                                     ! Total hydrometeor amounts [g m^-2]
  Real    (Kind=jprb), pointer :: cfrac_column(:), cfrac_column_tl(:)

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_HYDRO_TL',0_jpim,zhook_handle)

  de2mr =  1.0E+02_JPRB * mair / rgc

  !* Initialise cloud and rain properties of the cloudy/rainy column
  do ilayer=1,nlevels
    do iprof = 1, nprofiles  
      scatt_aux    % hydro (iprof,ilayer,:) = cld_profiles    (iprof) % hydro (ilayer,:)
      scatt_aux_tl % hydro (iprof,ilayer,:) = cld_profiles_tl (iprof) % hydro (ilayer,:)
    enddo
  enddo

  !* Change units
  do ilayer = 1, nlevels
    do iprof = 1, nprofiles

      mmr_to_density    = presf   (iprof,ilayer) * de2mr / profiles    (iprof) % t (ilayer)
      mmr_to_density_tl = presftl (iprof,ilayer) * de2mr / profiles    (iprof) % t (ilayer) &
                      & - presf   (iprof,ilayer) * de2mr * profiles_tl (iprof) % t (ilayer) &
                      & / (profiles (iprof) % t (ilayer) * profiles (iprof) % t (ilayer))

      do itype = 1, coef_scatt % nhydro
        if( cld_profiles (iprof) % flux_conversion (itype) == 0 ) then

          !* Condensate from g/g to g/m^3
          scatt_aux_tl % hydro (iprof,ilayer,itype) = &
            &   scatt_aux_tl % hydro (iprof,ilayer,itype) * mmr_to_density &
            & + scatt_aux    % hydro (iprof,ilayer,itype) * mmr_to_density_tl
          scatt_aux % hydro (iprof,ilayer,itype) = scatt_aux % hydro (iprof,ilayer,itype) * mmr_to_density

        else

          !* Rates from kg/m^2/s to g/m^3
          if( cld_profiles (iprof) % flux_conversion (itype) == 1) then
            rho  = rho_rain
            conv = conv_rain
          else
            rho  = rho_snow
            conv = conv_sp
          endif

          scatt_aux_tl % hydro (iprof,ilayer,itype) = scatt_aux_tl % hydro (iprof,ilayer,itype) / rho
          scatt_aux    % hydro (iprof,ilayer,itype) = scatt_aux    % hydro (iprof,ilayer,itype) / rho

          if (scatt_aux % hydro (iprof,ilayer,itype) > 0.0_JPRB) then
            scatt_aux_tl % hydro (iprof,ilayer,itype) =  scatt_aux_tl % hydro (iprof,ilayer,itype) &
              & * (conv(2)) * (scatt_aux % hydro (iprof,ilayer,itype) ** (conv(2) - 1.0_JPRB)) &
              & * (conv(1))**(conv(2))
            scatt_aux % hydro (iprof,ilayer,itype) = (scatt_aux % hydro (iprof,ilayer,itype) * conv(1))**(conv(2))
          else
            scatt_aux_tl % hydro (iprof,ilayer,itype) = 0.0_JPRB
          endif

        endif
      enddo
    enddo
  enddo

  scatt_aux    % cfrac (:) = 0.0_JPRB
  scatt_aux_tl % cfrac (:) = 0.0_JPRB

  do iprof = 1, nprofiles

    if (lreflectivity) then

      ! Cfrac has no meaning for radar simulator, but set it to 1 to pass checks elsewhere
      scatt_aux    % cfrac (iprof) = 1.0_JPRB
      scatt_aux_tl % cfrac (iprof) = 0.0_JPRB

      ! Each hydrometeor is weighted by its own subgrid fraction
      do itype = 1, coef_scatt % nhydro

        if (cld_profiles (iprof) % nhydro_frac == 1) then
          cfrac_column    => cld_profiles    (iprof) % hydro_frac (:,1)
          cfrac_column_tl => cld_profiles_tl (iprof) % hydro_frac (:,1)
        else
          cfrac_column    => cld_profiles    (iprof) % hydro_frac (:,itype)
          cfrac_column_tl => cld_profiles_tl (iprof) % hydro_frac (:,itype)
        endif

        where(cfrac_column > min_cfrac_radar)
          scatt_aux_tl % hydro (iprof,:,itype) = scatt_aux_tl % hydro (iprof,:,itype) / cfrac_column - &
            & cfrac_column_tl * scatt_aux % hydro (iprof,:,itype) / (cfrac_column**2)
          scatt_aux    % hydro (iprof,:,itype) = scatt_aux    % hydro (iprof,:,itype) / cfrac_column
        elsewhere
          scatt_aux    % hydro (iprof,:,itype) = 0.0_JPRB
          scatt_aux_tl % hydro (iprof,:,itype) = 0.0_JPRB
        endwhere

      enddo

    else

      if( usercfrac ) then

        !* User-supplied cloud fraction
        scatt_aux    % cfrac (iprof) = cld_profiles    (iprof) % cfrac
        scatt_aux_tl % cfrac (iprof) = cld_profiles_tl (iprof) % cfrac

      else

        !* Calculate a hydrometeor-weighted average cloudy sky fraction

        !* Partial column of hydrometeors in g/m^2
        do itype = 1, coef_scatt % nhydro
          hydro_weights (:,itype) = scatt_aux % hydro (iprof,:,itype) * scatt_aux % dz (iprof,:)
          if (opts_scatt%hydro_cfrac_tlad) then
            hydro_weights_tl(:,itype) = scatt_aux % hydro (iprof,:,itype) * scatt_aux_tl % dz (iprof,:) + &
              & scatt_aux_tl % hydro (iprof,:,itype) * scatt_aux % dz (iprof,:)
          else
            hydro_weights_tl(:,itype) = scatt_aux % hydro (iprof,:,itype) * scatt_aux_tl % dz (iprof,:)
          endif
        enddo

        hydro_column = sum(hydro_weights)
        hydro_column_tl = sum(hydro_weights_tl)

        !* Weighted mean cloud fraction
        if (hydro_column > 1e-10_JPRB) then

          if (cld_profiles (iprof) % nhydro_frac == 1) then

            ! One cloud fraction applies to all hydrometeors
            do itype = 1, coef_scatt % nhydro
              hydro_frac_weighted    (:,itype) = cld_profiles (iprof) % hydro_frac (:,1) * hydro_weights (:,itype)
              hydro_frac_weighted_tl (:,itype) = &
                cld_profiles_tl (iprof) % hydro_frac (:,1) * hydro_weights    (:,itype) + &
                cld_profiles    (iprof) % hydro_frac (:,1) * hydro_weights_tl (:,itype)
            enddo

          else

            ! One hydrometeor fraction per hydrometeor
            hydro_frac_weighted    (:,:) = cld_profiles (iprof) % hydro_frac (:,:) * hydro_weights (:,:)
            hydro_frac_weighted_tl (:,:) = &
              cld_profiles_tl (iprof) % hydro_frac (:,:) * hydro_weights    (:,:) + &
              cld_profiles    (iprof) % hydro_frac (:,:) * hydro_weights_tl (:,:)

          endif

          scatt_aux    % cfrac (iprof) = sum(hydro_frac_weighted)    / hydro_column
          scatt_aux_tl % cfrac (iprof) = sum(hydro_frac_weighted_tl) / hydro_column - &
            hydro_column_tl * sum(hydro_frac_weighted) / ( hydro_column**2 )

          if ( scatt_aux % cfrac (iprof) < 0.0_JPRB ) then
            scatt_aux    % cfrac (iprof) = 0.0_JPRB
            scatt_aux_tl % cfrac (iprof) = 0.0_JPRB
          endif

          if ( scatt_aux % cfrac (iprof) > 1.0_JPRB ) then
            scatt_aux    % cfrac (iprof) = 1.0_JPRB
            scatt_aux_tl % cfrac (iprof) = 0.0_JPRB
          endif

        else
          scatt_aux    % cfrac (iprof) = 0.0_JPRB
          scatt_aux_tl % cfrac (iprof) = 0.0_JPRB
        endif

      endif

      !* Partition all cloud and rain into the cloudy column
      if (scatt_aux % cfrac (iprof) > ccthres) Then

        scatt_aux_tl % hydro (iprof,:,:) = &
          &   scatt_aux_tl % hydro (iprof,:,:) / scatt_aux % cfrac (iprof) &
          & - scatt_aux_tl % cfrac (iprof) * scatt_aux % hydro (iprof,:,:) &
          & / (scatt_aux % cfrac (iprof)**2)

        scatt_aux % hydro (iprof,:,:) = scatt_aux % hydro (iprof,:,:) / scatt_aux % cfrac (iprof)

      else
        scatt_aux    % hydro (iprof,:,:) = 0.0_JPRB
        scatt_aux_tl % hydro (iprof,:,:) = 0.0_JPRB
      endif

    endif

  enddo

  if (lhook) call dr_hook('RTTOV_HYDRO_TL',1_jpim,zhook_handle)

end subroutine rttov_hydro_tl
