Subroutine rttov_iniedd_tl (&
     & lreflectivity, &! in
     & ccthres,       &! in
     & nlevels,       &! in
     & nchannels ,    &! in
     & nprofiles ,    &! in
     & chanprof  ,    &! in
     & angles    ,    &! in
     & coef      ,    &! in
     & scatt_aux,     &! inout
     & scatt_aux_tl)   ! inout 

  ! Description:
  ! to compute variables specific to Eddington approximation to RT
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
  !  1.2       03/2004   Added polarimetry   (R. Saunders)
  !  1.3       08/2004   Polarimetry fixes   (U. O'Keefe)
  !  1.4       11/2004   Clean-up            (P. Bauer)
  !  1.5       11/2007   RTTOV9 version      (A. Geer)
  !  1.6       07/2008   Consistent mininimum SSA (A. Geer)
  !  1.7       07/2008   Clear sky speed-ups (A. Geer)
  !  1.8       03/2010   Optimisation + don't delta scale outside Eddington (A. Geer)
  !  1.9       11/2017   R/T now done with radiances, not Tb (A. Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:

  Use rttov_types, Only :    &
       & rttov_coef           ,&
       & rttov_geometry        ,&
       & rttov_chanprof        ,&
       & rttov_profile_scatt_aux 

  Use parkind1, Only : jpim, jplm, jprb
!INTF_OFF
  Use rttov_const, Only : min_ssa, max_scatt_optical_depth

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
  use rttov_math_mod, only : planck, planck_tl
!INTF_ON
  Implicit None

!* Subroutine arguments:
  Logical (Kind=jplm), Intent (in) :: lreflectivity         ! Computation of radar reflectivity, not TB
  Real    (Kind=jprb), Intent (in) :: ccthres
  Integer (Kind=jpim), Intent (in) :: nlevels               ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles             ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels             ! Number of radiances
  Type(rttov_chanprof), Intent(in) :: chanprof(nchannels)   ! Channel and profile indices

  Type (rttov_geometry),          Intent (in)    :: angles (nprofiles) ! Zenith angles
  Type (rttov_coef),              Intent (in)    :: coef               ! RTTOV Coefficients
  Type (rttov_profile_scatt_aux), Intent (inout) :: scatt_aux          ! Auxiliary profile variables
  Type (rttov_profile_scatt_aux), Intent (inout) :: scatt_aux_tl       ! Auxiliary profile variables

!INTF_END

!* Local variables
  Integer(Kind=jpim) :: ilayer, iprof, ichan, ichanid

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_INIEDD_TL',0_jpim,zhook_handle)

  scatt_aux_tl % delta   = 0.0_JPRB
  scatt_aux_tl % lambda  = 0.0_JPRB
  scatt_aux_tl % h       = 0.0_JPRB
  scatt_aux_tl % tau     = 0.0_JPRB
  scatt_aux_tl % int_tau = 0.0_JPRB

  scatt_aux    % delta   = 0.0_JPRB
  scatt_aux    % lambda  = 0.0_JPRB
  scatt_aux    % h       = 0.0_JPRB
  scatt_aux    % tau     = 1.0_JPRB
  scatt_aux    % int_tau = 1.0_JPRB
  
  scatt_aux    % mclayer = 0

  ! Planck function top and bottom of layer, its gradient, top & bottom boundaries
  do ichan = 1, nchannels
    iprof   = chanprof(ichan)%prof
    ichanid = chanprof(ichan)%chan

    call planck   (coef%planck1(ichanid), coef%planck2(ichanid), scatt_aux%tbd(ichan,2:nlevels+1), &
                  & scatt_aux%b0(ichan,:))
    call planck_tl(coef%planck1(ichanid), coef%planck2(ichanid), &
                  & scatt_aux%tbd(ichan,2:nlevels+1), scatt_aux_tl%tbd(ichan,2:nlevels+1), &
                  & scatt_aux%b0(ichan,:), scatt_aux_tl%b0(ichan,:))
    call planck   (coef%planck1(ichanid), coef%planck2(ichanid), scatt_aux%tbd(ichan,1:nlevels  ), &
                  & scatt_aux%bn(ichan,:))
    call planck_tl(coef%planck1(ichanid), coef%planck2(ichanid), scatt_aux%tbd(ichan,1:nlevels  ), &
                  & scatt_aux_tl%tbd(ichan,1:nlevels  ), &
                  & scatt_aux%bn(ichan,:), scatt_aux_tl%bn(ichan,:))

    scatt_aux_tl%b1(ichan,:) = (scatt_aux_tl%bn(ichan,:)  -  scatt_aux_tl%b0(ichan,:)) &
                            & /  scatt_aux   %dz(iprof,:)  - (scatt_aux   %bn(ichan,:)  &
                            & -  scatt_aux   %b0(ichan,:)) *  scatt_aux_tl%dz(iprof,:)  &
                            & /  scatt_aux   %dz(iprof,:)  /  scatt_aux   %dz(iprof,:)
    scatt_aux%b1(ichan,:) = (scatt_aux%bn(ichan,:) - scatt_aux%b0(ichan,:)) / scatt_aux%dz(iprof,:)

    call planck   (coef%planck1(ichanid), coef%planck2(ichanid), scatt_aux%tsfc(ichan), &
                  & scatt_aux%bsfc(ichan))
    call planck_tl(coef%planck1(ichanid), coef%planck2(ichanid), scatt_aux%tsfc(ichan), scatt_aux_tl%tsfc(ichan), &
                  & scatt_aux%bsfc(ichan), scatt_aux_tl%bsfc(ichan))
    call planck   (coef%planck1(ichanid), coef%planck2(ichanid), scatt_aux%tcosmic(ichan), &
                  & scatt_aux%btop(ichan))
    call planck_tl(coef%planck1(ichanid), coef%planck2(ichanid), scatt_aux%tcosmic(ichan), scatt_aux_tl%tcosmic(ichan), &
                  & scatt_aux%btop(ichan), scatt_aux_tl%btop(ichan))

  end do

  do ichan = 1, nchannels
    !* Cloud top level index
    scatt_aux % mclayer (ichan) = nlevels+1 
    do ilayer = 1, nlevels 
      if (scatt_aux % ssa (ichan,ilayer) > min_ssa ) then 
        scatt_aux % mclayer (ichan) = ilayer
        exit
      endif
    end do
    if (scatt_aux % mclayer (ichan) > nlevels-2 .and. scatt_aux % mclayer (ichan) /= nlevels+1) &
      & scatt_aux % mclayer (ichan) = nlevels-2 !* DGBF imposes minimum number of layers in Eddington
  end do 

  !* Delta-scaling
  do ilayer = 1, nlevels
    do ichan = 1, nchannels
      iprof = chanprof(ichan)%prof
      if ( scatt_aux % cfrac(iprof) > ccthres ) then 

        if (ilayer >= scatt_aux % mclayer (ichan)) then
          if (.not. lreflectivity) then
            scatt_aux_tl % ext  (ichan,ilayer) = (1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer) &
                                    & * scatt_aux % asm (ichan,ilayer)) * scatt_aux_tl % ext  (ichan,ilayer) &
                                    & + scatt_aux % ext (ichan,ilayer) * scatt_aux % asm (ichan,ilayer) &
                                    & * (-1.0_JPRB * scatt_aux_tl % ssa  (ichan,ilayer) * scatt_aux % asm (ichan,ilayer) &
                                    & -2.0_JPRB * scatt_aux % ssa (ichan,ilayer) * scatt_aux_tl % asm  (ichan,ilayer))
            scatt_aux    % ext  (ichan,ilayer) = (1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer) &
                                    & * scatt_aux % asm (ichan,ilayer)) * scatt_aux % ext (ichan,ilayer)
            scatt_aux_tl % ssa  (ichan,ilayer) = (1.0_JPRB /  (1.0_JPRB - scatt_aux % asm (ichan,ilayer) &
                                    & * scatt_aux % asm (ichan,ilayer) &
                                    & * scatt_aux % ssa (ichan,ilayer))**2.0_JPRB ) * ( (1.0_JPRB - scatt_aux % asm (ichan,ilayer) &
                                    & * scatt_aux % asm (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer))         &
                                    & * (( 1.0_JPRB - scatt_aux % asm (ichan,ilayer) * scatt_aux % asm (ichan,ilayer)) *&
                                    & scatt_aux_tl % ssa (ichan,ilayer) &
                                    & +  (-2.0_JPRB * scatt_aux % asm (ichan,ilayer) * scatt_aux_tl % asm (ichan,ilayer) * &
                                    &scatt_aux % ssa (ichan,ilayer)) ) &
                                    & - (1.0_JPRB - scatt_aux % asm (ichan,ilayer) * scatt_aux % asm (ichan,ilayer)) * &
                                    &scatt_aux % ssa (ichan,ilayer) &
                                    & * (  (-1.0_JPRB * scatt_aux % asm (ichan,ilayer) * scatt_aux % asm (ichan,ilayer) * &
                                    &scatt_aux_tl % ssa  (ichan,ilayer)) &
                                    & +  (-2.0_JPRB * scatt_aux % asm (ichan,ilayer) * scatt_aux_tl % asm  (ichan,ilayer) * &
                                    &scatt_aux % ssa (ichan,ilayer)) ) )
            scatt_aux    % ssa  (ichan,ilayer) = (1.0_JPRB - scatt_aux % asm (ichan,ilayer) * scatt_aux % asm (ichan,ilayer)) &
                                    & * scatt_aux % ssa (ichan,ilayer) / (1.0_JPRB - scatt_aux % asm (ichan,ilayer)  &
                                    & * scatt_aux % asm (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer))
            scatt_aux_tl % asm  (ichan,ilayer) = scatt_aux_tl % asm (ichan,ilayer) / (1.0_JPRB + scatt_aux % asm (ichan,ilayer)) &
                                    & - scatt_aux % asm (ichan,ilayer) * scatt_aux_tl % asm (ichan,ilayer) / &
                                    &(1.0_JPRB + scatt_aux % asm (ichan,ilayer)) &
                                    & / (1.0_JPRB + scatt_aux % asm (ichan,ilayer))
            scatt_aux    % asm  (ichan,ilayer) = scatt_aux % asm (ichan,ilayer) / (1.0_JPRB + scatt_aux % asm (ichan,ilayer))
          endif

          scatt_aux_tl % lambda (ichan,ilayer) = (1.0_JPRB / ( 2._JPRB * sqrt (3.0_JPRB * scatt_aux % ext (ichan,ilayer) *&
                                     &scatt_aux % ext (ichan,ilayer)  &
                                     & * (1.0_JPRB - scatt_aux % ssa (ichan,ilayer)) &
                                     & * (1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer))))) &
                                     & * ( 6.0_JPRB * scatt_aux % ext (ichan,ilayer) * scatt_aux_tl % ext (ichan,ilayer) &
                                     & * (1.0_JPRB - scatt_aux % ssa (ichan,ilayer)) * &
                                     &(1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer)) &
                                     & - 3.0_JPRB * scatt_aux % ext (ichan,ilayer) * scatt_aux % ext(ichan,ilayer) &
                                     & * scatt_aux_tl % ssa (ichan,ilayer) *  &
                                     &(1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer)) &
                                     & - 3.0_JPRB * scatt_aux % ext (ichan,ilayer) * scatt_aux % ext (ichan,ilayer) * &
                                     &(1.0_JPRB - scatt_aux % ssa (ichan,ilayer)) &
                                     & * (scatt_aux % ssa (ichan,ilayer) * scatt_aux_tl % asm (ichan,ilayer)+ &
                                     &scatt_aux_tl % ssa (ichan,ilayer) &
                                     & * scatt_aux % asm (ichan,ilayer)) ) 

          scatt_aux    % lambda (ichan,ilayer) = sqrt (3.0_JPRB * scatt_aux % ext (ichan,ilayer) * scatt_aux % ext (ichan,ilayer) &
                                     & * (1.0_JPRB - scatt_aux % ssa (ichan,ilayer)) &
                                     & * (1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer))) 

          scatt_aux_tl % h      (ichan,ilayer) = 1.5_JPRB * scatt_aux_tl % ext (ichan,ilayer) &
                                     & * (1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer)) &
                                     & -1.5_JPRB * scatt_aux % ext (ichan,ilayer) &
                                     & * (scatt_aux_tl % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer) + &
                                     &scatt_aux % ssa (ichan,ilayer) * scatt_aux_tl % asm (ichan,ilayer)) 
          scatt_aux    % h      (ichan,ilayer) = 1.5_JPRB * scatt_aux % ext (ichan,ilayer) &
                                     & * (1.0_JPRB - scatt_aux % ssa (ichan,ilayer) * scatt_aux % asm (ichan,ilayer)) 

          if (scatt_aux    % h (ichan,ilayer) < 0.00001_JPRB) then
              scatt_aux_tl % h (ichan,ilayer) = 0.0_JPRB
              scatt_aux    % h (ichan,ilayer) = 0.00001_JPRB
          endif

        endif
      endif

      scatt_aux_tl % delta (ichan,ilayer) = (scatt_aux_tl % ext (ichan,ilayer) * scatt_aux    % dz (iprof,ilayer) &
                                & +  scatt_aux    % ext (ichan,ilayer) * scatt_aux_tl % dz (iprof,ilayer)) / &
                                &angles (iprof) % coszen
      scatt_aux    % delta (ichan,ilayer) = (scatt_aux    % ext (ichan,ilayer) * scatt_aux    % dz (iprof,ilayer)) /&
        & angles (iprof) % coszen

      if (scatt_aux    % delta (ichan,ilayer) >= max_scatt_optical_depth) then
          scatt_aux    % delta (ichan,ilayer)  = max_scatt_optical_depth
          scatt_aux_tl % delta (ichan,ilayer)  =  0.0_JPRB
      endif

      scatt_aux    % tau    (ichan,ilayer) =  1.0_JPRB / exp (scatt_aux % delta (ichan,ilayer))
      scatt_aux_tl % tau    (ichan,ilayer) = -1.0_JPRB * scatt_aux_tl % delta (ichan,ilayer) * scatt_aux % tau (ichan,ilayer)

      if( ilayer == 1) then
        scatt_aux    % int_tau (ichan,ilayer) = scatt_aux    % tau (ichan,ilayer)
        scatt_aux_tl % int_tau (ichan,ilayer) = scatt_aux_tl % tau (ichan,ilayer)
      else
        scatt_aux    % int_tau (ichan,ilayer) = scatt_aux    % int_tau (ichan,ilayer-1) &
                                            & * scatt_aux    % tau (ichan,ilayer)
        scatt_aux_tl % int_tau (ichan,ilayer) = scatt_aux    % int_tau (ichan,ilayer-1) &
                                            & * scatt_aux_tl % tau (ichan,ilayer)       &
                                            & + scatt_aux_tl % int_tau (ichan,ilayer-1) &
                                            & * scatt_aux    % tau (ichan,ilayer)
      endif

    end do
  end do

  if (lhook) call dr_hook('RTTOV_INIEDD_TL',1_jpim,zhook_handle)

End subroutine rttov_iniedd_tl
