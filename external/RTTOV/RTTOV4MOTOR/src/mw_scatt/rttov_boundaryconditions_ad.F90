!
Subroutine rttov_boundaryconditions_ad (& 
     & ccthres,       &! in
     & nlevels,       &! in
     & nchannels,     &! in
     & nprofiles,     &! in
     & nprofilesad,   &! in
     & lprofiles,     &! in
     & scatt_aux,     &! in
     & scatt_aux_ad,  &! inout
     & ftop,          &! in
     & ftop_ad,       &! inout
     & dp,            &! out
     & dp_ad,         &! inout
     & dm,            &! out
     & dm_ad)          ! inout 


  ! Description:
  ! to compute boundary conditions for Eddington approximation to RT
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
  !  1.0       09/2002   Initial version     (E. Moreau)
  !  1.1       05/2003   RTTOV7.3 compatible (F. Chevallier)
  !  1.2       03/2004   Added polarimetry   (R. Saunders)
  !  1.3       08/2004   Polarimetry fixes   (U. O'Keefe)
  !  1.4       11/2004   Clean-up            (P. Bauer)
  !  1.5       11/2007   RTTOV9 version      (A. Geer)
  !  1.6       07/2008   Clear sky speed-ups (A. Geer)
  !  1.7       03/2010   Optimisation        (A. Geer)
  !  1.8       11/2017   R/T now done with radiances, not Tb (A. Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarationsta:
  ! Modules used:
  ! Imported Type Definitions:

  Use rttov_types, Only :    &
       & rttov_profile_scatt_aux 

  Use parkind1, Only : jpim     ,jprb
!INTF_OFF
  Use rttov_const, Only: adk_adjoint, adk_k

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
#ifndef _RTTOV_ECMWF
  USE rttov_lapack_mod, ONLY : dgbtrf, dgbtrs
#endif
!INTF_ON
  Implicit none

!* Subroutine arguments:
  Real    (Kind=jprb), Intent (in) :: ccthres
  Integer (Kind=jpim), Intent (in) :: nlevels               ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles             ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nprofilesad           ! Number of profiles in adjoint variables
  Integer (Kind=jpim), Intent (in) :: nchannels             ! Number of radiances
  Integer (Kind=jpim), Intent (in) :: lprofiles (nchannels) ! Profile indices

  Type (rttov_profile_scatt_aux), Intent    (in) :: scatt_aux               ! Auxiliary profile variables for RTTOV_SCATT
  Type (rttov_profile_scatt_aux), Intent (inout) :: scatt_aux_ad            ! Auxiliary profile variables for RTTOV_SCATT

  Real (Kind=jprb), Intent    (in), dimension (nchannels)            :: ftop
  Real (Kind=jprb), Intent (inout), dimension (nchannels)            :: ftop_ad
  Real (Kind=jprb), Intent   (out), dimension (nchannels,nlevels) :: dp   , dm
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: dp_ad, dm_ad

!INTF_END

!* Local variables
  Real    (Kind=jprb), dimension (nchannels,nlevels) :: lh_p   , lh_m   , bh
  Real    (Kind=jprb), dimension (nchannels,nlevels) :: lh_p_ad, lh_m_ad, bh_ad
  Real    (Kind=jprb), allocatable :: b (:)
  Real    (Kind=jprb), allocatable :: dx_ad (:)
  Real    (Kind=jprb), allocatable :: ab_ad (:, :)
  Real    (Kind=jprb)              :: ztmp, ztmp_ad
  Integer (Kind=jpim)              :: ilayer, jlayer, klayer, ilin, icol, iband, uband
  Integer (Kind=jpim)              :: nmaxdim, iprof, ichan, jj, ii, mcly
  Integer (Kind=jpim)              :: iprofad, adk

!* Lapack/ESSL
#ifdef _RTTOV_ECMWF 
  Real    (Kind=jprb), allocatable :: ab (:, :), tab (:, :)
  Real    (Kind=jprb), allocatable :: dx (:)
  Real    (Kind=jprb), allocatable :: b_ad (:)
  Integer (Kind=jpim), allocatable :: ipiv(:)
  Integer (Kind=jpim)              :: ndim, kl, ku, ldab, info, nrhs
#else
  Double Precision,    allocatable :: ab (:, :), tab (:, :)
  Double Precision,    allocatable :: dx (:)
  Double Precision,    allocatable :: b_ad (:)
  Integer,             allocatable :: ipiv(:)
  Integer                          :: ndim, kl, ku, ldab, info, nrhs
#endif
  Character (len=1)                :: trans

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_BOUNDARYCONDITIONS_AD',0_jpim,zhook_handle)        

  if (nprofilesad == nprofiles) then 
    adk = adk_adjoint   ! Adjoint mode
  else if (nprofilesad == nchannels) then
    adk = adk_k         ! K mode
  endif 

  !* Indices for band matrix representation for lapack/essl
  kl = 2
  ku = 2
  iband = kl + ku 
  uband = 2 * kl + ku + 1
  ldab = uband 
  trans = 'N'
  nrhs  = 1
  info  = 0

  nmaxdim = 2 * (nlevels - minval(scatt_aux % mclayer (:)) + 1)

  allocate (b     (nmaxdim     ))
  allocate (b_ad  (nmaxdim     ))
  allocate (dx    (nmaxdim     ))
  allocate (dx_ad (nmaxdim     ))
  allocate (tab   (ldab,nmaxdim))
  allocate (ab    (ldab,nmaxdim))
  allocate (ab_ad (ldab,nmaxdim))
  allocate (ipiv  (nmaxdim     ))

!* FORWARD PART
!* Reset      
  dp (:,:) = 0.0_JPRB
  dm (:,:) = 0.0_JPRB

  lh_p_ad = 0.0_JPRB
  lh_m_ad = 0.0_JPRB
  bh_ad   = 0.0_JPRB

!* Channels * Profiles  
  do ilayer=1,nlevels
    do ichan = 1, nchannels
      iprof = lprofiles (ichan)
      if( scatt_aux % cfrac (iprof) > ccthres .and. ilayer >= scatt_aux % mclayer(ichan)) then 

        bh   (ichan,ilayer) = scatt_aux % b1 (ichan,ilayer) / scatt_aux % h (ichan,ilayer)
        lh_p (ichan,ilayer) = (1.0_JPRB + scatt_aux % lambda (ichan,ilayer) / scatt_aux % h (ichan,ilayer)) 
        lh_m (ichan,ilayer) = (1.0_JPRB - scatt_aux % lambda (ichan,ilayer) / scatt_aux % h (ichan,ilayer))

      endif
    enddo
  enddo
  do ichan = 1, nchannels
    iprof = lprofiles (ichan)
    if( scatt_aux % cfrac (iprof) > ccthres .and. scatt_aux % mclayer (ichan) <= nlevels ) then 

      if (adk == adk_adjoint) then
        iprofad = iprof  
      else if (adk == adk_k) then
        iprofad = ichan  
      endif
 
      mcly = scatt_aux % mclayer (ichan)

      ndim = 2 * (nlevels - mcly + 1)

      do ilayer = 2, ndim - 2, 2
        jlayer = nlevels - ilayer / 2 + 1
        klayer = jlayer - 1

        ilin =  ilayer
        icol = (ilayer - 1)

        ztmp = exp (scatt_aux % lambda (ichan,jlayer) * scatt_aux % dz (iprof,jlayer))

        !* From downward fluxes at i-th interface (@ level=dz for jlayer == level=0 for klayer)
        if(icol > 1) ab (iband+3,icol-1) = 0.0_JPRB
        ab (iband+2,icol  ) =             lh_p (ichan,jlayer) * ztmp
        ab (iband+1,icol+1) =             lh_m (ichan,jlayer) / ztmp
        ab (iband  ,icol+2) = -1.0_JPRB * lh_p (ichan,klayer) 
        ab (iband-1,icol+3) = -1.0_JPRB * lh_m (ichan,klayer) 
        
        b (ilin  ) = bh (ichan,klayer) - bh (ichan,jlayer)

        !* From upward fluxes at i-th interface (@ level=dz for jlayer == level=0 for klayer)
        ab (iband+3,icol  ) =             lh_m (ichan,jlayer) * ztmp
        ab (iband+2,icol+1) =             lh_p (ichan,jlayer) / ztmp
        ab (iband+1,icol+2) = -1.0_JPRB * lh_m (ichan,klayer) 
        ab (iband  ,icol+3) = -1.0_JPRB * lh_p (ichan,klayer) 
        if(icol < ndim-3) ab (iband-1,icol+4) = 0.0_JPRB        

        b (ilin+1) = bh (ichan,jlayer) - bh (ichan,klayer)
      end do

      !* From boundary conditions at bottom of the atmosphere with r_sfc=1-e_sfc
      ztmp = (2.0_JPRB - scatt_aux % ems_bnd (ichan)) &
        &   * scatt_aux % lambda (ichan,nlevels) / scatt_aux % h (ichan,nlevels)

      ab (iband+1,1) = scatt_aux % ems_bnd(ichan) - ztmp
      ab (iband  ,2) = scatt_aux % ems_bnd(ichan) + ztmp
      ab (iband-1,3) = 0.0_JPRB

      b (1)   = scatt_aux % ems_bnd(ichan) &
           & * (scatt_aux % bsfc (ichan) - scatt_aux % b0 (ichan,nlevels)) &
           & + (2.0_JPRB - scatt_aux % ems_bnd (ichan)) * bh (ichan,nlevels)  

      !* From boundary conditions at top of the atmosphere 
      ztmp = exp (scatt_aux % lambda (ichan,mcly) * scatt_aux % dz (iprof,mcly))
     
      ab (iband+3,ndim-2) = 0.0_JPRB
      ab (iband+2,ndim-1) = lh_p (ichan,mcly) * ztmp
      ab (iband+1,ndim  ) = lh_m (ichan,mcly) / ztmp
     
      b (ndim) = ftop (ichan) - scatt_aux % bn (ichan,mcly) - bh (ichan,mcly)

      ! Store a transposed copy of ab for the adjoint part
      do jj = 1, ndim
        do ii = max(1_jpim,jj-ku), min(ndim,jj+kl)
          tab (kl+ku+ii-jj+1,jj) = ab (kl+ku+jj-ii+1,ii)
        end do
      end do

      !* Solve equations A * DX = B          
      call dgbtrf (ndim, ndim, kl, ku, ab, ldab, ipiv, info)                  
      
      dx (1:ndim) = b (1:ndim)
     
      call dgbtrs (trans, ndim, kl, ku, nrhs, ab, ldab, ipiv, dx, ndim, info)

      !* Decompose D+ and D-
      do ilayer = 2, ndim, 2
        jlayer = nlevels - ilayer / 2 + 1
        
        dp (ichan,jlayer) = dx (ilayer-1)
        dm (ichan,jlayer) = dx (ilayer  )
      end do

!* ADJOINT PART
      dx_ad(1:ndim)              = 0.0_JPRB

      !* Decompose D+ and D-
      do ilayer = 2, ndim, 2
        jlayer = nlevels - ilayer / 2 + 1
        
        dx_ad (ilayer  ) = dx_ad (ilayer  ) + dm_ad (ichan,jlayer)
        dx_ad (ilayer-1) = dx_ad (ilayer-1) + dp_ad (ichan,jlayer)
        
        dm_ad (ichan,jlayer) = 0.0_JPRB
        dp_ad (ichan,jlayer) = 0.0_JPRB
      end do
     
      !* Solve equations A * DX = B     
      call dgbtrf (ndim, ndim, kl, ku, tab, ldab, ipiv, info)               

      b_ad (1:ndim) = dx_ad (1:ndim) 

      call dgbtrs (trans, ndim, kl, ku, nrhs, tab, ldab, ipiv, b_ad, ndim, info) 

!     Following is the band matrix equivalent of:
!      do ilayer = 1, ndim
!        do jlayer = 1, ndim
!          a_ad (ilayer,jlayer) = -1.0_JPRB * b_ad (ilayer) * dx (jlayer)
!        enddo
!      enddo
      do jj = 1, ndim
        do ii = max(1_jpim,jj-ku), min(ndim,jj+kl)
          ab_ad (kl+ku+ii-jj+1,jj) = -1.0_JPRB * b_ad (ii) * dx (jj)
        end do
      end do

      !* From boundary conditions at top of the atmosphere 
      ftop_ad (ichan) = ftop_ad (ichan) + b_ad (ndim)
      scatt_aux_ad % bn (ichan,mcly) = scatt_aux_ad % bn (ichan,mcly) - b_ad (ndim)
      bh_ad (ichan,mcly) = bh_ad (ichan,mcly) - b_ad (ndim)
      b_ad (ndim) = 0.0_JPRB

      ztmp_ad = -1.0_JPRB * ab_ad(iband+1,ndim) * lh_m (ichan,mcly) / ztmp / ztmp
      lh_m_ad (ichan,mcly) = lh_m_ad (ichan,mcly) + ab_ad (iband+1,ndim) / ztmp
      ab_ad (iband+1,ndim) = 0.0_JPRB
 
      ztmp_ad = ztmp_ad + ab_ad (iband+2,ndim-1) * lh_p (ichan,mcly) 
      lh_p_ad (ichan,mcly) = lh_p_ad (ichan,mcly) + ab_ad (iband+2,ndim-1) * ztmp
      ab_ad (iband+2,ndim-1) = 0.0_JPRB

      scatt_aux_ad % lambda (ichan,mcly)   = scatt_aux_ad % lambda (ichan,mcly)   &
       & + ztmp_ad * scatt_aux % dz (iprof,mcly) * ztmp
      scatt_aux_ad % dz     (iprofad,mcly) = scatt_aux_ad % dz     (iprofad,mcly) &
       & + ztmp_ad * scatt_aux % lambda (ichan,mcly) * ztmp
      ztmp_ad = 0.0_JPRB

      !* From boundary conditions at bottom of the atmosphere with r_sfc=1-e_sfc
      scatt_aux_ad % ems_bnd (ichan) = scatt_aux_ad % ems_bnd (ichan) + b_ad (1) &
       & * (scatt_aux % bsfc (ichan) - scatt_aux % b0 (ichan,nlevels)) &
       & - b_ad (1) * bh (ichan,nlevels) 
      scatt_aux_ad % bsfc (ichan) = scatt_aux_ad % bsfc (ichan) &
        & + b_ad (1) * scatt_aux % ems_bnd (ichan)
      scatt_aux_ad % b0 (ichan,nlevels) = scatt_aux_ad % b0 (ichan,nlevels) &
        & - b_ad (1) * scatt_aux % ems_bnd (ichan)
      bh_ad (ichan,nlevels) = bh_ad (ichan,nlevels) + b_ad (1) &
        & * (2.0_JPRB - scatt_aux % ems_bnd (ichan))
      b_ad (1) = 0.0_JPRB

      ztmp = (2.0_JPRB - scatt_aux % ems_bnd (ichan)) &
       &   * scatt_aux % lambda (ichan,nlevels) / scatt_aux % h (ichan,nlevels)

      scatt_aux_ad % ems_bnd (ichan) = scatt_aux_ad % ems_bnd (ichan) + ab_ad (iband,2)
      ztmp_ad = ab_ad (iband,2)
      ab_ad (iband,2) = 0.0_JPRB

      scatt_aux_ad % ems_bnd (ichan) = scatt_aux_ad % ems_bnd (ichan) + ab_ad (iband+1,1)
      ztmp_ad = ztmp_ad - ab_ad (iband+1,1)
      ab_ad (iband+1,1) = 0.0_JPRB

      scatt_aux_ad % ems_bnd (ichan) = scatt_aux_ad % ems_bnd (ichan) &
       & - ztmp_ad * scatt_aux % lambda (ichan,nlevels) / scatt_aux % h (ichan,nlevels) 
      scatt_aux_ad % lambda (ichan,nlevels) = scatt_aux_ad % lambda (ichan,nlevels) &
       & + ztmp_ad * (2.0_JPRB - scatt_aux % ems_bnd(ichan)) / scatt_aux % h (ichan,nlevels) 
      scatt_aux_ad % h (ichan,nlevels) = scatt_aux_ad % h (ichan,nlevels) &
       & - ztmp_ad * (2.0_JPRB - scatt_aux % ems_bnd (ichan)) &
       & * scatt_aux % lambda (ichan,nlevels) &
       & / scatt_aux % h (ichan,nlevels) / scatt_aux % h (ichan,nlevels) 
      ztmp_ad = 0.0_JPRB
     
      do ilayer = ndim - 2, 2, -2
        jlayer = nlevels - ilayer / 2 + 1
        klayer = jlayer - 1

        ilin =  ilayer
        icol = (ilayer - 1)

        ztmp = exp (scatt_aux % lambda (ichan,jlayer) * scatt_aux % dz (iprof,jlayer))

        !* From upward fluxes at i-th interface (@ level=dz for jlayer == level=0 for klayer)
        bh_ad (ichan,jlayer) = bh_ad (ichan,jlayer) + b_ad (ilin+1)
        bh_ad (ichan,klayer) = bh_ad (ichan,klayer) - b_ad (ilin+1)
        b_ad  (ilin+1) = 0.0_JPRB

        lh_p_ad (ichan,klayer) = lh_p_ad (ichan,klayer) - ab_ad (iband,icol+3)
        ab_ad (iband,icol+3) = 0.0_JPRB

        lh_m_ad (ichan,klayer) = lh_m_ad (ichan,klayer) - ab_ad (iband+1,icol+2)
        ab_ad (iband+1,icol+2) = 0.0_JPRB

        lh_p_ad (ichan,jlayer) = lh_p_ad (ichan,jlayer) + ab_ad (iband+2,icol+1) / ztmp
        ztmp_ad = -1.0_JPRB * ab_ad (iband+2,icol+1) * lh_p (ichan,jlayer) / ztmp / ztmp
        ab_ad (iband+2,icol+1) = 0.0_JPRB

        lh_m_ad (ichan,jlayer) = lh_m_ad (ichan,jlayer) + ab_ad (iband+3,icol) * ztmp
        ztmp_ad = ztmp_ad + ab_ad (iband+3,icol) * lh_m (ichan,jlayer) 
        ab_ad (iband+3,icol  ) = 0.0_JPRB

        !* From downward fluxes at i-th interface (@ level=dz for jlayer == level=0 for klayer)
        bh_ad (ichan,klayer) = bh_ad (ichan,klayer) + b_ad (ilin  )
        bh_ad (ichan,jlayer) = bh_ad (ichan,jlayer) - b_ad (ilin  )
        b_ad  (ilin  ) = 0.0_JPRB

        lh_m_ad (ichan,klayer) = lh_m_ad (ichan,klayer) - ab_ad (iband-1,icol+3)
        ab_ad (iband-1 ,icol+3) = 0.0_JPRB

        lh_p_ad (ichan,klayer) = lh_p_ad (ichan,klayer) - ab_ad (iband ,icol+2)
        ab_ad (iband ,icol+2) = 00._JPRB

        lh_m_ad (ichan,jlayer) = lh_m_ad (ichan,jlayer) + ab_ad (iband+1 ,icol+1) / ztmp
        ztmp_ad = ztmp_ad - ab_ad (iband+1 ,icol+1) * lh_m (ichan,jlayer) / ztmp / ztmp
        ab_ad (iband+1  ,icol+1) = 0.0_JPRB

        lh_p_ad (ichan,jlayer) = lh_p_ad (ichan,jlayer) + ab_ad (iband+2 ,icol) * ztmp
        ztmp_ad = ztmp_ad + ab_ad (iband+2 ,icol  ) * lh_p (ichan,jlayer) 
        ab_ad (iband+2 ,icol  ) = 0.0_JPRB

        scatt_aux_ad % lambda (ichan,jlayer) = scatt_aux_ad % lambda (ichan,jlayer) + &
         & ztmp_ad * ztmp * scatt_aux % dz     (iprof,jlayer)
        scatt_aux_ad % dz     (iprofad,jlayer) = scatt_aux_ad % dz (iprofad,jlayer) + &
         & ztmp_ad * ztmp * scatt_aux % lambda (ichan,jlayer)
        ztmp_ad = 0._JPRB
      enddo
    endif
  end do 

  do ilayer=1,nlevels
    do ichan = 1, nchannels
      iprof = lprofiles (ichan)
      if( scatt_aux % cfrac (iprof) > ccthres .and. ilayer >= scatt_aux % mclayer(ichan)) then 

        if (adk == adk_adjoint) then
          iprofad = iprof  
        else if (adk == adk_k) then
          iprofad = ichan  
        endif
        scatt_aux_ad % lambda (ichan,ilayer) = scatt_aux_ad % lambda (ichan,ilayer) &
         & - lh_m_ad (ichan,ilayer) / scatt_aux % h      (ichan,ilayer)
        scatt_aux_ad % h      (ichan,ilayer) = scatt_aux_ad % h      (ichan,ilayer) &
         & + lh_m_ad (ichan,ilayer) * scatt_aux % lambda (ichan,ilayer) &
         & / scatt_aux % h (ichan,ilayer) / scatt_aux % h (ichan,ilayer)
        lh_m_ad (ichan,ilayer) = 0.0_JPRB

        scatt_aux_ad % lambda (ichan,ilayer) = scatt_aux_ad % lambda (ichan,ilayer) &
         & + lh_p_ad (ichan,ilayer) / scatt_aux % h      (ichan,ilayer)
        scatt_aux_ad % h      (ichan,ilayer) = scatt_aux_ad % h      (ichan,ilayer) &
         & - lh_p_ad (ichan,ilayer) * scatt_aux % lambda (ichan,ilayer) &
         & / scatt_aux % h (ichan,ilayer) / scatt_aux % h (ichan,ilayer) 
        lh_p_ad (ichan,ilayer) = 0.0_JPRB

        scatt_aux_ad % b1 (ichan,ilayer) = scatt_aux_ad % b1 (ichan,ilayer) &
         & + bh_ad (ichan,ilayer) / scatt_aux % h  (ichan,ilayer)
        scatt_aux_ad % h  (ichan,ilayer)   = scatt_aux_ad % h  (ichan,ilayer)   &
         & - bh_ad (ichan,ilayer) * scatt_aux % b1 (ichan,ilayer) &
         & / scatt_aux % h (ichan,ilayer) / scatt_aux % h (ichan,ilayer)
        bh_ad   (ichan,ilayer) = 0.0_JPRB
      endif
    end do 
  end do 

  !* Reset      
  dp_ad (:,:) = 0.0_JPRB
  dm_ad (:,:) = 0.0_JPRB

  deallocate (b , dx , b_ad, dx_ad)
  deallocate (ab, ab_ad, tab, ipiv)

  if (lhook) call dr_hook('RTTOV_BOUNDARYCONDITIONS_AD',1_jpim,zhook_handle)

End subroutine rttov_boundaryconditions_ad
