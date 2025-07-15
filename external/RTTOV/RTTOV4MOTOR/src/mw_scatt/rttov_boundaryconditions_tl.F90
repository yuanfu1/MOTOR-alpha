!
Subroutine rttov_boundaryconditions_tl (&         
     & ccthres,       &! in
     & nlevels,       &! in
     & nchannels,     &! in
     & lprofiles,     &! in
     & scatt_aux,     &! in
     & scatt_aux_tl,  &! in
     & ftop,          &! in
     & ftop_tl,       &! in
     & dp,            &! out
     & dp_tl,         &! out
     & dm,            &! out
     & dm_tl)          ! out 


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
  !  1.0       09/2002   Initial version      (E. Moreau)
  !  1.1       05/2003   RTTOV7.3 compatible  (F. Chevallier)
  !  1.2       03/2004   Included polarimetry (R. Saunders)
  !  1.3       11/2004   Clean-up             (P. Bauer)
  !  1.4       11/2007   RTTOV9 version       (A. Geer)
  !  1.5       07/2008   Clear sky speed-ups  (A. Geer)
  !  1.6       03/2010   Optimisation         (A. Geer)
  !  1.7       11/2017   R/T now done with radiances, not Tb (A. Geer)
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
       & rttov_profile_scatt_aux 

  Use parkind1, Only : jpim     ,jprb
!INTF_OFF
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
#ifndef _RTTOV_ECMWF
  USE rttov_lapack_mod, ONLY : dgbtrf, dgbtrs
#endif
!INTF_ON
  Implicit none

!* Subroutine arguments:
  Real    (Kind=jprb), Intent (in) :: ccthres
  Integer (Kind=jpim), Intent (in) :: nlevels    ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nchannels  ! Number of radiances
  Integer (Kind=jpim), Intent (in) :: lprofiles (nchannels) ! Profile indices

  Type (rttov_profile_scatt_aux), Intent (in)    :: scatt_aux      ! Auxiliary profile variables for RTTOV_SCATT
  Type (rttov_profile_scatt_aux), Intent (in)    :: scatt_aux_tl   ! Auxiliary profile variables for RTTOV_SCATT

  Real (Kind=jprb), Intent  (in), dimension (nchannels)            :: ftop
  Real (Kind=jprb), Intent  (in), dimension (nchannels)            :: ftop_tl
  Real (Kind=jprb), Intent (out), dimension (nchannels,nlevels) :: dp   , dm
  Real (Kind=jprb), Intent (out), dimension (nchannels,nlevels) :: dp_tl, dm_tl

!INTF_END

!* Local variables
  Real    (Kind=jprb), dimension (nchannels,nlevels) :: lh_p   , lh_m   , bh
  Real    (Kind=jprb), dimension (nchannels,nlevels) :: lh_p_tl, lh_m_tl, bh_tl
  Real    (Kind=jprb), allocatable :: b (:)
  Real    (Kind=jprb), allocatable :: b_tl (:)
  Real    (Kind=jprb), allocatable :: ab_tl (:,:)
  Real    (Kind=jprb)              :: ztmp, ztmp_tl
  Integer (Kind=jpim)              :: ilayer, jlayer, klayer, ilin, icol, iband, uband
  Integer (Kind=jpim)              :: nmaxdim, iprof, ichan, jj, ii, mcly

!* Lapack/ESSL
#ifdef _RTTOV_ECMWF 
  Real    (Kind=jprb), allocatable :: ab (:,:)
  Real    (Kind=jprb), allocatable :: dx (:)
  Real    (Kind=jprb), allocatable :: dx_tl (:)
  Integer (Kind=jpim), allocatable :: ipiv (:)
  Integer (Kind=jpim)              :: ndim, kl, ku, ldab, info, nrhs
#else
  Double Precision,    allocatable :: ab (:,:)
  Double Precision,    allocatable :: dx (:)
  Double Precision,    allocatable :: dx_tl (:)
  Integer,             allocatable :: ipiv (:)
  Integer                          :: ndim, kl, ku, ldab, info, nrhs
#endif
  Character (len=1)                :: trans

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_BOUNDARYCONDITIONS_TL',0_jpim,zhook_handle)

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
  allocate (b_tl  (nmaxdim     ))
  allocate (dx    (nmaxdim     ))
  allocate (dx_tl (nmaxdim     ))
  allocate (ab    (ldab,nmaxdim))
  allocate (ab_tl (ldab,nmaxdim))
  allocate (ipiv  (nmaxdim     ))

!* Reset      
  dp_tl (:,:) = 0.0_JPRB
  dm_tl (:,:) = 0.0_JPRB
  dp    (:,:) = 0.0_JPRB
  dm    (:,:) = 0.0_JPRB

  !* Channels * Profiles      
  do ilayer=1,nlevels
    do ichan = 1, nchannels
      iprof = lprofiles (ichan)
      if( scatt_aux % cfrac (iprof) > ccthres .and. ilayer >= scatt_aux % mclayer(ichan)) then 

        bh    (ichan,ilayer) = scatt_aux    % b1 (ichan,ilayer) / scatt_aux    % h (ichan,ilayer)

        lh_p    (ichan,ilayer) = (1.0_JPRB + scatt_aux % lambda (ichan,ilayer) / scatt_aux % h (ichan,ilayer)) 

        lh_m    (ichan,ilayer) = (1.0_JPRB - scatt_aux % lambda (ichan,ilayer) / scatt_aux % h (ichan,ilayer))

      endif
    enddo
  enddo
  do ilayer=1,nlevels
    do ichan = 1, nchannels
      iprof = lprofiles (ichan)
      if( scatt_aux % cfrac (iprof) > ccthres .and. ilayer >= scatt_aux % mclayer(ichan)) then 

        bh_tl (ichan,ilayer) = scatt_aux_tl % b1 (ichan,ilayer) / scatt_aux    % h (ichan,ilayer) &
                     & - scatt_aux    % b1 (ichan,ilayer) * scatt_aux_tl % h (ichan,ilayer) &
                     & /( scatt_aux    % h  (ichan,ilayer) * scatt_aux    % h (ichan,ilayer) )
        lh_p_tl (ichan,ilayer) = scatt_aux_tl % lambda (ichan,ilayer) / scatt_aux    % h (ichan,ilayer) &
                     & - scatt_aux    % lambda (ichan,ilayer) * scatt_aux_tl % h (ichan,ilayer) &
                     & /( scatt_aux    % h      (ichan,ilayer) * scatt_aux    % h (ichan,ilayer) )
        lh_m_tl (ichan,ilayer) = - 1.0_JPRB * lh_p_tl (ichan,ilayer)
      endif
    enddo
  enddo

  do ichan = 1, nchannels
    iprof = lprofiles (ichan)
    if( scatt_aux % cfrac (iprof) > ccthres .and. scatt_aux % mclayer (ichan) <= nlevels ) then 

      mcly = scatt_aux % mclayer (ichan)
      ndim = 2 * (nlevels - mcly + 1)

      do ilayer = 2, ndim - 2, 2
        jlayer = nlevels - ilayer / 2 + 1
        klayer = jlayer - 1

        ilin =  ilayer
        icol = (ilayer - 1)

        ztmp    = exp (scatt_aux    % lambda (ichan,jlayer) * scatt_aux    % dz (iprof,jlayer))
        ztmp_tl =     (scatt_aux_tl % lambda (ichan,jlayer) * scatt_aux    % dz (iprof,jlayer) + &
                     & scatt_aux    % lambda (ichan,jlayer) * scatt_aux_tl % dz (iprof,jlayer)) * ztmp 

        !* From downward fluxes at i-th interface (@ level=dz for jlayer == level=0 for klayer)
        if(icol > 1) then
          ab_tl (iband+3,icol-1) = 0.0_JPRB
          ab    (iband+3,icol-1) = 0.0_JPRB
        endif

        ab_tl (iband+2 ,icol  ) = lh_p_tl (ichan,jlayer) * ztmp + lh_p (ichan,jlayer) * ztmp_tl
        ab    (iband+2 ,icol  ) = lh_p    (ichan,jlayer) * ztmp

        ab_tl (iband+1 ,icol+1) = lh_m_tl (ichan,jlayer) / ztmp - lh_m (ichan,jlayer) * ztmp_tl /(ztmp * ztmp)
        ab    (iband+1 ,icol+1) = lh_m    (ichan,jlayer) / ztmp

        ab_tl (iband   ,icol+2) = -1.0_JPRB * lh_p_tl (ichan,klayer) 
        ab    (iband   ,icol+2) = -1.0_JPRB * lh_p    (ichan,klayer) 

        ab_tl (iband-1 ,icol+3) = -1.0_JPRB * lh_m_tl (ichan,klayer) 
        ab    (iband-1 ,icol+3) = -1.0_JPRB * lh_m    (ichan,klayer) 

        b_tl (ilin  ) = bh_tl (ichan,klayer) - bh_tl (ichan,jlayer)
        b    (ilin  ) = bh    (ichan,klayer) - bh    (ichan,jlayer)

        !* From upward fluxes at i-th interface (@ level=dz for jlayer == level=0 for klayer)
        ab_tl (iband+3 ,icol  ) = lh_m_tl (ichan,jlayer) * ztmp + lh_m (ichan,jlayer) * ztmp_tl
        ab    (iband+3 ,icol  ) = lh_m    (ichan,jlayer) * ztmp

        ab_tl (iband+2 ,icol+1) = lh_p_tl (ichan,jlayer) / ztmp - lh_p (ichan,jlayer) * ztmp_tl /(ztmp * ztmp)
        ab    (iband+2 ,icol+1) = lh_p    (ichan,jlayer) / ztmp

        ab_tl (iband+1 ,icol+2) = -1.0_JPRB * lh_m_tl (ichan,klayer) 
        ab    (iband+1 ,icol+2) = -1.0_JPRB * lh_m    (ichan,klayer) 

        ab_tl (iband   ,icol+3) = -1.0_JPRB * lh_p_tl (ichan,klayer) 
        ab    (iband   ,icol+3) = -1.0_JPRB * lh_p    (ichan,klayer) 

        if(icol < ndim-3) then
          ab_tl (iband-1,icol+4) = 0.0_JPRB 
          ab    (iband-1,icol+4) = 0.0_JPRB 
        endif

        b_tl (ilin+1) = bh_tl (ichan,jlayer) - bh_tl (ichan,klayer)
        b    (ilin+1) = bh    (ichan,jlayer) - bh    (ichan,klayer)
      end do

      !* From boundary conditions at bottom of the atmosphere with r_sfc=1-e_sfc
      ztmp    = (2.0_JPRB - scatt_aux % ems_bnd (ichan)) * scatt_aux % lambda (ichan,nlevels) / &
        & scatt_aux % h (ichan,nlevels)
      ztmp_tl =   -1.0_JPRB * scatt_aux_tl % ems_bnd (ichan)  * scatt_aux    % lambda (ichan,nlevels) / &
        & scatt_aux    % h (ichan,nlevels) &
             & + (2.0_JPRB - scatt_aux    % ems_bnd (ichan)) * scatt_aux_tl % lambda (ichan,nlevels) / &
             &scatt_aux    % h (ichan,nlevels) &
             & - (2.0_JPRB - scatt_aux    % ems_bnd (ichan)) * scatt_aux    % lambda (ichan,nlevels) * &
             & scatt_aux_tl % h (ichan,nlevels) &
             & /( scatt_aux % h (ichan,nlevels) * scatt_aux % h (ichan,nlevels) )

      ab_tl (iband+1,1) = scatt_aux_tl % ems_bnd (ichan) - ztmp_tl
      ab    (iband+1,1) = scatt_aux    % ems_bnd (ichan) - ztmp

      ab_tl (iband,2) = scatt_aux_tl % ems_bnd (ichan) + ztmp_tl
      ab    (iband,2) = scatt_aux    % ems_bnd (ichan) + ztmp

      ab_tl (iband-1,3) = 0.0_JPRB
      ab    (iband-1,3) = 0.0_JPRB

      b_tl (1) = scatt_aux_tl % ems_bnd (ichan) * (scatt_aux % bsfc (ichan) - scatt_aux    % b0 (ichan,nlevels)) &
            & + scatt_aux    % ems_bnd (ichan) * (scatt_aux_tl % bsfc (ichan) - scatt_aux_tl % b0 (ichan,nlevels)) &
            & - scatt_aux_tl % ems_bnd (ichan) * bh (ichan,nlevels) &
            & + (2.0_JPRB - scatt_aux % ems_bnd (ichan)) * bh_tl (ichan,nlevels) 
      b (1)    = scatt_aux % ems_bnd (ichan) * (scatt_aux % bsfc (ichan) - scatt_aux % b0 (ichan,nlevels)) &
            & + (2.0_JPRB - scatt_aux % ems_bnd (ichan)) * bh    (ichan,nlevels)  

      !* From boundary conditions at top of the atmosphere 
      ztmp    = exp (scatt_aux % lambda (ichan,mcly) * scatt_aux % dz (iprof,mcly))
      ztmp_tl = (scatt_aux_tl % lambda (ichan,mcly) * scatt_aux    % dz (iprof,mcly) &
            & + scatt_aux    % lambda (ichan,mcly) * scatt_aux_tl % dz (iprof,mcly)) * ztmp 

      ab_tl (iband+3,ndim-2) = 0.0_JPRB  
      ab    (iband+3,ndim-2) = 0.0_JPRB  

      ab_tl (iband+2,ndim-1) = lh_p_tl (ichan,mcly) * ztmp + lh_p (ichan,mcly) * ztmp_tl
      ab    (iband+2,ndim-1) = lh_p    (ichan,mcly) * ztmp

      ab_tl (iband+1,ndim  ) = lh_m_tl (ichan,mcly) / ztmp - lh_m (ichan,mcly) * ztmp_tl /( ztmp * ztmp)
      ab    (iband+1,ndim  ) = lh_m    (ichan,mcly) / ztmp

      b_tl (ndim) = ftop_tl (ichan) - scatt_aux_tl % bn (ichan,mcly) - bh_tl (ichan,mcly)
      b    (ndim) = ftop    (ichan) - scatt_aux    % bn (ichan,mcly) - bh    (ichan,mcly)

      !* Solve equations A * DX = B, forward          
      call dgbtrf (ndim, ndim, kl, ku, ab, ldab, ipiv, info)                  

      dx (1:ndim) = b (1:ndim)
     
      call dgbtrs (trans, ndim, kl, ku, nrhs, ab, ldab, ipiv, dx, ndim, info)

      !* Solve equations A * DX = B, tangent-linear     

      ! Following code replaces this:     
      !  dx_tl (:) = b_tl (:) - matmul (a_tl, dx)
      ! since we are now working with "general band matrix representation" (see ESSL manuals) 
      ! and sparse a_tl has been replaced with packed ab_tl
      dx_tl(1:ndim) = b_tl(1:ndim)
      do jj = 1, ndim
        do ii = max(1_jpim,jj-ku), min(ndim,jj+kl)
          dx_tl(jj) = dx_tl(jj) - ab_tl(kl+ku+jj-ii+1,ii) * dx(ii)
        end do
      end do

      call dgbtrs (trans, ndim, kl, ku, nrhs, ab, ldab, ipiv, dx_tl, ndim, info)
     
      !* Decompose D+ and D-
      do ilayer = 2, ndim, 2
        jlayer = nlevels - ilayer / 2 + 1
        
        dp_tl (ichan,jlayer) = dx_tl (ilayer-1)
        dp    (ichan,jlayer) = dx    (ilayer-1)
        dm_tl (ichan,jlayer) = dx_tl (ilayer  )
        dm    (ichan,jlayer) = dx    (ilayer  )
      end do

    endif
  end do

  deallocate (b, dx, b_tl, dx_tl, ab, ab_tl, ipiv)
 
  if (lhook) call dr_hook('RTTOV_BOUNDARYCONDITIONS_TL',1_jpim,zhook_handle)

End subroutine rttov_boundaryconditions_tl
