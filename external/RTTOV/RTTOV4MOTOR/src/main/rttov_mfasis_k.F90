! Description:
!> @file
!!    K of MFASIS fast visible/near-IR scattering model.
!
!> @brief
!!    K of MFASIS fast visible/near-IR scattering model.
!!
!!
!! @param[out]    err               status on exit
!! @param[in]     chanprof          specifies channels and profiles to simulate
!! @param[in]     chanflag          flags to indicate which channels with LUT available
!! @param[in]     opts              options to configure the simulations
!! @param[in]     profiles          input atmospheric profiles and surface variables
!! @param[in,out] profiles_k        input profile increments
!! @param[in]     profiles_int      profiles in internal units
!! @param[in,out] profiles_int_k    profile increments in internal units
!! @param[in]     coefs             coefficients structure for instrument to simulate
!! @param[in]     ircld             information on cloud columns
!! @param[in,out] ircld_k           cloud column increments
!! @param[in]     aux               additional internal profile variables
!! @param[in,out] aux_k             additional internal profile variable increments
!! @param[in]     reflectance       surface BRDFs
!! @param[in,out] reflectance_k     surface BRDF increments
!! @param[in]     solar_spectrum    TOA solar irradiance for each channel
!! @param[in]     trans_scatt_ir    cloud/aerosol optical depths
!! @param[in,out] trans_scatt_ir_k  cloud/aerosol optical depth increments
!! @param[in]     mfasis_refl       quantities computed by rttov_mfasis used by TL/AD/K
!! @param[in,out] radiance_k        input gradient wrt radiances
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
!    Copyright 2017, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_mfasis_k( &
              err,              &
              chanprof,         &
              chanflag,         &
              opts,             &
              profiles,         &
              profiles_k,       &
              profiles_int,     &
              profiles_int_k,   &
              coefs,            &
              ircld,            &
              ircld_k,          &
              aux,              &
              aux_k,            &
              reflectance,      &
              reflectance_k,    &
              solar_spectrum,   &
              trans_scatt_ir,   &
              trans_scatt_ir_k, &
              mfasis_refl,      &
              radiance_k)

#include "throw.h"

  USE parkind1, ONLY : jpim, jprb, jplm

  USE rttov_types, ONLY :        &
    rttov_chanprof,              &
    rttov_options,               &
    rttov_profile,               &
    rttov_coefs,                 &
    rttov_ircld,                 &
    rttov_profile_aux,           &
    rttov_reflectance,           &
    rttov_transmission_scatt_ir, &
    rttov_radiance,              &
    rttov_mfasis_refl
!INTF_OFF
  USE rttov_types, ONLY : &
    rttov_coef_mfasis,           &
    rttov_mfasis_axis

  USE rttov_const, ONLY : &
    wcl_opac_deff,               &
    pi,                          &
    deg2rad,                     &
    pi_r,                        &
    clw_scheme_deff,             &
    mfasis_cld,                  &
!     mfasis_aer,                  &
    mfasis_dim_albedo,           &
    mfasis_dim_kfourier,         &
    mfasis_dim_lfourier,         &
    mfasis_dim_opdp,             &
    mfasis_dim_effdia,           &
    mfasis_dim_scaangle,         &
    gas_id_watervapour,          &
    gas_mass,                    &
    mair,                        &
    gravity,                     &
    nwcl_max
!INTF_ON
  IMPLICIT NONE
  INTEGER(jpim),                      INTENT(OUT)   :: err
  TYPE(rttov_chanprof),               INTENT(IN)    :: chanprof(:)
  LOGICAL(jplm),                      INTENT(IN)    :: chanflag(SIZE(chanprof))
  TYPE(rttov_options),                INTENT(IN)    :: opts
  TYPE(rttov_profile),                INTENT(IN)    :: profiles(:)
  TYPE(rttov_profile),                INTENT(INOUT) :: profiles_k(SIZE(chanprof))
  TYPE(rttov_profile),                INTENT(IN)    :: profiles_int(:)
  TYPE(rttov_profile),                INTENT(INOUT) :: profiles_int_k(SIZE(chanprof))
  TYPE(rttov_coefs),                  INTENT(IN)    :: coefs
  TYPE(rttov_ircld),                  INTENT(IN)    :: ircld
  TYPE(rttov_ircld),                  INTENT(INOUT) :: ircld_k
  TYPE(rttov_profile_aux),            INTENT(IN)    :: aux
  TYPE(rttov_profile_aux),            INTENT(INOUT) :: aux_k
  TYPE(rttov_reflectance),            INTENT(IN)    :: reflectance(SIZE(chanprof))
  TYPE(rttov_reflectance),            INTENT(INOUT) :: reflectance_k(SIZE(chanprof))
  REAL(jprb),                         INTENT(IN)    :: solar_spectrum(SIZE(chanprof))
  TYPE(rttov_transmission_scatt_ir),  INTENT(IN)    :: trans_scatt_ir
  TYPE(rttov_transmission_scatt_ir),  INTENT(INOUT) :: trans_scatt_ir_k
  TYPE(rttov_mfasis_refl),            INTENT(IN)    :: mfasis_refl(0:,:)
  TYPE(rttov_radiance),               INTENT(INOUT) :: radiance_k
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(jpim)              :: nchanprof, prof, ncolms, npar, chan
  INTEGER(jpim)              :: nlay, nlev, ndim, nid
  INTEGER(jpim)              :: d, n, i, j, k, cc, par, jj
  TYPE(rttov_coef_mfasis)    :: mfasis_coefs
  REAL(jprb),    ALLOCATABLE :: ip   (:)
  REAL(jprb),    ALLOCATABLE :: ip_k(:), iw_k(:)

  REAL(jprb),    ALLOCATABLE :: q_mxr(:),dvap(:), dvap_1(:)
  REAL(jprb),    ALLOCATABLE :: q_mxr_k(:),dvap_k(:)
  REAL(jprb),    ALLOCATABLE :: array1(:)                  
  REAL(jprb)                 :: wvint_bot, wvint_top 
  REAL(jprb)                 :: wvint_bot_k, wvint_top_k 
  REAL(jprb)                 :: iw_wv_k(3)
  INTEGER(jpim)              :: n_wvdim


  INTEGER(jpim), ALLOCATABLE :: di(:)
  REAL(jprb),    ALLOCATABLE :: od   (:,:), od_eff   (:), ed   (:,:), ed_eff(:), ed_aux(:)
  REAL(jprb),    ALLOCATABLE :: od_v1(:,:), od_eff_v1(:), ed_eff_v1(:) !versions of nonlinear quantities needed for TL/AD/K
  REAL(jprb),    ALLOCATABLE :: od_k(:,:), od_eff_k(:), ed_k(:,:), ed_eff_k(:)

  REAL(jprb)                 :: albedo
  REAL(jprb)                 :: albedo_k
  TYPE(rttov_mfasis_axis)    :: axis
  REAL(jprb)                 :: colwei_clr, colwei_cc
  REAL(jprb)                 :: refl_k, rflcol_k, refl_clr_k, colwei_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(jprb)                 :: od1_int,od1_mixed,od2_mixed
  REAL(jprb)                 :: od1_int_k,od1_mixed_k,od2_mixed_k, ooo1_k
  REAL(jprb),    ALLOCATABLE :: ooo1(:)
  REAL(jprb),    ALLOCATABLE :: fac_mx(:)
  REAL(jprb),    ALLOCATABLE :: fac_mx_k(:)
  REAL(jprb)                 :: fac_sw,xxx1,xxx2            ! variables    used for switching on/off mixed cloud correction
  REAL(jprb)                 :: fac_sw_k,xxx1_k,xxx2_k   ! variables    used for switching on/off mixed cloud correction
  ! --------------------------------------------------------------------------------------

  TRY

  ! TODO: Adapt to water vapor correction using three LUTs
  ! TODO: Test aerosol functionality

  ! --------------------------------------------------------------------------
  ! Initialisation
  ! --------------------------------------------------------------------------
  nchanprof = SIZE(chanprof)
  nlay = profiles(1)%nlayers
  nlev = nlay+1

  ! Check whether cloud or aerosol simulation (currently both simultaneously not supported)
  IF (opts%rt_ir%addclouds ) THEN
    mfasis_coefs = coefs%coef_mfasis_cld
  ELSE
    mfasis_coefs = coefs%coef_mfasis_aer
  ENDIF

  npar = mfasis_coefs%nparticles        ! # of particles: 2 for clouds (water and ice)
                                        !                 or # of aerosols considered
  ndim = mfasis_coefs%ndims             ! # of LUT dimensions

  ! Find number of dimensions that are not interpolated (fourier coefs and albedo)
  n = 0_jpim
  DO d = 1, ndim
    IF (mfasis_coefs%lut_axes(d)%dim_type /= mfasis_dim_kfourier .AND. &
        mfasis_coefs%lut_axes(d)%dim_type /= mfasis_dim_lfourier .AND. &
        mfasis_coefs%lut_axes(d)%dim_type /= mfasis_dim_albedo ) CYCLE
    n = n + 1
  END DO

  nid = ndim - n ! Total number of dimensions to be interpolated (total - albedo - fourier indices)

  ALLOCATE( ip     (nid) , &      ! interpolation point
            ip_k  (nid) , &      ! interpolation point
            iw_k  (nid) , &      ! interpolation weight
            di     (nid), STAT=err )  ! original dimension index in nid arrays
  THROWM(err.NE.0,"Allocation of memory for rttov_mfasis_k failed")


  ! Populate di and set some relevant dimensions indices
  di(:) = -1_jpim ! Indices of dimensions to interpolate

  n = 1_jpim
  DO d = 1, ndim
    SELECT CASE (mfasis_coefs%lut_axes(d)%dim_type)
    CASE (mfasis_dim_kfourier, mfasis_dim_lfourier, mfasis_dim_albedo)

    CASE DEFAULT
      di(n) = d
      n = n + 1
    END SELECT
  ENDDO


  ALLOCATE( od_eff    (npar)    , & ! Effective total optical depth per type of particle
            od_eff_v1 (npar)    , & ! version of od_eff for adjoint computations
            od_eff_k (npar)    , & ! Effective total optical depth per type of particle
            od     (npar, nlay) , & ! Effective optical depth per type of particle and layer
            od_v1  (npar, nlay), &     ! od saved for adjoint computations
            od_k  (npar, nlay), STAT=err)  ! Effective optical depth per type of particle and layer
  THROWM(err.NE.0,"Allocation of memory for rttov_mfasis_k failed")

  ! For the time being we assume effective diameters will be out of LUT for aerosols
  IF ( mfasis_coefs%file_type .EQ. 1 ) THEN
    ALLOCATE( ed_eff (npar)      , & ! Total effective diameter per type of particle
              ed_eff_v1(npar) ,  &      ! version of ed_eff for adjoint computations
              ed_aux (npar)      , & ! Effective diameters to be used with 0 opdp to avoid flagging
              ed_eff_k(npar)    , & ! Total effective diameter per type of particle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              fac_mx (nlay),      &         ! fraction of cloud considered as (potentially) mixed cloud
              fac_mx_k (nlay),   &         ! fraction of cloud considered as (potentially) mixed cloud
              ooo1 (nlay),        &         ! discriminator for switching on mixed layer summation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ed   (npar, nlay)  , & ! Effective diameter per type of particle and layer
              ed_k(npar, nlay),   &
              q_mxr(nlev),         &
              q_mxr_k(nlev),         &
              dvap (nlay),         &
              dvap_1(nlay),         &
              dvap_k (nlay),         &
              array1(0:nlay+1),       &
                              STAT=err) ! 
                           
    THROWM(err.NE.0,"Allocation of memory for rttov_mfasis_k failed")

    ! Effective diameters to be used with 0 opdp to avoid incorrect flagging for clear column
    n = 1_jpim
    DO d = 1, ndim
      IF (mfasis_coefs%lut_axes(d)%dim_type .NE. mfasis_dim_effdia) CYCLE
      ed_aux(n) = mfasis_coefs%lut_axes(d)%values(1)
      n = n + 1
    ENDDO
  ENDIF
  fac_mx (:) = 0
  ooo1 (:)   = 0

!=====================================================================================
! Make sure all lokal linear variables are initialized as zeros
!=====================================================================================
  ip_k(:)   = 0
  iw_k(:)   = 0
  iw_wv_k(:)= 0
  wvint_bot_k=0._jprb
  wvint_top_k=0._jprb


  od_k     (:,:)  = 0
  od_eff_k (:)    = 0
  ed_k     (:,:)  = 0
  ed_eff_k (:)    = 0
  fac_mx_k (:)    = 0

  albedo_k  = 0
  refl_k    = 0
  rflcol_k    = 0
  refl_clr_k  = 0
  colwei_k    = 0

  od1_int_k    = 0
  od1_mixed_k  = 0
  od2_mixed_k  = 0
  ooo1_k       = 0

  fac_sw_k   = 0
  xxx1_k   = 0
  xxx2_k   = 0
!=====================================================================================

  ! --------------------------------------------------------------------------
  ! Channel loop
  ! --------------------------------------------------------------------------
  DO i = 1, nchanprof ! channel loop
    IF (.NOT. chanflag(i)) CYCLE
    prof = chanprof(i)%prof
    chan = chanprof(i)%chan

    n_wvdim=size(mfasis_coefs%lut(chan)%qint,2)
    dvap(:)=0._jprb
    dvap_k(:)=0._jprb
    IF(n_wvdim==3) THEN
      q_mxr(:) =  profiles_int(prof)%q(:)* gas_mass(gas_id_watervapour)/  &
                  (mair * 1.E06_jprb +profiles_int(prof)%q(:)*gas_mass(gas_id_watervapour))
!     q_mxr_k(:)= profiles_int_k(prof)%q(:)*gas_mass(gas_id_watervapour) * (1._jprb - q_mxr(:) ) &
!                 /(mair * 1.E06_jprb +profiles_int(prof)%q(:)*gas_mass(gas_id_watervapour))

      dvap  (1:nlay)  = 100._jprb * (profiles(prof)%p(2:nlev)-profiles(prof)%p(1:nlev-1))* &
                                           0.5_jprb * (q_mxr(2:nlev)+q_mxr(1:nlev-1))/gravity
      dvap_1(1:nlay)  = dvap(1:nlay) 
!     dvap_k(1:nlay)=100._jprb *((profiles_k(prof)%p(2:nlev)-profiles_k(prof)%p(1:nlev-1))* &
!                                          0.5_jprb * (q_mxr(2:nlev)+q_mxr(1:nlev-1))/gravity  &
!                                 + (profiles(prof)%p(2:nlev)-profiles(prof)%p(1:nlev-1))* &
!                                    0.5_jprb * (q_mxr_k(2:nlev)+q_mxr_k(1:nlev-1))/gravity  &
!                                )  

      jj = aux%s(prof)%nearestlev_surf - 1
      IF(jj>=1 .AND. jj<=nlay) THEN
!       dvap_k(jj)=dvap_k(jj)*(1._jprb - aux   %s(prof)%pfraction_surf)           &
!                       - dvap   (jj) * aux_k%s(prof)%pfraction_surf
        dvap   (jj)=dvap   (jj)*(1._jprb - aux   %s(prof)%pfraction_surf)
      ENDIF
    ENDIF


    colwei_clr = ircld%xcolclr(prof) ! clear column

    ip(:) = -1

    IF ( mfasis_coefs%file_type == mfasis_cld ) THEN
      ncolms = ircld%ncolumn(prof)
    ELSE
      ncolms = 1_jpim
    ENDIF


!=====================================================================================
!  linear computations independent of column cc
!=====================================================================================
!================================
!    radiance_k%total(i) = refl_k * solar_spectrum(i) * pi_r  &
!      * COS(profiles(prof)%sunzenangle * deg2rad)

! To be consistent with other RTTOV calls we only use input AD/K increments in
! radiance%total and we assume all other radiance arrays are zero
    refl_k    = refl_k     + solar_spectrum(i) * pi_r   &
      * COS(profiles(prof)%sunzenangle * deg2rad) * radiance_k%total(i)



    colwei_k  =colwei_k   + mfasis_refl(0,i)%refl*refl_k
    refl_clr_k=refl_clr_k + colwei_clr     *refl_k

!    colwei_k = ircld_k%xcolclr(i) ! clear column
    ircld_k%xcolclr(i)=ircld_k%xcolclr(i) + colwei_k
    colwei_k=0._jprb

!   refl_clr_tl = refl_clr_tl+sum(iw_wv_tl(1:n_wvdim)*mfasis_refl(0,i)%refl_wv(1:n_wvdim))
    iw_wv_k(1:n_wvdim) = iw_wv_k(1:n_wvdim) + mfasis_refl(0,i)%refl_wv(1:n_wvdim)*refl_clr_k

    DO d = 1, nid
      IF ( mfasis_coefs%lut_axes(di(d))%dim_type .EQ. mfasis_dim_scaangle ) THEN
!       refl_clr_k = refl_clr_k+mfasis_refl(0,i)%refl_lin_coef(d)*albedo_k
        albedo_k = albedo_k + mfasis_refl(0,i)%refl_lin_coef(d)* refl_clr_k
      ELSE
        iw_k(d)    = iw_k(d)   +mfasis_refl(0,i)%refl_lin_coef(d)* refl_clr_k
      ENDIF
    ENDDO
    refl_clr_k = 0


    DO d = 1, nid
      IF ( mfasis_coefs%lut_axes(di(d))%dim_type .EQ. mfasis_dim_opdp ) THEN
        iw_k(d) = 0._jprb
      ENDIF
    ENDDO

!=============================================================================================================
!   IF(n_wvdim==3) THEN
!     wvint_bot_tl=0._jprb
!     wvint_top_tl=0._jprb
!     DO j = 1, nlay
!       wvint_top_tl = wvint_top_tl +  dvap_tl(j)
!     ENDDO
!   ENDIF
!   iw_wv_tl(:) = 0._jprb
!   if( n_wvdim==3) THEN
!       iw_wv_tl(2)=                  wvint_bot_tl  / &
!            (mfasis_coefs%lut(chan)%qint(1,2) - mfasis_coefs%lut(chan)%qint(1,1))
!       iw_wv_tl(3)=                  wvint_top_tl  / &
!            (mfasis_coefs%lut(chan)%qint(2,3) - mfasis_coefs%lut(chan)%qint(2,1))
!       iw_wv_tl(1) =       - iw_wv_tl(2) - iw_wv_tl(3)
!    ENDIF
!=============================================================================================================
     IF( n_wvdim==3) THEN
       iw_wv_k(2) = iw_wv_k(2) - iw_wv_k(1)
       iw_wv_k(3) = iw_wv_k(3) - iw_wv_k(1)
       iw_wv_k(1) = 0
       wvint_top_k = wvint_top_k + iw_wv_k(3) / &
              (mfasis_coefs%lut(chan)%qint(2,3) - mfasis_coefs%lut(chan)%qint(2,1))
       wvint_bot_k = wvint_bot_k + iw_wv_k(2)  / &
               (mfasis_coefs%lut(chan)%qint(1,2) - mfasis_coefs%lut(chan)%qint(1,1))
     ENDIF
     iw_wv_k(:) = 0._jprb
     IF(n_wvdim==3) THEN
       DO j = nlay, 1, -1
!        wvint_top_tl = wvint_top_tl +  dvap_tl(j)
         dvap_k(j) = dvap_k(j) + wvint_top_k
       ENDDO
       wvint_bot_k=0._jprb
       wvint_top_k=0._jprb
     ENDIF




!=====================================================================================

    DO cc = 1, ncolms
!==================================================================================
! First compute nonlinear quatities needed in linear computations (for given channel "i" and column "cc")
!==================================================================================
      ! --------------------------------------------------------------------------
      ! Set column weight, optical depth and effective diameters
      ! --------------------------------------------------------------------------
      IF ( mfasis_coefs%file_type .EQ. mfasis_cld ) THEN ! clouds
        colwei_cc = ircld%xcol(cc+1,prof) - ircld%xcol(cc,prof)
!       colwei_k = ircld_k%xcol(cc+1,i) - ircld_k%xcol(cc,i)
        od   (:,:) = 0._jprb
!       od_k(:,:) = 0._jprb
        ed   (:,:) = 0._jprb
!       ed_k(:,:) = 0._jprb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        od1_int = 0._jprb
        od1_mixed= 0._jprb
        od2_mixed= 0._jprb
        fac_mx(:)= 0._jprb
        wvint_bot=0._jprb
        wvint_top=0._jprb

!       od1_int_k = 0._jprb
!       od1_mixed_k= 0._jprb
!       od2_mixed_k= 0._jprb
!       fac_mx_k(:)= 0._jprb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO j = 1, nlay
          IF (j > aux%s(prof)%nearestlev_surf - 1) EXIT    ! Layer is entirely below surface pressure so nothing more to do
          IF ( ircld%icldarr(cc,j,prof) .EQ. 1 ) THEN
            od   (1,j) = SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i)) ! water clouds
!           od_k(1,j) = SUM(trans_scatt_ir_k%opdpext(1:nwcl_max,j,i)) ! water clouds
            od   (2,j) = trans_scatt_ir%opdpext(nwcl_max+1,j,i)        ! ice clouds
!           od_k(2,j) = trans_scatt_ir_k%opdpext(nwcl_max+1,j,i)        ! ice clouds
            IF (j == aux%s(prof)%nearestlev_surf - 1) THEN
              ! Modify optical depth in partial layer above surface
              od_v1(:,j) = od   (:,j)
!             od_k(:,j) = od_k(:,j) * (1 - aux%s(prof)%pfraction_surf)  &
!                              - od_v1(:,j) * aux_k%s(i)%pfraction_surf
              od   (:,j) = od   (:,j) * (1 - aux%s(prof)%pfraction_surf)
            ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            od1_int = od1_int + od(1,j)
!           od1_int_k = od1_int_k + od_k(1,j)
            ooo1(j)=od1_int-opts%dev%od1_thresh
!           ooo1_k=od1_int_k
            fac_mx(j)=(1._jprb + TANH(ooo1(j)/opts%dev%o_del1))/2._jprb

!           fac_mx_k(j)= 0._jprb
!           IF(abs(ooo1(j)/opts%dev%o_del1) < 5._jprb) &
!           fac_mx_k(j)=(ooo1_k/opts%dev%o_del1) / (COSH(ooo1(j)/opts%dev%o_del1)**2) /2._jprb

            od1_mixed=od1_mixed+fac_mx(j)*od(1,j)
            od2_mixed=od2_mixed+fac_mx(j)*od(2,j)
!           od1_mixed_k=od1_mixed_k+fac_mx_k(j)*od(1,j)+fac_mx(j)*od_k(1,j)
!           od2_mixed_k=od2_mixed_k+fac_mx_k(j)*od(2,j)+fac_mx(j)*od_k(2,j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF ( coefs%coef_mfasis_cld%clw_scheme == clw_scheme_deff ) THEN ! Deff water clouds scheme
              ed   (1,j) = aux%clw_dg(j,prof)
!           ed_k(1,j) = aux_k%clw_dg(j,i)
            ELSE ! OPAC water cloud parametrization
              ed   (1,j) = SUM(wcl_opac_deff(1:nwcl_max) * trans_scatt_ir%opdpext(1:nwcl_max,j,i))
!             ed_k(1,j) = SUM(wcl_opac_deff(1:nwcl_max) * trans_scatt_ir_k%opdpext(1:nwcl_max,j,i))
              IF ( od(1,j) .GT. 0. ) THEN
                ! TODO: take into account optical depth surface correction in eff diameter calculation
                ed   (1,j) =  ed   (1,j)/SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i))
!               ed_k(1,j) =  ed_k(1,j)/SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i)) - &
!                             ed(1,j)*SUM(trans_scatt_ir_k%opdpext(1:nwcl_max,j,i))/&
!                             SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i))
              ELSE
                ed   (1,j) = 0._jprb
!               ed_k(1,j) = 0._jprb
              ENDIF
            ENDIF
            ed   (2,j) = aux%ice_dg(j,prof)
!           ed_k(2,j) = aux_k%ice_dg(j,prof)
          ELSEIF(j>1) THEN
            fac_mx(j) = fac_mx(j-1)
          ENDIF
          IF(n_wvdim==3) THEN
            wvint_bot = wvint_bot + fac_mx(j)          * dvap(j)
            wvint_top = wvint_top + (1._jprb-fac_mx(j))* dvap(j)
          ENDIF
        ENDDO
        DO k = 1, npar
          od_eff   (k) = SUM( od   (k,:))
          od_eff_v1(k) = od_eff   (k)
!         od_eff_k(k) = SUM( od_k(k,:))
          IF ( od_eff_v1(k) .GT. 0. ) THEN
            ed_eff   (k) = SUM( od(k,:)*ed(k,:))/od_eff_v1(k)
            ed_eff_v1(k) = ed_eff(k)
!           ed_eff_k(k) = SUM( od_k(k,:)*ed(k,:) + od(k,:)*ed_k(k,:))/od_eff_v1(k)   &
!                          - ed_eff_v1(k) * od_eff_k(k)/od_eff_v1(k)

          ELSE
            ed_eff   (k) = ed_aux(k)
!           ed_eff_k(k) = 0._jprb
          ENDIF
        ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        xxx1= opts%dev%qq_1*od1_mixed-od2_mixed
        xxx2= opts%dev%qq_2*od_eff(1)-od_eff(2)
!       xxx1_k= opts%dev%qq_1*od1_mixed_k-od2_mixed_k
!       xxx2_k= opts%dev%qq_2*od_eff_k(1)-od_eff_k(2)
        fac_sw= (1._jprb + TANH(xxx1/opts%dev%x_del1))*(1._jprb + TANH(xxx2/opts%dev%x_del2))/4._jprb

!       fac_sw_k=0._jprb
!       IF(abs(xxx1/opts%dev%x_del1) < 5._jprb) &
!         fac_sw_k= (xxx1_k/opts%dev%x_del1) / &
!                    (COSH(xxx1/opts%dev%x_del1)**2) *(1._jprb + TANH(xxx2/opts%dev%x_del2)) /4._jprb
!       IF(abs(xxx2/opts%dev%x_del2) < 5._jprb) &
!         fac_sw_k= fac_sw_k + (xxx2_k/opts%dev%x_del2) / &
!                    (COSH(xxx2/opts%dev%x_del2)**2) *(1._jprb + TANH(xxx1/opts%dev%x_del1)) /4._jprb


        od_eff(1) = od_eff(1) + fac_sw*od2_mixed
        od_eff(2) = od_eff(2) - fac_sw*od2_mixed
        od_eff_v1(2) = od_eff(2)
        IF(od_eff_v1(2) < 0) od_eff(2) = 0._jprb ! Avoid negative optical depths (from mixed-phase correction)
!       od_eff_k(1) = od_eff_k(1) + fac_sw_k*od2_mixed + fac_sw*od2_mixed_k
!       od_eff_k(2) = od_eff_k(2) -(fac_sw_k*od2_mixed + fac_sw*od2_mixed_k)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ELSE  ! aerosols
        colwei_cc = 1._jprb     ! single column for aerosols
!       colwei_k = 0._jprb     ! single column for aerosols
        od   (:,:)= 0._jprb
!       od_k(:,:)= 0._jprb
        DO k = 1, npar
          par = mfasis_coefs%aer_types(k)
          DO j = 1, nlay
            IF (j > aux%s(prof)%nearestlev_surf - 1) EXIT    ! Layer is entirely below surface pressure so nothing more to do
            od(k,j)    = trans_scatt_ir%opdpext(par,j,i)
!           od_k(k,j) = trans_scatt_ir_k%opdpext(par,j,i)
            IF (j == aux%s(prof)%nearestlev_surf - 1) THEN
              ! Modify optical depth in partial layer above surface
              od_v1(:,j) = od   (:,j)
!             od_k(k,j) = od_k(k,j) * (1 - aux%s(prof)%pfraction_surf)  &
!                         - od_v1(k,j) * aux_k%s(i)%pfraction_surf
              od   (k,j) = od   (k,j) * (1 - aux%s(prof)%pfraction_surf)
            ENDIF
          ENDDO
          od_eff   (k) = SUM(od   (k,:))
!         od_eff_k(k) = SUM(od_k(k,:))
        ENDDO
      ENDIF
      ! --------------------------------------------------------------------------
      ! Fill into the interpolation point array the rest of dimensions
      ! --------------------------------------------------------------------------
      j = 1
      k = 1
      DO d = 1, nid
        SELECT CASE (mfasis_coefs%lut_axes(di(d))%dim_type)
        CASE (mfasis_dim_opdp)
          ip   (d) = od_eff   (j)
!         ip_k(d) = od_eff_k(j)
          j = j + 1
        CASE (mfasis_dim_effdia)
          ! Check IF stored as radius or diamater in LUT
          IF (mfasis_coefs%lut_axes(di(d))%name(1:1) .EQ. "R") THEN
            ip(d) = ed_eff(k)/2
!           ip_k(d) = ed_eff_k(k)/2
          ELSE
            ip(d) = ed_eff(k)
!           ip_k(d) = ed_eff_k(k)
          ENDIF
          k = k + 1
        END SELECT
      ENDDO

!==================================================================================
! Start linear computations for channel "i" and column "cc"
!==================================================================================
!     refl_k = refl_k + rflcol_k*colwei_cc + mfasis_refl(cc,i)%refl*colwei_k
      rflcol_k = rflcol_k+  colwei_cc        * refl_k
      colwei_k = colwei_k+   mfasis_refl(cc,i)%refl * refl_k                  ! rflcol !NNNNNNNNN

!     rflcol_tl = rflcol_tl+sum(iw_wv_tl(1:n_wvdim)*mfasis_refl(cc,i)%refl_wv(1:n_wvdim))
      iw_wv_k(1:n_wvdim) = iw_wv_k(1:n_wvdim) + mfasis_refl(cc,i)%refl_wv(1:n_wvdim) * rflcol_k

      DO d = 1, nid
        IF ( mfasis_coefs%lut_axes(di(d))%dim_type .EQ. mfasis_dim_scaangle ) THEN
!         rflcol_k = rflcol_k+mfasis_refl(cc,i)%refl_lin_coef(d)*albedo_k
          albedo_k = albedo_k +mfasis_refl(cc,i)%refl_lin_coef(d)* rflcol_k
        ELSE
!         rflcol_k = rflcol_k+mfasis_refl(cc,i)%refl_lin_coef(d)*iw_k(d)
          iw_k(d)  = iw_k(d) +mfasis_refl(cc,i)%refl_lin_coef(d)* rflcol_k
        ENDIF
      ENDDO
      rflcol_k = 0

!     iw_wv_tl(:) = 0._jprb
!     IF( n_wvdim==3) THEN
!         iw_wv_tl(2)=                  wvint_bot_tl  / &
!              (mfasis_coefs%lut(chan)%qint(1,2) - mfasis_coefs%lut(chan)%qint(1,1))
!         iw_wv_tl(3)=                  wvint_top_tl  / &
!              (mfasis_coefs%lut(chan)%qint(2,3) - mfasis_coefs%lut(chan)%qint(2,1))
!         iw_wv_tl(1) =       - iw_wv_tl(2) - iw_wv_tl(3)
!      ENDIF

      IF( n_wvdim==3) THEN
        iw_wv_k(2) = iw_wv_k(2) - iw_wv_k(1) 
        iw_wv_k(3) = iw_wv_k(3) - iw_wv_k(1) 
        iw_wv_k(1) = 0
        wvint_top_k = wvint_top_k + iw_wv_k(3) / &
               (mfasis_coefs%lut(chan)%qint(2,3) - mfasis_coefs%lut(chan)%qint(2,1))
        wvint_bot_k = wvint_bot_k + iw_wv_k(2)  / &
               (mfasis_coefs%lut(chan)%qint(1,2) - mfasis_coefs%lut(chan)%qint(1,1))
      ENDIF
      iw_wv_k(:) = 0._jprb




      DO d = 1, nid
        IF (mfasis_coefs%lut_axes(di(d))%dim_type .EQ. mfasis_dim_scaangle) CYCLE ! Skip scattering angle (already done)
        axis = mfasis_coefs%lut_axes(di(d))

        IF ( ip(d) .GE. axis%values(axis%nvalues) ) THEN
          iw_k(d)  = 0._jprb
        ELSE
!         DO j = 1, axis%nvalues - 1
          DO j = axis%nvalues - 1 ,1 ,-1
            IF ( ip(d) .GT. axis%values(j) ) THEN
              IF ( mfasis_coefs%lut_axes(di(d))%dim_type .EQ. mfasis_dim_opdp .AND. j .NE. 1 ) THEN
                ! optical depths ( lin. int. in log)
!               iw_k(d) =  ip_k(d)/ip(d)                   &
!                 / ( LOG(axis%values(j+1)) - LOG(axis%values(j)) )
                ip_k(d) = ip_k(d) + iw_k(d)/ip(d)                   &
                  / ( LOG(axis%values(j+1)) - LOG(axis%values(j)) )
                iw_k(d) = 0
              ELSE            ! effective radii and alpha, linear interpolation
!               iw_k(d) =  ip_k(d)               &
!                 / ( axis%values(j+1) - axis%values(j) )
                ip_k(d) = ip_k(d) + iw_k(d)               &
                  / ( axis%values(j+1) - axis%values(j) )
                iw_k(d) = 0
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        iw_k(d)  = 0._jprb
      ENDDO

!     ! --------------------------------------------------------------------------
!     ! Fill into the interpolation point array the rest of dimensions
!     ! --------------------------------------------------------------------------

      j = 1
      k = 1
      DO d = 1, nid
        SELECT CASE (mfasis_coefs%lut_axes(di(d))%dim_type)
        CASE (mfasis_dim_opdp)
!         ip_k(d) = od_eff_k(j)
          od_eff_k(j)= od_eff_k(j) + ip_k(d)
          ip_k(d) = 0
          j = j + 1
        CASE (mfasis_dim_effdia)
          ! Check if stored as radius or diamater in LUT
          IF (mfasis_coefs%lut_axes(di(d))%name(1:1) .EQ. "R") THEN
            ed_eff_k(k) = ed_eff_k(k) + ip_k(d)/2
          ELSE
            ed_eff_k(k) = ed_eff_k(k) + ip_k(d)
          ENDIF
          ip_k(d) = 0
          k = k + 1
        END SELECT
      ENDDO

!==================================================================================
!  IF CLOUD ELSE AEROSOL
!==================================================================================
      ! --------------------------------------------------------------------------
      ! Set column weight, optical depth and effective diameters
      ! --------------------------------------------------------------------------
      IF ( mfasis_coefs%file_type .EQ. mfasis_cld ) THEN ! clouds
!=========================!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
        IF(od_eff_v1(2) < 0) od_eff_k(2) = 0._jprb ! Avoid negative optical depths (from mixed-phase correction)


!       od_eff_k(2) = od_eff_k(2) -(fac_sw_k*od2_mixed + fac_sw*od2_mixed_k)
        fac_sw_k    = fac_sw_k    - od2_mixed * od_eff_k(2)
        od2_mixed_k = od2_mixed_k - fac_sw    * od_eff_k(2)
!       od_eff_k(1) = od_eff_k(1) + fac_sw_k*od2_mixed + fac_sw*od2_mixed_k
        fac_sw_k    = fac_sw_k    + od2_mixed * od_eff_k(1)
        od2_mixed_k = od2_mixed_k + fac_sw    * od_eff_k(1)

!       IF(abs(xxx2/opts%dev%x_del2) < 5._jprb) &
!         fac_sw_k= fac_sw_k + (xxx2_k/opts%dev%x_del2) / &
!                    (COSH(xxx2/opts%dev%x_del2)**2) *(1._jprb + TANH(xxx1/opts%dev%x_del1)) /4._jprb
        IF(abs(xxx2/opts%dev%x_del2) < 5._jprb) THEN
          xxx2_k   = xxx2_k + (fac_sw_k/opts%dev%x_del2) / &
                      (COSH(xxx2/opts%dev%x_del2)**2) *(1._jprb + TANH(xxx1/opts%dev%x_del1)) /4._jprb
        ENDIF

!       IF(abs(xxx1/opts%dev%x_del1) < 5._jprb) &
!         fac_sw_k= (xxx1_k/opts%dev%x_del1) / (COSH(xxx1/opts%dev%x_del1)**2) *(1._jprb + TANH(xxx2/opts%dev%x_del2)) /4._jprb
        IF(abs( xxx1/opts%dev%x_del1 ) < 5._jprb) THEN
          xxx1_k   = xxx1_k + (fac_sw_k/opts%dev%x_del1) / &
                      (COSH(xxx1/opts%dev%x_del1)**2) *(1._jprb + TANH(xxx2/opts%dev%x_del2)) /4._jprb
          fac_sw_k = 0
        ENDIF

        fac_sw_k=0._jprb

!       xxx2_k= opts%dev%qq_2*od_eff_k(1)-od_eff_k(2)
        od_eff_k(1) = od_eff_k(1) + opts%dev%qq_2*xxx2_k
        od_eff_k(2) = od_eff_k(2) -      xxx2_k
        xxx2_k      = 0

!       xxx1_k= opts%dev%qq_1*od1_mixed_k-od2_mixed_k
        od1_mixed_k = od1_mixed_k + opts%dev%qq_1*xxx1_k
        od2_mixed_k = od2_mixed_k -      xxx1_k
        xxx1_k= 0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        DO k = 1, npar
          IF ( od_eff_v1(k) .GT. 0. ) THEN
!           ed_eff_k(k) = SUM( od_k(k,:)*ed(k,:) + od(k,:)*ed_k(k,:))/od_eff_v1(k)   &
!                          - ed_eff_v1(k) * od_eff_k(k)/od_eff_v1(k)
            od_k(k,:) = od_k(k,:) + ed(k,:) * ed_eff_k(k) /od_eff_v1(k)
            ed_k(k,:) = ed_k(k,:) + od(k,:) * ed_eff_k(k) /od_eff_v1(k)
            od_eff_k(k) = od_eff_k(k) - ed_eff_v1(k) * ed_eff_k(k) /od_eff_v1(k)
            ed_eff_k(k) = 0._jprb
          ELSE
            ed_eff_k(k) = 0._jprb
          ENDIF
!         od_eff_k(k) = SUM( od_k(k,:))
          od_k(k,:) = od_k(k,:) + od_eff_k(k)
          od_eff_k(k) = 0._jprb
        ENDDO


        od1_int_k = 0
        DO j = nlay, 1, -1
          IF (j > aux%s(prof)%nearestlev_surf - 1) CYCLE    ! Layer is entirely below surface pressure so nothing more to do
!         IF(n_wvdim==3) THEN
!             wvint_bot_tl = wvint_bot_tl + fac_mx(j)          * dvap_tl(j) + fac_mx_tl(j) * dvap(j)
!             wvint_top_tl = wvint_top_tl + (1._jprb-fac_mx(j))* dvap_tl(j) - fac_mx_tl(j) * dvap(j)
!         ENDIF
          IF(n_wvdim==3) THEN
            dvap_k(j)  = dvap_k  (j) + fac_mx(j) * wvint_bot_k 
            fac_mx_k(j)= fac_mx_k(j) + dvap(j)   * wvint_bot_k

            dvap_k(j)  = dvap_k(j)   + (1._jprb-fac_mx(j))* wvint_top_k
            fac_mx_k(j)= fac_mx_k(j) - dvap(j)            * wvint_top_k
          ENDIF

          IF ( ircld%icldarr(cc,j,prof) .NE. 1 ) THEN
            IF(j>1) THEN
              fac_mx_k(j-1) = fac_mx_k(j)
            ELSE
              fac_mx_k(j)=0
            ENDIF
          ELSE !(ircld%icldarr(cc,j,prof) .EQ. 1 )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           ed_k(2,j) = aux_k%ice_dg(j,prof)
            aux_k%ice_dg(j,i) = aux_k%ice_dg(j,i) + ed_k(2,j)
            ed_k(2,j) = 0._jprb
            IF ( coefs%coef_mfasis_cld%clw_scheme == clw_scheme_deff ) THEN ! Deff water clouds scheme
!           ed_k(1,j) = aux_k%clw_dg(j,i)
              aux_k%clw_dg(j,i) = aux_k%clw_dg(j,i) + ed_k(1,j)
              ed_k(1,j) = 0._jprb
            ELSE ! OPAC water cloud parametrization
              ! TODO: take into account optical depth surface correction in eff diameter calculation
              IF ( od(1,j) .GT. 0._jprb ) THEN
!               ed_k(1,j) =  ed_k(1,j)/SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i)) - &
!                             ed(1,j)*SUM(trans_scatt_ir_k%opdpext(1:nwcl_max,j,i))/&
!                             SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i))
                trans_scatt_ir_k%opdpext(1:nwcl_max,j,i)=trans_scatt_ir_k%opdpext(1:nwcl_max,j,i) &
                                                          -ed(1,j)/SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i)) &
                                                           * ed_k(1,j)
                ed_k(1,j) =  ed_k(1,j)/SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i))
              ELSE
                ed_k(1,j) = 0._jprb
              ENDIF
!             ed_k(1,j) = SUM(wcl_opac_deff(1:nwcl_max) * trans_scatt_ir_k%opdpext(1:nwcl_max,j,i))
              trans_scatt_ir_k%opdpext(1:nwcl_max,j,i) = trans_scatt_ir_k%opdpext(1:nwcl_max,j,i) +  &
                                                          wcl_opac_deff(1:nwcl_max) * ed_k(1,j)
              ed_k(1,j) = 0
            ENDIF

!           od2_mixed_k=od2_mixed_k+fac_mx_k(j)*od(2,j)+fac_mx(j)*od_k(2,j)
            fac_mx_k(j) = fac_mx_k(j) + od2_mixed_k * od(2,j)
            od_k(2,j)   = od_k(2,j)   + od2_mixed_k * fac_mx(j)

!           od1_mixed_k=od1_mixed_k+fac_mx_k(j)*od(1,j)+fac_mx(j)*od_k(1,j)
            fac_mx_k(j) = fac_mx_k(j) + od1_mixed_k * od(1,j)
            od_k(1,j)   = od_k(1,j)   + od1_mixed_k * fac_mx(j)

!           IF(abs(ooo1(j)/opts%dev%o_del1) < 5._jprb) &
!           fac_mx_k(j)=(ooo1_k/opts%dev%o_del1) / (COSH(ooo1(j)/opts%dev%o_del1)**2) /2._jprb
            IF(abs( ooo1(j)/opts%dev%o_del1 ) < 5._jprb) THEN
               ooo1_k = ooo1_k + fac_mx_k(j)/opts%dev%o_del1 / (COSH(ooo1(j)/opts%dev%o_del1)**2) /2._jprb
               fac_mx_k(j)= 0
            ENDIF

            fac_mx_k(j)= 0._jprb

!           ooo1_k=od1_int_k
            od1_int_k = od1_int_k + ooo1_k
            ooo1_k    = 0

!           od1_int_k = od1_int_k + od_k(1,j)
            od_k(1,j) = od_k(1,j) + od1_int_k
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF (j == aux%s(prof)%nearestlev_surf - 1) THEN
              ! Modify optical depth in partial layer above surface
!             od_k(:,j) = od_k(:,j) * (1 - aux%s(prof)%pfraction_surf)  &
!                              - od_v1(:,j) * aux_k%s(i)%pfraction_surf
              aux_k%s(i)%pfraction_surf = aux_k%s(i)%pfraction_surf &
                               - SUM(od_v1(:,j) * od_k(:,j))
              od_k(:,j) = od_k(:,j) * (1 - aux%s(prof)%pfraction_surf)
            ENDIF
!           od_k(1,j) = SUM(trans_scatt_ir_k%opdpext(1:5,j,i)) ! water clouds
            trans_scatt_ir_k%opdpext(1:nwcl_max,j,i) = trans_scatt_ir_k%opdpext(1:nwcl_max,j,i) + od_k(1,j)
            od_k(1,j) = 0
!           od_k(2,j) = trans_scatt_ir_k%opdpext(6,j,i)        ! ice clouds
            trans_scatt_ir_k%opdpext(nwcl_max+1,j,i) = trans_scatt_ir_k%opdpext(nwcl_max+1,j,i) + od_k(2,j)
            od_k(2,j) = 0
          ENDIF
        ENDDO

!!!!!!!!!!!!!!!!!!!!!!!
        od1_int_k  = 0._jprb
        od1_mixed_k= 0._jprb
        od2_mixed_k= 0._jprb
        fac_mx_k(:)= 0._jprb
        wvint_bot_k= 0._jprb
        wvint_top_k= 0._jprb
!!!!!!!!!!!!!!!!!!!!!!!
        od_k(:,:) = 0._jprb
        ed_k(:,:) = 0._jprb
!       colwei_k = ircld_k%xcol(cc+1,i) - ircld_k%xcol(cc,i)
        ircld_k%xcol(cc+1,i) = ircld_k%xcol(cc+1,i) + colwei_k
        ircld_k%xcol(cc,i)   = ircld_k%xcol(cc,i)   - colwei_k
        colwei_k = 0
      ELSE  ! Aerosols (not yet tested!!!!)
      ! TODO: Test aerosol functionality

        colwei_cc = 1._jprb     ! single column for aerosols
        od   (:,:)= 0._jprb
        DO k = 1, npar
          par = mfasis_coefs%aer_types(k)
!         od_eff_k(k) = SUM(od_k(k,:))
          od_k(k,:) = od_k(k,:) + od_eff_k(k)
          od_eff_k(k) = 0
          DO j = 1, nlay
            IF (j > aux%s(prof)%nearestlev_surf - 1) EXIT    ! Layer is entirely below surface pressure so nothing more to do
            od(k,j)    = trans_scatt_ir%opdpext(par,j,i)
            IF (j == aux%s(prof)%nearestlev_surf - 1) THEN
              ! Modify optical depth in partial layer above surface
!             od_k(k,j) = od_k(k,j) * (1 - aux%s(prof)%pfraction_surf)  &
!                        - od_v1(k,j) * aux_k%s(i)%pfraction_surf
              aux_k%s(i)%pfraction_surf = aux_k%s(i)%pfraction_surf - od_v1(k,j) * od_k(k,j)
              od_k(k,j) = (1 - aux%s(prof)%pfraction_surf) * od_k(k,j)
            ENDIF
!           od_k(k,j) = trans_scatt_ir_k%opdpext(par,j,i)
            trans_scatt_ir_k%opdpext(par,j,i) = trans_scatt_ir_k%opdpext(par,j,i) + od_k(k,j)
            od_k(k,j) = 0
          ENDDO
        ENDDO
        colwei_k = 0._jprb     ! single column for aerosols A1
        od_k(:,:)= 0._jprb     !                            A2
      ENDIF !ENDIF CLOUD ELSE AEROSOL

    ENDDO ! END of cloud column loop

!====================================================================================
    refl_k     = 0._jprb
    refl_clr_k = 0._jprb

    iw_k(:)  = 0._jprb  ! interpolation weight
    ip_k(:) =  0._jprb         ! interpolation point
    ! --------------------------------------------------------------------------
    ! Albedo
    ! --------------------------------------------------------------------------
!   albedo = MIN(reflectance(i)%refl_out * pi, 1._jprb)
    albedo   =     reflectance(i)%refl_out * pi
    IF (albedo < 1._jprb) THEN
!     albedo_k= reflectance_k(i)%refl_out * pi
      reflectance_k(i)%refl_out = reflectance_k(i)%refl_out + albedo_k * pi
      albedo_k= 0
    ELSE
      albedo_k= 0._jprb
    ENDIF

!   dvap_tl(:)=0._jprb
!   IF(n_wvdim==3) THEN
!     q_mxr(:) =  profiles_int(prof)%q(:)* gas_mass(gas_id_watervapour)/  &
!                 (mair * 1.E06_jprb +profiles_int(prof)%q(:)*gas_mass(gas_id_watervapour))

!     q_mxr_tl(:)= profiles_int_tl(prof)%q(:)*gas_mass(gas_id_watervapour) * (1._jprb - q_mxr(:) ) &
!                 /(mair * 1.E06_jprb +profiles_int(prof)%q(:)*gas_mass(gas_id_watervapour))

!
!     dvap(1:nlay)  = 100._jprb * (profiles(prof)%p(2:nlev)-profiles(prof)%p(1:nlev-1))* &
!                                          0.5_jprb * (q_mxr(2:nlev)+q_mxr(1:nlev-1))/gravity
!     dvap_tl(1:nlay)=100._jprb *((profiles_tl(prof)%p(2:nlev)-profiles_tl(prof)%p(1:nlev-1))* &
!                                          0.5_jprb * (q_mxr(2:nlev)+q_mxr(1:nlev-1))/gravity  &
!                                 + (profiles(prof)%p(2:nlev)-profiles(prof)%p(1:nlev-1))* &
!                                    0.5_jprb * (q_mxr_tl(2:nlev)+q_mxr_tl(1:nlev-1))/gravity  &
!                                )  
!
!     jj = aux%s(prof)%nearestlev_surf - 1
!     IF(jj>=1 .AND. jj<=nlay) THEN
!       dvap_tl(jj)=dvap_tl(jj)*(1._jprb - aux   %s(prof)%pfraction_surf)           &
!                       - dvap   (jj) * aux_tl%s(prof)%pfraction_surf
!       dvap   (jj)=dvap   (jj)*(1._jprb - aux   %s(prof)%pfraction_surf)
!     ENDIF

    q_mxr_k(:) = 0
    IF(n_wvdim==3) THEN

      jj = aux%s(prof)%nearestlev_surf - 1
      IF(jj>=1 .AND. jj<=nlay) THEN
        aux_k%s(i)%pfraction_surf =               aux_k%s(i)%pfraction_surf  &
                                         - dvap_1(jj) * dvap_k(jj) 

        dvap_k(jj)=dvap_k(jj)*(1._jprb - aux   %s(prof)%pfraction_surf)  
      ENDIF

!------------------------------------------------------------------------------------------------------
!       IF (opts%interpolation%lgradp) THEN
!         dvap_tl(1:nlay)=100._jprb *((profiles_tl(prof)%p(2:nlev)-profiles_tl(prof)%p(1:nlev-1))* &
!                                              0.5_jprb * (q_mxr(2:nlev)+q_mxr(1:nlev-1))/gravity  &
!                                     + (profiles(prof)%p(2:nlev)-profiles(prof)%p(1:nlev-1))* &
!                                        0.5_jprb * (q_mxr_tl(2:nlev)+q_mxr_tl(1:nlev-1))/gravity  &
!                                    )
!      ELSE
!        dvap_tl(1:nlay)=100._jprb *((profiles(prof)%p(2:nlev)-profiles(prof)%p(1:nlev-1))* &
!                                        0.5_jprb * (q_mxr_tl(2:nlev)+q_mxr_tl(1:nlev-1))/gravity  &
!                                    )
!      ENDIF

      array1(0)      = 0
      array1(nlay+1) = 0

      IF (opts%interpolation%lgradp) THEN
        array1(1:nlay) = 100._jprb * 0.5_jprb * (q_mxr(2:nlev)+q_mxr(1:nlev-1))/gravity * dvap_k(1:nlay)
        profiles_k(i)%p(1:nlev)  = profiles_k(i)%p(1:nlev) + array1(0:nlay) - array1(1:nlay+1)
      ENDIF

      array1(1:nlay) = 100._jprb * 0.5_jprb * (profiles(prof)%p(2:nlev)-profiles(prof)%p(1:nlev-1))/gravity &
                          * dvap_k(1:nlay)
      q_mxr_k(1:nlev) = q_mxr_k(1:nlev) +  array1(0:nlay) + array1(1:nlay+1)
!------------------------------------------------------------------------------------------------------

      profiles_int_k(i)%q(:) = profiles_int_k(i)%q(:) +  &
                        q_mxr_k(:)* gas_mass(gas_id_watervapour) * (1._jprb - q_mxr(:) ) & 
                    /(mair * 1.E06_jprb +profiles_int(prof)%q(:)*gas_mass(gas_id_watervapour))
      q_mxr_k(:) = 0
    ENDIF
    dvap_k(:)=0._jprb

  ENDDO ! end of channel loop

  ! --------------------------------------------------------------------------
  ! Tidy up
  ! --------------------------------------------------------------------------

  DEALLOCATE( ip, ip_k, iw_k, di, od, od_k, od_eff, od_eff_k, od_v1, od_eff_v1, STAT=err)
  THROWM(err.NE.0,"Deallocation of memory for rttov_mfasis_k failed")
  IF (mfasis_coefs%file_type .EQ. mfasis_cld) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   DEALLOCATE(ed, ed_k, ed_eff, ed_eff_v1, ed_aux, ed_eff_k, fac_mx, fac_mx_k, STAT=err)
    DEALLOCATE(ed, ed_k, ed_eff, ed_eff_v1, ed_aux, ed_eff_k, fac_mx, fac_mx_k, &
               q_mxr, q_mxr_k, dvap, dvap_1, dvap_k, array1, STAT=err)
    THROWM(err.NE.0,"Deallocation of memory for rttov_mfasis_k failed")
  ENDIF


  CATCH

END SUBROUTINE rttov_mfasis_k
