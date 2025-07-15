! Description:
!> @file
!!   TL of MFASIS fast visible/near-IR scattering model.
!
!> @brief
!!   TL of MFASIS fast visible/near-IR scattering model.
!!
!!
!! @param[out]    err               status on exit
!! @param[in]     chanprof          specifies channels and profiles to simulate
!! @param[in]     chanflag          flags to indicate which channels with LUT available
!! @param[in]     opts              options to configure the simulations
!! @param[in]     profiles          input atmospheric profiles and surface variables
!! @param[in]     profiles_tl       input profile perturbations
!! @param[in]     profiles_int      profiles in internal units
!! @param[in]     profiles_int_tl   profile perturbations in internal units
!! @param[in]     coefs             coefficients structure for instrument to simulate
!! @param[in]     ircld             information on cloud columns
!! @param[in]     ircld_tl          cloud column perturbations
!! @param[in]     aux               additional internal profile variables
!! @param[in]     aux_tl            additional internal profile variable perturbations
!! @param[in]     reflectance       surface BRDFs
!! @param[in]     reflectance_tl    surface BRDF perturbations
!! @param[in]     solar_spectrum    TOA solar irradiance for each channel
!! @param[in]     trans_scatt_ir    cloud/aerosol optical depths
!! @param[in]     trans_scatt_ir_tl cloud/aerosol optical depth perturbations
!! @param[in]     mfasis_refl       quantities computed by rttov_mfasis used by TL/AD/K
!! @param[in,out] radiance_tl       output radiance and BRF perturbations
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
SUBROUTINE rttov_mfasis_tl( &
              err,              &
              chanprof,         &
              chanflag,         &
              opts,             &
              profiles,         &
              profiles_tl,      &
              profiles_int,     &
              profiles_int_tl,  &
              coefs,            &
              ircld,            &
              ircld_tl,         &
              aux,              &
              aux_tl,           &
              reflectance,      &
              reflectance_tl,   &
              solar_spectrum,   &
              trans_scatt_ir,   &
              trans_scatt_ir_tl,&
              mfasis_refl,      &
              radiance_tl)

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
  TYPE(rttov_profile),                INTENT(IN)    :: profiles_tl(:)
  TYPE(rttov_profile),                INTENT(IN)    :: profiles_int(:)
  TYPE(rttov_profile),                INTENT(IN)    :: profiles_int_tl(:)
  TYPE(rttov_coefs),                  INTENT(IN)    :: coefs
  TYPE(rttov_ircld),                  INTENT(IN)    :: ircld
  TYPE(rttov_ircld),                  INTENT(IN)    :: ircld_tl
  TYPE(rttov_profile_aux),            INTENT(IN)    :: aux
  TYPE(rttov_profile_aux),            INTENT(IN)    :: aux_tl
  TYPE(rttov_reflectance),            INTENT(IN)    :: reflectance(SIZE(chanprof))
  TYPE(rttov_reflectance),            INTENT(IN)    :: reflectance_tl(SIZE(chanprof))
  REAL(jprb),                         INTENT(IN)    :: solar_spectrum(SIZE(chanprof))
  TYPE(rttov_transmission_scatt_ir),  INTENT(IN)    :: trans_scatt_ir
  TYPE(rttov_transmission_scatt_ir),  INTENT(IN)    :: trans_scatt_ir_tl
  TYPE(rttov_mfasis_refl),            INTENT(IN)    :: mfasis_refl(0:,:)
  TYPE(rttov_radiance),               INTENT(INOUT) :: radiance_tl
!INTF_END

#include "rttov_errorreport.interface"

  INTEGER(jpim)              :: nchanprof, prof, ncolms, npar, chan
  INTEGER(jpim)              :: nlay, nlev, ndim, nid
  INTEGER(jpim)              :: d, n, i, j, k, cc, par, jj
  TYPE(rttov_coef_mfasis)    :: mfasis_coefs
  REAL(jprb),    ALLOCATABLE :: ip   (:)
  REAL(jprb),    ALLOCATABLE :: ip_tl(:), iw_tl(:)

  REAL(jprb),    ALLOCATABLE :: q_mxr(:),dvap(:)
  REAL(jprb),    ALLOCATABLE :: q_mxr_tl(:),dvap_tl(:)
  REAL(jprb)                 :: wvint_bot, wvint_top 
  REAL(jprb)                 :: wvint_bot_tl, wvint_top_tl 
! REAL(jprb)                 :: iw_wv(3)
  REAL(jprb)                 :: iw_wv_tl(3)
  INTEGER(jpim)              :: n_wvdim

  INTEGER(jpim), ALLOCATABLE :: di(:)
  REAL(jprb),    ALLOCATABLE :: od   (:,:), od_eff   (:), ed   (:,:), ed_eff(:), ed_aux(:)
  REAL(jprb),    ALLOCATABLE :: od_v1(:,:), od_eff_v1(:), ed_eff_v1(:) !versions of nonlinear quantities needed for TL/AD/K
  REAL(jprb),    ALLOCATABLE :: od_tl(:,:), od_eff_tl(:), ed_tl(:,:), ed_eff_tl(:)

  REAL(jprb)                 :: albedo
  REAL(jprb)                 :: albedo_tl
  TYPE(rttov_mfasis_axis)    :: axis
  REAL(jprb)                 :: colwei_clr, colwei_cc
  REAL(jprb)                 :: refl_tl, rflcol_tl, refl_clr_tl, colwei_tl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(jprb)                 :: od1_int,od1_mixed,od2_mixed
  REAL(jprb)                 :: od1_int_tl,od1_mixed_tl,od2_mixed_tl, ooo1_tl
  REAL(jprb),    ALLOCATABLE :: ooo1(:)
  REAL(jprb),    ALLOCATABLE :: fac_mx(:)
  REAL(jprb),    ALLOCATABLE :: fac_mx_tl(:)

  REAL(jprb)                 :: fac_sw,xxx1,xxx2            ! variables    used for switching on/off mixed cloud correction
  REAL(jprb)                 :: fac_sw_tl,xxx1_tl,xxx2_tl   ! variables    used for switching on/off mixed cloud correction
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
            ip_tl  (nid) , &      ! interpolation point
            iw_tl  (nid) , &      ! interpolation weight
            di     (nid), STAT=err )  ! original dimension index in nid arrays
  THROWM(err.NE.0,"Allocation of memory for rttov_mfasis_tl failed")


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
            od_eff_tl (npar)    , & ! Effective total optical depth per type of particle
            od     (npar, nlay) , & ! Effective optical depth per type of particle and layer
            od_v1  (npar, nlay), &     ! od saved for adjoint computations
            od_tl  (npar, nlay), STAT=err)  ! Effective optical depth per type of particle and layer
  THROWM(err.NE.0,"Allocation of memory for rttov_mfasis_tl failed")

  ! For the time being we assume effective diameters will be out of LUT for aerosols
  IF ( mfasis_coefs%file_type .EQ. 1 ) THEN
    ALLOCATE( ed_eff (npar)      , & ! Total effective diameter per type of particle
              ed_eff_v1(npar) ,  &      ! version of ed_eff for adjoint computations
              ed_aux (npar)      , & ! Effective diameters to be used with 0 opdp to avoid flagging
              ed_eff_tl(npar)    , & ! Total effective diameter per type of particle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              fac_mx (nlay),      &         ! fraction of cloud considered as (potentially) mixed cloud
              fac_mx_tl (nlay),   &         ! fraction of cloud considered as (potentially) mixed cloud
              ooo1 (nlay),        &         ! discriminator for switching on mixed layer summation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
              ed   (npar, nlay)  , & ! Effective diameter per type of particle and layer
              ed_tl(npar, nlay),   &
              q_mxr(nlev),         &
              q_mxr_tl(nlev),         &
              dvap (nlay),         &
              dvap_tl (nlay),         &
                           STAT=err) !
    THROWM(err.NE.0,"Allocation of memory for rttov_mfasis_tl failed")

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
  ! --------------------------------------------------------------------------
  ! Channel loop
  ! --------------------------------------------------------------------------
  DO i = 1, nchanprof
    IF (.NOT. chanflag(i)) CYCLE
    prof = chanprof(i)%prof
    chan = chanprof(i)%chan

    n_wvdim=size(mfasis_coefs%lut(chan)%qint,2)
    dvap(:)=0._jprb
    dvap_tl(:)=0._jprb
    IF(n_wvdim==3) THEN
      q_mxr(:) =  profiles_int(prof)%q(:)* gas_mass(gas_id_watervapour)/  &
                  (mair * 1.E06_jprb +profiles_int(prof)%q(:)*gas_mass(gas_id_watervapour))

      q_mxr_tl(:)= profiles_int_tl(prof)%q(:)*gas_mass(gas_id_watervapour) * (1._jprb - q_mxr(:) ) &
                  /(mair * 1.E06_jprb +profiles_int(prof)%q(:)*gas_mass(gas_id_watervapour))

      dvap(1:nlay)  = 100._jprb * (profiles(prof)%p(2:nlev)-profiles(prof)%p(1:nlev-1))* &
                                           0.5_jprb * (q_mxr(2:nlev)+q_mxr(1:nlev-1))/gravity


       IF (opts%interpolation%lgradp) THEN
         dvap_tl(1:nlay)=100._jprb *((profiles_tl(prof)%p(2:nlev)-profiles_tl(prof)%p(1:nlev-1))* &
                                              0.5_jprb * (q_mxr(2:nlev)+q_mxr(1:nlev-1))/gravity  &
                                     + (profiles(prof)%p(2:nlev)-profiles(prof)%p(1:nlev-1))* &
                                        0.5_jprb * (q_mxr_tl(2:nlev)+q_mxr_tl(1:nlev-1))/gravity  &
                                    )
      ELSE
        dvap_tl(1:nlay)=100._jprb *((profiles(prof)%p(2:nlev)-profiles(prof)%p(1:nlev-1))* &
                                        0.5_jprb * (q_mxr_tl(2:nlev)+q_mxr_tl(1:nlev-1))/gravity  &
                                    )
      ENDIF
 
      jj = aux%s(prof)%nearestlev_surf - 1
      IF(jj>=1 .AND. jj<=nlay) THEN
        dvap_tl(jj)=dvap_tl(jj)*(1._jprb - aux   %s(prof)%pfraction_surf)           &
                        - dvap   (jj) * aux_tl%s(prof)%pfraction_surf
        dvap   (jj)=dvap   (jj)*(1._jprb - aux   %s(prof)%pfraction_surf)
      ENDIF
    ENDIF

    colwei_clr = ircld%xcolclr(prof) ! clear column
    ! --------------------------------------------------------------------------
    ! Albedo
    ! --------------------------------------------------------------------------
!   albedo = MIN(reflectance(i)%refl_out * pi, 1._jprb)
    albedo   =     reflectance(i)%refl_out * pi
    IF(albedo < 1._jprb) THEN
      albedo_tl= reflectance_tl(i)%refl_out * pi
    ELSE
      albedo   = 1._jprb
      albedo_tl= 0._jprb
    ENDIF


    ip   (:) = -1._jprb          ! interpolation point
    ip_tl(:) =  0._jprb          ! interpolation point


!   iw   (:)  = 0._jprb  ! interpolation weight
    iw_tl(:)  = 0._jprb  ! interpolation weight


    ! --------------------------------------------------------------------------
    ! Loop over (cloud) columns
    ! --------------------------------------------------------------------------
    ! No cloud columns for aerosols calculations
    IF ( mfasis_coefs%file_type == mfasis_cld ) THEN
      ncolms = ircld%ncolumn(prof)
    ELSE
      ncolms = 1_jpim
    ENDIF


    refl_tl     = 0._jprb
    refl_clr_tl = 0._jprb

    DO cc = 1, ncolms
!==================================================================================
! First compute nonlinear quatities needed in linear computations (for given channel "i" and column "cc")
!==================================================================================
      ! --------------------------------------------------------------------------
      ! Set column weight, optical depth and effective diameters
      ! --------------------------------------------------------------------------
      IF ( mfasis_coefs%file_type .EQ. mfasis_cld ) THEN ! clouds
        colwei_cc = ircld%xcol(cc+1,prof) - ircld%xcol(cc,prof)
!       colwei_tl = ircld_tl%xcol(cc+1,prof) - ircld_tl%xcol(cc,prof)
        od   (:,:) = 0._jprb
!       od_tl(:,:) = 0._jprb
        ed   (:,:) = 0._jprb
!       ed_tl(:,:) = 0._jprb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        od1_int = 0._jprb
        od1_mixed= 0._jprb
        od2_mixed= 0._jprb
        fac_mx(:)= 0._jprb
        wvint_bot=0._jprb
        wvint_top=0._jprb


!       od1_int_tl = 0._jprb
!       od1_mixed_tl= 0._jprb
!       od2_mixed_tl= 0._jprb
!       fac_mx_tl(:)= 0._jprb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO j = 1, nlay
          IF (j > aux%s(prof)%nearestlev_surf - 1) EXIT    ! Layer is entirely below surface pressure so nothing more to do
          IF ( ircld%icldarr(cc,j,prof) .EQ. 1 ) THEN
            od   (1,j) = SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i)) ! water clouds
!           od_tl(1,j) = SUM(trans_scatt_ir_tl%opdpext(1:nwcl_max,j,i)) ! water clouds
            od   (2,j) = trans_scatt_ir%opdpext(nwcl_max+1,j,i)        ! ice clouds
!           od_tl(2,j) = trans_scatt_ir_tl%opdpext(nwcl_max+1,j,i)        ! ice clouds
            IF (j == aux%s(prof)%nearestlev_surf - 1) THEN
              ! Modify optical depth in partial layer above surface
              od_v1(:,j) = od   (:,j)
!             od_tl(:,j) = od_tl(:,j) * (1 - aux%s(prof)%pfraction_surf)  &
!                              - od_v1(:,j) * aux_tl%s(prof)%pfraction_surf
              od   (:,j) = od   (:,j) * (1 - aux%s(prof)%pfraction_surf)
            ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            od1_int = od1_int + od(1,j)
!           od1_int_tl = od1_int_tl + od_tl(1,j)
            ooo1(j)=od1_int-opts%dev%od1_thresh
!           ooo1_tl=od1_int_tl
            fac_mx(j)=(1._jprb + TANH(ooo1(j)/opts%dev%o_del1))/2._jprb

!           fac_mx_tl(j)= 0._jprb
!           IF(abs(ooo1(j)/opts%dev%o_del1) < 5._jprb) &
!           fac_mx_tl(j)=(ooo1_tl/opts%dev%o_del1) / (COSH(ooo1(j)/opts%dev%o_del1)**2) /2._jprb

            od1_mixed=od1_mixed+fac_mx(j)*od(1,j)
            od2_mixed=od2_mixed+fac_mx(j)*od(2,j)
!           od1_mixed_tl=od1_mixed_tl+fac_mx_tl(j)*od(1,j)+fac_mx(j)*od_tl(1,j)
!           od2_mixed_tl=od2_mixed_tl+fac_mx_tl(j)*od(2,j)+fac_mx(j)*od_tl(2,j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF ( coefs%coef_mfasis_cld%clw_scheme == clw_scheme_deff ) THEN ! Deff water clouds scheme
              ed   (1,j) = aux%clw_dg(j,prof)
!             ed_tl(1,j) = aux_tl%clw_dg(j,prof)
            ELSE ! OPAC water cloud parametrization
              ed   (1,j) = SUM(wcl_opac_deff(1:nwcl_max) * trans_scatt_ir%opdpext(1:nwcl_max,j,i))
!             ed_tl(1,j) = SUM(wcl_opac_deff(1:nwcl_max) * trans_scatt_ir_tl%opdpext(1:nwcl_max,j,i))
              IF ( od(1,j) .GT. 0. ) THEN
                ! TODO: take into account optical depth surface correction in eff diameter calculation
                ed   (1,j) =  ed   (1,j)/SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i))
!               ed_tl(1,j) =  ed_tl(1,j)/SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i)) - &
!                             ed(1,j)*SUM(trans_scatt_ir_tl%opdpext(1:nwcl_max,j,i))/&
!                             SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i))
              ELSE
                ed   (1,j) = 0._jprb
!               ed_tl(1,j) = 0._jprb
              ENDIF
            ENDIF
            ed   (2,j) = aux%ice_dg(j,prof)
!           ed_tl(2,j) = aux_tl%ice_dg(j,prof)
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
!         od_eff_tl(k) = SUM( od_tl(k,:))
          IF ( od_eff_v1(k) .GT. 0. ) THEN
            ed_eff   (k) = SUM( od(k,:)*ed(k,:))/od_eff_v1(k)
            ed_eff_v1(k) = ed_eff(k)
!           ed_eff_tl(k) = SUM( od_tl(k,:)*ed(k,:) + od(k,:)*ed_tl(k,:))/od_eff_v1(k)   &
!                          - ed_eff_v1(k) * od_eff_tl(k)/od_eff_v1(k)

          ELSE
            ed_eff   (k) = ed_aux(k)
!           ed_eff_tl(k) = 0._jprb
          ENDIF
        ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        xxx1= opts%dev%qq_1*od1_mixed-od2_mixed
        xxx2= opts%dev%qq_2*od_eff(1)-od_eff(2)
!       xxx1_tl= opts%dev%qq_1*od1_mixed_tl-od2_mixed_tl
!       xxx2_tl= opts%dev%qq_2*od_eff_tl(1)-od_eff_tl(2)
        fac_sw= (1._jprb + TANH(xxx1/opts%dev%x_del1))*(1._jprb + TANH(xxx2/opts%dev%x_del2))/4._jprb

!       fac_sw_tl=0._jprb
!       IF(abs(xxx1/opts%dev%x_del1) < 5._jprb) &
!         fac_sw_tl= (xxx1_tl/opts%dev%x_del1) / &
!                    (COSH(xxx1/opts%dev%x_del1)**2) *(1._jprb + TANH(xxx2/opts%dev%x_del2)) /4._jprb
!       IF(abs(xxx2/opts%dev%x_del2) < 5._jprb) &
!         fac_sw_tl= fac_sw_tl +       &
!                    (xxx2_tl/opts%dev%x_del2) / &
!                    (COSH(xxx2/opts%dev%x_del2)**2) *(1._jprb + TANH(xxx1/opts%dev%x_del1)) /4._jprb


        od_eff(1) = od_eff(1) + fac_sw*od2_mixed
        od_eff(2) = od_eff(2) - fac_sw*od2_mixed
        od_eff_v1(2) = od_eff(2)
        IF(od_eff_v1(2) < 0) od_eff(2) = 0._jprb ! Avoid negative optical depths (from mixed-phase correction above)
!       od_eff_tl(1) = od_eff_tl(1) + fac_sw_tl*od2_mixed + fac_sw*od2_mixed_tl
!       od_eff_tl(2) = od_eff_tl(2) -(fac_sw_tl*od2_mixed + fac_sw*od2_mixed_tl)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ELSE  ! aerosols
        colwei_cc = 1._jprb     ! single column for aerosols
!       colwei_tl = 0._jprb     ! single column for aerosols
        od   (:,:)= 0._jprb
!       od_tl(:,:)= 0._jprb
        DO k = 1, npar
          par = mfasis_coefs%aer_types(k)
          DO j = 1, nlay
            IF (j > aux%s(prof)%nearestlev_surf - 1) EXIT    ! Layer is entirely below surface pressure so nothing more to do
            od(k,j)    = trans_scatt_ir%opdpext(par,j,i)
!           od_tl(k,j) = trans_scatt_ir_tl%opdpext(par,j,i)
            IF (j == aux%s(prof)%nearestlev_surf - 1) THEN
              ! Modify optical depth in partial layer above surface
              od_v1(:,j) = od   (:,j)
!             od_tl(k,j) = od_tl(k,j) * (1 - aux%s(prof)%pfraction_surf)  &
!                         - od_v1(k,j) * aux_tl%s(prof)%pfraction_surf
              od   (k,j) = od   (k,j) * (1 - aux%s(prof)%pfraction_surf)
            ENDIF
          ENDDO
          od_eff   (k) = SUM(od   (k,:))
!         od_eff_tl(k) = SUM(od_tl(k,:))
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
!         ip_tl(d) = od_eff_tl(j)
          j = j + 1
        CASE (mfasis_dim_effdia)
          ! Check IF stored as radius or diamater in LUT
          IF (mfasis_coefs%lut_axes(di(d))%name(1:1) .EQ. "R") THEN
            ip(d) = ed_eff(k)/2
!           ip_tl(d) = ed_eff_tl(k)/2
          ELSE
            ip(d) = ed_eff(k)
!           ip_tl(d) = ed_eff_tl(k)
          ENDIF
          k = k + 1
        END SELECT
      ENDDO


!==================================================================================
! Start linear computations for channel "i" and column "cc"
!==================================================================================
      ! --------------------------------------------------------------------------
      ! Set column weight, optical depth and effective diameters
      ! --------------------------------------------------------------------------
      IF ( mfasis_coefs%file_type .EQ. mfasis_cld ) THEN ! clouds
!       colwei_cc = ircld%xcol(cc+1,prof) - ircld%xcol(cc,prof)
        colwei_tl = ircld_tl%xcol(cc+1,prof) - ircld_tl%xcol(cc,prof)
!       od   (:,:) = 0._jprb
        od_tl(:,:) = 0._jprb
!       ed   (:,:) = 0._jprb
        ed_tl(:,:) = 0._jprb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       od1_int = 0._jprb
!       od1_mixed= 0._jprb
!       od2_mixed= 0._jprb
!       fac_mx(:)= 0._jprb

        od1_int_tl = 0._jprb
        od1_mixed_tl= 0._jprb
        od2_mixed_tl= 0._jprb
        fac_mx_tl(:)= 0._jprb
        wvint_bot_tl=0._jprb
        wvint_top_tl=0._jprb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DO j = 1, nlay
          IF (j > aux%s(prof)%nearestlev_surf - 1) EXIT    ! Layer is entirely below surface pressure so nothing more to do
          IF ( ircld%icldarr(cc,j,prof) .EQ. 1 ) THEN
!           od   (1,j) = SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i)) ! water clouds
            od_tl(1,j) = SUM(trans_scatt_ir_tl%opdpext(1:nwcl_max,j,i)) ! water clouds
!           od   (2,j) = trans_scatt_ir%opdpext(nwcl_max+1,j,i)        ! ice clouds
            od_tl(2,j) = trans_scatt_ir_tl%opdpext(nwcl_max+1,j,i)        ! ice clouds
            IF (j == aux%s(prof)%nearestlev_surf - 1) THEN
              ! Modify optical depth in partial layer above surface
              od_tl(:,j) = od_tl(:,j) * (1 - aux%s(prof)%pfraction_surf)  &
                               - od_v1(:,j) * aux_tl%s(prof)%pfraction_surf
!             od   (:,j) = od   (:,j) * (1 - aux%s(prof)%pfraction_surf)
            ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           od1_int = od1_int + od(1,j)
            od1_int_tl = od1_int_tl + od_tl(1,j)
!           ooo1(j)=od1_int-opts%dev%od1_thresh
            ooo1_tl=od1_int_tl
!           fac_mx(j)=(1._jprb + TANH(ooo1(j)/opts%dev%o_del1))/2._jprb

            fac_mx_tl(j)= 0._jprb
            IF(abs(ooo1(j)/opts%dev%o_del1) < 5._jprb) &
            fac_mx_tl(j)=(ooo1_tl/opts%dev%o_del1) / (COSH(ooo1(j)/opts%dev%o_del1)**2) /2._jprb

!             od1_mixed=od1_mixed+fac_mx(j)*od(1,j)
!             od2_mixed=od2_mixed+fac_mx(j)*od(2,j)
              od1_mixed_tl=od1_mixed_tl+fac_mx_tl(j)*od(1,j)+fac_mx(j)*od_tl(1,j)
              od2_mixed_tl=od2_mixed_tl+fac_mx_tl(j)*od(2,j)+fac_mx(j)*od_tl(2,j)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            IF ( coefs%coef_mfasis_cld%clw_scheme == clw_scheme_deff ) THEN ! Deff water clouds scheme
!             ed   (1,j) = aux%clw_dg(j,prof)
              ed_tl(1,j) = aux_tl%clw_dg(j,prof)
            ELSE ! OPAC water cloud parametrization
!             ed   (1,j) = SUM(wcl_opac_deff(1:nwcl_max) * trans_scatt_ir%opdpext(1:nwcl_max,j,i))
              ed_tl(1,j) = SUM(wcl_opac_deff(1:nwcl_max) * trans_scatt_ir_tl%opdpext(1:nwcl_max,j,i))
              IF ( od(1,j) .GT. 0. ) THEN
                ! TODO: take into account optical depth surface correction in eff diameter calculation
!               ed   (1,j) =  ed   (1,j)/SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i))
                ed_tl(1,j) =  ed_tl(1,j)/SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i)) - &
                              ed(1,j)*SUM(trans_scatt_ir_tl%opdpext(1:nwcl_max,j,i))/&
                              SUM(trans_scatt_ir%opdpext(1:nwcl_max,j,i))
              ELSE
!               ed   (1,j) = 0._jprb
                ed_tl(1,j) = 0._jprb
              ENDIF
            ENDIF
!           ed   (2,j) = aux%ice_dg(j,prof)
            ed_tl(2,j) = aux_tl%ice_dg(j,prof)
          ELSEIF(j>1) THEN
            fac_mx_tl(j) = fac_mx_tl(j-1)
          ENDIF
          IF(n_wvdim==3) THEN
              wvint_bot_tl = wvint_bot_tl + fac_mx(j)          * dvap_tl(j) + fac_mx_tl(j) * dvap(j)
              wvint_top_tl = wvint_top_tl + (1._jprb-fac_mx(j))* dvap_tl(j) - fac_mx_tl(j) * dvap(j)
          ENDIF

        ENDDO

        DO k = 1, npar
!         od_eff   (k) = SUM( od   (k,:))
!         od_eff_v1(k) = od_eff   (k)
          od_eff_tl(k) = SUM( od_tl(k,:))
          IF ( od_eff_v1(k) .GT. 0. ) THEN
!           ed_eff   (k) = SUM( od(k,:)*ed(k,:))/od_eff_v1(k)
!           ed_eff_v1(k) = ed_eff(k)
            ed_eff_tl(k) = SUM( od_tl(k,:)*ed(k,:) + od(k,:)*ed_tl(k,:))/od_eff_v1(k)   &
                           - ed_eff_v1(k) * od_eff_tl(k)/od_eff_v1(k)
          ELSE
!           ed_eff   (k) = ed_aux(k)
            ed_eff_tl(k) = 0._jprb
          ENDIF
        ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       xxx1= opts%dev%qq_1*od1_mixed-od2_mixed
!       xxx2= opts%dev%qq_2*od_eff(1)-od_eff(2)
        xxx1_tl= opts%dev%qq_1*od1_mixed_tl-od2_mixed_tl
        xxx2_tl= opts%dev%qq_2*od_eff_tl(1)-od_eff_tl(2)
!       fac_sw= (1._jprb + TANH(xxx1/opts%dev%x_del1))*(1._jprb + TANH(xxx2/opts%dev%x_del2))/4._jprb

        fac_sw_tl=0._jprb
        IF(abs(xxx1/opts%dev%x_del1) < 5._jprb) &
          fac_sw_tl= (xxx1_tl/opts%dev%x_del1) / &
                     (COSH(xxx1/opts%dev%x_del1)**2) *(1._jprb + TANH(xxx2/opts%dev%x_del2)) /4._jprb
        IF(abs(xxx2/opts%dev%x_del2) < 5._jprb) &
          fac_sw_tl= fac_sw_tl +       &
                     (xxx2_tl/opts%dev%x_del2) / &
                     (COSH(xxx2/opts%dev%x_del2)**2) *(1._jprb + TANH(xxx1/opts%dev%x_del1)) /4._jprb

!       od_eff(1) = od_eff(1) + fac_sw*od2_mixed
!       od_eff(2) = od_eff(2) - fac_sw*od2_mixed
        od_eff_tl(1) = od_eff_tl(1) + fac_sw_tl*od2_mixed + fac_sw*od2_mixed_tl
        od_eff_tl(2) = od_eff_tl(2) -(fac_sw_tl*od2_mixed + fac_sw*od2_mixed_tl)
        IF(od_eff_v1(2) < 0) od_eff_tl(2) = 0._jprb ! Avoid negative optical depths (from mixed-phase correction above)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ELSE  ! aerosols
!       colwei_cc = 1._jprb     ! single column for aerosols
        colwei_tl = 0._jprb     ! single column for aerosols
!       od   (:,:)= 0._jprb
        od_tl(:,:)= 0._jprb
        DO k = 1, npar
          par = mfasis_coefs%aer_types(k)
          DO j = 1, nlay
            IF (j > aux%s(prof)%nearestlev_surf - 1) EXIT    ! Layer is entirely below surface pressure so nothing more to do
!           od(k,j)    = trans_scatt_ir%opdpext(par,j,i)
            od_tl(k,j) = trans_scatt_ir_tl%opdpext(par,j,i)
            IF (j == aux%s(prof)%nearestlev_surf - 1) THEN
              ! Modify optical depth in partial layer above surface
!             od_v1(k,j) = od   (k,j)
              od_tl(k,j) = od_tl(k,j) * (1 - aux%s(prof)%pfraction_surf)  &
                          - od_v1(k,j) * aux_tl%s(prof)%pfraction_surf
!             od   (k,j) = od   (k,j) * (1 - aux%s(prof)%pfraction_surf)
            ENDIF
          ENDDO
!         od_eff   (k) = SUM(od   (k,:))
          od_eff_tl(k) = SUM(od_tl(k,:))
        ENDDO
      ENDIF !IF CLOUD ELSE AEROSOL
      ! --------------------------------------------------------------------------
      ! Fill into the interpolation point array the rest of dimensions
      ! --------------------------------------------------------------------------
      j = 1
      k = 1
      DO d = 1, nid
        SELECT CASE (mfasis_coefs%lut_axes(di(d))%dim_type)
        CASE (mfasis_dim_opdp)
!         ip   (d) = od_eff   (j)
          ip_tl(d) = od_eff_tl(j)
          j = j + 1
        CASE (mfasis_dim_effdia)
          ! Check IF stored as radius or diamater in LUT
          IF (mfasis_coefs%lut_axes(di(d))%name(1:1) .EQ. "R") THEN
!           ip(d) = ed_eff(k)/2
            ip_tl(d) = ed_eff_tl(k)/2
          ELSE
!           ip(d) = ed_eff(k)
            ip_tl(d) = ed_eff_tl(k)
          ENDIF
          k = k + 1
        END SELECT
      ENDDO

      ! --------------------------------------------------------------------------
      ! Look for nearest entries in LUT and compute interpolation weights
      ! --------------------------------------------------------------------------
      ! Here we are assuming the LUT dimensions for optical depths (and effective
      ! diameter in the case of clouds) are always ordered as aer_types in the case
      ! of aerosols and first water and THEN ice in the case of clouds
      DO d = 1, nid
        IF (mfasis_coefs%lut_axes(di(d))%dim_type .EQ. mfasis_dim_scaangle) CYCLE ! Skip scattering angle (already DOne)
        axis = mfasis_coefs%lut_axes(di(d))
        ! Reset lc_aux and iw to make sure correct for every column calculation
!       iw(d)     = 0._jprb  ! interpolation weight
        iw_tl(d)  = 0._jprb

        IF ( ip(d) .GE. axis%values(axis%nvalues) ) THEN
!         iw   (d)  = 1._jprb
          iw_tl(d)  = 0._jprb
        ELSE
          DO j = 1, axis%nvalues - 1
            IF ( ip(d) .GT. axis%values(j) ) THEN
              IF ( mfasis_coefs%lut_axes(di(d))%dim_type .EQ. mfasis_dim_opdp .AND. j .NE. 1 ) THEN
                ! optical depths ( lin. int. in log)
!               iw(d) = ( LOG(ip(d)) - LOG(axis%values(j)) ) &
!                 / ( LOG(axis%values(j+1)) - LOG(axis%values(j)) )
                iw_tl(d) =  ip_tl(d)/ip(d)                   &
                  / ( LOG(axis%values(j+1)) - LOG(axis%values(j)) )
              ELSE            ! effective radii and alpha, linear interpolation
!               iw(d) = ( ip(d) - axis%values(j) ) &
!                 / ( axis%values(j+1) - axis%values(j) )
                iw_tl(d) =  ip_tl(d)               &
                  / ( axis%values(j+1) - axis%values(j) )
              ENDIF
            ENDIF
          ENDDO
        ENDIF

      ENDDO
      iw_wv_tl(:) = 0._jprb
      if( n_wvdim==3) THEN
          iw_wv_tl(2)=                  wvint_bot_tl  / &
               (mfasis_coefs%lut(chan)%qint(1,2) - mfasis_coefs%lut(chan)%qint(1,1))
          iw_wv_tl(3)=                  wvint_top_tl  / &
               (mfasis_coefs%lut(chan)%qint(2,3) - mfasis_coefs%lut(chan)%qint(2,1))
          iw_wv_tl(1) =       - iw_wv_tl(2) - iw_wv_tl(3)
       ENDIF


      ! --------------------------------------------------------------------------
      ! Linear interpolation in nid dimensions (replaces subroutine rttov_mfasis_interpolate)
      ! --------------------------------------------------------------------------

      rflcol_tl=0
      DO d = 1, nid
        IF ( mfasis_coefs%lut_axes(di(d))%dim_type .EQ. mfasis_dim_scaangle ) THEN
          rflcol_tl = rflcol_tl+mfasis_refl(cc,i)%refl_lin_coef(d)*albedo_tl
        ELSE
          rflcol_tl = rflcol_tl+mfasis_refl(cc,i)%refl_lin_coef(d)*iw_tl(d)
        ENDIF
      ENDDO
      rflcol_tl = rflcol_tl+sum(iw_wv_tl(1:n_wvdim)*mfasis_refl(cc,i)%refl_wv(1:n_wvdim))


      ! --------------------------------------------------------------------------
      ! Add column contribution to total reflectance
      ! --------------------------------------------------------------------------

      refl_tl = refl_tl + rflcol_tl*colwei_cc + mfasis_refl(cc,i)%refl*colwei_tl

    ENDDO ! END of cloud column loop

    ! --------------------------------------------------------------------------
    ! Compute clear reflectance and add clear column to total reflectabce
    ! --------------------------------------------------------------------------
    ! Rerun linear interpolation in nid dimensions setting optical depths
    ! to zero with weight 1
    ! This assumes zero optical depth is always the first value
    IF(n_wvdim==3) THEN
      wvint_bot_tl=0._jprb
      wvint_top_tl=0._jprb
      DO j = 1, nlay
        wvint_top_tl = wvint_top_tl +  dvap_tl(j)
      ENDDO
    ENDIF
    iw_wv_tl(:) = 0._jprb
    if( n_wvdim==3) THEN
        iw_wv_tl(2)=                  wvint_bot_tl  / &
             (mfasis_coefs%lut(chan)%qint(1,2) - mfasis_coefs%lut(chan)%qint(1,1))
        iw_wv_tl(3)=                  wvint_top_tl  / &
             (mfasis_coefs%lut(chan)%qint(2,3) - mfasis_coefs%lut(chan)%qint(2,1))
        iw_wv_tl(1) =       - iw_wv_tl(2) - iw_wv_tl(3)
     ENDIF

    DO d = 1, nid
      IF ( mfasis_coefs%lut_axes(di(d))%dim_type .EQ. mfasis_dim_opdp ) THEN
!       iw   (d) = 1._jprb
        iw_tl(d) = 0._jprb
      ENDIF
    ENDDO

    refl_clr_tl=0
    DO d = 1, nid
      IF ( mfasis_coefs%lut_axes(di(d))%dim_type .EQ. mfasis_dim_scaangle ) THEN
        refl_clr_tl = refl_clr_tl+mfasis_refl(0,i)%refl_lin_coef(d)*albedo_tl
      ELSE
        refl_clr_tl = refl_clr_tl+mfasis_refl(0,i)%refl_lin_coef(d)*iw_tl(d)
      ENDIF
    ENDDO
    refl_clr_tl = refl_clr_tl+sum(iw_wv_tl(1:n_wvdim)*mfasis_refl(0,i)%refl_wv(1:n_wvdim))


!   colwei_clr = ircld%xcolclr(prof) ! clear column
    colwei_tl = ircld_tl%xcolclr(prof) ! clear column

    refl_tl = refl_tl + refl_clr_tl*colwei_clr + mfasis_refl(0,i)%refl*colwei_tl

    ! --------------------------------------------------------------------------
    ! Fill radiance structure
    ! --------------------------------------------------------------------------
    radiance_tl%refl(i) = refl_tl
    radiance_tl%refl_clear(i) = refl_clr_tl

    radiance_tl%total(i) = refl_tl * solar_spectrum(i) * pi_r  &
      * COS(profiles(prof)%sunzenangle * deg2rad)
    radiance_tl%clear(i) = refl_clr_tl * solar_spectrum(i) * pi_r  &
      * COS(profiles(prof)%sunzenangle * deg2rad)
    radiance_tl%cloudy(i) = radiance_tl%total(i)

  ENDDO ! END of channel loop

  ! --------------------------------------------------------------------------
  ! Tidy up
  ! --------------------------------------------------------------------------

! DEALLOCATE( ip, ip_tl, iw_tl, di, od, od_tl, od_eff, od_eff_tl, STAT=err)
  DEALLOCATE( ip, ip_tl, iw_tl, di, od, od_tl, od_eff, od_eff_tl, od_v1, od_eff_v1, STAT=err)
  THROWM(err.NE.0,"Deallocation of memory for rttov_mfasis_tl failed")
  IF (mfasis_coefs%file_type .EQ. mfasis_cld) THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   DEALLOCATE(ed, ed_tl, ed_eff, ed_aux, ed_eff_tl, fac_mx, fac_mx_tl, STAT=err)
    DEALLOCATE(ed, ed_tl, ed_eff, ed_eff_v1, ed_aux, ed_eff_tl, fac_mx, fac_mx_tl, &
               q_mxr, q_mxr_tl, dvap, dvap_tl, STAT=err)

    THROWM(err.NE.0,"Deallocation of memory for rttov_mfasis_tl failed")
  ENDIF

  CATCH

END SUBROUTINE rttov_mfasis_tl
