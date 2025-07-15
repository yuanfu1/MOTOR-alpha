! Description:
!> @file
!!   Subroutines for interface to HTFRTC
!
!> @brief
!!   Subroutines for interface to HTFRTC
!!
!! @details
!!   The htfrtc_interface subroutine implements the HTFRTC forward and
!!   Jacobian models. The Jacobian model is run if one of the optional
!!   profiles_k_pc or profiles_k_rec arguments are present.
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
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
MODULE rttov_htfrtc_interface_mod

  USE parkind1, ONLY: jpim, jprb, jplm

  IMPLICIT NONE

#include "rttov_errorreport.interface"

  PRIVATE
  PUBLIC :: htfrtc_interface

CONTAINS

  SUBROUTINE htfrtc_interface(err, &
                              coefs, &
                              opts, &
                              profiles, &
                              pccomp, &
                              calcemis, &
                              emissivity, &
                              emissivity_k, &
                              profiles_k_pc, &
                              profiles_k_rec)
#include "throw.h"

  USE rttov_types, ONLY : &
      rttov_coefs, &
      rttov_options, &
      rttov_profile, &
      rttov_pccomp, &
      rttov_emissivity

  USE rttov_const, ONLY : &
      deg2rad, gravity, earthradius, mair, rgc, &
      planck_c1, planck_c2, &
      surftype_land, surftype_sea, surftype_seaice, &
      gas_unit_ppmvdry,   &
      gas_unit_specconc,  &
      gas_unit_ppmv,      &
      gas_id_watervapour, &
      gas_id_ozone,       &
      gas_id_co2,         &
      gas_id_n2o,         &
      gas_id_co,          &
      gas_id_ch4,         &
      gas_id_so2,         &
      gas_mass

    IMPLICIT NONE

! input/output
  INTEGER(jpim)         , INTENT(OUT)              :: err
  TYPE(rttov_coefs)     , INTENT(IN)               :: coefs
  TYPE(rttov_options)   , INTENT(IN)               :: opts
  TYPE(rttov_profile)   , INTENT(IN)               :: profiles(:)
  TYPE(rttov_pccomp)    , INTENT(INOUT)            :: pccomp
  LOGICAL(jplm)         , INTENT(IN)    , OPTIONAL :: calcemis(:)
  TYPE(rttov_emissivity), INTENT(INOUT) , OPTIONAL :: emissivity(:)
  TYPE(rttov_emissivity), INTENT(INOUT) , OPTIONAL :: emissivity_k(:)
  TYPE(rttov_profile)   , INTENT(INOUT) , OPTIONAL :: profiles_k_pc(:)
  TYPE(rttov_profile)   , INTENT(INOUT) , OPTIONAL :: profiles_k_rec(:)

  LOGICAL(jplm) :: do_k, do_lambertian, user_emis
  LOGICAL(jplm) :: laddf(coefs%coef_htfrtc%n_gas_l,coefs%coef_htfrtc%n_f)
  LOGICAL(jplm) :: laddch(coefs%coef_htfrtc%n_gas_l,coefs%coef_htfrtc%n_ch)

  INTEGER(jpim),PARAMETER :: ns=9
  INTEGER(jpim) :: n_p,n_f,n_b,n_lt,n_prof
  INTEGER(jpim) :: nlayers,nlevels
  INTEGER(jpim) :: n_pc,n_pc_oc,n_ch
  INTEGER(jpim) :: surf_level
  INTEGER(jpim) :: cloud_level
  INTEGER(jpim) :: i,i_f,i_p,i_pc,i_ch
  INTEGER(jpim) :: ja,jb,js
  INTEGER(jpim) :: t_b_surf_i
  INTEGER(jpim) :: sg
  INTEGER(jpim) :: igu,ngu,igc,ngc,is,j

  INTEGER(jpim),DIMENSION(profiles(1)%nlayers) :: p_opt_ly_i
  INTEGER(jpim),DIMENSION(profiles(1)%nlevels) :: t_b_i
  INTEGER(jpim),DIMENSION(coefs%coef_htfrtc%n_gas_l) :: gid,gmid

  REAL(jprb),DIMENSION(profiles(1)%nlevels) :: z
  REAL(jprb),DIMENSION(profiles(1)%nlevels) :: t_b_r,b
  REAL(jprb),DIMENSION(profiles(1)%nlevels) :: rdown,rup,rdownlb
  REAL(jprb),DIMENSION(profiles(1)%nlevels) :: rup_cld
  REAL(jprb),DIMENSION(profiles(1)%nlevels) :: dbdt
  REAL(jprb),DIMENSION(profiles(1)%nlevels) :: cv
  REAL(jprb),DIMENSION(profiles(1)%nlevels) :: totdown,totup,totdownlb

  REAL(jprb),DIMENSION(profiles(1)%nlayers) :: p_ly,dp,dz
  REAL(jprb),DIMENSION(profiles(1)%nlayers) :: t_ly,t2_ly,q_ly,q2_ly
  REAL(jprb),DIMENSION(profiles(1)%nlayers) :: p_opt_ly_r
  REAL(jprb),DIMENSION(profiles(1)%nlayers) :: blydown,blyup,blydownlb
  REAL(jprb),DIMENSION(profiles(1)%nlayers) :: lt,ltlb
  REAL(jprb),DIMENSION(profiles(1)%nlayers,coefs%coef_htfrtc%n_gas_l) :: x_ly
  REAL(jprb),DIMENSION(profiles(1)%nlayers) :: tau,tausl,tausllb,tr,trlb
  REAL(jprb),DIMENSION(profiles(1)%nlayers) :: gravm1,mpath,mft,mftlb
  REAL(jprb),DIMENSION(profiles(1)%nlayers) :: dblydowndt,dblyupdt,dblydowndtlb

  REAL(jprb),DIMENSION(profiles(1)%nlayers,coefs%coef_htfrtc%n_gas_l) :: dx
  REAL(jprb),DIMENSION(profiles(1)%nlayers) :: dttmp,dxtmp
  REAL(jprb),DIMENSION(profiles(1)%nlayers) :: dtrdt
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_gas_l,profiles(1)%nlayers) :: &
  dtrdx
  REAL(jprb),DIMENSION(profiles(1)%nlayers) :: dtrdtlb
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_gas_l,profiles(1)%nlayers) :: &
  dtrdxlb
  REAL(jprb),DIMENSION(profiles(1)%nlayers) :: drdowndt
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_gas_l,profiles(1)%nlayers) :: &
  drdowndx
  REAL(jprb),DIMENSION(profiles(1)%nlayers) :: drupdt
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_gas_l,profiles(1)%nlayers) :: &
  drupdx

  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f) :: rc,normm1
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f,profiles(1)%nlevels) :: drcdt
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f,coefs%coef_htfrtc%n_gas_l, &
  profiles(1)%nlevels) :: drcdx
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f,ns) :: drcds
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f) :: drcdsem
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f) :: &
  surf_em,aems,bems,cems,expf
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f) :: &
  daems,dbems,dcems,dexpf
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f) :: &
  dsurf_em,dsurf_emdu,dsurf_emdv,dsurf_emdskint
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f) :: sp
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f) :: swrad_f
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f) :: wn_f,wn3_f
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f) :: planck_c1_wn3_f,planck_c2_wn_f
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f) :: r_f,bt_f

  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_ch) :: wn,wn3
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_ch) :: planck_c1_wn3,planck_c2_wn
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_ch) :: swrad_tmp

  REAL(jprb),DIMENSION(profiles(1)%nlayers,coefs%coef_htfrtc%n_f) :: rc_oc
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f,profiles(1)%nlayers) :: rc_oc_t
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_f) :: rc_cld

  REAL(jprb),DIMENSION(profiles(1)%nlayers,coefs%coef_htfrtc%n_pc_oc) :: &
  overcast_pcscores
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_pc_oc,profiles(1)%nlayers) :: &
  overcast_pcscores_t

  REAL(jprb),DIMENSION(profiles(1)%nlevels,coefs%coef_htfrtc%n_pc) :: &
  profiles_k_pc_t
  REAL(jprb),DIMENSION(profiles(1)%nlevels,coefs%coef_htfrtc%n_gas_l, &
  coefs%coef_htfrtc%n_pc) :: profiles_k_pc_x
  REAL(jprb),DIMENSION(ns,coefs%coef_htfrtc%n_pc) :: profiles_k_pc_s

  REAL(jprb),DIMENSION(profiles(1)%nlevels,coefs%coef_htfrtc%n_ch) :: &
  profiles_k_rec_t
  REAL(jprb),DIMENSION(profiles(1)%nlevels,coefs%coef_htfrtc%n_gas_l, &
  coefs%coef_htfrtc%n_ch) :: profiles_k_rec_x
  REAL(jprb),DIMENSION(coefs%coef_htfrtc%n_ch,9) :: profiles_k_rec_s

  REAL(jprb) :: t_b_min,t_b_inc,t_b_inc_m1
  REAL(jprb) :: lt_min,lt_inc,lt_inc_m1
  REAL(jprb) :: surf_wind,surf_wind2,surf_windm1
  REAL(jprb) :: surf_layer_frac,t_b_surf_r
  REAL(jprb) :: cloud_layer_frac
  REAL(jprb) :: rdown_surf,b_surf,rdown_surflb
  REAL(jprb) :: db_surfdt
  REAL(jprb) :: tot
  REAL(jprb) :: cf
  REAL(jprb) :: qscv
  REAL(jprb) :: e_wvp

  TRY

  n_ch=coefs%coef_htfrtc%n_ch
  nlevels=profiles(1)%nlevels
  nlayers=profiles(1)%nlayers
  n_f=coefs%coef_htfrtc%n_f
  n_p=coefs%coef_htfrtc%n_p
  n_b=coefs%coef_htfrtc%n_b
  n_lt=coefs%coef_htfrtc%n_lt
  n_prof=size(profiles)
  do_k=(PRESENT(profiles_k_pc).OR.PRESENT(profiles_k_rec))

  IF (PRESENT(calcemis) .AND. .NOT. PRESENT(emissivity)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0,'emissivity arg is mandatory if calcemis arg is present')
  ENDIF
  IF (PRESENT(calcemis)) THEN
    IF (SIZE(calcemis) /= n_f * n_prof) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'calcemis arg must have size npredictors*nprofiles')
    ENDIF
  ENDIF
  IF (PRESENT(emissivity)) THEN
    IF (SIZE(emissivity) /= n_f * n_prof) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'emissivity arg must have size npredictors*nprofiles')
    ENDIF
  ENDIF
  IF (PRESENT(emissivity_k)) THEN
    IF (SIZE(emissivity_k) /= n_f * n_prof) THEN
      err = errorstatus_fatal
      THROWM(err.NE.0,'emissivity_k arg must have size npredictors*nprofiles')
    ENDIF
  ENDIF

  n_pc=opts%htfrtc_opts%n_pc_in
  IF (opts%htfrtc_opts%n_pc_in < 1 .OR. opts%htfrtc_opts%n_pc_in > coefs%coef_htfrtc%n_pc) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0,'invalid number of PC scores in opts%htfrtc_opts%n_pc_in')
  ENDIF
  IF (n_pc * n_prof > SIZE(pccomp%total_pcscores)) THEN
    err = errorstatus_fatal
    THROWM(err.NE.0,'pccomp structure is too small for specified opts%htfrtc_opts%n_pc_in')
  ENDIF
  n_pc_oc=min(n_pc,coefs%coef_htfrtc%n_pc_oc)

  x_ly=0.0_jprb

  laddf=(coefs%coef_htfrtc%addf==1)
  laddch=(coefs%coef_htfrtc%addch==1)

  t_b_min=coefs%coef_htfrtc%val_b(1)
  t_b_inc=coefs%coef_htfrtc%val_b(2)-coefs%coef_htfrtc%val_b(1)
  t_b_inc_m1=1.0/t_b_inc

  lt_min=coefs%coef_htfrtc%val_lt(1)
  lt_inc=coefs%coef_htfrtc%val_lt(2)-coefs%coef_htfrtc%val_lt(1)
  lt_inc_m1=1.0/lt_inc

  DO i_f=1,n_f
    normm1(i_f)=1.0/coefs%coef_htfrtc%val_norm(i_f)
  ENDDO

  IF(opts%htfrtc_opts%overcast) THEN
    pccomp%overcast_pcscores=0.0_jprb
  ENDIF

  !Reconstruct for radiances / brightness temperatures, precalc
  IF(opts%htfrtc_opts%reconstruct) THEN
    FORALL(i=1:n_ch)
      wn(i)=coefs%coef_htfrtc%sensor_freq(i)
      wn3(i)=wn(i)**3
      planck_c1_wn3(i)=planck_c1*wn3(i)
      planck_c2_wn(i)=planck_c2*wn(i)
    ENDFORALL
  ENDIF

  IF(PRESENT(emissivity_k).AND.opts%rt_all%switchrad) THEN
    FORALL(i=1:n_f)
      wn_f(i)=coefs%coef_htfrtc%freq(i)
      wn3_f(i)=wn_f(i)**3
      planck_c1_wn3_f(i)=planck_c1*wn3_f(i)
      planck_c2_wn_f(i)=planck_c2*wn_f(i)
    ENDFORALL
  ENDIF

  !Introduction of the RTTOV variable, optional, trace gases
  !"supergas" (sg) type = RTTOV predictor type

  ngu=0
    ngu=ngu+1
    gid(ngu)=1
    gmid(ngu)=gas_id_watervapour
  IF (opts%rt_all%ozone_data) THEN
    ngu=ngu+1
    gid(ngu)=2
    gmid(ngu)=gas_id_ozone
  ENDIF
  IF (opts%rt_all%co2_data) THEN
    ngu=ngu+1
    gid(ngu)=3
    gmid(ngu)=gas_id_co2
  ENDIF
  IF (opts%rt_all%n2o_data) THEN
    ngu=ngu+1
    gid(ngu)=4
    gmid(ngu)=gas_id_n2o
  ENDIF
  IF (opts%rt_all%co_data) THEN
    ngu=ngu+1
    gid(ngu)=5
    gmid(ngu)=gas_id_co
  ENDIF
  IF (opts%rt_all%ch4_data) THEN
    ngu=ngu+1
    gid(ngu)=6
    gmid(ngu)=gas_id_ch4
  ENDIF
  IF (opts%rt_all%so2_data) THEN
    ngu=ngu+1
    gid(ngu)=7
    gmid(ngu)=gas_id_so2
  ENDIF

  sg=9
  IF(.NOT.(opts%rt_all%n2o_data.OR.opts%rt_all%co_data.OR. &
  opts%rt_all%ch4_data.OR.opts%rt_all%so2_data)) THEN
    sg=8
    IF(.NOT.(opts%rt_all%co2_data)) THEN
      sg=7
    ENDIF
  ENDIF

  ngc=ngu
  IF (sg>=7) THEN
    IF (.NOT.(opts%rt_all%ozone_data)) THEN
      ngc=ngc+1
      gid(ngc)=2
    ENDIF
  ENDIF
  IF (sg>=8) THEN
    IF (.NOT.(opts%rt_all%co2_data)) THEN
      ngc=ngc+1
      gid(ngc)=3
    ENDIF
  ENDIF
  IF (sg>=9) THEN
    IF (.NOT.(opts%rt_all%n2o_data)) THEN
      ngc=ngc+1
      gid(ngc)=4
    ENDIF
    IF (.NOT.(opts%rt_all%co_data)) THEN
      ngc=ngc+1
      gid(ngc)=5
    ENDIF
    IF (.NOT.(opts%rt_all%ch4_data)) THEN
      ngc=ngc+1
      gid(ngc)=6
    ENDIF
    IF (.NOT.(opts%rt_all%so2_data)) THEN
      ngc=ngc+1
      gid(ngc)=7
    ENDIF
  ENDIF

  DO igc=8,10
    IF (sg+1==igc) THEN
      ngc=ngc+1
      gid(ngc)=igc
      x_ly(1:nlayers,ngc)=coefs%coef_htfrtc%mixed_ref_frac(1,igc)
    ENDIF
  ENDDO

  !Diffuse downwelling slant path
  mftlb=coefs%coef_htfrtc%mftlb(1)

  e_wvp=gas_mass(gas_id_watervapour)/mair

  !Main loop over profiles
  DO i_p=1,n_prof

    !User levels to layers
    FORALL(i=1:nlayers)
      p_ly(i)=0.5*(profiles(i_p)%p(i)+profiles(i_p)%p(i+1))
      dp(i)=profiles(i_p)%p(i+1)-profiles(i_p)%p(i)
      gravm1(i)=1.0/gravity
      mpath(i)=100.0*dp(i)*gravm1(i)
      t_ly(i)=0.5*(profiles(i_p)%t(i)+profiles(i_p)%t(i+1))
      dz(i)=dp(i)*gravm1(i)*rgc*t_ly(i)/mair/p_ly(i)
      mft(i)=1.0/cos(profiles(i_p)%zenangle*deg2rad)
    ENDFORALL

    !Map user levels and ref prof / coeff levels
    js=1
    DO jb=1,nlayers
      IF (p_ly(jb)<coefs%coef_htfrtc%p(1)) THEN
        p_opt_ly_i(jb)=1
        p_opt_ly_r(jb)=0.0
      ELSE IF (p_ly(jb)>=coefs%coef_htfrtc%p(n_p)) THEN
        p_opt_ly_i(jb)=n_p-1
        p_opt_ly_r(jb)=1.0
      ELSE
        DO ja=js,n_p-1
        IF (p_ly(jb)>=coefs%coef_htfrtc%p(ja) .AND. &
           p_ly(jb)<coefs%coef_htfrtc%p(ja+1)) THEN
          js=ja
          p_opt_ly_i(jb)=ja
          p_opt_ly_r(jb)=(p_ly(jb)-coefs%coef_htfrtc%p(ja))/ &
          (coefs%coef_htfrtc%p(ja+1)-coefs%coef_htfrtc%p(ja))
        ENDIF
        ENDDO
      ENDIF
    ENDDO

    !Find the level just above the surface
    surf_level=nlevels-1
    surf_layer_frac=1.0
    DO jb=1,nlevels-1
      IF ((profiles(i_p)%s2m%p.GT.profiles(i_p)%p(jb)).AND. &
         (profiles(i_p)%s2m%p.LE.profiles(i_p)%p(jb+1))) THEN
        surf_level=jb
        surf_layer_frac=(profiles(i_p)%s2m%p-profiles(i_p)%p(jb))/ &
        (profiles(i_p)%p(jb+1)-profiles(i_p)%p(jb))
      ENDIF
    ENDDO

    !Near surface layer
    p_ly(surf_level)=0.5*(profiles(i_p)%p(surf_level)+profiles(i_p)%s2m%p)
    dp(surf_level)=profiles(i_p)%s2m%p-profiles(i_p)%p(surf_level)
    gravm1(surf_level)=1.0/gravity
    mpath(surf_level)=100.0*dp(surf_level)*gravm1(surf_level)
    t_ly(surf_level)=0.5*(profiles(i_p)%t(surf_level)+profiles(i_p)%s2m%t)
    dz(surf_level)=dp(surf_level)*gravm1(surf_level)*rgc*t_ly(surf_level) &
    /mair/p_ly(surf_level)
    z(surf_level+1)=profiles(i_p)%elevation
    DO i=surf_level,1,-1
      z(i)=z(i+1)+dz(i)
    ENDDO
    IF(.NOT.opts%rt_all%plane_parallel) THEN
      FORALL(i=1:surf_level)
        mft(i)=1.0/sqrt(1.0-(sin(profiles(i_p)%zenangle*deg2rad)*earthradius &
        /(earthradius+0.5*(z(i)+z(i+1))))**2)
      ENDFORALL
    ENDIF

    !Planckians, precalc
    FORALL(i=1:nlevels)
      t_b_i(i)=min(n_b-1,max(1,1+int(t_b_inc_m1*(profiles(i_p)%t(i)-t_b_min))))
      t_b_r(i)=t_b_inc_m1*((profiles(i_p)%t(i)-t_b_min)-t_b_inc*(t_b_i(i)-1))
    ENDFORALL
    t_b_i(surf_level+1)=min(n_b-1, &
    max(1,1+int(t_b_inc_m1*(profiles(i_p)%s2m%t-t_b_min))))
    t_b_r(surf_level+1)=t_b_inc_m1 &
    *((profiles(i_p)%s2m%t-t_b_min)-t_b_inc*(t_b_i(surf_level+1)-1))
    t_b_surf_i=min(n_b-1,max(1,1+ &
    int(t_b_inc_m1*(profiles(i_p)%skin%t-t_b_min))))
    t_b_surf_r=t_b_inc_m1* &
    ((profiles(i_p)%skin%t-t_b_min)-t_b_inc*(t_b_surf_i-1))

    !Unit conversions
    SELECT CASE(profiles(1)%gas_units)
    CASE(gas_unit_ppmvdry)
      FORALL(i=1:nlevels)
        cv(i)=1.0_jprb/(Mair*1.E06_jprb+ &
        gas_mass(gas_id_watervapour)*profiles(i_p)%q(i))
      ENDFORALL
      FORALL(i=1:nlayers)
        x_ly(i,1)=gas_mass(gas_id_watervapour) &
        *0.5*(profiles(i_p)%q(i)*cv(i)+profiles(i_p)%q(i+1)*cv(i+1))
      ENDFORALL
      IF(opts%rt_all%use_q2m) THEN
        qscv=1.0_jprb/(Mair*1.E06_jprb+ &
        gas_mass(gas_id_watervapour)*profiles(i_p)%s2m%q)
        x_ly(surf_level,1)=gas_mass(gas_id_watervapour) &
        *0.5*(profiles(i_p)%q(surf_level)*cv(surf_level)+ &
        profiles(i_p)%s2m%q*qscv)
      ELSE
        x_ly(surf_level,1)=gas_mass(gas_id_watervapour)* &
        (profiles(i_p)%q(surf_level)*cv(surf_level) &
        +surf_layer_frac*(profiles(i_p)%q(surf_level+1)*cv(surf_level+1) &
        -profiles(i_p)%q(surf_level)*cv(surf_level)))
      ENDIF
    CASE(gas_unit_ppmv)
      FORALL(i=1:nlevels)
        cv(i)=1.0_jprb/(1.0e6_jprb*Mair+profiles(i_p)%q(i)* &
        (gas_mass(gas_id_watervapour)-Mair))
      ENDFORALL
      FORALL(i=1:nlayers)
        x_ly(i,1)=gas_mass(gas_id_watervapour) &
        *0.5*(profiles(i_p)%q(i)*cv(i)+profiles(i_p)%q(i+1)*cv(i+1))
      ENDFORALL
      IF(opts%rt_all%use_q2m) THEN
        qscv=1.0_jprb/(1.0e6_jprb*Mair+profiles(i_p)%s2m%q* &
        (gas_mass(gas_id_watervapour)-Mair))
        x_ly(surf_level,1)=gas_mass(gas_id_watervapour) &
        *0.5*(profiles(i_p)%q(surf_level)*cv(surf_level)+ &
        profiles(i_p)%s2m%q*qscv)
      ELSE
        x_ly(surf_level,1)=gas_mass(gas_id_watervapour)* &
        (profiles(i_p)%q(surf_level)*cv(surf_level) &
        +surf_layer_frac*(profiles(i_p)%q(surf_level+1)*cv(surf_level+1) &
        -profiles(i_p)%q(surf_level)*cv(surf_level)))
      ENDIF
    END SELECT

    SELECT CASE(profiles(1)%gas_units)
    CASE(gas_unit_ppmvdry,gas_unit_ppmv)
      DO igu=1,ngu
      SELECT CASE(gid(igu))
      CASE(2)
        FORALL(i=1:nlayers)
          x_ly(i,igu)=0.5*(profiles(i_p)%o3(i)*cv(i)+profiles(i_p)%o3(i+1) &
          *cv(i+1))*gas_mass(gas_id_ozone)
        ENDFORALL
        x_ly(surf_level,igu)=(profiles(i_p)%o3(surf_level)*cv(surf_level) &
        +surf_layer_frac*(profiles(i_p)%o3(surf_level+1)*cv(surf_level+1) &
        -profiles(i_p)%o3(surf_level)*cv(surf_level))) &
        *gas_mass(gas_id_ozone)
      CASE(3)
        FORALL(i=1:nlayers)
          x_ly(i,igu)=0.5*(profiles(i_p)%co2(i)*cv(i)+profiles(i_p)%co2(i+1) &
          *cv(i+1))*gas_mass(gas_id_co2)
        ENDFORALL
        x_ly(surf_level,igu)=(profiles(i_p)%co2(surf_level)*cv(surf_level) &
        +surf_layer_frac*(profiles(i_p)%co2(surf_level+1)*cv(surf_level+1) &
        -profiles(i_p)%co2(surf_level)*cv(surf_level))) &
        *gas_mass(gas_id_co2)
      CASE(4)
        FORALL(i=1:nlayers)
          x_ly(i,igu)=0.5*(profiles(i_p)%n2o(i)*cv(i)+profiles(i_p)%n2o(i+1) &
          *cv(i+1))*gas_mass(gas_id_n2o)
        ENDFORALL
        x_ly(surf_level,igu)=(profiles(i_p)%n2o(surf_level)*cv(surf_level) &
        +surf_layer_frac*(profiles(i_p)%n2o(surf_level+1)*cv(surf_level+1) &
        -profiles(i_p)%n2o(surf_level)*cv(surf_level))) &
        *gas_mass(gas_id_n2o)
      CASE(5)
        FORALL(i=1:nlayers)
          x_ly(i,igu)=0.5*(profiles(i_p)%co(i)*cv(i)+profiles(i_p)%co(i+1) &
          *cv(i+1))*gas_mass(gas_id_co)
        ENDFORALL
        x_ly(surf_level,igu)=(profiles(i_p)%co(surf_level)*cv(surf_level) &
        +surf_layer_frac*(profiles(i_p)%co(surf_level+1)*cv(surf_level+1) &
        -profiles(i_p)%co(surf_level)*cv(surf_level))) &
        *gas_mass(gas_id_co)
      CASE(6)
        FORALL(i=1:nlayers)
          x_ly(i,igu)=0.5*(profiles(i_p)%ch4(i)*cv(i)+profiles(i_p)%ch4(i+1) &
          *cv(i+1))*gas_mass(gas_id_ch4)
        ENDFORALL
        x_ly(surf_level,igu)=(profiles(i_p)%ch4(surf_level)*cv(surf_level) &
        +surf_layer_frac*(profiles(i_p)%ch4(surf_level+1)*cv(surf_level+1) &
        -profiles(i_p)%ch4(surf_level)*cv(surf_level))) &
        *gas_mass(gas_id_ch4)
      CASE(7)
        FORALL(i=1:nlayers)
          x_ly(i,igu)=0.5*(profiles(i_p)%so2(i)*cv(i)+profiles(i_p)%so2(i+1) &
          *cv(i+1))*gas_mass(gas_id_so2)
        ENDFORALL
        x_ly(surf_level,igu)=(profiles(i_p)%so2(surf_level)*cv(surf_level) &
        +surf_layer_frac*(profiles(i_p)%so2(surf_level+1)*cv(surf_level+1) &
        -profiles(i_p)%so2(surf_level)*cv(surf_level))) &
        *gas_mass(gas_id_so2)
      END SELECT
      ENDDO
    CASE(gas_unit_specconc)
      FORALL(i=1:nlayers)
        x_ly(i,1)=0.5*(profiles(i_p)%q(i)+profiles(i_p)%q(i+1))
      ENDFORALL
      IF(opts%rt_all%use_q2m) THEN
        x_ly(surf_level,1)=0.5*(profiles(i_p)%q(surf_level)+profiles(i_p)%s2m%q)
      ELSE
        x_ly(surf_level,1)=profiles(i_p)%q(surf_level) &
        +surf_layer_frac*(profiles(i_p)%q(surf_level+1) &
        -profiles(i_p)%q(surf_level))
      ENDIF
      DO igu=1,ngu
      SELECT CASE(gid(igu))
      CASE(2)
        FORALL(i=1:nlayers)
          x_ly(i,igu)=0.5*(profiles(i_p)%o3(i)+profiles(i_p)%o3(i+1))
        ENDFORALL
        x_ly(surf_level,igu)=profiles(i_p)%o3(surf_level) &
        +surf_layer_frac*(profiles(i_p)%o3(surf_level+1) &
        -profiles(i_p)%o3(surf_level))
      CASE(3)
        FORALL(i=1:nlayers)
          x_ly(i,igu)=0.5*(profiles(i_p)%co2(i)+profiles(i_p)%co2(i+1))
        ENDFORALL
        x_ly(surf_level,igu)=profiles(i_p)%co2(surf_level) &
        +surf_layer_frac*(profiles(i_p)%co2(surf_level+1) &
        -profiles(i_p)%co2(surf_level))
      CASE(4)
        FORALL(i=1:nlayers)
          x_ly(i,igu)=0.5*(profiles(i_p)%n2o(i)+profiles(i_p)%n2o(i+1))
        ENDFORALL
        x_ly(surf_level,igu)=profiles(i_p)%n2o(surf_level) &
        +surf_layer_frac*(profiles(i_p)%n2o(surf_level+1) &
        -profiles(i_p)%n2o(surf_level))
      CASE(5)
        FORALL(i=1:nlayers)
          x_ly(i,igu)=0.5*(profiles(i_p)%co(i)+profiles(i_p)%co(i+1))
        ENDFORALL
        x_ly(surf_level,igu)=profiles(i_p)%co(surf_level) &
        +surf_layer_frac*(profiles(i_p)%co(surf_level+1) &
        -profiles(i_p)%co(surf_level))
      CASE(6)
        FORALL(i=1:nlayers)
          x_ly(i,igu)=0.5*(profiles(i_p)%ch4(i)+profiles(i_p)%ch4(i+1))
        ENDFORALL
        x_ly(surf_level,igu)=profiles(i_p)%ch4(surf_level) &
        +surf_layer_frac*(profiles(i_p)%ch4(surf_level+1) &
        -profiles(i_p)%ch4(surf_level))
      CASE(7)
        FORALL(i=1:nlayers)
          x_ly(i,igu)=0.5*(profiles(i_p)%so2(i)+profiles(i_p)%so2(i+1))
        ENDFORALL
        x_ly(surf_level,igu)=profiles(i_p)%so2(surf_level) &
        +surf_layer_frac*(profiles(i_p)%so2(surf_level+1) &
        -profiles(i_p)%so2(surf_level))
      END SELECT
      ENDDO
    END SELECT

    !Ref prof where no user prof
    FORALL(i=1:nlayers)
      cv(i)=1.0_jprb/(Mair*1.E06_jprb+ &
      gas_mass(gas_id_watervapour) * profiles(i_p)%q(i))
    ENDFORALL
    DO igc=ngu+1,ngc
    SELECT CASE(gid(igc))
    CASE(2)
      FORALL(i=1:surf_level)
        x_ly(i,igc)=((1.0-p_opt_ly_r(i))* &
        coefs%coef_htfrtc%mixed_ref_frac(p_opt_ly_i(i),2)+ &
        p_opt_ly_r(i)*coefs%coef_htfrtc%mixed_ref_frac(p_opt_ly_i(i)+1,2)) &
        *cv(i)*gas_mass(gas_id_ozone)
      ENDFORALL
    CASE(3)
      FORALL(i=1:surf_level)
        x_ly(i,igc)=((1.0-p_opt_ly_r(i))* &
        coefs%coef_htfrtc%mixed_ref_frac(p_opt_ly_i(i),3)+ &
        p_opt_ly_r(i)*coefs%coef_htfrtc%mixed_ref_frac(p_opt_ly_i(i)+1,3)) &
        *cv(i)*gas_mass(gas_id_co2)
      ENDFORALL
    CASE(4)
      FORALL(i=1:surf_level)
        x_ly(i,igc)=((1.0-p_opt_ly_r(i))* &
        coefs%coef_htfrtc%mixed_ref_frac(p_opt_ly_i(i),4)+ &
        p_opt_ly_r(i)*coefs%coef_htfrtc%mixed_ref_frac(p_opt_ly_i(i)+1,4)) &
        *cv(i)*gas_mass(gas_id_n2o)
      ENDFORALL
    CASE(5)
      FORALL(i=1:surf_level)
        x_ly(i,igc)=((1.0-p_opt_ly_r(i))* &
        coefs%coef_htfrtc%mixed_ref_frac(p_opt_ly_i(i),5)+ &
        p_opt_ly_r(i)*coefs%coef_htfrtc%mixed_ref_frac(p_opt_ly_i(i)+1,5)) &
        *cv(i)*gas_mass(gas_id_co)
      ENDFORALL
    CASE(6)
      FORALL(i=1:surf_level)
        x_ly(i,igc)=((1.0-p_opt_ly_r(i))* &
        coefs%coef_htfrtc%mixed_ref_frac(p_opt_ly_i(i),6)+ &
        p_opt_ly_r(i)*coefs%coef_htfrtc%mixed_ref_frac(p_opt_ly_i(i)+1,6)) &
        *cv(i)*gas_mass(gas_id_ch4)
      ENDFORALL
    CASE(7)
      FORALL(i=1:surf_level)
        x_ly(i,igc)=((1.0-p_opt_ly_r(i))* &
        coefs%coef_htfrtc%mixed_ref_frac(p_opt_ly_i(i),7)+ &
        p_opt_ly_r(i)*coefs%coef_htfrtc%mixed_ref_frac(p_opt_ly_i(i)+1,7)) &
        *cv(i)*gas_mass(gas_id_so2)
      ENDFORALL
    END SELECT
    ENDDO

    !Gaseous optical properties, precalc
    FORALL(i=1:nlayers)
      t2_ly(i)=t_ly(i)*t_ly(i)
      q_ly(i)=x_ly(i,1)
      q2_ly(i)=q_ly(i)*q_ly(i)
    ENDFORALL

    do_lambertian=opts%rt_all%do_lambertian

    user_emis=.FALSE.
    IF (PRESENT(calcemis)) THEN
      user_emis=(.NOT. ANY(calcemis((i_p-1)*n_f+1:i_p*n_f)))
    ENDIF

    IF (user_emis) THEN

      DO i_f=1,n_f
        surf_em(i_f)=emissivity((i_p-1)*n_f+i_f)%emis_in
        emissivity((i_p-1)*n_f+i_f)%emis_out=emissivity((i_p-1)*n_f+i_f)%emis_in
      END DO

      IF(do_k) THEN

          dsurf_emdu=0.0_jprb
          dsurf_emdv=0.0_jprb
          dsurf_emdskint=0.0_jprb

      ENDIF

    ELSE

      !Non-directional horizontal wind
      IF(profiles(i_p)%skin%surftype==surftype_sea) THEN
        surf_wind=sqrt(profiles(i_p)%s2m%u**2+profiles(i_p)%s2m%v**2)
        surf_wind2=surf_wind**2
        surf_windm1=1.0_jprb/(surf_wind+1.0e-6_jprb)
      ENDIF

      IF(profiles(i_p)%skin%surftype==surftype_sea) THEN

        do_lambertian=.FALSE.

        IF (opts%rt_ir%ir_sea_emis_model==1) THEN
          DO i_f=1,n_f
            aems(i_f)=coefs%coef_htfrtc%coef_ssemp(1,i_f) &
            +surf_wind*coefs%coef_htfrtc%coef_ssemp(2,i_f) &
            +surf_wind2*coefs%coef_htfrtc%coef_ssemp(3,i_f)
            bems(i_f)=coefs%coef_htfrtc%coef_ssemp(4,i_f) &
            +surf_wind*coefs%coef_htfrtc%coef_ssemp(5,i_f) &
            +surf_wind2*coefs%coef_htfrtc%coef_ssemp(6,i_f)
            cems(i_f)=coefs%coef_htfrtc%coef_ssemp(7,i_f) &
            +surf_wind*coefs%coef_htfrtc%coef_ssemp(8,i_f)
            expf(i_f)= &
            exp( ( (coefs%coef_htfrtc%coef_ssemp(9,i_f)-60.0_jprb)**2 &
            - (profiles(i_p)%zenangle &
            -coefs%coef_htfrtc%coef_ssemp(9,i_f))**2 ) /cems(i_f) )
            surf_em(i_f)=aems(i_f)+(bems(i_f)-aems(i_f))*expf(i_f)
          ENDDO
        ELSE IF (opts%rt_ir%ir_sea_emis_model==2) THEN
          DO i_f=1,n_f
            aems(i_f)=coefs%coef_htfrtc%coef_iremis(1,i_f) &
            +surf_wind*coefs%coef_htfrtc%coef_iremis(2,i_f) &
            +surf_wind2*coefs%coef_htfrtc%coef_iremis(3,i_f)
            bems(i_f)=coefs%coef_htfrtc%coef_iremis(4,i_f) &
            +surf_wind*coefs%coef_htfrtc%coef_iremis(5,i_f) &
            +surf_wind2*coefs%coef_htfrtc%coef_iremis(6,i_f)
            cems(i_f)=coefs%coef_htfrtc%coef_iremis(7,i_f) &
            +surf_wind*coefs%coef_htfrtc%coef_iremis(8,i_f)
            expf(i_f)= &
            exp( ( (coefs%coef_htfrtc%coef_iremis(9,i_f)-60.0_jprb)**2 &
            - (profiles(i_p)%zenangle &
            -coefs%coef_htfrtc%coef_iremis(9,i_f))**2 ) /cems(i_f) )
            surf_em(i_f)=aems(i_f)+(bems(i_f)-aems(i_f))*expf(i_f)
            IF (ABS(coefs%coef_htfrtc%coef_iremis(10,i_f))>0.0 .AND. &
              ABS(coefs%coef_htfrtc%coef_iremis(11,i_f))>0.0) THEN
              surf_em(i_f)=surf_em(i_f)+ &
              (profiles(i_p)%skin%t- 301.2)* &
              (coefs%coef_htfrtc%coef_iremis(10,i_f)+ &
              exp(coefs%coef_htfrtc%coef_iremis(11,i_f)* &
              profiles(i_p)%zenangle**2 / 60.0_jprb**2))
            ENDIF
          ENDDO
        ENDIF

      ELSE IF (profiles(i_p)%skin%surftype==surftype_land) THEN
        surf_em=0.98_jprb
      ELSE IF (profiles(i_p)%skin%surftype==surftype_seaice) THEN
        surf_em=0.99_jprb
      ELSE
      ENDIF

      IF(PRESENT(emissivity)) THEN
        DO i_f=1,n_f
          emissivity((i_p-1)*n_f+i_f)%emis_out=surf_em(i_f)
        END DO
      END IF

      IF(do_k) THEN

        IF(profiles(i_p)%skin%surftype==surftype_sea) THEN

          IF (opts%rt_ir%ir_sea_emis_model==1) THEN
            DO i_f=1,n_f
              daems(i_f)=surf_windm1*(coefs%coef_htfrtc%coef_ssemp(2,i_f) &
              +2.0_jprb*surf_wind*coefs%coef_htfrtc%coef_ssemp(3,i_f))
              dbems(i_f)=surf_windm1*(coefs%coef_htfrtc%coef_ssemp(5,i_f) &
              +2.0_jprb*surf_wind*coefs%coef_htfrtc%coef_ssemp(6,i_f))
              dcems(i_f)=surf_windm1*coefs%coef_htfrtc%coef_ssemp(8,i_f)
              dexpf(i_f)= &
              -( (coefs%coef_htfrtc%coef_ssemp(9,i_f)-60.0_jprb)**2 &
              - (profiles(i_p)%zenangle &
              -coefs%coef_htfrtc%coef_ssemp(9,i_f))**2 ) &
              / (cems(i_f)*cems(i_f))
              dsurf_em(i_f)=(daems(i_f)+(dbems(i_f)-daems(i_f)) &
              +dcems(i_f)*(bems(i_f)-aems(i_f))*dexpf(i_f))*expf(i_f)
              dsurf_emdu(i_f)=dsurf_em(i_f)*profiles(i_p)%s2m%u
              dsurf_emdv(i_f)=dsurf_em(i_f)*profiles(i_p)%s2m%v
              dsurf_emdskint(i_f)=0.0
            ENDDO
          ELSE IF (opts%rt_ir%ir_sea_emis_model==2) THEN
            DO i_f=1,n_f
              daems(i_f)=surf_windm1*(coefs%coef_htfrtc%coef_iremis(2,i_f) &
              +2.0_jprb*surf_wind*coefs%coef_htfrtc%coef_iremis(3,i_f))
              dbems(i_f)=surf_windm1*(coefs%coef_htfrtc%coef_iremis(5,i_f) &
              +2.0_jprb*surf_wind*coefs%coef_htfrtc%coef_iremis(6,i_f))
              dcems(i_f)=surf_windm1*coefs%coef_htfrtc%coef_iremis(8,i_f)
              dexpf(i_f)= &
              -( (coefs%coef_htfrtc%coef_iremis(9,i_f)-60.0_jprb)**2 &
              - (profiles(i_p)%zenangle &
              -coefs%coef_htfrtc%coef_iremis(9,i_f))**2 ) &
              / (cems(i_f)*cems(i_f))
              dsurf_em(i_f)=(daems(i_f)+(dbems(i_f)-daems(i_f)) &
              +dcems(i_f)*(bems(i_f)-aems(i_f))*dexpf(i_f))*expf(i_f)
              dsurf_emdu(i_f)=dsurf_em(i_f)*profiles(i_p)%s2m%u
              dsurf_emdv(i_f)=dsurf_em(i_f)*profiles(i_p)%s2m%v
              IF (ABS(coefs%coef_htfrtc%coef_iremis(10,i_f))>0.0 .AND. &
                ABS(coefs%coef_htfrtc%coef_iremis(11,i_f))>0.0) THEN
                dsurf_emdskint(i_f)=profiles(i_p)%skin%t* &
                (coefs%coef_htfrtc%coef_iremis(10,i_f)+ &
                exp(coefs%coef_htfrtc%coef_iremis(11,i_f)* &
                profiles(i_p)%zenangle**2 / 60.0_jprb**2))
              ELSE
                dsurf_emdskint(i_f)=0.0
              ENDIF
           ENDDO
         ENDIF

        ELSE

          dsurf_emdu=0.0_jprb
          dsurf_emdv=0.0_jprb
          dsurf_emdskint=0.0_jprb

        ENDIF

      END IF ! do_k

    END IF ! useremis

    IF (do_lambertian) THEN
      IF (PRESENT(emissivity)) THEN
        DO i_f=1,n_f
          sp(i_f)=emissivity((i_p-1)*n_f+i_f)%specularity
        ENDDO
      ELSE
        sp=0.0
      ENDIF
    ELSE
      sp=1.0
    ENDIF

    !Find the level just above the cloud
    IF(opts%htfrtc_opts%simple_cloud) THEN
      cloud_level=nlevels-1
      cloud_layer_frac=1.0
      DO jb=1,nlevels-1
        IF ((profiles(i_p)%ctp.GT.profiles(i_p)%p(jb)).AND. &
          (profiles(i_p)%ctp.LE.profiles(i_p)%p(jb+1))) THEN
          cloud_level=jb
          cloud_layer_frac=(profiles(i_p)%ctp-profiles(i_p)%p(jb)) &
          /(profiles(i_p)%p(jb+1)-profiles(i_p)%p(jb))
        ENDIF
      ENDDO
      cf=profiles(i_p)%cfraction
    ELSE
      cf=0.0
    ENDIF

    !Main long loop over centroid frequencies
    DO i_f=1,n_f

      !Gaseous optical properties calc
      dx=0.0
      DO igc=1,ngc
        j=gid(igc)
        IF (laddf(j,i_f)) THEN
          DO i=1,nlayers
            dx(i,igc)= &
            ((1.0-p_opt_ly_r(i)) &
            *(coefs%coef_htfrtc%coef_l(1,p_opt_ly_i(i),j,i_f) &
            +t_ly(i)*coefs%coef_htfrtc%coef_l(2,p_opt_ly_i(i),j,i_f) &
            +t2_ly(i)*coefs%coef_htfrtc%coef_l(3,p_opt_ly_i(i),j,i_f)) &
            +p_opt_ly_r(i)*(coefs%coef_htfrtc%coef_l(1,p_opt_ly_i(i)+1,j,i_f) &
            +t_ly(i)*coefs%coef_htfrtc%coef_l(2,p_opt_ly_i(i)+1,j,i_f) &
            +t2_ly(i)*coefs%coef_htfrtc%coef_l(3,p_opt_ly_i(i)+1,j,i_f)))
          ENDDO
        ENDIF
      ENDDO

      FORALL(i=1:nlayers)
        tau(i)=max(0._jprb,mpath(i)* &
        (sum(x_ly(i,1:ngc)*dx(i,1:ngc))+p_ly(i) &
        *((coefs%coef_htfrtc%coef_ct(1,i_f) &
        *(coefs%coef_htfrtc%coef_ctt(1,t_b_i(i),i_f)+t_b_r(i)* &
        (coefs%coef_htfrtc%coef_ctt(1,t_b_i(i)+1,i_f) &
        -coefs%coef_htfrtc%coef_ctt(1,t_b_i(i),i_f))) &
        *q2_ly(i))+ &
        (coefs%coef_htfrtc%coef_ct(2,i_f) &
        *(coefs%coef_htfrtc%coef_ctt(2,t_b_i(i),i_f)+t_b_r(i)* &
        (coefs%coef_htfrtc%coef_ctt(2,t_b_i(i)+1,i_f) &
        -coefs%coef_htfrtc%coef_ctt(2,t_b_i(i),i_f))) &
        *(q_ly(i)*e_wvp-q2_ly(i))))))
      ENDFORALL
      tausl=tau*mft
      tr=exp(-tausl)

      IF(do_lambertian) THEN
        tausllb=tau*mftlb
        trlb=exp(-tausllb)
      ENDIF

      !Planckians calc
      FORALL(i=1:nlevels)
        b(i)=coefs%coef_htfrtc%coef_b(t_b_i(i),i_f) &
        +t_b_r(i)*(coefs%coef_htfrtc%coef_b(t_b_i(i)+1,i_f) &
        -coefs%coef_htfrtc%coef_b(t_b_i(i),i_f))
      ENDFORALL
      b_surf=coefs%coef_htfrtc%coef_b(t_b_surf_i,i_f) &
      +t_b_surf_r*(coefs%coef_htfrtc%coef_b(t_b_surf_i+1,i_f) &
      -coefs%coef_htfrtc%coef_b(t_b_surf_i,i_f))

      FORALL(i=1:nlayers)
        lt(i)=coefs%coef_htfrtc%coef_lt(min(n_lt, &
        max(1,1+int(lt_inc_m1*(tausl(i)-lt_min)))))
        blydown(i)=b(i)+lt(i)*(b(i+1)-b(i))
        blyup(i)=b(i+1)+lt(i)*(b(i)-b(i+1))
      ENDFORALL

      IF(do_lambertian) THEN
        FORALL(i=1:nlayers)
          ltlb(i)=coefs%coef_htfrtc%coef_lt(min(n_lt, &
          max(1,1+int(lt_inc_m1*(tausllb(i)-lt_min)))))
          blydownlb(i)=b(i)+ltlb(i)*(b(i+1)-b(i))
        ENDFORALL
      ENDIF

      !Radiative transfer, downwelling, then upwelling
      rdown(1)=0.0
      DO i=2,nlevels
        rdown(i)=blydown(i-1)+(rdown(i-1)-blydown(i-1))*tr(i-1)
      ENDDO
      rdown_surf=blydown(surf_level) &
      +(rdown(surf_level)-blydown(surf_level))*tr(surf_level)

      IF(do_lambertian) THEN
        rdownlb(1)=0.0
        DO i=2,nlevels
          rdownlb(i)=blydownlb(i-1)+(rdownlb(i-1)-blydown(i-1))*trlb(i-1)
        ENDDO
      rdown_surflb=blydownlb(surf_level) &
      +(rdownlb(surf_level)-blydownlb(surf_level))*trlb(surf_level)
      ENDIF

      DO i=nlevels,surf_level+2,-1
        rup(i)=0.0
      ENDDO
      IF(do_lambertian) THEN
        rup(surf_level+1)=b_surf*surf_em(i_f)+ &
        (rdown_surflb+sp(i_f)*(rdown_surf-rdown_surflb))*(1.0-surf_em(i_f))
      ELSE
        rup(surf_level+1)=b_surf*surf_em(i_f)+rdown_surf*(1.0-surf_em(i_f))
      ENDIF
      rup(surf_level)=blyup(surf_level)+(rup(surf_level+1) &
      -blyup(surf_level))*tr(surf_level)
      DO i=surf_level,2,-1
        rup(i-1)=blyup(i-1)+(rup(i)-blyup(i-1))*tr(i-1)
      ENDDO
      rc(i_f)=(rup(1)-coefs%coef_htfrtc%val_mean(i_f))*normm1(i_f)

      IF(opts%htfrtc_opts%overcast) THEN
        DO i=2,nlevels
          totup(i)=product(tr(1:i-1))
          rc_oc(i-1,i_f)= &
          (rup(1)+(b(i)-rup(i))*totup(i)-coefs%coef_htfrtc%val_mean(i_f)) &
          *normm1(i_f)
        ENDDO
      ENDIF

      rup_cld=0.0
      IF(opts%htfrtc_opts%simple_cloud) THEN
        rup_cld(cloud_level)=b(cloud_level)+ &
        0.5*(b(cloud_level+1)-b(cloud_level))* &
        cloud_layer_frac*(1.0+exp(-tausl(cloud_level)*cloud_layer_frac))
        DO i=cloud_level,2,-1
          rup_cld(i-1)=blyup(i-1)+(rup_cld(i)-blyup(i-1))*tr(i-1)
        ENDDO
        rc_cld(i_f)=(rup_cld(1)-coefs%coef_htfrtc%val_mean(i_f))*normm1(i_f)
      ENDIF

      r_f(i_f)=(1.0-cf)*rup(1)+cf*rup_cld(1)

      IF(do_k) THEN

        DO i=1,nlayers
          dttmp(i)= &
          sum(x_ly(i,1:ngc)* &
          ((1.0-p_opt_ly_r(i))* &
          (coefs%coef_htfrtc%coef_l(2,p_opt_ly_i(i),gid(1:ngc),i_f) &
          +2.0*t_ly(i) &
          *coefs%coef_htfrtc%coef_l(3,p_opt_ly_i(i),gid(1:ngc),i_f)) &
          +p_opt_ly_r(i)* &
          (coefs%coef_htfrtc%coef_l(2,p_opt_ly_i(i)+1,gid(1:ngc),i_f) &
          +2.0*t_ly(i)* &
          coefs%coef_htfrtc%coef_l(3,p_opt_ly_i(i)+1,gid(1:ngc),i_f)))) &
          +p_ly(i) &
          *((coefs%coef_htfrtc%coef_ct(1,i_f) &
          *t_b_inc_m1*(coefs%coef_htfrtc%coef_ctt(1,t_b_i(i)+1,i_f)- &
          coefs%coef_htfrtc%coef_ctt(1,t_b_i(i),i_f)) &
          *q2_ly(i))+ &
          (coefs%coef_htfrtc%coef_ct(2,i_f) &
          *t_b_inc_m1*(coefs%coef_htfrtc%coef_ctt(2,t_b_i(i)+1,i_f)- &
           coefs%coef_htfrtc%coef_ctt(2,t_b_i(i),i_f)) &
          *(q_ly(i)*e_wvp-q2_ly(i))))
        ENDDO

        DO i=1,nlayers
          dxtmp(i)=dx(i,1) &
          +p_ly(i) &
          *((coefs%coef_htfrtc%coef_ct(1,i_f) &
          *(coefs%coef_htfrtc%coef_ctt(1,t_b_i(i),i_f)+t_b_r(i)* &
          (coefs%coef_htfrtc%coef_ctt(1,t_b_i(i)+1,i_f) &
          -coefs%coef_htfrtc%coef_ctt(1,t_b_i(i),i_f))) &
          *2.0*q_ly(i))+ &
          (coefs%coef_htfrtc%coef_ct(2,i_f) &
          *(coefs%coef_htfrtc%coef_ctt(2,t_b_i(i),i_f)+t_b_r(i)* &
          (coefs%coef_htfrtc%coef_ctt(2,t_b_i(i)+1,i_f) &
          -coefs%coef_htfrtc%coef_ctt(2,t_b_i(i),i_f))) &
          *(e_wvp-2.0*q_ly(i))))
        ENDDO

        DO i=1,nlayers
          dtrdt(i)=-tr(i)*mpath(i)*mft(i)*dttmp(i)
          dtrdx(1,i)=-tr(i)*mpath(i)*mft(i)*dxtmp(i)
          DO igu=2,ngu
            dtrdx(igu,i)=-tr(i)*mpath(i)*mft(i)*dx(i,igu)
          ENDDO
        ENDDO

        IF(do_lambertian) THEN
          DO i=1,nlayers
            dtrdtlb(i)=-trlb(i)*mpath(i)*mftlb(i)*dttmp(i)
            dtrdxlb(1,i)=-tr(i)*mpath(i)*mft(i)*dxtmp(i)
            DO igu=2,ngu
              dtrdxlb(igu,i)=-trlb(i)*mpath(i)*mftlb(i)*dx(i,igu)
            ENDDO
          ENDDO
        ENDIF

        FORALL(i=1:nlevels)
          dbdt(i)=t_b_inc_m1*(coefs%coef_htfrtc%coef_b(t_b_i(i)+1,i_f) &
          -coefs%coef_htfrtc%coef_b(t_b_i(i),i_f))
        ENDFORALL
        db_surfdt=t_b_inc_m1*(coefs%coef_htfrtc%coef_b(t_b_surf_i+1,i_f) &
        -coefs%coef_htfrtc%coef_b(t_b_surf_i,i_f))

        FORALL(i=1:nlayers)
          dblydowndt(i)=dbdt(i)+lt(i)*(dbdt(i+1)-dbdt(i))
          dblyupdt(i)=dbdt(i+1)+lt(i)*(dbdt(i)-dbdt(i+1))
        ENDFORALL
        IF(do_lambertian) THEN
          FORALL(i=1:nlayers)
            dblydowndtlb(i)=dbdt(i)+ltlb(i)*(dbdt(i+1)-dbdt(i))
          ENDFORALL
        ENDIF

        tot=product(tr(1:surf_level))
        IF(do_lambertian) THEN
          DO i=1,surf_level
            totdown(i)=product(tr(i+1:surf_level))*tot*(1.0-surf_em(i_f))
            totdownlb(i)=product(trlb(i+1:surf_level))*tot*(1.0-surf_em(i_f))
            drdowndt(i)=(1.0-sp(i_f))*(dblydowndtlb(i)*(1.0-trlb(i))+ &
            (rdownlb(i)-blydownlb(i))*dtrdtlb(i))*totdownlb(i) &
            +sp(i_f)*(dblydowndt(i)*(1.0-tr(i))+ &
            (rdown(i)-blydown(i))*dtrdt(i))*totdown(i)
            DO igu=1,ngu
              drdowndx(igu,i)=(1.0-sp(i_f))* &
              (rdownlb(i)-blydownlb(i))*dtrdxlb(igu,i)*totdownlb(i) &
              +sp(i_f)*(rdown(i)-blydown(i))*dtrdx(igu,i)*totdown(i)
            ENDDO
          ENDDO
        ELSE
          DO i=1,surf_level
            totdown(i)=product(tr(i+1:surf_level))*tot*(1.0-surf_em(i_f))
            drdowndt(i)=(dblydowndt(i)*(1.0-tr(i))+ &
            (rdown(i)-blydown(i))*dtrdt(i))*totdown(i)
            DO igu=1,ngu
              drdowndx(igu,i)=(rdown(i)-blydown(i))*dtrdx(igu,i)*totdown(i)
            ENDDO
          ENDDO
        ENDIF

        IF(opts%htfrtc_opts%simple_cloud) THEN
          DO i=1,surf_level
            totup(i)=product(tr(1:i-1))
            DO igu=1,ngu
              drupdx(igu,i)=((1.0-cf)*rup(i+1)+ &
              cf*rup_cld(i+1)-blyup(i))*dtrdx(igu,i)*totup(i)
            ENDDO
          ENDDO
          DO i=1,cloud_level
            drupdt(i)=(dblyupdt(i)*(1.0-tr(i)) &
            +((1.0-cf)*rup(i+1)+cf*rup_cld(i+1)-blyup(i))*dtrdt(i))*totup(i)
          ENDDO
          DO i=cloud_level+1,surf_level
            drupdt(i)=((1.0-cf)*rup(i+1)-blyup(i))*dtrdt(i)*totup(i)
          ENDDO
        ELSE
          DO i=1,surf_level
            totup(i)=product(tr(1:i-1))
            drupdt(i)=(dblyupdt(i)*(1.0-tr(i)) &
            +(rup(i+1)-blyup(i))*dtrdt(i))*totup(i)
            DO igu=1,ngu
              drupdx(igu,i) =(rup(i+1)-blyup(i))*dtrdx(igu,i)*totup(i)
            ENDDO
          ENDDO
        ENDIF

        !Jacobians at centroids
        drcdt(i_f,1)=0.5*(drdowndt(1)+drupdt(1))*normm1(i_f)
        DO igu=1,ngu
          drcdx(i_f,igu,1)=0.5*(drdowndx(igu,1)+drupdx(igu,1))*normm1(i_f)
        ENDDO
        DO i=2,surf_level
        drcdt(i_f,i)=0.5*(drdowndt(i-1)+drdowndt(i)+drupdt(i-1)+drupdt(i)) &
        *normm1(i_f)
          DO igu=1,ngu
            drcdx(i_f,igu,i)=0.5*(drdowndx(igu,i-1)+drdowndx(igu,i) &
            +drupdx(igu,i-1)+drupdx(igu,i))*normm1(i_f)
          ENDDO
        ENDDO
        drcdt(i_f,surf_level+1)=0.5*(drdowndt(surf_level) &
        +drupdt(surf_level))*normm1(i_f)
        DO igu=1,ngu
          drcdx(i_f,igu,surf_level+1)=0.5*(drdowndx(igu,surf_level) &
          +drupdx(igu,surf_level))*normm1(i_f)
        ENDDO
        DO i=surf_level+2, nlevels
          drcdt(i_f,i)=0.0
          DO igu=1,ngu
            drcdx(i_f,igu,i)=0.0
          ENDDO
        ENDDO

        drcds(i_f,1)=(1.0-cf)* &
        (db_surfdt*surf_em(i_f)+(b_surf-rdown_surf)*dsurf_emdskint(i_f)) &
        *tot*normm1(i_f)
        drcds(i_f,3)=drcdt(i_f,surf_level+1)
        drcds(i_f,4)=drcdx(i_f,1,surf_level+1)
        IF(do_lambertian) THEN
          drcds(i_f,2)=(1.0-cf)*(rdown_surf-rdown_surflb)*(1.0-surf_em(i_f))
          drcds(i_f,5)=(1.0-cf)*totup(surf_level)*normm1(i_f)/dp(surf_level)* &
          ((1.0-sp(i_f))*(rdownlb(surf_level)-blydownlb(surf_level))* &
          (1.0-surf_em(i_f))*trlb(surf_level)*(-tausllb(surf_level)) &
          +((sp(i_f)*(rdown(surf_level)-blydown(surf_level))*(1.0-surf_em(i_f)))+ &
          (rup(surf_level+1)-blyup(surf_level))) &
          *tr(surf_level)*(-tausl(surf_level)))
          drcds(i_f,6)=0.0
          drcds(i_f,7)=0.0
          drcdsem(i_f)=(1.0-cf)* &
          (b_surf-(rdown_surflb+sp(i_f)*(rdown_surf-rdown_surflb)))*tot*normm1(i_f)
        ELSE
          drcds(i_f,2)=0.0
          drcds(i_f,5)=(1.0-cf)*totup(surf_level)*normm1(i_f)/dp(surf_level)* &
          ((rdown(surf_level)-blydown(surf_level))*(1.0-surf_em(i_f))+ &
          (rup(surf_level+1)-blyup(surf_level)))*tr(surf_level)*(-tausl(surf_level))
          drcds(i_f,6)=(1.0-cf)*(b_surf-rdown_surf)*dsurf_emdu(i_f)*tot*normm1(i_f)
          drcds(i_f,7)=(1.0-cf)*(b_surf-rdown_surf)*dsurf_emdv(i_f)*tot*normm1(i_f)
          drcdsem(i_f)=(1.0-cf)*(b_surf-rdown_surf)*tot*normm1(i_f)
        ENDIF
        IF(opts%htfrtc_opts%simple_cloud) THEN
          drcds(i_f,8)=rc_cld(i_f)-rc(i_f)
          drcds(i_f,9)=(1.0-cf)*totup(cloud_level)*normm1(i_f)/dp(cloud_level)* &
          0.5*(b(cloud_level+1)-b(cloud_level)) &
          *(1.0+(1.0-tausl(cloud_level)*cloud_layer_frac) &
          *exp(-tausl(cloud_level)*cloud_layer_frac))
        ELSE
          drcds(i_f,8)=0.0
          drcds(i_f,9)=0.0
        ENDIF

      ENDIF ! do_k

    ENDDO ! n_f

    IF (PRESENT(emissivity_k)) THEN

      DO i_f=1,n_f
        emissivity_k((i_p-1)*n_f+i_f)%emis_in=drcdsem(i_f)
        emissivity_k((i_p-1)*n_f+i_f)%specularity=drcds(i_f,2)
      END DO

      IF(opts%rt_all%switchrad) THEN

        WHERE (r_f(1:n_f)<=0.0_jprb)
          swrad_f(1:n_f)=0.0_jprb
        ELSEWHERE
          bt_f(1:n_f)=planck_c2_wn_f &
          /log(1.0+planck_c1_wn3_f(1:n_f)/r_f(1:n_f))
          swrad_f(1:n_f)=planck_c1_wn3_f(1:n_f) &
          *(bt_f(1:n_f)**2) &
          /(planck_c2_wn_f(1:n_f)* &
          r_f(1:n_f) &
          *(r_f(1:n_f) &
          +planck_c1_wn3_f(1:n_f)))
        ENDWHERE

        DO i_f=1,n_f
          emissivity_k((i_p-1)*n_f+i_f)%emis_in= &
          swrad_f(i_f)*emissivity_k((i_p-1)*n_f+i_f)%emis_in
          emissivity_k((i_p-1)*n_f+i_f)%specularity= &
          swrad_f(i_f)*emissivity_k((i_p-1)*n_f+i_f)%specularity
        ENDDO

      ENDIF

    END IF

    !PC scores
    FORALL(i_pc=1:n_pc)
      pccomp%clear_pcscores((i_p-1)*n_pc+i_pc)= &
      sum(rc(1:n_f)*coefs%coef_htfrtc%coef_pdt(1:n_f,i_pc))
    END FORALL ! n_pc

    IF(opts%htfrtc_opts%simple_cloud) THEN
      FORALL(i_pc=1:n_pc)
        pccomp%cloudy_pcscores((i_p-1)*n_pc+i_pc)= &
        sum(rc_cld(1:n_f)*coefs%coef_htfrtc%coef_pdt(1:n_f,i_pc))
        pccomp%total_pcscores((i_p-1)*n_pc+i_pc)= &
        (1.0-cf)*pccomp%clear_pcscores((i_p-1)*n_pc+i_pc) &
        +cf*pccomp%cloudy_pcscores((i_p-1)*n_pc+i_pc)
      END FORALL ! n_pc
    ELSE
      FORALL(i_pc=1:n_pc)
        pccomp%total_pcscores((i_p-1)*n_pc+i_pc)= &
        pccomp%clear_pcscores((i_p-1)*n_pc+i_pc)
      ENDFORALL
    ENDIF

    IF(opts%htfrtc_opts%overcast) THEN
    rc_oc_t=transpose(rc_oc)
      DO i=1,nlayers
        DO i_pc=1,n_pc_oc
          overcast_pcscores_t(i_pc,i)= &
          sum(rc_oc_t(1:n_f,i)*coefs%coef_htfrtc%coef_pdt(1:n_f,i_pc))
        ENDDO
      ENDDO
    overcast_pcscores=transpose(overcast_pcscores_t)
      DO i_pc=1,n_pc_oc
        DO i=1,nlayers
          pccomp%overcast_pcscores(i,(i_p-1)*n_pc+i_pc)= &
          overcast_pcscores(i,i_pc)
        ENDDO
      ENDDO
    ENDIF

    IF(do_k) THEN

      DO i_pc=1,n_pc
        FORALL(i=1:nlevels)
          profiles_k_pc_t(i,i_pc)= &
          sum(drcdt(1:n_f,i)*coefs%coef_htfrtc%coef_pdt(1:n_f,i_pc))
        ENDFORALL
        DO i=1,nlevels
          DO igu=1,ngu
            profiles_k_pc_x(i,igu,i_pc)= &
            sum(drcdx(1:n_f,igu,i)*coefs%coef_htfrtc%coef_pdt(1:n_f,i_pc))
          ENDDO
        ENDDO
        DO is=1,ns
          profiles_k_pc_s(is,i_pc)= &
          sum(drcds(1:n_f,is)*coefs%coef_htfrtc%coef_pdt(1:n_f,i_pc))
        ENDDO
      ENDDO ! n_pc

    ENDIF ! do_k

    IF (do_k) THEN
      SELECT CASE(profiles(1)%gas_units)
      CASE(gas_unit_ppmvdry,gas_unit_ppmv)
        DO igu=1,ngu
          DO i=1,nlevels
            profiles_k_pc_x(i,igu,1:n_pc)= &
            profiles_k_pc_x(i,igu,1:n_pc)*(cv(i)*gas_mass(gmid(igu)))
          ENDDO
        ENDDO
        profiles_k_pc_s(4,1:n_pc)= &
        profiles_k_pc_s(4,1:n_pc)*(cv(surf_level)*gas_mass(gas_id_watervapour))
      END SELECT
    ENDIF

    IF(PRESENT(profiles_k_pc)) THEN
      DO i_pc=1,n_pc
        DO i=1,nlevels
          profiles_k_pc((i_p-1)*n_pc+i_pc)%t(i)=profiles_k_pc_t(i,i_pc)
        ENDDO
        DO igu=1,ngu
          SELECT CASE (gid(igu))
          CASE(1)
            DO i=1,nlevels
              profiles_k_pc((i_p-1)*n_pc+i_pc)%q(i)= &
              profiles_k_pc_x(i,igu,i_pc)
            ENDDO
          CASE(2)
            DO i=1,nlevels
              profiles_k_pc((i_p-1)*n_pc+i_pc)%o3(i)= &
              profiles_k_pc_x(i,igu,i_pc)
            ENDDO
          CASE(3)
            DO i=1,nlevels
              profiles_k_pc((i_p-1)*n_pc+i_pc)%co2(i)= &
              profiles_k_pc_x(i,igu,i_pc)
            ENDDO
          CASE(4)
            DO i=1,nlevels
              profiles_k_pc((i_p-1)*n_pc+i_pc)%n2o(i)= &
              profiles_k_pc_x(i,igu,i_pc)
            ENDDO
          CASE(5)
            DO i=1,nlevels
              profiles_k_pc((i_p-1)*n_pc+i_pc)%co(i)= &
              profiles_k_pc_x(i,igu,i_pc)
            ENDDO
          CASE(6)
            DO i=1,nlevels
              profiles_k_pc((i_p-1)*n_pc+i_pc)%ch4(i)= &
              profiles_k_pc_x(i,igu,i_pc)
            ENDDO
          CASE(7)
            DO i=1,nlevels
              profiles_k_pc((i_p-1)*n_pc+i_pc)%so2(i)= &
              profiles_k_pc_x(i,igu,i_pc)
            ENDDO
          END SELECT
        ENDDO
        profiles_k_pc((i_p-1)*n_pc+i_pc)%skin%t= &
        profiles_k_pc_s(1,i_pc)
        profiles_k_pc((i_p-1)*n_pc+i_pc)%s2m%t= &
        profiles_k_pc_s(3,i_pc)
        profiles_k_pc((i_p-1)*n_pc+i_pc)%s2m%q= &
        profiles_k_pc_s(4,i_pc)
        profiles_k_pc((i_p-1)*n_pc+i_pc)%s2m%o=0.0_jprb
        profiles_k_pc((i_p-1)*n_pc+i_pc)%s2m%p= &
        profiles_k_pc_s(5,i_pc)
        profiles_k_pc((i_p-1)*n_pc+i_pc)%s2m%u= &
        profiles_k_pc_s(6,i_pc)
        profiles_k_pc((i_p-1)*n_pc+i_pc)%s2m%v= &
        profiles_k_pc_s(7,i_pc)
        profiles_k_pc((i_p-1)*n_pc+i_pc)%cfraction= &
        profiles_k_pc_s(8,i_pc)
        profiles_k_pc((i_p-1)*n_pc+i_pc)%ctp= &
        profiles_k_pc_s(9,i_pc)
      ENDDO ! i_pc
    ENDIF ! present profiles_k_pc

    !Reconstruct for radiances / brightness temperatures
    IF(opts%htfrtc_opts%reconstruct) THEN

      !Radiances
      FORALL(i_ch=1:n_ch)
        pccomp%clear_pccomp((i_p-1)*n_ch+i_ch)= &
        coefs%coef_htfrtc%ch_mean(i_ch)+ &
        sum(pccomp%clear_pcscores((i_p-1)*n_pc+1:i_p*n_pc) &
        *coefs%coef_htfrtc%pc(1:n_pc,i_ch))
      END FORALL
      IF(opts%htfrtc_opts%simple_cloud) THEN
        FORALL(i_ch=1:n_ch)
          pccomp%cloudy_pccomp((i_p-1)*n_ch+i_ch)= &
          coefs%coef_htfrtc%ch_mean(i_ch)+&
          sum(pccomp%cloudy_pcscores((i_p-1)*n_pc+1:i_p*n_pc) &
          *coefs%coef_htfrtc%pc(1:n_pc,i_ch))
          pccomp%total_pccomp((i_p-1)*n_ch+i_ch)= &
          coefs%coef_htfrtc%ch_mean(i_ch)+&
          sum(pccomp%total_pcscores((i_p-1)*n_pc+1:i_p*n_pc) &
          *coefs%coef_htfrtc%pc(1:n_pc,i_ch))
        END FORALL
      ELSE
        FORALL(i_ch=1:n_ch)
          pccomp%total_pccomp((i_p-1)*n_ch+i_ch)= &
          pccomp%clear_pccomp((i_p-1)*n_ch+i_ch)
        ENDFORALL
      ENDIF

      IF(opts%htfrtc_opts%overcast) THEN
        DO i_ch=1,n_ch
          DO i=1,nlayers
            pccomp%overcast_pccomp(i,(i_p-1)*n_ch+i_ch)= &
            coefs%coef_htfrtc%ch_mean(i_ch)+ &
            sum(overcast_pcscores_t(1:n_pc_oc,i)* &
            coefs%coef_htfrtc%pc(1:n_pc_oc,i_ch))
          ENDDO
        ENDDO
      ENDIF

      !Deal with unphysical values and calculate brightness temperatures
      WHERE (pccomp%clear_pccomp((i_p-1)*n_ch+1:i_p*n_ch)<0.0_jprb)
        pccomp%clear_pccomp((i_p-1)*n_ch+1:i_p*n_ch)=0.0_jprb
        pccomp%bt_clear_pccomp((i_p-1)*n_ch+1:i_p*n_ch)=0.0_jprb
      ELSEWHERE
        pccomp%bt_clear_pccomp((i_p-1)*n_ch+1:i_p*n_ch)= &
        planck_c2_wn(1:n_ch)/log(1.0+planck_c1_wn3(1:n_ch) &
        /pccomp%clear_pccomp((i_p-1)*n_ch+1:i_p*n_ch))
      ENDWHERE
      IF(opts%htfrtc_opts%simple_cloud) THEN
        WHERE (pccomp%cloudy_pccomp((i_p-1)*n_ch+1:i_p*n_ch)<0.0_jprb)
          pccomp%cloudy_pccomp((i_p-1)*n_ch+1:i_p*n_ch)=0.0_jprb
        ENDWHERE
        WHERE (pccomp%total_pccomp((i_p-1)*n_ch+1:i_p*n_ch)<0.0_jprb)
          pccomp%total_pccomp((i_p-1)*n_ch+1:i_p*n_ch)=0.0_jprb
          pccomp%bt_pccomp((i_p-1)*n_ch+1:i_p*n_ch)=0.0_jprb
        ELSEWHERE
          pccomp%bt_pccomp((i_p-1)*n_ch+1:i_p*n_ch)= &
          planck_c2_wn(1:n_ch)/log(1.0+planck_c1_wn3(1:n_ch) &
          /pccomp%total_pccomp((i_p-1)*n_ch+1:i_p*n_ch))
        ENDWHERE
      ELSE
        FORALL(i_ch=1:n_ch)
          pccomp%total_pccomp((i_p-1)*n_ch+i_ch)= &
          pccomp%clear_pccomp((i_p-1)*n_ch+i_ch)
          pccomp%bt_pccomp((i_p-1)*n_ch+i_ch)= &
          pccomp%bt_clear_pccomp((i_p-1)*n_ch+i_ch)
        ENDFORALL
      ENDIF

      IF(do_k) THEN

        profiles_k_rec_x=0._jprb
        DO i_ch = 1, n_ch
          FORALL(i=1:nlevels)
            profiles_k_rec_t(i,i_ch)= &
            sum(profiles_k_pc_t(i,1:n_pc)*coefs%coef_htfrtc%pc(1:n_pc,i_ch))
          ENDFORALL
          DO igu=1,ngu
            IF (laddch(gid(igu),i_ch)) THEN
              DO i=1,nlevels
                profiles_k_rec_x(i,igu,i_ch)= &
                sum(profiles_k_pc_x(i,igu,1:n_pc) &
                *coefs%coef_htfrtc%pc(1:n_pc,i_ch))
                !profiles_k_rec_x(i,igu,i_ch)= &
                !sum(profiles_k_pc_x_t(1:n_pc,i,igu)
                !*coefs%coef_htfrtc%pc(1:n_pc,i_ch))
              ENDDO
            ENDIF
          ENDDO
          DO is=1,ns
            profiles_k_rec_s(i_ch,is)= &
            sum(profiles_k_pc_s(is,1:n_pc)*coefs%coef_htfrtc%pc(1:n_pc,i_ch))
          ENDDO
        ENDDO

        IF(opts%rt_all%switchrad) THEN

          WHERE (pccomp%total_pccomp((i_p-1)*n_ch+1:i_p*n_ch)<=0.0_jprb)
            swrad_tmp(1:n_ch)=0.0_jprb
          ELSEWHERE
            swrad_tmp(1:n_ch)=planck_c1_wn3(1:n_ch) &
            *(pccomp%bt_pccomp((i_p-1)*n_ch+1:i_p*n_ch))**2 &
            /(planck_c2_wn(1:n_ch)* &
            pccomp%total_pccomp((i_p-1)*n_ch+1:i_p*n_ch) &
            *(pccomp%total_pccomp((i_p-1)*n_ch+1:i_p*n_ch) &
            +planck_c1_wn3(1:n_ch)))
          ENDWHERE

          DO i_ch=1,n_ch

            DO i=1,nlevels
              profiles_k_rec((i_p-1)*n_ch+i_ch)%t(i)= &
              swrad_tmp(i_ch)*profiles_k_rec_t(i,i_ch)
            ENDDO
            DO igu=1,ngu
              SELECT CASE (gid(igu))
              CASE(1)
                DO i=1,nlevels
                  profiles_k_rec((i_p-1)*n_ch+i_ch)%q(i)= &
                  swrad_tmp(i_ch)*profiles_k_rec_x(i,igu,i_ch)
                ENDDO
              CASE(2)
                DO i=1,nlevels
                  profiles_k_rec((i_p-1)*n_ch+i_ch)%o3(i)= &
                  swrad_tmp(i_ch)*profiles_k_rec_x(i,igu,i_ch)
                ENDDO
              CASE(3)
                DO i=1,nlevels
                  profiles_k_rec((i_p-1)*n_ch+i_ch)%co2(i)= &
                  swrad_tmp(i_ch)*profiles_k_rec_x(i,igu,i_ch)
                ENDDO
              CASE(4)
                DO i=1,nlevels
                  profiles_k_rec((i_p-1)*n_ch+i_ch)%n2o(i)= &
                  swrad_tmp(i_ch)*profiles_k_rec_x(i,igu,i_ch)
                ENDDO
              CASE(5)
                DO i=1,nlevels
                  profiles_k_rec((i_p-1)*n_ch+i_ch)%co(i)= &
                  swrad_tmp(i_ch)*profiles_k_rec_x(i,igu,i_ch)
                ENDDO
              CASE(6)
                 DO i=1,nlevels
                   profiles_k_rec((i_p-1)*n_ch+i_ch)%ch4(i)= &
                   swrad_tmp(i_ch)*profiles_k_rec_x(i,igu,i_ch)
                 ENDDO
              CASE(7)
                 DO i=1,nlevels
                   profiles_k_rec((i_p-1)*n_ch+i_ch)%so2(i)= &
                   swrad_tmp(i_ch)*profiles_k_rec_x(i,igu,i_ch)
                 ENDDO
              END SELECT
            ENDDO

            profiles_k_rec((i_p-1)*n_ch+i_ch)%skin%t= &
            swrad_tmp(i_ch)*profiles_k_rec_s(i_ch,1)
            profiles_k_rec((i_p-1)*n_ch+i_ch)%s2m%t= &
            swrad_tmp(i_ch)*profiles_k_rec_s(i_ch,3)
            profiles_k_rec((i_p-1)*n_ch+i_ch)%s2m%q= &
            swrad_tmp(i_ch)*profiles_k_rec_s(i_ch,4)
            profiles_k_rec((i_p-1)*n_ch+i_ch)%s2m%o=0.0_jprb
            profiles_k_rec((i_p-1)*n_ch+i_ch)%s2m%p= &
            swrad_tmp(i_ch)*profiles_k_rec_s(i_ch,5)
            profiles_k_rec((i_p-1)*n_ch+i_ch)%s2m%u= &
            swrad_tmp(i_ch)*profiles_k_rec_s(i_ch,6)
            profiles_k_rec((i_p-1)*n_ch+i_ch)%s2m%v= &
            swrad_tmp(i_ch)*profiles_k_rec_s(i_ch,7)
            profiles_k_rec((i_p-1)*n_ch+i_ch)%cfraction= &
            swrad_tmp(i_ch)*profiles_k_rec_s(i_ch,8)
            profiles_k_rec((i_p-1)*n_ch+i_ch)%ctp= &
            swrad_tmp(i_ch)*profiles_k_rec_s(i_ch,9)

          ENDDO ! i_ch

        ELSE

          DO i_ch=1,n_ch

            DO i=1,nlevels
               profiles_k_rec((i_p-1)*n_ch+i_ch)%t(i)=profiles_k_rec_t(i,i_ch)
            ENDDO
            DO igu=1,ngu
              SELECT CASE(gid(igu))
              CASE(1)
                DO i=1,nlevels
                  profiles_k_rec((i_p-1)*n_ch+i_ch)%q(i)= &
                  profiles_k_rec_x(i,igu,i_ch)
                ENDDO
              CASE(2)
                DO i=1,nlevels
                  profiles_k_rec((i_p-1)*n_ch+i_ch)%o3(i)= &
                  profiles_k_rec_x(i,igu,i_ch)
               ENDDO
              CASE(3)
                DO i=1,nlevels
                  profiles_k_rec((i_p-1)*n_ch+i_ch)%co2(i)= &
                  profiles_k_rec_x(i,igu,i_ch)
                ENDDO
              CASE(4)
                DO i=1,nlevels
                  profiles_k_rec((i_p-1)*n_ch+i_ch)%n2o(i)= &
                  profiles_k_rec_x(i,igu,i_ch)
                ENDDO
              CASE(5)
                DO i=1,nlevels
                  profiles_k_rec((i_p-1)*n_ch+i_ch)%co(i)= &
                  profiles_k_rec_x(i,igu,i_ch)
                ENDDO
              CASE(6)
                DO i=1,nlevels
                  profiles_k_rec((i_p-1)*n_ch+i_ch)%ch4(i)= &
                  profiles_k_rec_x(i,igu,i_ch)
                ENDDO
              CASE(7)
                DO i=1,nlevels
                  profiles_k_rec((i_p-1)*n_ch+i_ch)%so2(i)= &
                  profiles_k_rec_x(i,igu,i_ch)
                ENDDO
              END SELECT
            ENDDO

            profiles_k_rec((i_p-1)*n_ch+i_ch)%skin%t= &
            profiles_k_rec_s(i_ch,1)
            profiles_k_rec((i_p-1)*n_ch+i_ch)%s2m%t= &
            profiles_k_rec_s(i_ch,3)
            profiles_k_rec((i_p-1)*n_ch+i_ch)%s2m%q= &
            profiles_k_rec_s(i_ch,4)
            profiles_k_rec((i_p-1)*n_ch+i_ch)%s2m%o=0.0_jprb
            profiles_k_rec((i_p-1)*n_ch+i_ch)%s2m%p= &
            profiles_k_rec_s(i_ch,5)
            profiles_k_rec((i_p-1)*n_ch+i_ch)%s2m%u= &
            profiles_k_rec_s(i_ch,6)
            profiles_k_rec((i_p-1)*n_ch+i_ch)%s2m%v= &
            profiles_k_rec_s(i_ch,7)
            profiles_k_rec((i_p-1)*n_ch+i_ch)%cfraction= &
            profiles_k_rec_s(i_ch,8)
            profiles_k_rec((i_p-1)*n_ch+i_ch)%ctp= &
            profiles_k_rec_s(i_ch,9)

          ENDDO ! i_ch

        ENDIF ! switchrad

      ENDIF ! do_k

    ENDIF ! reconstruct

  ENDDO ! i_p (main loop over profiles)

  IF(opts%htfrtc_opts%overcast .AND. opts%htfrtc_opts%reconstruct) THEN
    WHERE(pccomp%overcast_pccomp<0.0_jprb) pccomp%overcast_pccomp=0.0_jprb
  ENDIF

  CATCH

  END SUBROUTINE htfrtc_interface

END MODULE rttov_htfrtc_interface_mod
