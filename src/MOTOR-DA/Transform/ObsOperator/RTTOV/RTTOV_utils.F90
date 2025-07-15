!!--------------------------------------------------------------------------------------------------
! PROJECT           : rttov_utils
! AFFILIATION       : GBA-MWF(SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu, 2021/8/27
! Information       : Add interfaces to call rttov N.L..
!                   : Only the clear-sky capability is enabled so far.
!   Modified by Yuanfu Xie, 2022/05/13
!                   : Add a routine for reading in RTTOV default profiles as backgrounds
!!--------------------------------------------------------------------------------------------------
!> @brief
!!
MODULE rttov_utils_m
  USE kinds_m, ONLY: r_kind, i_kind
  USE parameters_m
  USE State_m, ONLY: State_t
  USE Field_m, ONLY: Field_t
  USE ObsAttr_m, ONLY: ObsAttrSAT_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE rttov_typedefine_m

  USE rttov_TYPES, ONLY: rttov_options, rttov_profile, rttov_coefs, &
                         rttov_radiance, rttov_transmission, rttov_emissivity, &
                         rttov_reflectance, rttov_chanprof, rttov_predictors
  USE rttov_CONST
  USE parkind1, ONLY: jpim
  USE rttov_unix_env, ONLY: rttov_exit
  USE YAMLRead_m

  TYPE rttov_utils_t
    TYPE(rttov_coefs), ALLOCATABLE   :: coefs
    TYPE(rttov_options), ALLOCATABLE :: opts  ! use rttov defined options type
    ! Note: rttov_options includes rt_all, rt_ir, rt_mw, config, interpolation
    ! So: use this derived type (opts) in the following rttov_opts_all, rttov_opts_ir, etc
    TYPE(rttov_setup_input_t), ALLOCATABLE :: rttov_setup_input
    TYPE(rttov_in_calc_t), ALLOCATABLE     :: rttov_in_calc
    TYPE(rttov_diag_output), ALLOCATABLE   :: diag
    CHARACTER(LEN=1024)       :: configFile
    TYPE(State_t), POINTER :: X

  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: initialize
    FINAL :: destructor
    PROCEDURE, PUBLIC :: rttov_opts_all
    PROCEDURE, PUBLIC :: rttov_opts_ir
    PROCEDURE, PUBLIC :: rttov_opts_mw
    PROCEDURE, PUBLIC :: rttov_opts_vis_nir
    PROCEDURE, PUBLIC :: read_coefs   ! Read clr/aerosol/cld coefficient files
    PROCEDURE, PUBLIC :: rttov_GetMVars
    PROCEDURE, PUBLIC :: rttov_GetMVars_obsloc
    PROCEDURE, PUBLIC :: rttov_GetOption
    PROCEDURE, PUBLIC :: rttov_MVar2RttovVar ! Assign model (bkg) variables to rttov variables
    PROCEDURE, PUBLIC :: rttov_TLVars
    PROCEDURE, PUBLIC :: USStandardProfile2Bkgd
    PROCEDURE, PUBLIC :: rttovProfile2Bkgd

    PROCEDURE, PUBLIC :: alloc_direct
    PROCEDURE, PUBLIC :: alloc_k
    PROCEDURE, PUBLIC :: zero_k
    PROCEDURE, PUBLIC :: alloc_tl
    PROCEDURE, PUBLIC :: zero_tl
    PROCEDURE, PUBLIC :: alloc_ad
    PROCEDURE, PUBLIC :: zero_ad
    PROCEDURE, PUBLIC :: alloc_profs
    PROCEDURE, PUBLIC :: alloc_profs_k
    PROCEDURE, PUBLIC :: alloc_profs_tl
    PROCEDURE, PUBLIC :: alloc_profs_ad
    PROCEDURE, PUBLIC :: init_emiss
    PROCEDURE, PUBLIC :: Convert_to_diag_out
    PROCEDURE, PUBLIC :: rttov_dealloc_array
    PROCEDURE, PUBLIC :: rttov_dealloc_types

  END TYPE rttov_utils_t

CONTAINS

  SUBROUTINE initialize(this, configFile, X)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X

    this%configFile = configFile
    this%X => X
    IF (.NOT. ALLOCATED(this%coefs)) ALLOCATE (this%coefs)
    IF (.NOT. ALLOCATED(this%opts)) ALLOCATE (this%opts)
    IF (.NOT. ALLOCATED(this%rttov_setup_input)) THEN
      ALLOCATE (this%rttov_setup_input)
      IF (.NOT. ALLOCATED(this%rttov_setup_input%rttov_in_fixvars)) ALLOCATE (this%rttov_setup_input%rttov_in_fixvars)
      IF (.NOT. ALLOCATED(this%rttov_setup_input%rttov_in_cldvars)) ALLOCATE (this%rttov_setup_input%rttov_in_cldvars)
      IF (.NOT. ALLOCATED(this%rttov_setup_input%rttov_in_gasvars)) ALLOCATE (this%rttov_setup_input%rttov_in_gasvars)
      IF (.NOT. ALLOCATED(this%rttov_setup_input%rttov_in_userdefs)) ALLOCATE (this%rttov_setup_input%rttov_in_userdefs)
    END IF
    IF (.NOT. ALLOCATED(this%rttov_in_calc)) THEN
      ALLOCATE (this%rttov_in_calc)
!       ASSOCIATE(calc=>this%rttov_in_calc)
!         ALLOCATE(calc%radiance, calc%radiance_k, calc%radiance2, calc%transm, calc%transm_k)
! !        ALLOCATE(calc%profs(1), calc%profs_K(1))
!       END ASSOCIATE
    END IF
    IF (.NOT. ALLOCATED(this%diag)) ALLOCATE (this%diag)

  END SUBROUTINE initialize

  SUBROUTINE rttov_GetMVars(this, X, locX)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    INTEGER(i_kind), DIMENSION(2), INTENT(IN) :: locX

    ! local variables
    INTEGER(i_kind) :: i, istat, status, nlevs

    ASSOCIATE (dvars => this%rttov_setup_input%rttov_in_fixvars, &
               cldvars => this%rttov_setup_input%rttov_in_cldvars, &
               gasvars => this%rttov_setup_input%rttov_in_gasvars, &
               userdefs => this%rttov_setup_input%rttov_in_userdefs)

      ! dvars%pres = X%fields(X%getVarIdx('pres'))%data(:, locX(1), locX(2))
      dvars%pres = X%sg%FGPres(:, locX(1), locX(2))
      dvars%temp = X%fields(X%getVarIdx('temp'))%DATA(:, locX(1), locX(2))
      dvars%qvapor = X%fields(X%getVarIdx('qvapor'))%DATA(:, locX(1), locX(2))
      dvars%uwind = X%fields(X%getVarIdx('uwnd'))%DATA(:, locX(1), locX(2))
      dvars%vwind = X%fields(X%getVarIdx('vwnd'))%DATA(:, locX(1), locX(2))
      dvars%elevation = X%sg%topo(locX(1))
      dvars%tsfc = X%sg%tskin(locX(1), locX(2))
      dvars%tskin = X%sg%tskin(locX(1), locX(2))
      dvars%psfc = X%sg%psfc(locX(1), locX(2))
      dvars%u10m = X%sg%u10m(locX(1), locX(2))
      dvars%v10m = X%sg%v10m(locX(1), locX(2))
      dvars%landmask = X%sg%landmask(locX(1))
      dvars%soiltype = X%sg%soil_type(locX(1))
      dvars%snowc = X%sg%snowc(locX(1), locX(2))

      IF (this%opts%rt_ir%addclouds) THEN
        nlevs = SIZE(dvars%pres,1)
        ALLOCATE(cldvars%cldfrac(nlevs), cldvars%qcloud(nlevs), cldvars%qrain(nlevs), &
        cldvars%qice(nlevs), cldvars%qsnow(nlevs), cldvars%qgraup(nlevs))
        cldvars%qcloud = X%fields(X%getVarIdx('qcloud'))%data(:, locX(1), locX(2))
        cldvars%qrain = X%fields(X%getVarIdx('qrain'))%data(:, locX(1), locX(2))
        cldvars%qice = X%fields(X%getVarIdx('qice'))%data(:, locX(1), locX(2))
        cldvars%qsnow = X%fields(X%getVarIdx('qsnow'))%data(:, locX(1), locX(2))
        cldvars%qgraup = X%fields(X%getVarIdx('qgraupel'))%data(:, locX(1), locX(2))
        ! cldvars%cldfrac = 0.9999
        CALL Calc_cldfra(cldvars%cldfrac, dvars%qvapor, cldvars%qcloud, &
        cldvars%qice, cldvars%qsnow, dvars%temp, dvars%pres)
      END IF
  
    END ASSOCIATE

  END SUBROUTINE rttov_GetMVars

  SUBROUTINE rttov_GetMVars_obsloc(this, pres, t, q, u, v, tskin, landmask, elevation, soil_type, snowc, &
                                   qc, qr, qi, qs, qg)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    REAL(r_kind), INTENT(IN) :: pres(:), t(:), q(:), u(:), v(:)
    REAL(r_kind), INTENT(IN) :: tskin, landmask, soil_type, snowc, elevation
    REAL(r_kind), INTENT(IN), OPTIONAL :: qc(:), qr(:), qi(:), qs(:), qg(:)

    ! local variables
    INTEGER(i_kind) :: i, istat, status, nlevs

    nlevs = SIZE(pres, 1)
    ASSOCIATE (dvars => this%rttov_setup_input%rttov_in_fixvars, &
               cldvars => this%rttov_setup_input%rttov_in_cldvars, &
               gasvars => this%rttov_setup_input%rttov_in_gasvars, &
               userdefs => this%rttov_setup_input%rttov_in_userdefs)

      ALLOCATE (dvars%pres(nlevs), dvars%temp(nlevs), dvars%qvapor(nlevs), dvars%uwind(nlevs), dvars%vwind(nlevs))
      dvars%pres = pres
      dvars%temp = t
      dvars%qvapor = q
      dvars%uwind = u
      dvars%vwind = v
      dvars%elevation = elevation
      dvars%tsfc = tskin
      dvars%tskin = tskin
      IF (dvars%pres(1) > dvars%pres(2)) THEN
        dvars%psfc = pres(1) ! bottom is on the 1st level
        dvars%u10m = u(1)
        dvars%v10m = v(1)
      ELSE
        dvars%psfc = pres(nlevs) ! TOP is on the 1st level
        dvars%u10m = u(nlevs)
        dvars%v10m = v(nlevs)
      END IF
      ! PRINT *, 'Check bottom and top levels: ', pres(1), pres(nlevs)
      dvars%landmask = landmask
      dvars%soiltype = soil_type
      dvars%snowc = snowc
      IF (this%opts%rt_ir%addclouds .AND. PRESENT(qc)) THEN
        ALLOCATE(cldvars%cldfrac(nlevs), cldvars%qcloud(nlevs), cldvars%qrain(nlevs), &
        cldvars%qice(nlevs), cldvars%qsnow(nlevs), cldvars%qgraup(nlevs))
        cldvars%qcloud = qc
        cldvars%qrain = qr
        cldvars%qice = qi
        cldvars%qsnow = qs
        cldvars%qgraup = qg
        ! cldvars%cldfrac = 0.9999
        CALL Calc_cldfra(cldvars%cldfrac, dvars%qvapor, cldvars%qcloud, &
        cldvars%qice, cldvars%qsnow, dvars%temp, dvars%pres)
      END IF

    END ASSOCIATE

  END SUBROUTINE rttov_GetMVars_obsloc

  SUBROUTINE USStandardProfile2Bkgd(this, X)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X

    ! local variables
    INTEGER :: ios
    INTEGER(KIND=jpim), PARAMETER :: iup = 20
    INTEGER(KIND=jpim) :: errorstatus,gas_units
    INTEGER(KIND=jpim) :: iprof, nprof
    INTEGER(i_kind), PARAMETER :: tlevel01 = 101
    INTEGER(i_kind) :: i, ilevel, vLevel, ivar, surftype = 0, watertype = 1, it

    REAL(r_kind) :: t2m = 292.678, q2m = 5.53, u10m = 2.0, v10m = 5.0, tskin=292.678, psfc = 1100.0, wfetc = 0.1E6, salinity = 35.0, &
                    elevation=0.0, latitude=18.884, longitude=111.875, zenangle = 0.0, azangle = 0.0, &
                    sunzenangle = 45.0, sunazangle=30.0, ctp = 500.0, cfraction = 0.0
    REAL(r_kind), ALLOCATABLE :: pres(:), temp(:), qvapor(:)

    ALLOCATE(pres(tlevel01))
    ALLOCATE(temp(tlevel01))
    ALLOCATE(qvapor(tlevel01))
    pres = (/0.00500000,	0.01606500,	0.03838300,	0.07687900,	0.13695400, &
    0.22441200,	0.34540400,	0.50637400,	0.71402300,	0.97527400, &
    1.29724000,	1.68720000,	2.15257300,	2.70089700,	3.33981200, &
    4.07703800,	4.92036200,	5.87762300,	6.95669500,	8.16548000, &
    9.51188900,	11.00383500,	12.64922300,	14.45593600,	16.43183300, &
    18.58473200,	20.92240800,	23.45258300,	26.18291800,	29.12100900, &
    32.27437800,	35.65046700,	39.25663300,	43.10014400,	47.18817100, &
    51.52778600,	56.12595300,	60.98953000,	66.12525900,	71.53976800, &
    77.23956000,	83.23101600,	89.52039000,	96.11380300,	103.01724400, &
    110.23656500,	117.77748100,	125.64556200,	133.84624000,	142.38479600, &
    151.26636600,	160.49593900,	170.07834800,	180.01827900,	190.32026000, &
    200.98866500,	212.02771200,	223.44146100,	235.23381400,	247.40851400, &
    259.96914200,	272.91912000,	286.26170600,	300.00000000,	314.13693600, &
    328.67528600,	343.61765900,	358.96650300,	374.72409800,	390.89256600, &
    407.47386100,	424.46977600,	441.88194100,	459.71182100,	477.96072200, &
    496.62978500,	515.71998900,	535.23215300,	555.16693500,	575.52483200, &
    596.30618200,	617.51116300,	639.13979700,	661.19194600,	683.66731600, &
    706.56546000,	729.88577200,	753.62749400,	777.78971600,	802.37137600, &
    827.37125900,	852.78800300,	878.62009600,	904.86588000,	931.52354900, &
    958.59115400,	986.06660100,	1013.94765500,	1042.23194000,	1070.91694000, &
    1100.0				/)

    temp = (/190.195,	203.667,	215.175,	226.809,	237.789,	&
    247.507,	256.483,	263.555,	268.955,	270.636, &
    267.045,	261.561,	256.386,	251.670,	247.321, &
    243.261,	239.551,	236.057,	232.759,	229.832, &
    228.503,	227.233,	226.133,	225.245,	224.392, &
    223.595,	222.847,	222.126,	221.425,	220.733, &
    220.065,	219.420,	218.797,	218.196,	217.614, &
    217.106,	216.700,	216.700,	216.700,	216.700, &
    216.700,	216.700,	216.700,	216.700,	216.700, &
    216.700,	216.700,	216.700,	216.700,	216.700, &
    216.700,	216.700,	216.700,	216.700,	216.700, &
    216.723,	216.757,	216.790,	218.296,	220.415, &
    222.495,	224.553,	226.585,	228.580,	230.577, &
    232.588,	234.564,	236.514,	238.474,	240.401, &
    242.296,	244.201,	246.087,	247.943,	249.783, &
    251.627,	253.442,	255.229,	257.021,	258.798, &
    260.549,	262.275,	264.015,	265.730,	267.420, &
    269.095,	270.776,	272.433,	274.067,	275.689, &
    277.314,	278.917,	280.497,	282.066,	283.643, &
    285.200,	286.736,	288.251,	289.746,	291.221, &
    292.678				/)

    qvapor = (/8.86997040E-04,	1.52427987E-03,	2.01101378E-03,	2.40599752E-03,	2.72437587E-03,	&	
    2.96236227E-03,	3.10396312E-03,	3.19368814E-03,	3.23610094E-03,	3.25981037E-03,	&
    3.25671227E-03,	3.23196162E-03,	3.19285497E-03,	3.13971111E-03,	3.10625749E-03,	&
    3.08100578E-03,	3.06248483E-03,	3.04448043E-03,	3.02081059E-03,	2.99803553E-03,	&
    2.97440450E-03,	2.95184624E-03,	2.92511340E-03,	2.89196774E-03,	2.86016092E-03,	&
    2.82977471E-03,	2.80069281E-03,	2.77267144E-03,	2.73864263E-03,	2.68478764E-03,	&
    2.64118787E-03,	2.59714199E-03,	2.54483821E-03,	2.50568261E-03,	2.47311126E-03,	&
    2.44672175E-03,	2.42271064E-03,	2.40621873E-03,	2.39238370E-03,	2.38458004E-03,	&
    2.38110201E-03,	2.38850647E-03,	2.39914434E-03,	2.42737549E-03,	2.45492999E-03,	&
    2.71901083E-03,	2.99417069E-03,	3.24509817E-03,	3.47714919E-03,	3.78040842E-03,	&
    4.95999293E-03,	6.11455674E-03,	7.57663188E-03,	9.42311177E-03,	1.12322375E-02,	&
    1.42466434E-02,	1.78602949E-02,	2.14039606E-02,	2.73154771E-02,	3.41753036E-02,	&
    4.09074844E-02,	5.42724979E-02,	7.17133280E-02,	8.88417820E-02,	1.15945509E-01,	&
    1.56046471E-01,	1.95454032E-01,	2.34255117E-01,	2.72750760E-01,	3.10602956E-01,	&
    3.47828689E-01,	4.06535547E-01,	4.70306893E-01,	5.33055556E-01,	6.01901014E-01,	&
    6.85096335E-01,	7.66994753E-01,	8.47628708E-01,	9.65099457E-01,	1.09450413E+00,	&
    1.22196301E+00,	1.34952968E+00,	1.52006224E+00,	1.68809431E+00,	1.85367915E+00,	&
    2.03383153E+00,	2.26692831E+00,	2.49671002E+00,	2.72324444E+00,	2.94770972E+00,	&
    3.17163654E+00,	3.39246541E+00,	3.61025923E+00,	3.83453885E+00,	4.08727861E+00,	&
    4.33660869E+00,	4.58259770E+00,	4.82530792E+00,	5.06480362E+00,	5.30114363E+00,	&
    5.53E+00					/)
    qvapor = qvapor / 1000.0D0

    !========== Assign them to background fields ============
    ! DO ilevel = 1, SIZE(X%sg%FGPres(:,1,1),1)
    !   X%sg%SVapor(ilevel,:) = qvapor(ilevel) / q_mr_to_ppmv
    !   X%sg%STemp(ilevel,:) = 1.0D0
    ! END DO
    DO ilevel = 1, UBOUND(X%fields(1)%data,1)
      X%sg%FGPres(ilevel, :, :) = pres(ilevel)*hPa2Pa
    END DO

    DO ivar = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      vLevel = UBOUND(X%fields(ivar)%data,1)
        WRITE(*,1) ivar,vlevel,TRIM(X%fields(ivar)%Get_Name())
    1   FORMAT('Assign default profiles to bkgd: ',2I2,A)
      SELECT CASE (TRIM(X%fields(ivar)%Get_Name()))
      CASE ("pres")
        DO ilevel = 1, vLevel
          X%fields(ivar)%data(ilevel, :, :) = pres(ilevel)*hPa2Pa
        END DO
        PRINT *, 'bkg profile pres: ', X%fields(ivar)%data(:, 1, 1)
      CASE ("temp")
        DO ilevel = 1, vLevel
          X%fields(ivar)%data(ilevel, :, :) = temp(ilevel)
        END DO
      CASE ("qvapor")
        DO ilevel = 1, vLevel
          X%fields(ivar)%data(ilevel, :, :) = qvapor(ilevel)
        END DO
        PRINT *, 'bkg profile qv: ', X%fields(ivar)%data(:, 1, 1)
      CASE ("qcloud")
        X%fields(ivar)%data = 0.0
      CASE ("qice")
        X%fields(ivar)%data = 0.0
      CASE ("qrain")
        X%fields(ivar)%data = 0.0
      CASE ("qsnow")
        X%fields(ivar)%data = 0.0
      CASE ("qgraupel")
        X%fields(ivar)%data = 0.0
      CASE ("uwnd")
        X%fields(ivar)%data = 0.0
      CASE ("vwnd")
        X%fields(ivar)%data = 0.0
      CASE ("elevation")
        X%fields(ivar)%data(1, :, :) = elevation
      CASE ("ter")
        X%fields(ivar)%data(1, :, :) = elevation
      CASE ("psfc")
        X%fields(ivar)%data(1, :, :) = psfc*hPa2Pa
      CASE ("t2m")
        X%fields(ivar)%data(1, :, :) = t2m
      CASE ("q2m")
        X%fields(ivar)%data(1, :, :) = q2m
      CASE ("u10m")
        X%fields(ivar)%data(1, :, :) = u10m
      CASE ("v10m")
        X%fields(ivar)%data(1, :, :) = v10m
      CASE ("tskin")
        X%fields(ivar)%data(1, :, :) = tskin
      CASE ("tsfc")
        X%fields(ivar)%data(1, :, :) = t2m
      CASE ("surftype")
        X%fields(ivar)%data(1, :, :) = surftype
      CASE ("latitude")
        X%fields(ivar)%data(1, :, :) = latitude*degree2radian
      CASE ("longitude")
        X%fields(ivar)%data(1, :, :) = longitude*degree2radian
      CASE ("zenangle")
        X%fields(ivar)%data(1, :, :) = zenangle*degree2radian
      CASE ("azangle")
        X%fields(ivar)%data(1, :, :) = azangle*degree2radian
      CASE ("sunzenangle")
        X%fields(ivar)%data(1, :, :) = sunzenangle*degree2radian
      CASE ("sunazangle")
        X%fields(ivar)%data(1, :, :) = sunazangle*degree2radian
      CASE ("tbb")
        X%fields(ivar)%data(1, :, :) = 0.0
      CASE default
        PRINT *, 'No variables found'
      END SELECT
    END DO

    DEALLOCATE(pres, temp, qvapor)

  END SUBROUTINE USStandardProfile2Bkgd

  SUBROUTINE rttovProfile2Bkgd(this, profile, X)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: profile
    TYPE(State_t), INTENT(INOUT) :: X

    ! local variables

    INTEGER :: ios
    INTEGER(KIND=jpim), PARAMETER :: iup = 20
    INTEGER(KIND=jpim) :: errorstatus, gas_units
    INTEGER(KIND=jpim) :: iprof, nprof
    INTEGER(i_kind), PARAMETER :: tlevel01 = 51
    INTEGER(i_kind) :: i, ilevel, vLevel, ivar, surftype, watertype, it

    REAL(r_kind) :: t2m, q2m, u10m, v10m, tskin, psfc, wfetc, salinity, &
                    elevation, latitude, longitude, zenangle, azangle, &
                    sunzenangle, sunazangle, ctp, cfraction
    REAL(r_kind), ALLOCATABLE :: pres(:), temp(:), qvapor(:), rho(:), qcloud(:), qice(:), qrain(:), qsnow(:), qgraup(:)

    ALLOCATE(pres(tlevel01))
    ALLOCATE(temp(tlevel01))
    ALLOCATE(qvapor(tlevel01))
    ALLOCATE(rho(tlevel01))
    ALLOCATE(qcloud(tlevel01))
    ALLOCATE(qice(tlevel01))
    ALLOCATE(qrain(tlevel01), qsnow(tlevel01), qgraup(tlevel01))
    qrain = 0.0; qsnow = 0.0; qgraup = 0.0
    qcloud = (/ 0.000000, 0.000000,  0.000000,  0.000000,  0.000014, &
                0.000072, 0.000155,  0.000241,  0.000297,  0.000278, &
                0.000240, 0.000304,  0.000422,  0.000354,  0.000148, &
                0.000062, 0.000074,  0.000039,  0.000000,  0.000000, &
                0.000000, 0.000000,  0.000000,  0.000000,  0.000000, &
                0.000000, 0.000000,  0.000000,  0.000000,  0.000000, &
                0.000000, 0.000000,  0.000000,  0.000000,  0.000000, &
                0.000000, 0.000000,  0.000000,  0.000000,  0.000000, &
                0.000000, 0.000000,  0.000000,  0.000000,  0.000000, &
                0.000000, 0.000000,  0.000000,  0.000000,  0.000000, 0.000000 /)
    qice = (/ 0.000000, 0.000000,  0.000000,  0.000000,  0.000000, &
              0.000000, 0.000000,  0.000000,  0.000000,  0.000000, &
              0.000000, 0.000000,  0.000000,  0.000000,  0.000000, &
              0.000000, 0.000000,  0.000050,  0.000025,  0.000031, &
              0.000013, 0.000030,  0.000010,  0.000000,  0.000000, &
              0.000000, 0.000000,  0.000000,  0.000000,  0.000000, &
              0.000000, 0.000000,  0.000000,  0.000000,  0.000000, &
              0.000000, 0.000000,  0.000000,  0.000000,  0.000000, &
              0.000000, 0.000000,  0.000000,  0.000000,  0.000000, &
              0.000000, 0.000000,  0.000000,  0.000000,  0.000000, 0.000000 /)
    rho(:) = (/ 1.278577674731118, 1.258577674731118, &
              1.221116201857548, 1.06668013445105, 1.004988652329131, 0.984988652329131, &
              0.991116201857548, 0.96668013445105, 0.894988652329131, &
              0.765570918835984, 0.718989294358248, 0.675872676546365, &
              0.634722279059472, 0.595947730459753, 0.558596192446581, &
              0.523113755735217, 0.489837846936396, 0.458657813913866, &
              0.429345984663443, 0.401848861102272, 0.376047818716946, &
              0.351754883600606, 0.328752496432394, 0.306684770312954, &
              0.285136745177303, 0.264044251415438, 0.243626508867293, &
              0.224093877744312, 0.205570193422621, 0.188081657333435, &
              0.171263843998056, 0.154619340545843, 0.138237836964938, &
              0.122717282868665, 0.10883425566547, 0.0969172318085445, &
              0.0867089009509354, 0.0778268996918154, 0.0699852119712784, &
              0.0630524486776708, 0.0568675326879061, 0.0513378383385185, &
              0.0463859795747462, 0.0419483587686818, 0.0379669143046217, &
              0.0343882520830711, 0.0311748396505062, 0.028316365440844, &
              0.0257451999192356, 0.0234265445849186, 0.0213291583443029/) !, &
             ! 0.0194488258705759, 0.0176955442726988, 0.0126955442726988, &
             ! 0.0086955442726988, 0.0026955442726988/)

    !===============================================
    !========== Read profiles == start =============
    OPEN (iup, file=TRIM(Profile), status='old', iostat=ios)
    IF (ios /= 0) THEN
      WRITE (*, *) 'error opening profile file ios= ', ios
    END IF
    CALL rttov_skipcommentline(iup, errorstatus)
    ! Read gas units for profiles
    READ (iup, *) gas_units
    CALL rttov_skipcommentline(iup, errorstatus)
    nprof = 1
    ! Loop over all profiles and read data for each one
    DO iprof = 1, nprof
      ! Read pressure (hPa), temp (K), WV, O3 (gas units ppmv or kg/kg -
      ! as read above)
      READ (iup, *) pres(:)
      CALL rttov_skipcommentline(iup, errorstatus)
      READ (iup, *) temp(:)
      CALL rttov_skipcommentline(iup, errorstatus)
      READ (iup, *) qvapor(:)
      CALL rttov_skipcommentline(iup, errorstatus)
      ! Ozone profile is commented out in input profile data
      !     READ(iup,*) profiles(iprof) % o3(:)
      !     CALL rttov_skipcommentline(iup, errorstatus)
      ! 2 meter air variables
      READ (iup, *) t2m, q2m, psfc, u10m, v10m, wfetc
      CALL rttov_skipcommentline(iup, errorstatus)
      ! Skin variables
      READ (iup, *) tskin, salinity      ! FASTEM only applies to MW instruments
      CALL rttov_skipcommentline(iup, errorstatus)
      ! Surface type and water type
      READ (iup, *) surftype, watertype
      CALL rttov_skipcommentline(iup, errorstatus)
      ! Elevation, latitude and longitude
      READ (iup, *) elevation, &
        latitude, &
        longitude
      CALL rttov_skipcommentline(iup, errorstatus)
      ! Satellite and solar angles
      READ (iup, *) zenangle, azangle, sunzenangle, sunazangle
      CALL rttov_skipcommentline(iup, errorstatus)
      ! Cloud variables for simple cloud scheme, set cfraction to 0. to
      ! turn this off (VIS/IR only)
      READ (iup, *) ctp, cfraction
      CALL rttov_skipcommentline(iup, errorstatus)
    END DO
    CLOSE (iup)
    PRINT *, 'Read profile == end'

    ! Set pressure as a fixed variable for iteration
    DO ilevel = 1, UBOUND(X%fields(1)%DATA, 1)
      X%sg%FGPres(ilevel, :, :) = pres(ilevel) * hPa2Pa
    END DO

    !========== Read profiles == end =============
    !=============================================

    !========== Assign them to background fields ============
    DO ilevel = 1, SIZE(X%sg%FGPres(:, 1, 1), 1)
      X%sg%SVapor(ilevel, :) = qvapor(ilevel) / q_mr_to_ppmv
      X%sg%STemp(ilevel, :) = 1.0D0
    END DO
    DO ivar = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      vLevel = UBOUND(X%fields(ivar)%DATA, 1)
      WRITE (*, 1) ivar, vlevel, TRIM(X%fields(ivar)%Get_Name())
1     FORMAT('Assign default profiles to bkgd: ', 2I2, A)
      SELECT CASE (TRIM(X%fields(ivar)%Get_Name()))
      CASE ("pres")
        DO ilevel = 1, vLevel
          X%fields(ivar)%DATA(ilevel, :, :) = pres(ilevel) * hPa2Pa
        END DO
      CASE ("temp")
        DO ilevel = 1, vLevel
          X%fields(ivar)%DATA(ilevel, :, :) = temp(ilevel)
        END DO
      CASE ("temp_ctl")
        DO ilevel = 1, vLevel
          X%fields(ivar)%DATA(ilevel, :, :) = temp(ilevel) / X%sg%STemp(ilevel, 1)
        END DO
      CASE ("rho")
        DO ilevel = 1, vLevel
          X%fields(ivar)%DATA(vLevel + 1 - ilevel, :, :) = rho(ilevel)
        END DO
      CASE ("uwnd")
        DO ilevel = 1, vLevel
          X%fields(ivar)%DATA(ilevel, :, :) = 2.0D0
        END DO
      CASE ("vwnd")
        DO ilevel = 1, vLevel
          X%fields(ivar)%DATA(ilevel, :, :) = 5.0D0
        END DO
      CASE ("qvapor")
        DO ilevel = 1, vLevel
          X%fields(ivar)%DATA(ilevel, :, :) = qvapor(ilevel) / q_mr_to_ppmv
        END DO
        PRINT *, 'bkg profile qv: ', X%fields(ivar)%data(:, 1, 1)
        X%fields(ivar)%data = X%fields(ivar)%data * 1.0
      CASE ("qcloud")
        DO ilevel = 1, vLevel
          IF (qcloud(ilevel) < 1.0E-12) qcloud(ilevel) = 0.0D0
          X%fields(ivar)%data(vLevel+1-ilevel, :, :) = qcloud(ilevel)
        END DO
        X%fields(ivar)%data = X%fields(ivar)%data * 10.0
        PRINT *, 'bkg profile qcloud: ', X%fields(ivar)%data(:, 1, 1)
      CASE ("qice")
        DO ilevel = 1, vLevel
          IF (qice(ilevel) < 1.0E-12) qice(ilevel) = 0.0D0
          X%fields(ivar)%data(vLevel+1-ilevel, :, :) = qice(ilevel)
        END DO
        X%fields(ivar)%data = X%fields(ivar)%data * 10.0
        PRINT *, 'bkg profile qice: ', X%fields(ivar)%data(:, 1, 1)
      CASE ("qrain")
        DO ilevel = 1, vLevel
          IF (qrain(ilevel) < 1.0E-12) qrain(ilevel) = 0.0D0
          X%fields(ivar)%data(vLevel+1-ilevel, :, :) = qrain(ilevel)
        END DO
      CASE ("qsnow")
        DO ilevel = 1, vLevel
          IF (qsnow(ilevel) < 1.0E-12) qsnow(ilevel) = 0.0D0
          X%fields(ivar)%data(vLevel+1-ilevel, :, :) = qsnow(ilevel)
        END DO
      CASE ("qgraupel")
        DO ilevel = 1, vLevel
          IF (qgraup(ilevel) < 1.0E-12) qgraup(ilevel) = 0.0D0
          X%fields(ivar)%data(vLevel+1-ilevel, :, :) = qgraup(ilevel)
        END DO
      CASE ("qvapor_ctl")
        DO ilevel = 1, vLevel
          X%fields(ivar)%DATA(ilevel, :, :) = qvapor(ilevel) / (q_mr_to_ppmv * X%sg%SVapor(ilevel, 1))
        END DO
        PRINT *, 'bkg profile qv: ', X%fields(ivar)%DATA(:, 1, 1)
      CASE ("rhov")
        DO ilevel = 1, vLevel
          X%fields(ivar)%DATA(ilevel, :, :) = qvapor(ilevel) * rho(vLevel + 1 - ilevel) / q_mr_to_ppmv
        END DO
      CASE ("rhov_ctl")
        PRINT *, 'check sv: ', MAXVAL(X%sg%SVapor), MINVAL(X%sg%SVapor)
        DO ilevel = 1, vLevel
          X%fields(ivar)%DATA(ilevel, :, :) = qvapor(ilevel) * rho(vLevel + 1 - ilevel) / q_mr_to_ppmv
          DO it = 1, X%sg%tSlots
            X%fields(ivar)%DATA(ilevel, :, it) = X%fields(ivar)%DATA(ilevel, :, it) / X%sg%s1(ilevel, :)
          END DO
        END DO
      CASE ("elevation")
        X%fields(ivar)%DATA(1, :, :) = elevation
      CASE ("ter")
        X%fields(ivar)%DATA(1, :, :) = elevation
      CASE ("psfc")
        X%fields(ivar)%DATA(1, :, :) = psfc * hPa2Pa
      CASE ("t2m")
        X%fields(ivar)%DATA(1, :, :) = t2m
      CASE ("q2m")
        X%fields(ivar)%DATA(1, :, :) = q2m
      CASE ("u10m")
        X%fields(ivar)%DATA(1, :, :) = u10m
      CASE ("v10m")
        X%fields(ivar)%DATA(1, :, :) = v10m
      CASE ("tskin")
        X%fields(ivar)%DATA(1, :, :) = tskin
      CASE ("tsfc")
        X%fields(ivar)%DATA(1, :, :) = t2m
      CASE ("surftype")
        X%fields(ivar)%DATA(1, :, :) = surftype
      CASE ("latitude")
        X%fields(ivar)%DATA(1, :, :) = latitude * degree2radian
      CASE ("longitude")
        X%fields(ivar)%DATA(1, :, :) = longitude * degree2radian
      CASE ("zenangle")
        X%fields(ivar)%DATA(1, :, :) = zenangle * degree2radian
      CASE ("azangle")
        X%fields(ivar)%DATA(1, :, :) = azangle * degree2radian
      CASE ("sunzenangle")
        X%fields(ivar)%DATA(1, :, :) = sunzenangle * degree2radian
      CASE ("sunazangle")
        X%fields(ivar)%DATA(1, :, :) = sunazangle * degree2radian
      CASE ("tbb")
        X%fields(ivar)%DATA(1, :, :) = 0.0
      CASE default
        PRINT *, 'No variables found'
      END SELECT
    END DO

    DEALLOCATE(pres, temp, qvapor, rho, qcloud, qice, qrain, qsnow, qgraup)

  END SUBROUTINE rttovProfile2Bkgd

  SUBROUTINE rttov_MVar2RttovVar(this, vLevel, SatAngles, InLatLon)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    INTEGER(i_kind), INTENT(IN)   :: vLevel
    REAL(r_kind), INTENT(IN) :: SatAngles(4), InLatLon(2)
    INTEGER(i_kind) :: stride, rttov_p_start, rttov_p_end, layer_start, layer_end
    CHARACTER(10) :: tmp

    ASSOCIATE (profs => this%rttov_in_calc%profs, &
               dvars => this%rttov_setup_input%rttov_in_fixvars, &
               cldvars => this%rttov_setup_input%rttov_in_cldvars, &
               gasvars => this%rttov_setup_input%rttov_in_gasvars, &
               userdefs => this%rttov_setup_input%rttov_in_userdefs)

      IF (userdefs%debug_dev) THEN
        PRINT *, 'profs%nlevels=', profs(1)%nlevels
        PRINT *, 'profs%nlayers=', profs(1)%nlayers
      END IF
      IF (dvars%pres(1) > dvars%pres(2) .OR. dvars%pres(2) > dvars%pres(3)) THEN
        stride = -1
        rttov_p_start = vLevel
        rttov_p_end = 1
        layer_start = vLevel -1
        layer_end = 1
      ELSE
        stride = 1
        rttov_p_start = 1
        rttov_p_end = vLevel
        layer_start = 1
        layer_end = vLevel - 1
      END IF

      ! Get 3D Pressure for RTTOV profiles
      profs(1)%p = 0.0
      profs(1)%p(rttov_p_start:rttov_p_end:stride) = dvars%pres(1:vLevel) * Pa2hPa! Note: rttov uses hPa, MOTOR-DA uses Pa

      ! Note: if model variables are not interpolated to rttov levels, they will be
      ! assumed to be layer quantities (half pressure levels).
      ! userdefs%opts%interpolation % interp_mode will be called for interpolation.

      ! Get 3D Temperature for RTTOV profiles
      profs(1)%t = 0.0
      profs(1)%t(rttov_p_start:rttov_p_end:stride) = dvars%temp(1:vLevel)  ! unit: K

      ! Get 3D qvapor for RTTOV profiles
      profs(1)%q = 0.0
      IF (gasvars%ppmv_gasunit) THEN
        profs(1)%gas_units = gas_unit_ppmv ! ppmv or ppmvdry?
        profs(1)%q(rttov_p_start:rttov_p_end:stride) = dvars%qvapor(1:vLevel) * q_mr_to_ppmv
        !IF (input_gas_unit = 'g/kg') &
        !profs(1)%q = dvars%qvapor(:)*q_mr_to_ppmv
      ELSE
        profs(1)%gas_units = gas_unit_specconc ! kg/kg
        profs(1)%q(rttov_p_start:rttov_p_end:stride) = dvars%qvapor(1:vLevel)
      END IF
      ! PRINT *, 'check q2m ', profs(1)%s2m%q
      ! PRINT *, 'utils: check dvars%q ', dvars%qvapor(1:vLevel)
      ! PRINT *, 'utils: check profs(1)%q ', profs(1)%q

      IF (this%opts%rt_ir%addclouds) THEN
        profs(1)%clwde(:) = 0
        ! profs(1) % clwde(:) = re_cw(:)   !cloud liquid water effective diameter in nlayers
        profs(1) % cloud(2:5,:)      = 0.0 ! any of indices 1-5 is enough when clw_scheme=2
        profs(1) % cloud(1,layer_start:layer_end:stride)      = cldvars % qcloud(1:vLevel-1) ! any of indices 1-5 is enough when clw_scheme=2
        profs(1) % cloud(6,layer_start:layer_end:stride)      = cldvars % qice(1:vLevel-1) ! for ice cloud concentration
        profs(1) % cfrac(layer_start:layer_end:stride)        = cldvars % cldfrac(1:vLevel-1)
        profs(1) % clw_scheme      = 2
        profs(1) % ice_scheme      = 3  ! will ignore idg & icede etc. input when this option is set to 2 or 3
        profs(1) % clwde_param     = 1
        profs(1) % icede_param     = 2 
        profs(1) % mmr_cldaer      = .true. ! indicate cloud units (T => kg/kg; F => g/m^3)
      END IF

      ! Get 3D o3/co2/clw for RTTOV profiles, currently they are set to be zero
      IF (ASSOCIATED(profs(1)%o3)) profs(1)%o3 = 0.0
      IF (ASSOCIATED(profs(1)%co2)) profs(1)%co2 = 0.0
      IF (ASSOCIATED(profs(1)%clw)) profs(1)%clw = 0.0

      ! Get 2D surface variables. All surface variables are stored in 's2m', skin variables are stored in 'skin'.
      profs(1)%skin%fastem = (/3.0, 5.0, 15.0, 0.1, 0.3/)
     ! profs(1)%s2m%p = 1007.300
      ! IF (dvars%tskin < 1.0) THEN 
      !   profs(1)%s2m%t = 286.6682D0
      ! ELSE
      !   ! profs(1)%s2m%t = dvars%tskin
      IF (rttov_p_start > rttov_p_end) THEN
        profs(1)%s2m%t = profs(1)%t(rttov_p_start)
        profs(1)%s2m%q = profs(1)%q(rttov_p_start)
        profs(1)%s2m%p = MAX(400.0D0, profs(1)%p(rttov_p_start))
        profs(1)%s2m%u = dvars%uwind(1)
        profs(1)%s2m%v = dvars%vwind(1)
      ELSE
        profs(1)%s2m%t = profs(1)%t(rttov_p_end)
        profs(1)%s2m%q = profs(1)%q(rttov_p_end)
        profs(1)%s2m%p = MAX(400.0D0, profs(1)%p(rttov_p_end))
        profs(1)%s2m%u = dvars%uwind(vLevel)
        profs(1)%s2m%v = dvars%vwind(vLevel)
      END IF
      ! END IF
      ! profs(1)%s2m%t = 386.0 ! Test only
      ! profs(1)%s2m%t = profs(1)%s2m%t + 2.0 ! Test only

      profs(1)%s2m%wfetc = 100000.  !m
      IF (profs(1)%s2m%p > 1100.0 .OR. profs(1)%s2m%p < 400.0) THEN
        PRINT *, 'check surface variables: ', profs(1)%s2m%p, profs(1)%s2m%t, profs(1)%s2m%q, profs(1)%s2m%u, profs(1)%s2m%v
      END IF

      ! surface types
      ! RTTOV surftype: land = 0, sea = 1, seaice = 2
      ! GRAPES4MOTOR: sea = 0, land = 1
      profs(1)%skin%surftype = surftype_sea ! default is sea
      IF (dvars%landmask == 0) profs(1)%skin%surftype = surftype_sea
      IF (dvars%landmask == 1) profs(1)%skin%surftype = surftype_land
      ! RTTOV watertype: fresh = 0, ocean = 1
      profs(1)%skin%watertype = 1
      profs(1)%skin%snow_fraction = dvars%snowc

      IF (dvars%tskin < 1.0) THEN
        profs(1)%skin%t = 286.6682D0
      ELSE
        ! profs(1)%skin%t = dvars%tskin
        profs(1)%skin%t = profs(1)%t(rttov_p_start)
      END IF
      ! profs(1)%skin%t = 386.0 ！ Test only
      ! profs(1)%skin%t = profs(1)%skin%t + 2.0 ! Test only
      profs(1)%skin%salinity = 35.0

      ! Get geometry info from obs file
      profs(1)%elevation = dvars%elevation / km2m!0.162  ! Note: unit should be km
      ! profs(1)%elevation = 0.0 !0.162  ! Note: unit should be km
      profs(1)%latitude = InLatLon(1) / degree2radian
      profs(1)%longitude = InLatLon(2) / degree2radian
      profs(1)%zenangle = SatAngles(1) / degree2radian
      profs(1)%azangle = SatAngles(2) / degree2radian
      profs(1)%sunzenangle = SatAngles(3) / degree2radian
      profs(1)%sunazangle = SatAngles(4) / degree2radian

      ! Cloud variables for simple cloud scheme, set cfraction to 0. to turn this off (VIS/IR only)
      profs(1)%ctp = 500.0  !hPa
      profs(1)%cfraction = 0.0 !
      ! cloud top pressure 500 hPa
      ! cloud fraction 0.0
      IF (userdefs%debug_dev) &
        PRINT *, 'rttov_MVar2RttovVar OVER'
    END ASSOCIATE
  END SUBROUTINE rttov_MVar2RttovVar

  SUBROUTINE rttov_TLVars(this, pres, dX, locX, vLevel)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    REAL(r_kind), INTENT(IN) :: pres(:)
    TYPE(State_t), INTENT(IN) :: dX
    INTEGER(i_kind), DIMENSION(2), INTENT(IN) :: locX
    INTEGER(i_kind), INTENT(IN)   :: vLevel
    INTEGER(i_kind) :: stride, rttov_p_start, rttov_p_end, layer_start, layer_end
    CHARACTER(10) :: tmp

    ASSOCIATE (profs_tl => this%rttov_in_calc%profs_tl)

      IF (pres(1) > pres(2) .OR. pres(2) > pres(3)) THEN
        stride = -1
        rttov_p_start = vLevel
        rttov_p_end = 1
        layer_start = vLevel -1
        layer_end = 1
      ELSE
        stride = 1
        rttov_p_start = 1
        rttov_p_end = vLevel
        layer_start = 1
        layer_end = vLevel - 1
      END IF

      ! Get 3D Temperature for RTTOV profiles
      profs_tl(1)%t = 0.0
      profs_tl(1)%t(rttov_p_start:rttov_p_end:stride) = dX%fields(dX%getVarIdx('temp'))%data(1:vLevel, locX(1), locX(2))  ! unit: K

      ! Get 3D qvapor for RTTOV profiles
      profs_tl(1)%q = 0.0
      IF (this%rttov_setup_input%rttov_in_gasvars%ppmv_gasunit) THEN
        profs_tl(1)%gas_units = gas_unit_ppmv ! ppmv or ppmvdry?
        profs_tl(1)%q(rttov_p_start:rttov_p_end:stride) = dX%fields(dX%getVarIdx('qvapor'))%data(1:vLevel, locX(1), locX(2)) * q_mr_to_ppmv
      ELSE
        profs_tl(1)%gas_units = gas_unit_specconc ! kg/kg
        profs_tl(1)%q(rttov_p_start:rttov_p_end:stride) = dX%fields(dX%getVarIdx('qvapor'))%data(1:vLevel, locX(1), locX(2))
      END IF

      IF (this%opts%rt_ir%addclouds) THEN
        profs_tl(1) % cloud(2:5,:)      = 0.0 ! 
        profs_tl(1) % cloud(1,layer_start:layer_end:stride)    = dX%fields(dX%getVarIdx('qcloud'))%data(1:vLevel-1, locX(1), locX(2))
        profs_tl(1) % cloud(6,layer_start:layer_end:stride)    = dX%fields(dX%getVarIdx('qice'))%data(1:vLevel-1, locX(1), locX(2))
      END IF

      ! PRINT *, 'rttov_TLVars OVER'
    END ASSOCIATE
  END SUBROUTINE rttov_TLVars

  SUBROUTINE rttov_GetOption(this, inst_name, platform_name)
    USE YAMLRead_m
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    CHARACTER(len=*), INTENT(IN)   :: inst_name, platform_name
    CHARACTER(LEN=1024) :: static_file
    CHARACTER(LEN=1024) :: satinfo_file
    LOGICAL :: istat

    ! local
    LOGICAL :: haspara = .TRUE.
    INTEGER(i_kind) :: nplat, iplat, i, j, ifile

    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", static_file)

    ! The following input parameters are not mandatory. If not set, the default will be used.

    ! TODO: ishas function
    ! Default values have been set in the rttov_in_userdefs_t
    !haspara=ifhas(TRIM(configFile), 'platform_name')
    ASSOCIATE (dvars => this%rttov_setup_input%rttov_in_fixvars, &
               cldvars => this%rttov_setup_input%rttov_in_cldvars, &
               gasvars => this%rttov_setup_input%rttov_in_gasvars, &
               userdefs => this%rttov_setup_input%rttov_in_userdefs)

      userdefs%platform_name = platform_name
      userdefs%inst_name = inst_name

      !haspara=ifhas(TRIM(configFile), 'ppmv_gasunit')

      haspara = .TRUE.
      IF (haspara) THEN
        ifile = yaml_get_var(TRIM(this%configFile), 'RTTOV', 'ppmv_gasunit', gasvars%ppmv_gasunit)
        IF (userdefs%debug_dev) &
          PRINT *, 'ppmv_gasunit is: ', gasvars%ppmv_gasunit
      ELSE
        gasvars%ppmv_gasunit = .FALSE.
      END IF

      !haspara=ifhas(TRIM(configFile), 'gas_name')
      haspara = .TRUE.
      IF (haspara) THEN
        ifile = yaml_get_var(TRIM(this%configFile), 'RTTOV', 'gas_name', userdefs%rttov_insts(1)%gas_name)
        IF (userdefs%debug_dev) &
          PRINT *, 'gas_name is: ', userdefs%rttov_insts(1)%gas_name
      ELSE
        PRINT *, 'No gas is given'
      END IF
      userdefs%rttov_insts(1)%ngases = SIZE(userdefs%rttov_insts(1)%gas_name)
      IF (userdefs%debug_dev) &
        PRINT *, 'rttov_in_userdefs%rttov_insts(1)%ngases=', userdefs%rttov_insts(1)%ngases

      IF (userdefs%rttov_insts(1)%ngases == 0) THEN
        gasvars%nmixed = 0
        gasvars%nwater = 0
        gasvars%nozone = 0
        gasvars%nwvcont = 0
        gasvars%nco2 = 0
        gasvars%nn2o = 0
        gasvars%nco = 0
        gasvars%nch4 = 0
        gasvars%nso2 = 0
        gasvars%npmc = 0
      END IF

      !haspara=ifhas(TRIM(configFile), 'rttov_clouds')
      IF (haspara) THEN
        ifile = yaml_get_var(TRIM(this%configFile), 'RTTOV', 'rttov_clouds', userdefs%rttov_insts(1)%rttov_clouds)
        IF (userdefs%debug_dev) &
          PRINT *, 'rttov_clouds is: ', userdefs%rttov_insts(1)%rttov_clouds
      ELSE
        userdefs%rttov_insts(1)%rttov_clouds = .FALSE.
      END IF
      IF (userdefs%rttov_insts(1)%rttov_clouds) this%opts%rt_ir%addclouds = .TRUE.

      !haspara=ifhas(TRIM(configFile), 'rttov_aerosols')
      IF (haspara) THEN
        ifile = yaml_get_var(TRIM(this%configFile), 'RTTOV', 'rttov_aerosols', userdefs%rttov_insts(1)%rttov_aerosols)
        IF (userdefs%debug_dev) &
          PRINT *, 'rttov_aerosols is: ', userdefs%rttov_insts(1)%rttov_aerosols
      ELSE
        userdefs%rttov_insts(1)%rttov_aerosols = .FALSE.
      END IF

      !haspara=ifhas(TRIM(configFile), 'use_kmodel')
      IF (haspara) THEN
        ifile = yaml_get_var(TRIM(this%configFile), 'RTTOV', 'use_kmodel', userdefs%rttov_insts(1)%use_kmodel)
        IF (userdefs%debug_dev) &
          PRINT *, 'use_kmodel is: ', userdefs%rttov_insts(1)%use_kmodel
      ELSE
        userdefs%rttov_insts(1)%use_kmodel = .FALSE.
      END IF

      !haspara=ifhas(TRIM(configFile), 'only_over_sea')
      ! IF (haspara) THEN
      !   ifile = yaml_get_var(TRIM(this%configFile), 'FY4-AGRI', 'only_over_sea', userdefs%rttov_insts(1)%only_over_sea)
      !   IF (userdefs%debug_dev) &
      !     PRINT *, 'only_over_sea is: ', userdefs%rttov_insts(1)%only_over_sea
      ! ELSE
      !   userdefs%rttov_insts(1)%only_over_sea = .False.
      ! END IF

      !haspara=ifhas(TRIM(configFile), 'debug_dev')
      IF (haspara) THEN
        ifile = yaml_get_var(TRIM(this%configFile), 'RTTOV', 'debug_dev', userdefs%debug_dev)
        IF (userdefs%debug_dev) &
          PRINT *, 'debug_dev is: ', userdefs%debug_dev
      ELSE
        userdefs%debug_dev = .FALSE.
      END IF

      IF (haspara) THEN
        ifile = yaml_get_var(TRIM(this%configFile), 'RTTOV', 'interp_method', userdefs%interp_method)
        ! PRINT *, 'interp_method is: ', userdefs%interp_method
      ELSE
        userdefs%interp_method = 1
      END IF

    END ASSOCIATE

  END SUBROUTINE rttov_GetOption

  SUBROUTINE rttov_dealloc_array(this)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    LOGICAL :: alloc_status
    INTEGER(kind=jpim)  :: errstatus

    IF (ALLOCATED(this%rttov_setup_input)) THEN
      ASSOCIATE (dvars => this%rttov_setup_input%rttov_in_fixvars, &
                 cldvars => this%rttov_setup_input%rttov_in_cldvars, &
                 gasvars => this%rttov_setup_input%rttov_in_gasvars, &
                 userdefs => this%rttov_setup_input%rttov_in_userdefs)

        alloc_status = ALLOCATED(dvars%pres)
        IF (alloc_status) DEALLOCATE (dvars%pres)

        alloc_status = ALLOCATED(dvars%temp)
        IF (alloc_status) DEALLOCATE (dvars%temp)

        alloc_status = ALLOCATED(dvars%qvapor)
        IF (alloc_status) DEALLOCATE (dvars%qvapor)

        alloc_status = ALLOCATED(dvars%uwind)
        IF (alloc_status) DEALLOCATE (dvars%uwind)

        alloc_status = ALLOCATED(dvars%vwind)
        IF (alloc_status) DEALLOCATE (dvars%vwind)

        alloc_status = ALLOCATED(cldvars%qcloud)
        IF (alloc_status) DEALLOCATE (cldvars%qcloud)

        alloc_status = ALLOCATED(cldvars%qrain)
        IF (alloc_status) DEALLOCATE (cldvars%qrain)

        alloc_status = ALLOCATED(cldvars%qice)
        IF (alloc_status) DEALLOCATE (cldvars%qice)

        alloc_status = ALLOCATED(cldvars%qsnow)
        IF (alloc_status) DEALLOCATE (cldvars%qsnow)

        alloc_status = ALLOCATED(cldvars%qgraup)
        IF (alloc_status) DEALLOCATE (cldvars%qgraup)

        alloc_status = ALLOCATED(cldvars%re_cloud)
        IF (alloc_status) DEALLOCATE (cldvars%re_cloud)

        alloc_status = ALLOCATED(cldvars%re_rain)
        IF (alloc_status) DEALLOCATE (cldvars%re_rain)

        alloc_status = ALLOCATED(cldvars%re_ice)
        IF (alloc_status) DEALLOCATE (cldvars%re_ice)

        alloc_status = ALLOCATED(cldvars%re_snow)
        IF (alloc_status) DEALLOCATE (cldvars%re_snow)

        alloc_status = ALLOCATED(cldvars%re_graup)
        IF (alloc_status) DEALLOCATE (cldvars%re_graup)

        alloc_status = ALLOCATED(cldvars%cldfrac)
        IF (alloc_status) DEALLOCATE (cldvars%cldfrac)

        alloc_status = ALLOCATED(dvars%tsoil)
        IF (alloc_status) DEALLOCATE (dvars%tsoil)

        alloc_status = ALLOCATED(dvars%qsoil)
        IF (alloc_status) DEALLOCATE (dvars%qsoil)

        alloc_status = ALLOCATED(dvars%rttov_vlevel)
        IF (alloc_status) DEALLOCATE (dvars%rttov_vlevel)

        alloc_status = ALLOCATED(gasvars%h2o)
        IF (alloc_status) DEALLOCATE (gasvars%h2o)

        alloc_status = ALLOCATED(gasvars%co2)
        IF (alloc_status) DEALLOCATE (gasvars%co2)

        alloc_status = ALLOCATED(gasvars%o3)
        IF (alloc_status) DEALLOCATE (gasvars%o3)

        alloc_status = ALLOCATED(gasvars%ch4)
        IF (alloc_status) DEALLOCATE (gasvars%ch4)

        alloc_status = ALLOCATED(gasvars%co)
        IF (alloc_status) DEALLOCATE (gasvars%co)

        alloc_status = ALLOCATED(gasvars%so2)
        IF (alloc_status) DEALLOCATE (gasvars%so2)

        alloc_status = ALLOCATED(userdefs%rttov_insts(1)%gas_name)
        IF (alloc_status) DEALLOCATE (userdefs%rttov_insts(1)%gas_name)

      END ASSOCIATE

    END IF

    IF (ALLOCATED(this%diag)) THEN
      alloc_status = ALLOCATED(this%diag%pres)
      IF (alloc_status) DEALLOCATE (this%diag%pres)

      alloc_status = ALLOCATED(this%diag%tau_levels)
      IF (alloc_status) DEALLOCATE (this%diag%tau_levels)

      alloc_status = ALLOCATED(this%diag%tau_levels_cld)
      IF (alloc_status) DEALLOCATE (this%diag%tau_levels_cld)

      alloc_status = ALLOCATED(this%diag%up)
      IF (alloc_status) DEALLOCATE (this%diag%up)

      alloc_status = ALLOCATED(this%diag%down)
      IF (alloc_status) DEALLOCATE (this%diag%down)

      alloc_status = ALLOCATED(this%diag%surf)
      IF (alloc_status) DEALLOCATE (this%diag%surf)

      alloc_status = ALLOCATED(this%diag%t)
      IF (alloc_status) DEALLOCATE (this%diag%t)

      alloc_status = ALLOCATED(this%diag%q)
      IF (alloc_status) DEALLOCATE (this%diag%q)

      alloc_status = ALLOCATED(this%diag%jac_t)
      IF (alloc_status) DEALLOCATE (this%diag%jac_t)

      alloc_status = ALLOCATED(this%diag%jac_q)
      IF (alloc_status) DEALLOCATE (this%diag%jac_q)

      alloc_status = ALLOCATED(this%diag%cldfrac) 
      IF (alloc_status) DEALLOCATE(this%diag%cldfrac)

      alloc_status = ALLOCATED(this%diag%dtaudlnp) 
      IF (alloc_status) DEALLOCATE(this%diag%dtaudlnp)

      alloc_status = ALLOCATED(this%diag%outvar)
      IF (alloc_status) DEALLOCATE (this%diag%outvar)

    END IF

  END SUBROUTINE rttov_dealloc_array

  SUBROUTINE rttov_dealloc_types(this)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    LOGICAL :: alloc_status
    INTEGER(kind=jpim)  :: errstatus

    INCLUDE 'rttov_dealloc_coefs.interface'

    IF (ALLOCATED(this%rttov_in_calc)) THEN
!      alloc_status = ALLOCATED(this%rttov_in_calc%profs)
!      IF (alloc_status) DEALLOCATE(this%rttov_in_calc%profs)

!      alloc_status = ALLOCATED(this%rttov_in_calc%profs_K)
!      IF (alloc_status) DEALLOCATE(this%rttov_in_calc%profs_K)

      ! alloc_status = ALLOCATED(this%rttov_in_calc%radiance)
      ! IF (alloc_status) DEALLOCATE(this%rttov_in_calc%radiance)

      ! alloc_status = ALLOCATED(this%rttov_in_calc%radiance_k)
      ! IF (alloc_status) DEALLOCATE(this%rttov_in_calc%radiance_k)

      ! alloc_status = ALLOCATED(this%rttov_in_calc%radiance2)
      ! IF (alloc_status) DEALLOCATE(this%rttov_in_calc%radiance2)

      ! alloc_status = ALLOCATED(this%rttov_in_calc%transm)
      ! IF (alloc_status) DEALLOCATE(this%rttov_in_calc%transm)

      ! alloc_status = ALLOCATED(this%rttov_in_calc%transm_k)
      ! IF (alloc_status) DEALLOCATE(this%rttov_in_calc%transm_k)

    END IF

    IF (ALLOCATED(this%rttov_setup_input)) THEN
      IF (ALLOCATED(this%rttov_setup_input%rttov_in_fixvars)) DEALLOCATE (this%rttov_setup_input%rttov_in_fixvars)
      IF (ALLOCATED(this%rttov_setup_input%rttov_in_cldvars)) DEALLOCATE (this%rttov_setup_input%rttov_in_cldvars)
      IF (ALLOCATED(this%rttov_setup_input%rttov_in_gasvars)) DEALLOCATE (this%rttov_setup_input%rttov_in_gasvars)
      IF (ALLOCATED(this%rttov_setup_input%rttov_in_userdefs)) DEALLOCATE (this%rttov_setup_input%rttov_in_userdefs)
      DEALLOCATE (this%rttov_setup_input)
    END IF

    IF (ALLOCATED(this%diag)) DEALLOCATE (this%diag)

    IF (ALLOCATED(this%opts)) DEALLOCATE (this%opts)

    IF (ALLOCATED(this%rttov_in_calc)) DEALLOCATE (this%rttov_in_calc)

    IF (ALLOCATED(this%coefs)) THEN
      CALL rttov_dealloc_coefs(errstatus, this%coefs)
      DEALLOCATE (this%coefs)
    END IF

  END SUBROUTINE rttov_dealloc_types

  SUBROUTINE rttov_opts_all(this)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    ! For those frequently used options/configurations
    ! For all the rttov modules
    ! Revised to match the RTTOV 'run_example_fwd.sh' test

    this%opts%config%do_checkinput = .TRUE.
    this%opts%config%apply_reg_limits = .TRUE.
    this%opts%config%fix_hgpl = .TRUE.

    ! Note: Interpolation is required since the input vertical level is not the same
    ! with the RTTOV level
    this%opts%interpolation%addinterp = .TRUE.
    this%opts%interpolation%interp_mode = 1
    this%opts%interpolation%lgradp = .FALSE.
    this%opts%interpolation%spacetop = .TRUE.
    this%opts%interpolation%reg_limit_extrap = .TRUE.

    this%opts%rt_all%addrefrac = .TRUE. ! or False? How much does it matter?
    !this%opts%rt_all%switchrad = .False. ! Input K perturbation in radiance
    this%opts%rt_all%switchrad = .TRUE. ! Input K perturbation in BT
    this%opts%rt_all%use_t2m_opdep = .TRUE.
    this%opts%rt_all%use_q2m = .FALSE.
    this%opts%rt_all%do_lambertian = .FALSE.
    this%opts%rt_all%lambertian_fixed_angle = .TRUE.
    this%opts%rt_all%plane_parallel = .FALSE.
    this%opts%rt_all%rad_down_lin_tau = .TRUE.
    this%opts%rt_all%dtau_test = .FALSE. ! Deprecated option;

    this%opts%rt_all%ozone_data = .FALSE.
    this%opts%rt_all%co2_data = .FALSE.
    this%opts%rt_all%n2o_data = .FALSE.
    this%opts%rt_all%co_data = .FALSE.
    this%opts%rt_all%ch4_data = .FALSE.
    this%opts%rt_all%so2_data = .FALSE.

    this%opts%config%verbose = .TRUE. ! Enable priting of warnings

    !this%opts%dev%do_opdep_calc = .False.

  END SUBROUTINE rttov_opts_all

  SUBROUTINE rttov_opts_ir(this)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    ! For those frequently used options/configurations
    ! For IR channels only

    this%opts%rt_ir%solar_sea_brdf_model = 1
    this%opts%rt_ir%ir_sea_emis_model = 2
    this%opts%rt_ir%addsolar = .FALSE.
    this%opts%rt_ir%addaerosl = .FALSE.
    this%opts%rt_ir%addclouds = .FALSE.

    this%opts%rt_ir%do_nlte_correction = .FALSE.
    this%opts%rt_ir%rayleigh_single_scatt = .FALSE.
    this%opts%rt_ir%user_aer_opt_param = .FALSE.
    this%opts%rt_ir%user_cld_opt_param = .FALSE.
    this%opts%rt_ir%grid_box_avg_cloud = .FALSE.

    ! pc and htfrtc are for hi-res IR sounders only
    this%opts%rt_ir%pc%addpc = .FALSE.
    this%opts%rt_ir%pc%ipcbnd = -1
    this%opts%rt_ir%pc%ipcreg = -1
    this%opts%rt_ir%pc%npcscores = -1
    this%opts%rt_ir%pc%addradrec = .FALSE.

    this%opts%htfrtc_opts%htfrtc = .FALSE.
    this%opts%htfrtc_opts%simple_cloud = .FALSE.  ! NOTE: USE WITH CAUTION

    ! For rttov clouds
    this%opts%rt_ir%ir_scatt_model = 2    !  1 => DOM; 2 => Chou-scaling、
    this%opts%rt_ir%user_cld_opt_param = .False. ! use method 01 for cloud scattering. 
    ! Method 01: use pre-defined optical properties; 
    ! Method 02: provide optical properties explicitly.

  END SUBROUTINE rttov_opts_ir

  SUBROUTINE rttov_opts_mw(this)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    ! For those frequently used options/configurations
    ! For MW channels only

    this%opts%rt_mw%fastem_version = 6
    this%opts%rt_mw%clw_data = .False.
    this%opts%rt_mw%clw_scheme = 2

    IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
      PRINT *, "This subroutine has nothing for now"

  END SUBROUTINE rttov_opts_mw

  SUBROUTINE rttov_opts_vis_nir(this)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    ! For those frequently used options/configurations
    ! For VIS/NIR channels only

    IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
      PRINT *, "This subroutine has nothing for now"

  END SUBROUTINE rttov_opts_vis_nir

  SUBROUTINE read_coefs(this)
    ! Read coefficient files for each instrument
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    ! local
    INTEGER(i_kind) :: nlevels
    INTEGER(kind=jpim)  :: errstatus, access, status
    CHARACTER(LEN=1024) :: coef_path, clr_coef_file, aero_coef_file, cld_coef_file

    INCLUDE 'rttov_read_coefs.interface'

    !TODO: directly give a clr_coef_file name in the NL file
    ! inst_name, sensor_name, etc can be read in from the coef file
    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", coef_path)

    clr_coef_file = TRIM(coef_path)//"/Satellite/rttov_coeffs"//"/rtcoef_"//TRIM(this%rttov_setup_input%rttov_in_userdefs%platform_name) &
                    //"_"//TRIM(this%rttov_setup_input%rttov_in_userdefs%inst_name)//".dat"

    IF (this%opts%rt_ir%addaerosl) &
      aero_coef_file = TRIM(coef_path)//"/Satellite/rttov_coeffs"//"/scaercoef_"//TRIM(this%rttov_setup_input%rttov_in_userdefs%platform_name) &
                       //"_"//TRIM(this%rttov_setup_input%rttov_in_userdefs%inst_name)//".dat"

    IF (this%opts%rt_ir%addclouds) &
      cld_coef_file = TRIM(coef_path)//"/Satellite/rttov_coeffs"//"/sccldcoef_"//TRIM(this%rttov_setup_input%rttov_in_userdefs%platform_name) &
                      //"_"//TRIM(this%rttov_setup_input%rttov_in_userdefs%inst_name)//".dat"

    status = access(clr_coef_file, " ")
    IF (status .NE. 0) THEN
      PRINT *, "=======clr_coef_file is a mandatory file, but NOT existed======="
      PRINT *, 'clr_coef_file = ', TRIM(clr_coef_file)
      PRINT *, "STOP right now"
      STOP
    END IF
    status = access(aero_coef_file, " ")
    IF (status .NE. 0 .AND. this%opts%rt_ir%addaerosl) THEN
      PRINT *, "=======aero_coef_file is a mandatory file when aerosols are enabled, but NOT existed======="
      PRINT *, 'aero_coef_file = ', TRIM(aero_coef_file)
      PRINT *, "STOP right now"
      STOP
    END IF
    status = access(cld_coef_file, " ")
    IF (status .NE. 0 .AND. this%opts%rt_ir%addclouds) THEN
      PRINT *, "=======cld_coef_file is a mandatory file when clouds are enabled, but NOT existed======="
      PRINT *, 'cld_coef_file = ', TRIM(cld_coef_file)
      PRINT *, "STOP right now"
      STOP
    END IF

    IF (this%opts%rt_ir%addclouds) THEN
      CALL rttov_read_coefs(errstatus, &   ! out
                          this%coefs, &   ! out
                          this%opts, &   ! in
                          file_coef=clr_coef_file, &
                          file_sccld=cld_coef_file)
                          !file_scaer=aero_coef_file, &   ! in, optional
    ELSE
      CALL rttov_read_coefs(errstatus, &   ! out
                          this%coefs, &   ! out
                          this%opts, &   ! in
                          file_coef=clr_coef_file)     ! in, optional
    END IF

    IF (errstatus /= errorstatus_success) THEN
      PRINT *, "Read clr_coef_file failed"
      !  CALL rttov_exit(errstatus)
    END IF

    ! coefs%coef%nlevels = nlevels
    ! coefs%coef%nlayers = nlevels - 1
    ! IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) THEN
    IF (this%coefs%coef%nlevels == 0) THEN
      PRINT *, 'this%coefs%coef%nlevels=', this%coefs%coef%nlevels
      nlevels = this%coefs%coef%nlevels
      PRINT *, 'coef%fmv_model_ver=', this%coefs%coef%fmv_model_ver
      PRINT *, 'coef%id_sensor=', this%coefs%coef%id_sensor
      PRINT *, 'nmixed=', this%coefs%coef%nmixed
      PRINT *, 'nwater=', this%coefs%coef%nwater
      PRINT *, 'nozone=', this%coefs%coef%nozone
      PRINT *, 'nwvcont=', this%coefs%coef%nwvcont
      PRINT *, 'nco2=', this%coefs%coef%nco2
      PRINT *, 'coef%fmv_gas_id=', this%coefs%coef%fmv_gas_id
      PRINT *, 'coef%fmv_gas_pos=', this%coefs%coef%fmv_gas_pos
      PRINT *, 'coef%fmv_var=', this%coefs%coef%fmv_var
      PRINT *, 'coef%fmv_lvl=', this%coefs%coef%fmv_lvl
      PRINT *, 'coef%fmv_coe=', this%coefs%coef%fmv_coe
      PRINT *, 'coef%fmv_ncorr=', this%coefs%coef%fmv_ncorr
      PRINT *, 'coefs%coef%fmv_chn = ', this%coefs%coef%fmv_chn
    END IF
    ! this%coefs%coef%id_sensor = sensor_id_ir

  END SUBROUTINE read_coefs

  SUBROUTINE alloc_direct(this, asw, nprofs, nchans, nlevels, chans_rttov)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    INTEGER(i_kind), INTENT(IN) :: asw, nprofs, nlevels, nchans
    INTEGER(i_kind), INTENT(IN) :: chans_rttov
    ! local
    INTEGER(kind=jpim)  :: errstatus
    INTEGER(i_kind) :: nchanprof, nch, j, jch
    LOGICAL :: init

    INCLUDE 'rttov_alloc_direct.interface'

    nchanprof = nchans * nprofs
    this%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%chans = chans_rttov

    IF (asw > 0) THEN
      init = .TRUE.
    ELSE
      init = .FALSE.
    END IF

    !CALL mpddGlob%barrier
    IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
      PRINT *, 'rttov_utils alloc_direct start'
    ! Note to allocate profiles before this CALL
    CALL rttov_alloc_direct(errstatus, &
                            asw, &
                            nprofs, &
                            nchanprof, &
                            nlevels, &
                            this%rttov_in_calc%chanprof, &
                            this%opts, &
                            this%rttov_in_calc%profs, &        ! profs have already been allocated in alloc_profs
                            this%coefs, &
                            transmission=this%rttov_in_calc%transm, &
                            radiance=this%rttov_in_calc%radiance, &
                            radiance2=this%rttov_in_calc%radiance2, &
                            calcemis=this%rttov_in_calc%calcemis, &
                            emissivity=this%rttov_in_calc%emis, &
                            calcrefl=this%rttov_in_calc%calcrefl, &
                            reflectance=this%rttov_in_calc%reflectance, &
                            init=init)
    IF (errstatus /= errorstatus_success) THEN
      PRINT *, 'alloc_direct error = ', errstatus
      !  CALL rttov_exit(errstatus)
    END IF

    IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
      PRINT *, 'rttov_utils alloc_direct 1'
    IF (asw == 1) THEN
      this%rttov_in_calc%emis(:)%emis_in = -1.0
      this%rttov_in_calc%emis(:)%emis_out = -1.0
      this%rttov_in_calc%calcemis(:) = .FALSE.

      IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
        PRINT *, 'rttov_utils alloc_direct 2 nprofs,nchans:', nprofs, nchans
      ! Build the list of profiles/channel indices in chanprof
      nch = 0
      DO j = 1, nprofs
        DO jch = 1, nchans
          nch = nch + 1
          this%rttov_in_calc%chanprof(nch)%prof = j
          this%rttov_in_calc%chanprof(nch)%chan = this%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%chans(jch)
        END DO
      END DO
      ! Example:
      ! chanprof(1:6)%prof = (/1,1,1,2,2,2/)
      ! chanprof(1:6)%prof = (/1,2,3,1,2,3/)
      ! Interpretation: 2 different profiles, 3 different channels, 2x3 simulations
      ! The profiles must be specified in ascending order with all channels to simulate
      ! for profile 1 specified first, followed by all channels for profile 2, and so on.

    END IF

    IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
      PRINT *, 'rttov_utils alloc_direct over'

  END SUBROUTINE alloc_direct

  SUBROUTINE alloc_k(this, asw, nprofs, nchans, nlevels, chans_rttov)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    INTEGER(i_kind), INTENT(IN) :: asw, nprofs, nlevels, nchans
    INTEGER(i_kind), INTENT(IN) :: chans_rttov
    ! local
    INTEGER(i_kind) :: errstatus
    INTEGER(i_kind) :: nchanprof, nch, j, jch
    LOGICAL :: init

    INCLUDE 'rttov_alloc_k.interface'

    nchanprof = nchans * nprofs
    this%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%chans = chans_rttov

    IF (asw > 0) THEN
      init = .TRUE. !Important: adjoint variables must be initialized
    ELSE
      init = .FALSE.
    END IF

    IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
      PRINT *, 'rttov_utils alloc_k start'
    CALL rttov_alloc_k(errstatus, &
                       asw, &
                       nprofs, &
                       nchanprof, &
                       nlevels, &
                       this%rttov_in_calc%chanprof, &
                       this%opts, &
                       this%rttov_in_calc%profs, &
                       this%rttov_in_calc%profs_k, &            ! profs_k have already been allocated in alloc_profs_k
                       this%coefs, &
                       transmission=this%rttov_in_calc%transm, &
                       transmission_k=this%rttov_in_calc%transm_k, &
                       radiance=this%rttov_in_calc%radiance, &
                       radiance_k=this%rttov_in_calc%radiance_k, &
                       radiance2=this%rttov_in_calc%radiance2, &
                       calcemis=this%rttov_in_calc%calcemis, &
                       emissivity=this%rttov_in_calc%emis, &
                       emissivity_k=this%rttov_in_calc%emis_k, &
                       calcrefl=this%rttov_in_calc%calcrefl, &
                       reflectance=this%rttov_in_calc%reflectance, &
                       reflectance_k=this%rttov_in_calc%reflectance_k, &
                       init=init)

    IF (errstatus /= errorstatus_success) THEN
      PRINT *, ' alloc_k error = ', errstatus
      !  CALL rttov_exit(errstatus)
    END IF

    ! Build the list of profiles/channel indices in chanprof
    IF (asw > 0) THEN
      nch = 0
      DO j = 1, nprofs
        DO jch = 1, nchans
          nch = nch + 1
          this%rttov_in_calc%chanprof(nch)%prof = j
          this%rttov_in_calc%chanprof(nch)%chan = this%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%chans(jch)
        END DO
      END DO

      IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) THEN
        PRINT *, 'rttov_utils chanprof prof', this%rttov_in_calc%chanprof(1)%prof
        PRINT *, 'rttov_utils chanprof chan', this%rttov_in_calc%chanprof(1)%chan
        PRINT *, 'rttov_utils alloc_k over'
      END IF

      this%rttov_in_calc%emis_k(:)%emis_out = 0
      this%rttov_in_calc%emis_k(:)%emis_in = 0
      this%rttov_in_calc%emis(:)%emis_out = 0
      ! Input bt of K model is 1
      this%rttov_in_calc%radiance_k%bt = 1
      this%rttov_in_calc%radiance_k%total = 1

    END IF

    IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
      PRINT *, 'rttov_utils alloc_k over'

  END SUBROUTINE alloc_k

  SUBROUTINE zero_k(this)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this

    INCLUDE 'rttov_init_prof.interface'

    CALL rttov_init_prof(this%rttov_in_calc%profs_k)
    this%rttov_in_calc%emis_k(:)%emis_in = 0.0
    this%rttov_in_calc%emis_k(:)%emis_out = 0.0
    this%rttov_in_calc%emis(:)%emis_out = 0.0
    this%rttov_in_calc%radiance_k%bt(:) = 1.0
    this%rttov_in_calc%radiance_k%total(:) = 1.0

  END SUBROUTINE zero_k

  SUBROUTINE alloc_tl(this, asw, nprofs, nchans, nlevels, chans_rttov)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    INTEGER(i_kind), INTENT(IN) :: asw, nprofs, nlevels, nchans
    INTEGER(i_kind), INTENT(IN) :: chans_rttov
    ! local
    INTEGER(i_kind) :: errstatus
    INTEGER(i_kind) :: nchanprof, nch, j, jch
    LOGICAL :: init

    INCLUDE 'rttov_alloc_tl.interface'

    nchanprof = nchans*nprofs
    this%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%chans = chans_rttov

    IF (asw > 0) THEN
      init = .True. !Important: adjoint variables must be initialized
    ELSE
      init = .False.
    END IF

    IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
      PRINT *, 'rttov_utils alloc_tl start'
    CALL rttov_alloc_tl(errstatus, &
                       asw, &
                       nprofs, &
                       nchanprof, &
                       nlevels, &
                       this%rttov_in_calc%chanprof, &
                       this%opts, &
                       this%rttov_in_calc%profs, &
                       this%rttov_in_calc%profs_tl, &            ! profs_k have already been allocated in alloc_profs_k
                       this%coefs, &
                       transmission=this%rttov_in_calc%transm, &
                       transmission_tl=this%rttov_in_calc%transm_tl, &
                       radiance=this%rttov_in_calc%radiance, &
                       radiance_tl=this%rttov_in_calc%radiance_tl, &
                       radiance2=this%rttov_in_calc%radiance2, &
                       calcemis=this%rttov_in_calc%calcemis, &
                       emissivity=this%rttov_in_calc%emis, &
                       emissivity_tl=this%rttov_in_calc%emis_tl, &
                       calcrefl=this%rttov_in_calc%calcrefl, &
                       reflectance=this%rttov_in_calc%reflectance, &
                       reflectance_tl=this%rttov_in_calc%reflectance_tl, &
                       init=init)

    IF (errstatus /= errorstatus_success) THEN
      PRINT *, ' alloc_tl error = ', errstatus
      !  CALL rttov_exit(errstatus)
    END IF

    ! Build the list of profiles/channel indices in chanprof
    IF (asw > 0) THEN
      nch = 0
      DO j = 1, nprofs
        DO jch = 1, nchans
          nch = nch + 1
          this%rttov_in_calc%chanprof(nch)%prof = j
          this%rttov_in_calc%chanprof(nch)%chan = this%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%chans(jch)
        END DO
      END DO

      IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) THEN
        PRINT *, 'rttov_utils chanprof prof', this%rttov_in_calc%chanprof(1)%prof
        PRINT *, 'rttov_utils chanprof chan', this%rttov_in_calc%chanprof(1)%chan
        PRINT *, 'rttov_utils alloc_tl over'
      END IF

      this%rttov_in_calc%emis_tl(:)%emis_out = 0
      this%rttov_in_calc%emis_tl(:)%emis_in = 0
      this%rttov_in_calc%emis(:)%emis_out = 0
      ! Input bt of K model is 1
      this%rttov_in_calc%radiance_tl%bt = 1
      this%rttov_in_calc%radiance_tl%total = 1

    END IF

    IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
      PRINT *, 'rttov_utils alloc_tl over'

  END SUBROUTINE alloc_tl

  SUBROUTINE zero_tl(this)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this

    include 'rttov_init_prof.interface'

    call rttov_init_prof(this%rttov_in_calc%profs_tl)
    this%rttov_in_calc%emis_tl(:)%emis_in = 0.0
    this%rttov_in_calc%emis_tl(:)%emis_out = 0.0
    this%rttov_in_calc%emis(:)%emis_out = 0.0
    this%rttov_in_calc%radiance_tl%bt(:) = 1.0
    this%rttov_in_calc%radiance_tl%total(:) = 1.0

  END SUBROUTINE zero_tl

  SUBROUTINE alloc_ad(this, asw, norm_d_sp, nprofs, nchans, nlevels, chans_rttov)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    REAL(r_kind), INTENT(IN) :: norm_d_sp
    INTEGER(i_kind), INTENT(IN) :: asw, nprofs, nlevels, nchans
    INTEGER(i_kind), INTENT(IN) :: chans_rttov
    ! local
    INTEGER(i_kind) :: errstatus
    INTEGER(i_kind) :: nchanprof, nch, j, jch
    LOGICAL :: init

    INCLUDE 'rttov_alloc_ad.interface'

    nchanprof = nchans*nprofs
    this%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%chans = chans_rttov

    IF (asw > 0) THEN
      init = .True. !Important: adjoint variables must be initialized
    ELSE
      init = .False.
    END IF

    IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
      PRINT *, 'rttov_utils alloc_ad start'
    CALL rttov_alloc_ad(errstatus, &
                       asw, &
                       nprofs, &
                       nchanprof, &
                       nlevels, &
                       this%rttov_in_calc%chanprof, &
                       this%opts, &
                       this%rttov_in_calc%profs, &
                       this%rttov_in_calc%profs_ad, &            ! profs_k have already been allocated in alloc_profs_k
                       this%coefs, &
                       transmission=this%rttov_in_calc%transm, &
                       transmission_ad=this%rttov_in_calc%transm_ad, &
                       radiance=this%rttov_in_calc%radiance, &
                       radiance_ad=this%rttov_in_calc%radiance_ad, &
                       radiance2=this%rttov_in_calc%radiance2, &
                       calcemis=this%rttov_in_calc%calcemis, &
                       emissivity=this%rttov_in_calc%emis, &
                       emissivity_ad=this%rttov_in_calc%emis_ad, &
                       calcrefl=this%rttov_in_calc%calcrefl, &
                       reflectance=this%rttov_in_calc%reflectance, &
                       reflectance_ad=this%rttov_in_calc%reflectance_ad, &
                       init=init)

    IF (errstatus /= errorstatus_success) THEN
      PRINT *, ' alloc_ad error = ', errstatus
      !  CALL rttov_exit(errstatus)
    END IF

    ! Build the list of profiles/channel indices in chanprof
    IF (asw > 0) THEN
      nch = 0
      DO j = 1, nprofs
        DO jch = 1, nchans
          nch = nch + 1
          this%rttov_in_calc%chanprof(nch)%prof = j
          this%rttov_in_calc%chanprof(nch)%chan = this%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%chans(jch)
        END DO
      END DO

      IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) THEN
        PRINT *, 'rttov_utils chanprof prof', this%rttov_in_calc%chanprof(1)%prof
        PRINT *, 'rttov_utils chanprof chan', this%rttov_in_calc%chanprof(1)%chan
        PRINT *, 'rttov_utils alloc_ad over'
      END IF

      this%rttov_in_calc%emis_ad(:)%emis_out = 0
      this%rttov_in_calc%emis_ad(:)%emis_in = 0
      this%rttov_in_calc%emis(:)%emis_out = 0
      ! Input bt of K model is 1
      this%rttov_in_calc%radiance_ad%bt = norm_d_sp
      this%rttov_in_calc%radiance_ad%total = 0.0

    END IF

    IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
      PRINT *, 'rttov_utils alloc_ad over'

  END SUBROUTINE alloc_ad

  SUBROUTINE zero_ad(this)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this

    include 'rttov_init_prof.interface'

    call rttov_init_prof(this%rttov_in_calc%profs_ad)
    this%rttov_in_calc%emis_ad(:)%emis_in = 0.0
    this%rttov_in_calc%emis_ad(:)%emis_out = 0.0
    this%rttov_in_calc%emis(:)%emis_out = 0.0
    this%rttov_in_calc%radiance_ad%bt(:) = 1.0
    this%rttov_in_calc%radiance_ad%total(:) = 1.0

  END SUBROUTINE zero_ad

  SUBROUTINE alloc_profs(this, asw, nprofs, nchans, nlevels)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    INTEGER(i_kind), INTENT(IN) :: asw, nprofs, nlevels, nchans
    ! local
    INTEGER(kind=jpim)  :: errstatus
    LOGICAL :: init

    INCLUDE 'rttov_alloc_prof.interface'

    IF (asw > 0) THEN
      init = .TRUE.
    ELSE
      init = .FALSE.
    END IF

    IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
      PRINT *, ' rttov_utils alloc_profs start'
!    IF (asw == 1 .AND. (.NOT. ALLOCATED(this%rttov_in_calc%profs))) &
!    ALLOCATE (this%rttov_in_calc%profs(nprofs))

    CALL rttov_alloc_prof(errstatus, &       ! out
                          nprofs, &        ! in
                          this%rttov_in_calc%profs, &  ! inout
                          nlevels, &  !in
                          this%opts, &    ! in. Note that opts is mandatory in this CALL
                          asw, &  ! in: 1=allocate, 0=deallocate
                          coefs=this%coefs, &
                          init=init)   !in, optional. Additionally initialize profs sturcture.

    IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
      PRINT *, ' rttov_utils alloc_profs '
    IF (errstatus /= errorstatus_success) THEN
      PRINT *, 'alloc_profs error = ', errstatus
      !  CALL rttov_exit(errstatus)
    END IF

!    IF (asw == 0 .AND. ALLOCATED(this%rttov_in_calc%profS)) &
!    DEALLOCATE (this%rttov_in_calc%profs)

    IF (asw > 0) THEN
      this%rttov_in_calc%profs%skin%surftype = -1  !???
    END IF

    IF (this%rttov_setup_input%rttov_in_userdefs%debug_dev) &
      PRINT *, ' rttov_utils alloc_profs over'
  END SUBROUTINE alloc_profs

  SUBROUTINE alloc_profs_k(this, asw, nprofs, nchans, nlevels)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    INTEGER(i_kind), INTENT(IN)  :: asw, nprofs, nlevels, nchans
    ! locall
    INTEGER(kind=jpim)  :: errstatus
    LOGICAL :: init

    INCLUDE 'rttov_alloc_prof.interface'

    IF (asw > 0) THEN
      init = .TRUE.
    ELSE
      init = .FALSE.
    END IF

!    IF (asw == 1 .AND. (.NOT. ALLOCATED(this%rttov_in_calc%profs_k)) ) &
!    ALLOCATE (this%rttov_in_calc%profs_k(nprofs))

    CALL rttov_alloc_prof(errstatus, &       ! out
                          nprofs, &        ! in
                          this%rttov_in_calc%profs_k, &  ! inout
                          nlevels, &  !in
                          this%opts, &    ! in. Note that opts is mandatory in this CALL
                          asw, &  ! in: 1=allocate, 0=deallocate
                          coefs=this%coefs, &
                          init=init)   !in, optional. Additionally initialize profs sturcture.

    IF (errstatus /= errorstatus_success) THEN
      PRINT *, 'alloc_profs error = ', errstatus
      !  CALL rttov_exit(errstatus)
    END IF

!    IF (asw == 0 .AND. ALLOCATED(this%rttov_in_calc%profs_K)) &
!    DEALLOCATE (this%rttov_in_calc%profs_k)

    !TODO: how to retain rttov_in_calc%profs_k(nprofs)?

  END SUBROUTINE alloc_profs_k

  SUBROUTINE alloc_profs_tl(this, asw, nprofs, nchans, nlevels)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    INTEGER(i_kind), INTENT(IN)  :: asw, nprofs, nlevels, nchans
    ! locall
    INTEGER(kind=jpim)  :: errstatus
    LOGICAL :: init

    INCLUDE 'rttov_alloc_prof.interface'

    IF (asw > 0) THEN
      init = .True.
    ELSE
      init = .False.
    END IF

    CALL rttov_alloc_prof(errstatus, &       ! out
                          nprofs, &        ! in
                          this%rttov_in_calc%profs_tl, &  ! inout
                          nlevels, &  !in
                          this%opts, &    ! in. Note that opts is mandatory in this CALL
                          asw, &  ! in: 1=allocate, 0=deallocate
                          coefs=this%coefs, &
                          init=init)   !in, optional. Additionally initialize profs sturcture.

    IF (errstatus /= errorstatus_success) THEN
      PRINT *, 'alloc_profs error = ', errstatus
      !  CALL rttov_exit(errstatus)
    END IF

  END SUBROUTINE alloc_profs_tl

  SUBROUTINE alloc_profs_ad(this, asw, nprofs, nchans, nlevels)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    INTEGER(i_kind), INTENT(IN)  :: asw, nprofs, nlevels, nchans
    ! locall
    INTEGER(kind=jpim)  :: errstatus
    LOGICAL :: init

    INCLUDE 'rttov_alloc_prof.interface'

    IF (asw > 0) THEN
      init = .True.
    ELSE
      init = .False.
    END IF

    CALL rttov_alloc_prof(errstatus, &       ! out
                          nprofs, &        ! in
                          this%rttov_in_calc%profs_ad, &  ! inout
                          nlevels, &  !in
                          this%opts, &    ! in. Note that opts is mandatory in this CALL
                          asw, &  ! in: 1=allocate, 0=deallocate
                          coefs=this%coefs, &
                          init=init)   !in, optional. Additionally initialize profs sturcture.

    IF (errstatus /= errorstatus_success) THEN
      PRINT *, 'alloc_profs error = ', errstatus
      !  CALL rttov_exit(errstatus)
    END IF

  END SUBROUTINE alloc_profs_ad

  SUBROUTINE init_emiss(this, nchanprof)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    INTEGER(i_kind), INTENT(IN) ::  nchanprof
    ! local
    INTEGER(kind=jpim)  :: errstatus
    INTEGER(i_kind) :: ichan, i_prof

    IF (this%coefs%coef%id_sensor == sensor_id_mw) THEN
      ! PRINT *, "sensor_id_mw is activated"
      DO ichan = 1, nchanprof
        i_prof = this%rttov_in_calc%chanprof(ichan)%prof
        IF (this%rttov_in_calc%profs(i_prof)%skin%surftype == surftype_land) THEN
          this%rttov_in_calc%emis(ichan)%emis_in = 0.95
        ELSE IF (this%rttov_in_calc%profs(i_prof)%skin%surftype == surftype_seaice) THEN
          this%rttov_in_calc%emis(ichan)%emis_in = 0.92
        END IF
        this%rttov_in_calc%calcemis(ichan) = .False.
        this%rttov_in_calc%calcrefl(ichan) = .False.
      END DO
    ELSE IF (this%coefs%coef%id_sensor == sensor_id_ir .or. &
             this%coefs%coef%id_sensor == sensor_id_hi) then
      ! PRINT *, "sensor_id_ir/sensor_id_hi is activated"
      DO ichan = 1, nchanprof
        i_prof = this%rttov_in_calc%chanprof(ichan)%prof
        IF (this%rttov_in_calc%profs(i_prof)%skin%surftype == surftype_sea) THEN
          this%rttov_in_calc%emis(ichan)%emis_in = 0.0
          this%rttov_in_calc%calcemis(ichan) = .TRUE.
          this%rttov_in_calc%calcrefl(ichan) = .TRUE.
        ELSE
          IF (this%rttov_in_calc%profs(i_prof)%skin%surftype == surftype_land) THEN
            this%rttov_in_calc%emis(ichan)%emis_in = 0.95
          ELSE IF (this%rttov_in_calc%profs(i_prof)%skin%surftype == surftype_seaice) THEN
            this%rttov_in_calc%emis(ichan)%emis_in = 0.92
          END IF
          this%rttov_in_calc%calcemis(ichan) = .FALSE.
          this%rttov_in_calc%calcrefl(ichan) = .FALSE.
        END IF
      END DO
    END IF

  END SUBROUTINE init_emiss

  SUBROUTINE Calc_WeightFunc(this, diag)
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    TYPE(rttov_diag_output), INTENT(INOUT) :: diag
    INTEGER(i_kind) :: npres, ipres

    ! in RTTOV, pressure increased with level number. e.g.
    ! level-1: 0.1 hPa
    ! level-2: 1 hPa
    ! level-3: 100 hPa
    ! level-4: 1000 hPa
    npres = SIZE(diag%pres,1)
    DO ipres = npres, 2, -1
      IF (.NOT. diag%pres(ipres) - diag%pres(ipres - 1) > 0) PRINT *, 'check dtaudlnp: Negative log values'
      diag%dtaudlnp(ipres) = ABS((diag%tau_levels(ipres) - diag%tau_levels(ipres - 1)) / (LOG(diag%pres(ipres) - diag%pres(ipres - 1))))
    END DO
    diag%dtaudlnp(1) = ABS((diag%tau_levels(2) - diag%tau_levels(1)) / (LOG(diag%pres(2) - diag%pres(1))))

  END SUBROUTINE Calc_WeightFunc

  SUBROUTINE Convert_to_diag_out(this, diag)
    USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup
    USE NMLRead_m
    IMPLICIT NONE
    CLASS(rttov_utils_t) :: this
    INTEGER(i_kind) :: ifile, nlevs
    TYPE(rttov_diag_output), INTENT(INOUT) :: diag
    LOGICAL  :: write_diag
    INTEGER(i_kind) :: ilev
    CHARACTER*10 :: char_ih, char_it

    ifile = yaml_get_var(TRIM(this%configFile), 'FY4-AGRI', 'write_diag', write_diag)
    !ifile = yaml_get_var(TRIM(this%configFile), 'FY4-GIIRS', 'write_diag', write_diag)

    diag%pres = this%rttov_in_calc%profs(1)%p
    diag%bt = this%rttov_in_calc%radiance%bt(1)
    diag%btclr = this%rttov_in_calc%radiance%bt_clear(1)
    diag%tau_total = this%rttov_in_calc%transm%tau_total(1)
    diag%tau_levels = this%rttov_in_calc%transm%tau_levels(:, 1)
    diag%tau_total_cld = this%rttov_in_calc%transm%tau_total_cld(1)
    diag%tau_levels_cld = this%rttov_in_calc%transm%tau_levels_cld(:, 1)
    diag%upclr = this%rttov_in_calc%radiance2%upclear(1)
    diag%dnclr = this%rttov_in_calc%radiance2%dnclear(1)
    diag%up = this%rttov_in_calc%radiance2%up(:, 1)
    diag%down = this%rttov_in_calc%radiance2%down(:, 1)
    diag%surf = this%rttov_in_calc%radiance2%surf(:, 1)
    diag%t = this%rttov_in_calc%profs(1)%t
    IF (this%rttov_setup_input%rttov_in_gasvars%ppmv_gasunit) THEN
      diag%q = this%rttov_in_calc%profs(1)%q / q_mr_to_ppmv
    ELSE
      diag%q = this%rttov_in_calc%profs(1)%q
    END IF

    IF (this%opts%rt_ir%addclouds) THEN
      diag%cldfrac = this%rttov_in_calc%profs(1)%cfrac
      diag%qc = this%rttov_in_calc%profs(1)%cloud(1,:)
      diag%qi = this%rttov_in_calc%profs(1)%cloud(6,:)
      diag%jac_qc = this%rttov_in_calc%profs_k(1)%cloud(1,:)
      diag%jac_qi = this%rttov_in_calc%profs_k(1)%cloud(6,:)
    ELSE
      IF (.NOT. ALLOCATED(diag%cldfrac)) THEN
        ALLOCATE(diag%cldfrac(SIZE(this%rttov_in_calc%profs(1)%q)))
        ALLOCATE(diag%qc(SIZE(this%rttov_in_calc%profs(1)%q)))
        ALLOCATE(diag%qi(SIZE(this%rttov_in_calc%profs(1)%q)))
        ALLOCATE(diag%jac_qc(SIZE(this%rttov_in_calc%profs(1)%q)))
        ALLOCATE(diag%jac_qi(SIZE(this%rttov_in_calc%profs(1)%q)))
      END IF
      diag%cldfrac = 0.0
      diag%qc = 0.0
      diag%qi = 0.0
      diag%jac_qc = 0.0
      diag%jac_qi = 0.0
    END IF

    diag%jac_t = this%rttov_in_calc%profs_k(1)%t
    IF (this%rttov_setup_input%rttov_in_gasvars%ppmv_gasunit) THEN
      diag%jac_q = this%rttov_in_calc%profs_k(1)%q * q_mr_to_ppmv
    ELSE
      diag%jac_q = this%rttov_in_calc%profs_k(1)%q
    END IF

    diag%emis_out = this%rttov_in_calc%emis(1)%emis_out
    IF (.NOT. ALLOCATED(diag%dtaudlnp)) ALLOCATE (diag%dtaudlnp(SIZE(diag%tau_levels, 1)))
    CALL Calc_WeightFunc(this, diag)
    IF ( write_diag ) THEN 
      nlevs = SIZE(diag%pres,1)
      diag%nvars = 14
      IF (.NOT. ALLOCATED(diag%outvar)) ALLOCATE(diag%outvar(diag%nvars, nlevs))
      diag%outvar = 0.0D0
      ! PRINT *, 'write var with ', nlevs, 'levels'
      diag%outvar(1, 1:nlevs) = diag%pres(nlevs:1:-1)
      diag%outvar(2, 1:nlevs) = diag%t(nlevs:1:-1) 
      diag%outvar(3, 1:nlevs) = diag%q(nlevs:1:-1) 
      IF (this%opts%rt_ir%addclouds) THEN
        diag%outvar(4, 1:nlevs-1) = diag%qc(nlevs-1:1:-1) 
        diag%outvar(5, 1:nlevs-1) = diag%qi(nlevs-1:1:-1) 
        diag%outvar(6, 1:nlevs-1) = diag%cldfrac(nlevs-1:1:-1)
        diag%outvar(9, 1:nlevs-1) = diag%jac_qc(nlevs-1:1:-1) 
        diag%outvar(10, 1:nlevs-1) = diag%jac_qi(nlevs-1:1:-1) 
      ELSE
        diag%outvar(4, 1:nlevs-1) = 0.0
        diag%outvar(5, 1:nlevs-1) = 0.0
        diag%outvar(6, 1:nlevs-1) = 0.0
        diag%outvar(9, 1:nlevs-1) = 0.0
        diag%outvar(10, 1:nlevs-1) = 0.0
      END IF
      diag%outvar(7, 1:nlevs) = diag%jac_t(nlevs:1:-1) 
      diag%outvar(8, 1:nlevs) = diag%jac_q(nlevs:1:-1) 
      diag%outvar(11, 1:nlevs) = diag%tau_levels(nlevs:1:-1) 
      diag%outvar(12, 1:nlevs) = diag%dtaudlnp(nlevs:1:-1) 
      diag%outvar(13, 1:nlevs) = diag%bt
      diag%outvar(14,1:nlevs) = diag%btclr
      ! PRINT *, 'diag outvar: bt ', diag%bt
      ! PRINT *, 'pres out = ', diag%outvar(1, 1:nlevs)
      ! PRINT *, 'pres diag = ', diag%pres(nlevs:1:-1)
    END IF

    ! PRINT *, 'Write out diag data'

  END SUBROUTINE Convert_to_diag_out

  SUBROUTINE Calc_cldfra(cldfra, qv, qc, qi, qs, t_phy, p_phy)

    ! refer to Xu and Randall, 1996 
    
    IMPLICIT NONE

    REAL(r_kind), INTENT(INOUT) :: CLDFRA(:)
    REAL(r_kind), INTENT(IN) ::  QV(:), QI(:),QC(:), QS(:), t_phy(:), p_phy(:)

    INTEGER:: i,j,k, nlevs
    REAL    :: RHUM, tc, esw, esi, weight, qvsw, qvsi, qvs_weight, QIMID, QWMID, QCLD, DENOM, ARG, SUBSAT

    REAL    ,PARAMETER :: ALPHA0=100., GAMMA=0.49, QCLDMIN=1.E-12,    &
                                        PEXP=0.25, RHGRID=1.0
    REAL    , PARAMETER ::  SVP1=0.61078
    REAL    , PARAMETER ::  SVP2=17.2693882
    REAL    , PARAMETER ::  SVPI2=21.8745584
    REAL    , PARAMETER ::  SVP3=35.86
    REAL    , PARAMETER ::  SVPI3=7.66
    REAL    , PARAMETER ::  SVPT0=273.15
    REAL    , PARAMETER ::  r_d = 287.
    REAL    , PARAMETER ::  r_v = 461.6
    REAL    , PARAMETER ::  ep_2=r_d/r_v

    nlevs = SIZE(QV, 1)
    CLDFRA = 0.0
    !print *, 'check cldfra, nlevs =', nlevs
    DO k=1, nlevs
      tc   = t_phy(k) - SVPT0
      esw  = 1000.0 * SVP1 * EXP( SVP2  * tc / ( t_phy(k) - SVP3  ) )
      esi  = 1000.0 * SVP1 * EXP( SVPI2 * tc / ( t_phy(k) - SVPI3 ) )
      QVSW = EP_2 * esw / ( p_phy(k) - esw )
      QVSI = EP_2 * esi / ( p_phy(k) - esi )
      
      QCLD = QI(k)+QC(k)+QS(k)
      ! print *, 'check esw esi ', esw, esi
      ! print *, 'check QCLD ', QCLD

      IF (QCLD .LT. QCLDMIN) THEN
        weight = 0.
      ELSE
        weight = (QI(k)+QS(k)) / QCLD
      ENDIF
  
      QVS_WEIGHT = (1-weight)*QVSW + weight*QVSI
      RHUM=QV(k)/QVS_WEIGHT   !--- Relative humidity

      IF (QCLD .LT. QCLDMIN) THEN
        CLDFRA(k)=0.
      ELSEIF(RHUM.GE.RHGRID)THEN
        CLDFRA(k)=1.
      ELSE
        SUBSAT=MAX(1.E-10,RHGRID*QVS_WEIGHT-QV(k))
        DENOM=(SUBSAT)**GAMMA
        ARG=MAX(-6.9, -ALPHA0*QCLD/DENOM)    ! <-- EXP(-6.9)=.001

        RHUM=MAX(1.E-10, RHUM)
        CLDFRA(k)=(RHUM/RHGRID)**PEXP*(1.-EXP(ARG))
        IF (CLDFRA(k) .LT. .01) CLDFRA(k)=0.  
      ENDIF    
    ENDDO        
    
  END SUBROUTINE Calc_cldfra  

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(rttov_utils_t), INTENT(INOUT) :: this

    CALL rttov_dealloc_array(this)
    CALL rttov_dealloc_types(this)

  END SUBROUTINE destructor

END MODULE rttov_utils_m
