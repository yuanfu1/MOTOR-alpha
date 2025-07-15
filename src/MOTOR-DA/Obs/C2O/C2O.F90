!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE C2O_m
  USE M2ODirect_m, ONLY: M2ODirect_t
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE UV2DirSpd_m, ONLY: UV2DirSpd_t
  USE Theta2Temp_m, ONLY: Theta2Temp_t
  USE Qvapor2Rhov_m, ONLY: Qvapor2Rhov_t
  USE OprRadarRef_m, ONLY: OprRadarRef_t
  USE rttov_nl_m, ONLY: rttov_nl_t
  USE rttov_tlad_m, ONLY: rttov_tlad_t
  USE OprRadarVel_m, ONLY: OprRadarVel_t
  USE OprGNSSRefrac_m, ONLY: OprGNSSRefrac_t
  USE UV2W_m, ONLY: UV2W_t
  USE CumInterface_m, ONLY: CumInterface_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE Ctl2State_m, ONLY: Ctl2State_t
  USE GeosBal_m, ONLY: GeosBal_t
  USE parameters_m, ONLY: qv_limit, rhov_ctl_limit
  USE YAMLRead_m
  USE Satellite_utils_m
  USE UV2W_Poisson_m, ONLY: UV2W_Poisson_t

!#define TRACE_PRESSURE

  TYPE C2O_t
    TYPE(M2ODirect_t) :: M2O
    TYPE(OprRadarRef_t) :: OprRadarRef
    TYPE(rttov_nl_t), ALLOCATABLE :: OprRTTOVfwd(:)
    TYPE(rttov_tlad_t), ALLOCATABLE :: OprRTTOVtlad(:)
    TYPE(OprRadarVel_t) :: OprRadarVel
    CLASS(UV2W_t), ALLOCATABLE :: UV2W
    ! TYPE(CumInterface_t) :: CumInterface
    TYPE(Ctl2State_t) :: Ctl2State
    TYPE(CumInterface_t) :: CumInterface
    TYPE(GeosBal_t) :: GeosBal

    ! Added by Yuanfu Xie 2022-08-02:
    TYPE(UV2DirSpd_t) :: UV2DirSpd
    TYPE(OprGNSSRefrac_t):: GNSSRefrac
    TYPE(State_t), POINTER :: XbRef
    CHARACTER(len=30) :: framework = 'FullState'
    INTEGER :: n_insts
    CHARACTER(len=50), ALLOCATABLE :: platform_name(:), inst_name(:)
    LOGICAL, ALLOCATABLE :: turnOn(:)
    LOGICAL :: useUV2W = .FALSE.

  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: initialize
    FINAL :: destructor
    PROCEDURE, PUBLIC :: transFwdNonLinear_opr
    PROCEDURE, PUBLIC :: transFwdTanLinear_opr
    PROCEDURE, PUBLIC :: transAdjMultiply_opr

    GENERIC, PUBLIC :: fwdNL_opr => transFwdNonLinear_opr
    GENERIC, PUBLIC :: fwdTL_opr => transFwdTanLinear_opr
    GENERIC, PUBLIC :: adjMul_opr => transAdjMultiply_opr

    ! GENERIC :: OPERATOR(.TL.) => transFwdNonLinear_opr
    ! GENERIC :: OPERATOR(.ADJ.) => transAdjoint
    ! GENERIC :: OPERATOR(.yaml.) => transFwdNonLinear

  END TYPE C2O_t

CONTAINS

  SUBROUTINE initialize(this, configFile, X, Y, XbRef)
    IMPLICIT NONE
    CLASS(C2O_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(ObsSet_t), TARGET, INTENT(IN) :: Y
    TYPE(State_t), TARGET, INTENT(IN) :: X
    TYPE(State_t), TARGET, OPTIONAL :: XbRef
    CHARACTER(LEN=1024) :: satinfo_file
    INTEGER, ALLOCATABLE :: chan_lists(:), rttov_chan_lists(:)
    INTEGER, ALLOCATABLE :: chan_list_used(:), numsout(:)
    INTEGER :: ifile, n_plats, i_inst, nchans, ichan, nchan_used, NCH, ifuse, II
    INTEGER :: IdxStart
    LOGICAL :: istat

    IF (PRESENT(XbRef)) this%XbRef => XbRef
    ifile = yaml_get_var(TRIM(configFile), 'RunMode', 'Framework', this%framework)

    CALL this%Ctl2State%initialize(configFile)

    ! Conventional observations
    this%M2O = M2ODirect_t(configFile, X, Y)

    ! ! Control to model
    ! this%RhoCtl2Rho = RhoCtl2Rho_t(configFile, X)
    ! this%RhoRCtl2RhoR = RhoRCtl2RhoR_t(configFile, X)

    ! For option of assimilating wind speed and directions:
    IF (Y%getObsIdx('cdir') .NE. 0 .AND. &
        Y%getObsIdx('sdir') .NE. 0 .AND. &
        Y%getObsIdx('wspd') .NE. 0) THEN
      this%UV2DirSpd = UV2DirSpd_t(configFile, X, Y)
    END IF

    IF (Y%getObsIdx('refractivity') .NE. 0) THEN
      this%GNSSRefrac = OprGNSSRefrac_t(configFile, XbRef, Y)
      ! BLOCK
      ! INTEGER(i_kind) :: i
      ! DO i = LBOUND(X%Fields, 1), UBOUND(X%Fields, 1)
      !   PRINT*, 'Hyj+++, X name', i, X%Fields(i)%Get_Name()
      ! END DO
      ! END BLOCK
      ! stop
    END IF

    ! Radar operator
    IF (Y%getObsIdx('ref') .NE. 0 .AND. X%getVarIdx(TRIM('rhor_ctl')) .NE. 0) THEN
      this%OprRadarRef = OprRadarRef_t(configFile, X, Y)
    END IF

    ! Radar radial wind forward
    IF (Y%getObsIdx('rwnd') .NE. 0 &
        .AND. X%getVarIdx(TRIM('uwnd')) .NE. 0 &
        .AND. X%getVarIdx(TRIM('vwnd')) .NE. 0) THEN
      this%OprRadarVel = OprRadarVel_t(configFile, X, Y)
    END IF

    ifile = yaml_get_var(TRIM(configFile), 'C2O', 'useUV2W', this%useUV2W)
    IF (ifile /= 0) this%useUV2W = .FALSE.

    IF (this%useUV2W .AND. X%getVarIdx('wwnd') /= 0) THEN
      PRINT *, 'STOP: has wwnd before get into UV2W!'
      STOP
    END IF

    ! Radar radial wind strong constraint forward
    IF (this%useUV2W) THEN
      this%UV2W = UV2W_Poisson_t(configFile, X)
    END IF

    ! Satellite operator
    ! NOTE: for satellite DA,
    ! the instrument cycle is done in C2O, instead of in RTTOV operators. Otherwise, the required memory can be very large
    ! Each instrument has its own Hx, with the obstype name as 'fy4-1-agri' or 'fy4-2-agri'
    ! channels are treated as different ObsFields, e.g. Hx%ObsFileds(nchans)

    IdxStart = Y%getObsIdx('tbb')
    IF (Y%getObsIdx('tbb') .NE. 0) THEN
      ! Get parameters from the YAML file
      ifile = yaml_get_var(TRIM(configFile), 'RTTOV', 'inst_name', this%inst_name)
      ifile = yaml_get_var(TRIM(configFile), 'RTTOV', 'platform_name', this%platform_name)
      ifile = yaml_get_var(TRIM(configFile), 'RTTOV', 'turnOn', this%turnOn)

      n_plats = SIZE(this%platform_name, 1)
      this%n_insts = SIZE(this%inst_name, 1)
      PRINT *, 'C2O inst_name:', this%inst_name, '  ', this%n_insts, n_plats
      IF (n_plats .NE. this%n_insts) THEN
        PRINT *, 'STOP: Numbers of platform_name and inst_name are inconsistent'
        RETURN
      END IF

      ! ALLOCATE(this%OprRTTOVfwd(COUNT(this%turnOn)))
      ! ALLOCATE(this%OprRTTOVtlad(COUNT(this%turnOn)))
      ! PRINT *, COUNT(this%turnOn), ' instruments will be assimilated'

      ALLOCATE (this%OprRTTOVfwd(SIZE(this%turnOn)))
      ALLOCATE (this%OprRTTOVtlad(SIZE(this%turnOn)))
      PRINT *, COUNT(this%turnOn), 'out of ', SIZE(this%turnOn), ' instruments will be assimilated'

      ! For each instrument, channel info is determined in the 'satinfo' file, which
      ! tells us which channels will be assimilated
      DO i_inst = 1, this%n_insts

        IF (this%turnOn(i_inst)) THEN
          CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", satinfo_file)
          satinfo_file = TRIM(satinfo_file)//"/Satellite/satinfo_"// &
                         TRIM(this%platform_name(i_inst))//'-'//TRIM(this%inst_name(i_inst))//".txt"

          PRINT *, 'C2O: satinfo_file = ', TRIM(satinfo_file)

          IF (TRIM(this%inst_name(i_inst)) .NE. 'giirs') THEN
            CALL Get_rttov_chan_info(TRIM(this%platform_name(i_inst))//'-'//TRIM(this%inst_name(i_inst)), nchans)

            ALLOCATE (chan_lists(nchans), rttov_chan_lists(nchans))

            CALL Get_rttov_chan_info(TRIM(this%platform_name(i_inst))//'-'//TRIM(this%inst_name(i_inst)), nchans, chan_lists, rttov_chan_lists)
            nchan_used = nchans
          ELSE
            INQUIRE (file=satinfo_file, exist=istat)
            IF (istat) THEN
              OPEN (unit=121, file=satinfo_file, status='old')
              READ (121, *)
              READ (121, *) nchans
              ALLOCATE (chan_lists(nchans))
              ALLOCATE (chan_list_used(nchans))
              NCH = 0
              DO ichan = 1, nchans
                READ (121, *) chan_lists(ichan), II, ifuse
                IF (ifuse .GT. 0) THEN                      !! =1: use, / 0: not use
                  NCH = NCH + 1
                  chan_list_used(NCH) = chan_lists(ichan)
                  ! print *,ichan,NCH,"C20 chan_lists:",chan_lists(ichan),chan_list_used(NCH),ifuse
                END IF
              END DO
              nchan_used = NCH
              ALLOCATE (rttov_chan_lists(nchan_used))
              rttov_chan_lists(1:nchan_used) = chan_list_used(1:nchan_used)
              CLOSE (121)
            ELSE
              PRINT *, '----------satinfo file was not found in C2O----------'
              STOP
            END IF
          END IF

          ! Initialization for NL/TL/AD
          ! For a specific instrument only
          ALLOCATE (numsout(SIZE(Y%ObsFields, 1)))
          DO ii = 1, SIZE(Y%ObsFields, 1)
            numsout(ii) = SIZE(Y%ObsFields(ii)%values, 1)
          END DO
          ! PRINT *, 'myrank = ', X%sg%mpddInfo_sg%myrank, 'total number of obs into DA is: ', numsout
          ! PRINT *, 'IdxStart = ', IdxStart, ' nchan_used = ', nchan_used

          CALL this%OprRTTOVfwd(i_inst)%initialize(configFile, X, Y, IdxStart, nchan_used, rttov_chan_lists, TRIM(this%inst_name(i_inst)), TRIM(this%platform_name(i_inst)))
          CALL this%OprRTTOVtlad(i_inst)%initialize(configFile, X, Y, IdxStart, nchan_used, rttov_chan_lists, TRIM(this%inst_name(i_inst)), TRIM(this%platform_name(i_inst)))

          IdxStart = IdxStart + nchan_used

          IF (ALLOCATED(chan_lists)) DEALLOCATE (chan_lists)
          IF (ALLOCATED(rttov_chan_lists)) DEALLOCATE (rttov_chan_lists)
          DEALLOCATE (numsout)
        END IF
      END DO
    END IF

    ! Currently hard coded to test geostrophic balance requirement
    ! Consider build a tester inside GeosBal.F90
    IF (X%getVarIdx('uwnd') .GT. 0 .AND. &
        X%getVarIdx('vwnd') .GT. 0 .AND. &
        X%getVarIdx('temp') .GT. 0 .AND. &
        (X%getVarIdx('rho') .GT. 0 .OR. &
         X%getVarIdx('pres_ctl') .GT. 0 .OR. &
         X%getVarIdx('pres_log') .GT. 0)) THEN
      CALL this%GeosBal%initial(configFile, X)
      PRINT *, 'GeosBal initialized'
    END IF

    PRINT *, 'C2O is initialized...!'

  END SUBROUTINE initialize

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(C2O_t), INTENT(INOUT) :: this

    ! PRINT *, 'End of destructor of C2O.'
    IF (ALLOCATED(this%OprRTTOVfwd)) DEALLOCATE (this%OprRTTOVfwd)
    IF (ALLOCATED(this%OprRTTOVtlad)) DEALLOCATE (this%OprRTTOVtlad)
    IF (ALLOCATED(this%platform_name)) DEALLOCATE (this%platform_name)
    IF (ALLOCATED(this%inst_name)) DEALLOCATE (this%inst_name)
    IF (ALLOCATED(this%turnOn)) DEALLOCATE (this%turnOn)

  END SUBROUTINE destructor

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dY)
    IMPLICIT NONE
    CLASS(C2O_t), INTENT(IN) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t), INTENT(IN), OPTIONAL :: X
    TYPE(ObsSet_t) :: dY
    INTEGER :: i_inst, i

    TYPE(State_t) :: dXTemp, XTemp

    dXTemp = dX
    XTemp = X
    ! IF (dXTemp%getVarIdx('qvapor') /= 0) &
    !   dXTemp%Fields(dXTemp%getVarIdx('qvapor'))%data = dXTemp%Fields(dXTemp%getVarIdx('qvapor'))%data/1000.0D0
    ! IF (dXTemp%getVarIdx('pres') /= 0) &
    !   dXTemp%Fields(dXTemp%getVarIdx('pres'))%data = dXTemp%Fields(dXTemp%getVarIdx('pres'))%data*100.0D0
    ! CALL this%Ctl2State%transFwdTanLinear(XTemp, X)
    CALL this%Ctl2State%transFwdNonLinear(XTemp)
    IF (TRIM(this%framework) .EQ. 'Incremental') THEN
      CALL this%Ctl2State%SumState(XTemp, this%XbRef)
    END IF
    ! Yuanfu Xie added X to the call of transFwdTanLinear as Ctl2State has nonlinear operators:
    CALL this%Ctl2State%transFwdTanLinear(dXTemp, X)

    ! Convention observations
    dY = this%M2O%transFwdNonLinear_opr(dXTemp)

    ! Radar reflectivity forward
    ! IF (dY%getObsIdx('ref') .NE. 0 .AND. dX%getVarIdx(TRIM('rhor_ctl')) .NE. 0) THEN
    !   CALL this%RhoRCtl2RhoR%transFwdNonLinear(dXTemp)
    !   CALL this%OprRadarRef%transFwdNonLinear(dXTemp, dY)
    ! END IF
    ! GNSS RO Refractivity
    IF (dY%getObsIdx('refractivity') .NE. 0 &
        .AND. (dXTemp%getVarIdx(TRIM('pres')) .NE. 0) & !.OR. (dX%getVarIdx(TRIM('pres_ctl')) .NE. 0))  &
        .AND. (dXTemp%getVarIdx(TRIM('qvapor')) .NE. 0) &!.OR. (dX%getVarIdx(TRIM('qvapor_ctl')) .NE. 0)) &
        .AND. (dXTemp%getVarIdx(TRIM('temp')) .NE. 0) &
        ) THEN

      CALL this%GNSSRefrac%transFwdTanLinear(dXTemp, dY, XTemp)
    END IF

    IF (this%useUV2W) THEN
      CALL this%UV2W%transFwdNonLinear(dXTemp)
    END IF

    ! Radar radial wind forward
    IF (dY%getObsIdx('rwnd') .NE. 0 &
        .AND. dX%getVarIdx(TRIM('uwnd')) .NE. 0 &
        .AND. dX%getVarIdx(TRIM('vwnd')) .NE. 0) THEN
      CALL this%OprRadarVel%transFwdNonLinear(dXTemp, dY)
    END IF

    ! Geostrophic balance:
    IF (.FALSE. .AND. X%getVarIdx('uwnd') .GT. 0 .AND. &
        X%getVarIdx('vwnd') .GT. 0 .AND. &
        X%getVarIdx('temp') .GT. 0 .AND. &
        (X%getVarIdx('rho') .GT. 0 .OR. X%getVarIdx('pres') .GT. 0)) THEN

#ifdef TRACE_PRESSURE
      WRITE (*, 1) X%mpddSub%myrank
1     FORMAT('transFwdTanLinear_opr calling Geo-T:', I2)
#endif
      IF (X%sg%gLevel .GE. this%GeosBal%geobalBeg .AND. &
          X%sg%gLevel .LE. this%GeosBal%geobalEnd) &
        CALL this%GeosBal%transFwdTanLinear(dXTemp, X)
    END IF

    ! Satellite TL forward
    IF (dY%getObsIdx('tbb') .NE. 0) THEN

      WHERE (XTemp%fields(XTemp%getVarIdx('qvapor'))%DATA < qv_limit) XTemp%fields(XTemp%getVarIdx('qvapor'))%DATA = qv_limit
      IF (XTemp%getVarIdx('qice') .GT. 0) THEN
        WHERE (XTemp%fields(XTemp%getVarIdx('qcloud'))%DATA < 0.0D0) XTemp%fields(XTemp%getVarIdx('qcloud'))%DATA = 0.0D0
        WHERE (XTemp%fields(XTemp%getVarIdx('qice'))%DATA < 0.0D0) XTemp%fields(XTemp%getVarIdx('qice'))%DATA = 0.0D0
      END IF

      DO i_inst = 1, this%n_insts
        IF (this%turnOn(i_inst) .AND. &
            dY%getObsTypeIdx(TRIM(this%platform_name(i_inst))) .NE. 0 .AND. dY%getObsTypeIdx(TRIM(this%inst_name(i_inst))) .NE. 0) &
          CALL this%OprRTTOVtlad(i_inst)%rttov_tl_simobs(XTemp, dXTemp, dY)
      END DO

    END IF
  END FUNCTION transFwdTanLinear_opr

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(Y)
    IMPLICIT NONE
    CLASS(C2O_t) :: this
    TYPE(State_t), INTENT(in) :: X
    TYPE(ObsSet_t) :: Y
    INTEGER :: i_inst, i

    TYPE(State_t) :: XTemp

    XTemp = X

    ! IF (XTemp%getVarIdx('qvapor') /= 0) &
    !   XTemp%Fields(XTemp%getVarIdx('qvapor'))%data = XTemp%Fields(XTemp%getVarIdx('qvapor'))%data/1000.0D0
    ! IF (XTemp%getVarIdx('pres') /= 0) &
    !   XTemp%Fields(XTemp%getVarIdx('pres'))%data = XTemp%Fields(XTemp%getVarIdx('pres'))%data*100.0D0

    CALL this%Ctl2State%transFwdNonLinear(XTemp)
    IF (TRIM(this%framework) .EQ. 'Incremental') THEN
      CALL this%Ctl2State%SumState(XTemp, this%XbRef)
    END IF

    ! CALL this%Pres2Rho%fwdNL(XTemp)
    ! CALL this%Qvapor2Rhov%fwdNL(XTemp)
    ! CALL this%Theta2Temp%fwdNL(XTemp)

    ! ! Cumlus forward
    ! IF (XTemp%getVarIdx(TRIM('uwnd')) .NE. 0 &
    !     .AND. XTemp%getVarIdx(TRIM('vwnd')) .NE. 0 &
    !     .AND. XTemp%getVarIdx(TRIM('wwnd')) .NE. 0 &
    !     .AND. XTemp%getVarIdx(TRIM('theta')) .NE. 0 &
    !     .AND. XTemp%getVarIdx(TRIM('qvapor')) .NE. 0 &
    !     .AND. XTemp%getVarIdx(TRIM('pres')) .NE. 0 &
    !     ) THEN
    !   CALL this%CumInterface%CumForward(XTemp)
    ! END IF

    ! Convention observations
    Y = this%M2O%transFwdNonLinear_opr(XTemp)

    ! ! For directly assimilating wind direction and speed:
    ! IF (Y%getObsIdx('cdir') .NE. 0 .AND. &
    !     Y%getObsIdx('sdir') .NE. 0 .AND. &
    !     Y%getObsIdx('wspd') .NE. 0 .AND. &
    !     X%getVarIdx('uwnd') .NE. 0 .AND. &
    !     X%getVarIdx('vwnd') .NE. 0) THEN
    !   CALL this%UV2DirSpd%transFwdNonLinear(XTemp, Y)
    ! END IF

    ! Radar reflectivity forward
    ! IF (Y%getObsIdx('ref') .NE. 0 .AND. X%getVarIdx(TRIM('rhor_ctl')) .NE. 0) THEN
    !   CALL this%RhoRCtl2RhoR%transFwdNonLinear(XTemp)
    !   CALL this%OprRadarRef%transFwdNonLinear(XTemp, Y)
    ! END IF

    IF (Y%getObsIdx('refractivity') .NE. 0 &
        .AND. (XTemp%getVarIdx(TRIM('pres')) .NE. 0) & !.OR. (X%getVarIdx(TRIM('pres_ctl')) .NE. 0))  &
        .AND. (XTemp%getVarIdx(TRIM('qvapor')) .NE. 0) & !.OR. (X%getVarIdx(TRIM('qvapor_ctl')) .NE. 0)) &
        .AND. (XTemp%getVarIdx(TRIM('temp')) .NE. 0) &
        ) THEN
      CALL this%GNSSRefrac%transFwdNonLinear(XTemp, Y)
    END IF

    ! Radar radial wind strong constraint forward
    IF (this%useUV2W) THEN
      CALL this%UV2W%transFwdNonLinear(XTemp)
    END IF

    ! Radar radial wind forward
    IF (Y%getObsIdx('rwnd') .NE. 0 &
        .AND. X%getVarIdx(TRIM('uwnd')) .NE. 0 &
        .AND. X%getVarIdx(TRIM('vwnd')) .NE. 0) THEN
      CALL this%OprRadarVel%transFwdNonLinear(XTemp, Y)
    END IF

    ! ! print*,'I am out of transFwdNonlinear_opr....',X%mpddGlob%myrank

    ! Satellite NL forward
    IF (Y%getObsIdx('tbb') .NE. 0 ) THEN
    
      WHERE (XTemp%fields(XTemp%getVarIdx('qvapor'))%data < qv_limit) XTemp%fields(XTemp%getVarIdx('qvapor'))%data = qv_limit
      IF (XTemp%getVarIdx('qice') .GT. 0) THEN
        WHERE (XTemp%fields(XTemp%getVarIdx('qcloud'))%DATA < 0.0D0) XTemp%fields(XTemp%getVarIdx('qcloud'))%DATA = 0.0D0
        WHERE (XTemp%fields(XTemp%getVarIdx('qice'))%DATA < 0.0D0) XTemp%fields(XTemp%getVarIdx('qice'))%DATA = 0.0D0
      END IF
      
      IF (MAXVAL(XTemp%sg%FGPres) .GT. 115000.0 .OR. MINVAL(XTemp%sg%FGPres) .LT. 0.1)  THEN
        PRINT *, '--- STOP: pressure values exceed max/min check --- ' 
        PRINT *, 'MAX/MIN: ', X%sg%mpddInfo_sg%myrank, MAXVAL(XTemp%sg%FGPres), MINVAL(XTemp%sg%FGPres)
        STOP
      END IF

      IF (MAXVAL(XTemp%fields(XTemp%getVarIdx('qvapor'))%DATA) .GT. 0.055) THEN
        PRINT *, '--- STOP: qvapor values after minimization exceed max check ---'
        PRINT *, 'MAX: ', MAXVAL(XTemp%fields(XTemp%getVarIdx('qvapor'))%DATA)
        STOP
      END IF

      IF (MAXVAL(XTemp%fields(XTemp%getVarIdx('temp'))%DATA) .GT. 350.0 .OR. &
          MINVAL(XTemp%fields(XTemp%getVarIdx('temp'))%DATA) .LT. 150.0) THEN
        PRINT *, '--- STOP: Temp values after minimization exceed max/min check ---'
        PRINT *, 'MAX/MIN: ', MAXVAL(XTemp%fields(XTemp%getVarIdx('temp'))%DATA), MINVAL(XTemp%fields(XTemp%getVarIdx('temp'))%DATA)
        STOP
      END IF
      IF (XTemp%getVarIdx('qcloud') .GT. 0) THEN
        IF (MINVAL(XTemp%fields(XTemp%getVarIdx('qcloud'))%DATA) .LT. 0.0) THEN
          PRINT *, '--- STOP: qcloud values after minimization are negative  ---'
          PRINT *, 'MIN: ', MINVAL(XTemp%fields(XTemp%getVarIdx('qcloud'))%DATA)
          STOP
        END IF
      END IF
      IF (XTemp%getVarIdx('qice') .GT. 0) THEN
        IF (MINVAL(XTemp%fields(XTemp%getVarIdx('qice'))%DATA) .LT. 0.0) THEN
          PRINT *, '--- STOP: qice values after minimization are negative  ---'
          PRINT *, 'MIN: ', MINVAL(XTemp%fields(XTemp%getVarIdx('qice'))%DATA)
          STOP
        END IF
      END IF

      ! PRINT *, 'check NL temp: ', MAXVAL(XTemp%fields(XTemp%getVarIdx('temp'))%data), &
      ! MINVAL(XTemp%fields(XTemp%getVarIdx('temp'))%data)
      ! PRINT *, 'check NL qvapor: ', MAXVAL(XTemp%fields(XTemp%getVarIdx('qvapor'))%data), &
      ! MINVAL(XTemp%fields(XTemp%getVarIdx('qvapor'))%data)

      DO i_inst = 1, this%n_insts
        IF (this%turnOn(i_inst) .AND. &
            Y%getObsTypeIdx(TRIM(this%platform_name(i_inst))) .NE. 0 .AND. Y%getObsTypeIdx(TRIM(this%inst_name(i_inst))) .NE. 0) &
          CALL this%OprRTTOVfwd(i_inst)%rttov_nl_simobs(XTemp, Y)
      END DO
    END IF

    ! Geostrophic balance:
    IF (.FALSE. .AND. X%getVarIdx('uwnd') .GT. 0 .AND. &
        X%getVarIdx('vwnd') .GT. 0 .AND. &
        X%getVarIdx('temp') .GT. 0 .AND. &
        (X%getVarIdx('rho') .GT. 0 .OR. X%getVarIdx('pres') .GT. 0)) THEN
      IF (X%sg%gLevel .GE. this%GeosBal%geobalBeg .AND. &
          X%sg%gLevel .LE. this%GeosBal%geobalEnd) &
        CALL this%GeosBal%transFwdNonLinear(XTemp)
    END IF

    ! print*,'I am out of transFwdNonlinear_opr....',X%mpddGlob%myrank

  END FUNCTION transFwdNonLinear_opr

  FUNCTION transAdjMultiply_opr(this, D, X) RESULT(dX)
    IMPLICIT NONE
    CLASS(C2O_t) :: this
    TYPE(State_t) :: dX
    TYPE(ObsSet_t), INTENT(IN) :: D
    TYPE(State_t), INTENT(in) :: X
    TYPE(State_t) :: XTemp
    INTEGER :: i_inst

    XTemp = X

    CALL this%Ctl2State%transFwdNonLinear(XTemp)
    IF (TRIM(this%framework) .EQ. 'Incremental') THEN
      CALL this%Ctl2State%SumState(XTemp, this%XbRef)
    END IF
    dX = XTemp%zeroCopy()

    ! Conventional observation
    CALL this%M2O%transAdjMultiply(D, dX)
#ifdef TRACE_PRESSURE
    PRINT *, 'MMdX1 = ', MINVAL(dX%fields(5)%DATA(1, :, :)), MAXVAL(dX%fields(5)%DATA(1, :, :))
#endif

    ! ! For directly assimilating wind direction and speed:
    ! IF (D%getObsIdx('cdir') .NE. 0 .AND. &
    !     D%getObsIdx('sdir') .NE. 0 .AND. &
    !     D%getObsIdx('wspd') .NE. 0 .AND. &
    !     X%getVarIdx('uwnd') .NE. 0 .AND. &
    !     X%getVarIdx('vwnd') .NE. 0) THEN
    !   BLOCK
    !     TYPE(State_t) :: DUV
    !     dUV = this%UV2DirSpd%transAdjMultiply_opr(D, X)

    !     dX = dX + dUV
    !   END BLOCK
    ! END IF

    ! Radar reflectivity
    ! IF (D%getObsIdx('ref') .NE. 0 .AND. X%getVarIdx(TRIM('rhor_ctl')) .NE. 0) THEN
    !   BLOCK
    !     TYPE(State_t) :: dXR
    !     TYPE(State_t) :: xTemp

    !     XTemp = X
    !     CALL this%RhoRCtl2RhoR%transFwdNonLinear(XTemp)

    !     dXR = this%OprRadarRef%transAdjMultiply_opr(D, XTemp)

    !     ! Rho To Rho ctl
    !     CALL this%RhoRCtl2RhoR%transAdjMultiply(dXR)

    !     ! Add the increment
    !     dX = dX + dXR
    !   END BLOCK
    ! END IF
    ! GNSS RO refractivity
    IF (D%getObsIdx('refractivity') .NE. 0 &
        .AND. (XTemp%getVarIdx(TRIM('temp')) .NE. 0) &
        .AND. (XTemp%getVarIdx(TRIM('pres')) .NE. 0) & !.OR. (X%getVarIdx(TRIM('pres_ctl')) .NE. 0))  &
        .AND. (XTemp%getVarIdx(TRIM('qvapor')) .NE. 0) & !.OR. (X%getVarIdx(TRIM('qvapor_ctl')) .NE. 0)) &
        ) THEN
      CALL this%GNSSRefrac%transAdjMultiply(D, dX, XTemp)
    END IF

    ! Radar radial wind
    IF (D%getObsIdx('rwnd') .NE. 0 &
        .AND. X%getVarIdx(TRIM('uwnd')) .NE. 0 &
        .AND. X%getVarIdx(TRIM('vwnd')) .NE. 0) THEN
      CALL this%OprRadarVel%transAdjMultiply(D, dX)
    END IF

    ! Radar radial wind strong constraint forward
    IF (this%useUV2W) THEN
      ! w to u, v
      CALL this%UV2W%transAdjMultiply(dX)
      CALL dX%rmVar('wwnd')
    END IF

#ifdef TRACE_PRESSURE
    PRINT *, 'MMdX2 = ', MINVAL(dX%fields(5)%DATA(1, :, :)), MAXVAL(dX%fields(5)%DATA(1, :, :))
#endif

    ! ! Cumulus Adjoint
    ! IF (dX%getVarIdx(TRIM('uwnd')) .NE. 0 &
    !     .AND. dX%getVarIdx(TRIM('vwnd')) .NE. 0 &
    !     .AND. dX%getVarIdx(TRIM('wwnd')) .NE. 0 &
    !     .AND. dX%getVarIdx(TRIM('theta')) .NE. 0 &
    !     .AND. dX%getVarIdx(TRIM('qvapor')) .NE. 0 &
    !     .AND. dX%getVarIdx(TRIM('pres')) .NE. 0 &
    !     ) THEN
    !       BLOCK
    !         TYPE(State_t) :: xTemp
    !       XTemp = X
    !   CALL this%CumInterface%CumForward(XTemp)
    !   CALL this%CumInterface%CumAdjoint(XTemp, dX)

    !   END BLOCK
    ! END IF

    ! IF (dX%getVarIdx('qvapor') /= 0) &
    !   dX%Fields(dX%getVarIdx('qvapor'))%data = dX%Fields(dX%getVarIdx('qvapor'))%data/1000.0D0
    ! IF (dX%getVarIdx('pres') /= 0) &
    !   dX%Fields(dX%getVarIdx('pres'))%data = dX%Fields(dX%getVarIdx('pres'))%data*100.0D0

    ! Satellite
    IF (D%getObsIdx('tbb') .NE. 0) THEN

      BLOCK

        TYPE(State_t) :: gradx

        WHERE (XTemp%fields(XTemp%getVarIdx('qvapor'))%DATA < qv_limit) XTemp%fields(XTemp%getVarIdx('qvapor'))%DATA = qv_limit
        IF (XTemp%getVarIdx('qice') .GT. 0) THEN
          WHERE (XTemp%fields(XTemp%getVarIdx('qcloud'))%DATA < 0.0D0) XTemp%fields(XTemp%getVarIdx('qcloud'))%DATA = 0.0D0
          WHERE (XTemp%fields(XTemp%getVarIdx('qice'))%DATA < 0.0D0) XTemp%fields(XTemp%getVarIdx('qice'))%DATA = 0.0D0
        END IF

        gradx = XTemp%zeroCopy()

        ! CALL RTTOV AD for gradient calculation

        DO i_inst = 1, this%n_insts

          IF (this%turnOn(i_inst) .AND. D%getObsTypeIdx(TRIM(this%platform_name(i_inst))) .NE. 0 &
              .AND. D%getObsTypeIdx(TRIM(this%inst_name(i_inst))) .NE. 0) &
            CALL this%OprRTTOVtlad(i_inst)%rttov_ad_simobs(XTemp, D, gradx)

          ! PRINT *, 'check ad qvapor: ', gradx%sg%mpddInfo_sg%myrank, MAXVAL(gradx%fields(gradx%getVarIdx('qvapor'))%data), MINVAL(gradx%fields(gradx%getVarIdx('qvapor'))%data)
          ! PRINT *, 'check ad temp: ', gradx%sg%mpddInfo_sg%myrank, MAXVAL(gradx%fields(gradx%getVarIdx('temp'))%data), MINVAL(gradx%fields(gradx%getVarIdx('temp'))%data)

          dX = dX + gradx
        END DO
      END BLOCK
    END IF

    IF (.FALSE. .AND. X%getVarIdx('uwnd') .GT. 0 .AND. &
        X%getVarIdx('vwnd') .GT. 0 .AND. &
        X%getVarIdx('temp') .GT. 0 .AND. &
        (X%getVarIdx('rho') .GT. 0 .OR. X%getVarIdx('pres') .GT. 0)) THEN
      IF (X%sg%gLevel .GE. this%GeosBal%geobalBeg .AND. &
          X%sg%gLevel .LE. this%GeosBal%geobalEnd) &
        dX = this%GeosBal%transAdjMultiply_opr(dX, X)
    END IF

    CALL this%Ctl2State%transAdjMultiply(dX, X)

    ! CALL dX%exHalo()
  END FUNCTION transAdjMultiply_opr

END MODULE C2O_m
