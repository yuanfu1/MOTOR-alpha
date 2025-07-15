!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/05/30, @GBA-MWF, Shenzhen for adding obs
!     forward, tangent and adjoint operators.
!!--------------------------------------------------------------------------------------------------

!> @brief
MODULE ObsRadarVelEach_m
  USE ObsBase_m, ONLY: ObsBase_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  !USE NMLRead_m
  USE YAMLRead_m
  USE State_m, ONLY: State_t
  USE OprRadarVel_m, ONLY: OprRadarVel_t, calRadarVelFactorElem
  USE parameters_m, ONLY: degree2radian
  IMPLICIT NONE

  TYPE, EXTENDS(ObsBase_t) :: ObsRadarVelEach_t
    CHARACTER(LEN=20), ALLOCATABLE :: obsVars(:)    ! Obs variable names
    TYPE(OprRadarVel_t) :: OprRadarVel
    REAL(r_kind) :: locRadar(3)
    CHARACTER(len=5) :: staName
    CHARACTER(LEN=1024) :: inputFileDir

  CONTAINS
    PROCEDURE :: ObsInitial => radarVelInitialInferface
    PROCEDURE :: ObsIngest => radarVelcIngest
    PROCEDURE :: ObsForward => radarVelcForward
    PROCEDURE :: ObsTangent => radarVelcTangent
    PROCEDURE :: ObsAdjoint => radarVelcAdjoint
    PROCEDURE :: ObsQC => radarVelQC

    PROCEDURE, PUBLIC :: radarVelInitial
    PROCEDURE, PUBLIC :: radarVelInitialForMegerFile
    PROCEDURE, PUBLIC :: GetForwardValue  ! Used by superObs to calculate the increment of Obs
    PROCEDURE, PUBLIC :: ObsPrepaVelorSg  ! Used for observations which has to be prepared before thinning, like radar and satellite
    PROCEDURE, NOPASS, PUBLIC :: Field2ObsIsExisted

    PROCEDURE, PUBLIC :: dumpVADToFile
    PROCEDURE, PUBLIC :: readVADFromFile
  END TYPE ObsRadarVelEach_t

CONTAINS

  SUBROUTINE ObsPrepaVelorSg(this, X)
    CLASS(ObsRadarVelEach_t) :: this
    TYPE(State_t) :: X
    ! This is an optional holder for data prepare polymorphic subroutine

    this%OprRadarVel = OprRadarVel_t(this%configfile, X)
    IF (.NOT. this%Field2ObsIsExisted(X, 'rwnd')) THEN
      PRINT *, 'Error using ObsRadarVelEach_m in ObsPrepaVelorSg. STOP.'
      STOP
    END IF
  END SUBROUTINE

  FUNCTION GetForwardValue(this, X, varnametmp, vIdx, hIdx, tIdx, iv) RESULT(Vel)
    CLASS(ObsRadarVelEach_t) :: this
    TYPE(State_t) :: X
    INTEGER(i_kind), INTENT(IN) :: vIdx, hIdx, tIdx
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: iv
    REAL(r_kind) :: Vel
    CHARACTER(LEN=*) :: varnametmp

    ! PRINT *, 'aaa'
    IF (X%getVarIdx('rwnd') .NE. 0) THEN
      Vel = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(vIdx, hIdx, tIdx)
    ELSE
      BLOCK
        REAL(r_kind) :: uwnd, vwnd, wwnd, locGrid(3)

        uwnd = X%fields(X%getVarIdx('uwnd'))%DATA(vIdx, hIdx, tIdx)
        vwnd = X%fields(X%getVarIdx('vwnd'))%DATA(vIdx, hIdx, tIdx)

        IF (X%getVarIdx('wwnd') .NE. 0) THEN
          wwnd = X%fields(X%getVarIdx('wwnd'))%DATA(vIdx, hIdx, tIdx)
        ELSE
          wwnd = 0.0D0
        END IF

        locGrid = (/X%sg%cell_cntr(1, hIdx), X%sg%cell_cntr(2, hIdx), X%sg%zHght(vIdx, hIdx)/)

        Vel = this%OprRadarVel%transElemForwardForThinning(uwnd, vwnd, wwnd, this%locRadar, locGrid)
      END BLOCK
    END IF
  END FUNCTION

  FUNCTION Field2ObsIsExisted(X, varnametmp)
    TYPE(State_t) :: X
    CHARACTER(LEN=*) :: varnametmp
    LOGICAL :: Field2ObsIsExisted

    Field2ObsIsExisted = .FALSE.

    ! PRINT*, 'In Field2ObsIsExisted of Radar Rwnd'
    IF ((X%getVarIdx(TRIM('uwnd')) .NE. 0) .AND. &
        (X%getVarIdx(TRIM('vwnd')) .NE. 0)) Field2ObsIsExisted = .TRUE.

    IF ((X%getVarIdx(TRIM('rwnd')) .NE. 0)) Field2ObsIsExisted = .TRUE.
  END FUNCTION

  SUBROUTINE radarVelInitialInferface(this, configFile, idxFile)
    IMPLICIT NONE
    CLASS(ObsRadarVelEach_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, OPTIONAL :: idxFile
  END SUBROUTINE

  SUBROUTINE radarVelInitialForMegerFile(this, configFile, mergeFileName)
    IMPLICIT NONE
    CLASS(ObsRadarVelEach_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    CHARACTER(LEN=1024) :: inputFileDir, radarFileDir
    CHARACTER(LEN=1024), ALLOCATABLE :: obsFileName(:)
    CHARACTER(LEN=20), ALLOCATABLE :: obsVarList(:)
    CHARACTER(LEN=*) :: mergeFileName
    INTEGER :: i, istatus

    ! read the obs surface files names
    istatus = yaml_get_var(configFile, 'IO', 'input_dir_Radar', inputFileDir)
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsFileList_Radar', obsFileName)

    ! Got the inputFileDir from the yaml file
    this%inputFileDir = inputFileDir
    this%obsType = mergeFileName(1:9)
    PRINT *, 'this%obsType: ', this%obsType

    ! set the var names of obs files into this%obsVars(:)
    this%numVars = 1
    ALLOCATE (this%obsVars(this%numVars))
    this%obsVars(1) = "Velocity"

    ! to match the var names between obs file and analysis fields.
    ALLOCATE (this%varNames(this%numVars))
    DO i = 1, this%numVars
      this%varNames(i) = "rwnd"
    END DO

    ! Read in an interpolation option:
    ! CALL namelist_read(TRIM(configFile), "interpolation", this%interpolation)
    istatus = yaml_get_var(TRIM(configFile), 'obs_thinning', 'interpolation', this%interpolation)

    ALLOCATE (this%radius(4, this%numVars), this%sizeInc(this%numVars), &
              this%qcThreshold(this%numVars))

    this%radius(1, :) = 5.0D2   ! Temporary holders, (1:2)=5.0D2 is the optimal one right now
    this%radius(2, :) = 5.0D1   ! Temporal influence radius, 5.0D1 is the optimal one right now
    this%radius(3, :) = 5.0D2   ! Topography influence radius, 5.0D2 is the optimal one right now
    this%radius(4, :) = 5.0D1   ! Temporal influence radius, 5.0D1 is the optimal one right now
    this%sizeInc = 10.0D0      ! Temporary holders
    this%configFile = configFile

    this%correlation_threshold = 0.6D0 ! No influence under this threshold value

    ! Threshold values for QC:
    this%qcThreshold(1) = 5.0D0 ! Temporary holders

  END SUBROUTINE

  SUBROUTINE radarVelInitial(this, configFile, idxFile)
    IMPLICIT NONE
    CLASS(ObsRadarVelEach_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, ALLOCATABLE, OPTIONAL :: idxFile(:)

    CHARACTER(LEN=1024) :: inputFileDir, radarFileDir
    CHARACTER(LEN=1024), ALLOCATABLE :: obsFileName(:)
    CHARACTER(LEN=20), ALLOCATABLE :: obsVarList(:)
    INTEGER :: i, istatus

    IF ((.NOT. PRESENT(idxFile)) .OR. SIZE(idxFile) .EQ. 0) THEN
      PRINT *, 'Error using ObsRadarVelEach_m in radarVelInitial. STOP.'
      STOP
    END IF

    ! read the obs surface files names
    istatus = yaml_get_var(configFile, 'IO', 'input_dir_Radar', inputFileDir)
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsFileList_Radar', obsFileName)

    ! Got the inputFileDir from the yaml file
    this%inputFileDir = inputFileDir

    this%numFiles = SIZE(idxFile)

    ALLOCATE (this%fileNames(SIZE(idxFile)))
    DO i = 1, SIZE(idxFile)
      this%fileNames(i) = TRIM(inputFileDir)//"/"//TRIM(obsFileName(idxFile(i)))
      PRINT *, 'this%fileNames(i): ', TRIM(this%fileNames(i))
      ! PRINT *, TRIM(obsFileName(idxFile(i)))
      this%obsType = obsFileName(idxFile(i))
      this%obsType = "RAD_"//this%obsType(10:14)
      PRINT *, 'this%obsType: ', this%obsType
    END DO

    ! set the var names of obs files into this%obsVars(:)
    this%numVars = 1
    ALLOCATE (this%obsVars(this%numVars))
    this%obsVars(1) = "Velocity"

    ! to match the var names between obs file and analysis fields.
    ALLOCATE (this%varNames(this%numVars))
    DO i = 1, this%numVars
      this%varNames(i) = "rwnd"
    END DO

    ! Read in an interpolation option:
    ! CALL namelist_read(TRIM(configFile), "interpolation", this%interpolation)
    istatus = yaml_get_var(TRIM(configFile), 'obs_thinning', 'interpolation', this%interpolation)

    ALLOCATE (this%radius(4, this%numVars), this%sizeInc(this%numVars), &
              this%qcThreshold(this%numVars))

    this%radius(1, :) = 5.0D2   ! Temporary holders, (1:2)=5.0D2 is the optimal one right now
    this%radius(2, :) = 5.0D1   ! Temporal influence radius, 5.0D1 is the optimal one right now
    this%radius(3, :) = 5.0D2   ! Topography influence radius, 5.0D2 is the optimal one right now
    this%radius(4, :) = 5.0D1   ! Temporal influence radius, 5.0D1 is the optimal one right now
    this%sizeInc = 10.0D0      ! Temporary holders
    this%configFile = configFile

    this%correlation_threshold = 0.6D0 ! No influence under this threshold value

    ! Threshold values for QC:
    this%qcThreshold(1) = 5.0D0 ! Temporary holders

  END SUBROUTINE radarVelInitial

  SUBROUTINE dumpVADToFile(this)
    USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup
    USE parameters_m, ONLY: degree2radian
    CLASS(ObsRadarVelEach_t) :: this

    CHARACTER(len=1024) :: fileName

    TYPE(NcDataset)   :: nc
    TYPE(NcDimension) :: dim1
    TYPE(NcVariable)  :: var, scalar
    TYPE(NcGroup)     :: grp

    INTEGER(i_kind) :: recNum
    fileName = TRIM(this%inputFileDir)//"/"//TRIM(this%obsType)//'_VAD.nc'

    nc = NcDataset(filename, "w")

    recNum = this%numObs
    dim1 = nc%setDimension("recNum", recNum)

    var = nc%setVariable("latVel", "f32", (/dim1/)); CALL var%setData(this%olatlon(1, :) / degree2radian)
    var = nc%setVariable("lonVel", "f32", (/dim1/)); CALL var%setData(this%olatlon(2, :) / degree2radian)
    var = nc%setVariable("timVel", "f32", (/dim1/)); CALL var%setData(this%obsTime)
    var = nc%setVariable("altVel", "f32", (/dim1/)); CALL var%setData(this%obsHght)
    var = nc%setVariable("valVel", "f32", (/dim1/)); CALL var%setData(this%obsData(:, 1))
    var = nc%setVariable("errVel", "f32", (/dim1/)); CALL var%setData(this%obsErrs(:, 1))

    CALL nc%setAttribute("types", this%types(1))
    CALL nc%setAttribute("obsType", this%obsType)
    CALL nc%setAttribute("SiteLat", this%locRadar(1) / degree2radian)
    CALL nc%setAttribute("SiteLon", this%locRadar(2) / degree2radian)
    CALL nc%setAttribute("SiteAlt", this%locRadar(3))
    CALL nc%setAttribute("MaxLat", MAXVAL(this%olatlon(1, :)) / degree2radian)
    CALL nc%setAttribute("MinLat", MINVAL(this%olatlon(1, :)) / degree2radian)
    CALL nc%setAttribute("MaxLon", MAXVAL(this%olatlon(2, :)) / degree2radian)
    CALL nc%setAttribute("MinLon", MINVAL(this%olatlon(2, :)) / degree2radian)

    PRINT *, 'Output radar files.'

    ! close the dataset
    CALL nc%CLOSE()

  END SUBROUTINE

  SUBROUTINE readVADFromFile(this, X)
    USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable, NcGroup
    USE parameters_m, ONLY: degree2radian
    CLASS(ObsRadarVelEach_t) :: this

    TYPE(State_t) :: X
    TYPE(NcDataset)   :: nc
    TYPE(NcDimension) :: dim1
    TYPE(NcVariable)  :: var, scalar
    TYPE(NcGroup)     :: grp

    CHARACTER(len=1024) :: fileName
    INTEGER(i_kind) :: recNum
    REAL(r_kind), ALLOCATABLE :: tmpData(:)

    fileName = TRIM(this%inputFileDir)//"/"//TRIM(this%obsType)//'_VAD.nc'
    nc = NcDataset(filename, "r")

    dim1 = nc%getDimension("recNum"); recNum = dim1%getLength()

    ALLOCATE (this%obsData(recNum, this%numVars))
    ALLOCATE (this%obsErrs(recNum, this%numVars))
    ALLOCATE (this%olatlon(2, recNum))
    ALLOCATE (this%obsTime(recNum))
    ALLOCATE (this%obsHght(recNum))
    ALLOCATE (this%land(recNum))
    ! ALLOCATE (this%StNames(sizeVel))
    ALLOCATE (this%types(this%numVars))
    this%land = 0

    var = nc%getVariable("latVel"); CALL var%getData(tmpData); this%olatlon(1, :) = tmpData * degree2radian
    var = nc%getVariable("lonVel"); CALL var%getData(tmpData); this%olatlon(2, :) = tmpData * degree2radian
    var = nc%getVariable("timVel"); CALL var%getData(this%obsTime)
    var = nc%getVariable("altVel"); CALL var%getData(this%obsHght)
    var = nc%getVariable("valVel"); CALL var%getData(tmpData); this%obsData(:, 1) = tmpData
    var = nc%getVariable("errVel"); CALL var%getData(tmpData); this%obsErrs(:, 1) = tmpData

    CALL nc%getAttribute("types", this%types(1))
    CALL nc%getAttribute("obsType", this%obsType)
    CALL nc%getAttribute("SiteLat", this%locRadar(1))
    CALL nc%getAttribute("SiteLon", this%locRadar(2))
    CALL nc%getAttribute("SiteAlt", this%locRadar(3))

    this%locRadar(1:2) = this%locRadar(1:2) * degree2radian
    ! CALL nc%getAttribute("MaxLat", MAXVAL(this%olatlon(1, :))/ degree2radian)
    ! CALL nc%getAttribute("MinLat", MINVAL(this%olatlon(1, :))/ degree2radian)
    ! CALL nc%getAttribute("MaxLon", MAXVAL(this%olatlon(2, :))/ degree2radian)
    ! CALL nc%getAttribute("MinLon", MINVAL(this%olatlon(2, :))/ degree2radian)

    BLOCK
      USE geoTools_m, ONLY: GeoBox_t
      TYPE(GeoBox_t) :: geoBox
      LOGICAL, ALLOCATABLE :: idxMask(:)
      INTEGER(i_kind), ALLOCATABLE :: idxSelected(:)
      INTEGER(i_kind) :: i

      GeoBox = GeoBox_t(MAXVAL(X%sg%cell_cntr(1, :)), MINVAL(X%sg%cell_cntr(1, :)), &
                        MAXVAL(X%sg%cell_cntr(2, :)), MINVAL(X%sg%cell_cntr(2, :)))
      ALLOCATE (idxMask(recNum)); idxMask = .FALSE.

      idxMask = geoBox%inBoxElem(this%olatlon(1, :), this%olatlon(2, :), GeoBox%maxLat, GeoBox%minLat, GeoBox%maxLon, GeoBox%minLon)

      FORALL (i=1:recNum, this%obsTime(i) < X%sg%tt(X%sg%tSlots) - 900)
        idxMask(i) = .FALSE.
      END FORALL

      ! DO i = 1, recNum
      !   IF (this%obsHght(i) - this%locRadar(3) > 1700) PRINT*, this%obsHght(i)
      ! END DO
      ! STOP

      FORALL (i=1:recNum, this%obsHght(i) - this%locRadar(3) < 1000)
        idxMask(i) = .FALSE.
      END FORALL

      FORALL (i=1:recNum, mod(i, 2) == 1)
        idxMask(i) = .FALSE.
      END FORALL

      idxSelected = PACK((/(i, i=1, recNum)/), idxMask) !将所需的格点索引，按次序排成一列
      PRINT *, 'SIZE(idxSelected)', SIZE(idxSelected)

      ! Setect the inbox data of this%olatlon, put it into new arrays
      this%obsData = this%obsData(idxSelected, :)
      this%obsErrs = this%obsErrs(idxSelected, :)
      this%olatlon = this%olatlon(:, idxSelected)
      this%obsTime = this%obsTime(idxSelected)
      this%obsHght = this%obsHght(idxSelected)

      this%numObs = SIZE(idxSelected)

      DEALLOCATE (idxMask)
      DEALLOCATE (idxSelected)

    END BLOCK

    IF (ALLOCATED(tmpData)) DEALLOCATE (tmpData)

    PRINT *, 'In reading VAD from file.', recNum
    CALL nc%CLOSE()

  END SUBROUTINE

  SUBROUTINE radarVelcIngest(this, X)
    USE RadarRAW_m, ONLY: RadarRAW_t
    IMPLICIT NONE
    CLASS(ObsRadarVelEach_t) :: this
    TYPE(State_t) :: X
    CHARACTER(LEN=10) :: obsFileType_Radar

    TYPE(RadarRAW_t), ALLOCATABLE :: RadarRAW(:)
    INTEGER(i_kind) :: sizeVel
    INTEGER(i_kind) :: i, j, k

    ! Whether exist the pre-processed file ----------------------
    CHARACTER(len=1024) :: fileName
    LOGICAL :: file_is_exists

    fileName = TRIM(this%inputFileDir)//"/"//TRIM(this%obsType)//'_VAD.nc'

    INQUIRE (file=TRIM(fileName), exist=file_is_exists)
    PRINT *, 'file_is_exists: ', file_is_exists, TRIM(fileName)

    IF (file_is_exists) THEN
      PRINT *, 'Get pre processed VAD data from file: ', TRIM(fileName)
      CALL this%readVADFromFile(X)
      RETURN
    END IF
    ! -----------------------------------------------------------

    sizeVel = 0
    ALLOCATE (RadarRAW(this%numFiles))
    ! i = 1
    DO i = 1, this%numFiles
      ! Constructer for RadarRAW
      RadarRAW(i) = RadarRAW_t(this%configFile)

      CALL X%sg%mpddInfo_sg%barrier
      PRINT *, "Reading: ", TRIM(this%fileNames(i)), ' ...'
      ! Read the radar data
      CALL RadarRAW(i)%ReadRadarData(this%fileNames(i), 'vel')
      sizeVel = sizeVel + RadarRAW(i)%sizeVel

      CALL X%sg%mpddInfo_sg%barrier
      PRINT *, 'Done reading, sizeVel: ', sizeVel
    END DO

    ! get the unique elements
    ALLOCATE (this%obsData(sizeVel, this%numVars))
    ALLOCATE (this%obsErrs(sizeVel, this%numVars))
    ALLOCATE (this%olatlon(2, sizeVel))
    ALLOCATE (this%obsTime(sizeVel))
    ALLOCATE (this%obsHght(sizeVel))
    ALLOCATE (this%land(sizeVel))
    ! ALLOCATE (this%StNames(sizeVel))
    ALLOCATE (this%types(this%numVars))
    this%obsType = 'RAD_'//TRIM(RadarRAW(1)%Station)
    this%locRadar = RadarRAW(1)%locRadar

    ! this%obsType = "RADAR"
    DO j = 1, this%numVars
      this%types(j) = this%obsType//"_"//TRIM(this%varNames(j))
    END DO

    BLOCK
      INTEGER(i_kind) :: countIdxStart, countIdxEnd
      this%numObs = sizeVel
      countIdxStart = 1
      DO i = 1, this%numFiles
        countIdxEnd = countIdxStart + RadarRAW(i)%sizeVel - 1

        ! this%obsData(countIdxStart:countIdxEnd, 1) = RadarRAW(i)%valVel
        this%obsData(countIdxStart:countIdxEnd, 1) = RadarRAW(i)%valVel    ! Convert from dBZ to mm6/m3
        this%obsErrs(countIdxStart:countIdxEnd, 1) = 1.0D0
        this%olatlon(1, countIdxStart:countIdxEnd) = RadarRAW(i)%latVel
        this%olatlon(2, countIdxStart:countIdxEnd) = RadarRAW(i)%lonVel
        this%obsTime(countIdxStart:countIdxEnd) = RadarRAW(i)%timVel
        this%obsHght(countIdxStart:countIdxEnd) = RadarRAW(i)%altVel
        this%land = 0.0D0

        ! ------------------- For unitest
        ! DO j = countIdxStart, countIdxEnd
        !   BLOCK
        !     REAL(r_kind) :: locGrid(3), factors(3), Velo(3), xyt(3,1)

        !     locGrid = (/this%olatlon(1, j), this%olatlon(2, j), this%obsHght(j)/)

        !     factors = calRadarVelFactorElem(this%locRadar, locGrid)
        !     Velo = (/3, 0, 0/)
        !     xyt(1,1) = (this%olatlon(1, j) - 18*degree2radian)/(9*degree2radian);
        !     xyt(2,1) = (this%olatlon(2, j) - 107*degree2radian)/(12*degree2radian);
        !     xyt(3,1) = 0.5D0;

        !     CALL analytic(1, xyt, Velo(1))

        !     xyt(2,1) = (12*degree2radian-(this%olatlon(2, j) - 107*degree2radian))/(12*degree2radian);

        !     CALL analytic(1, xyt, Velo(2))

        !     this%obsData(j, 1) = sum(Velo*factors)
        !   END BLOCK
        ! END DO
        ! ------------------- End of unitest

        countIdxStart = countIdxEnd + 1
      END DO
    END BLOCK

    DEALLOCATE (RadarRAW)
  END SUBROUTINE

  FUNCTION radarVelcForward(this, X, O) RESULT(Y)
    USE RadarRAW_m, ONLY: RadarRAW_t
    USE ObsSet_m, ONLY: ObsSet_t
    IMPLICIT NONE
    CLASS(ObsRadarVelEach_t) :: this
    TYPE(State_t) :: X
    TYPE(ObsSet_t), INTENT(IN) :: O ! Telling the operator to map X on to this obs
    TYPE(ObsSet_t) :: Y

    Y = O
    PRINT *, 'radarVelcThis forward operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION radarVelcForward

  FUNCTION radarVelcTangent(this, dX, X) RESULT(Y)
    USE RadarRAW_m, ONLY: RadarRAW_t
    USE ObsSet_m, ONLY: ObsSet_t
    IMPLICIT NONE
    CLASS(ObsRadarVelEach_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: Y

    PRINT *, 'radarVelcThis tangent operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION radarVelcTangent

  FUNCTION radarVelcAdjoint(this, dY, X) RESULT(dX)
    USE RadarRAW_m, ONLY: RadarRAW_t
    USE ObsSet_m, ONLY: ObsSet_t
    IMPLICIT NONE
    CLASS(ObsRadarVelEach_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: dY

    PRINT *, 'radarVelcThis adjoint operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION radarVelcAdjoint

  SUBROUTINE analytic(npts, xyt, func)

    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN) :: npts
    REAL(r_kind), INTENT(IN) :: xyt(3, npts)
    REAL(r_kind), INTENT(OUT) :: func(npts)

    ! Local variables:
    INTEGER(i_kind) :: i
    REAL(r_kind) :: speed, pi, st, rg, fr, x, y, r, r2

    !speed = 10.0  ! About 74mph  if st=-14 rg=28 fr=20
    !speed = 7.6   ! About 90km/h if st=-14 rg=28 fr=20
    speed = 4.0   ! About 90km/h

    pi = 4.0 * ATAN(1.0)

    st = -3.5
    rg = 7.0
    fr = 5.0

    ! For all observation locations, replace the obs values with the analytic:
    DO i = 1, npts

      x = st + rg * xyt(1, i)
      y = st + rg * xyt(2, i)
      r = ABS(x - y + fr + speed * xyt(3, i))

      r2 = SQRT((x + st)**2 + (y - st)**2)
      r = r - 0.0 - (0.0 + 0.8 * r2) * SQRT(2.0)

      IF (ABS(r) .GT. 1.0E-8) THEN
        func(i) = 0.8 * (1.0 + 2.0 * (1.0 + TANH(r)) * (1.0 + 1.0 / r * SIN(r)**2))
      ELSE
        func(i) = 0.8
      END IF

    END DO

  END SUBROUTINE analytic

  SUBROUTINE radarVelQC(this)
    CLASS(ObsRadarVelEach_t) :: this
  END SUBROUTINE radarVelQC

END MODULE ObsRadarVelEach_m
