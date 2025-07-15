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
MODULE ObsRadarRef_m
  USE ObsBase_m, ONLY: ObsBase_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  !USE NMLRead_m
  USE YAMLRead_m
  USE State_m, ONLY: State_t
  USE RhoRCtl2RhoR_m, ONLY: RhoRCtl2RhoR_t
  USE OprRadarRef_m, ONLY: OprRadarRef_t
  USE parameters_m, ONLY: missing
  IMPLICIT NONE

  TYPE, EXTENDS(ObsBase_t) :: ObsRadarRef_t
    CHARACTER(LEN=20), ALLOCATABLE :: obsVars(:)    ! Obs variable names
    TYPE(RhoRCtl2RhoR_t) :: RhoRCtl2RhoR
    TYPE(OprRadarRef_t) :: OprRadarRef

  CONTAINS
    PROCEDURE :: ObsInitial => radarRefInitial
    PROCEDURE :: ObsIngest => radarRefcIngest
    PROCEDURE :: ObsForward => radarRefcForward
    PROCEDURE :: ObsTangent => radarRefcTangent
    PROCEDURE :: ObsAdjoint => radarRefcADjoint
    PROCEDURE :: ObsQC => radarRefQC

    PROCEDURE, PUBLIC :: GetForwardValue  ! Used by superObs to calculate the increment of Obs
    PROCEDURE, PUBLIC :: ObsPrepareForSg  ! Used for observations which has to be prepared before thinning, like radar and satellite
    PROCEDURE, NOPASS, PUBLIC :: Field2ObsIsExisted
  END TYPE ObsRadarRef_t

CONTAINS

  SUBROUTINE ObsPrepareForSg(this, X)
    CLASS(ObsRadarRef_t) :: this
    TYPE(State_t) :: X
    ! This is an optional holder for data prepare polymorphic subroutine

    this%RhoRCtl2RhoR = RhoRCtl2RhoR_t(this%configfile, X)
    this%OprRadarRef = OprRadarRef_t(this%configfile, X)

  END SUBROUTINE

  FUNCTION GetForwardValue(this, X, varnametmp, vIdx, hIdx, tIdx, iv) RESULT(ref)
    CLASS(ObsRadarRef_t) :: this
    TYPE(State_t) :: X
    INTEGER(i_kind), INTENT(IN) :: vIdx, hIdx, tIdx
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: iv
    REAL(r_kind) :: ref
    CHARACTER(LEN=*) :: varnametmp

    ! PRINT *, 'aaa'
    IF (X%getVarIdx('ref') .NE. 0) THEN
      ref = X%Fields(X%getVarIdx(TRIM(varnametmp)))%DATA(vIdx, hIdx, tIdx)
      ! PRINT *, 'ref', ref
    ELSE IF (X%getVarIdx('rhor_ctl') .NE. 0) THEN
      BLOCK
        REAL(r_kind) :: rhor
        rhor = this%RhoRCtl2RhoR%transElemForward(X%fields(X%getVarIdx('rhor_ctl'))%DATA(vIdx, hIdx, tIdx), X%sg%s1(vIdx, hIdx))
        ref = this%OprRadarRef%transElemForward(rhor)
      END BLOCK
    ELSE IF (X%getVarIdx('rhor') .NE. 0) THEN
      ref = X%fields(X%getVarIdx('rhor'))%DATA(vIdx, hIdx, tIdx)
    END IF
  END FUNCTION

  FUNCTION Field2ObsIsExisted(X, varnametmp)
    TYPE(State_t) :: X
    CHARACTER(LEN=*) :: varnametmp    ! Yuanfu Xie changed the string length to * from 20 as some codes complains about the length
    LOGICAL :: Field2ObsIsExisted

    Field2ObsIsExisted = .FALSE.

    IF ((X%getVarIdx(TRIM('rhor_ctl')) .NE. 0)) Field2ObsIsExisted = .TRUE.
    IF ((X%getVarIdx(TRIM('ref')) .NE. 0)) Field2ObsIsExisted = .TRUE.
    IF ((X%getVarIdx(TRIM('rhor')) .NE. 0)) Field2ObsIsExisted = .TRUE.
  END FUNCTION

  SUBROUTINE radarRefInitial(this, configFile, idxFile)
    IMPLICIT NONE
    CLASS(ObsRadarRef_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER, OPTIONAL :: idxFile

    CHARACTER(LEN=1024) :: inputFileDir, radarFileDir
    CHARACTER(LEN=1024), ALLOCATABLE :: obsFileName(:)
    CHARACTER(LEN=20), ALLOCATABLE :: obsVarList(:)
    INTEGER :: i, istatus

    ! read the obs surface files names
    istatus = yaml_get_var(configFile, 'IO', 'input_dir_Radar', inputFileDir)

    !CALL namelist_read(TRIM(configFile), "obsFileList_Radar", obsFileName)
    !CALL namelist_read(TRIM(configFile), "obsFilePath_Radar", radarFileDir)
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsFileList_Radar', obsFileName)
    ! istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsFilePath_Radar', radarFileDir)

    this%numFiles = UBOUND(obsFileName, 1)
    PRINT *, "amount of radar files: ", this%numFiles

    ALLOCATE (this%fileNames(this%numfiles))
    DO i = 1, this%numFiles
      this%fileNames(i) = TRIM(inputFileDir)//"/"//TRIM(obsFileName(i))
      PRINT *, "name of radar file ", i, ": ", TRIM(this%fileNames(i))
    END DO

    ! set the var names of obs files into this%obsVars(:)
    this%numVars = 1
    ALLOCATE (this%obsVars(this%numVars))
    this%obsVars(1) = "reflectivity"
    ! this%obsVars(1) = "densityOfRain"

    ! to match the var names between obs file and analysis fields.
    ALLOCATE (this%varNames(this%numVars))
    DO i = 1, this%numVars
      IF (TRIM(this%obsVars(i)) .EQ. "reflectivity") this%varNames(i) = "ref"
      IF (TRIM(this%obsVars(i)) .EQ. "densityOfRain") this%varNames(i) = "rhor"
      PRINT *, "varName: ", TRIM(this%varNames(i))
    END DO

    ! Read in an interpolation option:
    !CALL namelist_read(TRIM(configFile), "interpolation",this%interpolation)
    istatus = yaml_get_var(TRIM(configFile), 'obs_thinning', 'interpolation', this%interpolation)

    this%obsType = "RADAR"
    ALLOCATE (this%types(this%numVars))
    DO i = 1, this%numVars
      this%types(i) = this%obsType//"_"//TRIM(this%varNames(i))
    END DO

    ALLOCATE (this%radius(4, this%numVars), this%sizeInc(this%numVars), &
              this%qcThreshold(this%numVars))

    this%radius(1, :) = 5.0D2   ! Temporary holders, (1:2)=5.0D2 is the optimal one right now
    this%radius(2, :) = 5.0D1   ! Temporal influence radius, 5.0D1 is the optimal one right now
    this%radius(3, :) = 5.0D2   ! Topography influence radius, 5.0D2 is the optimal one right now
    this%radius(4, :) = 5.0D1   ! Temporal influence radius, 5.0D1 is the optimal one right now
    this%sizeInc = 10.0D0      ! Temporary holders
    this%configFile = configFile

    this%correlation_threshold = 0.01D0 ! No influence under this threshold value

    ! Threshold values for QC:
    this%qcThreshold(1) = 500.0D0 ! Temporary holders

  END SUBROUTINE radarRefInitial

  SUBROUTINE radarRefcIngest(this, X)
    USE RadarRAW_m, ONLY: RadarRAW_t
    IMPLICIT NONE
    CLASS(ObsRadarRef_t) :: this
    TYPE(State_t) :: X

    ! Local variables:
    LOGICAL :: fileExist(this%numFiles)
    TYPE(RadarRAW_t), ALLOCATABLE :: RadarRAW(:)
    INTEGER(i_kind) :: sizeRef
    INTEGER(i_kind) :: i, j, k

    sizeRef = 0
    ALLOCATE (RadarRAW(this%numFiles))
    DO i = 1, this%numFiles
      ! Yuanfu Xie added a safeguard for missing file 2022-05-31
      INQUIRE (FILE=TRIM(this%fileNames(i)), EXIST=fileExist(i))
      IF (fileExist(i)) THEN
        WRITE (*, 1) TRIM(this%fileNames(i))
1       FORMAT('ObsRadarRef: Radar reflectivity data file: ', A, ' is found and reading...')
      ELSE
        WRITE (*, 3) TRIM(this%fileNames(i))
3       FORMAT('ObsRadarRef: This radar reflectivity data file does not exist but the run continues, CHECK! ', A)
        CYCLE
      END IF
      ! Constructer for RadarRAW
      RadarRAW(i) = RadarRAW_t(this%configFile)
      PRINT *, "Reading: ", TRIM(this%fileNames(i)), ' ...'
      ! Read the radar data
      CALL RadarRAW(i)%ReadRadarData(this%fileNames(i), 'ref')
      sizeRef = sizeRef + RadarRAW(i)%sizeRef
      PRINT *, 'Done reading file: ', i, ' sizeRef: ', sizeRef
    END DO

    ! get the unique elements
    ALLOCATE (this%obsData(sizeRef, this%numVars))
    ALLOCATE (this%obsErrs(sizeRef, this%numVars))
    ALLOCATE (this%olatlon(2, sizeRef))
    ALLOCATE (this%obsTime(sizeRef))
    ALLOCATE (this%obsHght(sizeRef))
    ALLOCATE (this%land(sizeRef))
    ! ALLOCATE (this%StNames(sizeRef))

    BLOCK
      INTEGER(i_kind) :: countIdxStart, countIdxEnd
      this%numObs = sizeRef
      countIdxStart = 1
      DO i = 1, this%numFiles

        IF (.NOT. fileExist(i)) CYCLE  ! Yuanfu Xie added a safeguard for missing file 2022-05-31

        countIdxEnd = countIdxStart + RadarRAW(i)%sizeRef - 1

        IF (TRIM(this%varNames(1)) == 'ref') this%obsData(countIdxStart:countIdxEnd, 1) = RadarRAW(i)%valRef
        ! Yuanfu Xie converted the rho_r to the standard density unit: m3/m3
        ! Note: Z-R: ref (dBZ) = C_1 + C_2 * log(rhor), where rho_r is in g/m^3 = 1.0D-3kg/(1.0D3 kg),
        ! where C_1 = 43.1, and C_2 = 17.5. Thus, rho_r (kg/kg) = 10**[(ref-C_1)/C_2]*1.0D-6
        IF (TRIM(this%varNames(1)) == 'rhor') this%obsData(countIdxStart:countIdxEnd, 1) = &
          10**((RadarRAW(i)%valRef - 43.1) / 17.5) * 1.0D-6    ! Convert from dBZ to kg/kg

        ! Check the radar data:
#ifdef DEBUG
        WRITE (*, 4) MAXVAL(RadarRAW(i)%valRef), MINVAL(RadarRAW(i)%valRef), &
          MAXVAL(this%obsData(countIdxStart:countIdxEnd, 1)), &
          MINVAL(this%obsData(countIdxStart:countIdxEnd, 1)), X%mpddGlob%myrank
4       FORMAT('RefIngest - MAX/MIN RawRef: ', 2D10.2, ' Converted: ', 2D10.2, ' proc/file: ', 2I3)
#endif
        this%obsErrs(countIdxStart:countIdxEnd, 1) = 1.0D0
        this%olatlon(1, countIdxStart:countIdxEnd) = RadarRAW(i)%latRef
        this%olatlon(2, countIdxStart:countIdxEnd) = RadarRAW(i)%lonRef
        this%obsTime(countIdxStart:countIdxEnd) = RadarRAW(i)%timRef
        this%obsHght(countIdxStart:countIdxEnd) = RadarRAW(i)%altRef
        this%land = 0.0D0

        countIdxStart = countIdxEnd + 1
      END DO
    END BLOCK
    DEALLOCATE (RadarRAW)
  END SUBROUTINE

  FUNCTION radarRefcForward(this, X, O) RESULT(Y)
    USE RadarRAW_m, ONLY: RadarRAW_t
    USE ObsSet_m, ONLY: ObsSet_t
    IMPLICIT NONE
    CLASS(ObsRadarRef_t) :: this
    TYPE(State_t) :: X
    TYPE(ObsSet_t), INTENT(IN) :: O ! Telling the operator to map X on to this obs
    TYPE(ObsSet_t) :: Y

    Y = O
    PRINT *, 'radarRefc: This forward operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION radarRefcForward

  FUNCTION radarRefcTangent(this, dX, X) RESULT(Y)
    USE RadarRAW_m, ONLY: RadarRAW_t
    USE ObsSet_m, ONLY: ObsSet_t
    IMPLICIT NONE
    CLASS(ObsRadarRef_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: Y

    PRINT *, 'radarRefc: This tangent operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION radarRefcTangent

  FUNCTION radarRefcAdjoint(this, dY, X) RESULT(dX)
    USE RadarRAW_m, ONLY: RadarRAW_t
    USE ObsSet_m, ONLY: ObsSet_t
    IMPLICIT NONE
    CLASS(ObsRadarRef_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: dY

    PRINT *, 'radarRefc: This adjoint operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION radarRefcAdjoint

  SUBROUTINE radarRefQC(this)
    CLASS(ObsRadarRef_t) :: this
  END SUBROUTINE radarRefQC

END MODULE ObsRadarRef_m
