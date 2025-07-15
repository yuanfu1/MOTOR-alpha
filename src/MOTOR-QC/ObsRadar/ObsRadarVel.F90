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
MODULE ObsRadarVel_m
  USE ObsBase_m, ONLY: ObsBase_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  !USE NMLRead_m
  USE YAMLRead_m
  USE State_m, ONLY: State_t
  USE ObsRadarVelEach_m, ONLY: ObsRadarVelEach_t
  USE RadarRAW_m, ONLY: getStaName

  TYPE, EXTENDS(ObsBase_t) :: ObsRadarVel_t
    CHARACTER(LEN=20), ALLOCATABLE :: obsVars(:)    ! Obs variable names
    TYPE(ObsRadarVelEach_t), ALLOCATABLE :: radars(:)
    INTEGER(i_kind) :: nStations
    CHARACTER(LEN=10) :: obsFileType_Radar

  CONTAINS
    PROCEDURE :: ObsInitial => radarVelInitial
    PROCEDURE :: ObsIngest => radarVelcIngest
    PROCEDURE :: ObsForward => radarVelcForward
    PROCEDURE :: ObsTangent => radarVelcTangent
    PROCEDURE :: ObsAdjoint => radarVelcAdjoint
    PROCEDURE :: ObsQC => radarVelQC

    PROCEDURE, PUBLIC :: ObsPrepareForSg  ! Used for observations which has to be prepared before thinning, like radar and satellite
    PROCEDURE, PUBLIC :: ObsSuper         ! mapping raw observation to super observations
    PROCEDURE, PUBLIC :: ObsThinning      ! Modified ObsSuper by Xie
    FINAL :: destructor                      ! destroy arrays of observations
  END TYPE ObsRadarVel_t

  TYPE idxArray_t
    INTEGER(i_kind), ALLOCATABLE :: idx(:)
  CONTAINS
    FINAL :: destructorOfidxArray_t   ! destroy arrays
  END TYPE

CONTAINS

  IMPURE ELEMENTAL SUBROUTINE destructorOfidxArray_t(this)
    TYPE(idxArray_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%idx)) DEALLOCATE (this%idx)
  END SUBROUTINE

  SUBROUTINE ObsPrepareForSg(this, X)
    CLASS(ObsRadarVel_t) :: this
    TYPE(State_t) :: X
    ! This is an optional holder for data prepare polymorphic subroutine
    INTEGER :: i

    DO i = 1, this%nStations
      CALL this%radars(i)%ObsPrepareForSg(X)
    END DO
  END SUBROUTINE

  SUBROUTINE ObsSuper(this, state, thinObs, mpObs)
    USE State_m, ONLY: State_t
    USE ObsSet_m, ONLY: ObsSet_t
    USE mpObs_m, ONLY: mpObs_t
    IMPLICIT NONE

    CLASS(ObsRadarVel_t) :: this
    TYPE(State_t), INTENT(IN)  :: state
    TYPE(ObsSet_t), INTENT(OUT) :: thinObs
    TYPE(mpObs_t), INTENT(IN)  :: mpObs
    INTEGER :: i

    IF (mpObs%isActiveProc()) THEN
      thinObs = ObsSet_t(this%configFile, mpObs)
      ALLOCATE (thinObs%ObsFields(this%nStations))
      DO i = 1, this%nStations
        PRINT *, "Finish thinning radar radial wind in ", i, ": ", TRIM(this%fileNames(i))
        BLOCK
          TYPE(ObsSet_t) :: thinObsTemp
          CALL this%radars(i)%ObsThinning(state, thinObsTemp, mpObs, .TRUE., .FALSE.)
          thinObs%ObsFields(i) = thinObsTemp%ObsFields(1)
          thinObs%ObsFields(i)%locObs = this%radars(i)%locRadar
        END BLOCK
        PRINT *, "Finish thinning radar radial wind in ", i, ": ", TRIM(this%fileNames(i))
        CALL state%sg%mpddInfo_sg%barrier()
      END DO
    END IF

  END SUBROUTINE ObsSuper

  SUBROUTINE ObsThinning(this, state, thinObs, mpObs, useBkg, vector)
    USE State_m, ONLY: State_t
    USE ObsSet_m, ONLY: ObsSet_t
    USE mpObs_m, ONLY: mpObs_t
    IMPLICIT NONE

    CLASS(ObsRadarVel_t) :: this
    TYPE(State_t), INTENT(IN)  :: state
    TYPE(ObsSet_t), INTENT(INOUT) :: thinObs
    TYPE(mpObs_t), INTENT(IN)  :: mpObs
    LOGICAL, INTENT(IN) :: useBkg, vector
    INTEGER :: i

    PRINT *, "In RADAR-VEL thinning."

    IF (mpObs%isActiveProc()) THEN
      thinObs = ObsSet_t(this%configFile, mpObs)
      ALLOCATE (thinObs%ObsFields(this%nStations))
      DO i = 1, this%nStations
        PRINT *, "Start thinning radar radial wind in station: ", i
        BLOCK
          TYPE(ObsSet_t) :: thinObsTemp
          CALL this%radars(i)%ObsThinning(state, thinObsTemp, mpObs, useBkg, vector)
          thinObs%ObsFields(i) = thinObsTemp%ObsFields(1)
          thinObs%ObsFields(i)%locObs = this%radars(i)%locRadar
        END BLOCK
        ! CALL state%sg%mpddInfo_sg%barrier()
        ! PRINT *, "Finish thinning radar radial wind in ", i, ": ", TRIM(this%fileNames(i))
        ! CALL state%sg%mpddInfo_sg%barrier()
      END DO
    END IF

  END SUBROUTINE ObsThinning

  SUBROUTINE radarVelInitial(this, configFile, idxFile)
    IMPLICIT NONE
    CLASS(ObsRadarVel_t) :: this
    INTEGER, OPTIONAL :: idxFile
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    CHARACTER(LEN=14), ALLOCATABLE :: staNames(:)
    TYPE(idxArray_t), ALLOCATABLE ::idxArray(:)
    CHARACTER(LEN=1024) :: inputFileDir

    INTEGER :: i, istatus, j

    ! read the obs surface files names
    istatus = yaml_get_var(configFile, 'IO', 'input_dir_Radar', inputFileDir)
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsFileList_Radar', this%fileNames)
    istatus = yaml_get_var(TRIM(configFile), 'IO', 'obsFileType_Radar', this%obsFileType_Radar)
    PRINT *, "amount of radar files: ", UBOUND(this%fileNames, 1)

    this%nStations = 0

    IF (this%obsFileType_Radar == 'raw') THEN
      DO i = 1, UBOUND(this%fileNames, 1)
        BLOCK
          CHARACTER(LEN=14) :: staName
          LOGICAL :: isNewSta
          INTEGER(i_kind) :: idxSta
          TYPE(idxArray_t) :: temp

          isNewSta = .TRUE.
          PRINT *, "FileName: ", TRIM(TRIM(inputFileDir)//"/"//this%fileNames(i))
          ! staName = getStaName(TRIM(TRIM(inputFileDir)//"/"//this%fileNames(i)))
          staName = this%fileNames(i) (1:14)
          PRINT *, "staName: ", staName, i, SIZE(staNames)

          IF (ALLOCATED(staNames)) THEN
          DO j = 1, SIZE(staNames)
            IF (TRIM(staNames(j)) == TRIM(staName)) THEN
              isNewSta = .FALSE.
              idxSta = j
              EXIT
            END IF
          END DO
          END IF

          IF (isNewSta) THEN
            this%nStations = this%nStations + 1
            IF (.NOT. ALLOCATED(staNames)) THEN
              ALLOCATE (staNames(1), idxArray(1))
              staNames(1) = staName
              ALLOCATE (idxArray(1)%idx(1))
              idxArray(1)%idx(1) = i
            ELSE
              staNames = (/staNames, staName/)
              ALLOCATE (temp%idx(1))
              temp%idx(1) = i
              BLOCK
                TYPE(idxArray_t), ALLOCATABLE :: tmp(:)
                tmp = idxArray
                DEALLOCATE (idxArray)
                ALLOCATE (idxArray(SIZE(tmp) + 1))
                idxArray(1:SIZE(tmp)) = tmp
                idxArray(SIZE(tmp) + 1) = temp
                DEALLOCATE (tmp)
              END BLOCK
            END IF
          ELSE
            BLOCK
              INTEGER(i_kind), ALLOCATABLE :: tmp(:)
              tmp = idxArray(idxSta)%idx
              DEALLOCATE (idxArray(idxSta)%idx)
              ALLOCATE (idxArray(idxSta)%idx(SIZE(tmp) + 1))
              idxArray(idxSta)%idx(1:SIZE(tmp)) = tmp
              idxArray(idxSta)%idx(SIZE(tmp) + 1) = i
              DEALLOCATE (tmp)
              PRINT *, 'idxArray(i)%idx: ', idxArray(idxSta)%idx
            END BLOCK
          END IF
        END BLOCK
      END DO
      ALLOCATE (this%radars(this%nStations))
      DO i = 1, this%nStations
        CALL this%radars(i)%radarVelInitial(configFile, idxArray(i)%idx)
      END DO
    ELSEIF (this%obsFileType_Radar == 'merge') THEN
      this%nStations = SIZE(this%fileNames)
      ALLOCATE (this%radars(this%nStations))
      DO i = 1, this%nStations
        CALL this%radars(i)%radarVelInitialForMegerFile(configFile, this%fileNames(i))
      END DO
    END IF

    IF (ALLOCATED(staNames)) DEALLOCATE (staNames)

    IF (ALLOCATED(idxArray)) THEN
      DO i = 1, SIZE(idxArray)
        IF (ALLOCATED(idxArray(i)%idx)) DEALLOCATE (idxArray(i)%idx)
      END DO
      ! Warining intel compiler has some stranger errors here.
      IF (SIZE(idxArray) /= 1) DEALLOCATE (idxArray)
    END IF

    PRINT *, 'Done Radar vel initial - 3 !'
  END SUBROUTINE radarVelInitial

  SUBROUTINE radarVelcIngest(this, X)
    CLASS(ObsRadarVel_t) :: this
    TYPE(State_t) :: X

    DO i = 1, this%nStations
      PRINT *, "Ingest radar station: ", i
      CALL this%radars(i)%ObsIngest(X)
    END DO
  END SUBROUTINE radarVelcIngest

  FUNCTION radarVelcForward(this, X, O) RESULT(Y)
    USE ObsSet_m, ONLY: ObsSet_t
    CLASS(ObsRadarVel_t) :: this
    TYPE(State_t) :: X
    TYPE(ObsSet_t), INTENT(IN) :: O ! Telling the operator to map X on to this obs
    TYPE(ObsSet_t) :: Y

    Y = O
    PRINT *, 'radarVelc: This forward operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION radarVelcForward

  FUNCTION radarVelcTangent(this, dX, X) RESULT(Y)
    USE ObsSet_m, ONLY: ObsSet_t
    CLASS(ObsRadarVel_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: Y

    PRINT *, 'radarVelc: This tangent operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION radarVelcTangent

  FUNCTION radarVelcAdjoint(this, dY, X) RESULT(dX)
    USE ObsSet_m, ONLY: ObsSet_t
    CLASS(ObsRadarVel_t) :: this
    TYPE(State_t) :: dX, X
    TYPE(ObsSet_t) :: dY

    PRINT *, 'radarVelc: This adjoint operator has not implemented from ctl2State code, the code stops here'
    STOP
  END FUNCTION radarVelcAdjoint

  SUBROUTINE radarVelQC(this)
    CLASS(ObsRadarVel_t) :: this

  END SUBROUTINE radarVelQC

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(ObsRadarVel_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%radars)) DEALLOCATE (this%radars)
  END SUBROUTINE destructor
END MODULE ObsRadarVel_m
