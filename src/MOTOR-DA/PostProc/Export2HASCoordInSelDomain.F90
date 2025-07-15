!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2022/10/30, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module export the model state from the coordinate from MOTOR to coordinate of Height Above Surface (HAS)
MODULE Export2HASCoordInSelDomain_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE YAMLRead_m
  USE State_m, ONLY: State_t
  USE parameters_m, ONLY: degree2radian
  USE NCOutput_m, ONLY: NCOutput_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE RectangleGridUtility_m, ONLY: RectangleGridUtility_t
  USE Interp1D_m, ONLY: interp1d
  IMPLICIT NONE

  TYPE Export2HASCoordInSelDomain_t
    CHARACTER(len=1024) :: configFile
    REAL(r_kind), ALLOCATABLE :: latEx(:)
    REAL(r_kind), ALLOCATABLE :: lonEx(:)
    REAL(r_kind), ALLOCATABLE :: llEx(:, :)
    INTEGER(i_kind) :: NX, NY, vLevelOutput
    INTEGER(i_kind), ALLOCATABLE :: idxtgt(:, :)
    REAL(r_kind), ALLOCATABLE :: coetgt(:, :)
    REAL(r_kind), ALLOCATABLE :: HgtAbvMSL(:, :)
    REAL(r_kind), ALLOCATABLE :: HgtAbv(:)

  CONTAINS
    FINAL :: destructor
    PROCEDURE, PUBLIC :: Export_State_AVST
    PROCEDURE, PUBLIC :: Export_State_AVLNST
    PROCEDURE, PUBLIC :: Export_State_SVST
    PROCEDURE, PUBLIC :: Export_State_SVLNST
  END TYPE

  INTERFACE Export2HASCoordInSelDomain_t
    PROCEDURE :: constructor
  END INTERFACE

CONTAINS

  FUNCTION constructor(configFile, sg, X) RESULT(this)
    IMPLICIT NONE
    TYPE(Export2HASCoordInSelDomain_t) :: this
    CHARACTER(len=1024), INTENT(IN) :: configFile
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    TYPE(State_t) :: X

    ! INTEGER(i_kind) ::
    REAL(r_kind) ::  lat1, lat2, lon1, lon2
    REAL(r_kind), ALLOCATABLE :: ll(:)
    REAL(r_kind), ALLOCATABLE :: res(:)
    INTEGER :: ifile
    INTEGER(i_kind) :: i, j, k
    TYPE(RectangleGridUtility_t) :: regular
    REAL(r_kind), ALLOCATABLE :: cell_cntr(:, :)

    this%configFile = configFile
    ifile = yaml_get_var(configFile, 'PostProc', 'ExportDomain', ll)
    ifile = yaml_get_var(configFile, 'PostProc', 'ExportRes', res)
    lat1 = ll(1); lat2 = ll(2); lon1 = ll(3); lon2 = ll(4)

    PRINT *, "Domain is: ", ll
    PRINT *, "Resolution is: ", res

    this%NY = NINT((lat2 - lat1) / res(1)) + 1
    this%NX = NINT((lon2 - lon1) / res(2)) + 1

    ALLOCATE (this%latEx(this%NY))
    ALLOCATE (this%lonEx(this%NX))
    ALLOCATE (this%llEx(2, this%NX * this%NY))

    FORALL (j=1:this%NY) this%latEx(j) = (j - 1) * res(1) + lat1
    FORALL (i=1:this%NX) this%lonEx(i) = (i - 1) * res(2) + lon1

    DO i = 1, this%NX
      DO j = 1, this%NY
        ! FORALL (i=1:this%NX)
        !   FORALL (j=1:this%NY)
        this%llEx(:, i + (j - 1) * this%NX) = (/this%latEx(j), this%lonEx(i)/) * degree2radian
        !   END FORALL
        ! END FORALL
      END DO
    END DO

    IF (sg%isBaseProc()) THEN
      ALLOCATE (cell_cntr(2, sg%num_icell_global))
      ALLOCATE (this%idxtgt(4, this%NX * this%NY), this%coetgt(4, this%NX * this%NY))
    END IF
    CALL sg%aggrGridReal2D(sg%cell_cntr, cell_cntr, [2, sg%num_icell_global])
    IF (sg%isBaseProc()) CALL regular%RectangleHorizontalIntp(TRANSPOSE(cell_cntr), sg%num_icell_global, &
                                                              TRANSPOSE(this%llEx), this%NX * this%NY, &
                                                              this%idxtgt, this%coetgt, sg%mpddInfo_sg%myrank)

    IF (ALLOCATED(cell_cntr)) DEALLOCATE (cell_cntr)

    ! PRINT *, 'lonEx is: ', this%lonEx
    ! PRINT *, 'latEx is: ', this%latEx

    BLOCK
      CHARACTER(len=20) :: Coord
      REAL(r_kind), ALLOCATABLE :: PresLevel(:)

      ifile = yaml_get_var(configFile, 'PostProc', 'Coordinate', Coord)
      PRINT *, 'Coorinate is: ', Coord

      IF (TRIM(Coord) == "HAS") THEN
        ! Construct the coordinates of height above surface
        ifile = yaml_get_var(configFile, 'PostProc', 'HeightLevelAboveSurface', this%HgtAbv)
        this%vLevelOutput = SIZE(this%HgtAbv)
        PRINT *, "HeightLevelAboveSurface: ", this%HgtAbv

        ALLOCATE (this%HgtAbvMSL(this%vLevelOutput, sg%num_cell))
        FORALL (i=1:sg%num_cell)
          this%HgtAbvMSL(:, i) = this%HgtAbv + sg%topo(i)
        END FORALL
      ELSE IF (TRIM(Coord) == "MSL") THEN
        ! Construct the coordinates of height above surface
        ! ifile = yaml_get_var(configFile, 'PostProc', 'HeightLevelAboveSurface', this%HgtAbv)
        ALLOCATE (this%HgtAbv(INT(sg%ztop / 50.0 + 1)))
        DO i = 1, INT(sg%ztop / 50.0 + 1)
          this%HgtAbv(i) = (i - 1) * 50.0
        END DO

        this%vLevelOutput = SIZE(this%HgtAbv)
        PRINT *, "HeightLevelAboveSurface: ", this%HgtAbv

        ALLOCATE (this%HgtAbvMSL(this%vLevelOutput, sg%num_cell))
        FORALL (i=1:sg%num_cell)
          this%HgtAbvMSL(:, i) = this%HgtAbv
        END FORALL

      ELSE IF (TRIM(Coord) == "PRES") THEN
        ! Construct the coordinates of pressure surface
        ifile = yaml_get_var(configFile, 'PostProc', 'PressureLevel', PresLevel)

        this%vLevelOutput = SIZE(PresLevel)
        PRINT *, "PressureLevel: ", PresLevel

        ALLOCATE (this%HgtAbvMSL(this%vLevelOutput, sg%num_cell))

        DO j = 1, sg%num_cell
          CALL interp1d(-1.0D0 * X%fields(X%getVarIdx('pres'))%DATA(:, j, sg%tSlots), sg%zHght(:, j), &
                        -1.0D0 * PresLevel, this%HgtAbvMSL(:, j))
        END DO
        PRINT *, 'this%HgtAbvMSL(:, 10)', this%HgtAbvMSL(:, 10)

        this%HgtAbv = PresLevel
      END IF
    END BLOCK

  END

!> @brief
!! Output state to nc file for single time slot (all variables).
  SUBROUTINE Export_State_AVST(this, X, pathName, FNPrefix, tIdx, hasGLabel, hasMOTOR)
    USE AdvanceTime_m

    CLASS(Export2HASCoordInSelDomain_t) :: this
    TYPE(State_t), INTENT(IN) :: X             !< State.
    CHARACTER(*), INTENT(IN) :: pathName       !< Path for output this state.
    CHARACTER(*), INTENT(IN) :: FNPrefix       !< Prefix of filename.
    INTEGER(i_kind), INTENT(IN) :: tIdx        !< Index of time slot ouptut.

    REAL(r_kind), ALLOCATABLE :: valueGlob(:, :)
    CHARACTER(len=1024) :: fileName
    INTEGER(i_kind) :: gTime(6), uTime_sec, i, j, k
    TYPE(NCOutput_t) :: NCOutput
    LOGICAL, OPTIONAL :: hasGLabel, hasMOTOR
    CHARACTER(len=20) :: timeStr
    REAL(r_kind), ALLOCATABLE :: valueSel(:, :)
    REAL(r_kind), ALLOCATABLE :: valueHAS_MP(:, :)

    uTime_sec = X%sg%tt(tIdx)
    CALL Time_Unix_to_GMT(uTime_sec, gTime)
    ! PRINT *, 'gTime for this time slot is: ', gTime

1001 FORMAT(I4.4, I2.2, I2.2, A, I2.2, I2.2, A)
    WRITE (timeStr, 1001) gTime(1), gTime(2), gTime(3), '_', gTime(4), gTime(5), '_'

    fileName = TRIM(FNPrefix)

    IF (PRESENT(hasGLabel)) THEN
      IF (hasGLabel) THEN
        BLOCK
          CHARACTER(LEN=2) :: iGrid
          WRITE (iGrid, "(I2.2)") X%sg%gLevel
          fileName = TRIM(fileName)//'_G'//iGrid
        END BLOCK
      END IF
    END IF

    IF (PRESENT(hasMOTOR)) THEN
      IF (hasMOTOR) THEN
        fileName = 'MOTOR-'//TRIM(fileName)
      END IF
    END IF

    fileName = TRIM(pathName)//'/'//TRIM(timeStr)//TRIM(fileName)//'.nc'

    ASSOCIATE (sg => X%sg)
      IF (sg%isBaseProc()) THEN
        NCOutput = NCOutput_t(fileName, this%lonEx, this%latEx, this%HgtAbv, sg%tt(tIdx:tIdx), &
                              this%NX, this%NY, this%vLevelOutput, 1)
      END IF

      IF (sg%isBaseProc()) THEN
        ALLOCATE (valueGlob(this%vLevelOutput, sg%num_icell_global))
        ALLOCATE (valueSel(this%vLevelOutput, this%NX * this%NY))
      END IF

      ALLOCATE (valueHAS_MP(this%vLevelOutput, sg%num_cell))

      DO i = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
        DO j = 1, sg%num_cell
          CALL interp1d(sg%zHght(:, j), X%fields(i)%DATA(:, j, tIdx), &
                        this%HgtAbvMSL(:, j), valueHAS_MP(:, j))

          DO k = 1, this%vLevelOutput
            IF (this%HgtAbvMSL(k, j) .LT. sg%topo(j)) THEN
              valueHAS_MP(k, j) = 0.0D0
            END IF
          END DO
        END DO

        CALL sg%aggrGridReal2D(valueHAS_MP, &
                               valueGlob, [this%vLevelOutput, sg%num_icell_global])

        IF (sg%isBaseProc()) THEN
          FORALL (i=1:this%NX * this%NY)
            valueSel(:, i) = valueGlob(:, this%idxtgt(1, i)) * this%coetgt(1, i) + &
                             valueGlob(:, this%idxtgt(2, i)) * this%coetgt(2, i) + &
                             valueGlob(:, this%idxtgt(3, i)) * this%coetgt(3, i) + &
                             valueGlob(:, this%idxtgt(4, i)) * this%coetgt(4, i)
          END FORALL
          CALL NCOutput%addVar(X%fields(i)%Get_Name(), valueSel, &
                               this%vLevelOutput, this%NX, this%NY, 1)
        END IF
      END DO

      ! Export height
      CALL sg%aggrGridReal2D(this%HgtAbvMSL, &
                             valueGlob, [this%vLevelOutput, sg%num_icell_global])

      IF (sg%isBaseProc()) THEN
        FORALL (i=1:this%NX * this%NY)
          valueSel(:, i) = valueGlob(:, this%idxtgt(1, i)) * this%coetgt(1, i) + &
                           valueGlob(:, this%idxtgt(2, i)) * this%coetgt(2, i) + &
                           valueGlob(:, this%idxtgt(3, i)) * this%coetgt(3, i) + &
                           valueGlob(:, this%idxtgt(4, i)) * this%coetgt(4, i)
        END FORALL
        CALL NCOutput%addVar('height', valueSel, &
                             this%vLevelOutput, this%NX, this%NY, 1)
      END IF

      IF (ALLOCATED(valueGlob)) DEALLOCATE (valueGlob)
      IF (ALLOCATED(valueSel)) DEALLOCATE (valueSel)
      IF (ALLOCATED(valueHAS_MP)) DEALLOCATE (valueHAS_MP)

      IF (sg%isBaseProc()) CALL NCOutput%CLOSE()
    END ASSOCIATE

  END SUBROUTINE

!> @brief
!! Output state to nc file for single time slot (single variable).
  SUBROUTINE Export_State_SVST(this, X, pathName, FNPrefix, varName, tIdx, hasGLabel, hasMOTOR)
    USE AdvanceTime_m

    CLASS(Export2HASCoordInSelDomain_t) :: this
    TYPE(State_t), INTENT(IN) :: X             !< State.
    CHARACTER(*), INTENT(IN) :: pathName       !< Path for output this state.
    CHARACTER(*), INTENT(IN) :: FNPrefix       !< Prefix of filename.
    INTEGER(i_kind), INTENT(IN) :: tIdx        !< Index of time slot ouptut.

    REAL(r_kind), ALLOCATABLE :: valueGlob(:, :)
    INTEGER :: i, j, k, l
    CHARACTER(len=1024) :: fileName
    INTEGER(i_kind) :: gTime(6), uTime_sec
    CHARACTER(*), INTENT(IN) :: varName        !< varName for output.
    TYPE(NCOutput_t) :: NCOutput
    CHARACTER(len=20) :: timeStr
    LOGICAL, OPTIONAL :: hasGLabel, hasMOTOR
    REAL(r_kind), ALLOCATABLE :: valueSel(:, :)
    REAL(r_kind), ALLOCATABLE :: valueHAS_MP(:, :)

    uTime_sec = X%sg%tt(tIdx)
    CALL Time_Unix_to_GMT(uTime_sec, gTime)
    ! PRINT *, 'gTime for this time slot is: ', gTime

1001 FORMAT(I4.4, I2.2, I2.2, A, I2.2, I2.2, A)
    WRITE (timeStr, 1001) gTime(1), gTime(2), gTime(3), '_', gTime(4), gTime(5), '_'

    fileName = TRIM(FNPrefix)//'_'//TRIM(varName)

    IF (PRESENT(hasGLabel)) THEN
      IF (hasGLabel) THEN
        BLOCK
          CHARACTER(LEN=2) :: iGrid
          WRITE (iGrid, "(I2.2)") X%sg%gLevel
          fileName = TRIM(fileName)//'_G'//iGrid
        END BLOCK
      END IF
    END IF

    IF (PRESENT(hasMOTOR)) THEN
      IF (hasMOTOR) THEN
        fileName = 'MOTOR-'//TRIM(fileName)
      END IF
    END IF

    fileName = TRIM(pathName)//'/'//TRIM(timeStr)//TRIM(fileName)//'.nc'

    ! PRINT *, 'fileName is: ', TRIM(fileName)
    ASSOCIATE (sg => X%sg)
      IF (X%getVarIdx(TRIM(varName)) .NE. 0) THEN

        IF (sg%isBaseProc()) NCOutput = NCOutput_t(fileName, this%lonEx, this%latEx, this%HgtAbv, sg%tt(tIdx:tIdx), &
                                                   this%NX, this%NY, this%vLevelOutput, 1)

        i = X%getVarIdx(TRIM(varName))
        IF (sg%isBaseProc()) THEN
          ALLOCATE (valueGlob(this%vLevelOutput, sg%num_icell_global))
          ALLOCATE (valueSel(this%vLevelOutput, this%NX * this%NY))
        END IF

        ALLOCATE (valueHAS_MP(this%vLevelOutput, sg%num_cell))

        DO j = 1, sg%num_cell
          CALL interp1d(sg%zHght(:, j), X%fields(i)%DATA(:, j, tIdx), &
                        this%HgtAbvMSL(:, j), valueHAS_MP(:, j))
          DO k = 1, this%vLevelOutput
            IF (this%HgtAbvMSL(k, j) .LT. sg%topo(j)) THEN
              valueHAS_MP(k, j) = 0.0D0
            END IF
          END DO
        END DO

        CALL sg%aggrGridReal2D(valueHAS_MP, &
                               valueGlob, [this%vLevelOutput, sg%num_icell_global])

        IF (sg%isBaseProc()) THEN
          FORALL (i=1:this%NX * this%NY)
            valueSel(:, i) = valueGlob(:, this%idxtgt(1, i)) * this%coetgt(1, i) + &
                             valueGlob(:, this%idxtgt(2, i)) * this%coetgt(2, i) + &
                             valueGlob(:, this%idxtgt(3, i)) * this%coetgt(3, i) + &
                             valueGlob(:, this%idxtgt(4, i)) * this%coetgt(4, i)
          END FORALL
          CALL NCOutput%addVar(X%fields(i)%Get_Name(), valueSel, &
                               this%vLevelOutput, this%NX, this%NY, 1)
        END IF

        ! Export height
        CALL sg%aggrGridReal2D(this%HgtAbvMSL, &
                               valueGlob, [this%vLevelOutput, sg%num_icell_global])

        IF (sg%isBaseProc()) THEN
          FORALL (i=1:this%NX * this%NY)
            valueSel(:, i) = valueGlob(:, this%idxtgt(1, i)) * this%coetgt(1, i) + &
                             valueGlob(:, this%idxtgt(2, i)) * this%coetgt(2, i) + &
                             valueGlob(:, this%idxtgt(3, i)) * this%coetgt(3, i) + &
                             valueGlob(:, this%idxtgt(4, i)) * this%coetgt(4, i)
          END FORALL
          CALL NCOutput%addVar('height', valueSel, &
                               this%vLevelOutput, this%NX, this%NY, 1)
        END IF

        IF (ALLOCATED(valueGlob)) DEALLOCATE (valueGlob)
        IF (ALLOCATED(valueSel)) DEALLOCATE (valueSel)
        IF (ALLOCATED(valueHAS_MP)) DEALLOCATE (valueHAS_MP)

        IF (sg%isBaseProc()) CALL NCOutput%CLOSE()
      END IF
    END ASSOCIATE
  END SUBROUTINE

!> @brief
!! Output state to nc file for last N time slot (all variables).
  SUBROUTINE Export_State_AVLNST(this, X, pathName, FNPrefix, ntSlots, hasGLabel, hasMOTOR)
    USE AdvanceTime_m

    CLASS(Export2HASCoordInSelDomain_t) :: this
    TYPE(State_t), INTENT(IN) :: X             !< State.
    CHARACTER(*), INTENT(IN) :: pathName       !< Path for output this state.
    CHARACTER(*), INTENT(IN) :: FNPrefix       !< Prefix of filename.
    INTEGER(i_kind), INTENT(IN) :: ntSlots     !< Index of time slot ouptut.
    LOGICAL, OPTIONAL :: hasGLabel, hasMOTOR

    INTEGER(i_kind) :: i

    IF (.NOT. X%sg%isActiveProc()) RETURN ! Return if it is not on the active procs

    IF (ntSlots > X%sg%tSlots) ERROR STOP "Output ntSlots error! "

    DO i = 1, ntSlots
      CALL this%Export_State_AVST(X, pathName, FNPrefix, X%sg%tSlots - (i - 1), hasGLabel, hasMOTOR)
    END DO

  END SUBROUTINE

!> @brief
!! Output state to nc file for last N time slot (single variable).
  SUBROUTINE Export_State_SVLNST(this, X, pathName, FNPrefix, varName, ntSlots, hasGLabel, hasMOTOR)
    USE AdvanceTime_m

    CLASS(Export2HASCoordInSelDomain_t) :: this
    TYPE(State_t), INTENT(IN) :: X             !< State.
    CHARACTER(*), INTENT(IN) :: pathName       !< Path for output this state.
    CHARACTER(*), INTENT(IN) :: FNPrefix       !< Prefix of filename.
    INTEGER(i_kind), INTENT(IN) :: ntSlots     !< Index of time slot ouptut.
    CHARACTER(*), INTENT(IN) :: varName        !< varName for output.
    LOGICAL, OPTIONAL :: hasGLabel, hasMOTOR

    INTEGER(i_kind) :: i

    IF (.NOT. X%sg%isActiveProc()) RETURN ! Return if it is not on the active procs

    IF (ntSlots > X%sg%tSlots) ERROR STOP "Output ntSlots error! "

    DO i = 1, ntSlots
      CALL this%Export_State_SVST(X, pathName, FNPrefix, varName, X%sg%tSlots - (i - 1), hasGLabel, hasMOTOR)
    END DO

  END SUBROUTINE

  SUBROUTINE destructor(this)
    TYPE(Export2HASCoordInSelDomain_t) :: this

    IF (ALLOCATED(this%latEx)) DEALLOCATE (this%latEx)
    IF (ALLOCATED(this%lonEx)) DEALLOCATE (this%lonEx)
    IF (ALLOCATED(this%llEx)) DEALLOCATE (this%llEx)
    IF (ALLOCATED(this%idxtgt)) DEALLOCATE (this%idxtgt)
    IF (ALLOCATED(this%coetgt)) DEALLOCATE (this%coetgt)
    IF (ALLOCATED(this%HgtAbvMSL)) DEALLOCATE (this%HgtAbvMSL)
    IF (ALLOCATED(this%HgtAbv)) DEALLOCATE (this%HgtAbv)
  END SUBROUTINE

END MODULE
