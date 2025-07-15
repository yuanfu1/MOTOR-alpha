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
MODULE Export2SelDomain_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE YAMLRead_m
  USE State_m, ONLY: State_t
  USE parameters_m, ONLY: degree2radian
  USE NCOutput_m, ONLY: NCOutput_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE RectangleGridUtility_m, ONLY: RectangleGridUtility_t
  IMPLICIT NONE

  TYPE Export2SelDomain_t
    CHARACTER(len=1024) :: configFile
    REAL(r_kind), ALLOCATABLE :: latEx(:)
    REAL(r_kind), ALLOCATABLE :: lonEx(:)
    REAL(r_kind), ALLOCATABLE :: llEx(:, :)
    INTEGER(i_kind) :: NX, NY
    INTEGER(i_kind), ALLOCATABLE :: idxtgt(:, :)
    REAL(r_kind), ALLOCATABLE :: coetgt(:, :)

  CONTAINS
    FINAL :: destructor
    PROCEDURE, PUBLIC :: Export_State_AVST
    PROCEDURE, PUBLIC :: Export_State_AVLNST
    PROCEDURE, PUBLIC :: Export_State_SVST
    PROCEDURE, PUBLIC :: Export_State_SVLNST
  END TYPE

  INTERFACE Export2SelDomain_t
    PROCEDURE :: constructor
  END INTERFACE

CONTAINS

  FUNCTION constructor(configFile, sg) RESULT(this)
    IMPLICIT NONE
    TYPE(Export2SelDomain_t) :: this
    CHARACTER(len=1024), INTENT(IN) :: configFile
    TYPE(SingleGrid_t), INTENT(IN) :: sg
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
  END

!> @brief
!! Output state to nc file for single time slot (all variables).
  SUBROUTINE Export_State_AVST(this, X, pathName, FNPrefix, tIdx, hasGLabel, hasMOTOR)
    USE AdvanceTime_m

    CLASS(Export2SelDomain_t) :: this
    TYPE(State_t), INTENT(IN) :: X             !< State.
    CHARACTER(*), INTENT(IN) :: pathName       !< Path for output this state.
    CHARACTER(*), INTENT(IN) :: FNPrefix       !< Prefix of filename.
    INTEGER(i_kind), INTENT(IN) :: tIdx        !< Index of time slot ouptut.

    REAL(r_kind), ALLOCATABLE :: valueGlob(:, :, :)
    CHARACTER(len=1024) :: fileName
    INTEGER(i_kind) :: gTime(6), uTime_sec, i
    TYPE(NCOutput_t) :: NCOutput
    LOGICAL, OPTIONAL :: hasGLabel, hasMOTOR
    CHARACTER(len=20) :: timeStr
    REAL(r_kind), ALLOCATABLE :: valueSel(:, :)

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
        NCOutput = NCOutput_t(fileName, this%lonEx, this%latEx, sg%sigma, sg%tt(tIdx:tIdx), &
                              this%NX, this%NY, sg%vLevel, 1)
      END IF

      DO i = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
        IF (sg%isBaseProc()) THEN
          ALLOCATE (valueGlob(sg%vLevel, sg%num_icell_global, 1))
          ALLOCATE (valueSel(sg%vLevel, this%NX * this%NY))
        END IF

        CALL sg%aggrGridRealForFieldGrid(X%fields(i)%DATA(:, :, tIdx:tIdx), &
                                         valueGlob, [sg%vLevel, sg%num_icell_global, 1])

        IF (sg%isBaseProc()) THEN
          FORALL (i=1:this%NX * this%NY)
            valueSel(:, i) = valueGlob(:, this%idxtgt(1, i), 1) * this%coetgt(1, i) + &
                             valueGlob(:, this%idxtgt(2, i), 1) * this%coetgt(2, i) + &
                             valueGlob(:, this%idxtgt(3, i), 1) * this%coetgt(3, i) + &
                             valueGlob(:, this%idxtgt(4, i), 1) * this%coetgt(4, i)
          END FORALL
          CALL NCOutput%addVar(X%fields(i)%Get_Name(), valueSel, &
                               sg%vLevel, this%NX, this%NY, 1)
        END IF

        IF (ALLOCATED(valueGlob)) DEALLOCATE (valueGlob)
        IF (ALLOCATED(valueSel)) DEALLOCATE (valueSel)

      END DO
      IF (sg%isBaseProc()) CALL NCOutput%CLOSE()
    END ASSOCIATE

  END SUBROUTINE

!> @brief
!! Output state to nc file for single time slot (single variable).
  SUBROUTINE Export_State_SVST(this, X, pathName, FNPrefix, varName, tIdx, hasGLabel, hasMOTOR)
    USE AdvanceTime_m

    CLASS(Export2SelDomain_t) :: this
    TYPE(State_t), INTENT(IN) :: X             !< State.
    CHARACTER(*), INTENT(IN) :: pathName       !< Path for output this state.
    CHARACTER(*), INTENT(IN) :: FNPrefix       !< Prefix of filename.
    INTEGER(i_kind), INTENT(IN) :: tIdx        !< Index of time slot ouptut.

    REAL(r_kind), ALLOCATABLE :: valueGlob(:, :, :)
    INTEGER :: i, j, k, l
    CHARACTER(len=1024) :: fileName
    INTEGER(i_kind) :: gTime(6), uTime_sec
    CHARACTER(*), INTENT(IN) :: varName        !< varName for output.
    TYPE(NCOutput_t) :: NCOutput
    CHARACTER(len=20) :: timeStr
    LOGICAL, OPTIONAL :: hasGLabel, hasMOTOR
    REAL(r_kind), ALLOCATABLE :: valueSel(:, :)

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

        IF (sg%isBaseProc()) NCOutput = NCOutput_t(fileName, this%lonEx, this%latEx, sg%sigma, sg%tt(tIdx:tIdx), &
                                                   this%NX, this%NY, sg%vLevel, 1)

        i = X%getVarIdx(TRIM(varName))
        IF (sg%isBaseProc()) THEN
          ALLOCATE (valueGlob(sg%vLevel, sg%num_icell_global, 1))
          ALLOCATE (valueSel(sg%vLevel, this%NX * this%NY))
        END IF

        CALL sg%aggrGridRealForFieldGrid(X%fields(i)%DATA(:, :, tIdx:tIdx), &
                                         valueGlob, [sg%vLevel, sg%num_icell_global, 1])

        IF (sg%isBaseProc()) THEN
          FORALL (i=1:this%NX * this%NY)
            valueSel(:, i) = valueGlob(:, this%idxtgt(1, i), 1) * this%coetgt(1, i) + &
                             valueGlob(:, this%idxtgt(2, i), 1) * this%coetgt(2, i) + &
                             valueGlob(:, this%idxtgt(3, i), 1) * this%coetgt(3, i) + &
                             valueGlob(:, this%idxtgt(4, i), 1) * this%coetgt(4, i)
          END FORALL
          CALL NCOutput%addVar(X%fields(i)%Get_Name(), valueSel, &
                               sg%vLevel, this%NX, this%NY, 1)
        END IF

        IF (ALLOCATED(valueGlob)) DEALLOCATE (valueGlob)
        IF (ALLOCATED(valueSel)) DEALLOCATE (valueSel)

        IF (sg%isBaseProc()) CALL NCOutput%CLOSE()
      END IF
    END ASSOCIATE
  END SUBROUTINE

!> @brief
!! Output state to nc file for last N time slot (all variables).
  SUBROUTINE Export_State_AVLNST(this, X, pathName, FNPrefix, ntSlots, hasGLabel, hasMOTOR)
    USE AdvanceTime_m

    CLASS(Export2SelDomain_t) :: this
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

    CLASS(Export2SelDomain_t) :: this
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
    TYPE(Export2SelDomain_t) :: this

    IF (ALLOCATED(this%latEx)) DEALLOCATE (this%latEx)
    IF (ALLOCATED(this%lonEx)) DEALLOCATE (this%lonEx)
    IF (ALLOCATED(this%llEx)) DEALLOCATE (this%llEx)
    IF (ALLOCATED(this%idxtgt)) DEALLOCATE (this%idxtgt)
    IF (ALLOCATED(this%coetgt)) DEALLOCATE (this%coetgt)

  END SUBROUTINE

END MODULE
