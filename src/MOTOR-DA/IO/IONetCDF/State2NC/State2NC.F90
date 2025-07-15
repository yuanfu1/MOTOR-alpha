!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE State2NC_m
  USE State_m, ONLY: State_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE NCOutput_m, ONLY: NCOutput_t
  USE parameters_m, ONLY: degree2radian

CONTAINS

!> @brief
!! Output a single variable in state to NC file.
  SUBROUTINE Output_NC_State_SV(X, pathName, FNPrefix, varName, hasGLabel, hasMOTOR)
    TYPE(State_t), INTENT(IN) :: X             !< State.
    CHARACTER(*), INTENT(IN) :: pathName       !< Path for output this state.
    CHARACTER(*), INTENT(IN) :: FNPrefix       !< Prefix of filename.
    CHARACTER(*), INTENT(IN) :: varName        !< varName for output.
    LOGICAL, OPTIONAL :: hasGLabel, hasMOTOR
    CHARACTER(len=1024) :: fileName
    REAL(r_kind), ALLOCATABLE :: valueState(:, :, :, :), valueGlob(:, :, :), valueGlob2D(:, :)
    TYPE(NCOutput_t) :: NCOutput

    fileName = TRIM(FNPrefix)//'_'//TRIM(varName)

    IF (PRESENT(hasGLabel)) THEN
      IF (hasGLabel) THEN
        BLOCK
          CHARACTER(LEN=2) :: iGrid
          WRITE (iGrid, "(I2.2)") X%sg%gLevel
          fileName = TRIM(fileName)//'_G'//iGrid//'.nc'
        END BLOCK
      END IF
    END IF

    IF (PRESENT(hasMOTOR)) THEN
      IF (hasMOTOR) THEN
        fileName = 'MOTOR-'//TRIM(fileName)
      END IF
    END IF

    fileName = TRIM(pathName)//'/'//TRIM(fileName)//'.nc'

    ASSOCIATE (sg => X%sg)
      IF (X%getVarIdx(TRIM(varName)) .NE. 0) THEN

        ALLOCATE (valueGlob(sg%vLevel, sg%num_icell_global, sg%tSlots))
        CALL sg%aggrGridRealForFieldGrid(X%fields(X%getVarIdx(TRIM(varName)))%DATA, &
                                         valueGlob, [sg%vLevel, sg%num_icell_global, sg%tSlots])

        IF (sg%isBaseProc()) THEN
          NCOutput = NCOutput_t(fileName, sg%lon1DAtBase, sg%lat1DAtBase, sg%sigma, sg%tt, &
                                sg%dimCell_global(2), sg%dimCell_global(1), sg%vLevel, sg%tSlots)

          ! valueState = reshape(valueGlob, (/sg%vLevel, sg%dimCell_global(2), sg%dimCell_global(1), sg%tSlots/));
          PRINT *, 'MAXVAL(valueGlob)-------- : ', MAXVAL(valueGlob)
          CALL NCOutput%addVar(varName, valueGlob, sg%vLevel, sg%dimCell_global(2), sg%dimCell_global(1), sg%tSlots)
        END IF
        DEALLOCATE (valueGlob)

        ! Export height
        ALLOCATE (valueGlob2D(sg%vLevel, sg%num_icell_global))
        CALL sg%aggrGridReal2D(sg%zHght, &
                               valueGlob2D, [sg%vLevel, sg%num_icell_global])

        IF (sg%isBaseProc()) THEN
          CALL NCOutput%addVar('height', valueGlob2D, &
                               sg%vLevel, sg%dimCell_global(2), sg%dimCell_global(1), 1)
          CALL NCOutput%CLOSE()
        END IF
        DEALLOCATE (valueGlob2D)
      END IF
    END ASSOCIATE
  END SUBROUTINE

  !> @brief
!! Output a single variable in state to NC file.
  SUBROUTINE Output_NC_State_AV(X, pathName, FNPrefix, hasGLabel, hasMOTOR)
    USE parameters_m, ONLY: degree2radian
    USE ncWriteGrid_m

    TYPE(State_t), INTENT(IN) :: X             !< State.
    CHARACTER(*), INTENT(IN) :: pathName       !< Path for output this state.
    CHARACTER(*), INTENT(IN) :: FNPrefix       !< Prefix of filename.
    LOGICAL, OPTIONAL :: hasGLabel, hasMOTOR
    CHARACTER(len=1024) :: fileName
    REAL(r_kind), ALLOCATABLE :: valueGlob(:, :, :), valueGlob2D(:, :)
    INTEGER :: i
    TYPE(NCOutput_t) :: NCOutput

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

    fileName = TRIM(pathName)//'/'//TRIM(fileName)//'.nc'

    ASSOCIATE (sg => X%sg)
      IF (sg%isBaseProc()) NCOutput = NCOutput_t(fileName, sg%lon1DAtBase, sg%lat1DAtBase, sg%sigma, sg%tt, &
                                                 sg%dimCell_global(2), sg%dimCell_global(1), sg%vLevel, sg%tSlots)
      ALLOCATE (valueGlob(sg%vLevel, sg%num_icell_global, sg%tSlots))

      DO i = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
        CALL sg%aggrGridRealForFieldGrid(X%fields(i)%DATA, &
                                         valueGlob, [sg%vLevel, sg%num_icell_global, sg%tSlots])

        IF (sg%isBaseProc()) THEN
          CALL NCOutput%addVar(X%fields(i)%Get_Name(), valueGlob, &
                               sg%vLevel, sg%dimCell_global(2), sg%dimCell_global(1), sg%tSlots)
        END IF

      END DO

      ! Export height
      ALLOCATE (valueGlob2D(sg%vLevel, sg%num_icell_global))
      CALL sg%aggrGridReal2D(sg%zHght, &
                             valueGlob2D, [sg%vLevel, sg%num_icell_global])

      IF (sg%isBaseProc()) THEN
        CALL NCOutput%addVar('height', valueGlob2D, &
                             sg%vLevel, sg%dimCell_global(2), sg%dimCell_global(1), 1)
      END IF
      DEALLOCATE (valueGlob2D)

      DEALLOCATE (valueGlob)
      IF (sg%isBaseProc()) CALL NCOutput%CLOSE()
    END ASSOCIATE
  END SUBROUTINE

!> @brief
!! Output state to nc file for single time slot (all variables).
  SUBROUTINE Output_NC_State_AVST(X, pathName, FNPrefix, tIdx, hasGLabel, hasMOTOR)
    USE parameters_m, ONLY: degree2radian
    USE ncWriteGrid_m
    USE AdvanceTime_m

    TYPE(State_t), INTENT(IN) :: X             !< State.
    CHARACTER(*), INTENT(IN) :: pathName       !< Path for output this state.
    CHARACTER(*), INTENT(IN) :: FNPrefix       !< Prefix of filename.
    INTEGER(i_kind), INTENT(IN) :: tIdx        !< Index of time slot ouptut.

    REAL(r_kind), ALLOCATABLE :: valueGlob(:, :, :), valueGlob2D(:, :)
    CHARACTER(len=1024) :: fileName
    INTEGER(i_kind) :: gTime(6), uTime_sec, i
    TYPE(NCOutput_t) :: NCOutput
    LOGICAL, OPTIONAL :: hasGLabel, hasMOTOR
    CHARACTER(len=20) :: timeStr

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
      IF (sg%isBaseProc()) NCOutput = NCOutput_t(fileName, sg%lon1DAtBase, sg%lat1DAtBase, sg%sigma, sg%tt(tIdx:tIdx), &
                                                 sg%dimCell_global(2), sg%dimCell_global(1), sg%vLevel, 1)

      ALLOCATE (valueGlob(sg%vLevel, sg%num_icell_global, 1))
      DO i = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
        CALL sg%aggrGridRealForFieldGrid(X%fields(i)%DATA(:, :, tIdx:tIdx), &
                                         valueGlob, [sg%vLevel, sg%num_icell_global, 1])

        IF (sg%isBaseProc()) THEN
          CALL NCOutput%addVar(X%fields(i)%Get_Name(), valueGlob, &
                               sg%vLevel, sg%dimCell_global(2), sg%dimCell_global(1), 1)
        END IF

      END DO

      ! Export height
      ALLOCATE (valueGlob2D(sg%vLevel, sg%num_icell_global))
      CALL sg%aggrGridReal2D(sg%zHght, &
                             valueGlob2D, [sg%vLevel, sg%num_icell_global])

      IF (sg%isBaseProc()) THEN
        CALL NCOutput%addVar('height', valueGlob2D, &
                             sg%vLevel, sg%dimCell_global(2), sg%dimCell_global(1), 1)
      END IF
      DEALLOCATE (valueGlob2D)

      DEALLOCATE (valueGlob)
      IF (sg%isBaseProc()) CALL NCOutput%CLOSE()
    END ASSOCIATE

  END SUBROUTINE

!> @brief
!! Output state to nc file for single time slot (single variable).
  SUBROUTINE Output_NC_State_SVST(X, pathName, FNPrefix, varName, tIdx, hasGLabel, hasMOTOR)
    USE parameters_m, ONLY: degree2radian
    USE ncWriteGrid_m
    USE AdvanceTime_m

    TYPE(State_t), INTENT(IN) :: X             !< State.
    CHARACTER(*), INTENT(IN) :: pathName       !< Path for output this state.
    CHARACTER(*), INTENT(IN) :: FNPrefix       !< Prefix of filename.
    INTEGER(i_kind), INTENT(IN) :: tIdx        !< Index of time slot ouptut.

    REAL(r_kind), ALLOCATABLE :: valueGlob(:, :, :), valueGlob2D(:, :)
    INTEGER :: i, j, k, l
    CHARACTER(len=1024) :: fileName
    INTEGER(i_kind) :: gTime(6), uTime_sec
    CHARACTER(*), INTENT(IN) :: varName        !< varName for output.
    TYPE(NCOutput_t) :: NCOutput
    CHARACTER(len=20) :: timeStr
    LOGICAL, OPTIONAL :: hasGLabel, hasMOTOR

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

        IF (sg%isBaseProc()) NCOutput = NCOutput_t(fileName, sg%lon1DAtBase, sg%lat1DAtBase, sg%sigma, sg%tt(tIdx:tIdx), &
                                                   sg%dimCell_global(2), sg%dimCell_global(1), sg%vLevel, 1)

        i = X%getVarIdx(TRIM(varName))
        ALLOCATE (valueGlob(sg%vLevel, sg%num_icell_global, 1))
        CALL sg%aggrGridRealForFieldGrid(X%fields(i)%DATA(:, :, tIdx:tIdx), &
                                         valueGlob, [sg%vLevel, sg%num_icell_global, 1])

        IF (sg%isBaseProc()) THEN
          CALL NCOutput%addVar(X%fields(i)%Get_Name(), valueGlob, &
                               sg%vLevel, sg%dimCell_global(2), sg%dimCell_global(1), 1)
        END IF

        ! Export height
        ALLOCATE (valueGlob2D(sg%vLevel, sg%num_icell_global))
        CALL sg%aggrGridReal2D(sg%zHght, &
                               valueGlob2D, [sg%vLevel, sg%num_icell_global])

        IF (sg%isBaseProc()) THEN
          CALL NCOutput%addVar('height', valueGlob2D, &
                               sg%vLevel, sg%dimCell_global(2), sg%dimCell_global(1), 1)
        END IF
        DEALLOCATE (valueGlob2D)

        DEALLOCATE (valueGlob)

        IF (sg%isBaseProc()) CALL NCOutput%CLOSE()
      END IF
    END ASSOCIATE
  END SUBROUTINE

!> @brief
!! Output state to nc file for last N time slot (all variables).
  SUBROUTINE Output_NC_State_AVLNST(X, pathName, FNPrefix, ntSlots, hasGLabel, hasMOTOR)
    USE parameters_m, ONLY: degree2radian
    USE ncWriteGrid_m
    USE AdvanceTime_m

    TYPE(State_t), INTENT(IN) :: X             !< State.
    CHARACTER(*), INTENT(IN) :: pathName       !< Path for output this state.
    CHARACTER(*), INTENT(IN) :: FNPrefix       !< Prefix of filename.
    INTEGER(i_kind), INTENT(IN) :: ntSlots     !< Index of time slot ouptut.
    LOGICAL, OPTIONAL :: hasGLabel, hasMOTOR

    INTEGER(i_kind) :: i

    IF (.NOT. X%sg%isActiveProc()) RETURN ! Return if it is not on the active procs

    IF (ntSlots > X%sg%tSlots) ERROR STOP "Output ntSlots error! "

    DO i = 1, ntSlots
      CALL Output_NC_State_AVST(X, pathName, FNPrefix, X%sg%tSlots - (i - 1), hasGLabel, hasMOTOR)
    END DO

  END SUBROUTINE

  !> @brief
!! Output state to nc file for last N time slot (single variable).
  SUBROUTINE Output_NC_State_SVLNST(X, pathName, FNPrefix, varName, ntSlots, hasGLabel, hasMOTOR)
    USE parameters_m, ONLY: degree2radian
    USE ncWriteGrid_m
    USE AdvanceTime_m

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
      CALL Output_NC_State_SVST(X, pathName, FNPrefix, varName, X%sg%tSlots - (i - 1), hasGLabel, hasMOTOR)
    END DO

  END SUBROUTINE
END MODULE State2NC_m
