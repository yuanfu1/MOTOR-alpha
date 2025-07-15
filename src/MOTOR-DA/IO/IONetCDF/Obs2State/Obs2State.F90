!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2022/3/4, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE Obs2State_m
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t

CONTAINS
  FUNCTION Obs2State_BaseTypeName(sg, Y) RESULT(X)
    TYPE(SingleGrid_t) :: sg
    TYPE(State_t) :: X
    TYPE(ObsSet_t) :: Y
    CHARACTER(len=25), ALLOCATABLE :: varList(:)
    INTEGER :: i, j, k

    ALLOCATE (varList(UBOUND(Y%ObsFields, 1)))

    DO i = 1, UBOUND(Y%ObsFields, 1)
      varList(i) = TRIM(Y%ObsFields(i)%Get_Id_Name())
      ! PRINT *,'Obs2State_BaseTypeName - varList(i): ', Y%ObsFields(i)%Get_Id_Name()
    END DO

    CALL X%ctorVarList(varList, sg)

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      DO i = LBOUND(Y%ObsFields, 1), UBOUND(Y%ObsFields, 1)
        IF ((TRIM(Y%ObsFields(i)%Get_Id_Name())) .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN
          DO k = LBOUND(Y%ObsFields(i)%idx, 1), UBOUND(Y%ObsFields(i)%idx, 1)
            CALL X%fields(j)%Set_Value(Y%ObsFields(i)%idx(k), Y%ObsFields(i)%values(k))
          END DO
        END IF
      END DO
    END DO

    DEALLOCATE (varList)
  END FUNCTION

  FUNCTION Obs2State_BaseName(sg, Y) RESULT(X)
    TYPE(SingleGrid_t) :: sg
    TYPE(State_t) :: X
    TYPE(ObsSet_t) :: Y
    CHARACTER(len=25), ALLOCATABLE :: varList(:)
    INTEGER :: i, j, k

    ALLOCATE (varList(UBOUND(Y%ObsFields, 1)))

    DO i = 1, UBOUND(Y%ObsFields, 1)
      varList(i) = TRIM(Y%ObsFields(i)%Get_Name())
    END DO

    CALL X%ctorVarList(varList, sg)

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      DO i = LBOUND(Y%ObsFields, 1), UBOUND(Y%ObsFields, 1)
        IF ((TRIM(Y%ObsFields(i)%Get_Name())) .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN
          DO k = LBOUND(Y%ObsFields(i)%idx, 1), UBOUND(Y%ObsFields(i)%idx, 1)

            !  PRINT *,'X%fields(j)%Get_Name(): ',X%fields(j)%Get_Name(), sg%gLevel
            CALL X%fields(j)%Set_Value(Y%ObsFields(i)%idx(k), Y%ObsFields(i)%values(k))
          END DO
        END IF
      END DO
    END DO

    DEALLOCATE (varList)
  END FUNCTION

END MODULE Obs2State_m
