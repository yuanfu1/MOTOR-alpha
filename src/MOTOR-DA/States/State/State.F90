!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.State_Xm
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : Beta 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2020/11/19, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/11/01, @GBA-MWF, Shenzhen
!     for adding state filter operator.
!!--------------------------------------------------------------------------------------------------

! Contains states of Xm, etc.
MODULE State_m
  ! IMPLICIT NONE
  USE kinds_m, ONLY: i_kind, r_kind
  USE singleGrid_m, ONLY: singleGrid_t
  USE field_m, ONLY: field_t !, field2D_t
  USE AuxTypeSG_m, ONLY: AuxTypeSG_t
  USE Filter_m, ONLY: smoothField
  USE conversions_m, ONLY: Saturated_Vapor
  USE parameters_m
  

  TYPE, EXTENDS(AuxTypeSG_t) :: State_t
    TYPE(field_t), ALLOCATABLE :: fields(:)
  CONTAINS
    FINAL :: destructor
    PROCEDURE, PUBLIC :: multiply
    PROCEDURE, PUBLIC :: dot_multiply
    PROCEDURE, PUBLIC :: divideByValue
    PROCEDURE, PUBLIC :: divideByVec
    PROCEDURE, PUBLIC :: add
    PROCEDURE, PUBLIC :: subtract
    PROCEDURE, PUBLIC :: getVarIdx
    PROCEDURE, PUBLIC :: exHalo
    PROCEDURE, PUBLIC :: exHaloRevSum
    PROCEDURE, PUBLIC :: cleanHalo
    PROCEDURE, PUBLIC :: zeroCopy
    PROCEDURE, PUBLIC :: setAllFieldData
    PROCEDURE, PUBLIC :: getVar
    PROCEDURE, PUBLIC :: addVar
    PROCEDURE, PUBLIC :: rmVar
    PROCEDURE, PUBLIC :: showInfo
    PROCEDURE, PUBLIC :: getMeanSeaLevelTempaure
    PROCEDURE, PUBLIC :: clearHalo
    PROCEDURE, PUBLIC :: fillPresWithHydrostatic

    PROCEDURE, PUBLIC :: hStateFilter   ! Horizontal filter calling a SingleGrid operator

    GENERIC :: OPERATOR(.DOT.) => dot_multiply
    GENERIC :: OPERATOR(*) => multiply
    GENERIC :: OPERATOR(+) => add
    GENERIC :: OPERATOR(-) => subtract
    GENERIC :: OPERATOR(/) => divideByValue, divideByVec

    PROCEDURE, PUBLIC :: initialize, ctorVarList

  END TYPE State_t

  ! INTERFACE State_t
  !   PROCEDURE :: constructor
  !   PROCEDURE :: ctorVarList
  ! END INTERFACE State_t

CONTAINS

  !> @brief
  !! Constructor 1
  !! For a model initial conditions, it requires a few time frames so Yuanfu Xie added an
  !! option to specify a state with far less time frames.
  SUBROUTINE initialize(this, configFile, sg, numTimes)
    USE namelist_State_m, ONLY: namelist_State

    CLASS(State_t) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    INTEGER(i_kind), INTENT(IN), OPTIONAL :: numTimes !< Yuanfu Xie added this for model IC

    ! Local variables:
    CHARACTER(LEN=10), ALLOCATABLE :: varList(:)
    INTEGER(i_kind) :: istatus
    INTEGER(i_kind), ALLOCATABLE :: smoother(:, :)   ! smoothness of var, multigrid level

    ! this%AuxTypeSG_t = AuxTypeSG_t(sg)
    CALL this%AuxTypeSG_t%aux_initialize(sg)
    CALL namelist_State(configFile, varList)

    BLOCK
      INTEGER(i_kind) :: i
      ALLOCATE (this%fields(UBOUND(varList, 1)))

      WRITE (*, 1) UBOUND(varList, 1), (TRIM(varList(i)), i=1, UBOUND(varList, 1))
1     FORMAT('Model Vars: ', I2, ' VarList ', 20A8)
      ! The fields input from the config files
      DO i = 1, UBOUND(varList, 1)
        IF (TRIM(varList(i)) .EQ. 'psfc0') THEN
          CALL this%fields(i)%initialize(this%sg, varList(i), '3DS')
        ELSE
          IF (PRESENT(numTimes)) THEN
            CALL this%fields(i)%initialize(this%sg, varList(i), '4DU', numTimes)
          ELSE
            CALL this%fields(i)%initialize(this%sg, varList(i), '4DU')
          END IF
        END IF
      END DO
    END BLOCK

    IF (ALLOCATED(varList)) DEALLOCATE (varList)
  END SUBROUTINE

  SUBROUTINE clearHalo(this)
    CLASS(State_t) :: this

    DO i = 1, UBOUND(this%fields, 1)
      this%fields(i)%DATA(:, this%sg%num_icell + 1:this%sg%num_cell, :) = 0.0D0
    END DO
  END SUBROUTINE

  SUBROUTINE ctorVarList(this, varList, sg)
    USE namelist_State_m, ONLY: namelist_State

    CLASS(State_t) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    CHARACTER(LEN=*), ALLOCATABLE :: varList(:)

    !this%AuxTypeSG_t = AuxTypeSG_t(sg)
    CALL this%AuxTypeSG_t%aux_initialize(sg)

    BLOCK
      INTEGER(i_kind) :: i
      ALLOCATE (this%fields(UBOUND(varList, 1)))

      ! The fields input from the config files
      DO i = 1, UBOUND(varList, 1)
        IF (TRIM(varList(i)) .EQ. 'psfc0') THEN
          CALL this%fields(i)%initialize(this%sg, varList(i), '3DS')
        ELSE
          CALL this%fields(i)%initialize(this%sg, varList(i), '4DU')
        END IF
      END DO
    END BLOCK
  END SUBROUTINE

  FUNCTION getVar(this, varName) RESULT(var)
    CLASS(State_t), TARGET :: this
    CHARACTER(*) :: varName
    TYPE(Field_t), POINTER :: var

    DO i = LBOUND(this%fields, 1), UBOUND(this%fields, 1)
      IF (TRIM(this%fields(i)%Get_Name()) .EQ. TRIM(varName)) THEN
        var => this%fields(i)
        RETURN
      END IF
    END DO
  END FUNCTION getVar

  SUBROUTINE rmVar(this, varName)
    CLASS(State_t) :: this
    CHARACTER(*) :: varName
    INTEGER(i_kind) :: i, j, varIdx
    TYPE(field_t), ALLOCATABLE :: temp(:)

    varIdx = this%getVarIdx(varName)
    IF (varIdx /= 0) THEN
      temp = this%fields
      DEALLOCATE (this%fields)
      ALLOCATE (this%fields(SIZE(temp) - 1))
      j = 1
      DO i = 1, SIZE(temp)
        IF (i .NE. varIdx) THEN
          this%fields(j) = temp(i)
          j = j + 1
        END IF
      END DO
      DEALLOCATE (temp)
    END IF

    ! &
    ! this%fields = pack(this%fields, (/(i, i=1, SIZE(this%fields))/) /= this%getVarIdx(varName))
  END SUBROUTINE

  SUBROUTINE setAllFieldData(this, val)
    CLASS(State_t) :: this
    REAL(r_kind) :: val
    INTEGER(i_kind) :: i

    DO i = LBOUND(this%Fields, 1), UBOUND(this%Fields, 1)
      this%fields(i)%DATA = val
    END DO
  END SUBROUTINE

  SUBROUTINE addVar(this, varName)
    CLASS(State_t) :: this
    CHARACTER(*) :: varName
    TYPE(field_t), ALLOCATABLE :: temp(:)

    IF (this%getVarIdx(varName) == 0) THEN
      temp = this%fields
      DEALLOCATE (this%fields)
      ALLOCATE (this%fields(SIZE(temp) + 1))
      this%fields(1:SIZE(temp)) = temp
      this%fields(SIZE(temp) + 1) = field_t(this%sg, varName, '4DU')
      DEALLOCATE (temp)
    END IF

  END SUBROUTINE

  SUBROUTINE exHalo(this, varName)
    CLASS(State_t) :: this
    INTEGER(i_kind) :: i
    CHARACTER(len=*), OPTIONAL :: varName

    IF (PRESENT(varName)) THEN
      CALL this%sg%ExchangeMatOnHaloForFieldGrid(this%sg%tSlots, this%sg%vLevel, this%Fields(this%getVarIdx(TRIM(varName)))%DATA)
    ELSE
      DO k = LBOUND(this%Fields, 1), UBOUND(this%Fields, 1)
        CALL this%sg%ExchangeMatOnHaloForFieldGrid(this%sg%tSlots, this%sg%vLevel, this%Fields(k)%DATA)
      END DO
    END IF
  END SUBROUTINE

  SUBROUTINE cleanHalo(this, varName)
    CLASS(State_t) :: this
    INTEGER(i_kind) :: i
    CHARACTER(len=*), OPTIONAL :: varName

    IF (PRESENT(varName)) THEN
      DO i = this%sg%num_icell + 1, this%sg%num_cell
        this%Fields(this%getVarIdx(TRIM(varName)))%DATA(:, i, :) = 0.0D0
      END DO
    ELSE
      DO k = LBOUND(this%Fields, 1), UBOUND(this%Fields, 1)
        DO i = this%sg%num_icell + 1, this%sg%num_cell
          this%Fields(k)%DATA(:, i, :) = 0.0D0
        END DO
      END DO
    END IF
  END SUBROUTINE
  SUBROUTINE exHaloRevSum(this)
    CLASS(State_t) :: this
    INTEGER(i_kind) :: i

    DO k = LBOUND(this%Fields, 1), UBOUND(this%Fields, 1)
      CALL this%sg%ExchangeMatOnHaloReverseSumForFieldGrid(this%sg%tSlots, this%sg%vLevel, this%Fields(k)%DATA)
    END DO
  END SUBROUTINE

  ! @brief Yuanfu Xie modified this function to allow wind speed and direction as obs:
  FUNCTION getVarIdx(this, varName) RESULT(idx)
    CLASS(State_t) :: this
    CHARACTER(*) :: varName
    INTEGER(i_kind) :: idx
    INTEGER(i_kind) :: i

    idx = 0
    ! Yuanfu Xie added wspd and wdir option, 2022-08-01, as state_t uses uwnd and vwnd for wind states:
    IF (TRIM(varName) .EQ. 'wspd' .OR. TRIM(varName) .EQ. 'cdir' .OR. TRIM(varName) .EQ. 'sdir') THEN
      DO i = LBOUND(this%fields, 1), UBOUND(this%fields, 1)
        IF ((TRIM(this%fields(i)%Get_Name()) .EQ. 'uwnd' .AND. &
             TRIM(varName) .EQ. 'wspd') .OR. &
            (TRIM(this%fields(i)%Get_Name()) .EQ. 'vwnd' .AND. &
             (TRIM(varName) .EQ. 'cdir') .OR. TRIM(varName) .EQ. 'sdir')) THEN
          idx = i
          RETURN
        END IF
      END DO
    ELSE IF (TRIM(varName) .EQ. 'refractivity') THEN
      !! for refractivity
      DO i = LBOUND(this%fields, 1), UBOUND(this%fields, 1)
        IF (TRIM(this%fields(i)%Get_Name()) .EQ. 'temp' .OR. &
            TRIM(this%fields(i)%Get_Name()) .EQ. 'qvapor' .OR. &
            TRIM(this%fields(i)%Get_Name()) .EQ. 'pres' &
            ) THEN
          idx = i
          RETURN
        END IF
      END DO
    ELSE
      DO i = LBOUND(this%fields, 1), UBOUND(this%fields, 1)
        IF (this%fields(i)%Get_Name() .EQ. (varName)) THEN
          idx = i
          RETURN
        END IF
      END DO
    END IF
  END FUNCTION

  FUNCTION add(X1, X2) RESULT(X3)
    CLASS(State_t), INTENT(IN) :: X1
    TYPE(State_t), INTENT(IN) :: X2
    TYPE(State_t) :: X3

    X3 = X1
    DO i = LBOUND(X1%fields, 1), UBOUND(X1%fields, 1)
      X3%fields(i) = X1%fields(i) + X2%fields(i)
    END DO

  END FUNCTION add

  FUNCTION multiply(X1, mul) RESULT(X2)
    CLASS(State_t), INTENT(IN) :: X1
    REAL(r_kind), INTENT(IN) :: mul
    TYPE(State_t) :: X2

    X2 = X1
    DO i = LBOUND(X1%fields, 1), UBOUND(X1%fields, 1)
      X2%fields(i) = X1%fields(i) * mul
    END DO

  END FUNCTION

  FUNCTION divideByVec(X1, X0) RESULT(X2)
    CLASS(State_t), INTENT(IN) :: X1
    TYPE(State_t), INTENT(IN) :: X0
    TYPE(State_t) :: X2

    X2 = X1
    DO i = LBOUND(X1%fields, 1), UBOUND(X1%fields, 1)
      X2%fields(i) = X1%fields(i) / X0%fields(i)
    END DO
  END FUNCTION

  FUNCTION divideByValue(X1, divider) RESULT(X2)
    CLASS(State_t), INTENT(IN) :: X1
    REAL(r_kind), INTENT(IN) :: divider
    TYPE(State_t) :: X2

    X2 = X1
    DO i = LBOUND(X1%fields, 1), UBOUND(X1%fields, 1)
      X2%fields(i) = X1%fields(i) / divider
    END DO
  END FUNCTION

  FUNCTION subtract(X1, X2) RESULT(X3)
    CLASS(State_t), INTENT(IN) :: X1
    TYPE(State_t), INTENT(IN) :: X2
    TYPE(State_t) ::X3

    X3 = X1
    DO i = LBOUND(X3%Fields, 1), UBOUND(X3%Fields, 1)
      X3%Fields(i) = X2%Fields(i) - X1%Fields(i)
    END DO

  END FUNCTION subtract

  REAL(r_kind) FUNCTION dot_multiply(X1, X2)
    CLASS(State_t), INTENT(IN) :: X1
    TYPE(State_t), INTENT(IN) :: X2
    INTEGER(i_kind) :: i, j

    dot_multiply = 0.0D0
    DO i = LBOUND(X1%fields, 1), UBOUND(X1%fields, 1)
      dot_multiply = (X1%fields(i) .DOT.X2%fields(i)) + dot_multiply
    END DO
  END FUNCTION

  FUNCTION zeroCopy(this) RESULT(Xout)
    CLASS(State_t) :: this
    TYPE(State_t) :: Xout
    INTEGER :: i

    Xout = this
    DO i = LBOUND(Xout%fields, 1), UBOUND(Xout%fields, 1)
      Xout%fields(i)%DATA = 0.0D0
    END DO
  END FUNCTION

  SUBROUTINE hStateFilter(this, smoother)
    CLASS(State_t) :: this
    INTEGER(i_kind), INTENT(IN) :: smoother(*)

    ! Local variables:
    INTEGER(i_kind) :: i

    ! For each field applying filter:
    DO i = 1, UBOUND(this%fields, 1)
      CALL smoothField(this%sg, this%sg%vLevel, this%sg%tSlots, smoother(i), &
                       this%fields(i)%DATA, this%fields(i)%DATA)
    END DO
  END SUBROUTINE hStateFilter

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(State_t), INTENT(INOUT) :: this

    ! Release the fields in Xm
    IF (ALLOCATED(this%fields)) DEALLOCATE (this%fields)
  END SUBROUTINE destructor

  SUBROUTINE showInfo(this)
    CLASS(State_t) :: this
    PRINT *, 'State var List:'
    DO i = LBOUND(this%Fields, 1), UBOUND(this%Fields, 1)
      PRINT *, this%Fields(i)%Get_Name()
    END DO
  END SUBROUTINE

  FUNCTION getMeanSeaLevelTempaure(this) RESULT(TMSL)
    CLASS(State_t) :: this
    REAL(r_kind) :: sumAllAtEachProc, TMSL, sumAll
    INTEGER(i_kind) :: i, j, it
    INTEGER(i_kind) :: countAll, countAllAtEachProc

    TMSL = 0.0D0
    IF (.NOT. this%mpddSub%isActiveProc()) RETURN

    IF (this%getVarIdx('temp') == 0) THEN
      PRINT *, 'No temperature field in the state'
      STOP
    END IF
    it = this%getVarIdx('temp')
    ! Note: without these precalculated indices, some times the compiler complains about the getVarIdx Yuanfu Xie 2025-05-19
    ASSOCIATE (temp => this%fields(it)%DATA)
      sumAllAtEachProc = 0.0D0
      countAllAtEachProc = 0

      DO i = 1, this%sg%num_icell
        DO j = 1, this%sg%tSlots
          IF (ABS(this%sg%zHght(1, i)) < 0.000000001) THEN
            sumAllAtEachProc = sumAllAtEachProc + temp(1, i, j)
            ! sumAllAtEachProc = sumAllAtEachProc + this%fields(this%getVarIdx('temp'))%DATA(1, i, j)
            countAllAtEachProc = countAllAtEachProc + 1
          END IF
        END DO
      END DO

      CALL this%mpddSub%AllReduceSumReal(sumAllAtEachProc, sumAll)
      CALL this%mpddSub%AllReduceSumInt(countAllAtEachProc, countAll)

      PRINT *, 'sumAll:', sumAll
      PRINT *, 'countAll:', countAll

      TMSL = sumAll / countAll; 
    END ASSOCIATE

  END FUNCTION

  SUBROUTINE fillPresWithHydrostatic(this)
    CLASS(State_t) :: this
    REAL(r_kind), ALLOCATABLE :: swapDA(:,:,:)
    INTEGER(i_kind) :: j, k, ip, iq, it     ! Yuanfu Xie added the ip, iq, it for the indices of the fields on 2025-05-19

    IF (this%getVarIdx('pres') == 0 .OR. this%getVarIdx('temp') == 0 .OR. this%getVarIdx('qvapor') == 0 ) THEN
      RETURN
    END IF

    ALLOCATE (swapDA(this%sg%vLevel, this%sg%num_cell, this%sg%tSlots))
    swapDA = this%Fields(this%getVarIdx('pres'))%DATA

    ASSOCIATE(sg => this%sg)
    ! Calculate the pressure at the model level with surface pressure and hydrosatic balance
    DO j = 2, sg%vLevel
      DO k = 1, sg%tSlots
      swapDA(j, :, k) = swapDA(j - 1, :, k) * EXP(-(sg%zHght(j, :) - sg%zHght(j - 1, :)) * g &
                                                                        / (dry_air_gas_const * &
                                                                           (this%Fields(this%getVarIdx('temp'))%DATA(j, :, k) * &
                                                                            (this%Fields(this%getVarIdx('qvapor'))%DATA(j, :, k) * 0.608 + 1) + &
                                                                            this%Fields(this%getVarIdx('temp'))%DATA(j - 1, :, k) * &
                                                                            (this%Fields(this%getVarIdx('qvapor'))%DATA(j - 1, :, k) * 0.608 + 1)) / 2.0D0))
      END DO
    END DO
    this%Fields(this%getVarIdx('pres'))%DATA = swapDA

    ! Calculate the Saturated vapor at each level
    ip = this%getVarIdx('pres')
    iq = this%getVarIdx('qvapor')
    it = this%getVarIdx('temp')
    ! Note: without these precalculated indices, some times the compiler complains about the getVarIdx Yuanfu Xie 2025-05-19
    ASSOCIATE (temp => this%Fields(it)%DATA, pres => this%Fields(ip)%DATA, &
               qvapor => this%fields(iq)%DATA)
      DO j = 1, sg%vLevel
        DO i = 1, sg%num_cell
          DO k = 1, sg%tSlots
          IF (qvapor(j, i, k) > Saturated_Vapor(temp(j, i, k), pres(j, i, k))) THEN
            qvapor(j, i, k) = Saturated_Vapor(temp(j, i, k), pres(j, i, k))
          END IF
          END DO
        END DO
      END DO
    END ASSOCIATE

    ! Recalculate the pressure at the model level with surface pressure and hydrosatic balance
    DO j = 2, sg%vLevel
      DO k = 1, sg%tSlots
      swapDA(j, :, k) = swapDA(j - 1, :, k) * EXP(-(sg%zHght(j, :) - sg%zHght(j - 1, :)) * g &
                                                                        / (dry_air_gas_const * &
                                                                           (this%Fields(this%getVarIdx('temp'))%DATA(j, :, k) * &
                                                                            (this%Fields(this%getVarIdx('qvapor'))%DATA(j, :, k) * 0.608 + 1) + &
                                                                            this%Fields(this%getVarIdx('temp'))%DATA(j - 1, :, k) * &
                                                                            (this%Fields(this%getVarIdx('qvapor'))%DATA(j - 1, :, k) * 0.608 + 1)) / 2.0D0))
      END DO
    END DO
    this%Fields(this%getVarIdx('pres'))%DATA = swapDA
    
    DEALLOCATE (swapDA)
    END ASSOCIATE

  END SUBROUTINE

END MODULE State_m
