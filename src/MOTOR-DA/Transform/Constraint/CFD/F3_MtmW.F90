!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2022/3/20, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the F3_MtmW Object which is tranform the variables from control variable space to
!! model space.
MODULE F3_MtmW_m
  USE State_m, ONLY: State_t
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE M2MBase_m, ONLY: M2MBase_t

  TYPE, EXTENDS(M2MBase_t) :: F3_MtmW_t
    TYPE(State_t), POINTER :: X

  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: transFwdNonLinear_opr
    PROCEDURE, PUBLIC, PASS(this) :: transFwdTanLinear_opr
    PROCEDURE, PUBLIC, PASS(this) :: transAdjMultiply_opr

    PROCEDURE, PUBLIC, PASS(this) :: transFwdNonLinear
    PROCEDURE, PUBLIC, PASS(this) :: transFwdTanLinear
    PROCEDURE, PUBLIC, PASS(this) :: transAdjMultiply

    PROCEDURE :: fwdNL_opr => transFwdNonLinear_opr
    PROCEDURE :: fwdTL_opr => transFwdTanLinear_opr
    PROCEDURE :: adjMul_opr => transAdjMultiply_opr

    PROCEDURE :: fwdNL => transFwdNonLinear
    PROCEDURE :: fwdTL => transFwdTanLinear
    PROCEDURE :: adjMul => transAdjMultiply

    PROCEDURE, PUBLIC, NOPASS :: transElemForward
    PROCEDURE, PUBLIC, NOPASS :: transElemAdjointMultiply
    PROCEDURE :: calParameters

    FINAL :: destructor
  END TYPE F3_MtmW_t

  INTERFACE F3_MtmW_t
    PROCEDURE :: constructor
  END INTERFACE F3_MtmW_t

CONTAINS

  SUBROUTINE calParameters(this)
    IMPLICIT NONE
    CLASS(F3_MtmW_t) :: this
    INTEGER(i_kind) :: i, j, k
    REAL(r_kind) :: t1, t2, t3
    REAL(r_kind), ALLOCATABLE :: psMatrix(:, :)

    ASSOCIATE (sg => this%X%sg)

      ALLOCATE (this%c1(this%X%sg%vLevel, this%X%sg%num_cell), this%c2(this%X%sg%vLevel, this%X%sg%num_cell), &
                this%c3(this%X%sg%num_cell))

      ! Calculate the C term
      BLOCK
        REAL(r_kind), ALLOCATABLE :: pZspx(:), &
                                     pZspy(:)

        ALLOCATE (pZspx(sg%num_cell), pZspy(sg%num_cell))

        pZspx = 0.0D0
        pZspy = 0.0D0

        DO i = 1, sg%num_icell
          IF ((sg%cell_stcl(4, i) .NE. 0) .AND. (sg%cell_stcl(6, i) .NE. 0)) THEN
            pZspx(i) = (sg%zHght(1, sg%cell_stcl(6, i)) - sg%zHght(1, sg%cell_stcl(4, i))) / (2 * sg%cell_dist(4, i))
          END IF

          IF ((sg%cell_stcl(2, i) .NE. 0) .AND. (sg%cell_stcl(8, i) .NE. 0)) THEN
            pZspy(i) = (sg%zHght(1, sg%cell_stcl(8, i)) - sg%zHght(1, sg%cell_stcl(2, i))) / (2 * sg%cell_dist(2, i))
          END IF
        END DO

        this%c1 = 0.0D0
        this%c2 = 0.0D0
        this%c3 = 0.0D0

        DO i = 1, sg%vLevel
          this%c1(i, :) = sg%zTop * (sg%ztop - sg%zHght(i, :)) / ((sg%ztop - sg%zHght(1, :))**2) * pZspx
          this%c2(i, :) = sg%zTop * (sg%ztop - sg%zHght(i, :)) / ((sg%ztop - sg%zHght(1, :))**2) * pZspy
          this%c3 = sg%zTop / (sg%zTop - sg%zHght(1, :))
        END DO
        ! CALL sg%ExchangeMatOnHalo2D(sg%vLevel, this%C)! Exchange hale

        DEALLOCATE (pZspx, pZspy)

      END BLOCK

      ! Calculate the stencils and coefficients
      BLOCK
        REAL(r_kind) :: h1, h2, h3, a1, a2, a3

        REAL(r_kind), ALLOCATABLE :: coef_p_pX(:, :)      ! 3 stcl at horizontal / num_icell
        REAL(r_kind), ALLOCATABLE :: coef_p_pY(:, :)      ! 3 stcl at horizontal / num_icell
        REAL(r_kind), ALLOCATABLE :: coef_p_pSigma(:, :)  ! 3 stcl at verticl/ 3 stcl at horizontal / vLevel / num_icell

        INTEGER(i_kind), ALLOCATABLE :: stcl_p_pX(:, :)       ! 3 stcl / num_icell
        INTEGER(i_kind), ALLOCATABLE :: stcl_p_pY(:, :)       ! 3 stcl / num_icell
        INTEGER(i_kind), ALLOCATABLE :: stcl_p_pSigma(:, :)   ! 3 stcl / vLevel

        ALLOCATE (this%coef_p_pX(3, sg%num_cell), &      ! 3 stcl at verticl/ 3 stcl at horizontal/ num_icell
                  this%coef_p_pY(3, sg%num_cell), &      ! 3 stcl at verticl/ 3 stcl at horizontal/ num_icell
                  this%coef_p_pSigma(3, sg%vLevel), &    ! 3 stcl / num_icell
                  this%stcl_p_pX(3, sg%num_cell), &    ! 3 stcl / num_icell
                  this%stcl_p_pY(3, sg%num_cell), &
                  this%stcl_p_pSigma(3, sg%vLevel))   ! 3 stcl / vLevel

        this%coef_p_pX = 0.0D0
        this%coef_p_pY = 0.0D0
        this%coef_p_pSigma = 0.0D0
        this%stcl_p_pX = 0
        this%stcl_p_pY = 0
        this%stcl_p_pSigma = 0

        DO i = 1, sg%num_icell
          IF ((sg%cell_stcl(4, i) .NE. 0) .AND. (sg%cell_stcl(6, i) .NE. 0)) THEN
            t1 = -1.0D0 / (2.0D0 * sg%cell_dist(4, i))
            t2 = 0.0D0
            t3 = 1.0D0 / (2.0D0 * sg%cell_dist(4, i))

            this%stcl_p_pX(1, i) = sg%cell_stcl(4, i)
            this%coef_p_pX(1, i) = t1
            this%stcl_p_pX(2, i) = sg%cell_stcl(5, i)
            this%coef_p_pX(2, i) = t2
            this%stcl_p_pX(3, i) = sg%cell_stcl(6, i)
            this%coef_p_pX(3, i) = t3
          END IF

          IF ((sg%cell_stcl(2, i) .NE. 0) .AND. (sg%cell_stcl(8, i) .NE. 0)) THEN

            t1 = -1.0D0 / (2.0D0 * sg%cell_dist(2, i))
            t2 = 0.0D0
            t3 = 1.0D0 / (2.0D0 * sg%cell_dist(2, i))

            this%stcl_p_pY(1, i) = sg%cell_stcl(2, i)
            this%coef_p_pY(1, i) = t1
            this%stcl_p_pY(2, i) = sg%cell_stcl(5, i)
            this%coef_p_pY(2, i) = t2
            this%stcl_p_pY(3, i) = sg%cell_stcl(8, i)
            this%coef_p_pY(3, i) = t3
          END IF
        END DO

        DO i = 2, sg%vLevel - 1
          h1 = sg%sigma(i - 1) - sg%sigma(i)
          h2 = 0
          h3 = sg%sigma(i + 1) - sg%sigma(i)

          this%stcl_p_pSigma(1, i) = i - 1
          this%stcl_p_pSigma(2, i) = i
          this%stcl_p_pSigma(3, i) = i + 1

          a1 = -(h2 + h3) / (h1 - h2) / (h1 - h3)
          a2 = (h1 + h3) / (h1 - h2) / (h2 - h3)
          a3 = (h1 + h2) / (h1 - h3) / (h3 - h2)

          this%coef_p_pSigma(1, i) = a1
          this%coef_p_pSigma(2, i) = a2
          this%coef_p_pSigma(3, i) = a3

        END DO
      END BLOCK
    END ASSOCIATE
  END SUBROUTINE

  FUNCTION constructor(configFile, X) RESULT(this)
    IMPLICIT NONE
    TYPE(F3_MtmW_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X

    this%X => X
    IF (X%sg%vLevel < 3) RETURN
    CALL this%calParameters
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(F3_MtmW_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%coef_p_pX)) DEALLOCATE (this%coef_p_pX)
    IF (ALLOCATED(this%coef_p_pY)) DEALLOCATE (this%coef_p_pY)
    IF (ALLOCATED(this%coef_p_pSigma)) DEALLOCATE (this%coef_p_pSigma)

    IF (ALLOCATED(this%stcl_p_pX)) DEALLOCATE (this%stcl_p_pX)
    IF (ALLOCATED(this%stcl_p_pY)) DEALLOCATE (this%stcl_p_pY)
    IF (ALLOCATED(this%stcl_p_pSigma)) DEALLOCATE (this%stcl_p_pSigma)

    IF (ALLOCATED(this%c1)) DEALLOCATE (this%c1)
    IF (ALLOCATED(this%c2)) DEALLOCATE (this%c2)
    IF (ALLOCATED(this%c3)) DEALLOCATE (this%c3)
  END SUBROUTINE destructor

  SUBROUTINE transFwdNonLinear(this, X)
    IMPLICIT NONE
    CLASS(F3_MtmW_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER :: i, j, k, hStcl, vStcl

    IF (X%sg%vLevel < 3) RETURN

    IF ((X%getVarIdx(TRIM('uwnd')) .EQ. 0) .OR. &
        (X%getVarIdx(TRIM('vwnd')) .EQ. 0) .OR. &
        (X%getVarIdx(TRIM('wwnd')) .EQ. 0)) THEN
      PRINT *, 'Error use the transFwdNonLinear in F3_MtmW_t, STOP! '
      STOP
    END IF

    CALL X%AddVar('divg')

    ASSOCIATE (uwnd => X%Fields(X%getVarIdx(TRIM('uwnd')))%DATA, &
               vwnd => X%Fields(X%getVarIdx(TRIM('vwnd')))%DATA, &
               wwnd => X%Fields(X%getVarIdx(TRIM('wwnd')))%DATA, &
               divg => X%Fields(X%getVarIdx(TRIM('divg')))%DATA)

      BLOCK
        divg = 0.0D0
        DO j = 1, X%sg%num_icell
          DO i = 2, X%sg%vLevel - 1

            IF (X%sg%cell_type(j) /= 0) CYCLE
            !First term
            DO hStcl = 1, 3
              divg(i, j, :) = divg(i, j, :) + &
                              uwnd(i, this%stcl_p_pX(hStcl, j), :) * this%coef_p_pX(hStcl, j) + &
                              vwnd(i, this%stcl_p_pY(hStcl, j), :) * this%coef_p_pY(hStcl, j)
            END DO

            ! Second term
            DO vStcl = 1, 3
              divg(i, j, :) = divg(i, j, :) &
                              - 1.0D0 * this%c1(i, j) * uwnd(this%stcl_p_pSigma(vStcl, i), j, :) * this%coef_p_pSigma(vStcl, i) &
                              - 1.0D0 * this%c2(i, j) * vwnd(this%stcl_p_pSigma(vStcl, i), j, :) * this%coef_p_pSigma(vStcl, i) + &
                              this%c3(i) * wwnd(this%stcl_p_pSigma(vStcl, i), j, :) * this%coef_p_pSigma(vStcl, i)
            END DO
          END DO
        END DO

      END BLOCK
    END ASSOCIATE

    ! CALL X%exHalo()
    CALL X%sg%ExchangeMatOnHaloForFieldGrid(X%sg%tSlots, X%sg%vLevel, X%Fields(X%getVarIdx(TRIM('divg')))%DATA)

  END SUBROUTINE

  SUBROUTINE transAdjMultiply(this, dX, X)
    IMPLICIT NONE
    CLASS(F3_MtmW_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(State_t) :: dXi
    INTEGER(i_kind)  :: i, j, hStcl, vStcl, k

    IF (dX%sg%vLevel < 3) RETURN

    IF ((dX%getVarIdx(TRIM('uwnd')) .EQ. 0) .OR. &
        (dX%getVarIdx(TRIM('vwnd')) .EQ. 0) .OR. &
        (dX%getVarIdx(TRIM('wwnd')) .EQ. 0) .OR. &
        (dX%getVarIdx(TRIM('divg')) .EQ. 0)) THEN
      PRINT *, 'Error use the transFwdNonLinear in F3_MtmW_t, STOP! '
      STOP
    END IF

    ! Construct an empty state to contain increment
    dXi = dX%zeroCopy()
    dXi%fields(dXi%getVarIdx('divg')) = dX%fields(dX%getVarIdx('divg'))

    ASSOCIATE (uwnd => dXi%Fields(dXi%getVarIdx(TRIM('uwnd')))%DATA, &
               vwnd => dXi%Fields(dXi%getVarIdx(TRIM('vwnd')))%DATA, &
               wwnd => dXi%Fields(dXi%getVarIdx(TRIM('wwnd')))%DATA, &
               divg => dXi%Fields(dXi%getVarIdx(TRIM('divg')))%DATA)
      ! PRINT*, 'Going here 4'

      BLOCK
        DO j = 1, dXi%sg%num_icell
          DO i = 2, dXi%sg%vLevel - 1
            IF (dX%sg%cell_type(j) /= 0) CYCLE

            ! First term
            DO hStcl = 1, 3
              uwnd(i, this%stcl_p_pX(hStcl, j), :) = uwnd(i, this%stcl_p_pX(hStcl, j), :) + &
                                                     divg(i, j, :) * this%coef_p_pX(hStcl, j)
              vwnd(i, this%stcl_p_pY(hStcl, j), :) = vwnd(i, this%stcl_p_pY(hStcl, j), :) + &
                                                     divg(i, j, :) * this%coef_p_pY(hStcl, j)
            END DO

            ! Second term
            DO vStcl = 1, 3
              uwnd(this%stcl_p_pSigma(vStcl, i), j, :) = uwnd(this%stcl_p_pSigma(vStcl, i), j, :) &
                                                         - this%c1(i, j) * divg(i, j, :) * this%coef_p_pSigma(vStcl, i)
              vwnd(this%stcl_p_pSigma(vStcl, i), j, :) = vwnd(this%stcl_p_pSigma(vStcl, i), j, :) &
                                                         - this%c2(i, j) * divg(i, j, :) * this%coef_p_pSigma(vStcl, i)
              wwnd(this%stcl_p_pSigma(vStcl, i), j, :) = wwnd(this%stcl_p_pSigma(vStcl, i), j, :) &
                                                         + this%c3(j) * divg(i, j, :) * this%coef_p_pSigma(vStcl, i)
            END DO
          END DO
        END DO
        wwnd(1, :, :) = 0.0D0
        IF (dXi%sg%ztop > 25000.0) wwnd(dXi%sg%vLevel, :, :) = 0.0D0

        ! ===================================
        DO j = 1, dXi%sg%num_cell
          DO k = 1, dXi%sg%vLevel
            IF (dXi%sg%zHght(k, j) < dXi%sg%topo(j)) THEN
              uwnd(k, j, :) = 0.0D0
              vwnd(k, j, :) = 0.0D0
              wwnd(k, j, :) = 0.0D0
            END IF
          END DO
        END DO

        ! wwnd = 0.0D0
!=================================
      END BLOCK
    END ASSOCIATE

    CALL dXi%rmVar('divg')
    CALL dXi%exHaloRevSum()
    CALL dXi%exHalo()

    ! Add the increment with dX
    CALL dx%rmVar('divg')
    dX = dX + dXi
  END SUBROUTINE transAdjMultiply

  SUBROUTINE transFwdTanLinear(this, dX, X)
    CLASS(F3_MtmW_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    CALL this%transFwdNonLinear(dX)
  END SUBROUTINE

  FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
    IMPLICIT NONE
    CLASS(F3_MtmW_t) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t) :: dX1
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    dX1 = dX
    CALL this%transAdjMultiply(dX1)
  END FUNCTION

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(F3_MtmW_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1
    INTEGER(i_kind) :: i, j, k

    X1 = X
    CALL this%transFwdNonLinear(X1)

    DO i = 1, SIZE(X1%fields)
      IF (X1%fields(i)%Get_Name() /= 'divg') THEN
        X1%fields(i)%DATA = 0.0D0
      END IF
    END DO
  END FUNCTION

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
    CLASS(F3_MtmW_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(State_t):: dX1

    dX1 = this%transFwdNonLinear_opr(dX)
  END FUNCTION

  FUNCTION transElemForward(valIn, s1) RESULT(valOut)
    REAL(r_kind) :: valIn, valOut, s1

  END FUNCTION

  FUNCTION transElemAdjointMultiply(valIn, s1) RESULT(valOut)
    REAL(r_kind) :: valIn, valOut, s1

  END FUNCTION
END MODULE F3_MtmW_m

