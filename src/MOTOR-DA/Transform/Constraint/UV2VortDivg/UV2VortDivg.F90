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
!! This module contains the UV2VortDivg Object which is tranform the variables from control variable space to
!! model space.
MODULE UV2VortDivg_m
  USE State_m, ONLY: State_t
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE M2MBase_m, ONLY: M2MBase_t

  TYPE, EXTENDS(M2MBase_t) :: UV2VortDivg_t
    TYPE(State_t), POINTER :: X

    ! First term
    REAL(r_kind), ALLOCATABLE :: coef_u_1st(:, :, :, :)      ! 3 stcl at verticl/ 3 stcl at horizontal/ num_icell
    REAL(r_kind), ALLOCATABLE :: coef_v_1st(:, :, :, :)      ! 3 stcl at verticl/ 3 stcl at horizontal/ num_icell
    INTEGER(i_kind), ALLOCATABLE :: stcl_u_1st_h(:, :)    ! 3 stcl / num_icell
    INTEGER(i_kind), ALLOCATABLE :: stcl_v_1st_h(:, :)    ! 3 stcl / num_icell
    INTEGER(i_kind), ALLOCATABLE :: stcl_uv_1st_v(:, :)   ! 3 stcl / vLevel

    ! Second term
    REAL(r_kind), ALLOCATABLE :: coef_u_2nd(:, :, :)         ! 4 stcl at verticl/ vLevel
    REAL(r_kind), ALLOCATABLE :: coef_v_2nd(:, :, :)         ! 4 stcl at verticl/ vLevel
    INTEGER(i_kind), ALLOCATABLE :: stcl_uv_2nd_v(:, :)    ! 4 stcl /vLevel

    REAL(r_kind), ALLOCATABLE :: qrMatrix(:, :)
    REAL(r_kind), ALLOCATABLE :: tauMatrix(:, :)

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

    ! GENERIC :: OPERATOR(.yaml.) => transFwdNonLinear_opr
    ! GENERIC :: OPERATOR(.ADJ.) => transAdjMultiply_opr
    FINAL :: destructor
  END TYPE UV2VortDivg_t

  INTERFACE UV2VortDivg_t
    PROCEDURE :: constructor
  END INTERFACE UV2VortDivg_t

CONTAINS

  SUBROUTINE calParameters(this)
    IMPLICIT NONE
    CLASS(UV2VortDivg_t) :: this
    INTEGER(i_kind) :: i, j, k
    REAL(r_kind) :: t1, t2, t3
    REAL(r_kind), ALLOCATABLE :: psMatrix(:, :)
    REAL(r_kind), ALLOCATABLE :: c1(:, :), c2(:, :), c3(:)

    ASSOCIATE (sg => this%X%sg)

      ALLOCATE (psMatrix(this%X%sg%vLevel, this%X%sg%vLevel))
      ALLOCATE (c1(this%X%sg%vLevel, this%X%sg%num_cell), c2(this%X%sg%vLevel, this%X%sg%num_cell), c3(this%X%sg%num_cell))

      ! Calculate the C term
      BLOCK
        REAL(r_kind), ALLOCATABLE :: pZspx(:), &
                                     pZspy(:)

        ALLOCATE (pZspx(sg%num_cell), pZspy(sg%num_cell))

        pZspx = 0.0D0
        pZspy = 0.0D0

        DO i = 1, sg%num_icell
          IF ((sg%cell_stcl(4, i) .NE. 0) .AND. (sg%cell_stcl(6, i) .NE. 0)) THEN
            pZspx(i) = (sg%topo(sg%cell_stcl(6, i)) - sg%topo(sg%cell_stcl(4, i))) / (2 * sg%cell_dist(4, i))
          END IF

          IF ((sg%cell_stcl(4, i) .EQ. 0)) THEN
            t1 = -3.0D0 / 2.0D0 / sg%cell_dist(6, i)
            t2 = 2.0D0 / sg%cell_dist(6, i)
            t3 = -1.0D0 / 2.0D0 / sg%cell_dist(6, i)
            pZspx(i) = sg%topo(sg%cell_stcl(5, i)) * t1 + &
                       sg%topo(sg%cell_stcl(6, i)) * t2 + &
                       sg%topo(sg%cell_stcl(6, sg%cell_stcl(6, i))) * t3
          END IF

          IF ((sg%cell_stcl(6, i) .EQ. 0)) THEN
            t1 = 1.0D0 / 2.0D0 / sg%cell_dist(4, i)
            t2 = -2.0D0 / sg%cell_dist(4, i)
            t3 = 3.0D0 / 2.0D0 / sg%cell_dist(4, i)
            pZspx(i) = sg%topo(sg%cell_stcl(4, sg%cell_stcl(4, i))) * t1 + &
                       sg%topo(sg%cell_stcl(4, i)) * t2 + &
                       sg%topo(sg%cell_stcl(5, i)) * t3
          END IF

          IF ((sg%cell_stcl(2, i) .NE. 0) .AND. (sg%cell_stcl(8, i) .NE. 0)) THEN
            pZspy(i) = (sg%topo(sg%cell_stcl(8, i)) - sg%topo(sg%cell_stcl(2, i))) / (2 * sg%cell_dist(2, i))
          END IF

          IF ((sg%cell_stcl(2, i) .EQ. 0)) THEN
            t1 = -3.0D0 / 2.0D0 / sg%cell_dist(8, i)
            t2 = 2.0D0 / sg%cell_dist(8, i)
            t3 = -1.0D0 / 2.0D0 / sg%cell_dist(8, i)
            pZspy(i) = sg%topo(sg%cell_stcl(5, i)) * t1 + &
                       sg%topo(sg%cell_stcl(8, i)) * t2 + &
                       sg%topo(sg%cell_stcl(8, sg%cell_stcl(8, i))) * t3
          END IF

          IF ((sg%cell_stcl(8, i) .EQ. 0)) THEN
            t1 = 1.0D0 / 2.0D0 / sg%cell_dist(2, i)
            t2 = -2.0D0 / sg%cell_dist(2, i)
            t3 = 3.0D0 / 2.0D0 / sg%cell_dist(2, i)
            pZspy(i) = sg%topo(sg%cell_stcl(2, sg%cell_stcl(2, i))) * t1 + &
                       sg%topo(sg%cell_stcl(2, i)) * t2 + &
                       sg%topo(sg%cell_stcl(5, i)) * t3
          END IF
        END DO

        DO i = 1, sg%vLevel
          c1(i, :) = sg%zTop * (sg%ztop - sg%zHght(i, :)) / ((sg%ztop - sg%topo)**2) * pZspx
          c2(i, :) = sg%zTop * (sg%ztop - sg%zHght(i, :)) / ((sg%ztop - sg%topo)**2) * pZspy
          c3 = sg%zTop / (sg%zTop - sg%topo)
        END DO
        ! CALL sg%ExchangeMatOnHalo2D(sg%vLevel, this%C)! Exchange hale

        DEALLOCATE (pZspx, pZspy)

      END BLOCK

      ! Calculate the psMatrix
      BLOCK
        ! Middle points, use i-1, i, i+1 and i+2
        REAL(r_kind) :: h1, h2, h3, h4
        INTEGER(i_kind) :: LWORK, INFO
        REAL(r_kind), ALLOCATABLE :: ceof_2nd_v(:, :), WORK(:)

        ALLOCATE (WORK(sg%vLevel * sg%vLevel), ceof_2nd_v(4, sg%vLevel))

        ALLOCATE (this%coef_u_2nd(4, sg%vLevel, sg%num_cell), &         ! 4 stcl at verticl/ vLevel
                  this%coef_v_2nd(4, sg%vLevel, sg%num_cell), &         ! 4 stcl at verticl/ vLevel
                  this%stcl_uv_2nd_v(4, sg%vLevel))    ! 4 stcl /vLevel

        this%coef_u_2nd = 0.0D0
        this%coef_v_2nd = 0.0D0
        this%stcl_uv_2nd_v = 0.0D0

        LWORK = sg%vLevel * sg%vLevel
        psMatrix = 0.0D0

        psMatrix(1, 1) = 1.0D0
        psMatrix(sg%vLevel, sg%vLevel) = 1.0D0

        DO i = 2, sg%vLevel - 2
          h1 = sg%sigma(i - 1) - sg%sigma(i)
          h2 = 0.0D0
          h3 = sg%sigma(i + 1) - sg%sigma(i)
          h4 = sg%sigma(i + 2) - sg%sigma(i)

          psMatrix(i, i - 1) = (-2 * (h2 + h3 + h4) / ((h1 - h2) * (h1 - h3) * (h1 - h4)))
          psMatrix(i, i) = (2 * (h1 + h3 + h4) / ((h1 - h2) * (h2 - h3) * (h2 - h4)))
          psMatrix(i, i + 1) = (2 * (h1 + h2 + h4) / ((h1 - h3) * (h3 - h2) * (h3 - h4)))
          psMatrix(i, i + 2) = (2 * (h1 + h2 + h3) / ((h1 - h4) * (h4 - h2) * (h4 - h3)))

          ceof_2nd_v(:, i) = psMatrix(i, i - 1:i + 2)
          this%stcl_uv_2nd_v(:, i) = (/i - 1, i, i + 1, i + 2/)
        END DO

        i = sg%vLevel - 1
        h1 = sg%sigma(i - 2) - sg%sigma(i)
        h2 = sg%sigma(i - 1) - sg%sigma(i)
        h3 = 0.0D0
        h4 = sg%sigma(i + 1) - sg%sigma(i)

        psMatrix(i, i - 2) = (-2 * (h2 + h3 + h4) / ((h1 - h2) * (h1 - h3) * (h1 - h4)))
        psMatrix(i, i - 1) = (2 * (h1 + h3 + h4) / ((h1 - h2) * (h2 - h3) * (h2 - h4)))
        psMatrix(i, i) = (2 * (h1 + h2 + h4) / ((h1 - h3) * (h3 - h2) * (h3 - h4)))
        psMatrix(i, i + 1) = (2 * (h1 + h2 + h3) / ((h1 - h4) * (h4 - h2) * (h4 - h3)))

        ceof_2nd_v(:, i) = psMatrix(i, i - 2:i + 1)
        this%stcl_uv_2nd_v(:, i) = (/i - 2, i - 1, i, i + 1/)

        this%qrMatrix = psMatrix

        ! QR Decomposition
        ALLOCATE (this%tauMatrix(sg%vLevel, sg%vLevel))
        CALL DGEQRF(sg%vLevel, sg%vLevel, this%qrMatrix, sg%vLevel, this%tauMatrix, WORK, LWORK, INFO)

        DO i = 1, sg%num_icell
          DO j = 1, 4
            this%coef_u_2nd(j, :, i) = -1.0D0 * ceof_2nd_v(j, :) * c1(:, i) / c3(i)
            this%coef_v_2nd(j, :, i) = -1.0D0 * ceof_2nd_v(j, :) * c2(:, i) / c3(i)
          END DO
        END DO

        DEALLOCATE (WORK, ceof_2nd_v)
      END BLOCK

      ! Calculate the stencils and coefficients
      BLOCK
        REAL(r_kind) :: h1, h2, h3, a1, a2, a3

        ALLOCATE (this%coef_u_1st(3, 3, sg%vLevel, sg%num_cell), &      ! 3 stcl at verticl/ 3 stcl at horizontal/ num_icell
                  this%coef_v_1st(3, 3, sg%vLevel, sg%num_cell), &      ! 3 stcl at verticl/ 3 stcl at horizontal/ num_icell
                  this%stcl_u_1st_h(3, sg%num_cell), &    ! 3 stcl / num_icell
                  this%stcl_v_1st_h(3, sg%num_cell), &    ! 3 stcl / num_icell
                  this%stcl_uv_1st_v(3, sg%vLevel))   ! 3 stcl / vLevel

        this%coef_u_1st = 0.0D0
        this%coef_v_1st = 0.0D0
        this%stcl_u_1st_h = 0
        this%stcl_v_1st_h = 0
        this%stcl_uv_1st_v = 0

        DO i = 1, sg%num_icell
          IF ((sg%cell_stcl(4, i) .NE. 0) .AND. (sg%cell_stcl(6, i) .NE. 0)) THEN
            t1 = -1.0D0 / (2.0D0 * sg%cell_dist(4, i))
            t2 = 0.0D0
            t3 = 1.0D0 / (2.0D0 * sg%cell_dist(4, i))

            this%stcl_u_1st_h(1, i) = sg%cell_stcl(4, i)
            this%coef_u_1st(:, 1, :, i) = t1
            this%stcl_u_1st_h(2, i) = sg%cell_stcl(5, i)
            this%coef_u_1st(:, 2, :, i) = t2
            this%stcl_u_1st_h(3, i) = sg%cell_stcl(6, i)
            this%coef_u_1st(:, 3, :, i) = t3
          END IF

          IF ((sg%cell_stcl(4, i) .EQ. 0)) THEN
            t1 = -3.0D0 / 2.0D0 / sg%cell_dist(6, i)
            t2 = 2.0D0 / sg%cell_dist(6, i)
            t3 = -1.0D0 / 2.0D0 / sg%cell_dist(6, i)

            this%stcl_u_1st_h(1, i) = sg%cell_stcl(5, i)
            this%coef_u_1st(:, 1, :, i) = t1
            this%stcl_u_1st_h(2, i) = sg%cell_stcl(6, i)
            this%coef_u_1st(:, 2, :, i) = t2
            this%stcl_u_1st_h(3, i) = sg%cell_stcl(6, sg%cell_stcl(6, i))
            this%coef_u_1st(:, 3, :, i) = t3
          END IF

          IF ((sg%cell_stcl(6, i) .EQ. 0)) THEN
            t1 = 1.0D0 / 2.0D0 / sg%cell_dist(4, i)
            t2 = -2.0D0 / sg%cell_dist(4, i)
            t3 = 3.0D0 / 2.0D0 / sg%cell_dist(4, i)

            this%stcl_u_1st_h(1, i) = sg%cell_stcl(4, sg%cell_stcl(4, i))
            this%coef_u_1st(:, 1, :, i) = t1
            this%stcl_u_1st_h(2, i) = sg%cell_stcl(4, i)
            this%coef_u_1st(:, 2, :, i) = t2
            this%stcl_u_1st_h(3, i) = sg%cell_stcl(5, i)
            this%coef_u_1st(:, 3, :, i) = t3
          END IF

          IF ((sg%cell_stcl(2, i) .NE. 0) .AND. (sg%cell_stcl(8, i) .NE. 0)) THEN

            t1 = -1.0D0 / (2.0D0 * sg%cell_dist(2, i))
            t2 = 0.0D0
            t3 = 1.0D0 / (2.0D0 * sg%cell_dist(2, i))

            this%stcl_v_1st_h(1, i) = sg%cell_stcl(2, i)
            this%coef_v_1st(:, 1, :, i) = t1
            this%stcl_v_1st_h(2, i) = sg%cell_stcl(5, i)
            this%coef_v_1st(:, 2, :, i) = t2
            this%stcl_v_1st_h(3, i) = sg%cell_stcl(8, i)
            this%coef_v_1st(:, 3, :, i) = t3
          END IF

          IF ((sg%cell_stcl(2, i) .EQ. 0)) THEN
            t1 = -3.0D0 / 2.0D0 / sg%cell_dist(8, i)
            t2 = 2.0D0 / sg%cell_dist(8, i)
            t3 = -1.0D0 / 2.0D0 / sg%cell_dist(8, i)

            this%stcl_v_1st_h(1, i) = sg%cell_stcl(5, i)
            this%coef_v_1st(:, 1, :, i) = t1
            this%stcl_v_1st_h(2, i) = sg%cell_stcl(8, i)
            this%coef_v_1st(:, 2, :, i) = t2
            this%stcl_v_1st_h(3, i) = sg%cell_stcl(8, sg%cell_stcl(8, i))
            this%coef_v_1st(:, 3, :, i) = t3
          END IF

          IF ((sg%cell_stcl(8, i) .EQ. 0)) THEN
            t1 = 1.0D0 / 2.0D0 / sg%cell_dist(2, i)
            t2 = -2.0D0 / sg%cell_dist(2, i)
            t3 = 3.0D0 / 2.0D0 / sg%cell_dist(2, i)

            this%stcl_v_1st_h(1, i) = sg%cell_stcl(2, sg%cell_stcl(2, i))
            this%coef_v_1st(:, 1, :, i) = t1
            this%stcl_v_1st_h(2, i) = sg%cell_stcl(2, i)
            this%coef_v_1st(:, 2, :, i) = t2
            this%stcl_v_1st_h(3, i) = sg%cell_stcl(5, i)
            this%coef_v_1st(:, 3, :, i) = t3
          END IF
        END DO

        DO i = 2, sg%vLevel - 1
          h1 = sg%sigma(i - 1) - sg%sigma(i)
          h2 = 0
          h3 = sg%sigma(i + 1) - sg%sigma(i)

          this%stcl_uv_1st_v(1, i) = i - 1
          this%stcl_uv_1st_v(2, i) = i
          this%stcl_uv_1st_v(3, i) = i + 1

          a1 = -(h2 + h3) / (h1 - h2) / (h1 - h3)
          a2 = (h1 + h3) / (h1 - h2) / (h2 - h3)
          a3 = (h1 + h2) / (h1 - h3) / (h3 - h2)

          this%coef_u_1st(1, :, i, :) = this%coef_u_1st(1, :, i, :) * a1
          this%coef_u_1st(2, :, i, :) = this%coef_u_1st(2, :, i, :) * a2
          this%coef_u_1st(3, :, i, :) = this%coef_u_1st(3, :, i, :) * a3

          this%coef_v_1st(1, :, i, :) = this%coef_v_1st(1, :, i, :) * a1
          this%coef_v_1st(2, :, i, :) = this%coef_v_1st(2, :, i, :) * a2
          this%coef_v_1st(3, :, i, :) = this%coef_v_1st(3, :, i, :) * a3
        END DO
      END BLOCK

      DO i = 1, 3
        DO j = 1, 3
          DO k = 1, sg%vLevel
            this%coef_u_1st(i, j, k, :) = this%coef_u_1st(i, j, k, :) / c3
            this%coef_v_1st(i, j, k, :) = this%coef_v_1st(i, j, k, :) / c3
          END DO
        END DO
      END DO

      DEALLOCATE (psMatrix, c1, c2, c3)

    END ASSOCIATE
  END SUBROUTINE

  FUNCTION constructor(configFile, X) RESULT(this)
    IMPLICIT NONE
    TYPE(UV2VortDivg_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X

    this%X => X
    IF (X%sg%vLevel < 3) RETURN
    CALL this%calParameters
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(UV2VortDivg_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%qrMatrix)) DEALLOCATE (this%qrMatrix)
    IF (ALLOCATED(this%tauMatrix)) DEALLOCATE (this%tauMatrix)

    IF (ALLOCATED(this%coef_u_1st)) DEALLOCATE (this%coef_u_1st)
    IF (ALLOCATED(this%coef_v_1st)) DEALLOCATE (this%coef_v_1st)
    IF (ALLOCATED(this%stcl_u_1st_h)) DEALLOCATE (this%stcl_u_1st_h)
    IF (ALLOCATED(this%stcl_v_1st_h)) DEALLOCATE (this%stcl_v_1st_h)
    IF (ALLOCATED(this%stcl_uv_1st_v)) DEALLOCATE (this%stcl_uv_1st_v)

    IF (ALLOCATED(this%coef_u_2nd)) DEALLOCATE (this%coef_u_2nd)
    IF (ALLOCATED(this%coef_v_2nd)) DEALLOCATE (this%coef_v_2nd)
    IF (ALLOCATED(this%stcl_uv_2nd_v)) DEALLOCATE (this%stcl_uv_2nd_v)

  END SUBROUTINE destructor

  SUBROUTINE transFwdNonLinear(this, X)
    IMPLICIT NONE
    CLASS(UV2VortDivg_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER :: i, j, k, hStcl, vStcl

    IF (X%sg%vLevel < 3) RETURN

    IF ((X%getVarIdx(TRIM('uwnd')) .EQ. 0) .OR. &
        (X%getVarIdx(TRIM('vwnd')) .EQ. 0)) THEN
      PRINT *, 'Error use the transFwdNonLinear in UV2VortDivg_t, STOP! '
      STOP
    END IF

    IF (X%getVarIdx(TRIM('wwnd')) .EQ. 0) CALL X%addVar('wwnd')

    ASSOCIATE (uwnd => X%Fields(X%getVarIdx(TRIM('uwnd')))%DATA, &
               vwnd => X%Fields(X%getVarIdx(TRIM('vwnd')))%DATA, &
               wwnd => X%Fields(X%getVarIdx(TRIM('wwnd')))%DATA)

      BLOCK
        wwnd = 0.0D0

        DO i = 2, X%sg%vLevel - 1
          DO j = 1, X%sg%num_icell
            !First term
            DO vStcl = 1, 3
              DO hStcl = 1, 3
                wwnd(i, j, :) = uwnd(this%stcl_uv_1st_v(vStcl, i), this%stcl_u_1st_h(hStcl, j), :) &
                                * this%coef_u_1st(vStcl, hStcl, i, j) + &
                                vwnd(this%stcl_uv_1st_v(vStcl, i), this%stcl_v_1st_h(hStcl, j), :) &
                                * this%coef_v_1st(vStcl, hStcl, i, j) + &
                                wwnd(i, j, :)
              END DO
            END DO

            ! Second term
            DO vStcl = 1, 4
              ! wwnd(i, j, :) = this%coef_u_2nd(vStcl, i, j)
              wwnd(i, j, :) = uwnd(this%stcl_uv_2nd_v(vStcl, i), j, :) * this%coef_u_2nd(vStcl, i, j) + &
                              vwnd(this%stcl_uv_2nd_v(vStcl, i), j, :) * this%coef_v_2nd(vStcl, i, j) + &
                              wwnd(i, j, :)
            END DO
          END DO
        END DO
      END BLOCK

      BLOCK
        DO i = 1, X%sg%num_icell
          DO j = 1, X%sg%tSlots
            BLOCK
              REAL(r_kind) :: WORK(X%sg%vLevel * X%sg%vLevel)
              INTEGER(i_kind) :: LWORK, INFO
              LWORK = X%sg%vLevel * X%sg%vLevel

              CALL DORMQR('L', 'T', X%sg%vLevel, 1, X%sg%vLevel, this%qrMatrix, X%sg%vLevel, &
                          this%tauMatrix, wwnd(:, i, j), X%sg%vLevel, WORK, LWORK, INFO)
              CALL DTRTRS('U', 'N', 'N', X%sg%vLevel, 1, this%qrMatrix, X%sg%vLevel, wwnd(:, i, j), X%sg%vLevel, INFO)
            END BLOCK
          END DO
        END DO
      END BLOCK
    END ASSOCIATE

    CALL X%exHalo()
  END SUBROUTINE

  SUBROUTINE transAdjMultiply(this, dX, X)
    IMPLICIT NONE
    CLASS(UV2VortDivg_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    INTEGER(i_kind)  :: i, j, k, hStcl, vStcl
    TYPE(State_t) :: dXi

    IF ((dX%getVarIdx(TRIM('uwnd')) .EQ. 0) .OR. &
        (dX%getVarIdx(TRIM('wwnd')) .EQ. 0) .OR. &
        (dX%getVarIdx(TRIM('vwnd')) .EQ. 0)) THEN
      PRINT *, 'Error use the transAdjMultiply in UV2VortDivg_t, STOP! '
      STOP
    END IF

    ! Construct an empty state to contain increment
    dXi = dX%zeroCopy()
    dXi%fields(dXi%getVarIdx('wwnd')) = dX%fields(dX%getVarIdx('wwnd'))

    ASSOCIATE (uwnd => dXi%Fields(dXi%getVarIdx(TRIM('uwnd')))%DATA, &
               vwnd => dXi%Fields(dXi%getVarIdx(TRIM('vwnd')))%DATA, &
               wwnd => dXi%Fields(dXi%getVarIdx(TRIM('wwnd')))%DATA)

      BLOCK
        REAL(r_kind) :: WORK(dXi%sg%vLevel * dXi%sg%vLevel)
        INTEGER(i_kind) :: LWORK, INFO
        LWORK = dXi%sg%vLevel * dXi%sg%vLevel

        wwnd(1, :, :) = 0.0D0
        wwnd(dXi%sg%vLevel, :, :) = 0.0D0

        DO i = 1, dXi%sg%num_icell
          DO j = 1, dXi%sg%tSlots
            CALL DTRTRS('U', 'T', 'N', dXi%sg%vLevel, 1, this%qrMatrix, dXi%sg%vLevel, wwnd(:, i, j), dXi%sg%vLevel, INFO)
            CALL DORMQR('L', 'N', dXi%sg%vLevel, 1, dXi%sg%vLevel, this%qrMatrix, dXi%sg%vLevel, &
                        this%tauMatrix, wwnd(:, i, j), dXi%sg%vLevel, WORK, LWORK, INFO)
          END DO
        END DO
      END BLOCK

      BLOCK
        uwnd = 0.0D0
        vwnd = 0.0D0
        DO i = 2, dXi%sg%vLevel - 1
          DO j = 1, dXi%sg%num_icell
            ! First term
            DO vStcl = 1, 3
              DO hStcl = 1, 3
                uwnd(this%stcl_uv_1st_v(vStcl, i), this%stcl_u_1st_h(hStcl, j), :) = &
                  uwnd(this%stcl_uv_1st_v(vStcl, i), this%stcl_u_1st_h(hStcl, j), :) + &
                  this%coef_u_1st(vStcl, hStcl, i, j) * wwnd(i, j, :)

                vwnd(this%stcl_uv_1st_v(vStcl, i), this%stcl_v_1st_h(hStcl, j), :) = &
                  vwnd(this%stcl_uv_1st_v(vStcl, i), this%stcl_v_1st_h(hStcl, j), :) + &
                  this%coef_v_1st(vStcl, hStcl, i, j) * wwnd(i, j, :)

              END DO
            END DO

            ! Second term
            DO vStcl = 1, 4
              uwnd(this%stcl_uv_2nd_v(vStcl, i), j, :) = &
                uwnd(this%stcl_uv_2nd_v(vStcl, i), j, :) + this%coef_u_2nd(vStcl, i, j) * wwnd(i, j, :)
              vwnd(this%stcl_uv_2nd_v(vStcl, i), j, :) = &
                vwnd(this%stcl_uv_2nd_v(vStcl, i), j, :) + this%coef_v_2nd(vStcl, i, j) * wwnd(i, j, :)
            END DO

          END DO
        END DO
      END BLOCK

    END ASSOCIATE

    CALL dXi%rmVar('wwnd')
    CALL dXi%exHaloRevSum()
    CALL dXi%exHalo()

    ! Add the increment with dX
    CALL dx%rmVar('wwnd')
    dX = dX + dXi
  END SUBROUTINE transAdjMultiply

  SUBROUTINE transFwdTanLinear(this, dX, X)
    CLASS(UV2VortDivg_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    CALL this%transFwdNonLinear(dX)
  END SUBROUTINE

  FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
    IMPLICIT NONE
    CLASS(UV2VortDivg_t) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t) :: dX1
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    dX1 = dX
    CALL this%transAdjMultiply(dX1)
  END FUNCTION

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(UV2VortDivg_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1
    INTEGER(i_kind) :: i, j, k

    X1 = X
    CALL this%transFwdNonLinear(X1)
  END FUNCTION

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
    CLASS(UV2VortDivg_t) :: this
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
END MODULE UV2VortDivg_m
