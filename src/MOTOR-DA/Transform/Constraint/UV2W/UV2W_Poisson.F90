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
!! This module contains the UV2W Object which is tranform the variables from control variable space to
!! model space.
MODULE UV2W_Poisson_m
  USE State_m, ONLY: State_t
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t, GridIdx_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE M2MBase_m, ONLY: M2MBase_t
  USE UV2W_m, ONLY: UV2W_t

  TYPE, EXTENDS(UV2W_t) :: UV2W_Poisson_t
    TYPE(State_t), POINTER :: X

    ! First term
    REAL(r_kind), ALLOCATABLE :: coef_p_pXpSigma(:, :, :, :)      ! 3 stcl at verticl / 3 stcl at horizontal / num_vLevel/ num_icell
    REAL(r_kind), ALLOCATABLE :: coef_p_pYpSigma(:, :, :, :)      ! 3 stcl at verticl / 3 stcl at horizontal / num_vLevel/ num_icell
    INTEGER(i_kind), ALLOCATABLE :: stcl_p_pX(:, :)    ! 3 stcl / num_icell
    INTEGER(i_kind), ALLOCATABLE :: stcl_p_pY(:, :)    ! 3 stcl / num_icell
    INTEGER(i_kind), ALLOCATABLE :: stcl_p_pSigma(:, :)   ! 3 stcl / vLevel

    ! Second term
    REAL(r_kind), ALLOCATABLE :: coef_pU_pSigma2(:, :, :)         ! 4 stcl at verticl/ vLevel
    REAL(r_kind), ALLOCATABLE :: coef_pV_pSigma2(:, :, :)         ! 4 stcl at verticl/ vLevel
    INTEGER(i_kind), ALLOCATABLE :: stcl_p_pSigma2(:, :)    ! 4 stcl /vLevel

    ! Third term
    REAL(r_kind), ALLOCATABLE :: coef_p_pX(:, :)      ! 3 stcl at horizontal / num_icell
    REAL(r_kind), ALLOCATABLE :: coef_p_pY(:, :)      ! 3 stcl at horizontal / num_icell
    REAL(r_kind), ALLOCATABLE :: pc1pSigma(:, :)
    REAL(r_kind), ALLOCATABLE :: pc2pSigma(:, :)
    REAL(r_kind), ALLOCATABLE :: coef_p_pSigma(:, :)  ! 3 stcl at verticl/ 3 stcl at horizontal / vLevel / num_icell
    REAL(r_kind), ALLOCATABLE :: c1(:, :), c2(:, :), c3(:)

    ! Possion equation
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
  END TYPE UV2W_Poisson_t

  INTERFACE UV2W_Poisson_t
    PROCEDURE :: constructor
  END INTERFACE UV2W_Poisson_t

CONTAINS

  SUBROUTINE calParameters(this)
    IMPLICIT NONE
    CLASS(UV2W_Poisson_t) :: this
    INTEGER(i_kind) :: i, j, k
    REAL(r_kind) :: t1, t2, t3
    REAL(r_kind), ALLOCATABLE :: psMatrix(:, :)

    ASSOCIATE (sg => this%X%sg)

      ALLOCATE (psMatrix(this%X%sg%vLevel, this%X%sg%vLevel))
      ALLOCATE (this%c1(this%X%sg%vLevel, this%X%sg%num_cell), this%c2(this%X%sg%vLevel, this%X%sg%num_cell), this%c3(this%X%sg%num_cell))
      ALLOCATE (this%pc1pSigma(this%X%sg%vLevel, this%X%sg%num_cell), this%pc2pSigma(this%X%sg%vLevel, this%X%sg%num_cell))

      ! Calculate the C term
      BLOCK
        REAL(r_kind), ALLOCATABLE :: pZspx(:), pZspy(:)

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
          this%c1(i, :) = (sg%ztop - sg%sigma(i)) / (sg%ztop - sg%topo) * pZspx
          this%c2(i, :) = (sg%ztop - sg%sigma(i)) / (sg%ztop - sg%topo) * pZspy
          ! c1(i, :) = sg%zTop*(sg%ztop - sg%zHght(i, :))/((sg%ztop - sg%topo)**2)*pZspx
          ! c2(i, :) = sg%zTop*(sg%ztop - sg%zHght(i, :))/((sg%ztop - sg%topo)**2)*pZspy
          this%c3 = sg%zTop / (sg%zTop - sg%topo)
          this%pc1pSigma(i, :) = -1.0D0 / (sg%ztop - sg%topo) * pZspx / this%c3
          this%pc2pSigma(i, :) = -1.0D0 / (sg%ztop - sg%topo) * pZspy / this%c3
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

        ALLOCATE (this%coef_pU_pSigma2(4, sg%vLevel, sg%num_cell), &         ! 4 stcl at verticl/ vLevel
                  this%coef_pV_pSigma2(4, sg%vLevel, sg%num_cell), &         ! 4 stcl at verticl/ vLevel
                  this%stcl_p_pSigma2(4, sg%vLevel), &    ! 4 stcl /vLevel
                  this%coef_p_pSigma(3, sg%vLevel))

        this%coef_pU_pSigma2 = 0.0D0
        this%coef_pV_pSigma2 = 0.0D0
        this%stcl_p_pSigma2 = 0.0D0

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
          this%stcl_p_pSigma2(:, i) = (/i - 1, i, i + 1, i + 2/)
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
        this%stcl_p_pSigma2(:, i) = (/i - 2, i - 1, i, i + 1/)

        this%qrMatrix = psMatrix

        ! QR Decomposition
        ALLOCATE (this%tauMatrix(sg%vLevel, sg%vLevel))
        CALL DGEQRF(sg%vLevel, sg%vLevel, this%qrMatrix, sg%vLevel, this%tauMatrix, WORK, LWORK, INFO)

        DO i = 1, sg%num_icell
          DO j = 1, 4
            this%coef_pU_pSigma2(j, :, i) = ceof_2nd_v(j, :) * this%c1(:, i) / this%c3(i)
            this%coef_pV_pSigma2(j, :, i) = ceof_2nd_v(j, :) * this%c2(:, i) / this%c3(i)
          END DO
        END DO

        DEALLOCATE (WORK, ceof_2nd_v)
      END BLOCK

      ! Calculate the stencils and coefficients
      BLOCK
        REAL(r_kind) :: h1, h2, h3, a1, a2, a3

        ALLOCATE (this%coef_p_pXpSigma(3, 3, sg%vLevel, sg%num_cell), &      ! 3 stcl at verticl/ 3 stcl at horizontal/ num_icell
                  this%coef_p_pYpSigma(3, 3, sg%vLevel, sg%num_cell), &      ! 3 stcl at verticl/ 3 stcl at horizontal/ num_icell
                  this%stcl_p_pX(3, sg%num_cell), &    ! 3 stcl / num_icell
                  this%stcl_p_pY(3, sg%num_cell), &    ! 3 stcl / num_icell
                  this%stcl_p_pSigma(3, sg%vLevel), &
                  this%coef_p_pX(3, sg%num_cell), &      ! 3 stcl at verticl/ 3 stcl at horizontal/ num_icell
                  this%coef_p_pY(3, sg%num_cell) &      ! 3 stcl at verticl/ 3 stcl at horizontal/ num_icell
                  )   ! 3 stcl / vLevel

        this%coef_p_pXpSigma = 0.0D0
        this%coef_p_pYpSigma = 0.0D0
        this%stcl_p_pX = 0
        this%stcl_p_pY = 0
        this%stcl_p_pSigma = 0

        DO i = 1, sg%num_icell
          IF ((sg%cell_stcl(4, i) .NE. 0) .AND. (sg%cell_stcl(6, i) .NE. 0)) THEN
            t1 = -1.0D0 / (2.0D0 * sg%cell_dist(4, i))
            t2 = 0.0D0
            t3 = 1.0D0 / (2.0D0 * sg%cell_dist(4, i))

            this%stcl_p_pX(1, i) = sg%cell_stcl(4, i)
            this%coef_p_pXpSigma(:, 1, :, i) = t1
            this%stcl_p_pX(2, i) = sg%cell_stcl(5, i)
            this%coef_p_pXpSigma(:, 2, :, i) = t2
            this%stcl_p_pX(3, i) = sg%cell_stcl(6, i)
            this%coef_p_pXpSigma(:, 3, :, i) = t3
            this%coef_p_pX(1, i) = t1
            this%coef_p_pX(2, i) = t2
            this%coef_p_pX(3, i) = t3
          END IF

          IF ((sg%cell_stcl(4, i) .EQ. 0)) THEN
            t1 = -3.0D0 / 2.0D0 / sg%cell_dist(6, i)
            t2 = 2.0D0 / sg%cell_dist(6, i)
            t3 = -1.0D0 / 2.0D0 / sg%cell_dist(6, i)

            this%stcl_p_pX(1, i) = sg%cell_stcl(5, i)
            this%coef_p_pXpSigma(:, 1, :, i) = t1
            this%stcl_p_pX(2, i) = sg%cell_stcl(6, i)
            this%coef_p_pXpSigma(:, 2, :, i) = t2
            this%stcl_p_pX(3, i) = sg%cell_stcl(6, sg%cell_stcl(6, i))
            this%coef_p_pXpSigma(:, 3, :, i) = t3
            this%coef_p_pX(1, i) = t1
            this%coef_p_pX(2, i) = t2
            this%coef_p_pX(3, i) = t3
          END IF

          IF ((sg%cell_stcl(6, i) .EQ. 0)) THEN
            t1 = 1.0D0 / 2.0D0 / sg%cell_dist(4, i)
            t2 = -2.0D0 / sg%cell_dist(4, i)
            t3 = 3.0D0 / 2.0D0 / sg%cell_dist(4, i)

            this%stcl_p_pX(1, i) = sg%cell_stcl(4, sg%cell_stcl(4, i))
            this%coef_p_pXpSigma(:, 1, :, i) = t1
            this%stcl_p_pX(2, i) = sg%cell_stcl(4, i)
            this%coef_p_pXpSigma(:, 2, :, i) = t2
            this%stcl_p_pX(3, i) = sg%cell_stcl(5, i)
            this%coef_p_pXpSigma(:, 3, :, i) = t3
            this%coef_p_pX(1, i) = t1
            this%coef_p_pX(2, i) = t2
            this%coef_p_pX(3, i) = t3
          END IF

          IF ((sg%cell_stcl(2, i) .NE. 0) .AND. (sg%cell_stcl(8, i) .NE. 0)) THEN

            t1 = -1.0D0 / (2.0D0 * sg%cell_dist(2, i))
            t2 = 0.0D0
            t3 = 1.0D0 / (2.0D0 * sg%cell_dist(2, i))

            this%stcl_p_pY(1, i) = sg%cell_stcl(2, i)
            this%coef_p_pYpSigma(:, 1, :, i) = t1
            this%stcl_p_pY(2, i) = sg%cell_stcl(5, i)
            this%coef_p_pYpSigma(:, 2, :, i) = t2
            this%stcl_p_pY(3, i) = sg%cell_stcl(8, i)
            this%coef_p_pYpSigma(:, 3, :, i) = t3
            this%coef_p_pY(1, i) = t1
            this%coef_p_pY(2, i) = t2
            this%coef_p_pY(3, i) = t3
          END IF

          IF ((sg%cell_stcl(2, i) .EQ. 0)) THEN
            t1 = -3.0D0 / 2.0D0 / sg%cell_dist(8, i)
            t2 = 2.0D0 / sg%cell_dist(8, i)
            t3 = -1.0D0 / 2.0D0 / sg%cell_dist(8, i)

            this%stcl_p_pY(1, i) = sg%cell_stcl(5, i)
            this%coef_p_pYpSigma(:, 1, :, i) = t1
            this%stcl_p_pY(2, i) = sg%cell_stcl(8, i)
            this%coef_p_pYpSigma(:, 2, :, i) = t2
            this%stcl_p_pY(3, i) = sg%cell_stcl(8, sg%cell_stcl(8, i))
            this%coef_p_pYpSigma(:, 3, :, i) = t3
            this%coef_p_pY(1, i) = t1
            this%coef_p_pY(2, i) = t2
            this%coef_p_pY(3, i) = t3
          END IF

          IF ((sg%cell_stcl(8, i) .EQ. 0)) THEN
            t1 = 1.0D0 / 2.0D0 / sg%cell_dist(2, i)
            t2 = -2.0D0 / sg%cell_dist(2, i)
            t3 = 3.0D0 / 2.0D0 / sg%cell_dist(2, i)

            this%stcl_p_pY(1, i) = sg%cell_stcl(2, sg%cell_stcl(2, i))
            this%coef_p_pYpSigma(:, 1, :, i) = t1
            this%stcl_p_pY(2, i) = sg%cell_stcl(2, i)
            this%coef_p_pYpSigma(:, 2, :, i) = t2
            this%stcl_p_pY(3, i) = sg%cell_stcl(5, i)
            this%coef_p_pYpSigma(:, 3, :, i) = t3
            this%coef_p_pY(1, i) = t1
            this%coef_p_pY(2, i) = t2
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

          this%coef_p_pXpSigma(1, :, i, :) = this%coef_p_pXpSigma(1, :, i, :) * a1
          this%coef_p_pXpSigma(2, :, i, :) = this%coef_p_pXpSigma(2, :, i, :) * a2
          this%coef_p_pXpSigma(3, :, i, :) = this%coef_p_pXpSigma(3, :, i, :) * a3

          this%coef_p_pYpSigma(1, :, i, :) = this%coef_p_pYpSigma(1, :, i, :) * a1
          this%coef_p_pYpSigma(2, :, i, :) = this%coef_p_pYpSigma(2, :, i, :) * a2
          this%coef_p_pYpSigma(3, :, i, :) = this%coef_p_pYpSigma(3, :, i, :) * a3

          this%coef_p_pSigma(1, i) = a1
          this%coef_p_pSigma(2, i) = a2
          this%coef_p_pSigma(3, i) = a3
        END DO
      END BLOCK

      DO i = 1, 3
        DO j = 1, 3
          DO k = 1, sg%vLevel
            this%coef_p_pXpSigma(i, j, k, :) = -1.0D0 * this%coef_p_pXpSigma(i, j, k, :) / this%c3
            this%coef_p_pYpSigma(i, j, k, :) = -1.0D0 * this%coef_p_pYpSigma(i, j, k, :) / this%c3
          END DO
        END DO
      END DO

      DEALLOCATE (psMatrix)
    END ASSOCIATE

  END SUBROUTINE

  FUNCTION constructor(configFile, X) RESULT(this)
    IMPLICIT NONE
    TYPE(UV2W_Poisson_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X

    this%X => X
    IF (X%sg%vLevel < 3) RETURN
    CALL this%calParameters
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(UV2W_Poisson_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%qrMatrix)) DEALLOCATE (this%qrMatrix)
    IF (ALLOCATED(this%tauMatrix)) DEALLOCATE (this%tauMatrix)

    IF (ALLOCATED(this%coef_p_pXpSigma)) DEALLOCATE (this%coef_p_pXpSigma)
    IF (ALLOCATED(this%coef_p_pYpSigma)) DEALLOCATE (this%coef_p_pYpSigma)
    IF (ALLOCATED(this%stcl_p_pX)) DEALLOCATE (this%stcl_p_pX)
    IF (ALLOCATED(this%stcl_p_pY)) DEALLOCATE (this%stcl_p_pY)
    IF (ALLOCATED(this%stcl_p_pSigma)) DEALLOCATE (this%stcl_p_pSigma)

    IF (ALLOCATED(this%coef_pU_pSigma2)) DEALLOCATE (this%coef_pU_pSigma2)
    IF (ALLOCATED(this%coef_pV_pSigma2)) DEALLOCATE (this%coef_pV_pSigma2)
    IF (ALLOCATED(this%stcl_p_pSigma2)) DEALLOCATE (this%stcl_p_pSigma2)

    IF (ALLOCATED(this%coef_p_pX)) DEALLOCATE (this%coef_p_pX)
    IF (ALLOCATED(this%coef_p_pY)) DEALLOCATE (this%coef_p_pY)

    IF (ALLOCATED(this%pc1pSigma)) DEALLOCATE (this%pc1pSigma)
    IF (ALLOCATED(this%pc2pSigma)) DEALLOCATE (this%pc2pSigma)

    IF (ALLOCATED(this%c1)) DEALLOCATE (this%c1, this%c2, this%c3)

  END SUBROUTINE destructor

  SUBROUTINE transFwdNonLinear(this, X)
    IMPLICIT NONE
    CLASS(UV2W_Poisson_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER :: i, j, k, hStcl, vStcl

    IF (X%sg%vLevel < 3) RETURN

    IF ((X%getVarIdx(TRIM('uwnd')) .EQ. 0) .OR. &
        (X%getVarIdx(TRIM('vwnd')) .EQ. 0)) THEN
      PRINT *, 'Error use the transFwdNonLinear in UV2W_Poisson_t, STOP! '
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
            IF (X%sg%cell_type(j) /= 0) CYCLE

            !First term
            DO vStcl = 1, 3
              DO hStcl = 1, 3
                wwnd(i, j, :) = uwnd(this%stcl_p_pSigma(vStcl, i), this%stcl_p_pX(hStcl, j), :) &
                                * this%coef_p_pXpSigma(vStcl, hStcl, i, j) + &
                                vwnd(this%stcl_p_pSigma(vStcl, i), this%stcl_p_pY(hStcl, j), :) &
                                * this%coef_p_pYpSigma(vStcl, hStcl, i, j) + &
                                wwnd(i, j, :)
              END DO
            END DO

            ! Second term
            DO vStcl = 1, 4
              wwnd(i, j, :) = uwnd(this%stcl_p_pSigma2(vStcl, i), j, :) * this%coef_pU_pSigma2(vStcl, i, j) + &
                              vwnd(this%stcl_p_pSigma2(vStcl, i), j, :) * this%coef_pV_pSigma2(vStcl, i, j) + &
                              wwnd(i, j, :)
            END DO

            ! Third term
            DO vStcl = 1, 3
              wwnd(i, j, :) = uwnd(this%stcl_p_pSigma(vStcl, i), j, :) * this%coef_p_pSigma(vStcl, i) * this%pc1pSigma(i, j) + &
                              vwnd(this%stcl_p_pSigma(vStcl, i), j, :) * this%coef_p_pSigma(vStcl, i) * this%pc2pSigma(i, j) + &
                              wwnd(i, j, :)
            END DO

          END DO
        END DO
      END BLOCK

      BLOCK
        DO i = 1, X%sg%num_icell
          IF (X%sg%cell_type(i) /= 0) CYCLE

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

    ! CALL X%exHalo()
    CALL X%sg%ExchangeMatOnHaloForFieldGrid(X%sg%tSlots, X%sg%vLevel, X%Fields(X%getVarIdx(TRIM('wwnd')))%DATA)

  END SUBROUTINE

  SUBROUTINE transAdjMultiply(this, dX, X)
    IMPLICIT NONE
    CLASS(UV2W_Poisson_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    INTEGER(i_kind)  :: i, j, k, hStcl, vStcl
    TYPE(State_t) :: dXi

    IF ((dX%getVarIdx(TRIM('uwnd')) .EQ. 0) .OR. &
        (dX%getVarIdx(TRIM('wwnd')) .EQ. 0) .OR. &
        (dX%getVarIdx(TRIM('vwnd')) .EQ. 0)) THEN
      PRINT *, 'Error use the transAdjMultiply in UV2W_Poisson_t, STOP! '
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
          IF (dXi%sg%cell_type(i) /= 0) CYCLE
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
            IF (dXi%sg%cell_type(j) /= 0) CYCLE
            ! First term
            DO vStcl = 1, 3
              DO hStcl = 1, 3
                uwnd(this%stcl_p_pSigma(vStcl, i), this%stcl_p_pX(hStcl, j), :) = &
                  uwnd(this%stcl_p_pSigma(vStcl, i), this%stcl_p_pX(hStcl, j), :) + &
                  this%coef_p_pXpSigma(vStcl, hStcl, i, j) * wwnd(i, j, :)

                vwnd(this%stcl_p_pSigma(vStcl, i), this%stcl_p_pY(hStcl, j), :) = &
                  vwnd(this%stcl_p_pSigma(vStcl, i), this%stcl_p_pY(hStcl, j), :) + &
                  this%coef_p_pYpSigma(vStcl, hStcl, i, j) * wwnd(i, j, :)

              END DO
            END DO

            ! Second term
            DO vStcl = 1, 4
              uwnd(this%stcl_p_pSigma2(vStcl, i), j, :) = &
                uwnd(this%stcl_p_pSigma2(vStcl, i), j, :) + this%coef_pU_pSigma2(vStcl, i, j) * wwnd(i, j, :)
              vwnd(this%stcl_p_pSigma2(vStcl, i), j, :) = &
                vwnd(this%stcl_p_pSigma2(vStcl, i), j, :) + this%coef_pV_pSigma2(vStcl, i, j) * wwnd(i, j, :)
            END DO

            ! Third term
            DO vStcl = 1, 3
              uwnd(this%stcl_p_pSigma(vStcl, i), j, :) = &
                uwnd(this%stcl_p_pSigma(vStcl, i), j, :) + &
                this%coef_p_pSigma(vStcl, i) * wwnd(i, j, :) * this%pc1pSigma(i, j)

              vwnd(this%stcl_p_pSigma(vStcl, i), j, :) = &
                vwnd(this%stcl_p_pSigma(vStcl, i), j, :) + &
                this%coef_p_pSigma(vStcl, i) * wwnd(i, j, :) * this%pc2pSigma(i, j)
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
    CLASS(UV2W_Poisson_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    CALL this%transFwdNonLinear(dX)
  END SUBROUTINE

  FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
    IMPLICIT NONE
    CLASS(UV2W_Poisson_t) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t) :: dX1
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    dX1 = dX
    CALL this%transAdjMultiply(dX1)
  END FUNCTION

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(UV2W_Poisson_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1
    INTEGER(i_kind) :: i, j, k

    X1 = X
    CALL this%transFwdNonLinear(X1)
  END FUNCTION

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
    CLASS(UV2W_Poisson_t) :: this
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
END MODULE
