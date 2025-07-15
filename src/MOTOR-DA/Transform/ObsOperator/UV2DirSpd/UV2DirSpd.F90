!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie, Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/08/03, @GBA-MWF, Shenzhen
!     Copying Zilong's OprRadarVel.F90 and modified it for use of wind
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module forwards wind u and v components of a model state
!! toward COS(wdir) and SIN(wdir) and wspd in observation space
MODULE UV2DirSpd_m
  USE State_m, ONLY: State_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE ObsSet_m, ONLY: ObsSet_t
  USE M2OBase_m, ONLY: M2OBase_t
  USE parameters_m, ONLY: degree2radian

  TYPE, EXTENDS(M2OBase_t) :: UV2DirSpd_t
    TYPE(ObsSet_t), POINTER :: Y
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
    FINAL :: destructor
  END TYPE UV2DirSpd_t

  INTERFACE UV2DirSpd_t
    PROCEDURE :: constructor
  END INTERFACE UV2DirSpd_t

CONTAINS

  FUNCTION constructor(configFile, X, Y) RESULT(this)
    IMPLICIT NONE
    TYPE(UV2DirSpd_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X
    TYPE(ObsSet_t), OPTIONAL, TARGET, INTENT(IN) :: Y

    this%X => X
    IF (PRESENT(Y)) THEN
      this%Y => Y
    ELSE
      WRITE (*, 1)
1     FORMAT('UV2DirSpd - transFwdTanLinear_opr requires UV2DirSpd is constructed by option Y but missing! STOP')
      STOP
    END IF
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(UV2DirSpd_t), INTENT(INOUT) :: this
    INTEGER(i_kind) :: i
  END SUBROUTINE destructor

  SUBROUTINE transFwdNonLinear(this, X, Y)
    IMPLICIT NONE
    CLASS(UV2DirSpd_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(ObsSet_t), INTENT(INOUT) :: Y

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, ic, is, iw, iu, iv
    REAL(r_kind) :: wdir, u, v, r

    ic = Y%getObsIdx('cdir')
    is = Y%getObsIdx('sdir')
    iw = Y%getObsIdx('wspd')
    iu = X%getVarIdx('uwnd')
    iv = X%getVarIdx('vwnd')

    ! If needed variables are not present:
    IF ((iu .EQ. 0) .OR. &
        (iv .EQ. 0) .OR. &
        (ic .EQ. 0) .OR. &
        (is .EQ. 0) .OR. &
        (iv .EQ. 0)) RETURN

    ASSOCIATE (uwnd => X%Fields(iu), &
               vwnd => X%Fields(iv))

      ASSOCIATE (cdir => Y%ObsFields(ic), &
                 sdir => Y%ObsFields(is), &
                 wspd => Y%ObsFields(iw))
        DO k = LBOUND(cdir%idx, 1), UBOUND(cdir%idx, 1)
          u = uwnd%Get_Value(cdir%idx(k))
          v = vwnd%Get_Value(cdir%idx(k))
          r = DSQRT(u * u + v * v)
          cdir%values(k) = -v / r
          sdir%values(k) = -u / r
          wspd%values(k) = r * r
        END DO
      END ASSOCIATE
    END ASSOCIATE
  END SUBROUTINE transFwdNonLinear

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(Y)
    IMPLICIT NONE
    CLASS(UV2DirSpd_t) :: this
    TYPE(State_t), INTENT(in) :: X
    TYPE(ObsSet_t) :: Y
    INTEGER(i_kind) :: i, j, k

    Y = this%Y%zeroCopy()
    CALL this%transFwdNonLinear(X, Y)
  END FUNCTION

  SUBROUTINE transAdjMultiply(this, D, dX, X)
    !> @brief
    !! This multiplies D by the adjoint operator to dx:
    !! Assume the Jacobian is L, dX = L^T D:
    !! Note: the wind direction is defined as pointing to the north
    !! where the wind originated.
    IMPLICIT NONE
    CLASS(UV2DirSpd_t) :: this
    TYPE(ObsSet_t), INTENT(IN) :: D
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, iu, iv, ic, is, iw
    REAL(r_kind) :: u, v, r, dp, dq, dw

    ! This is a nonlinear operator and the current state is required:
    IF (.NOT. PRESENT(X)) THEN
      WRITE (*, 2)
2     FORMAT('UV2DirSpd - transAdjMultiply requires a state but it is missing, check and rerun!')
      STOP
    END IF

    iu = dX%getVarIdx(TRIM('uwnd'))
    iv = dX%getVarIdx(TRIM('vwnd'))
    ic = D%getObsIdx(TRIM('cdir'))
    is = D%getObsIdx(TRIM('sdir'))
    iw = D%getObsIdx(TRIM('wspd'))
    IF (iu .EQ. 0 .OR. iv .EQ. 0 .OR. &
        ic .EQ. 0 .OR. is .EQ. 0 .OR. iw .EQ. 0) THEN
      WRITE (*, 1)
1     FORMAT('Error use the UV2DirSpd_t in transAdjMultiply', /, &
             ' check your yaml for uwnd and vwnd for model state and cdir,sdir and wspd!')
      STOP
    END IF

    ASSOCIATE (uwnd => X%Fields(iu), &
               vwnd => X%Fields(iv), &
               duwnd => dX%Fields(iu), &
               dvwnd => dX%Fields(iv))

      duwnd%DATA = 0.0D0
      dvwnd%DATA = 0.0D0

      ASSOCIATE (cdir => D%ObsFields(ic), &
                 sdir => D%ObsFields(is), &
                 wspd => D%ObsFields(iw))
        DO k = LBOUND(cdir%idx, 1), UBOUND(cdir%idx, 1)
          u = uwnd%Get_Value(cdir%idx(k))
          v = vwnd%Get_Value(cdir%idx(k))
          dp = cdir%values(k)
          dq = sdir%values(k)
          dw = wspd%values(k)
          r = DSQRT(u * u + v * v)**3
          CALL duwnd%Add_Value(cdir%idx(k), (u * v * dp - v * v * dq) / r + 2.0D0 * u * dw)
          CALL dvwnd%Add_Value(cdir%idx(k), (-u * u * dp + u * v * dq) / r + 2.0D0 * v * dw)
        END DO
      END ASSOCIATE

    END ASSOCIATE

    CALL dX%exHalo()
  END SUBROUTINE transAdjMultiply

  FUNCTION transAdjMultiply_opr(this, D, X) RESULT(dX)
    IMPLICIT NONE
    CLASS(UV2DirSpd_t) :: this
    TYPE(ObsSet_t), INTENT(IN) :: D
    TYPE(State_t) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    dX = X%zeroCopy()
    CALL this%transAdjMultiply(D, dX, X)
  END FUNCTION

  ! UV2DirSpd is a nonlinear operator:
  SUBROUTINE transFwdTanLinear(this, dX, dY, X)
    IMPLICIT NONE
    CLASS(UV2DirSpd_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(ObsSet_t), INTENT(INOUT) :: dY
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    ! Local variables:
    INTEGER(i_kind) :: idu, idv, iu, iv, ic, is, iw, k
    REAL(r_kind) :: u, v, r, du, dv, dp, dq, dw

    ! For a nonlinear operator, X is required:
    IF (.NOT. PRESENT(X)) THEN
      WRITE (*, 2)
2     FORMAT('UV2DirSpd - transFwdNonLinear requires a state but it is missing, check and rerun!')
      STOP
    END IF

    idu = dX%getVarIdx(TRIM('uwnd'))
    idv = dX%getVarIdx(TRIM('vwnd'))
    iu = X%getVarIdx(TRIM('uwnd'))
    iv = X%getVarIdx(TRIM('vwnd'))
    ic = dY%getObsIdx(TRIM('cdir'))
    is = dY%getObsIdx(TRIM('sdir'))
    iw = dY%getObsIdx(TRIM('wspd'))
    IF (iu .EQ. 0 .OR. iv .EQ. 0 .OR. idu .EQ. 0 .OR. idv .EQ. 0 .OR. &
        ic .EQ. 0 .OR. is .EQ. 0 .OR. iw .EQ. 0) THEN
      WRITE (*, 1)
1     FORMAT('Error in transFwdNonLinear of UV2DirSpd_t ', /, &
             ' check your yaml for uwnd and vwnd for model state and cdir,sdir and wspd!')
      STOP
    END IF

    ! F(X):
    CALL this%transFwdNonLinear(X, this%Y)

    ! H*dX:
    ASSOCIATE (uwnd => X%Fields(iu), &
               vwnd => X%Fields(iv), &
               duwnd => dX%Fields(iu), &
               dvwnd => dX%Fields(iv), &
               dcdir => dY%ObsFields(ic), &
               dsdir => dY%ObsFields(is), &
               dwspd => dY%ObsFields(iw))

      dcdir%values = 0.0D0
      dsdir%values = 0.0D0
      dwspd%values = 0.0D0

      DO k = LBOUND(dcdir%idx, 1), UBOUND(dcdir%idx, 1)
        u = uwnd%Get_Value(dcdir%idx(k))
        v = vwnd%Get_Value(dcdir%idx(k))
        du = duwnd%Get_Value(dcdir%idx(k))
        dv = dvwnd%Get_Value(dcdir%idx(k))
        dp = dcdir%values(k)
        dq = dsdir%values(k)
        dw = dwspd%values(k)
        r = DSQRT(u * u + v * v)**3

        dcdir%values(k) = (u * v * du - u * u * dv) / r
        dsdir%values(k) = (u * v * dv - v * v * du) / r
        dwspd%values(k) = 2.0D0 * (u * du + v * dv)
      END DO
    END ASSOCIATE

    ! Tangent linear: F(X)+H*dX:
    dY = dY + this%Y  ! Note for ObsSet_t addition, it has to use dY + this%Y; otherwise it ends 2*this%Y

  END SUBROUTINE transFwdTanLinear

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dY)
    CLASS(UV2DirSpd_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(ObsSet_t) :: dY

    CALL this%transFwdTanLinear(dX, dY, X)
  END FUNCTION

END MODULE UV2DirSpd_m
