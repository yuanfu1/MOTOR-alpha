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
!!
!! method.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention

MODULE Ctl2State_m
  USE M2MBase_m, ONLY: M2MBase_t
  USE State_m, ONLY: State_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE RhoVCtl2RhoV_Power_m, ONLY: RhoVCtl2RhoV_Power_t
  USE TempCtl2Temp_Power_m, ONLY: TempCtl2Temp_Power_t
  USE QVCtl2QV_sv_m, ONLY: QVCtl2QV_sv_t

  USE YAMLRead_m
  IMPLICIT NONE

! #define TRACK_CTL

  TYPE, EXTENDS(M2MBase_t) :: Ctl2State_t
    CHARACTER(LEN=20), ALLOCATABLE :: ctlVarNames(:),ctlVarHead(:)
    CHARACTER(LEN=20), ALLOCATABLE :: stateVarNames(:), CtlStateVarNames(:)
    INTEGER(i_kind) :: numCtlVars, numStateVars
    INTEGER(i_kind), ALLOCATABLE :: nameHeadLen(:)
  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: initialize
    PROCEDURE, PUBLIC, PASS(this) :: transBackward

    PROCEDURE, PUBLIC, PASS(this) :: transFwdNonLinear_opr
    PROCEDURE, PUBLIC, PASS(this) :: transFwdTanLinear_opr
    PROCEDURE, PUBLIC, PASS(this) :: transAdjMultiply_opr

    PROCEDURE, PUBLIC, PASS(this) :: transFwdNonLinear
    PROCEDURE, PUBLIC, PASS(this) :: transFwdTanLinear
    PROCEDURE, PUBLIC, PASS(this) :: transAdjMultiply

    PROCEDURE, PUBLIC, PASS(this) :: SumState
    PROCEDURE, PUBLIC, PASS(this) :: SumBothCtlState
    PROCEDURE, PUBLIC, PASS(this) :: SubState

    PROCEDURE :: fwdNL_opr => transFwdNonLinear_opr
    PROCEDURE :: fwdTL_opr => transFwdTanLinear_opr
    PROCEDURE :: adjMul_opr => transAdjMultiply_opr

    PROCEDURE :: fwdNL => transFwdNonLinear
    PROCEDURE :: fwdTL => transFwdTanLinear
    PROCEDURE :: adjMul => transAdjMultiply

    FINAL :: destructor
  END TYPE Ctl2State_t

CONTAINS

  SUBROUTINE initialize(this, configFile)
    IMPLICIT NONE
    CLASS(Ctl2State_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    ! Local variables:
    CHARACTER(LEN=20) name
    INTEGER(i_kind) :: i,istatus,j

    istatus = yaml_get_var(TRIM(configFile), 'modelState', 'varList', this%stateVarNames)
    this%numStateVars = UBOUND(this%stateVarNames, 1)

    istatus = yaml_get_var(TRIM(configFile), 'modelState', 'ctlVarNames', this%ctlVarNames)
    PRINT *, 'ctlVarName read status: ', istatus
    IF (istatus .NE. 0) THEN
      this%numCtlVars = 0
    ELSE
      this%numCtlVars = UBOUND(this%ctlVarNames, 1)
      WRITE (*, 11) this%numCtlVars
11    FORMAT("Ctl2State: Number of control vars: ", I2)
      WRITE (*, 12) (this%ctlVarNames(istatus), istatus=1, this%numCtlVars)
12    FORMAT('Ctl2State: control names: ', 50(1X, A8))

      ! Set the nameHeadLen:
      ALLOCATE (this%ctlVarHead(this%numCtlVars), this%nameHeadLen(this%numCtlVars))
      DO istatus = 1, this%numCtlVars
        this%nameHeadLen(istatus) = 0
        name = this%ctlVarNames(istatus)
        DO i = 1, 8
          IF (name(i:i) .EQ. '_') THEN
            this%nameHeadLen(istatus) = i - 1
            IF (this%nameHeadLen(istatus) .LT. 1) THEN
              WRITE (*, 14) istatus, name
14            FORMAT('Ctl2State - invalid control var name ', I2, A)
              STOP
            END IF
            this%ctlVarHead(istatus) = name(1:this%nameHeadLen(istatus))

            EXIT
          END IF
        END DO
      END DO
      WRITE (*, 13) (this%nameHeadLen(istatus), istatus=1, this%numCtlVars)
13    FORMAT('Ctl2State: control name head len: ', 50I2)
    END IF
    WRITE (*, 15) (this%ctlVarNames(istatus) (1:this%nameHeadLen(istatus)), istatus=1, this%numCtlVars)
15  FORMAT('Ctl2State: control origs: ', 50(1X, A8))

    ALLOCATE(this%CtlStateVarNames(this%numStateVars))
    DO i = 1, this%numStateVars
      this%CtlStateVarNames(i) = TRIM(this%stateVarNames(i))
      DO j = 1, this%numCtlVars
        IF (TRIM(this%stateVarNames(i)) .EQ. TRIM(this%ctlVarHead(j))) THEN
          this%CtlStateVarNames(i) = ''
          this%CtlStateVarNames(i) = TRIM(this%ctlVarNames(j))
        END IF
      END DO
    END DO
    
    PRINT *, 'this%ctlVarNames: ', this%ctlVarNames
    PRINT *, 'this%ctlVarHead: ', this%ctlVarHead
    PRINT *, 'this%stateVarNames: ', this%stateVarNames
    PRINT *, 'this%CtlStateVarNames: ', this%CtlStateVarNames

  END SUBROUTINE initialize

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(Ctl2State_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%ctlVarHead)) DEALLOCATE(this%ctlVarHead)
    IF (ALLOCATED(this%nameHeadLen)) DEALLOCATE(this%nameHeadLen)
    IF (ALLOCATED(this%stateVarNames)) DEALLOCATE(this%stateVarNames)
    IF (ALLOCATED(this%CtlStateVarNames)) DEALLOCATE(this%CtlStateVarNames)

  END SUBROUTINE destructor

  SUBROUTINE SumState(this, X_Inc2FS, XbRef)
    ! X_Inc2FS: State of Increment to FullState
    ! XbRef: Background State
    IMPLICIT NONE
    CLASS(Ctl2State_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X_Inc2FS
    TYPE(State_t), INTENT(IN) :: XbRef

    ! Local variables:
    CHARACTER(LEN=20) :: name
    INTEGER(i_kind) :: i,j,k

    DO i=1,this%numStateVars
      name = this%stateVarNames(i)

      DO j = LBOUND(X_Inc2FS%fields, 1), UBOUND(X_Inc2FS%fields, 1)
        IF (TRIM(name) .EQ. TRIM(X_Inc2FS%fields(j)%Get_Name())) THEN
          ! PRINT *, 'Check SumState: ', TRIM(X_Inc2FS%fields(j)%Get_Name())
          X_Inc2FS%fields(j)%data = X_Inc2FS%fields(j)%data + XbRef%fields(j)%data
        END IF
      END DO
    END DO

  END SUBROUTINE SumState

  SUBROUTINE SumBothCtlState(this, X_Inc2FS, XbRef)
    ! X_Inc2FS: State of Increment to FullState
    ! XbRef: Background State
    IMPLICIT NONE
    CLASS(Ctl2State_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X_Inc2FS
    TYPE(State_t), INTENT(IN) :: XbRef

    ! Local variables:
    CHARACTER(LEN=20) :: name
    INTEGER(i_kind) :: i,j,k

    DO i=1,this%numStateVars
      name = this%CtlStateVarNames(i)

      DO j = LBOUND(X_Inc2FS%fields, 1), UBOUND(X_Inc2FS%fields, 1)
        IF (TRIM(name) .EQ. TRIM(X_Inc2FS%fields(j)%Get_Name())) THEN
          ! PRINT *, 'Check SumBothCtlState: ', TRIM(name), ' is implemented'
          X_Inc2FS%fields(j)%data = X_Inc2FS%fields(j)%data + XbRef%fields(j)%data
        END IF
      END DO
    END DO

  END SUBROUTINE SumBothCtlState

  SUBROUTINE SubState(this, X_FS2Inc, XbRef)
    ! X_FS2Inc: State of FullState to Increment
    ! XbRef: Background State
    IMPLICIT NONE
    CLASS(Ctl2State_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X_FS2Inc
    TYPE(State_t), INTENT(IN) :: XbRef

    ! Local variables:
    CHARACTER(LEN=20) :: name
    INTEGER(i_kind) :: i,j,k

    DO i=1,this%numStateVars

      name = this%stateVarNames(i)

      DO j = LBOUND(X_FS2Inc%fields, 1), UBOUND(X_FS2Inc%fields, 1)
        IF (TRIM(name) .EQ. TRIM(X_FS2Inc%fields(j)%Get_Name())) THEN
          ! PRINT *, 'Check SubState: ', TRIM(name), ' is implemented'
          X_FS2Inc%fields(j)%data = X_FS2Inc%fields(j)%data - XbRef%fields(j)%data
        END IF
      END DO
    END DO

  END SUBROUTINE SubState
  
  !> @brief
  !! Transfer the model state back to controls:
  SUBROUTINE transBackward(this, X)
    IMPLICIT NONE
    CLASS(Ctl2State_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X

    ! Local variables:
    CHARACTER(LEN=20) :: name
    INTEGER(i_kind) :: i, j, k

    ! For all controls, converting:
    DO i = 1, this%numCtlVars
      name = this%ctlVarNames(i)
#ifdef TRACK_CTL
      WRITE (*, 1) i, TRIM(name), X%sg%mpddInfo_sg%myrank
1     FORMAT('transBackward: ', I2, 1X, A, ' pc', I2)
#endif

      DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
        IF (name(1:this%nameHeadLen(i)) .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN

          ! Possible controls:
          SELECT CASE (TRIM(name))
          CASE ('pres_ctl')
            CALL X%fields(j)%Set_Name('pres_ctl')
            X%fields(j)%DATA = X%fields(j)%DATA / 100.0D0
#ifdef TRACK_CTL
            PRINT *, 'transBackward pres is converted to pres_ctl', X%sg%mpddInfo_sg%myrank
#endif

          CASE ('psl_ctl')
            CALL X%fields(j)%Set_Name('psl_ctl')
            X%fields(j)%DATA = X%fields(j)%DATA / 100.0D0
#ifdef TRACK_CTL
            PRINT *, 'transBackward psl is converted to psl_ctl', X%sg%mpddInfo_sg%myrank
#endif

          CASE ('qvapor_ctl')
            CALL X%fields(j)%Set_Name('qvapor_ctl')
            ! X%fields(j)%data = X%fields(j)%data/0.001D0
            DO k = 1, X%sg%tSlots
              X%fields(j)%DATA(:, :, k) = X%fields(j)%DATA(:, :, k) / X%sg%SVapor
            END DO
#ifdef TRACK_CTL
            PRINT *, 'transBackward qvapor is converted to qvapor_ctl with min SVapor', MINVAL(X%sg%SVapor), X%sg%mpddInfo_sg%myrank
#endif
          CASE ('qcloud_ctl')
            CALL X%fields(j)%Set_Name('qcloud_ctl')
            ! X%fields(j)%data = X%fields(j)%data/0.001D0
            DO k=1, X%sg%tSlots
              X%fields(j)%data(:, :, k) = X%fields(j)%data(:, :, k)/X%sg%SCloud
            ENDDO
#ifdef TRACK_CTL
  print*,'transBackward qcloud is converted to qcloud_ctl',X%sg%mpddInfo_sg%myrank
#endif
          CASE ('qice_ctl')
            CALL X%fields(j)%Set_Name('qice_ctl')
            ! X%fields(j)%data = X%fields(j)%data/0.001D0
            DO k=1, X%sg%tSlots
              X%fields(j)%data(:, :, k) = X%fields(j)%data(:, :, k)/X%sg%SIce
            ENDDO
#ifdef TRACK_CTL
print*,'transBackward qice is converted to qice_ctl',X%sg%mpddInfo_sg%myrank
#endif
          CASE ('qrain_ctl')
            CALL X%fields(j)%Set_Name('qrain_ctl')
            ! X%fields(j)%data = X%fields(j)%data/0.001D0
            DO k=1, X%sg%tSlots
              X%fields(j)%data(:, :, k) = X%fields(j)%data(:, :, k)/X%sg%SRain
            ENDDO
#ifdef TRACK_CTL
print*,'transBackward qrain is converted to qrain_ctl',X%sg%mpddInfo_sg%myrank
#endif
          CASE ('qsnow_ctl')
            CALL X%fields(j)%Set_Name('qsnow_ctl')
            ! X%fields(j)%data = X%fields(j)%data/0.001D0
            DO k=1, X%sg%tSlots
              X%fields(j)%data(:, :, k) = X%fields(j)%data(:, :, k)/X%sg%SSnow
            ENDDO
#ifdef TRACK_CTL
print*,'transBackward qsnow is converted to qsnow_ctl',X%sg%mpddInfo_sg%myrank
#endif
          CASE ('qgraupel_ctl')
            CALL X%fields(j)%Set_Name('qgraupel_ctl')
            ! X%fields(j)%data = X%fields(j)%data/0.001D0
            DO k=1, X%sg%tSlots
              X%fields(j)%data(:, :, k) = X%fields(j)%data(:, :, k)/X%sg%SGraupel
            ENDDO
#ifdef TRACK_CTL
print*,'transBackward qraupel is converted to qgraupel_ctl',X%sg%mpddInfo_sg%myrank
#endif
          END SELECT
        END IF
      END DO
    END DO
  END SUBROUTINE transBackward

  !> @brief
  !! Transfer the controls to model states:
  SUBROUTINE transFwdNonLinear(this, X)
    IMPLICIT NONE
    CLASS(Ctl2State_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X

    ! Local variables:
    CHARACTER(LEN=20) :: name
    INTEGER(i_kind) :: i, j, k

    ! For all controls:
    DO i = 1, this%numCtlVars
      name = this%ctlVarNames(i)
#ifdef TRACK_CTL
      WRITE (*, 1) i, TRIM(name), X%sg%mpddInfo_sg%myrank
1     FORMAT('transFwdNonLinear: ', I2, 1X, A, ' pc', I2)
#endif

      DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
        IF (TRIM(name) .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN

          ! Possible controls:
          SELECT CASE (TRIM(name))
          CASE ('pres_ctl')
            CALL X%fields(j)%Set_Name('pres')
            X%fields(j)%DATA = X%fields(j)%DATA * 100.0D0
#ifdef TRACK_CTL
            PRINT *, 'transFwdNonLinear pres_ctl is converted to pres', &
              X%sg%mpddInfo_sg%myrank, X%sg%gLevel
#endif

          CASE ('psl_ctl')
            CALL X%fields(j)%Set_Name('psl')
            X%fields(j)%DATA = X%fields(j)%DATA * 100.0D0
#ifdef TRACK_CTL
            PRINT *, 'transFwdNonLinear psl_ctl is converted to psl', &
              X%sg%mpddInfo_sg%myrank, X%sg%gLevel
#endif

          CASE ('qvapor_ctl')
            CALL X%fields(j)%Set_Name('qvapor')
            ! X%fields(j)%data = X%fields(j)%data*0.001D0
            DO k = 1, X%sg%tSlots
              X%fields(j)%DATA(:, :, k) = X%fields(j)%DATA(:, :, k) * X%sg%SVapor
            END DO
#ifdef TRACK_CTL
            PRINT *, 'transFwdNonLinear qvapor_ctl is converted to qvapor', &
              X%sg%mpddInfo_sg%myrank, X%sg%gLevel
#endif
          CASE ('qcloud_ctl')
            CALL X%fields(j)%Set_Name('qcloud')
            ! X%fields(j)%data = X%fields(j)%data*0.001D0
            DO k=1, X%sg%tSlots
              X%fields(j)%data(:, :, k) = X%fields(j)%data(:, :, k)*X%sg%SCloud
            ENDDO
#ifdef TRACK_CTL
  print*,'transFwdNonLinear qcloud_ctl is converted to qcloud',X%sg%mpddInfo_sg%myrank
#endif
          CASE ('qice_ctl')
            CALL X%fields(j)%Set_Name('qice')
            ! X%fields(j)%data = X%fields(j)%data*0.001D0
            DO k=1, X%sg%tSlots
              X%fields(j)%data(:, :, k) = X%fields(j)%data(:, :, k)*X%sg%SIce
            ENDDO
#ifdef TRACK_CTL
  print*,'transFwdNonLinear qice_ctl is converted to qice',X%sg%mpddInfo_sg%myrank
#endif
          CASE ('qrain_ctl')
            CALL X%fields(j)%Set_Name('qrain')
            ! X%fields(j)%data = X%fields(j)%data*0.001D0
            DO k=1, X%sg%tSlots
              X%fields(j)%data(:, :, k) = X%fields(j)%data(:, :, k)*X%sg%SRain
            ENDDO
#ifdef TRACK_CTL
  print*,'transFwdNonLinear qrain_ctl is converted to qrain',X%sg%mpddInfo_sg%myrank
#endif
          CASE ('qsnow_ctl')
            CALL X%fields(j)%Set_Name('qsnow')
            ! X%fields(j)%data = X%fields(j)%data*0.001D0
            DO k=1, X%sg%tSlots
              X%fields(j)%data(:, :, k) = X%fields(j)%data(:, :, k)*X%sg%SSnow
            ENDDO
#ifdef TRACK_CTL
  print*,'transFwdNonLinear qsnow_ctl is converted to qsnow',X%sg%mpddInfo_sg%myrank
#endif
          CASE ('qgraupel_ctl')
            CALL X%fields(j)%Set_Name('qgraupel')
            ! X%fields(j)%data = X%fields(j)%data*0.001D0
            DO k=1, X%sg%tSlots
              X%fields(j)%data(:, :, k) = X%fields(j)%data(:, :, k)*X%sg%SGraupel
            ENDDO
#ifdef TRACK_CTL
  print*,'transFwdNonLinear qgraupel_ctl is converted to qgraupel',X%sg%mpddInfo_sg%myrank
#endif
          END SELECT

        END IF
      END DO
    END DO
  END SUBROUTINE transFwdNonLinear

  SUBROUTINE transFwdTanLinear(this, dX, X)
    CLASS(Ctl2State_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    ! CALL this%transFwdNonLinear(dX)

    ! Local variables:
    CHARACTER(LEN=20) :: name
    INTEGER :: i, j, k

    ! As there is nonlinear operators involved, X must be present:
    IF (.NOT. PRESENT(X)) THEN
      WRITE (*, 1) dX%sg%mpddInfo_sg%myrank
1     FORMAT('Ctl2State transFwdTanLinear: X must be present! pc', I2)
      CALL dX%sg%mpddInfo_sg%barrier()
      STOP
    END IF

    ! For all controls:
    DO i = 1, this%numCtlVars
      name = this%ctlVarNames(i)
#ifdef TRACK_CTL
      WRITE (*, 2) i, TRIM(name), X%sg%mpddInfo_sg%myrank, &
        (TRIM(X%fields(j)%Get_Name()), j=LBOUND(X%fields, 1), UBOUND(X%fields, 1))

2     FORMAT('transFwdTanLinear: ', I2, 1X, A, ' pc', I2, ' ', 100A)
#endif

      DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
        IF (TRIM(name) .EQ. (TRIM(X%fields(j)%Get_Name()))) THEN

          ! Possible controls:
          SELECT CASE (TRIM(name))
          CASE ('pres_ctl')
            CALL dX%fields(j)%Set_Name('pres')
            dX%fields(j)%DATA = dX%fields(j)%DATA * 100.0D0
#ifdef TRACK_CTL
            PRINT *, 'transFwdTanLinear pres_ctl is converted to pres', X%sg%mpddInfo_sg%myrank
#endif

          CASE ('psl_ctl')
            CALL dX%fields(j)%Set_Name('psl')
            dX%fields(j)%DATA = dX%fields(j)%DATA * 100.0D0
#ifdef TRACK_CTL
            PRINT *, 'transFwdTanLinear psl_ctl is converted to psl', X%sg%mpddInfo_sg%myrank
#endif

          CASE ('qvapor_ctl')
            CALL dX%fields(j)%Set_Name('qvapor')
            DO k = 1, X%sg%tSlots
              dX%fields(j)%DATA(:, :, k) = dX%fields(j)%DATA(:, :, k) * X%sg%SVapor
            END DO
#ifdef TRACK_CTL
            PRINT *, 'transFwdTanLinear qvapor_ctl is converted to qvapor', X%sg%mpddInfo_sg%myrank
#endif
          CASE ('qcloud_ctl')
            CALL dX%fields(j)%Set_Name('qcloud')
            DO k=1, X%sg%tSlots
              dX%fields(j)%data(:, :, k) = dX%fields(j)%data(:, :, k) * X%sg%SCloud
            ENDDO
#ifdef TRACK_CTL
  print*,'transFwdTanLinear qcloud_ctl is converted to qcloud',X%sg%mpddInfo_sg%myrank
#endif
          CASE ('qice_ctl')
            CALL dX%fields(j)%Set_Name('qice')
            DO k=1, X%sg%tSlots
              dX%fields(j)%data(:, :, k) = dX%fields(j)%data(:, :, k) * X%sg%SIce
            ENDDO
#ifdef TRACK_CTL
  print*,'transFwdTanLinear qice_ctl is converted to qice',X%sg%mpddInfo_sg%myrank
#endif
          CASE ('qrain_ctl')
            CALL dX%fields(j)%Set_Name('qrain')
            DO k=1, X%sg%tSlots
              dX%fields(j)%data(:, :, k) = dX%fields(j)%data(:, :, k) * X%sg%SRain
            ENDDO
#ifdef TRACK_CTL
  print*,'transFwdTanLinear qrain_ctl is converted to qrain',X%sg%mpddInfo_sg%myrank
#endif
          CASE ('qsnow_ctl')
            CALL dX%fields(j)%Set_Name('qsnow')
            DO k=1, X%sg%tSlots
              dX%fields(j)%data(:, :, k) = dX%fields(j)%data(:, :, k) * X%sg%SSnow
            ENDDO
#ifdef TRACK_CTL
  print*,'transFwdTanLinear qsnow_ctl is converted to qsnow',X%sg%mpddInfo_sg%myrank
#endif
          CASE ('qgraupel_ctl')
            CALL dX%fields(j)%Set_Name('qgraupel')
            DO k=1, X%sg%tSlots
              dX%fields(j)%data(:, :, k) = dX%fields(j)%data(:, :, k) * X%sg%SGraupel
            ENDDO
#ifdef TRACK_CTL
  print*,'transFwdTanLinear qgraupel_ctl is converted to qgraupel',X%sg%mpddInfo_sg%myrank
#endif
          END SELECT

        END IF
      END DO
    END DO
  END SUBROUTINE transFwdTanLinear

  SUBROUTINE transAdjMultiply(this, dX, X)
    IMPLICIT NONE
    CLASS(Ctl2State_t) :: this
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    ! Local variables:
    CHARACTER(LEN=20) :: name
    INTEGER :: i, j, k

#ifdef TRACK_CTL
    PRINT *, 'Entering transAdjMultiply....', this%numCtlVars, PRESENT(X)
#endif

    ! Check current state is in as required by a nonlinear operator:
    IF (.NOT. PRESENT(X)) THEN
      WRITE (*, 1) dX%sg%mpddInfo_sg%myrank
1     FORMAT('Ctl2State transAdjMultiply: X must be present for a nonlinear operator ! pc', I2)
      CALL dX%sg%mpddInfo_sg%barrier()
      STOP
    END IF

    ! For all controls:
    DO i = 1, this%numCtlVars
      name = this%ctlVarNames(i)
#ifdef TRACK_CTL
      WRITE (*, 2) i, TRIM(name), dX%sg%mpddInfo_sg%myrank
2     FORMAT('transAdjMultiply: ', I2, 1X, A, ' pc', I2)
#endif

      DO j = LBOUND(dX%fields, 1), UBOUND(dX%fields, 1)
        IF (name(1:this%nameHeadLen(i)) .EQ. (TRIM(dX%fields(j)%Get_Name()))) THEN

          ! Possible controls:
          SELECT CASE (TRIM(name))
          CASE ('pres_ctl')
            CALL dX%fields(j)%Set_Name('pres_ctl')
            dX%fields(j)%DATA = dX%fields(j)%DATA * 100.0D0
#ifdef TRACK_CTL
            PRINT *, 'transAdjMultiply pres is converted to pres_ctl', dX%sg%mpddInfo_sg%myrank
#endif

          CASE ('psl_ctl')
            CALL dX%fields(j)%Set_Name('psl_ctl')
            dX%fields(j)%DATA = dX%fields(j)%DATA * 100.0D0
#ifdef TRACK_CTL
            PRINT *, 'transAdjMultiply psl is converted to psl_ctl', dX%sg%mpddInfo_sg%myrank
#endif

          CASE ('qvapor_ctl')
            CALL dX%fields(j)%Set_Name('qvapor_ctl')
            ! dX%fields(j)%data = dX%fields(j)%data*0.001D0
            DO k = 1, dX%sg%tSlots
              dX%fields(j)%DATA(:, :, k) = dX%fields(j)%DATA(:, :, k) * dX%sg%SVapor
            END DO
#ifdef TRACK_CTL
            PRINT *, 'transAdjMultiply qvapor is converted to qvapor_ctl', dX%sg%mpddInfo_sg%myrank
#endif
          CASE ('qcloud_ctl')
            CALL dX%fields(j)%Set_Name('qcloud_ctl')
            ! dX%fields(j)%data = dX%fields(j)%data*0.001D0
            DO k=1, dX%sg%tSlots
              dX%fields(j)%data(:, :, k) = dX%fields(j)%data(:, :, k)*dX%sg%SCloud
            ENDDO
#ifdef TRACK_CTL
  print*,'transAdjMultiply qcloud is converted to qcloud_ctl',dX%sg%mpddInfo_sg%myrank
#endif
          CASE ('qice_ctl')
            CALL dX%fields(j)%Set_Name('qice_ctl')
            ! dX%fields(j)%data = dX%fields(j)%data*0.001D0
            DO k=1, dX%sg%tSlots
              dX%fields(j)%data(:, :, k) = dX%fields(j)%data(:, :, k)*dX%sg%SIce
            ENDDO
#ifdef TRACK_CTL
  print*,'transAdjMultiply qice is converted to qice_ctl',dX%sg%mpddInfo_sg%myrank
#endif
          CASE ('qrain_ctl')
            CALL dX%fields(j)%Set_Name('qrain_ctl')
            ! dX%fields(j)%data = dX%fields(j)%data*0.001D0
            DO k=1, dX%sg%tSlots
              dX%fields(j)%data(:, :, k) = dX%fields(j)%data(:, :, k)*dX%sg%SRain
            ENDDO
#ifdef TRACK_CTL
  print*,'transAdjMultiply qrain is converted to qrain_ctl',dX%sg%mpddInfo_sg%myrank
#endif
          CASE ('qsnow_ctl')
            CALL dX%fields(j)%Set_Name('qsnow_ctl')
            ! dX%fields(j)%data = dX%fields(j)%data*0.001D0
            DO k=1, dX%sg%tSlots
              dX%fields(j)%data(:, :, k) = dX%fields(j)%data(:, :, k)*dX%sg%SSnow
            ENDDO
#ifdef TRACK_CTL
  print*,'transAdjMultiply qsnow is converted to qsnow_ctl',dX%sg%mpddInfo_sg%myrank
#endif
          CASE ('qgraupel_ctl')
            CALL dX%fields(j)%Set_Name('qgraupel_ctl')
            ! dX%fields(j)%data = dX%fields(j)%data*0.001D0
            DO k=1, dX%sg%tSlots
              dX%fields(j)%data(:, :, k) = dX%fields(j)%data(:, :, k)*dX%sg%SGraupel
            ENDDO
#ifdef TRACK_CTL
  print*,'transAdjMultiply qgraupel is converted to qgraupel_ctl',dX%sg%mpddInfo_sg%myrank
#endif
          END SELECT

          ! IF ('pres' .EQ. (TRIM(dX%fields(j)%Get_Name()))) THEN
          !   CALL dX%fields(j)%Set_Name('pres_ctl')
          !   dX%fields(j)%data = dX%fields(j)%data*100.0D0
          ! END IF
          ! IF ('qvapor' .EQ. (TRIM(dX%fields(j)%Get_Name()))) THEN
          !   CALL dX%fields(j)%Set_Name('qvapor_ctl')
          !   dX%fields(j)%data = dX%fields(j)%data*0.001D0

          ! DO i = 1, dX%sg%tSlots
          !   dX%fields(j)%data(:, :, i) = dX%fields(j)%data(: ,:, i)*dX%sg%s1
          ! END DO

        END IF
      END DO
    END DO
  END SUBROUTINE transAdjMultiply

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(X1)
    IMPLICIT NONE
    CLASS(Ctl2State_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: X1
    INTEGER(i_kind) :: i, j, k

    X1 = X
    CALL this%transFwdNonLinear(X1)
  END FUNCTION transFwdNonLinear_opr

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dX1)
    CLASS(Ctl2State_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(State_t):: dX1

    dX1 = this%transFwdNonLinear_opr(dX)
  END FUNCTION transFwdTanLinear_opr

  FUNCTION transAdjMultiply_opr(this, dX, X) RESULT(dX1)
    IMPLICIT NONE
    CLASS(Ctl2State_t) :: this
    TYPE(State_t), INTENT(IN) :: dX
    TYPE(State_t) :: dX1
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    dX1 = dX
    CALL this%transAdjMultiply(dX1)
  END FUNCTION transAdjMultiply_opr
END MODULE Ctl2State_m
