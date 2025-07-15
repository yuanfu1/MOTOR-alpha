!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.1
! HISTORY           :
!                     V 0.1 modified for X updated in time scale
!                           by Jiongming Pang (pang.j.m@hotmail.com), 2021-11-10, GBA-MWF, Shenzhen
!
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie: 2022/05/18 modified for vertical multigrid restriction
!!--------------------------------------------------------------------------------------------------

!> @brief
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
MODULE MGOpts_m
  USE State_m, ONLY: State_t
  USE Multigrid_m, ONLY: Multigrid_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpolateTime_m   ! For X updated in time scale. -Jiongming Pang
  USE YAMLRead_m
  IMPLICIT NONE

CONTAINS

  SUBROUTINE prolongationMG(XCoarser, XFiner, XbFiner, mg, IncFlag)
    TYPE(State_t) :: XCoarser, XFiner
    TYPE(State_t) :: XbFiner
    TYPE(Multigrid_t), INTENT(IN) :: mg
    LOGICAL, INTENT(IN) :: IncFlag
    INTEGER(i_kind) :: i, j, k, indCoarser

    ! Yuanfu Xie: 2022/05/18 added a multigrid in vertical and temporal directions:
    INTEGER(i_kind) :: vInc

    ! Setting vertical increment:
    IF (XFiner%sg%vLevel .EQ. XCoarser%sg%vLevel) THEN
      vInc = 1
    ELSE IF ((XFiner%sg%vLevel - 1) / 2 + 1 .EQ. XCoarser%sg%vLevel) THEN
      vInc = 2
    ELSE
      WRITE (*, 1) XFiner%sg%vLevel, XCoarser%sg%vLevel
1     FORMAT('prolongationMG -- error: vertical levels not suitable for a multigrid: ', 2I3)
      STOP
    END IF

    DO i = LBOUND(XFiner%fields, 1), UBOUND(XFiner%fields, 1)
      ! Horizontal and vertical interpretation, and time assigment
      DO k = 1, XFiner%sg%tSlots
        DO j = 1, XCoarser%sg%tSlots
          IF (ABS(XFiner%sg%tt(k) - XCoarser%sg%tt(j)) .LT. 1E-6) THEN
            PRINT *, 'In prolongation.'
            CALL mg%prolongation(XFiner%sg, XCoarser%sg, XFiner%Fields(i)%DATA(:, :, k), &
                                 XCoarser%Fields(i)%DATA(:, :, j), XCoarser%sg%vLevel, vInc, .TRUE.)
          END IF
        END DO
      END DO

    END DO

    ! Analysis interpolation for X in time scale.
    ! whole-interp way was added by Zilong Qin
    ! increment-interp way was added by Jiongming Pang

    ! whole-interp way
    IF (.NOT. IncFlag) THEN

      DO i = LBOUND(XFiner%fields, 1), UBOUND(XFiner%fields, 1)

        DO k = 1, XFiner%sg%tSlots
          BLOCK
            LOGICAL :: isPresent
            isPresent = .FALSE.
            indCoarser = 0 !PJM
            DO j = 1, XCoarser%sg%tSlots
              IF (ABS(XFiner%sg%tt(k) - XCoarser%sg%tt(j)) .LT. 1E-6) THEN
                isPresent = .TRUE.
                ! indCoarser = j !PJM
                EXIT !PJM
              END IF
            END DO
            IF (.NOT. isPresent) THEN
              XFiner%Fields(i)%DATA(:, :, k) = (XFiner%Fields(i)%DATA(:, :, k + 1) + XFiner%Fields(i)%DATA(:, :, k - 1)) / 2
            END IF
          END BLOCK
        END DO

      END DO

    END IF

    ! increment-interp way
    IF (IncFlag) THEN

      IF (XFiner%sg%tSlots .GT. XCoarser%sg%tSlots) THEN

        BLOCK
          TYPE(State_t) :: DX
          INTEGER(i_kind), ALLOCATABLE :: intp_ind_coarser(:)
          INTEGER(i_kind), ALLOCATABLE :: intp_ind_finer(:)
          TYPE(InterpolateData_Time_t), ALLOCATABLE :: DXIntpTime

          ! can be interpolated from 2 slots to 5 slots
          ALLOCATE (intp_ind_coarser(XCoarser%sg%tSlots))
          ALLOCATE (intp_ind_finer(XFiner%sg%tSlots))
          CALL IntpLinearIndex(XCoarser%sg%tSlots, XFiner%sg%tSlots, &
                               intp_ind_coarser, intp_ind_finer)

          DX = XFiner - XbFiner
          DO i = LBOUND(XFiner%fields, 1), UBOUND(XFiner%fields, 1)
            ALLOCATE (DXIntpTime)
            CALL DXIntpTime%IntpInitialize(DX%Fields(i)%DATA(:, :, intp_ind_coarser), intp_ind_coarser, intp_ind_finer)
            CALL DXIntpTime%IntpLinearTime(DX%Fields(i)%DATA(:, :, intp_ind_coarser))
            XFiner%Fields(i)%DATA = XbFiner%Fields(i)%DATA + DXIntpTime%data_intperpolate
            DEALLOCATE (DXIntpTime)
          END DO

          IF (ALLOCATED(intp_ind_coarser)) DEALLOCATE (intp_ind_coarser)
          IF (ALLOCATED(intp_ind_finer)) DEALLOCATE (intp_ind_finer)
        END BLOCK

        ! ! the original one, for reference
        ! DO i = LBOUND(XFiner%fields, 1), UBOUND(XFiner%fields, 1)

        !   DO k = 1, XFiner%sg%tSlots
        !     BLOCK
        !       LOGICAL :: isPresent
        !       isPresent  = .FALSE.
        !       DO j = 1, XCoarser%sg%tSlots
        !         IF (ABS(XFiner%sg%tt(k) - XCoarser%sg%tt(j)) .lt. 1e-6) THEN
        !           isPresent = .TRUE.
        !           EXIT !PJM
        !         END IF
        !       END DO
        !       IF (.not. isPresent) THEN
        !         XFiner%Fields(i)%data(:, :, k) = XbFiner%Fields(i)%data(:, :, k) + &
        !         (XFiner%Fields(i)%data(:, :, k + 1) - XbFiner%Fields(i)%data(:, :, k + 1)) / 2 + &
        !         (XFiner%Fields(i)%data(:, :, k - 1) - XbFiner%Fields(i)%data(:, :, k - 1)) / 2
        !       END IF
        !     END BLOCK
        !   END DO

        ! END DO

      END IF
    END IF

  END SUBROUTINE prolongationMG

  SUBROUTINE prolongationMGWhole(XCoarser, XFiner, sgCoarser, sgFiner, mg)
    TYPE(State_t) :: XCoarser, XFiner
    TYPE(State_t) :: XbFiner
    TYPE(Multigrid_t), INTENT(IN) :: mg
    TYPE(SingleGrid_t), INTENT(IN) :: sgCoarser, sgFiner
    INTEGER(i_kind) :: i, j, k
    REAL(r_kind) :: fake(1, 1)

    ! Yuanfu Xie: 2022/05/18 added a multigrid in vertical and temporal directions:
    INTEGER(i_kind) :: vInc

    ! Setting vertical increment:
    IF (sgFiner%vLevel .EQ. sgCoarser%vLevel) THEN
      vInc = 1
    ELSE IF ((sgFiner%vLevel - 1) / 2 + 1 .EQ. sgCoarser%vLevel) THEN
      vInc = 2
    ELSE
      WRITE (*, 1) sgFiner%vLevel, sgCoarser%vLevel
1     FORMAT('prolongationMG -- error: vertical levels not suitable for a multigrid: ', 2I3)
      STOP
    END IF

    IF (.NOT. sgFiner%isActiveProc()) RETURN

    DO i = LBOUND(XFiner%fields, 1), UBOUND(XFiner%fields, 1)
      ! Horizontal and vertical interpretation, and time assigment
      DO k = 1, sgFiner%tSlots
        DO j = 1, sgCoarser%tSlots
          IF (ABS(sgFiner%tt(k) - sgCoarser%tt(j)) .LT. 1E-6) THEN
            IF (sgCoarser%isActiveProc()) THEN
              CALL mg%prolongation(sgFiner, sgCoarser, XFiner%Fields(i)%DATA(:, :, k), &
                                   XCoarser%Fields(i)%DATA(:, :, j), sgCoarser%vLevel, vInc, .TRUE.)
            ELSE
              CALL mg%prolongation(sgFiner, sgCoarser, XFiner%Fields(i)%DATA(:, :, k), &
                                   fake, sgCoarser%vLevel, vInc, .TRUE.)
            END IF
          END IF
        END DO
      END DO
    END DO

    DO i = LBOUND(XFiner%fields, 1), UBOUND(XFiner%fields, 1)
      DO k = 1, sgFiner%tSlots
        BLOCK
          LOGICAL :: isPresent
          isPresent = .FALSE.
          DO j = 1, sgCoarser%tSlots
            IF (ABS(sgFiner%tt(k) - sgCoarser%tt(j)) .LT. 1E-6) THEN
              isPresent = .TRUE.
              EXIT
            END IF
          END DO
          IF (.NOT. isPresent) THEN
            XFiner%Fields(i)%DATA(:, :, k) = (XFiner%Fields(i)%DATA(:, :, k + 1) + XFiner%Fields(i)%DATA(:, :, k - 1)) / 2.0D0
          END IF
        END BLOCK
      END DO
    END DO

  END SUBROUTINE prolongationMGWhole

  ! The new propagation submodule - Zilong 20220526
  SUBROUTINE prolongationMGInc(XCoarser, XFiner, XbCoarser, XbFiner, sgCoarser, sgFiner, mg)
    TYPE(State_t) :: XCoarser, XFiner
    TYPE(State_t) :: XbFiner, XbCoarser
    TYPE(Multigrid_t), INTENT(IN) :: mg
    TYPE(SingleGrid_t), INTENT(IN) :: sgCoarser, sgFiner
    TYPE(State_t) :: DXCoarser, DXFiner

    IF (sgCoarser%isActiveProc()) DXCoarser = XCoarser - XbCoarser
    IF (sgFiner%isActiveProc()) DXFiner = XbFiner

    CALL prolongationMGWhole(DXCoarser, DXFiner, sgCoarser, sgFiner, mg)

    IF (sgFiner%isActiveProc()) XFiner = XbFiner + DXFiner

    PRINT *, 'In prolongation from g', XCoarser%sg%gLevel, ' to ', XFiner%sg%gLevel
  END SUBROUTINE prolongationMGInc

  ! The new propagation submodule - Zilong 20220526
  SUBROUTINE restrictionMGInc(XCoarser, XFiner, XbCoarser, XbFiner, sgCoarser, sgFiner, mg)
    TYPE(State_t) :: XCoarser, XFiner
    TYPE(State_t) :: XbFiner, XbCoarser
    TYPE(Multigrid_t), INTENT(IN) :: mg
    TYPE(SingleGrid_t), INTENT(IN) :: sgCoarser, sgFiner
    TYPE(State_t) :: DXCoarser, DXFiner

    IF (sgFiner%isActiveProc()) DXFiner = XFiner - XbFiner
    IF (sgCoarser%isActiveProc()) DXCoarser = XbCoarser%zeroCopy()

    CALL restrictionMGWhole(DXCoarser, DXFiner, sgCoarser, sgFiner, mg)

    IF (sgCoarser%isActiveProc()) XCoarser = DXCoarser + XbCoarser

    PRINT *, 'In restriction from g', XCoarser%sg%gLevel, ' to ', XFiner%sg%gLevel
  END SUBROUTINE restrictionMGInc

  SUBROUTINE restrictionMGWhole(XCoarser, XFiner, sgCoarser, sgFiner, mg)
    TYPE(State_t) :: XFiner
    TYPE(State_t) :: XCoarser
    TYPE(Multigrid_t), INTENT(IN) :: mg
    TYPE(SingleGrid_t), INTENT(IN) :: sgCoarser, sgFiner
    INTEGER(i_kind) :: i, j, k
    REAL(r_kind) :: fake(1, 1)

    ! Yuanfu Xie: 2022/05/18 modified for vertical multigrid restriction:
    INTEGER(i_kind) :: vInc

    IF (.NOT. sgFiner%isActiveProc()) RETURN

    DO i = LBOUND(XFiner%fields, 1), UBOUND(XFiner%fields, 1)

      ! Horizontal and vertical interpretation
      DO j = 1, sgCoarser%tSlots
        DO k = 1, sgFiner%tSlots

          ! This following print statement is a bad example for people hard to track it down!!!
          ! PRINT *, XFiner%sg%tt(k), XFiner%sg%tt(j), j, k, ABS(XFiner%sg%tt(k) - XCoarser%sg%tt(j))
          IF (ABS(sgFiner%tt(k) - sgCoarser%tt(j)) .LT. 1E-6) THEN
            IF (sgFiner%vLevel .EQ. sgCoarser%vLevel) THEN
              vInc = 1
            ELSE IF ((sgFiner%vLevel - 1) / 2 + 1 .EQ. sgCoarser%vLevel) THEN
              vInc = 2
            ELSE
              WRITE (*, 1) sgFiner%vLevel, sgCoarser%vLevel
1             FORMAT('restrictionMG -- error: vertical level is not suitable for a multigrid: ', 2I3)
              STOP
            END IF
            ! CALL mg%restriction(sgFiner, sgCoarser, XFiner%Fields(i)%data(:, :, k), &
            !                  XCoarser%Fields(i)%data(:, :, j), XCoarser%sg%vLevel, vInc, .TRUE.)

            IF (sgCoarser%isActiveProc()) THEN
              CALL mg%restriction(sgFiner, sgCoarser, XFiner%Fields(i)%DATA(:, :, k), &
                                  XCoarser%Fields(i)%DATA(:, :, j), sgCoarser%vLevel, vInc, .TRUE.)
            ELSE
              CALL mg%restriction(sgFiner, sgCoarser, XFiner%Fields(i)%DATA(:, :, k), &
                                  fake, sgCoarser%vLevel, vInc, .TRUE.)
            END IF
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE restrictionMGWhole

  SUBROUTINE restrictionMG(XCoarser, XFiner, mg)
    TYPE(State_t) :: XFiner
    TYPE(State_t) :: XCoarser
    TYPE(Multigrid_t), INTENT(IN) :: mg
    INTEGER(i_kind) :: i, j, k
    REAL(r_kind) :: fake(1, 1)

    ! Yuanfu Xie: 2022/05/18 modified for vertical multigrid restriction:
    INTEGER(i_kind) :: vInc

    DO i = LBOUND(XFiner%fields, 1), UBOUND(XFiner%fields, 1)

      ! Horizontal and vertical interpretation
      DO j = 1, XCoarser%sg%tSlots
        DO k = 1, XFiner%sg%tSlots
          ! This following print statement is a bad example for people hard to track it down!!!
          IF (ABS(XFiner%sg%tt(k) - XCoarser%sg%tt(j)) .LT. 1E-6) THEN
            IF (XFiner%sg%vLevel .EQ. XCoarser%sg%vLevel) THEN
              vInc = 1
            ELSE IF ((XFiner%sg%vLevel - 1) / 2 + 1 .EQ. XCoarser%sg%vLevel) THEN
              vInc = 2
            ELSE
              WRITE (*, 1) XFiner%sg%vLevel, XCoarser%sg%vLevel
1             FORMAT('restrictionMG -- error: vertical level is not suitable for a multigrid: ', 2I3)
              STOP
            END IF
            IF (XCoarser%sg%isActiveProc()) THEN
              CALL mg%restriction(XFiner%sg, XCoarser%sg, XFiner%Fields(i)%DATA(:, :, k), &
                                  XCoarser%Fields(i)%DATA(:, :, j), XCoarser%sg%vLevel, vInc, .TRUE.)
            ELSE
              CALL mg%restriction(XFiner%sg, XCoarser%sg, XFiner%Fields(i)%DATA(:, :, k), &
                                  fake, XCoarser%sg%vLevel, vInc, .TRUE.)
            END IF
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE restrictionMG

  SUBROUTINE update_restrictionOfStatics(configFile, X, sg)
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), INTENT(IN) :: X
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg

    CALL Calc_Scale_TProf(configFile, X, sg)
    ! CALL Calc_Scale_RhovProf(configFile, X, sg)
    CALL Calc_Scale_QvProf(configFile, X, sg)
    ! CALL Calc_sg_vars(configFile, X, sgCoarser)

  END SUBROUTINE update_restrictionOfStatics

  SUBROUTINE Calc_sg_vars(configFile, X, sg)
    ! Prepare pres and saved in sg
    USE Pres2Rho_m, ONLY: Pres2Rho_t
    TYPE(State_t), INTENT(IN) :: X
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    TYPE(State_t) :: XTemp
    TYPE(Pres2Rho_t) :: Pres2Rho
    INTEGER :: j

    XTemp = X
    !Pres2Rho = Pres2Rho_t(configFile, X)

    IF (XTemp%getVarIdx(TRIM('pres')) .NE. 0) THEN
      sg%FGPres = XTemp%fields(XTemp%getVarIdx(TRIM('pres')))%DATA
      sg%psfc = sg%FGPres(1, :, :)
      !ELSE
      !  IF (XTemp%getVarIdx(TRIM('rho')) .NE. 0 .AND. XTemp%getVarIdx(TRIM('pres')) .EQ. 0 ) &
      !  CALL Pres2Rho%fwdNL(XTemp) ! Get pres
    END IF
    IF (XTemp%getVarIdx(TRIM('temp')) .NE. 0) THEN
      sg%tskin = XTemp%fields(XTemp%getVarIdx(TRIM('temp')))%DATA(1, :, :)
    END IF
    IF (XTemp%getVarIdx(TRIM('uwnd')) .NE. 0) THEN
      sg%u10m = XTemp%fields(XTemp%getVarIdx(TRIM('uwnd')))%DATA(1, :, :)
    END IF
    IF (XTemp%getVarIdx(TRIM('vwnd')) .NE. 0) THEN
      sg%v10m = XTemp%fields(XTemp%getVarIdx(TRIM('vwnd')))%DATA(1, :, :)
    END IF

  END SUBROUTINE Calc_sg_vars

  SUBROUTINE Calc_Scale_RhovProf(configFile, X, sg)
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), INTENT(IN) :: X
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    CHARACTER(len=20) :: rhov_scale_scheme
    INTEGER :: ifile

    ifile = yaml_get_var(configFile, 'CV_Transform', 'rhov_scale_scheme', rhov_scale_scheme)

    IF (TRIM(rhov_scale_scheme) .NE. 'None') THEN

      SELECT CASE (TRIM(rhov_scale_scheme))
      CASE ("DomainAvg")
        !PRINT *, 'USE SCHEME ', TRIM(rhov_scale_scheme), 'max rhov_ctl = ', MAXVAL(X%Fields(X%getVarIdx('rhov_ctl'))%data)
        ! IF(MAXVAL(X%Fields(X%getVarIdx('rhov_ctl'))%data) > 0.1) THEN
        CALL RhovProf_DomainAvg(configFile, X, sg)
      CASE ("EXPO")
        !PRINT *, 'USE SCHEME ', TRIM(rhov_scale_scheme)
        CALL RhovProf_EXPO(configFile, X, sg)
      CASE ("Power")
        !PRINT *, 'USE SCHEME ', TRIM(rhov_scale_scheme)
        sg%SVapor = 1.0D0
      CASE default
        PRINT *, 'No rhov transform scheme found'
      END SELECT
      !IF (PRESENT(sgFiner)) PRINT *, 'sgFiner%SVapor = ', sgFiner%SVapor1D
      !PRINT *, 'rhov: sg%SVapor = ', sg%SVapor(::5,1)

    ELSE
      sg%SVapor = 1.0D0
    END IF

  END SUBROUTINE Calc_Scale_RhovProf

  SUBROUTINE Calc_Scale_QvProf(configFile, X, sg)
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), INTENT(IN) :: X
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    CHARACTER(len=20) :: qvapor_scale_scheme
    INTEGER :: ifile, i

    ifile = yaml_get_var(configFile, 'CV_Transform', 'qvapor_scale_scheme', qvapor_scale_scheme)

    IF (TRIM(qvapor_scale_scheme) .NE. 'None') THEN

      SELECT CASE (TRIM(qvapor_scale_scheme))
      CASE ("DomainAvg")
        !PRINT *, 'USE SCHEME ', TRIM(qvapor_scale_scheme), 'max qvapor = ', MAXVAL(X%Fields(X%getVarIdx('qvapor'))%data)
        CALL QvProf_DomainAvg(configFile, X, sg)
      CASE ("Power")
        !PRINT *, 'USE SCHEME ', TRIM(qvapor_scale_scheme)
        DO i = 1, SIZE(sg%SVapor, 2)
          sg%SVapor(:, i) = sg%s_qvapor(:)
        END DO
        PRINT *, "In power of Calc_Scale_QvProf"
      CASE default
        PRINT *, 'No rhov transform scheme found'
      END SELECT
      !IF (PRESENT(sgFiner)) PRINT *, 'sgFiner%SVapor = ', sgFiner%SVapor1D
      !PRINT *, 'qvapor: sg%SVapor = ', sg%SVapor(::5,1)

    ELSE
      sg%SVapor = 1.0D0
    END IF

  END SUBROUTINE Calc_Scale_QvProf

  SUBROUTINE Calc_Scale_TProf(configFile, X, sg)
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), INTENT(IN) :: X
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    CHARACTER(len=20) :: temp_scale_scheme
    INTEGER :: ifile

    ifile = yaml_get_var(configFile, 'CV_Transform', 'temp_scale_scheme', temp_scale_scheme)

    IF (TRIM(temp_scale_scheme) .NE. 'None') THEN

      SELECT CASE (TRIM(temp_scale_scheme))
      CASE ("DomainAvg")
        !PRINT *, 'USE SCHEME ', TRIM(temp_scale_scheme), 'max temp = ', MAXVAL(X%Fields(X%getVarIdx('temp'))%data)
        ! IF(MAXVAL(X%Fields(X%getVarIdx('temp_ctl'))%data) < 100.0) THEN
        CALL TProf_DomainAvg(configFile, X, sg)
          !PRINT *, 'sg%STemp of xm ', sg%STemp(:, 1)
      CASE ("Power")
        !PRINT *, 'USE SCHEME ', TRIM(temp_scale_scheme)
        sg%STemp = 1.0D0
        Print*, "In power of Calc_Scale_QvProf"
      CASE default
        PRINT *, 'No temp transform scheme found'
      END SELECT
      !IF (PRESENT(sgFiner))  PRINT *, 'sgFiner%STemp = ', sgFiner%STemp1D
      !PRINT *, 'sg%STemp = ', sg%STemp(::5,1)

    ELSE
      sg%STemp = 1.0D0
    END IF

  END SUBROUTINE Calc_Scale_TProf

  SUBROUTINE RhovProf_EXPO(configFile, X, sg)
    USE Pres2Rho_m, ONLY: Pres2Rho_t
    USE Qvapor2Rhov_m, ONLY: Qvapor2Rhov_t
    USE RhoVCtl2RhoV_sv_m, ONLY: RhoVCtl2RhoV_sv_t
    USE QVCtl2QV_sv_m, ONLY: QVCtl2QV_sv_t
    TYPE(State_t), INTENT(IN) :: X
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    TYPE(State_t) :: XTemp
    TYPE(Pres2Rho_t) :: Pres2Rho
    TYPE(QVCtl2QV_sv_t) :: QVCtl2QV
    TYPE(Qvapor2Rhov_t) :: Qvapor2Rhov
    TYPE(RhoVCtl2RhoV_sv_t) :: RhoVCtl2RhoV_sv

    INTEGER(i_kind) :: ih, iv
    REAL(r_kind) :: temp_s, beta
    REAL(r_kind), ALLOCATABLE :: temp1D(:)

    XTemp = X
    RhoVCtl2RhoV_sv = RhoVCtl2RhoV_sv_t(configFile, X)
    QVCtl2QV = QVCtl2QV_sv_t(configFile, X)
    Qvapor2Rhov = Qvapor2Rhov_t(configFile, X)
    Pres2Rho = Pres2Rho_t(configFile, X)

    IF (XTemp%sg%isActiveProc() .AND. XTemp%getVarIdx(TRIM('rho')) .NE. 0) THEN
      ! For sg < meEnd, rhov_ctl is the real rhov_ctl, of which the value is always greater than 0.1
      ! For sg = mgEnd, rhov_ctl is actually qvapor, of which the value is always smaller than 0.1
      IF (XTemp%getVarIdx(TRIM('rhov')) .NE. 0) THEN
        IF (MAXVAL(XTemp%Fields(XTemp%getVarIdx('rhov'))%DATA) > 0.1) THEN
          ! IF(MAXVAL(XTemp%Fields(XTemp%getVarIdx('rhov_ctl'))%data) > 30.0) THEN
          CALL Qvapor2Rhov%fwdNL(XTemp)  ! Get d qvapor
          temp1D = XTemp%Fields(XTemp%getVarIdx('qvapor'))%DATA(1, :, 1)
        ELSE
          temp1D = X%Fields(X%getVarIdx('rhov'))%DATA(1, :, 1)
        END IF
      END IF

      IF (XTemp%getVarIdx(TRIM('qvapor')) .NE. 0) THEN
        IF (MAXVAL(XTemp%Fields(XTemp%getVarIdx('qvapor'))%DATA) > 0.1) THEN
          CALL QVCtl2QV%fwdNL(XTemp) ! get d rhov
          temp1D = XTemp%Fields(XTemp%getVarIdx('qvapor'))%DATA(1, :, 1)
        ELSE
          temp1D = X%Fields(X%getVarIdx('qvapor_ctl'))%DATA(1, :, 1)
        END IF
      END IF

      DO ih = 1, SIZE(temp1D, 1)

        temp_s = temp1D(ih) * XTemp%Fields(XTemp%getVarIdx('rho'))%DATA(1, ih, 1)

        ! beta=1.0/(sg%ztop-sg%topo(ih))
        ! sg%SVapor(:, ih) = temp_s * 1.5D0**(-beta * (sg%zHght(:, ih)-sg%topo(ih)))
        beta = 0.0005
        sg%SVapor(:, ih) = temp_s * EXP(-beta * (sg%zHght(:, ih) - sg%topo(ih)))
        ! print *, 'check scaling: ', temp_s, sg%zHght(1, ih), sg%topo(ih), sg%SVapor(1, ih)

      END DO

      DEALLOCATE (temp1D)
      WHERE (sg%SVapor < 1.0E-12) sg%SVapor = 1.0E-12
      ! CALL sg%mpddInfo_sg%bCast(sg%SVapor)
    END IF

  END SUBROUTINE RhovProf_EXPO

  SUBROUTINE RhovProf_DomainAvg(configFile, X, sg)

    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), INTENT(IN) :: X
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg

    INTEGER(i_kind) :: ih, iv, it, j
    REAL(r_kind), ALLOCATABLE :: temp1D(:), temp3D(:, :, :)
    REAL(r_kind) :: sum1, sum2
    INTEGER(i_kind) :: num1, num2

    DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      !PRINT *, 'Check RhovProf_DomainAvg varnames: ', TRIM(X%fields(j)%Get_Name())
      !PRINT *, TRIM(X%fields(j)%Get_Name()), MAXVAL(X%fields(j)%data), MINVAL(X%fields(j)%data)
    END DO

    IF (X%getVarIdx(TRIM('rhov')) .NE. 0) THEN
      ! ALLOCATE(temp3D(sg%vLevel, SIZE(X%Fields(X%getVarIdx('rhov'))%data,2), sg%tSlots))
      ! temp3D = X%Fields(X%getVarIdx('rhov'))%data
      ! ALLOCATE(temp1D(SIZE(temp3D,2)))

      ! IF (X%sg%isBaseProc()) THEN

      ALLOCATE (temp3D(sg%vLevel, SIZE(X%Fields(X%getVarIdx('rhov'))%DATA, 2), sg%tSlots))
      temp3D = X%Fields(X%getVarIdx('rhov'))%DATA
      ALLOCATE (temp1D(SIZE(temp3D, 2)))

      DO iv = 1, sg%vLevel

        ! Actually rhov is qvapor here
        ! Use the fields of 1st time for all times
        temp1D(:) = temp3D(iv, :, 1) * X%Fields(X%getVarIdx('rho'))%DATA(iv, :, 1)
        sum1 = SUM(temp1D)
        num1 = SIZE(temp1D, 1)
        CALL X%mpddSub%AllReduceSumReal(sum1, sum2)
        CALL X%mpddSub%AllReduceSumInt(num1, num2)
        sg%SVapor1D(iv) = sum2 / REAL(num2)
        sg%SVapor1D(iv) = sg%SVapor1D(iv) / 100.0D0
        sg%SVapor(iv, :) = sg%SVapor1D(iv)

      END DO
      DEALLOCATE (temp1D, temp3D)

      ! END IF

      ! DEALLOCATE(temp1D, temp3D)
    END IF
    ! CALL sg%mpddInfo_sg%bCast(sg%SVapor)
    ! CALL sg%mpddInfo_sg%bCast(sg%SVapor1D)

    ! DO it = 1, sg%tSlots
    !   DO iv = 1, sg%vLevel

    !     ! Actually rhov is qvapor here
    !     temp1D(:) = temp3D(iv,:,it) * X%Fields(X%getVarIdx('rho'))%data(iv,:,it)
    !     sg%SVapor(iv, :) = sqrt((temp1D(:) - sg%SVapor1D(iv))**2)
    !     WHERE (sg%SVapor(iv, :) .LE. 1.0E-5) sg%SVapor(iv, :) = 1.0E-5

    !   END DO
    ! END DO
    ! DEALLOCATE(temp1D, temp3D)

    ! PRINT *, 'check sg%SVapor: ', sg%SVapor(:,1)
    ! PRINT *, 'check max/min of sg%SVapor: ', MAXVAL(sg%SVapor), MINVAL(sg%SVapor)

  END SUBROUTINE RhovProf_DomainAvg

  SUBROUTINE QvProf_DomainAvg(configFile, X, sg)

    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), INTENT(IN) :: X
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    REAL(r_kind) :: sum1, sum2
    INTEGER(i_kind) :: num1, num2

    INTEGER(i_kind) :: ih, iv, it, j
    REAL(r_kind), ALLOCATABLE :: temp1D(:), temp3D(:, :, :)

    !DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
    !PRINT *, 'Check QvProf_DomainAvg varnames: ', TRIM(X%fields(j)%Get_Name())
    !PRINT *, TRIM(X%fields(j)%Get_Name()), MAXVAL(X%fields(j)%data), MINVAL(X%fields(j)%data)
    !END DO

    IF (X%getVarIdx(TRIM('qvapor')) .NE. 0) THEN
      ! ALLOCATE(temp3D(sg%vLevel, SIZE(X%Fields(X%getVarIdx('qvapor'))%data,2), sg%tSlots))
      ! temp3D = X%Fields(X%getVarIdx('qvapor'))%data
      ! ALLOCATE(temp1D(SIZE(temp3D,2)))

      ! IF (X%sg%isBaseProc()) THEN

      ALLOCATE (temp3D(sg%vLevel, SIZE(X%Fields(X%getVarIdx('qvapor'))%DATA, 2), sg%tSlots))
      temp3D = X%Fields(X%getVarIdx('qvapor'))%DATA
      ALLOCATE (temp1D(SIZE(temp3D, 2)))

      DO iv = 1, sg%vLevel

        temp1D(:) = temp3D(iv, :, 1)  ! Use the fields of 1st time for all times
        sum1 = SUM(temp1D)
        num1 = SIZE(temp1D, 1)
        CALL X%mpddSub%AllReduceSumReal(sum1, sum2)
        CALL X%mpddSub%AllReduceSumInt(num1, num2)
        sg%SVapor1D(iv) = sum2 / REAL(num2)
        sg%SVapor1D(iv) = 0.001D0
        ! sg%SVapor1D(iv) = sg%s_qvapor(iv)
        ! sg%SVapor1D(iv) = sg%SVapor1D(iv) / 10.0D0 ! for qvapor=0.02kg/kg, variance~=0.001kg/kg, 1/20 times
        sg%SVapor(iv, :) = sg%SVapor1D(iv)
        ! sg%SVapor(iv, :) = sg%s1(iv, :)

      END DO

      DEALLOCATE (temp1D, temp3D)

      ! END IF

      ! DEALLOCATE(temp1D, temp3D)
    END IF
    ! CALL sg%mpddInfo_sg%bCast(sg%SVapor)
    ! CALL sg%mpddInfo_sg%bCast(sg%SVapor1D)

    ! DO it = 1, sg%tSlots
    !   DO iv = 1, sg%vLevel

    !     temp1D(:) = temp3D(iv,:,it)
    !     sg%SVapor(iv, :) = sqrt((temp1D(:) - sg%SVapor1D(iv))**2)
    !     WHERE (sg%SVapor(iv, :) .LE. 1.0E-5) sg%SVapor(iv, :) = 1.0E-5

    !   END DO
    ! END DO
    ! DEALLOCATE(temp1D, temp3D)

    ! PRINT *, 'check sg%SVapor: ', sg%SVapor(:,1)
    ! PRINT *, 'check max/min of sg%SVapor: ', MAXVAL(sg%SVapor), MINVAL(sg%SVapor)

  END SUBROUTINE QvProf_DomainAvg

  SUBROUTINE TProf_DomainAvg(configFile, X, sg)
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), INTENT(IN) :: X
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    REAL(r_kind) :: sum1, sum2
    INTEGER(i_kind) :: num1, num2

    INTEGER(i_kind) :: ih, iv, it
    REAL(r_kind), ALLOCATABLE :: temp1D(:), temp3D(:, :, :)

    !DO j = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
    !PRINT *, 'Check TProf_DomainAvg varnames: ', TRIM(X%fields(j)%Get_Name())
    !PRINT *, TRIM(X%fields(j)%Get_Name()), MAXVAL(X%fields(j)%data), MINVAL(X%fields(j)%data)
    !END DO

    IF (X%getVarIdx(TRIM('temp')) .NE. 0) THEN
      ! ALLOCATE(temp3D(sg%vLevel, SIZE(X%Fields(X%getVarIdx('temp'))%data,2), sg%tSlots))
      ! temp3D = X%Fields(X%getVarIdx('temp'))%data
      ! ALLOCATE(temp1D(SIZE(temp3D,2)))

      ! IF (X%sg%isBaseProc()) THEN

      ALLOCATE (temp3D(sg%vLevel, SIZE(X%Fields(X%getVarIdx('temp'))%DATA, 2), sg%tSlots))
      temp3D = X%Fields(X%getVarIdx('temp'))%DATA
      ALLOCATE (temp1D(SIZE(temp3D, 2)))
      !PRINT *, 'check temp ', maxval(temp3D), minval(temp3D)

      DO iv = 1, sg%vLevel

        temp1D(:) = temp3D(iv, :, 1) ! Use the fields of 1st time for all times
        sum1 = SUM(temp1D)
        num1 = SIZE(temp1D, 1)
        CALL X%mpddSub%AllReduceSumReal(sum1, sum2)
        CALL X%mpddSub%AllReduceSumInt(num1, num2)
        sg%STemp1D(iv) = sum2 / REAL(num2)
        sg%STemp1D(iv) = 1.0D0
        sg%STemp1D(iv) = sg%STemp1D(iv) / 1.0D0
        sg%STemp(iv, :) = sg%STemp1D(iv)
        ! PRINT *, 'sg%STemp 1 ', sg%STemp(iv, 1)

      END DO
      DEALLOCATE (temp1D, temp3D)

      ! END IF

      ! DEALLOCATE(temp1D, temp3D)
    END IF
    ! CALL sg%mpddInfo_sg%bCast(sg%STemp)
    ! CALL sg%mpddInfo_sg%bCast(sg%STemp1D)
    ! DO it = 1, sg%tSlots
    !   DO iv = 1, sg%vLevel
    !     temp1D(:) = temp3D(iv,:,it)
    !     sg%STemp(iv, :) = sqrt((temp1D(:) - sg%STemp1D(iv))**2)
    !     WHERE (sg%STemp(iv, :) .LE. 1.0E-5) sg%STemp(iv, :) = 1.0E-5
    !   END DO
    ! END DO
    ! PRINT *, 'sg%STemp 1 ', sg%STemp(:, 1)

    ! DEALLOCATE(temp1D, temp3D)
    ! PRINT *, 'sg%STemp 2 ', sg%STemp(::5, 1)

  END SUBROUTINE TProf_DomainAvg

  SUBROUTINE update_FGPres(XA, Xctl)
    USE parameters_m, ONLY: dry_air_gas_const
    TYPE(State_t), INTENT(IN) :: XA
    TYPE(State_t), INTENT(INOUT) :: Xctl
    REAL(r_kind), ALLOCATABLE :: swapDA(:, :, :)
    INTEGER :: i, j

    ALLOCATE (swapDA(XA%sg%vLevel, XA%sg%num_cell, XA%sg%tSlots))

    ! IF (XA%getVarIdx(TRIM('pres')) .NE. 0) THEN
    !   XA%sg%FGPres = XA%fields(XA%getVarIdx(TRIM('pres')))%data
    ! END IF

    IF (XA%getVarIdx('rho') .NE. 0) THEN
      swapDA(1, :, :) = XA%Fields(XA%getVarIdx('rho'))%DATA(1, :, :) * &
                        XA%Fields(XA%getVarIdx('temp'))%DATA(1, :, :) * dry_air_gas_const
    ELSE IF (XA%getVarIdx('pres') .NE. 0) THEN
      swapDA(1, :, :) = XA%Fields(XA%getVarIdx('pres'))%DATA(1, :, :)
    END IF

    DO i = 1, XA%sg%tSlots
      DO j = 2, XA%sg%vLevel
        swapDA(j, :, i) = swapDA(j - 1, :, i) * EXP(-(XA%sg%zHght(j, :) - XA%sg%zHght(j - 1, :)) * 9.80665D0 &
                                                    / (dry_air_gas_const * &
                                                       (XA%Fields(XA%getVarIdx('temp'))%DATA(j, :, i) * &
                                                        (XA%Fields(XA%getVarIdx('qvapor'))%DATA(j, :, i) * 0.608 + 1) + &
                                                        XA%Fields(XA%getVarIdx('temp'))%DATA(j - 1, :, i) * &
                                                        (XA%Fields(XA%getVarIdx('qvapor'))%DATA(j - 1, :, i) * 0.608 + 1)) / 2.0D0))
      END DO
    END DO
    Xctl%sg%FGPres = swapDA

    DEALLOCATE (swapDA)

  END SUBROUTINE

END MODULE MGOpts_m
