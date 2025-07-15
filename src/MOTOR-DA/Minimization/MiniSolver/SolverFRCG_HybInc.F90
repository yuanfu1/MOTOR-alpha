!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Tu Sanshan, 2021/1/26, @GBA-MWF, Shenzhen
!   Modified by Zilong Qin (zilong.qin@gmail.com) and Yuanfu Xie, 2022/3/13, @GBA-MWF, Shenzhen
!     Add constraints in the Minimization, change the search strategy of step length, from recursively searching to only one try (1/2).
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022-10-08, @GBA-MWF, Shenzhen
!     Add output of components of Jo for debugging purpose.
!!--------------------------------------------------------------------------------------------------

!> @brief
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
MODULE SolverFRCG_HybInc_m
  USE State_m, ONLY: State_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE JFunc_HybInc_m, ONLY: JFunc_HybInc_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE MiniSolver_m, ONLY: MiniSolver_t
  USE State2NC_m
  ! #define TRACK_DEBUG_INFO

  TYPE, EXTENDS(MiniSolver_t) :: SolverFRCG_HybInc_t
    LOGICAL :: HybInc
    INTEGER(i_kind) :: ensNum

  CONTAINS
    PROCEDURE, PUBLIC :: initialize
    FINAL :: destructor
    ! PROCEDURE, PUBLIC, PASS :: run
    PROCEDURE, PUBLIC, PASS :: run_HyyInc

  END TYPE

CONTAINS

  SUBROUTINE initialize(this, configFile)
    !USE NMLRead_m
    USE YAMLRead_m
    IMPLICIT NONE
    CLASS(SolverFRCG_HybInc_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    INTEGER(i_kind) :: istatus

    CALL this%MiniSolver_t%initialize(configFile)

    istatus = yaml_get_var(TRIM(configFile), 'BMat', 'HybInc', this%HybInc)
    IF (istatus .NE. 0) THEN
      this%HybInc = .FALSE.
    END IF

    istatus = yaml_get_var(TRIM(configFile), 'BMat', 'ensNum', this%ensNum)
    IF (istatus .NE. 0) THEN
      this%ensNum = 0
    END IF

  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(SolverFRCG_HybInc_t), INTENT(INOUT) :: this

  END SUBROUTINE destructor

  SUBROUTINE run_HyyInc(this, X, JFunc, sg, iters)
    IMPLICIT NONE
    CLASS(SolverFRCG_HybInc_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    TYPE(JFunc_HybInc_t), INTENT(INOUT) :: JFunc
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg
    INTEGER(i_kind), OPTIONAL :: iters
    INTEGER    :: i, j, k, nSwap, loopIdx, tIdx, hIdx, vIdx, lineSearches
    INTEGER    :: n = 25, m = 5, iprint = 1, MaxOptStep, maxSearch = 20 ! Yuanfu Xie added maxSearch controlling the number of lin searches.

    REAL(r_kind), PARAMETER    :: factr = 1.0D+7, pgtol = 1.0D-5
    REAL(r_kind) :: maxdRng
    REAL(r_kind)               :: f_min, rng_a, rng_b, lamda, miu, f_rng_b, f_lamda, f_miu, f_rng_b_last
    REAL(r_kind)               :: beta, dotg0g0, dotg1g1, fd

    TYPE(State_t) :: X_Temp, g0, d0, X_fullstate
    LOGICAL :: isTouchBound, isReCalJAndGradJ, isTouchBoundAll

    TYPE(State_t), ALLOCATABLE :: dXe(:), d0e(:), dXe_Temp(:), g0e(:), dXe2(:)
    TYPE(State_t) :: dXs, dXs_Temp, Aw_tmp, DAw_tmp, Aw_tmp2, DAw_tmp2
    INTEGER       :: ii, t, idmiddle, idmiddle2

    CHARACTER(LEN=1)  :: ensID

    ! Normalizing factors:
    REAL(r_kind), ALLOCATABLE :: costsJo(:,:)
    LOGICAL :: ifdebug = .TRUE.

    !     Declare a few additional variables for this sample problem
    IF (.NOT. sg%isActiveProc()) RETURN

    PRINT *, 'Start FR-CG minimization...', sg%gLevel, sg%mpddInfo_sg%myrank

    idmiddle = sg%dimCell_global(1) * sg%dimCell_global(2) / 2 + sg%dimCell_global(2) / 2
    idmiddle2 = sg%dimCell_global(1) * sg%dimCell_global(2) / 2 + sg%dimCell_global(2) / 2 + sg%dimCell_global(2) / 6

    n = 0
    DO i = LBOUND(X%fields, 1), UBOUND(X%fields, 1)
      n = n + sg%num_icell_global * sg%vLevel * sg%tSlots
    END DO

    WRITE (*, 1) n, sg%mpddInfo_sg%myrank, sg%gLevel
1   FORMAT('Minimization - n is: ', I12, ' myrank is: ', I2, ' Glevel: ', I2)

    IF (sg%isBaseProc()) PRINT *, 'Solving the simple variational analysis by FR-CG ... at ', sg%gLevel

    IF (this%HybInc .AND. this%ensNum > 0) THEN
      IF (ALLOCATED(dXe)) DEALLOCATE (dXe)
      ALLOCATE (dXe(this%ensNum), d0e(this%ensNum))
      IF (ALLOCATED(dXe_Temp)) DEALLOCATE (dXe_Temp)
      ALLOCATE (dXe_Temp(this%ensNum))
      IF (ALLOCATED(g0e)) DEALLOCATE (g0e)
      ALLOCATE (g0e(this%ensNum))

      DO i = 1, this%ensNum
        dXe(i) = X%zeroCopy()
        d0e(i) = X%zeroCopy()
        g0e(i) = X%zeroCopy()
        dXe_Temp(i) = X%zeroCopy()
      END DO
      dXs = X%zeroCopy()
    END IF

    IF (this%HybInc) THEN
      CALL JFunc%get_J_and_grad_J_vec(X, f_min, g0, dXs, dXe, g0e)
    ELSE
      CALL JFunc%get_J_and_grad_J_vec(X, f_min, g0)
    END IF

    IF (this%HybInc) THEN
      dotg0g0 = (g0.DOT.g0)
      d0 = g0 * (-1.0D0)
      DO i = 1, this%ensNum
        dotg0g0 = dotg0g0 + (g0e(i) .DOT.g0e(i))
        d0e(i) = g0e(i) * (-1.0D0)
      END DO
    ELSE
      dotg0g0 = (g0.DOT.g0)
      d0 = g0 * (-1.0D0)
    END IF

    IF (sg%isBaseProc()) PRINT *, 'step', 0, '|g|: ', SQRT(dotg0g0), 'f_min', f_min, 'g', &
      SIZE(g0%fields(1)%DATA)

    loopIdx = 1

    MaxOptStep = this%MaxOptStep
    IF (PRESENT(iters)) THEN
      IF (iters > MaxOptStep) MaxOptStep = iters
    END IF

    ! Initialize cost components of Jo:
    ALLOCATE (costsJo(MaxOptStep, 50))
    costsJo = 0.0D0

    IF (sg%isBaseProc()) PRINT *, "MaxOptStep: ", MaxOptStep
!   The beginning of the loop
    DO WHILE (dotg0g0 > pgtol .AND. loopIdx < MaxOptStep)
      isReCalJAndGradJ = .FALSE.

!=========== compute lamda begin =============
      IF (this%HybInc) THEN
        CALL JFunc%get_J0_value(d0, X, fd, d0e)
      ELSE
        CALL JFunc%get_J0_value(d0, X, fd)
      END IF
      ! PRINT *, 'PJM-DEBUG --- JFunc%get_J0_value --- fd:', fd
      rng_b = dotg0g0 / fd * 0.5D0
!=========== compute lamda end  =============

      IF (this%HybInc) THEN
        dXs_Temp = dXs + (d0 * rng_b)
        DO i = 1, this%ensNum
          dXe_Temp(i) = dXe(i) + (d0e(i) * rng_b)
        END DO
      ELSE
        X_Temp = X + (d0*rng_b)
        X_fullstate = X_Temp
      END IF
      CALL JFunc%get_fullstate_CtlVar(X_fullstate)

      ! =========== compute the boundary =============
      ! Add constraints for qrain > 0, pres >0 and rho >0
      ! isTouchBound = .FALSE.
      IF (.NOT. this%HybInc) THEN

      DO i = 1, SIZE(X_fullstate%Fields, 1)
        IF (TRIM(X_fullstate%Fields(i)%Get_Name()) == 'rho_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'rho' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'rhov' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'rhov_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qrain' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qvapor' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qvapor_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'rhor' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'pcpa' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'pcpa5min' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'ref' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'temp' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'rhor_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qcloud_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qice' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qice_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qrain' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qrain_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qsnow' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qsnow_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qgraupel' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qgraupel_ctl') THEN
            FORALL (tIdx=1:sg%tSlots, vIdx=1:sg%vLevel, hIdx=1:sg%num_cell, &
            X_fullstate%Fields(i)%data(vIdx, hIdx, tIdx) <= 0.0D0)
              d0%Fields(i)%data(vIdx, hIdx, tIdx) = -X%Fields(i)%data(vIdx, hIdx, tIdx)/rng_b/2.0D0
              X_Temp%Fields(i)%data(vIdx, hIdx, tIdx) = X%Fields(i)%data(vIdx, hIdx, tIdx)/2.0D0
              ! isTouchBound = .TRUE.
            END FORALL    
        END IF
      END DO

      END IF
      ! CALL sg%mpddInfo_sg%AllRedLOR(isTouchBound, isTouchBoundAll)
      ! isReCalJAndGradJ = isReCalJAndGradJ .OR. isTouchBoundAll
      ! IF(isTouchBound) X_Temp = X + (d0*rng_b)
!=========== compute the boundary =============

      IF (this%HybInc) THEN
        CALL JFunc%get_J_and_grad_J_vec(X, f_rng_b, g0, dXs_Temp, dXe_Temp, g0e)
      ELSE
        CALL JFunc%get_J_and_grad_J_vec(X_Temp, f_rng_b, g0)
      END IF
      IF (sg%isBaseProc()) PRINT *, 'rng_b: ', rng_b, ' f_rng_b:', f_rng_b


      IF (f_rng_b > f_min) THEN
        f_rng_b_last = f_rng_b
        lineSearches = 0
        DO WHILE (f_rng_b > f_min)
          f_rng_b_last = f_rng_b
          rng_b = rng_b * 0.618D0

          IF (this%HybInc) THEN
            dXs_Temp = dXs + (d0 * (rng_b))
            DO i = 1, this%ensNum
              dXe_Temp(i) = dXe(i) + (d0e(i) * (rng_b))
            END DO
            CALL JFunc%get_J_value(X_Temp, f_rng_b, dXs_Temp, dXe_Temp)
          ELSE
            X_Temp = X + (d0 * (rng_b))
            CALL JFunc%get_J_value(X_Temp, f_rng_b)
          END IF

          ! Yuanfu Xie added Glevel for for clarity: Cover it up for saving time 2023/10/22
!           IF (sg%isBaseProc()) WRITE(*,3) rng_b, f_rng_b, f_rng_b_last,sg%gLevel
! 3         FORMAT('LineSearch - New rng_b: ',E14.6,' f_rng_b: ',E14.6,' f_rng_b_last: ',E14.6,' Glvl:',I3)
          lineSearches = lineSearches + 1
          IF (lineSearches .GT. maxSearch) EXIT ! Take the current iteration as the approximate solution
        END DO
        isReCalJAndGradJ = .TRUE.
      END IF

      ! IF (f_rng_b > f_min) THEN
      !   rng_a = 0.0
      !   miu = rng_a + 0.618*(rng_b - rng_a)
      !   X_Temp = X + (d0*miu)
      !   CALL JFunc%get_J_value(X_Temp, f_miu)

      !   lamda = rng_a + 0.382*(rng_b - rng_a)
      !   X_Temp = X + (d0*lamda)
      !   CALL JFunc%get_J_value(X_Temp, f_lamda)

      !   IF (sg%isBaseProc()) PRINT *, 'rng_a', rng_a, 'rng_b', rng_b, 'f_lamda', f_lamda, 'f_miu', f_miu
      !   if (f_lamda > f_miu) then
      !     rng_a = lamda
      !   else
      !     rng_b = miu
      !   end if

      !   maxdRng = rng_b - rng_a
      !   do while (rng_b - rng_a > 0.00001*(maxdRng) .and. abs(f_lamda - f_miu) > 0.001)

      !     if (f_lamda > f_miu) then
      !       rng_a = lamda
      !       lamda = miu
      !       miu = rng_a + 0.618*(rng_b - rng_a)
      !       f_lamda = f_miu
      !       X_Temp = X + (d0*miu)
      !       CALL JFunc%get_J_value(X_Temp, f_miu)
      !     else
      !       rng_b = miu
      !       miu = lamda
      !       lamda = rng_a + 0.382*(rng_b - rng_a)
      !       f_miu = f_lamda
      !       X_Temp = X + (d0*lamda)

      !       CALL JFunc%get_J_value(X_Temp, f_lamda)
      !     end if

      !     IF (sg%isBaseProc()) PRINT *, 'rng_a', rng_a, 'rng_b', rng_b, 'f_lamda', f_lamda, 'f_miu', f_miu

      !   end do

      !   X_Temp = X + (d0*((rng_b + rng_a)/2.0))
      !   CALL JFunc%get_J_and_grad_J_vec(X_Temp, f_rng_b, g0)
      ! end if

      IF (isReCalJAndGradJ) THEN
        IF (sg%isBaseProc()) PRINT *, 'In re-calculation of J and Gradient J.'
        IF (this%HybInc) THEN
          CALL JFunc%get_J_and_grad_J_vec(X, f_rng_b, g0, dXs_Temp, dXe_Temp, g0e)
        ELSE
          CALL JFunc%get_J_and_grad_J_vec(X_Temp, f_rng_b, g0)
        END IF

!=========== compute the boundary =============
        ! Add constraints for qrain > 0, pres >0 and rho >0
        ! isTouchBound = .FALSE.
        DO i = 1, SIZE(X_Temp%Fields, 1)
          IF (TRIM(X_Temp%Fields(i)%Get_Name()) == 'rho_ctl' .OR. &
              TRIM(X_Temp%Fields(i)%Get_Name()) == 'rho' .OR. &
              TRIM(X_Temp%Fields(i)%Get_Name()) == 'rhov' .OR. &
              TRIM(X_Temp%Fields(i)%Get_Name()) == 'rhov_ctl' .OR. &
              TRIM(X_Temp%Fields(i)%Get_Name()) == 'qrain' .OR. &
              TRIM(X_Temp%Fields(i)%Get_Name()) == 'qvapor' .OR. &
              TRIM(X_Temp%Fields(i)%Get_Name()) == 'qvapor_ctl' .OR. &
              TRIM(X_Temp%Fields(i)%Get_Name()) == 'rhor' .OR. &
              TRIM(X_Temp%Fields(i)%Get_Name()) == 'pcpa' .OR. &
              TRIM(X_Temp%Fields(i)%Get_Name()) == 'pcpa5min' .OR. &
              TRIM(X_Temp%Fields(i)%Get_Name()) == 'ref' .OR. &
              TRIM(X_Temp%Fields(i)%Get_Name()) == 'temp' .OR. &
              TRIM(X_Temp%Fields(i)%Get_Name()) == 'rhor_ctl') THEN
            FORALL (tIdx=1:sg%tSlots, vIdx=1:sg%vLevel, hIdx=1:sg%num_cell, &
                    X_Temp%Fields(i)%DATA(vIdx, hIdx, tIdx) <= 0.0D0)
              d0%Fields(i)%DATA(vIdx, hIdx, tIdx) = -X%Fields(i)%DATA(vIdx, hIdx, tIdx) / rng_b / 2.0D0
              X_Temp%Fields(i)%DATA(vIdx, hIdx, tIdx) = X%Fields(i)%DATA(vIdx, hIdx, tIdx) / 2.0D0
              ! isTouchBound = .TRUE.
            END FORALL
          END IF
        END DO
        ! CALL sg%mpddInfo_sg%AllRedLOR(isTouchBound, isTouchBoundAll)
        ! isReCalJAndGradJ = isReCalJAndGradJ .OR. isTouchBoundAll
        ! IF(isTouchBound) X_Temp = X + (d0*rng_b)
!=========== compute the boundary =============
      END IF

      !=========== compute the boundary =============
      ! Add constraints for qrain > 0, pres >0 and rho >0
      ! isTouchBound = .FALSE.
      X_fullstate = X_Temp
      CALL JFunc%get_fullstate_CtlVar(X_fullstate)
      DO i = 1, SIZE(X_fullstate%Fields, 1)
        IF (TRIM(X_fullstate%Fields(i)%Get_Name()) == 'rho_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'rho' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'rhov' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'rhov_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qrain' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qvapor' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qvapor_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'rhor' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'pcpa' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'pcpa5min' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'ref' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'temp' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'rhor_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qcloud' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qcloud_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qice' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qice_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qrain' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qrain_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qsnow' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qsnow_ctl' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qgraupel' .OR. &
            TRIM(X_fullstate%Fields(i)%Get_Name()) == 'qgraupel_ctl') THEN
            FORALL (tIdx=1:sg%tSlots, vIdx=1:sg%vLevel, hIdx=1:sg%num_cell, &
            X_fullstate%Fields(i)%data(vIdx, hIdx, tIdx) <= 0.0D0)
              d0%Fields(i)%data(vIdx, hIdx, tIdx) = -X%Fields(i)%data(vIdx, hIdx, tIdx)/rng_b/2.0D0
              X_Temp%Fields(i)%data(vIdx, hIdx, tIdx) = X%Fields(i)%data(vIdx, hIdx, tIdx)/2.0D0
              ! isTouchBound = .TRUE.
            END FORALL
          
        END IF
      END DO
      ! CALL sg%mpddInfo_sg%AllRedLOR(isTouchBound, isTouchBoundAll)
      ! isReCalJAndGradJ = isReCalJAndGradJ .OR. isTouchBoundAll
      ! IF(isTouchBound) X_Temp = X + (d0*rng_b)
!=========== compute the boundary =============

      ! Add constraints for qrain > 0, pres >0 and rho >0
      ! BLOCK
      !   USE Field_m, ONLY: Field_t
      !   TYPE(Field_t) :: rngVec
      !   REAL(r_kind) :: rngTemp, rngMin
      !   rngMin = 1.0D100

      !   IF (X_Temp%getVarIdx('rho_ctl') .NE. 0) THEN
      !     IF (sg%mpddInfo_sg%AllRedMinReal(MINVAL(X_Temp%Fields(X_Temp%getVarIdx('rho_ctl'))%data)) < 0 .OR. isnan(rng_b)) THEN
      !       rngVec = X%Fields(X%getVarIdx('rho_ctl'))/d0%Fields(d0%getVarIdx('rho_ctl'))*(-1.0D0)
      !       rngTemp = sg%mpddInfo_sg%AllRedMinReal(MINVAL(rngVec%data, MASK=rngVec%data .GT. 0.0D0))
      !       PRINT *, 'rngTemp: ', rngTemp
      !       IF (rngTemp < rngMin) rngMin = rngTemp
      !     END IF
      !   END IF

      !   IF (X_Temp%getVarIdx('qrain_ctl') .NE. 0) THEN
      !     IF (sg%mpddInfo_sg%AllRedMinReal(MINVAL(X_Temp%Fields(X_Temp%getVarIdx('qrain_ctl'))%data)) < 0 .OR. isnan(rng_b)) THEN
      !       rngVec = X%Fields(X%getVarIdx('qrain_ctl'))/d0%Fields(d0%getVarIdx('qrain_ctl'))*(-1.0D0)
      !       rngTemp = sg%mpddInfo_sg%AllRedMinReal(MINVAL(rngVec%data, MASK=rngVec%data .GT. 0.0D0))
      !       PRINT *, 'rngTemp: ', rngTemp
      !       IF (rngTemp < rngMin) rngMin = rngTemp
      !     END IF
      !   END IF

      !   IF (X_Temp%getVarIdx('pres') .NE. 0) THEN
      !     IF (sg%mpddInfo_sg%AllRedMinReal(MINVAL(X_Temp%Fields(X_Temp%getVarIdx('pres'))%data)) < 0 .OR. isnan(rng_b)) THEN
      !       rngVec = X%Fields(X%getVarIdx('pres'))/d0%Fields(d0%getVarIdx('pres'))*(-1.0D0)
      !       rngTemp = sg%mpddInfo_sg%AllRedMinReal(MINVAL(rngVec%data, MASK=rngVec%data .GT. 0.0D0))
      !       PRINT *, 'rngTemp: ', rngTemp
      !       IF (rngTemp < rngMin) rngMin = rngTemp
      !     END IF
      !   END IF

      !   IF (X_Temp%getVarIdx('qvapor') .NE. 0) THEN
      !     IF (sg%mpddInfo_sg%AllRedMinReal(MINVAL(X_Temp%Fields(X_Temp%getVarIdx('qvapor'))%data)) < 0 .OR. isnan(rng_b)) THEN
      !       rngVec = X%Fields(X%getVarIdx('qvapor'))/d0%Fields(d0%getVarIdx('qvapor'))*(-1.0D0)
      !       rngTemp = sg%mpddInfo_sg%AllRedMinReal(MINVAL(rngVec%data, MASK=rngVec%data .GT. 0.0D0))
      !       PRINT *, 'rngTemp: ', rngTemp
      !       IF (rngTemp < rngMin) rngMin = rngTemp
      !     END IF
      !   END IF

      !   IF (X_Temp%getVarIdx('rhor_ctl') .NE. 0) THEN
      !     IF (sg%mpddInfo_sg%AllRedMinReal(MINVAL(X_Temp%Fields(X_Temp%getVarIdx('rhor_ctl'))%data)) < 0 .OR. isnan(rng_b)) THEN
      !       rngVec = X%Fields(X%getVarIdx('rhor_ctl'))/d0%Fields(d0%getVarIdx('rhor_ctl'))*(-1.0D0)
      !       rngTemp = sg%mpddInfo_sg%AllRedMinReal(MINVAL(rngVec%data, MASK=rngVec%data .GT. 0.0D0))
      !       PRINT *, 'rngTemp: ', rngTemp
      !       IF (rngTemp < rngMin) rngMin = rngTemp
      !     END IF
      !   END IF

      !   ! IF (X_Temp%getVarIdx('rhor') .NE. 0) THEN
      !   !   IF (sg%mpddInfo_sg%AllRedMinReal(MINVAL(X_Temp%Fields(X_Temp%getVarIdx('rhor'))%data)) < 0 .OR. isnan(rng_b)) THEN
      !   !     rngVec = X%Fields(X%getVarIdx('rhor'))/d0%Fields(d0%getVarIdx('rhor'))*(-1.0D0)
      !   !     rngTemp = sg%mpddInfo_sg%AllRedMinReal(MINVAL(rngVec%data, MASK=rngVec%data .GT. 0.0D0))
      !   !     PRINT *, 'rngTemp: ', rngTemp
      !   !     IF (rngTemp < rngMin) rngMin = rngTemp
      !   !   END IF
      !   ! END IF

      !   IF (rngMin < rng_b .OR. (isnan(rng_b) .AND. (rngMin < 1.0D3))) THEN
      !     X_Temp = X + d0*(rngMin*0.5)
      !     ! PRINT *, 'MINVAL(X_Temp%FIELD(1)%DATA)', MINVAL(X_Temp%FIELDs(1)%DATA)
      !     CALL JFunc%get_J_and_grad_J_vec(X_Temp, f_rng_b, g0)

      !     rng_b = rngMin*0.5
      !     IF (f_rng_b > f_min) THEN
      !       f_rng_b_last = f_rng_b
      !       DO WHILE (f_rng_b > f_min .OR. (f_rng_b_last - f_rng_b) > 1e-3)
      !         f_rng_b_last = f_rng_b
      !         rng_b = rng_b/2.0D0
      !         X_Temp = X + (d0*(rng_b))
      !         CALL JFunc%get_J_value(X_Temp, f_rng_b)
      !         IF (sg%isBaseProc()) PRINT *, 'New rng_b: ', rng_b, ' f_rng_b:', f_rng_b
      !       END DO

      !       rng_b = rng_b*2
      !       X_Temp = X + (d0*(rng_b))
      !       CALL JFunc%get_J_and_grad_J_vec(X_Temp, f_rng_b, g0)
      !       IF (sg%isBaseProc()) PRINT *, 'New rng_b: ', rng_b, ' f_rng_b:', f_rng_b
      !     END IF

      !     ! IF (sg%isBaseProc()) PRINT *, 'New rng_b from constraint: ', rngMin, ' rng_b:', rng_b
      !   END IF
      ! END BLOCK

      ! Final step length
      f_min = f_rng_b

      ! Jiongming Pang added HybInc (Hybrid in incremental analysis) referred to GSI
      IF (this%HybInc) THEN
        dXs = dXs_Temp
        dXe = dXe_Temp
      ELSE
        X = X_Temp
      END IF

#ifdef TRACK_DEBUG_INFO
      ! Yuanfu added a print information for checking Jo components:
      CALL JFunc%obsResidue(X_Temp, loopIdx, MaxOptStep, costsJo(1, 2))
#endif

      ! Update d0 to CG direction
      dotg1g1 = (g0.DOT.g0)
      IF (this%HybInc) THEN
        ! dotg1g1 = (g0.DOT.g0)
        ! IF ( JFunc%beta_e > 0.0D0 ) THEN
        DO i = 1, this%ensNum
          dotg1g1 = dotg1g1 + (g0e(i) .DOT.g0e(i))
        END DO
        ! END IF
        ! ELSE
        !   dotg1g1 = (g0.DOT.g0)
      END IF
      beta = dotg1g1 / dotg0g0

      d0 = (g0 * (-1.0D0)) + (d0 * beta)
      IF (this%HybInc) THEN
        DO i = 1, this%ensNum
          d0e(i) = (g0e(i) * (-1.0D0)) + (d0e(i) * beta)
        END DO
      END IF

      dotg0g0 = dotg1g1


      ! Print infos
      costsJo(loopIdx, 1) = f_min ! Use the first iteration function as normalizing factor
      IF (sg%isBaseProc()) &
        WRITE (*, 2) X%sg%gLevel, loopIdx, f_min, SQRT(dotg0g0), beta
2     FORMAT('Minimizing G', I1, ' step ', I3, ' f', E14.6, ' |g|', E14.6, ' beta', E14.6)

      loopIdx = loopIdx + 1

      ! CALL X%exHalo()
      ! CALL d0%exHalo()

      ! Output state file to NC for view.
      ! IF (loopIdx .eq. 100) THEN
      !   CALL Output_NC_State_SV(JFunc%Xb, '/Users/qzl/source/MOTOR-DA/grid/ncOutputFileBeforeDA1.nc', "ua")
      !   CALL Output_NC_State_SV(X, '/Users/qzl/source/MOTOR-DA/grid/ncOutputFileAfterDA1.nc', "ua")
      !   CALL Output_NC_State_SV(X-JFunc%Xb, '/Users/qzl/source/MOTOR-DA/grid/ncOutputFileDX1.nc', "ua")
      ! END IF

      ! Jiongming Pang added HybInc (Hybrid in incremental analysis) referred to GSI
      ! Updated the final analysis field X
      IF (this%HybInc) THEN
        X = JFunc%Xb_inc + dXs
        
        DO i = 1, this%ensNum
          Aw_tmp = JFunc%B_e.SQRTINVMUL. (dXe(i))
          CALL JFunc%D_opr(Aw_tmp, JFunc%Xb_enperts(i), DAw_tmp)
          CALL JFunc%H%Ctl2State%transBackward(DAw_tmp)
          X = X + DAw_tmp

        END DO
      END IF

    END DO

    ! Save the f value at the last.
    this%J = f_min
#ifdef DEBUG
    ! CALL JFunc%normalized(X,loopIdx-1,MaxOptStep,costsJo)
#endif
    IF (ALLOCATED(costsJo)) DEALLOCATE (costsJo)

  END SUBROUTINE run_HyyInc

END MODULE SolverFRCG_HybInc_m
