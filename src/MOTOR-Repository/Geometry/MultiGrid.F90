!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.MultiGrid.MultiGrid
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie, Zilong Qin
! VERSION           : V 0.1
! HISTORY           :
!   Created by Yuanfu Xie, 2019/04, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie 2020-4 for adding boundary conditions
!   Modified by Zilong Qin (zilong.qin@gmail.com), 2020/11/24, @GBA-MWF, Shenzhen, for parallelization
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/05/18, @GBA-MWF, Shenzhen, for extension of
!     the multigrid into vertical and temporal directions.
!!--------------------------------------------------------------------------------------------------

!> @brief
!! Modules for Data structures for multi grid structures.
!! @author Zilong Qin, Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
!! @note The structure is reforged from the original poisson solver.
! @warning
! @attention
MODULE multiGrid_m
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE singleGrid_m, ONLY: singleGrid_t
  USE mpddGlob_m, ONLY: mpddGlob_t

  !  USE mpi_f08
  INCLUDE "mpif.h"
  PUBLIC  ::  multiGrid_t

!> @brief
!! Data structures for multi grid structures.
! @see
! @note
! @warning
! @attention
  TYPE :: multiGrid_t
    INTEGER(i_kind) :: mg_coarsest, &              !<  Level of coarsest grid
                       mg_finest                   !<  Level of finest grid
    TYPE(mpddGlob_t), POINTER :: mpddGlob        !< Pointer of global multiprocessing domain decomposition class

    TYPE(singleGrid_t), ALLOCATABLE :: sg(:)       !<  Array of multiple single grid structs
    CHARACTER(len=1024) :: configFile
  CONTAINS

    FINAL :: destructor
    PROCEDURE, PUBLIC :: initialize

    PROCEDURE, PUBLIC :: restrictionAtGLevel
    PROCEDURE, PUBLIC :: prolongationAtGLevel

    PROCEDURE, PUBLIC, NOPASS :: restriction
    PROCEDURE, PUBLIC, NOPASS :: prolongation

    PROCEDURE, PUBLIC :: restrictionOfStatics

    PROCEDURE, PUBLIC, NOPASS :: prolongation3D, restriction3D

  END TYPE multiGrid_t

CONTAINS

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(multiGrid_t), INTENT(INOUT) :: this

    ! Deallocate memory
    IF (ALLOCATED(this%sg)) DEALLOCATE (this%sg)
  END SUBROUTINE destructor

  SUBROUTINE initialize(this, mpddGlob, group, rank, proc_layout, t_steps_mg, start_time, end_time, &
                        mg_coarsest, mg_finest, vLevel, configFile)
    IMPLICIT NONE
    CLASS(multiGrid_t) this
    TYPE(mpddGlob_t), TARGET, INTENT(IN) :: mpddGlob        !< Global multiprocessing domain decomposition class pointer
    INTEGER(i_kind), INTENT(IN) :: group, vLevel
    INTEGER(i_kind), INTENT(IN) :: mg_coarsest, mg_finest
    INTEGER(i_kind), INTENT(IN) :: rank, proc_layout(mg_coarsest:mg_finest, 2), t_steps_mg(mg_coarsest:mg_finest)
    REAL(r_kind), INTENT(IN) :: start_time, end_time
    CHARACTER(len=1024) :: configFile

    !===============================================================================
    ! Yuanfu Xie: 2022/05/18: for the multigrid in vertical and temporal directions:
    ! The strategy is to use vLevel and time_steps given in the namelist &geometry
    ! to construct the multigrid in vertical and temporal directions.
    !
    ! The current implementation is even distances in both directions!
    ! Local variables:
    INTEGER(i_kind) :: v_level_mg(mg_coarsest:mg_finest), ilvl

    v_level_mg(mg_finest) = vLevel
    DO ilvl = mg_finest - 1, mg_coarsest, -1
      ! Keep the coarsest grid in the coarser levels, i.e., only apply possible multigrid
      ! at the finest levels:
      v_level_mg(ilvl) = v_level_mg(ilvl + 1)
      IF (MOD(v_level_mg(ilvl) - 1, 2) .EQ. 0 .AND. (v_level_mg(ilvl) > 300)) &
        v_level_mg(ilvl) = (v_level_mg(ilvl) - 1) / 2 + 1
    END DO
    DO ilvl = mg_coarsest, mg_finest
      WRITE (*, 1) ilvl, v_level_mg(ilvl), t_steps_mg(ilvl), proc_layout(ilvl, :)
1     FORMAT('Multigrid levels at level: ', I3, ' v_levels: ', I3, ' t_step: ', I3, ' processor layout: ', 2I2)
    END DO

    this%mpddGlob => mpddGlob
    this%mg_finest = mg_finest
    this%mg_coarsest = mg_coarsest
    ALLOCATE (this%sg(this%mg_coarsest:this%mg_finest))

    DO ilvl = this%mg_coarsest, this%mg_finest
      ! this%sg(i) = singleGrid_t(mpddGlob, ilvl, vLevel, group, rank, &
      !proc_layout(ilvl, :), t_steps_mg(ilvl), start_time, end_time, configFile)
      ! Yuanfu Xie: 2022/05/18 changed to multigrid in vertical and temporal directions:
      CALL this%sg(ilvl)%initialize(mpddGlob, ilvl, v_level_mg(ilvl), group, rank, &
                                    proc_layout(ilvl, :), t_steps_mg(ilvl), start_time, end_time, configFile)

      CALL mpddGlob%barrier
    END DO

    this%configFile = configFile

  END SUBROUTINE

!> @brief
!! Prolongate the grid values from this level to a finer level
! @see
! @note
! @warning
! @attention
  SUBROUTINE prolongationAtGLevel(this, gLevel, vLevel, bufFiner, bufCoarser, isUpdateBndy)
    IMPLICIT NONE
    CLASS(multiGrid_t) this
    INTEGER(i_kind), INTENT(IN) :: gLevel, &                   !< Current multigrid level of grid
                                   vLevel                      !< Vertical level of the grid
    REAL(r_kind), INTENT(OUT) :: bufFiner(vLevel, this%sg(gLevel + 1)%num_cell)   !< Buffer address at finer level
    REAL(r_kind), INTENT(IN) :: bufCoarser(vLevel, this%sg(gLevel)%num_cell)      !< Buffer address at coarser level
    LOGICAL, INTENT(IN), OPTIONAL :: isUpdateBndy              !< Whether prolongate the cells on the boundary

    ! Yuanfu Xie: 2022/05/18 added a vertical increment for a multigrid:
    INTEGER(i_kind) :: vInc

    IF (gLevel .EQ. this%mg_finest) THEN
      PRINT *, 'Data at finest level cannot be interpolated.'
      STOP
    END IF

    ! Yuanfu Xie: 2022/05/18 changed this for multigrid in vertical:
    IF (this%sg(gLevel)%vLevel .EQ. this%sg(gLevel + 1)%vLevel) THEN
      vInc = 1
    ELSE IF ((this%sg(gLevel)%vLevel - 1) / 2 .EQ. this%sg(gLevel + 1)%vLevel) THEN
      vInc = 2
    ELSE
      WRITE (*, 1) this%sg(gLevel)%vLevel, this%sg(gLevel - 1)%vLevel
1     FORMAT('restrictionAtGLevel -- error: vertical levels are not suitable for a multigrid: ', 2I3)
    END IF

    CALL prolongation(this%sg(gLevel + 1), this%sg(gLevel), bufFiner, bufCoarser, vLevel, vInc, isUpdateBndy)
  END SUBROUTINE prolongationAtGLevel

  !> @brief
!! Restrict the grid values from this level to a finer level
! @see
! @note
! @warning
! @attention
  SUBROUTINE restrictionAtGLevel(this, gLevel, vLevel, bufFiner, bufCoarser, isUpdateBndy)
    IMPLICIT NONE
    CLASS(multiGrid_t) this
    INTEGER(i_kind), INTENT(IN) :: gLevel, &                   !< Current multigrid level of grid
                                   vLevel                      !< Vertical level of the grid
    REAL(r_kind), INTENT(IN) :: bufFiner(vLevel, this%sg(gLevel)%num_cell)          !< Buffer address at finer level
    REAL(r_kind), INTENT(OUT) :: bufCoarser(vLevel, this%sg(gLevel - 1)%num_cell)   !< Buffer address at coarser level
    LOGICAL, INTENT(IN), OPTIONAL :: isUpdateBndy              !< Whether prolongate the cells on the boundary

    ! Yuanfu Xie: 2022/05/18 added a vertical increment for a multigrid:
    INTEGER(i_kind) :: vInc

    IF (gLevel .EQ. this%mg_coarsest) THEN
      PRINT *, 'Data at coarsest level cannot be restricted.'
      STOP
    END IF

    ! Yuanfu Xie: 2022/05/18 changed this for multigrid in vertical:
    IF (this%sg(gLevel)%vLevel .EQ. this%sg(gLevel - 1)%vLevel) THEN
      vInc = 1
    ELSE IF ((this%sg(gLevel)%vLevel - 1) / 2 .EQ. this%sg(gLevel - 1)%vLevel) THEN
      vInc = 2
    ELSE
      WRITE (*, 1) this%sg(gLevel)%vLevel, this%sg(gLevel - 1)%vLevel
1     FORMAT('restrictionAtGLevel -- error: vertical levels are not suitable for a multigrid: ', 2I3)
    END IF
    CALL restriction(this%sg(gLevel), this%sg(gLevel - 1), bufFiner, bufCoarser, vLevel, vInc, isUpdateBndy)
  END SUBROUTINE restrictionAtGLevel

  SUBROUTINE restriction(sgFiner, sgCoarser, bufFiner, bufCoarser, vLevel, vInc, isUpdateBndy)
    IMPLICIT NONE
    TYPE(singleGrid_t), INTENT(IN) :: sgFiner, sgCoarser
    REAL(r_kind), INTENT(IN) :: bufFiner(:, :)
    REAL(r_kind), INTENT(OUT) :: bufCoarser(:, :)
    INTEGER(i_kind), INTENT(IN) :: vLevel, vInc ! Yuanfu Xie 2022/05/18 added vInc for vertical restriction
    REAL(r_kind), ALLOCATABLE :: bufTemp(:, :)
    LOGICAL, INTENT(IN), OPTIONAL :: isUpdateBndy

    INTEGER(i_kind) :: i, j, k, norm
    LOGICAL :: isUpdBndy

    IF (.NOT. PRESENT(isUpdateBndy)) THEN
      isUpdBndy = .TRUE.
    ELSE
      isUpdBndy = isUpdateBndy
    END IF

    IF (.NOT. sgFiner%mpddInfo_sg%isActiveProc()) RETURN
    ALLOCATE (bufTemp(vLevel, sgCoarser%num_icell_toFiner))

    ! print *,'Myrank is ',sgCoarser%mpddInfo_sg%myrank,'Index is', sgCoarser%c_tc_idx
    DO i = LBOUND(sgCoarser%c_tc_idx, 2), UBOUND(sgCoarser%c_tc_idx, 2)
      bufTemp(:, i) = 0.0D0
      norm = 0.0D0
      DO j = LBOUND(sgCoarser%c_tc_idx, 1), UBOUND(sgCoarser%c_tc_idx, 1)
        IF (sgCoarser%c_tc_idx(j, i) .EQ. 0) CYCLE
        bufTemp(:, i) = bufTemp(:, i) + bufFiner(::vInc, sgCoarser%c_tc_idx(j, i))
        norm = norm + 1.0D0
      END DO
      bufTemp(:, i) = bufTemp(:, i) / norm
    END DO

    IF (sgFiner%mpddInfo_sg%nProc .NE. sgCoarser%mpddInfo_sg%nProc) THEN
      CALL sgFiner%mpddInfo_sg%gather_to_coarser_threads(bufCoarser, &
                                                         bufTemp, &
                                                         sgFiner%mpddInfo_sg%proc_layout, &
                                                         sgCoarser%mpddInfo_sg%proc_layout, &
                                                         vlevel, &
                                                         sgCoarser%num_icell_toFiner, &
                                                         sgCoarser%sp_t_g_idx_toFiner, &
                                                         sgCoarser%g_t_sp_idx)
    ELSE
      ! PRINT *, 'layout is equal.', shape(bufTemp), shape(bufCoarser)
      bufCoarser(:, 1:sgCoarser%num_icell) = bufTemp
    END IF

    ! Set the boundary cell to 0.0 if not update boundaries.
    IF (.NOT. isUpdBndy) THEN
      DO i = 1, sgCoarser%num_icell
        IF ((sgCoarser%cell_type(i) .EQ. 1) .OR. (sgCoarser%cell_type(i) .EQ. 2)) THEN              ! On the inner cells, rhs is the vorticity
          bufCoarser(:, i) = 0.0D0
        END IF
      END DO
    END IF

    CALL sgCoarser%mpddInfo_sg%ExchangeMatOnHalo(sgCoarser%num_icell, &
                                                 vlevel, bufCoarser, sgCoarser%g_t_sp_idx)
    DEALLOCATE (bufTemp)

  END SUBROUTINE restriction

  SUBROUTINE prolongation(sgFiner, sgCoarser, bufFiner, bufCoarser, vLevel, vInc, isUpdateBndy)
    IMPLICIT NONE
    TYPE(singleGrid_t), INTENT(IN) :: sgFiner, sgCoarser
    REAL(r_kind), INTENT(OUT) :: bufFiner(:, :)
    REAL(r_kind), INTENT(IN) :: bufCoarser(:, :)
    INTEGER(i_kind), INTENT(IN) :: vLevel, vInc   ! Yuanfu Xie: 2022/05/18 added vInc for a multigrid in vertical
    LOGICAL, INTENT(IN), OPTIONAL :: isUpdateBndy

    REAL(r_kind), ALLOCATABLE :: bufTemp(:, :)
    INTEGER(i_kind) :: i, j, k

    ! Yuanfu Xie: 2022/05/18 added a fvLevel for a multigrid in the vertical:
    INTEGER(i_kind) :: fvLevel
    REAL(r_kind) :: c1, c2

    LOGICAL :: isUpdBndy

    ! Fine grid vertical levels:
    fvLevel = (vLevel - 1) * vInc + 1

    IF (.NOT. PRESENT(isUpdateBndy)) THEN
      isUpdBndy = .TRUE.
    ELSE
      isUpdBndy = isUpdateBndy
    END IF

    IF (.NOT. sgFiner%mpddInfo_sg%isActiveProc()) RETURN

    IF (sgCoarser%mpddInfo_sg%isActiveProc()) THEN
      ALLOCATE (bufTemp(fvLevel, sgFiner%num_icell_toCoarser))
      bufTemp = 0
      DO i = LBOUND(sgFiner%c2f_intp_stcl, 2), UBOUND(sgFiner%c2f_intp_stcl, 2)
        DO j = LBOUND(sgFiner%c2f_intp_stcl, 1), UBOUND(sgFiner%c2f_intp_stcl, 1)
          IF (sgFiner%c2f_intp_stcl(j, i) .EQ. 0) CYCLE
          IF (vInc .EQ. 1) THEN
            bufTemp(:, i) = bufTemp(:, i) + sgFiner%c2f_intp_cf(j, i) * bufCoarser(:, sgFiner%c2f_intp_stcl(j, i))
          ELSE
            DO k = 1, fvLevel, 2
              bufTemp(k, i) = bufTemp(k, i) + sgFiner%c2f_intp_cf(j, i) * &
                              bufCoarser((k - 1) / 2 + 1, sgFiner%c2f_intp_stcl(j, i))
            END DO
            DO k = 2, fvLevel, 2
              ! Change the interpolation coeifficients by Zilong
              c1 = (sgFiner%sigma(k + 1) - sgFiner%sigma(k)) / (sgFiner%sigma(k + 1) - sgFiner%sigma(k - 1))
              c2 = (sgFiner%sigma(k) - sgFiner%sigma(k - 1)) / (sgFiner%sigma(k + 1) - sgFiner%sigma(k - 1))
              bufTemp(k, i) = bufTemp(k, i) + sgFiner%c2f_intp_cf(j, i) * &
                              (bufCoarser(k / 2, sgFiner%c2f_intp_stcl(j, i)) * c1 + &
                               bufCoarser(k / 2 + 1, sgFiner%c2f_intp_stcl(j, i)) * c2)
            END DO
          END IF
        END DO
      END DO
      ! print *, sgFiner%c2f_intp_stcl
    END IF

    IF (sgFiner%mpddInfo_sg%nProc .NE. sgCoarser%mpddInfo_sg%nProc) THEN
      ! IF (sgCoarser%mpddInfo_sg%isActiveProc()) PRINT *, 'layout is not equal.', shape(bufTemp), shape(bufFiner)
      CALL sgFiner%mpddInfo_sg%scatter_to_finer_threads(bufFiner, &
                                                        bufTemp, &
                                                        sgFiner%mpddInfo_sg%proc_layout, &
                                                        sgCoarser%mpddInfo_sg%proc_layout, &
                                                        fvlevel, &
                                                        sgFiner%num_icell, &
                                                        sgFiner%g_t_sp_idx_toCoarser, &
                                                        sgFiner%sp_t_g_idx)
    ELSE
      ! PRINT *, 'layout is equal.', shape(bufTemp), shape(bufFiner)
      bufFiner(:, 1:sgFiner%num_icell) = bufTemp
    END IF

    ! Set the boundary cell to 0.0 if not update boundaries.
    IF (.NOT. isUpdBndy) THEN
      DO i = 1, sgFiner%num_icell
        IF ((sgFiner%cell_type(i) .EQ. 1) .OR. (sgFiner%cell_type(i) .EQ. 2)) THEN              ! On the inner cells, rhs is the vorticity
          bufFiner(:, i) = 0.0D0
        END IF
      END DO
    END IF

    CALL sgFiner%mpddInfo_sg%ExchangeMatOnHalo(sgFiner%num_icell, &
                                               fvlevel, bufFiner, sgFiner%g_t_sp_idx)

    IF (sgCoarser%mpddInfo_sg%isActiveProc()) DEALLOCATE (bufTemp)
  END SUBROUTINE prolongation

  SUBROUTINE restriction3D(sgFiner, sgCoarser, bufFiner, bufCoarser, vLevel, isUpdateBndy)
    IMPLICIT NONE
    TYPE(singleGrid_t), INTENT(IN) :: sgFiner, sgCoarser
    REAL(r_kind), INTENT(IN) :: bufFiner(:, :, :)
    REAL(r_kind), INTENT(OUT) :: bufCoarser(:, :, :)
    INTEGER(i_kind), INTENT(IN) :: vLevel
    REAL(r_kind), ALLOCATABLE :: bufTemp(:, :, :)
    LOGICAL, INTENT(IN), OPTIONAL :: isUpdateBndy

    INTEGER(i_kind) :: i, j, k, norm
    LOGICAL :: isUpdBndy

    IF (.NOT. PRESENT(isUpdateBndy)) THEN
      isUpdBndy = .TRUE.
    ELSE
      isUpdBndy = isUpdateBndy
    END IF

    IF (.NOT. sgFiner%mpddInfo_sg%isActiveProc()) RETURN
    ALLOCATE (bufTemp(vLevel, sgCoarser%num_icell_toFiner, sgCoarser%tSlots))

    ! print *,'Myrank is ',sgCoarser%mpddInfo_sg%myrank,'Index is', sgCoarser%c_tc_idx
    DO i = LBOUND(sgCoarser%c_tc_idx, 2), UBOUND(sgCoarser%c_tc_idx, 2)
      bufTemp(:, i, :) = 0.0D0
      norm = 0.0D0
      DO j = LBOUND(sgCoarser%c_tc_idx, 1), UBOUND(sgCoarser%c_tc_idx, 1)
        IF (sgCoarser%c_tc_idx(j, i) .EQ. 0) CYCLE
        bufTemp(:, i, :) = bufTemp(:, i, :) + bufFiner(:, sgCoarser%c_tc_idx(j, i), :)
        norm = norm + 1.0D0
      END DO
      bufTemp(:, i, :) = bufTemp(:, i, :) / norm
    END DO

    IF (sgFiner%mpddInfo_sg%nProc .NE. sgCoarser%mpddInfo_sg%nProc) THEN
      DO i = 1, sgFiner%tSlots
        CALL sgFiner%mpddInfo_sg%gather_to_coarser_threads(bufCoarser(:, :, i), &
                                                           bufTemp(:, :, i), &
                                                           sgFiner%mpddInfo_sg%proc_layout, &
                                                           sgCoarser%mpddInfo_sg%proc_layout, &
                                                           vlevel, &
                                                           sgCoarser%num_icell_toFiner, &
                                                           sgCoarser%sp_t_g_idx_toFiner, &
                                                           sgCoarser%g_t_sp_idx)
      END DO
    ELSE
      ! PRINT *, 'layout is equal.', shape(bufTemp), shape(bufCoarser)
      bufCoarser(:, 1:sgCoarser%num_icell, :) = bufTemp
    END IF

    ! Set the boundary cell to 0.0 if not update boundaries.
    IF (.NOT. isUpdBndy) THEN
      DO i = 1, sgCoarser%num_icell
        IF ((sgCoarser%cell_type(i) .EQ. 1) .OR. (sgCoarser%cell_type(i) .EQ. 2)) THEN              ! On the inner cells, rhs is the vorticity
          bufCoarser(:, i, :) = 0.0D0
        END IF
      END DO
    END IF

    CALL sgCoarser%ExchangeMatOnHaloForFieldGrid(sgCoarser%tSlots, vlevel, bufCoarser)

    DEALLOCATE (bufTemp)

  END SUBROUTINE restriction3D

  SUBROUTINE prolongation3D(sgFiner, sgCoarser, bufFiner, bufCoarser, vLevel, isUpdateBndy)
    IMPLICIT NONE
    TYPE(singleGrid_t), INTENT(IN) :: sgFiner, sgCoarser
    REAL(r_kind), INTENT(OUT) :: bufFiner(:, :, :)
    REAL(r_kind), INTENT(IN) :: bufCoarser(:, :, :)
    INTEGER(i_kind), INTENT(IN) :: vLevel
    LOGICAL, INTENT(IN), OPTIONAL :: isUpdateBndy

    REAL(r_kind), ALLOCATABLE :: bufTemp(:, :, :)
    INTEGER(i_kind) :: i, j, k

    LOGICAL :: isUpdBndy

    IF (.NOT. PRESENT(isUpdateBndy)) THEN
      isUpdBndy = .TRUE.
    ELSE
      isUpdBndy = isUpdateBndy
    END IF

    IF (.NOT. sgFiner%mpddInfo_sg%isActiveProc()) RETURN

    IF (sgCoarser%mpddInfo_sg%isActiveProc()) THEN
      ALLOCATE (bufTemp(vLevel, sgFiner%num_icell_toCoarser, sgFiner%tSlots))
      bufTemp = 0
      DO i = LBOUND(sgFiner%c2f_intp_stcl, 2), UBOUND(sgFiner%c2f_intp_stcl, 2)
        DO j = LBOUND(sgFiner%c2f_intp_stcl, 1), UBOUND(sgFiner%c2f_intp_stcl, 1)
          IF (sgFiner%c2f_intp_stcl(j, i) .EQ. 0) CYCLE
          bufTemp(:, i, :) = bufTemp(:, i, :) + sgFiner%c2f_intp_cf(j, i) * bufCoarser(:, sgFiner%c2f_intp_stcl(j, i), :)
        END DO
      END DO
      ! print *, sgFiner%c2f_intp_stcl
    END IF

    IF (sgFiner%mpddInfo_sg%nProc .NE. sgCoarser%mpddInfo_sg%nProc) THEN
      ! IF (sgCoarser%mpddInfo_sg%isActiveProc()) PRINT *, 'layout is not equal.', shape(bufTemp), shape(bufFiner)
      DO i = 1, sgFiner%tSlots
        CALL sgFiner%mpddInfo_sg%scatter_to_finer_threads(bufFiner(:, :, i), &
                                                          bufTemp(:, :, i), &
                                                          sgFiner%mpddInfo_sg%proc_layout, &
                                                          sgCoarser%mpddInfo_sg%proc_layout, &
                                                          vlevel, &
                                                          sgFiner%num_icell, &
                                                          sgFiner%g_t_sp_idx_toCoarser, &
                                                          sgFiner%sp_t_g_idx)
      END DO
    ELSE
      ! PRINT *, 'layout is equal.', shape(bufTemp), shape(bufFiner)
      bufFiner(:, 1:sgFiner%num_icell, :) = bufTemp
    END IF

    ! Set the boundary cell to 0.0 if not update boundaries.
    IF (.NOT. isUpdBndy) THEN
      DO i = 1, sgFiner%num_icell
        IF ((sgFiner%cell_type(i) .EQ. 1) .OR. (sgFiner%cell_type(i) .EQ. 2)) THEN              ! On the inner cells, rhs is the vorticity
          bufFiner(:, i, :) = 0.0D0
        END IF
      END DO
    END IF

    CALL sgFiner%ExchangeMatOnHaloForFieldGrid(sgFiner%tSlots, vlevel, bufFiner)

    IF (sgCoarser%mpddInfo_sg%isActiveProc()) DEALLOCATE (bufTemp)
  END SUBROUTINE prolongation3D

  SUBROUTINE restrictionOfStatics(this, sgFiner, sgCoarser)
    CLASS(multiGrid_t) this
    TYPE(singleGrid_t), INTENT(INOUT) :: sgFiner, sgCoarser
    REAL(r_kind), ALLOCATABLE :: bufCoarser(:, :)
    INTEGER(i_kind) :: vInc, i, j, k

    ! IF (.not. sgFiner%isActiveProc() .OR. .not. sgCoarser%isActiveProc()) RETURN ! This is incorrect and will cause hang.
    IF (.NOT. sgFiner%isActiveProc()) RETURN

    IF (sgFiner%vLevel .EQ. sgCoarser%vLevel) THEN
      vInc = 1
    ELSE IF ((sgFiner%vLevel - 1) / 2 + 1 .EQ. sgCoarser%vLevel) THEN
      vInc = 2
    END IF

    ALLOCATE (bufCoarser(1, sgCoarser%num_cell))

    DO j = 1, sgCoarser%tSlots
      DO k = 1, sgFiner%tSlots
        ! This following print statement is a bad example for people hard to track it down!!!
        IF (ABS(sgFiner%tt(k) - sgCoarser%tt(j)) .LT. 1E-6) THEN
          CALL this%restriction(sgFiner, sgCoarser, RESHAPE(sgFiner%tskin(:, k), (/1, sgFiner%num_cell/)), &
                                bufCoarser, 1, 1, .TRUE.)
          sgCoarser%tskin(:, j) = bufCoarser(1, :)
        END IF
      END DO
    END DO

    DO j = 1, sgCoarser%tSlots
      DO k = 1, sgFiner%tSlots
        ! This following print statement is a bad example for people hard to track it down!!!
        IF (ABS(sgFiner%tt(k) - sgCoarser%tt(j)) .LT. 1E-6) THEN
          CALL this%restriction(sgFiner, sgCoarser, RESHAPE(sgFiner%psfc(:, k), (/1, sgFiner%num_cell/)), &
                                bufCoarser, 1, 1, .TRUE.)
          sgCoarser%psfc(:, j) = bufCoarser(1, :)
        END IF
      END DO
    END DO

    DO j = 1, sgCoarser%tSlots
      DO k = 1, sgFiner%tSlots
        ! This following print statement is a bad example for people hard to track it down!!!
        IF (ABS(sgFiner%tt(k) - sgCoarser%tt(j)) .LT. 1E-6) THEN
          CALL this%restriction(sgFiner, sgCoarser, RESHAPE(sgFiner%u10m(:, k), (/1, sgFiner%num_cell/)), &
                                bufCoarser, 1, 1, .TRUE.)
          sgCoarser%u10m(:, j) = bufCoarser(1, :)
        END IF
      END DO
    END DO

    DO j = 1, sgCoarser%tSlots
      DO k = 1, sgFiner%tSlots
        ! This following print statement is a bad example for people hard to track it down!!!
        IF (ABS(sgFiner%tt(k) - sgCoarser%tt(j)) .LT. 1E-6) THEN
          CALL this%restriction(sgFiner, sgCoarser, RESHAPE(sgFiner%v10m(:, k), (/1, sgFiner%num_cell/)), &
                                bufCoarser, 1, 1, .TRUE.)
          sgCoarser%v10m(:, j) = bufCoarser(1, :)
        END IF
      END DO
    END DO

    DO j = 1, sgCoarser%tSlots
      DO k = 1, sgFiner%tSlots
        ! This following print statement is a bad example for people hard to track it down!!!
        IF (ABS(sgFiner%tt(k) - sgCoarser%tt(j)) .LT. 1E-6) THEN
          CALL this%restriction(sgFiner, sgCoarser, sgFiner%FGPres(:, :, k), &
                                sgCoarser%FGPres(:, :, j), sgCoarser%vLevel, vInc, .TRUE.)
        END IF
      END DO
    END DO

    CALL this%restriction(sgFiner, sgCoarser, RESHAPE(sgFiner%topo, (/1, sgFiner%num_cell/)), &
                          bufCoarser, 1, 1, .TRUE.)
    sgCoarser%topo = bufCoarser(1, :)
    CALL this%restriction(sgFiner, sgCoarser, RESHAPE(sgFiner%landmask, (/1, sgFiner%num_cell/)), &
                          bufCoarser, 1, 1, .TRUE.)
    sgCoarser%landmask = bufCoarser(1, :)

    ! Initialize horizonSimilarity: Yuanfu Xie 2023/10/18:
    CALL this%restriction(sgFiner, sgCoarser, RESHAPE(sgFiner%horizonSimilarity, (/1, sgFiner%num_cell/)), &
                          bufCoarser, 1, 1, .TRUE.)
    sgCoarser%horizonSimilarity = bufCoarser(1, :)

    DEALLOCATE (bufCoarser)

    CALL this%restriction(sgFiner, sgCoarser, sgFiner%pres, &
                          sgCoarser%pres, sgCoarser%vLevel, vInc, .TRUE.)
    CALL this%restriction(sgFiner, sgCoarser, sgFiner%s1, &
                          sgCoarser%s1, sgCoarser%vLevel, vInc, .TRUE.)
    CALL this%restriction(sgFiner, sgCoarser, sgFiner%SVapor, &
                          sgCoarser%SVapor, sgCoarser%vLevel, vInc, .TRUE.)
    CALL this%restriction(sgFiner, sgCoarser, sgFiner%STemp, &
                          sgCoarser%STemp, sgCoarser%vLevel, vInc, .TRUE.)
    IF (sgCoarser%vLevel == sgFiner%vLevel) THEN
      sgCoarser%sigma = sgFiner%sigma
      sgCoarser%SVapor1D = sgFiner%SVapor1D
      sgCoarser%STemp1D = sgFiner%STemp1D
      sgCoarser%s_qvapor = sgFiner%s_qvapor
    ELSE
      sgCoarser%sigma = sgFiner%sigma(::2)
      sgCoarser%SVapor1D = sgFiner%SVapor1D(::2)
      sgCoarser%STemp1D = sgFiner%STemp1D(::2)
      sgCoarser%s_qvapor = sgFiner%s_qvapor(::2)
    END IF
    CALL sgCoarser%update_zHght_from_topo_and_sigma

  END SUBROUTINE

END MODULE multiGrid_m
