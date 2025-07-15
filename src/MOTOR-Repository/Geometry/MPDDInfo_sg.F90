!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.mpdd.mpddInfo_sg
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Ning Wang, Jinrong Jiang
! VERSION           : V 0.0
! HISTORY           :
! The code is originally written by Ning Wang and Jinrong Jiang in mpi_meta_init.F90
!   Reforged by Zilong Qin (zilong.qin@gmail.com), 2020/11/25, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/09/07, @GBA-MWF, Shenzhen
!     Modified RemapAllGridIdxToSP2D so that it allows negative indices for adnb array's
!     grid distribution over multiple processors.
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module contains the data type for model_states method.
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
MODULE mpddInfo_sg_m
  USE kinds_m, ONLY: i_kind, r_kind, r_Double, i_llong
  USE mpddSub_m, ONLY: mpddSub_t
  USE omp_lib
  IMPLICIT NONE

  INCLUDE "mpif.h"

  TYPE NBlock
    INTEGER(i_kind), ALLOCATABLE :: idx(:)
    INTEGER(i_kind) :: sz
    INTEGER(i_kind) :: blk_num
  END TYPE NBLOCK

  TYPE HaloBlock
    TYPE(NBlock), ALLOCATABLE :: nb_blk(:)
    INTEGER(i_kind) :: nb_blk_sz ! Zilong added
    ! INTEGER(i_kind) :: my_start_idx, my_END_idx
    INTEGER(i_kind) :: num_icell
    INTEGER(i_kind) :: halo_sz
  END TYPE HaloBlock

  TYPE, EXTENDS(mpddSub_t) :: mpddInfo_sg_t
    ! INTEGER(i_kind) ::  start, end
    INTEGER(i_kind) :: halo_size
    TYPE(HaloBlock) :: recv_info, send_info
    INTEGER(i_kind), ALLOCATABLE :: remap(:)
    REAL(r_kind), ALLOCATABLE :: sendbuf(:, :), recvbuf(:, :)
    REAL(r_kind), ALLOCATABLE :: swapSendbuf(:, :), swapRecvbuf(:, :)
    INTEGER(i_kind) :: proc_layout(2)
    INTEGER(i_kind) :: gLevel

    ! INTEGER(i_kind) :: num_icell_global
    ! INTEGER(i_kind) :: num_icell
  CONTAINS

    FINAL  :: destructor

    PROCEDURE, PUBLIC :: initializemMPDDSub
    PROCEDURE, PUBLIC, PASS   :: buf_alloc
    PROCEDURE, PUBLIC, PASS   :: buf_dealloc

    PROCEDURE, PUBLIC, PASS   :: DistGridIntSeq1D
    PROCEDURE, PUBLIC, PASS   :: DistGridIntSeq2D
    PROCEDURE, PUBLIC, PASS   :: DistGridIntSeq3D
    PROCEDURE, PUBLIC, PASS   :: DistGridIntSeq4D

    GENERIC, PUBLIC :: DistGridIntSeq => DistGridIntSeq1D, DistGridIntSeq2D, DistGridIntSeq3D, DistGridIntSeq4D

    PROCEDURE, PUBLIC, PASS   :: DistGridRealSeq1D
    PROCEDURE, PUBLIC, PASS   :: DistGridRealSeq2D
    PROCEDURE, PUBLIC, PASS   :: DistGridRealSeq3D
    PROCEDURE, PUBLIC, PASS   :: DistGridRealSeq4D

    GENERIC, PUBLIC :: DistGridRealSeq => DistGridRealSeq1D, DistGridRealSeq2D, DistGridRealSeq3D, DistGridRealSeq4D

    PROCEDURE, PUBLIC, PASS   :: aggrGridReal2D, aggrGridReal1D
    PROCEDURE, PUBLIC, PASS   :: aggrGridReal3D

    PROCEDURE, PUBLIC, PASS   :: genMapIdx
    PROCEDURE, PRIVATE, NOPASS   :: calRowColRangeSP

    PROCEDURE, PUBLIC, PASS   :: generate_halo_info
    PROCEDURE, PUBLIC, PASS   :: generate_recv_info
    PROCEDURE, PUBLIC, PASS   :: generate_send_info
    PROCEDURE, PUBLIC, PASS   :: cal_num_cell_stcl

    PROCEDURE, PRIVATE, PASS  :: i_global2ProcRank
    PROCEDURE, PUBLIC, PASS   :: trim_idx_all_halos
    PROCEDURE, PRIVATE, PASS  :: halo_stcl_idx_remap

    PROCEDURE, PUBLIC, PASS   :: RemapInnerGridIdxToSP3D
    PROCEDURE, PUBLIC, PASS   :: RemapAllGridIdxToSP2D
    PROCEDURE, PUBLIC, PASS   :: RemapAllGridIdxToSP3D

    PROCEDURE, PUBLIC, PASS   :: genCellType
    PROCEDURE, PUBLIC, PASS   :: ExchangeMatOnHalo
    PROCEDURE, PUBLIC, PASS   :: ExchangeMatOnHaloReverseSum

    PROCEDURE, PUBLIC, PASS   :: calCellSize
    PROCEDURE, PUBLIC, PASS   :: gather_to_coarser_threads
    PROCEDURE, PUBLIC, PASS   :: scatter_to_finer_threads

  END TYPE mpddInfo_sg_t

  INTEGER, PARAMETER :: max_nb = 10

CONTAINS

  SUBROUTINE aggrGridReal1D(this, nCell_perProc, buf_src, buf_dist, buf_dist_shape, seq_store_idx)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) this
    INTEGER(i_kind), INTENT(IN) :: buf_dist_shape, nCell_perProc, seq_store_idx(:)
    REAL(r_kind), INTENT(IN)    :: buf_src(:)
    REAL(r_kind), INTENT(INOUT) :: buf_dist(:)
    INTEGER, ALLOCATABLE :: gcount_group(:), disp_group(:)
    INTEGER(i_kind) :: k1
    REAL(r_kind), ALLOCATABLE :: bufTemp(:)

    ! Return if the proc is not active at this grid level
    IF (.NOT. this%isActiveProc()) RETURN

    IF (this%isBaseProc()) THEN
      ALLOCATE (bufTemp(buf_dist_shape))
    END IF

    ALLOCATE (gcount_group(this%nProc), disp_group(this%nProc))
    CALL MPI_ALLGATHER(nCell_perProc, 1, MPI_INTEGER4, gcount_group, 1, MPI_INTEGER4, this%comm, this%IERR)

    gcount_group = gcount_group

    disp_group = 0
    DO k1 = 2, this%nProc
      disp_group(k1) = gcount_group(k1 - 1) + disp_group(k1 - 1)
    END DO

    CALL this%barrier
    CALL mpi_gatherv(buf_src, gcount_group(this%myrank), MPI_DOUBLE_PRECISION, &
                     bufTemp, gcount_group, disp_group, MPI_DOUBLE_PRECISION, this%rankBase, this%comm, this%ierr)

    IF (this%isBaseProc()) THEN
      buf_dist(seq_store_idx) = bufTemp
    END IF

    DEALLOCATE (gcount_group, disp_group)
    IF (this%isBaseProc()) DEALLOCATE (bufTemp)
  END SUBROUTINE aggrGridReal1D

  SUBROUTINE aggrGridReal2D(this, nCell_perProc, buf_src, buf_dist, buf_dist_shape, seq_store_idx)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) this
    INTEGER(i_kind), INTENT(IN) :: buf_dist_shape(2), nCell_perProc, seq_store_idx(:)
    REAL(r_kind), INTENT(IN)    :: buf_src(:, :)
    REAL(r_kind), INTENT(INOUT) :: buf_dist(:, :)
    INTEGER, ALLOCATABLE :: gcount_group(:), disp_group(:)
    INTEGER(i_kind) :: k1
    REAL(r_kind), ALLOCATABLE :: bufTemp(:, :)

    ! Return if the proc is not active at this grid level
    IF (.NOT. this%isActiveProc()) RETURN

    IF (this%isBaseProc()) THEN
      ALLOCATE (bufTemp(buf_dist_shape(1), buf_dist_shape(2)))
    END IF

    ALLOCATE (gcount_group(this%nProc), disp_group(this%nProc))
    CALL MPI_ALLGATHER(nCell_perProc, 1, MPI_INTEGER4, gcount_group, 1, MPI_INTEGER4, this%comm, this%IERR)

    gcount_group = gcount_group * buf_dist_shape(1)

    disp_group = 0
    DO k1 = 2, this%nProc
      disp_group(k1) = gcount_group(k1 - 1) + disp_group(k1 - 1)
    END DO

    CALL this%barrier
    CALL mpi_gatherv(buf_src, gcount_group(this%myrank), MPI_DOUBLE_PRECISION, &
                     bufTemp, gcount_group, disp_group, MPI_DOUBLE_PRECISION, this%rankBase, this%comm, this%ierr)

    IF (this%isBaseProc()) THEN
      buf_dist(:, seq_store_idx) = bufTemp
    END IF

    DEALLOCATE (gcount_group, disp_group)
    IF (this%isBaseProc()) DEALLOCATE (bufTemp)
  END SUBROUTINE aggrGridReal2D

  SUBROUTINE aggrGridReal3D(this, nCell_perProc, buf_src, buf_dist, buf_dist_shape, seq_store_idx)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) this
    INTEGER(i_kind), INTENT(IN) :: buf_dist_shape(3), nCell_perProc, seq_store_idx(:)
    REAL(r_kind), INTENT(IN)    :: buf_src(:, :, :)
    REAL(r_kind), INTENT(INOUT) :: buf_dist(:, :, :)
    INTEGER, ALLOCATABLE :: gcount_group(:), disp_group(:)
    INTEGER(i_kind) :: k1
    REAL(r_kind), ALLOCATABLE :: bufTemp(:, :, :)

    ! Return if the proc is not active at this grid level
    IF (.NOT. this%isActiveProc()) RETURN

    IF (this%isBaseProc()) THEN
      ALLOCATE (bufTemp(buf_dist_shape(1), buf_dist_shape(2), buf_dist_shape(3)))
    END IF

    ALLOCATE (gcount_group(this%nProc), disp_group(this%nProc))
    CALL MPI_ALLGATHER(nCell_perProc, 1, MPI_INTEGER4, gcount_group, 1, MPI_INTEGER4, this%comm, this%IERR)

    gcount_group = gcount_group * buf_dist_shape(1) * buf_dist_shape(2)

    disp_group = 0
    DO k1 = 2, this%nProc
      disp_group(k1) = gcount_group(k1 - 1) + disp_group(k1 - 1)
    END DO

    CALL this%barrier
    CALL mpi_gatherv(buf_src, gcount_group(this%myrank), MPI_DOUBLE_PRECISION, &
                     bufTemp, gcount_group, disp_group, MPI_DOUBLE_PRECISION, this%rankBase, this%comm, this%ierr)

    IF (this%isBaseProc()) THEN
      buf_dist(:, :, seq_store_idx) = bufTemp
    END IF

    DEALLOCATE (gcount_group, disp_group)
    IF (this%isBaseProc()) DEALLOCATE (bufTemp)
  END SUBROUTINE aggrGridReal3D

  SUBROUTINE initializemMPDDSub(this, gLevel, proc_layout, group, rank)
    CLASS(mpddInfo_sg_t) this
    INTEGER(i_kind), INTENT(IN) :: proc_layout(2), rank, gLevel
    INTEGER(i_kind), INTENT(IN) :: group

    this%gLevel = gLevel
    this%proc_layout = proc_layout
    this%halo_size = 0

    PRINT *, 'Begin construction of mpddSub on G', gLevel
    CALL this%initialize(proc_layout(1) * proc_layout(2), group, rank)
    PRINT *, 'End construction of mpddSub on G', gLevel

    ! Yuanfu Xie added more explanation of the output:
    WRITE (*, 1) this%gLevel, rank, this%rank, this%comm, MPI_COMM_NULL
1   FORMAT('MPDDInfo_sg - Grid Level:', I2, ' My rank', 2I4, ' comm: ', I2, ' comm_null: ', I2)

  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(mpddInfo_sg_t), INTENT(INOUT) ::this
    INTEGER(i_kind) :: j

    CALL this%buf_dealloc

    IF (ALLOCATED(this%send_info%nb_blk)) THEN
      DO j = LBOUND(this%send_info%nb_blk, 1), UBOUND(this%send_info%nb_blk, 1)   ! The glvl starts from 2
        IF (ALLOCATED(this%send_info%nb_blk(j)%idx)) DEALLOCATE (this%send_info%nb_blk(j)%idx)
      END DO
      DEALLOCATE (this%send_info%nb_blk)
    END IF

    IF (ALLOCATED(this%recv_info%nb_blk)) THEN
      DO j = LBOUND(this%recv_info%nb_blk, 1), UBOUND(this%recv_info%nb_blk, 1)   ! The glvl starts from 2
        IF (ALLOCATED(this%recv_info%nb_blk(j)%idx)) DEALLOCATE (this%recv_info%nb_blk(j)%idx)
      END DO
      DEALLOCATE (this%recv_info%nb_blk)
    END IF

    ! CALL this%dest_proc_layout
    PRINT *, 'initialize destructor is invoked!'
  END SUBROUTINE destructor

  SUBROUTINE scatter_to_finer_threads(this, toFiner, toFinerSwap, nproc_layoutFiner, nproc_layoutCoarser, &
                                      vlevel, ncells_SgFiner, g_t_sp_idx_toCoarser, sp_t_g_idx_atFiner)
    CLASS(mpddInfo_sg_t) :: this
    INTEGER(i_kind), INTENT(IN) :: nproc_layoutFiner(2), nproc_layoutCoarser(2), vlevel, g_t_sp_idx_toCoarser(:), &
                                   sp_t_g_idx_atFiner(:), ncells_SgFiner
    REAL(r_kind), INTENT(OUT) :: toFiner(vlevel, *)
    REAL(r_kind), INTENT(INOUT) :: toFinerSwap(vlevel, *)
    INTEGER(i_kind) :: tRatio, i, dimRatio(2)
    INTEGER(i_kind), ALLOCATABLE :: req1(:), req2(:)
    INTEGER(i_kind), ALLOCATABLE :: ncells_SgFinerGroup(:), disp_SgFinerGroup(:), ncells_sgFinerTotal
    INTEGER(i_kind), ALLOCATABLE :: sp_t_g_idx_SgFinerGroup(:), nProcFiner, nProcCoarser
    REAL(r_kind), ALLOCATABLE :: toFinerSwap_SgGroup(:, :)
    INTEGER(i_kind) :: srcProcRank, rankRowIdxF, rankColIdxF, rankRowIdxC, rankColIdxC, idx, j
    INTEGER(i_kind), ALLOCATABLE :: distProcRank(:)

    IF (.NOT. this%isActiveProc()) RETURN

    CALL mpi_bcast(nproc_layoutFiner, 2, MPI_INTEGER4, this%rankBase, this%comm, this%ierr)
    CALL mpi_bcast(nproc_layoutCoarser, 2, MPI_INTEGER4, this%rankBase, this%comm, this%ierr)

    nProcFiner = nproc_layoutFiner(1) * nproc_layoutFiner(2)
    dimRatio(1) = nproc_layoutFiner(1) / nproc_layoutCoarser(1)
    dimRatio(2) = nproc_layoutFiner(2) / nproc_layoutCoarser(2)
    ALLOCATE (req1(nProcFiner), req2(nProcFiner))

    nProcCoarser = nproc_layoutCoarser(1) * nproc_layoutCoarser(2)
    tRatio = nProcFiner / nProcCoarser

    rankRowIdxF = (this%myrank - 1) / nproc_layoutFiner(2) + 1
    rankColIdxF = MOD(this%myrank - 1, nproc_layoutFiner(2)) + 1

    rankRowIdxC = (rankRowIdxF - 1) / dimRatio(1) + 1
    rankColIdxC = (rankColIdxF - 1) / dimRatio(2) + 1

    srcProcRank = (rankRowIdxC - 1) * nproc_layoutCoarser(2) + rankColIdxC
    srcProcRank = srcProcRank - 1

    IF (this%myrank <= nProcCoarser) THEN
      ALLOCATE (ncells_SgFinerGroup(tRatio), &
                disp_SgFinerGroup(tRatio), &
                distProcRank(tRatio))

      idx = 0
      rankRowIdxC = (this%myrank - 1) / nproc_layoutCoarser(2) + 1
      rankColIdxC = MOD(this%myrank - 1, nproc_layoutCoarser(2)) + 1

      DO i = 1, nproc_layoutFiner(1) / nproc_layoutCoarser(1)
        DO j = 1, nproc_layoutFiner(2) / nproc_layoutCoarser(2)
          idx = idx + 1
          rankRowIdxF = (rankRowIdxC - 1) * dimRatio(1) + i
          rankColIdxF = (rankColIdxC - 1) * dimRatio(2) + j
          distProcRank(idx) = (rankRowIdxF - 1) * nproc_layoutFiner(2) + rankColIdxF
        END DO
      END DO
      distProcRank = distProcRank - 1

      ! print *, 'myrank is', this%myrank, 'I send to', distProcRank
    END IF

    ! print *, 'myrank is', this%myrank, 'I receive from', srcProcRank

!!
    CALL mpi_isend(ncells_SgFiner, 1, MPI_INTEGER4, &
                   srcProcRank, 997, this%comm, req1(1), this%ierr)

    IF (this%myrank <= nProcCoarser) THEN
      DO i = 1, tRatio
        CALL mpi_irecv(ncells_SgFinerGroup(i), 1, MPI_INTEGER4, &
                       distProcRank(i), 997, this%comm, req2(i), this%ierr)
      END DO

      DO i = 1, tRatio
        CALL mpi_wait(req2(i), this%Status, this%ierr)
      END DO

      disp_SgFinerGroup = 0
      ncells_sgFinerTotal = ncells_SgFinerGroup(1)
      DO i = 2, tRatio
        disp_SgFinerGroup(i) = disp_SgFinerGroup(i - 1) + ncells_SgFinerGroup(i - 1)
        ncells_sgFinerTotal = ncells_sgFinerTotal + ncells_SgFinerGroup(i)
      END DO
      ALLOCATE (sp_t_g_idx_SgFinerGroup(ncells_sgFinerTotal))
      ALLOCATE (toFinerSwap_SgGroup(vlevel, ncells_sgFinerTotal))
    END IF

    CALL mpi_wait(req1(1), this%Status, this%ierr)
! !!
!         !!!!!!!!!!!!!!!!!!!!!!
    ! CALL MPI_BARRIER(this%comm, this%ierr)

    CALL mpi_isend(sp_t_g_idx_atFiner, ncells_SgFiner, &
                   MPI_INTEGER4, &
                   srcProcRank, 997, this%comm, req1(1), this%ierr)

    IF (this%myrank <= nProcCoarser) THEN
      DO i = 1, tRatio
        CALL mpi_irecv(sp_t_g_idx_SgFinerGroup(disp_SgFinerGroup(i) + 1), ncells_SgFinerGroup(i), &
                       MPI_INTEGER4, &
                       distProcRank(i), 997, this%comm, req2(i), this%ierr)
      END DO

      DO i = 1, tRatio
        CALL mpi_wait(req2(i), this%Status, this%ierr)
      END DO

    END IF
    CALL mpi_wait(req1(1), this%Status, this%ierr)

    IF (this%myrank <= nProcCoarser) THEN
      toFinerSwap_SgGroup = toFinerSwap(:, g_t_sp_idx_toCoarser(sp_t_g_idx_SgFinerGroup))
      ! print *, 'toFinerSwap_SgGroup', toFinerSwap_SgGroup
    END IF

    IF (this%myrank <= nProcCoarser) THEN
      DO i = 1, tRatio
        CALL mpi_isend(toFinerSwap_SgGroup(1, disp_SgFinerGroup(i) + 1), ncells_SgFinerGroup(i) * vlevel, &
                       MPI_DOUBLE_PRECISION, &
                       distProcRank(i), 997, this%comm, req1(i), this%ierr)
      END DO
    END IF

    CALL mpi_irecv(toFiner, ncells_SgFiner * vlevel, &
                   MPI_DOUBLE_PRECISION, &
                   srcProcRank, 997, this%comm, req2(1), this%ierr)

    CALL mpi_wait(req2(1), this%Status, this%ierr)

    IF (this%myrank <= nProcCoarser) THEN
      DO i = 1, tRatio
        CALL mpi_wait(req1(i), this%Status, this%ierr)
      END DO
    END IF

    !!!!!!!!!!!!!!!!!!!!!!
    DEALLOCATE (req1, req2)
    IF (this%myrank <= nProcCoarser) DEALLOCATE (ncells_SgFinerGroup, &
                                                 disp_SgFinerGroup, &
                                                 sp_t_g_idx_SgFinerGroup, &
                                                 toFinerSwap_SgGroup, &
                                                 distProcRank)
    CALL MPI_BARRIER(this%comm, this%ierr)
  END SUBROUTINE scatter_to_finer_threads

  SUBROUTINE gather_to_coarser_threads(this, toCoarse, toCoarseSwap, nproc_layoutFiner, nproc_layoutCoarser, &
                                       vlevel, ncells_SgCoarser, sp_t_g_idx_toFiner, g_t_sp_idx_atCoarser)
    CLASS(mpddInfo_sg_t) :: this
    INTEGER(i_kind), INTENT(IN) :: nproc_layoutFiner(2), nproc_layoutCoarser(2), vlevel, sp_t_g_idx_toFiner(:), &
                                   g_t_sp_idx_atCoarser(:), ncells_SgCoarser
    REAL(r_kind), INTENT(OUT) :: toCoarse(vlevel, *)
    REAL(r_kind), INTENT(INOUT) :: toCoarseSwap(vlevel, *)
    INTEGER(i_kind) :: tRatio, i, dimRatio(2)
    INTEGER(i_kind), ALLOCATABLE :: req1(:), req2(:)
    INTEGER(i_kind), ALLOCATABLE :: ncells_SgCoarserGroup(:), disp_SgCoarserGroup(:), ncells_sgCoarserTotal
    INTEGER(i_kind), ALLOCATABLE :: sp_t_g_idx_SgCoarserGroup(:), nProcFiner, nProcCoarser
    REAL(r_kind), ALLOCATABLE :: toCoarserSwap_SgGroup(:, :)
    INTEGER(i_kind) :: distProcRank, rankRowIdxF, rankColIdxF, rankRowIdxC, rankColIdxC, idx, j
    INTEGER(i_kind), ALLOCATABLE :: srcProcRank(:)

    IF (.NOT. this%isActiveProc()) RETURN

    CALL mpi_bcast(nproc_layoutFiner, 2, MPI_INTEGER4, this%rankBase, this%comm, this%ierr)
    CALL mpi_bcast(nproc_layoutCoarser, 2, MPI_INTEGER4, this%rankBase, this%comm, this%ierr)

    nProcFiner = nproc_layoutFiner(1) * nproc_layoutFiner(2)
    dimRatio(1) = nproc_layoutFiner(1) / nproc_layoutCoarser(1)
    dimRatio(2) = nproc_layoutFiner(2) / nproc_layoutCoarser(2)
    ALLOCATE (req1(nProcFiner), req2(nProcFiner))

    nProcCoarser = nproc_layoutCoarser(1) * nproc_layoutCoarser(2)
    tRatio = nProcFiner / nProcCoarser

    rankRowIdxF = (this%myrank - 1) / nproc_layoutFiner(2) + 1
    rankColIdxF = MOD(this%myrank - 1, nproc_layoutFiner(2)) + 1

    rankRowIdxC = (rankRowIdxF - 1) / dimRatio(1) + 1
    rankColIdxC = (rankColIdxF - 1) / dimRatio(2) + 1

    distProcRank = (rankRowIdxC - 1) * nproc_layoutCoarser(2) + rankColIdxC
    distProcRank = distProcRank - 1

    IF (this%myrank <= nProcCoarser) THEN
      ALLOCATE (ncells_SgCoarserGroup(tRatio), &
                disp_SgCoarserGroup(tRatio), &
                srcProcRank(tRatio))

      idx = 0
      rankRowIdxC = (this%myrank - 1) / nproc_layoutCoarser(2) + 1
      rankColIdxC = MOD(this%myrank - 1, nproc_layoutCoarser(2)) + 1

      DO i = 1, nproc_layoutFiner(1) / nproc_layoutCoarser(1)
        DO j = 1, nproc_layoutFiner(2) / nproc_layoutCoarser(2)
          idx = idx + 1
          rankRowIdxF = (rankRowIdxC - 1) * dimRatio(1) + i
          rankColIdxF = (rankColIdxC - 1) * dimRatio(2) + j
          srcProcRank(idx) = (rankRowIdxF - 1) * nproc_layoutFiner(2) + rankColIdxF
        END DO
      END DO
      srcProcRank = srcProcRank - 1

      ! print *, 'myrank is', this%myrank, 'I receive from', srcProcRank
    END IF

    ! print *, 'myrank is', this%myrank, 'I send to', distProcRank

    CALL mpi_isend(ncells_SgCoarser, 1, MPI_INTEGER4, &
                   distProcRank, 997, this%comm, req1(1), this%ierr)

    IF (this%myrank <= nProcCoarser) THEN
      DO i = 1, tRatio
        CALL mpi_irecv(ncells_SgCoarserGroup(i), 1, MPI_INTEGER4, &
                       srcProcRank(i), 997, this%comm, req2(i), this%ierr)
      END DO

      DO i = 1, tRatio
        CALL mpi_wait(req2(i), this%Status, this%ierr)
      END DO

      disp_SgCoarserGroup = 0
      ncells_sgCoarserTotal = ncells_SgCoarserGroup(1)
      DO i = 2, tRatio
        disp_SgCoarserGroup(i) = disp_SgCoarserGroup(i - 1) + ncells_SgCoarserGroup(i - 1)
        ncells_sgCoarserTotal = ncells_sgCoarserTotal + ncells_SgCoarserGroup(i)
      END DO
      ALLOCATE (sp_t_g_idx_SgCoarserGroup(ncells_sgCoarserTotal))
      ALLOCATE (toCoarserSwap_SgGroup(vlevel, ncells_sgCoarserTotal))

    END IF

    CALL mpi_wait(req1(1), this%Status, this%ierr)

    ! !!!!!!!!!!!!!!!!!!!!!!
    ! CALL MPI_BARRIER(this%comm, this%ierr)

    CALL mpi_isend(sp_t_g_idx_toFiner, ncells_SgCoarser, &
                   MPI_INTEGER4, &
                   distProcRank, 997, this%comm, req1(1), this%ierr)

    IF (this%myrank <= nProcCoarser) THEN
      DO i = 1, tRatio
        CALL mpi_irecv(sp_t_g_idx_SgCoarserGroup(disp_SgCoarserGroup(i) + 1), ncells_SgCoarserGroup(i), &
                       MPI_INTEGER4, &
                       srcProcRank(i), 997, this%comm, req2(i), this%ierr)
      END DO

      DO i = 1, tRatio
        CALL mpi_wait(req2(i), this%Status, this%ierr)
      END DO

    END IF
    CALL mpi_wait(req1(1), this%Status, this%ierr)

    ! !!!!!!!!!!!!!!!!!!!!!!
    ! CALL MPI_BARRIER(this%comm, this%ierr)

    CALL mpi_isend(toCoarseSwap, ncells_SgCoarser * vlevel, &
                   MPI_DOUBLE_PRECISION, &
                   distProcRank, 997, this%comm, req1(1), this%ierr)

    IF (this%myrank <= nProcCoarser) THEN
      DO i = 1, tRatio
        CALL mpi_irecv(toCoarserSwap_SgGroup(1, disp_SgCoarserGroup(i) + 1), ncells_SgCoarserGroup(i) * vlevel, &
                       MPI_DOUBLE_PRECISION, &
                       srcProcRank(i), 997, this%comm, req2(i), this%ierr)
      END DO

      DO i = 1, tRatio
        CALL mpi_wait(req2(i), this%Status, this%ierr)
      END DO

    END IF
    CALL mpi_wait(req1(1), this%Status, this%ierr)

    ! CALL MPI_BARRIER(this%comm, this%ierr)

    IF (this%myrank <= nProcCoarser) THEN

      ! PRINT *, 'myrank', this%myrank, g_t_sp_idx_atCoarser

      DO i = 1, ncells_sgCoarserTotal
        ! PRINT *, i, sp_t_g_idx_SgCoarserGroup(i), g_t_sp_idx_atCoarser(sp_t_g_idx_SgCoarserGroup(i))
        toCoarse(:, g_t_sp_idx_atCoarser(sp_t_g_idx_SgCoarserGroup(i))) = toCoarserSwap_SgGroup(:, i)
      END DO
    END IF

    DEALLOCATE (req1, req2)
    IF (this%myrank <= nProcCoarser) DEALLOCATE (ncells_SgCoarserGroup, &
                                                 disp_SgCoarserGroup, &
                                                 sp_t_g_idx_SgCoarserGroup, &
                                                 toCoarserSwap_SgGroup, &
                                                 srcProcRank)
    CALL MPI_BARRIER(this%comm, this%ierr)
  END SUBROUTINE gather_to_coarser_threads

  SUBROUTINE ExchangeMatOnHalo(this, num_icell, vlevel, exBuf, g_t_sp_idx)
    IMPLICIT NONE
    ! Multigrid  and vertical level
    CLASS(mpddInfo_sg_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vlevel, num_icell, g_t_sp_idx(:)
    ! REAL(r_kind), ALLOCATABLE, INTENT(INOUT) :: swapSENDbuf(:, :), swapRecvbuf(:, :)
    REAL(r_kind), INTENT(INOUT) :: exBuf(vlevel, *)
    INTEGER(i_kind) :: i, j, k
    INTEGER(i_kind), ALLOCATABLE :: req1(:), req2(:)
    REAL(r_kind), ALLOCATABLE :: swapSendbuf(:, :)
    ! REAL(r_kind), ALLOCATABLE :: swapRecvbuf(:, :)

    IF (.NOT. this%isActiveProc()) RETURN

    ALLOCATE (req1(this%nProc), req2(this%nProc))
    ALLOCATE (swapSendbuf(vlevel, this%halo_size))
    ! ALLOCATE (swapRecvbuf(vlevel, this%halo_size))

    k = 0
    swapSENDbuf = 0

    DO i = 1, this%send_info%nb_blk_sz
      DO j = 1, this%send_info%nb_blk(i)%sz
        k = k + 1
        swapSENDbuf(:, k) = exBuf(:, g_t_sp_idx(this%send_info%nb_blk(i)%idx(j)))
      END DO
    END DO

    k = 0
    DO i = 1, this%send_info%nb_blk_sz
      CALL mpi_isend(swapSENDbuf(:, k + 1), this%send_info%nb_blk(i)%sz * vlevel, &
                     MPI_DOUBLE_PRECISION, &
                     this%send_info%nb_blk(i)%blk_num - 1, 997, this%comm, req1(i), this%ierr)
      k = k + this%send_info%nb_blk(i)%sz
    END DO

    k = 0
    DO i = 1, this%recv_info%nb_blk_sz
      CALL mpi_irecv(exBuf(:, num_icell + k + 1), this%recv_info%nb_blk(i)%sz * vlevel, &
                     MPI_DOUBLE_PRECISION, &
                     this%recv_info%nb_blk(i)%blk_num - 1, 997, this%comm, req2(i), this%ierr)
      k = k + this%recv_info%nb_blk(i)%sz
    END DO

    DO i = 1, this%recv_info%nb_blk_sz
      CALL mpi_wait(req2(i), this%Status, this%ierr)
    END DO

    DO i = 1, this%send_info%nb_blk_sz
      CALL mpi_wait(req1(i), this%Status, this%ierr)
    END DO

    CALL MPI_BARRIER(this%comm, this%ierr)
    ! exBuf(:, num_icell + 1:num_icell + this%halo_size) = swapRecvbuf(:, 1:this%halo_size)
    DEALLOCATE (req1, req2)
    ! DEALLOCATE (swapRecvbuf)

    DEALLOCATE (swapSendbuf)
    ! CALL MPI_BARRIER(this%comm, this%ierr)

  END SUBROUTINE ExchangeMatOnHalo

  SUBROUTINE ExchangeMatOnHaloReverseSum(this, num_icell, vlevel, exBuf, g_t_sp_idx)
    IMPLICIT NONE
    ! Multigrid  and vertical level
    CLASS(mpddInfo_sg_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vlevel, num_icell, g_t_sp_idx(:)
    ! REAL(r_kind), ALLOCATABLE, INTENT(INOUT) :: swapSENDbuf(:, :), swapRecvbuf(:, :)
    REAL(r_kind), INTENT(INOUT) :: exBuf(vlevel, *)
    INTEGER(i_kind) :: i, j, k
    INTEGER(i_kind), ALLOCATABLE :: req1(:), req2(:)
    REAL(r_kind), ALLOCATABLE :: swapSendbuf(:, :)
    ! REAL(r_kind), ALLOCATABLE :: swapRecvbuf(:, :)

    IF (.NOT. this%isActiveProc()) RETURN

    ALLOCATE (req1(this%nProc), req2(this%nProc))
    ALLOCATE (swapSendbuf(vlevel, this%halo_size))
    ! ALLOCATE (swapRecvbuf(vlevel, this%halo_size))

    k = 0
    DO i = 1, this%recv_info%nb_blk_sz
      CALL mpi_isend(exBuf(:, num_icell + k + 1), this%recv_info%nb_blk(i)%sz * vlevel, &
                     MPI_DOUBLE_PRECISION, &
                     this%recv_info%nb_blk(i)%blk_num - 1, 997, this%comm, req2(i), this%ierr)
      k = k + this%recv_info%nb_blk(i)%sz
    END DO

    k = 0
    DO i = 1, this%send_info%nb_blk_sz
      CALL mpi_irecv(swapSENDbuf(:, k + 1), this%send_info%nb_blk(i)%sz * vlevel, &
                     MPI_DOUBLE_PRECISION, &
                     this%send_info%nb_blk(i)%blk_num - 1, 997, this%comm, req1(i), this%ierr)
      k = k + this%send_info%nb_blk(i)%sz
    END DO

    DO i = 1, this%send_info%nb_blk_sz
      CALL mpi_wait(req1(i), this%Status, this%ierr)
    END DO

    DO i = 1, this%recv_info%nb_blk_sz
      CALL mpi_wait(req2(i), this%Status, this%ierr)
    END DO

    CALL MPI_BARRIER(this%comm, this%ierr)

    k = 0
    DO i = 1, this%send_info%nb_blk_sz
      DO j = 1, this%send_info%nb_blk(i)%sz
        k = k + 1
        exBuf(:, g_t_sp_idx(this%send_info%nb_blk(i)%idx(j))) = exBuf(:, g_t_sp_idx(this%send_info%nb_blk(i)%idx(j))) + swapSENDbuf(:, k)
      END DO
    END DO

    DEALLOCATE (req1, req2)
    ! DEALLOCATE (swapRecvbuf)

    DEALLOCATE (swapSendbuf)
    !

  END SUBROUTINE ExchangeMatOnHaloReverseSum

  SUBROUTINE genCellType(this, cell_type, edge_stcl, num_icell, sp_t_g_idx, dimCell_global)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) :: this
    INTEGER(i_kind), INTENT(OUT) :: cell_type(:)
    INTEGER(i_kind), INTENT(IN) :: edge_stcl(:, :, :), num_icell, sp_t_g_idx(:), dimCell_global(2)
    INTEGER(i_kind) :: i, idx_global, i_rowIdx, i_ColIdx

    IF (.NOT. this%isActiveProc()) RETURN

    PRINT *, UBOUND(cell_type, 1)
    cell_type = 0
    IF (num_icell .LE. UBOUND(cell_type, 1)) cell_type(num_icell + 1:UBOUND(cell_type, 1)) = 3

    DO i = 1, num_icell
      idx_global = sp_t_g_idx(i)
      i_rowIdx = (idx_global - 1) / dimCell_global(2) + 1
      i_ColIdx = MOD(idx_global - 1, dimcell_global(2)) + 1

      IF ((i_rowIdx .EQ. 1) .OR. (i_ColIdx .EQ. 1) .OR. &
          (i_rowIdx .EQ. dimcell_global(1)) .OR. (i_ColIdx .EQ. dimcell_global(2))) THEN
        cell_type(i) = 1
      END IF

      IF (((i_rowIdx .EQ. 1) .AND. (i_ColIdx .EQ. 1)) .OR. &
          ((i_rowIdx .EQ. 1) .AND. (i_ColIdx .EQ. dimCell_global(2))) .OR. &
          ((i_rowIdx .EQ. dimCell_global(1)) .AND. (i_ColIdx .EQ. 1)) .OR. &
          ((i_rowIdx .EQ. dimCell_global(1)) .AND. (i_ColIdx .EQ. dimCell_global(2)))) THEN
        cell_type(i) = 2
        CYCLE
      END IF
    END DO
  END SUBROUTINE genCellType

  !> @brief
!! This function remap the both the inner cell and boundary cell indices.
! @see
! @note
! @warning
! @attention
  SUBROUTINE RemapAllGridIdxToSP2D(this, matrix, g_t_sp_idx, dimCell_global)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) :: this
    INTEGER(i_kind), INTENT(IN) :: g_t_sp_idx(:), dimCell_global(2)
    INTEGER(i_kind), INTENT(INOUT) :: matrix(:, :)
    INTEGER :: i, j, idx_sp

    IF (.NOT. this%isActiveProc()) RETURN
    DO j = LBOUND(matrix, 2), UBOUND(matrix, 2)
      DO i = LBOUND(matrix, 1), UBOUND(matrix, 1)
        ! IF (matrix(i, j) .eq. 0) THEN ! Yuanfu Xie changed this to eliminate negative and o
        IF (matrix(i, j) .LE. 0) THEN
          idx_sp = 0
        ELSE
          idx_sp = g_t_sp_idx(matrix(i, j))
          IF (idx_sp .EQ. 0) THEN
            idx_sp = this%halo_stcl_idx_remap(matrix(i, j), this%recv_info, dimCell_global)
          END IF
        END IF
        matrix(i, j) = idx_sp
      END DO
    END DO
  END SUBROUTINE RemapAllGridIdxToSP2D

  !> @brief
!! This function remap the both the inner cell and boundary cell indices.
! @see
! @note
! @warning
! @attention
! Yuanfu Xie added a new routine for remapping all grid indices including the boundary 2025-02-07:
  SUBROUTINE RemapAllGridIdxToSP3D(this, matrix, g_t_sp_idx, dimCell_global)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) :: this
    INTEGER(i_kind), INTENT(IN) :: g_t_sp_idx(:), dimCell_global(2)
    INTEGER(i_kind), INTENT(INOUT) :: matrix(:, :, :)
    INTEGER :: i, j, k, idx_sp

    IF (.NOT. this%isActiveProc()) RETURN
    DO k = LBOUND(matrix, 3), UBOUND(matrix, 3)
      DO j = LBOUND(matrix, 2), UBOUND(matrix, 2)
        DO i = LBOUND(matrix, 1), UBOUND(matrix, 1)
          ! IF (matrix(i, j) .eq. 0) THEN ! Yuanfu Xie changed this to eliminate negative and o
          IF (matrix(i, j, k) .LE. 0) THEN
            idx_sp = 0
          ELSE
            idx_sp = g_t_sp_idx(matrix(i, j, k))
            IF (idx_sp .EQ. 0) THEN
              idx_sp = this%halo_stcl_idx_remap(matrix(i, j, k), this%recv_info, dimCell_global)
            END IF
          END IF
          matrix(i, j, k) = idx_sp
        END DO
      END DO
    END DO
  END SUBROUTINE RemapAllGridIdxToSP3D

!> @brief
!! This function remap the inner cells only.
! @see
! @note
! @warning
! @attention
  SUBROUTINE RemapInnerGridIdxToSP3D(this, matrix, g_t_sp_idx, cell_type, dimCell_global)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) :: this
    INTEGER(i_kind), INTENT(IN) :: g_t_sp_idx(:), cell_type(:), dimCell_global(2)
    INTEGER(i_kind), INTENT(INOUT) :: matrix(:, :, :)
    INTEGER :: i, j, k, idx_sp

    IF (.NOT. this%isActiveProc()) RETURN

    DO k = LBOUND(matrix, 3), UBOUND(matrix, 3)
      IF (cell_type(k) .NE. 0) CYCLE
      DO j = LBOUND(matrix, 2), UBOUND(matrix, 2)
        DO i = LBOUND(matrix, 1), UBOUND(matrix, 1)
          idx_sp = g_t_sp_idx(matrix(i, j, k))
          IF (idx_sp .EQ. 0) THEN
            idx_sp = this%halo_stcl_idx_remap(matrix(i, j, k), this%recv_info, dimCell_global)
          END IF
          matrix(i, j, k) = idx_sp
        END DO
      END DO
    END DO
  END SUBROUTINE RemapInnerGridIdxToSP3D

  INTEGER FUNCTION halo_stcl_idx_remap(this, idx, hBlock, dimCell_global)
    CLASS(mpddInfo_sg_t) :: this
    INTEGER(i_kind), INTENT(IN) :: idx, dimCell_global(2)
    TYPE(HaloBlock), INTENT(in)::hBlock
    INTEGER :: idxCounter, j, k, blk

    IF (.NOT. this%isActiveProc()) RETURN

    blk = this%i_global2ProcRank(idx, dimCell_global)
    IF (blk == this%myrank) THEN
      PRINT *, 'error indx in stcl 1', idx, this%myrank
      STOP
      RETURN
    END IF

    halo_stcl_idx_remap = 0
    idxCounter = 0
    DO j = 1, hBlock%nb_blk_sz
      IF (blk == hBlock%nb_blk(j)%blk_num) THEN
        DO k = 1, hBlock%nb_blk(j)%sz
          IF (hBlock%nb_blk(j)%idx(k) == idx) THEN
            halo_stcl_idx_remap = idxCounter + k + hBlock%num_icell
            RETURN
          END IF
        END DO
      END IF
      idxCounter = idxCounter + hBlock%nb_blk(j)%sz
    END DO

    IF (halo_stcl_idx_remap == 0) THEN
      PRINT *, 'error index in stcl 2', idx, ' pc: ',this%myrank
      STOP
    END IF
  END FUNCTION halo_stcl_idx_remap

  SUBROUTINE generate_halo_info(this, cell_stcl, num_iCell, num_cell_stcl, dimCell_global) ! ncp is limited at this point 4^p*20, p = 0, 1, ...
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) :: this
    ! INTEGER, intent(in) :: adjnb(:, :)
    INTEGER(i_kind), INTENT(in) :: cell_stcl(:, :), num_iCell, num_cell_stcl, dimCell_global(2)

    IF (.NOT. this%isActiveProc()) RETURN

    CALL this%generate_recv_info(cell_stcl, num_icell, num_cell_stcl, dimCell_global)
    CALL this%generate_send_info
    CALL this%trim_idx_all_halos

    IF (this%recv_info%nb_blk_sz > max_nb) &
      PRINT *, this%myrank, this%gLevel, 'recv nb_blk_sz reALLOCATEd', this%recv_info%nb_blk_sz
    IF (this%send_info%nb_blk_sz > max_nb) &
      PRINT *, this%myrank, this%gLevel, 'send nb_blk_sz reALLOCATEd', this%send_info%nb_blk_sz

  END SUBROUTINE generate_halo_info

  INTEGER FUNCTION i_global2ProcRank(this, i, dimCell_global)
    CLASS(mpddInfo_sg_t) :: this
    INTEGER(i_kind), INTENT(IN) :: i, dimCell_global(2)
    INTEGER(i_kind) :: i_rowIdx, i_ColIdx, procRank_rowIdx, procRank_colIdx

    i_rowIdx = (i - 1) / dimCell_global(2) + 1
    i_ColIdx = MOD(i - 1, dimcell_global(2)) + 1

    IF (i_rowIdx .EQ. 1) i_rowIdx = i_rowIdx + 1
    IF (i_ColIdx .EQ. 1) i_ColIdx = i_ColIdx + 1
    IF (i_rowIdx .EQ. (dimCell_global(1))) i_rowIdx = i_rowIdx - 1
    IF (i_ColIdx .EQ. (dimCell_global(2))) i_ColIdx = i_ColIdx - 1

    procRank_rowIdx = (i_rowIdx - 1 - 1) / ((dimCell_global(1) - 2) / this%proc_layout(1)) + 1
    procRank_colIdx = (i_ColIdx - 1 - 1) / ((dimCell_global(2) - 2) / this%proc_layout(2)) + 1

    ! PRINT *, 'procRank', i_rowIdx, i_ColIdx, procRank_rowIdx, procRank_colIdx, this%dimCell_global
    i_global2ProcRank = (procRank_rowIdx - 1) * this%proc_layout(2) + procRank_colIdx
  END FUNCTION i_global2ProcRank

  SUBROUTINE calRowColRangeSP(rowRange, colRange, myrank, proc_layout, dimCell_global)
    INTEGER(i_kind), INTENT(in) :: myrank, proc_layout(2), dimCell_global(2)
    INTEGER(i_kind), INTENT(out):: rowRange(2), colRange(2)
    INTEGER(i_kind) :: i, j, k, myrank_i, myrank_j

    rowRange = 0
    colRange = 0

    myrank_i = (myrank - 1) / proc_layout(2) + 1
    myrank_j = MOD(myrank - 1, proc_layout(2)) + 1

    rowRange(1) = 1 + (dimCell_global(1) - 2) / proc_layout(1) * (myrank_i - 1) + 1
    rowRange(2) = 1 + (dimCell_global(1) - 2) / proc_layout(1) * (myrank_i)
    colRange(1) = 1 + (dimCell_global(2) - 2) / proc_layout(2) * (myrank_j - 1) + 1
    colRange(2) = 1 + (dimCell_global(2) - 2) / proc_layout(2) * (myrank_j)

    IF (rowRange(1) .EQ. 2) rowRange(1) = rowRange(1) - 1
    IF (colRange(1) .EQ. 2) colRange(1) = colRange(1) - 1
    IF (rowRange(2) .EQ. (dimCell_global(1) - 1)) rowRange(2) = rowRange(2) + 1
    IF (colRange(2) .EQ. (dimCell_global(2) - 1)) colRange(2) = colRange(2) + 1

    ! Check the domain decomposition:
    WRITE (*, 1) myrank, myrank_i, myrank_j, rowRange, colRange, dimCell_global
1   FORMAT('MPDDInfo_sg - Row/Col proc: ', I2, ' procIJ: ', 2I2, ' RowColRange: ', 4I4, ' cellGlobal: ', 2I4)

  END SUBROUTINE calRowColRangeSP

  SUBROUTINE genMapIdx(this, sp_t_g_idx, g_t_sp_idx, seq_store_idx, dimCell_global, dimCell_sg)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) this
    INTEGER(i_kind), INTENT(IN) :: dimCell_global(2)
    INTEGER(i_kind), INTENT(out) :: sp_t_g_idx(:), g_t_sp_idx(:), seq_store_idx(:)
    INTEGER(i_kind) :: i, j, rowRange(2), colRange(2), idxSP, idxWorld
    INTEGER(i_kind), ALLOCATABLE :: gcount_group(:), disp_group(:)
    INTEGER(i_kind), OPTIONAL, INTENT(OUT) :: dimCell_sg(2)
    ! Return if the proc is not active at this grid level
    IF (.NOT. this%isActiveProc()) RETURN

    CALL this%calRowColRangeSP(rowRange, colRange, this%myrank, this%proc_layout, dimCell_global)
    PRINT *, 'End the calRowColRangeSP on grid: ', this%gLevel, 'on proc: ', this%myrank

    idxSp = 0
    idxWorld = 0
    DO i = 1, dimCell_global(1)
      DO j = 1, dimCell_global(2)
        idxWorld = idxWorld + 1
        g_t_sp_idx(idxWorld) = 0
        IF ((i .GE. rowRange(1)) .AND. &
            (i .LE. rowRange(2)) .AND. &
            (j .GE. colRange(1)) .AND. &
            (j .LE. colRange(2))) THEN
          idxSp = idxSp + 1
          sp_t_g_idx(idxSp) = idxWorld
          g_t_sp_idx(idxWorld) = idxSp
        END IF
      END DO
    END DO

    IF (PRESENT(dimCell_sg)) THEN
      dimCell_sg(1) = rowRange(2) - rowRange(1) + 1
      dimCell_sg(2) = colRange(2) - colRange(1) + 1
    END IF

    PRINT *, 'End the idxCal on grid: ', this%gLevel, 'on proc: ', this%myrank

    ! CALL this%barrier
    ALLOCATE (gcount_group(this%nProc), disp_group(this%nProc))
    CALL MPI_ALLGATHER(idxSp, 1, MPI_INTEGER4, gcount_group, 1, MPI_INTEGER4, this%comm, this%IERR)

    ! CALL this%bCast(gcount_group)
    disp_group = 0
    DO i = 2, this%nProc
      disp_group(i) = gcount_group(i - 1) + disp_group(i - 1)
    END DO

    ! CALL this%barrier
    PRINT *, 'Before the MPI_GATHER in genMapIdx on grid: ', this%gLevel, 'on proc: ', &
      this%myrank, SUM(gcount_group), this%rankBase, SIZE(sp_t_g_idx), idxSp, SIZE(seq_store_idx)

    ! Gen seq_store_idx
    CALL MPI_GATHERV(sp_t_g_idx, idxSp, MPI_INTEGER4, seq_store_idx, &
                     gcount_group, disp_group, MPI_INTEGER4, this%rankBase, this%comm, this%IERR)

    DEALLOCATE (gcount_group, disp_group)
  END SUBROUTINE genMapIdx

  SUBROUTINE buf_dealloc(this)
    CLASS(mpddInfo_sg_t) this

    IF (ALLOCATED(this%sendbuf)) DEALLOCATE (this%sendbuf)
    IF (ALLOCATED(this%recvbuf)) DEALLOCATE (this%recvbuf)
    IF (ALLOCATED(this%swapSendbuf)) DEALLOCATE (this%swapSENDbuf)
    IF (ALLOCATED(this%swapRecvbuf)) DEALLOCATE (this%swapRecvbuf)

  END SUBROUTINE buf_dealloc

  SUBROUTINE buf_alloc(this, num_vlevels)
    CLASS(mpddInfo_sg_t) this
    INTEGER, INTENT(IN) :: num_vlevels
    INTEGER :: i

    ! Return if the proc is not active at this grid level
    IF (.NOT. this%isActiveProc()) RETURN

    CALL this%buf_dealloc

    ALLOCATE (this%sendbuf(num_vlevels, this%send_info%halo_sz))
    ALLOCATE (this%recvbuf(num_vlevels, this%recv_info%halo_sz))
  END SUBROUTINE buf_alloc

  SUBROUTINE calCellSize(this, dimCell_global, num_iCell)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) :: this
    INTEGER(i_kind), INTENT(IN) :: dimCell_global(2)
    INTEGER(i_kind), INTENT(OUT) :: num_iCell
    INTEGER(i_kind) :: i, j, rowRange(2), colRange(2)

    num_iCell = 0
    ! Return if the proc is not active at this grid level
    IF (.NOT. this%isActiveProc()) RETURN

    CALL this%calRowColRangeSP(rowRange, colRange, this%myrank, this%proc_layout, dimCell_global)

    DO i = 1, dimCell_global(1)
      DO j = 1, dimCell_global(2)
        IF ((i .GE. rowRange(1)) .AND. &
            (i .LE. rowRange(2)) .AND. &
            (j .GE. colRange(1)) .AND. &
            (j .LE. colRange(2))) THEN
          num_iCell = num_iCell + 1
        END IF
      END DO
    END DO

  END SUBROUTINE calCellSize

  SUBROUTINE DistGridRealSeq4D(this, nCell_perProc, buf_src, buf_dist, buf_src_shape, seq_store_idx)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) this
    INTEGER(i_kind), INTENT(IN) :: buf_src_shape(4), nCell_perProc, seq_store_idx(:)
    REAL(r_kind), INTENT(IN)    :: buf_src(:, :, :, :)
    REAL(r_kind), INTENT(INOUT) :: buf_dist(:, :, :, :)
    INTEGER(i_kind), ALLOCATABLE :: gcount_group(:), disp_group(:)
    INTEGER(i_kind) :: k1
    REAL(r_kind), ALLOCATABLE :: bufTemp(:, :, :, :)

    ! Return if the proc is not active at this grid level
    IF (.NOT. this%isActiveProc()) RETURN

    IF (this%isBaseProc()) THEN
      ALLOCATE (bufTemp(buf_src_shape(1), buf_src_shape(2), buf_src_shape(3), buf_src_shape(4)))
      bufTemp = buf_src(:, :, :, seq_store_idx)
    END IF

    ALLOCATE (gcount_group(this%nProc), disp_group(this%nProc))
    CALL MPI_ALLGATHER(nCell_perProc, 1, MPI_INTEGER4, gcount_group, 1, MPI_INTEGER4, this%comm, this%IERR)

    gcount_group = gcount_group * buf_src_shape(1) * buf_src_shape(2) * buf_src_shape(3)

    disp_group = 0

    DO k1 = 2, this%nProc
      disp_group(k1) = gcount_group(k1 - 1) + disp_group(k1 - 1)
    END DO

    CALL this%barrier
    CALL mpi_scatterv(bufTemp, gcount_group, disp_group, MPI_DOUBLE_PRECISION, &
                      buf_dist, gcount_group(this%myrank), MPI_DOUBLE_PRECISION, this%rankBase, this%comm, this%ierr)

    DEALLOCATE (gcount_group, disp_group)
    IF (this%isBaseProc()) DEALLOCATE (bufTemp)
  END SUBROUTINE DistGridRealSeq4D

  SUBROUTINE DistGridRealSeq3D(this, nCell_perProc, buf_src, buf_dist, buf_src_shape, seq_store_idx)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) this
    INTEGER(i_kind), INTENT(IN) :: buf_src_shape(3), nCell_perProc, seq_store_idx(:)
    REAL(r_kind), INTENT(IN)    :: buf_src(:, :, :)
    REAL(r_kind), INTENT(INOUT) :: buf_dist(:, :, :)
    INTEGER, ALLOCATABLE :: gcount_group(:), disp_group(:)
    INTEGER(i_kind) :: k1
    REAL(r_kind), ALLOCATABLE :: bufTemp(:, :, :)

    ! Return if the proc is not active at this grid level
    IF (.NOT. this%isActiveProc()) RETURN

    IF (this%isBaseProc()) THEN
      ALLOCATE (bufTemp(buf_src_shape(1), buf_src_shape(2), buf_src_shape(3)))
      bufTemp = buf_src(:, :, seq_store_idx)
    END IF

    ALLOCATE (gcount_group(this%nProc), disp_group(this%nProc))
    CALL MPI_ALLGATHER(nCell_perProc, 1, MPI_INTEGER4, gcount_group, 1, MPI_INTEGER4, this%comm, this%IERR)

    gcount_group = gcount_group * buf_src_shape(1) * buf_src_shape(2)

    disp_group = 0

    DO k1 = 2, this%nProc
      disp_group(k1) = gcount_group(k1 - 1) + disp_group(k1 - 1)
    END DO

    CALL this%barrier
    CALL mpi_scatterv(bufTemp, gcount_group, disp_group, MPI_DOUBLE_PRECISION, &
                      buf_dist, gcount_group(this%myrank), MPI_DOUBLE_PRECISION, this%rankBase, this%comm, this%ierr)

    DEALLOCATE (gcount_group, disp_group)
    IF (this%isBaseProc()) DEALLOCATE (bufTemp)
  END SUBROUTINE DistGridRealSeq3D

  SUBROUTINE DistGridRealSeq2D(this, nCell_perProc, buf_src, buf_dist, buf_src_shape, seq_store_idx)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) this
    INTEGER(i_kind), INTENT(IN) :: buf_src_shape(2), nCell_perProc, seq_store_idx(:)
    REAL(r_kind), INTENT(IN)    :: buf_src(:, :)
    REAL(r_kind), INTENT(INOUT) :: buf_dist(:, :)
    INTEGER, ALLOCATABLE :: gcount_group(:), disp_group(:)
    INTEGER(i_kind) :: k1
    REAL(r_kind), ALLOCATABLE :: bufTemp(:, :)

    ! Return if the proc is not active at this grid level
    IF (.NOT. this%isActiveProc()) RETURN

    IF (this%isBaseProc()) THEN
      ALLOCATE (bufTemp(buf_src_shape(1), buf_src_shape(2)))
      bufTemp = buf_src(:, seq_store_idx)
    END IF

    ALLOCATE (gcount_group(this%nProc), disp_group(this%nProc))
    CALL MPI_ALLGATHER(nCell_perProc, 1, MPI_INTEGER4, gcount_group, 1, MPI_INTEGER4, this%comm, this%IERR)

    gcount_group = gcount_group * buf_src_shape(1)

    disp_group = 0

    DO k1 = 2, this%nProc
      disp_group(k1) = gcount_group(k1 - 1) + disp_group(k1 - 1)
    END DO

    CALL this%barrier
    CALL mpi_scatterv(bufTemp, gcount_group, disp_group, MPI_DOUBLE_PRECISION, &
                      buf_dist, gcount_group(this%myrank), MPI_DOUBLE_PRECISION, this%rankBase, this%comm, this%ierr)

    DEALLOCATE (gcount_group, disp_group)
    IF (this%isBaseProc()) DEALLOCATE (bufTemp)

  END SUBROUTINE DistGridRealSeq2D

  SUBROUTINE DistGridRealSeq1D(this, nCell_perProc, buf_src, buf_dist, buf_src_shape, seq_store_idx)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) this
    INTEGER(i_kind), INTENT(IN) :: buf_src_shape(1), nCell_perProc, seq_store_idx(:)
    REAL(r_kind), INTENT(IN)    :: buf_src(:)
    REAL(r_kind), INTENT(INOUT) :: buf_dist(:)
    INTEGER, ALLOCATABLE :: gcount_group(:), disp_group(:)
    INTEGER(i_kind) :: k1
    REAL(r_kind), ALLOCATABLE :: bufTemp(:)

    ! Return if the proc is not active at this grid level
    IF (.NOT. this%isActiveProc()) RETURN

    IF (this%isBaseProc()) THEN
      ALLOCATE (bufTemp(buf_src_shape(1)))
      bufTemp = buf_src(seq_store_idx)
    END IF

    ALLOCATE (gcount_group(this%nProc), disp_group(this%nProc))
    CALL MPI_ALLGATHER(nCell_perProc, 1, MPI_INTEGER4, gcount_group, 1, MPI_INTEGER4, this%comm, this%IERR)

    gcount_group = gcount_group

    disp_group = 0

    DO k1 = 2, this%nProc
      disp_group(k1) = gcount_group(k1 - 1) + disp_group(k1 - 1)
    END DO

    CALL this%barrier
    CALL mpi_scatterv(bufTemp, gcount_group, disp_group, MPI_DOUBLE_PRECISION, &
                      buf_dist, gcount_group(this%myrank), MPI_DOUBLE_PRECISION, this%rankBase, this%comm, this%ierr)

    DEALLOCATE (gcount_group, disp_group)
    IF (this%isBaseProc()) DEALLOCATE (bufTemp)
  END SUBROUTINE DistGridRealSeq1D

  SUBROUTINE DistGridIntSeq4D(this, nCell_perProc, buf_src, buf_dist, buf_src_shape, seq_store_idx)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) this
    INTEGER(i_kind), INTENT(IN) :: buf_src_shape(4), nCell_perProc, buf_src(:, :, :, :), seq_store_idx(:)
    INTEGER(i_kind), INTENT(INOUT) :: buf_dist(:, :, :, :)
    INTEGER, ALLOCATABLE :: gcount_group(:), disp_group(:)
    INTEGER(i_kind) :: k1
    INTEGER(i_kind), ALLOCATABLE :: bufTemp(:, :, :, :)

    ! Return if the proc is not active at this grid level
    IF (.NOT. this%isActiveProc()) RETURN

    IF (this%isBaseProc()) THEN
      ALLOCATE (bufTemp(buf_src_shape(1), buf_src_shape(2), buf_src_shape(3), buf_src_shape(4)))
      bufTemp = buf_src(:, :, :, seq_store_idx)
    END IF

    ALLOCATE (gcount_group(this%nProc), disp_group(this%nProc))
    CALL MPI_ALLGATHER(nCell_perProc, 1, MPI_INTEGER4, gcount_group, 1, MPI_INTEGER4, this%comm, this%IERR)

    gcount_group = gcount_group * buf_src_shape(1) * buf_src_shape(2) * buf_src_shape(3)

    disp_group = 0

    DO k1 = 2, this%nProc
      disp_group(k1) = gcount_group(k1 - 1) + disp_group(k1 - 1)
    END DO

    CALL this%barrier
    CALL mpi_scatterv(bufTemp, gcount_group, disp_group, MPI_INTEGER4, &
                      buf_dist, gcount_group(this%myrank), MPI_INTEGER4, this%rankBase, this%comm, this%ierr)

    DEALLOCATE (gcount_group, disp_group)
    IF (this%isBaseProc()) DEALLOCATE (bufTemp)
  END SUBROUTINE DistGridIntSeq4D

  SUBROUTINE DistGridIntSeq3D(this, nCell_perProc, buf_src, buf_dist, buf_src_shape, seq_store_idx)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) this
    INTEGER(i_kind), INTENT(IN) :: buf_src_shape(3), nCell_perProc, buf_src(:, :, :), seq_store_idx(:)
    INTEGER(i_kind), INTENT(INOUT) :: buf_dist(:, :, :)
    INTEGER, ALLOCATABLE :: gcount_group(:), disp_group(:)
    INTEGER(i_kind) :: k1
    INTEGER(i_kind), ALLOCATABLE :: bufTemp(:, :, :)

    ! Return if the proc is not active at this grid level
    IF (.NOT. this%isActiveProc()) RETURN

    IF (this%isBaseProc()) THEN
      ALLOCATE (bufTemp(buf_src_shape(1), buf_src_shape(2), buf_src_shape(3)))
      bufTemp = buf_src(:, :, seq_store_idx)
    END IF

    ALLOCATE (gcount_group(this%nProc), disp_group(this%nProc))
    CALL MPI_ALLGATHER(nCell_perProc, 1, MPI_INTEGER4, gcount_group, 1, MPI_INTEGER4, this%comm, this%IERR)

    gcount_group = gcount_group * buf_src_shape(1) * buf_src_shape(2)

    disp_group = 0

    DO k1 = 2, this%nProc
      disp_group(k1) = gcount_group(k1 - 1) + disp_group(k1 - 1)
    END DO

    CALL this%barrier
    CALL mpi_scatterv(bufTemp, gcount_group, disp_group, MPI_INTEGER4, &
                      buf_dist, gcount_group(this%myrank), MPI_INTEGER4, this%rankBase, this%comm, this%ierr)

    DEALLOCATE (gcount_group, disp_group)
    IF (this%isBaseProc()) DEALLOCATE (bufTemp)

  END SUBROUTINE DistGridIntSeq3D

  SUBROUTINE DistGridIntSeq2D(this, nCell_perProc, buf_src, buf_dist, buf_src_shape, seq_store_idx)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) this
    INTEGER(i_kind), INTENT(IN) :: buf_src_shape(2), nCell_perProc, buf_src(:, :), seq_store_idx(:)
    INTEGER(i_kind), INTENT(INOUT) :: buf_dist(:, :)
    INTEGER, ALLOCATABLE :: gcount_group(:), disp_group(:)
    INTEGER(i_kind) :: k1
    INTEGER(i_kind), ALLOCATABLE :: bufTemp(:, :)

    ! Return if the proc is not active at this grid level
    IF (.NOT. this%isActiveProc()) RETURN

    IF (this%isBaseProc()) THEN
      ALLOCATE (bufTemp(buf_src_shape(1), buf_src_shape(2)))
      bufTemp = buf_src(:, seq_store_idx)
    END IF

    ALLOCATE (gcount_group(this%nProc), disp_group(this%nProc))
    CALL MPI_ALLGATHER(nCell_perProc, 1, MPI_INTEGER4, gcount_group, 1, MPI_INTEGER4, this%comm, this%IERR)

    gcount_group = gcount_group * buf_src_shape(1)

    disp_group = 0

    DO k1 = 2, this%nProc
      disp_group(k1) = gcount_group(k1 - 1) + disp_group(k1 - 1)
    END DO

    CALL this%barrier
    CALL mpi_scatterv(bufTemp, gcount_group, disp_group, MPI_INTEGER4, &
                      buf_dist, gcount_group(this%myrank), MPI_INTEGER4, this%rankBase, this%comm, this%ierr)

    DEALLOCATE (gcount_group, disp_group)
    IF (this%isBaseProc()) DEALLOCATE (bufTemp)
  END SUBROUTINE DistGridIntSeq2D

  SUBROUTINE DistGridIntSeq1D(this, nCell_perProc, buf_src, buf_dist, buf_src_shape, seq_store_idx)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) this
    INTEGER(i_kind), INTENT(IN) :: buf_src_shape(1), nCell_perProc, buf_src(:), seq_store_idx(:)
    INTEGER(i_kind), INTENT(INOUT) :: buf_dist(:)
    INTEGER, ALLOCATABLE :: gcount_group(:), disp_group(:)
    INTEGER(i_kind) :: k1
    INTEGER(i_kind), ALLOCATABLE :: bufTemp(:)

    ! Return if the proc is not active at this grid level
    IF (.NOT. this%isActiveProc()) RETURN

    IF (this%isBaseProc()) THEN
      ALLOCATE (bufTemp(buf_src_shape(1)))
      bufTemp = buf_src(seq_store_idx)
    END IF

    ALLOCATE (gcount_group(this%nProc), disp_group(this%nProc))
    CALL MPI_ALLGATHER(nCell_perProc, 1, MPI_INTEGER4, gcount_group, 1, MPI_INTEGER4, this%comm, this%IERR)

    gcount_group = gcount_group

    disp_group = 0

    DO k1 = 2, this%nProc
      disp_group(k1) = gcount_group(k1 - 1) + disp_group(k1 - 1)
    END DO

    CALL this%barrier
    CALL mpi_scatterv(bufTemp, gcount_group, disp_group, MPI_INTEGER4, &
                      buf_dist, gcount_group(this%myrank), MPI_INTEGER4, this%rankBase, this%comm, this%ierr)

    DEALLOCATE (gcount_group, disp_group)
    IF (this%isBaseProc()) DEALLOCATE (bufTemp)
  END SUBROUTINE DistGridIntSeq1D

  SUBROUTINE generate_recv_info(this, cell_stcl, num_iCell, num_cell_stcl, dimCell_global) ! ncp is limited at this point 4^p*20, p = 0, 1, ...
    IMPLICIT NONE

    CLASS(mpddInfo_sg_t) :: this
    ! INTEGER, intent(in) :: adjnb(:, :)
    INTEGER(i_kind), INTENT(in) :: cell_stcl(:, :), num_iCell, num_cell_stcl, dimCell_global(2)

    INTEGER :: size_shape(1)
    INTEGER :: i, j, k, k1
    INTEGER, ALLOCATABLE :: halo_mask(:)
    CHARACTER(len=2) :: gls
    INTEGER :: cur_blk, blk_num, jp1, blk_sz, blk_side_len, num_iCell_global

    this%halo_size = 0

    IF (.NOT. this%isActiveProc()) RETURN

    num_iCell_global = dimCell_global(1) * dimCell_global(2)

    blk_sz = num_iCell ! size of the block for a computer task
    blk_side_len = 2 * 2**this%gLevel + 2!ntmp/nproc_mg !Zilong amENDed here

    ! recv_info%my_start_idx = (this%myrank - 1)*blk_sz + 1
    ! recv_info%my_END_idx = this%myrank*blk_sz
    this%recv_info%num_icell = num_icell
    this%recv_info%nb_blk_sz = 0  !Zilong added

    ALLOCATE (this%recv_info%nb_blk(max_nb))
    ALLOCATE (halo_mask(num_iCell_global))

    halo_mask = 0
    !$omp PARALLEL DO PRIVATE(i,j)
    DO i = 1, blk_sz
      DO j = 1, num_cell_stcl
        ! blk_num = i2blk_gl(glvl, edge_stcl(k, j, i), nproc_mg)
        IF (cell_stcl(j, i) .EQ. 0) CYCLE
        blk_num = this%i_global2ProcRank(cell_stcl(j, i), dimCell_global)
        IF (blk_num == this%myrank) THEN
          CYCLE
        ELSE
          halo_mask(cell_stcl(j, i)) = 1
        END IF
      END DO
    END DO
    !$omp END PARALLEL DO

! #ifdef DIAG
!     PRINT *, 'Current block =', this%myrank
!     PRINT *, 'Block start and END =', recv_info%my_start_idx, recv_info%my_END_idx
! #endif

! ! packing halo into idx array for originating block
    DO j = 1, max_nb
      this%recv_info%nb_blk(j)%sz = 0
      this%recv_info%nb_blk(j)%blk_num = 0
      ALLOCATE (this%recv_info%nb_blk(j)%idx(blk_side_len))
      this%recv_info%nb_blk(j)%idx = 0
    END DO

    DO i = 1, num_iCell_global
      IF (halo_mask(i) == 0) THEN
        CYCLE
      ELSE
        blk_num = this%i_global2ProcRank(i, dimCell_global)
        DO j = 1, this%nProc
          IF (this%recv_info%nb_blk(j)%blk_num == blk_num) THEN
            this%recv_info%nb_blk(j)%sz = this%recv_info%nb_blk(j)%sz + 1
            this%recv_info%nb_blk(j)%idx(this%recv_info%nb_blk(j)%sz) = i
            !Dynamic memeory Zilong added
            size_shape = SHAPE(this%recv_info%nb_blk(j)%idx)
            IF (size_shape(1) == this%recv_info%nb_blk(j)%sz) THEN
              CALL re_ALLOCATE_idx_size_increase1D(this%recv_info%nb_blk(j)%idx, blk_side_len)
            END IF

            EXIT
          ELSE IF (this%recv_info%nb_blk(j)%blk_num == 0) THEN
            this%recv_info%nb_blk(j)%blk_num = blk_num
            this%recv_info%nb_blk(j)%sz = this%recv_info%nb_blk(j)%sz + 1
            this%recv_info%nb_blk_sz = this%recv_info%nb_blk_sz + 1  !Zilong added here
            !Dynamic memeory Zilong added
            size_shape = SHAPE(this%recv_info%nb_blk)
            IF (size_shape(1) == this%recv_info%nb_blk_sz) THEN
              CALL re_ALLOCATE_nbblk(this%recv_info%nb_blk, max_nb, blk_side_len)
            END IF

            this%recv_info%nb_blk(j)%idx(this%recv_info%nb_blk(j)%sz) = i
            EXIT
          END IF
        END DO
      END IF
    END DO

    this%recv_info%halo_sz = 0
    DO j = 1, this%recv_info%nb_blk_sz
      this%recv_info%halo_sz = this%recv_info%halo_sz + this%recv_info%nb_blk(j)%sz
    END DO

    this%halo_size = this%recv_info%halo_sz

    CALL mpi_barrier(this%comm, this%ierr)

    PRINT *, 'halo_size is: ', this%halo_size

! #ifdef DIAG
!     DO j = 1, recv_info%nb_blk_sz
!       PRINT *, 'task num = ', recv_info%nb_blk(j)%blk_num, 'task halo size = ', recv_info%nb_blk(j)%sz
!       PRINT *, recv_info%nb_blk(j)%idx
!     ENDDO
!     PRINT *, 'halo size = ', recv_info%halo_sz
!     PRINT *, '=================================================================='
! #endif
!     ! ENDDO

    DEALLOCATE (halo_mask)
  END SUBROUTINE generate_recv_info

  SUBROUTINE re_ALLOCATE_idx_size_increase1D(cache, isize)
    IMPLICIT NONE
    INTEGER, ALLOCATABLE, INTENT(INOUT):: cache(:)
    INTEGER, INTENT(IN) :: isize
    INTEGER :: size_shape(1)
    INTEGER :: osize
    INTEGER, ALLOCATABLE :: swap(:)

    size_shape = SHAPE(cache)
    osize = size_shape(1)
    ALLOCATE (swap(osize))
    swap(1:osize) = cache(1:osize)
    DEALLOCATE (cache)
    ALLOCATE (cache(osize + isize))
    cache(1:osize) = swap(1:osize)
    DEALLOCATE (swap)
  END SUBROUTINE re_ALLOCATE_idx_size_increase1D

  SUBROUTINE re_ALLOCATE_nbblk(nb_blk, isize, blk_side_len)
    IMPLICIT NONE
    TYPE(NBlock), ALLOCATABLE, INTENT(INOUT):: nb_blk(:)
    INTEGER, INTENT(IN) :: isize, blk_side_len
    INTEGER :: size_shape(1)
    INTEGER :: osize, k1
    TYPE(NBlock), ALLOCATABLE :: swap_nbblk(:)

    size_shape = SHAPE(nb_blk)
    osize = size_shape(1)

    ALLOCATE (swap_nbblk(osize))
    DO k1 = 1, osize
      ALLOCATE (swap_nbblk(k1)%idx(nb_blk(k1)%sz))
      swap_nbblk(k1) = nb_blk(k1)
      DEALLOCATE (nb_blk(k1)%idx)
    END DO

    DEALLOCATE (nb_blk)
    ALLOCATE (nb_blk(osize + isize))

    DO k1 = 1, osize
      ALLOCATE (nb_blk(k1)%idx(swap_nbblk(k1)%sz))
      nb_blk(k1) = swap_nbblk(k1)
      DEALLOCATE (swap_nbblk(k1)%idx)
    END DO
    DEALLOCATE (swap_nbblk)

    DO k1 = osize + 1, osize + max_nb
      nb_blk(k1)%sz = 0
      nb_blk(k1)%blk_num = 0
      ALLOCATE (nb_blk(k1)%idx(blk_side_len))
      nb_blk(k1)%idx = 0
    END DO
  END SUBROUTINE

  SUBROUTINE trim_idx_size_1D(cache, rsize)
    IMPLICIT NONE
    INTEGER, ALLOCATABLE, INTENT(INOUT):: cache(:)
    INTEGER, INTENT(IN) :: rsize
    INTEGER, ALLOCATABLE :: swap(:)

    ALLOCATE (swap(rsize))
    swap(1:rsize) = cache(1:rsize)
    DEALLOCATE (cache)
    ALLOCATE (cache(rsize))
    cache(1:rsize) = swap(1:rsize)
    DEALLOCATE (swap)
  END SUBROUTINE trim_idx_size_1D

  !! for cell_stcl
  INTEGER FUNCTION cal_num_cell_stcl(this, halo_width)
    CLASS(mpddInfo_sg_t) :: this
    INTEGER(i_kind), INTENT(IN) :: halo_width

    cal_num_cell_stcl = (1 + 2 * halo_width)**2

  END FUNCTION cal_num_cell_stcl

  SUBROUTINE generate_send_info(this)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) :: this

    ! INTEGER, INTENT(IN) :: mycomp_task
    ! TYPE(HaloBlock), INTENT(INOUT) :: send_info
    ! TYPE(HaloBlock), INTENT(IN) :: recv_info
    TYPE(HaloBlock) :: swap_info
    ! INTEGER, INTENT(IN) :: nproc_mg

    INTEGER :: i, j, k, blk, sinfo_size
    INTEGER :: length_halo(1)
    INTEGER :: size_shape(1)
    INTEGER(i_kind), ALLOCATABLE :: req2(:)
    INTEGER(i_kind), ALLOCATABLE :: req1(:)

    IF (.NOT. this%isActiveProc()) RETURN

    ALLOCATE (req1(this%nProc), req2(this%nProc))

    req1 = 0
    ALLOCATE (swap_info%nb_blk(this%nProc))

    DO i = 1, this%nProc
      CALL mpi_irecv(swap_info%nb_blk(i)%sz, 1, MPI_INTEGER4, i - 1, 99, this%comm, req2(i), this%ierr)
      swap_info%nb_blk(i)%blk_num = i
    END DO

    DO j = 1, this%nProc
      DO i = 1, this%recv_info%nb_blk_sz
        IF (j == this%recv_info%nb_blk(i)%blk_num) THEN
          CALL MPI_ssend(this%recv_info%nb_blk(i)%sz, 1, MPI_INTEGER4, &
                         this%recv_info%nb_blk(i)%blk_num - 1, 99, this%comm, this%ierr)
          req1(j) = 1
        END IF
      END DO
      IF (req1(j) /= 1) THEN
        CALL MPI_ssend(0, 1, MPI_INTEGER4, j - 1, 99, this%comm, this%ierr)
      END IF
    END DO

    DO i = 1, this%nProc
      CALL MPI_wait(req2(i), this%Status, this%ierr)
    END DO

    CALL mpi_barrier(this%comm, this%ierr)
    PRINT *, 'End of sending halo szie.', this%myrank
    ! CALL mpi_barrier(this%comm, this%ierr)

    sinfo_size = 0
    DO i = 1, this%nProc
      IF (swap_info%nb_blk(i)%sz /= 0) THEN
        sinfo_size = sinfo_size + 1
      END IF
    END DO

    ALLOCATE (this%send_info%nb_blk(sinfo_size))

    k = 1
    this%send_info%halo_sz = 0
    this%send_info%nb_blk_sz = sinfo_size
    DO i = 1, this%nProc
      IF (swap_info%nb_blk(i)%sz /= 0) THEN
        this%send_info%nb_blk(k)%blk_num = swap_info%nb_blk(i)%blk_num
        this%send_info%nb_blk(k)%sz = swap_info%nb_blk(i)%sz
        ALLOCATE (this%send_info%nb_blk(k)%idx(swap_info%nb_blk(i)%sz))  ! this is a fake ALLOCATE
        this%send_info%nb_blk(k)%idx = 0
        this%send_info%halo_sz = this%send_info%halo_sz + this%send_info%nb_blk(k)%sz
        k = k + 1
      END IF
    END DO

    DO i = 1, sinfo_size
      CALL mpi_irecv(this%send_info%nb_blk(i)%idx, this%send_info%nb_blk(i)%sz, MPI_INTEGER4, &
                     this%send_info%nb_blk(i)%blk_num - 1, 97, this%comm, req2(i), this%ierr)
    END DO

    DO i = 1, this%recv_info%nb_blk_sz
      CALL MPI_ssend(this%recv_info%nb_blk(i)%idx, this%recv_info%nb_blk(i)%sz, MPI_INTEGER4, &
                     this%recv_info%nb_blk(i)%blk_num - 1, 97, this%comm, this%ierr)
    END DO

    DO i = 1, sinfo_size
      CALL MPI_wait(req2(i), this%Status, this%ierr)
    END DO

    CALL mpi_barrier(this%comm, this%ierr)
    PRINT *, 'End of sending halo swap info.', this%myrank
    ! CALL mpi_barrier(this%comm, this%ierr)

    DEALLOCATE (swap_info%nb_blk)

    DEALLOCATE (req1)
    DEALLOCATE (req2)
  END SUBROUTINE generate_send_info

  SUBROUTINE trim_idx_all_halos(this)
    IMPLICIT NONE
    CLASS(mpddInfo_sg_t) :: this
    INTEGER :: length_nb(1)
    INTEGER :: i

    length_nb = SHAPE(this%recv_info%nb_blk)
    DO i = 1, length_nb(1)
      CALL trim_idx_size_1D(this%recv_info%nb_blk(i)%idx, this%recv_info%nb_blk(i)%sz)
    END DO
    length_nb = SHAPE(this%send_info%nb_blk)
    DO i = 1, length_nb(1)
      CALL trim_idx_size_1D(this%send_info%nb_blk(i)%idx, this%send_info%nb_blk(i)%sz)
    END DO
  END SUBROUTINE trim_idx_all_halos

END MODULE mpddInfo_sg_m
