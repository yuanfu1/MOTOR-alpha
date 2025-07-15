MODULE InterpHPNew_m
  USE slint, ONLY: slint_init, bilinear_interp, nn_interp, tgt_grid
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE mpddSub_m, ONLY: mpddSub_t
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE parameters_m, ONLY: degree2radian

  INCLUDE "mpif.h"
  TYPE, EXTENDS(mpddSub_t) :: InterpHPNew_t
    TYPE(mpddGlob_t), POINTER :: mpdd
    INTEGER(i_kind), ALLOCATABLE :: ncount_group(:), disp_group(:)
    INTEGER(i_kind), ALLOCATABLE :: idxArray(:)
    INTEGER(i_kind) :: nTgtEach, nSrcAll

  CONTAINS
    PROCEDURE, PUBLIC :: interp_multiLevel_d
    PROCEDURE, PUBLIC :: interp_multiLevel_s
    PROCEDURE, PUBLIC :: interp_multiLevel_nearest_s
    PROCEDURE, PUBLIC :: interp_multiLevel_nearest_d

    PROCEDURE, PUBLIC :: interp_singleLevel_d
    PROCEDURE, PUBLIC :: interp_singleLevel_s
    PROCEDURE, PUBLIC :: interp_singleLevel_nearest_s
    PROCEDURE, PUBLIC :: interp_singleLevel_nearest_d

    GENERIC, PUBLIC :: interpN => interp_multiLevel_d, &
      interp_multiLevel_s, interp_singleLevel_d, interp_singleLevel_s
    FINAL :: destructor

    PROCEDURE, PUBLIC :: initializeInterpHPNew
  END TYPE InterpHPNew_t

CONTAINS

  SUBROUTINE initializeInterpHPNew(this, mpdd, sg, llsrc, nsrc)
    CLASS(InterpHPNew_t) :: this
    TYPE(mpddGlob_t), TARGET, INTENT(IN) :: mpdd
    TYPE(SingleGrid_t), INTENT(IN) :: sg

    INTEGER(i_kind), INTENT(IN) :: nsrc
    REAL(r_kind) :: llsrc(2, nsrc)
    INTEGER(i_kind) :: nt = 0, i
    REAL(r_kind), ALLOCATABLE :: lltgtEach(:, :)
    REAL(r_kind), ALLOCATABLE :: llsrcEach(:, :)
    LOGICAL, ALLOCATABLE :: idxMask(:)
    INTEGER(i_kind), ALLOCATABLE :: idxSelected(:)
    INTEGER(i_kind) :: nStcl, nsrcAllSwap
    REAL(r_kind), ALLOCATABLE :: llsrcAllSwap(:, :)
    REAL(r_kind) :: lltgtAllSwap(2, sg%num_cell), cnrtPtTemp(2)

    ! Yuanfu Xie added safeguard check the interpolation: 2024-07-17:
    INTEGER(i_kind) :: imn, j
    REAL(r_kind) :: amn

    this%mpdd => mpdd
    CALL this%initialize(sg%mpddInfo_sg%nproc, mpdd%group, mpdd%rank)

    ! CALL slint_init 获取 sg%cell_cntr 所对应llsrc中的插值格点，索引返回在tgt_grid%nn中
    IF (this%mpddSub_t%isActiveProc()) THEN

      IF (this%mpddSub_t%isBaseProc()) nsrcAllSwap = nsrc
      CALL this%bcast(nsrcAllSwap)
      ALLOCATE (llsrcAllSwap(2, nsrcAllSwap))
      IF (this%mpddSub_t%isBaseProc()) llsrcAllSwap = llsrc
      CALL this%bcast(llsrcAllSwap)

      ! Set all boundary point to be a known point
      PRINT *, "Before beginning of slint_init on ", this%myrank, nsrcAllSwap, llsrcAllSwap(:, 1)
      lltgtAllSwap = sg%cell_cntr

      cnrtPtTemp(1) = MAXVAL(lltgtAllSwap(1, :)) / 2 + MINVAL(lltgtAllSwap(1, :)) / 2
      cnrtPtTemp(2) = MAXVAL(lltgtAllSwap(2, :)) / 2 + MINVAL(lltgtAllSwap(2, :)) / 2
      PRINT *, 'cnrtPtTemp: ', cnrtPtTemp / degree2radian

      FORALL (i=1:sg%num_cell, sg%cell_type(i) .GT. 0) lltgtAllSwap(:, i) = cnrtPtTemp

      ! Calculate the stencil
      ! print*, "(MAXVAL(lltgtAllSwap(1,:))): ", (MAXVAL( lltgtAllSwap(1,:)))/degree2radian, sg%mpddInfo_sg%myrank, sg%gLevel
      ! print*, "(MAXVAL(lltgtAllSwap(2,:))): ", (MAXVAL( lltgtAllSwap(2,:)))/degree2radian, sg%mpddInfo_sg%myrank, sg%gLevel
      ! print*, "(MINVAL(lltgtAllSwap(1,:))): ", (MINVAL( lltgtAllSwap(1,:)))/degree2radian, sg%mpddInfo_sg%myrank, sg%gLevel
      ! print*, "(MINVAL(lltgtAllSwap(2,:))): ", (MINVAL( lltgtAllSwap(2,:)))/degree2radian, sg%mpddInfo_sg%myrank, sg%gLevel

      CALL sg%mpddInfo_sg%barrier
      ! Calculate the stencil
      ! print*, "(MAXVAL(llsrcAllSwap(1,:))): ", (MAXVAL( llsrcAllSwap(1,:)))/degree2radian, sg%mpddInfo_sg%myrank, sg%gLevel
      ! print*, "(MAXVAL(llsrcAllSwap(2,:))): ", (MAXVAL( llsrcAllSwap(2,:)))/degree2radian, sg%mpddInfo_sg%myrank, sg%gLevel
      ! print*, "(MINVAL(llsrcAllSwap(1,:))): ", (MINVAL( llsrcAllSwap(1,:)))/degree2radian, sg%mpddInfo_sg%myrank, sg%gLevel
      ! print*, "(MINVAL(llsrcAllSwap(2,:))): ", (MINVAL( llsrcAllSwap(2,:)))/degree2radian, sg%mpddInfo_sg%myrank, sg%gLevel

      CALL slint_init(TRANSPOSE(llsrcAllSwap), nsrcAllSwap, TRANSPOSE(lltgtAllSwap), sg%num_cell)
      PRINT *, '|-> SHAPE of sg%cell_cntr: ', SHAPE(lltgtAllSwap), 'rank :', this%myrank

      ! Check the interpolation coefficients:
      ! Yuanfu Xie added safeguard check the interpolation: 2024-07-17:
      imn = 0
      amn = 1.0D8
      DO j = 1, UBOUND(tgt_grid%coeffs, 2)
        IF (SUM(tgt_grid%coeffs(:, j)) .LT. amn) THEN
          imn = j
          amn = SUM(tgt_grid%coeffs(:, j))
        END IF
      END DO
      IF (amn .LT. 0.9D0) THEN
        WRITE (*, 3)
        WRITE (*, 1) tgt_grid%coeffs(:, imn)
1       FORMAT('| The source grid cannot cover the target grid with min coef:    |', /, &
               '|   ', 3E12.4, '                         |', /, &
               '| Please check and rerun!                                        |')
        WRITE (*, 2) MINVAL(llsrcAllSwap(1, :)) / degree2radian, MAXVAL(llsrcAllSwap(1, :)) / degree2radian, &
          MINVAL(llsrcAllSwap(2, :)) / degree2radian, MAXVAL(llsrcAllSwap(2, :)) / degree2radian, &
          MINVAL(lltgtAllSwap(1, :)) / degree2radian, MAXVAL(lltgtAllSwap(1, :)) / degree2radian, &
          MINVAL(lltgtAllSwap(2, :)) / degree2radian, MAXVAL(lltgtAllSwap(2, :)) / degree2radian
2       FORMAT('| Source LL: ', 4E12.4, '    |', /, '| Target LL: ', 4E12.4, '    |')
        WRITE (*, 3)
3       FORMAT('+----------------------------------------------------------------+')
        STOP
      END IF

      DEALLOCATE (llsrcAllSwap)

      nStcl = 3 !MOTOR-DA中，一个格点，所需llsrc插值点的个数
      ALLOCATE (idxMask(nsrcAllSwap)); 
      idxMask = .FALSE.

      idxMask(RESHAPE(tgt_grid%nn, (/sg%num_cell * nStcl/))) = .TRUE.!将所需的网格点标记为真
      idxSelected = PACK((/(i, i=1, nsrcAllSwap)/), idxMask) !将所需的格点索引，按次序排成一列
      PRINT *, this%myrank, MINVAL(tgt_grid%nn), SHAPE(tgt_grid%nn), "SIZE(idxSelected): ", SIZE(idxSelected)

      ALLOCATE (this%ncount_group(this%nproc))
      ALLOCATE (this%disp_group(this%nproc))

      ! 将所需的插值格点的个数，汇总到主线程 0
      CALL MPI_GATHER(SIZE(idxSelected), 1, MPI_INTEGER4, &
                      this%ncount_group, 1, MPI_INTEGER4, 0, &
                      this%comm, this%ierr)

      IF (this%isBaseProc()) THEN
        ! PRINT *, "ncount_group is ", this%ncount_group
        ALLOCATE (this%idxArray(SUM(this%ncount_group)))
        this%disp_group = 0
        FORALL (i=2:this%nproc) this%disp_group(i) = SUM(this%ncount_group(1:i - 1)) !Gatherv 时，每个进程所需的位移
        PRINT *, this%myrank, this%disp_group, this%ncount_group
      END IF

      CALL this%bcast(this%ncount_group)
      CALL this%bcast(this%disp_group)

      ! 将所有进程的，所需插值格点的索引汇总，存在this%idxArray中
      CALL MPI_GATHERV(idxSelected, SIZE(idxSelected), MPI_INTEGER4, &
                       this%idxArray, this%ncount_group, this%disp_group, MPI_INTEGER4, 0, &
                       this%comm, this%ierr)

      BLOCK
        REAL(r_kind) :: llsrcEach(2, SIZE(idxSelected))
        REAL(r_kind), ALLOCATABLE :: llsrcAll(:, :)

        ! 按照汇总的索引，得到汇总的索引对应的，经纬度的值
        IF (this%isBaseProc()) THEN
          ALLOCATE (llsrcAll(2, SUM(this%ncount_group)))
          FORALL (i=1:SUM(this%ncount_group)) llsrcAll(:, i) = llsrc(:, this%idxArray(i))
        END IF

        ! 将经纬度分发给各个线程
        CALL MPI_SCATTERV(llsrcAll, this%ncount_group * 2, this%disp_group * 2, MPI_DOUBLE_PRECISION, &
                          llsrcEach, SIZE(idxSelected) * 2, MPI_DOUBLE_PRECISION, 0, &
                          this%comm, this%ierr)

        ! IF (this%isBaseProc()) THEN
        !   PRINT *, llsrcAll(:, this%disp_group + 1)
        ! END IF
        ! CALL this%barrier
        ! PRINT *, this%myrank, llsrcEach(:, 1)

        ! 对每个线程，按照分发的格点，重新插值
        CALL slint_init(TRANSPOSE(llsrcEach), SIZE(idxSelected), TRANSPOSE(lltgtAllSwap), sg%num_cell)
        this%nTgtEach = sg%num_cell
        this%nSrcAll = nsrcAllSwap
      END BLOCK
    END IF
    ! CALL this%mpdd%barrier
  END SUBROUTINE

  SUBROUTINE interp_multiLevel_d(this, vLevel, datasrc, datatgt)
    CLASS(InterpHPNew_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vLevel
    REAL(r_kind), INTENT(IN) :: datasrc(:, :)
    REAL(r_kind), INTENT(INOUT) :: datatgt(vLevel, this%nTgtEach)
    REAL(r_kind), ALLOCATABLE :: datatgtEach(:, :)
    INTEGER(i_kind) :: i, disp_group(this%nproc)
    REAL(r_kind) :: dataSrcEach(vLevel, this%ncount_group(this%myrank))
    REAL(r_kind), ALLOCATABLE :: dataSrcAll(:, :)

    IF (this%isActiveProc()) THEN
      PRINT *, 'ENTERING interp_multiLevel_d...'

      IF (this%isBaseProc()) THEN
        ALLOCATE (dataSrcAll(vLevel, SUM(this%ncount_group)))
        FORALL (i=1:SUM(this%ncount_group)) dataSrcAll(:, i) = datasrc(:, this%idxArray(i))! 数组重排
      END IF

      ! PRINT *, this%myrank, "this%ncount_group*vLevel: ", this%ncount_group*vLevel
      ! PRINT *, this%myrank, "this%disp_group*vLevel: ", this%disp_group*vLevel

      ! 分发dataSrcAll
      CALL MPI_SCATTERV(dataSrcAll, this%ncount_group * vLevel, this%disp_group * vLevel, MPI_DOUBLE_PRECISION, &
                        dataSrcEach, this%ncount_group(this%myrank) * vLevel, MPI_DOUBLE_PRECISION, 0, &
                        this%comm, this%ierr)

      CALL bilinear_interp(dataSrcEach, datatgt)

      IF (ALLOCATED(dataSrcAll)) DEALLOCATE (dataSrcAll)
      CALL this%barrier
    END IF
  END SUBROUTINE interp_multiLevel_d

  SUBROUTINE interp_multiLevel_s(this, vLevel, datasrc, datatgt)
    CLASS(InterpHPNew_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vLevel
    REAL(r_single), INTENT(IN) :: datasrc(:, :)
    REAL(r_kind), INTENT(INOUT) :: datatgt(vLevel, this%nTgtEach)
    REAL(r_single), ALLOCATABLE :: datatgtEach(:, :)
    INTEGER(i_kind) :: i, disp_group(this%nproc)
    REAL(r_single) :: dataSrcEach(vLevel, this%ncount_group(this%myrank))
    REAL(r_single), ALLOCATABLE :: dataSrcAll(:, :)
    REAL(r_single) :: datatgt_s(vLevel, this%nTgtEach)

    datatgt_s = 0.0D0

    IF (this%isActiveProc()) THEN
      PRINT *, 'ENTERING interp_multiLevel_s...'

      IF (this%isBaseProc()) THEN
        ALLOCATE (dataSrcAll(vLevel, SUM(this%ncount_group)))
        FORALL (i=1:SUM(this%ncount_group)) dataSrcAll(:, i) = datasrc(:, this%idxArray(i))! 数组重排
      END IF

      ! 分发dataSrcAll
      CALL MPI_SCATTERV(dataSrcAll, this%ncount_group * vLevel, this%disp_group * vLevel, MPI_REAL4, &
                        dataSrcEach, this%ncount_group(this%myrank) * vLevel, MPI_REAL4, 0, &
                        this%comm, this%ierr)

      CALL bilinear_interp(dataSrcEach, datatgt_s)

      datatgt = datatgt_s
      IF (ALLOCATED(dataSrcAll)) DEALLOCATE (dataSrcAll)
      CALL this%barrier
    END IF
  END SUBROUTINE interp_multiLevel_s

  SUBROUTINE interp_multiLevel_nearest_d(this, vLevel, datasrc, datatgt)
    CLASS(InterpHPNew_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vLevel
    REAL(r_kind), INTENT(IN) :: datasrc(:, :)
    REAL(r_kind), INTENT(INOUT) :: datatgt(vLevel, this%nTgtEach)
    REAL(r_kind), ALLOCATABLE :: datatgtEach(:, :)
    INTEGER(i_kind) :: i, disp_group(this%nproc)
    REAL(r_kind) :: dataSrcEach(vLevel, this%ncount_group(this%myrank))
    REAL(r_kind), ALLOCATABLE :: dataSrcAll(:, :)

    IF (this%isActiveProc()) THEN
      IF (this%isBaseProc()) THEN
        ALLOCATE (dataSrcAll(vLevel, SUM(this%ncount_group)))
        FORALL (i=1:SUM(this%ncount_group)) dataSrcAll(:, i) = datasrc(:, this%idxArray(i))! 数组重排
      END IF

      ! PRINT *, this%myrank, "this%ncount_group*vLevel: ", this%ncount_group*vLevel
      ! PRINT *, this%myrank, "this%disp_group*vLevel: ", this%disp_group*vLevel

      ! 分发dataSrcAll
      CALL MPI_SCATTERV(dataSrcAll, this%ncount_group * vLevel, this%disp_group * vLevel, MPI_DOUBLE_PRECISION, &
                        dataSrcEach, this%ncount_group(this%myrank) * vLevel, MPI_DOUBLE_PRECISION, 0, &
                        this%comm, this%ierr)

      FORALL (i=1:this%nTgtEach) datatgt(:, i) = dataSrcEach(:, tgt_grid%nn(1, i))

      IF (ALLOCATED(dataSrcAll)) DEALLOCATE (dataSrcAll)
      CALL this%barrier
    END IF
  END SUBROUTINE

  SUBROUTINE interp_multiLevel_nearest_s(this, vLevel, datasrc, datatgt)
    CLASS(InterpHPNew_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vLevel
    REAL(r_single), INTENT(IN) :: datasrc(:, :)
    REAL(r_kind), INTENT(INOUT) :: datatgt(vLevel, this%nTgtEach)
    REAL(r_single), ALLOCATABLE :: datatgtEach(:, :)
    INTEGER(i_kind) :: i, disp_group(this%nproc)
    REAL(r_single) :: dataSrcEach(vLevel, this%ncount_group(this%myrank))
    REAL(r_single), ALLOCATABLE :: dataSrcAll(:, :)
    REAL(r_single) :: datatgt_s(vLevel, this%nTgtEach)

    IF (this%isActiveProc()) THEN
      IF (this%isBaseProc()) THEN
        ALLOCATE (dataSrcAll(vLevel, SUM(this%ncount_group)))
        FORALL (i=1:SUM(this%ncount_group)) dataSrcAll(:, i) = datasrc(:, this%idxArray(i))! 数组重排
      END IF

      ! PRINT *, this%myrank, "this%ncount_group*vLevel: ", this%ncount_group*vLevel
      ! PRINT *, this%myrank, "this%disp_group*vLevel: ", this%disp_group*vLevel

      ! 分发dataSrcAll
      CALL MPI_SCATTERV(dataSrcAll, this%ncount_group * vLevel, this%disp_group * vLevel, MPI_REAL4, &
                        dataSrcEach, this%ncount_group(this%myrank) * vLevel, MPI_REAL4, 0, &
                        this%comm, this%ierr)

      FORALL (i=1:this%nTgtEach) datatgt_s(:, i) = dataSrcEach(:, tgt_grid%nn(1, i))

      datatgt = datatgt_s
      IF (ALLOCATED(dataSrcAll)) DEALLOCATE (dataSrcAll)
      CALL this%barrier
    END IF
  END SUBROUTINE

  SUBROUTINE interp_singleLevel_d(this, datasrc, datatgt)
    CLASS(InterpHPNew_t) :: this
    REAL(r_kind), INTENT(IN) :: datasrc(:)
    REAL(r_kind), INTENT(INOUT) :: datatgt(this%nTgtEach)
    REAL(r_kind), ALLOCATABLE :: datatgtEach(:)
    INTEGER(i_kind) :: i, disp_group(this%nproc)
    REAL(r_kind) :: dataSrcEach(this%ncount_group(this%myrank))
    REAL(r_kind), ALLOCATABLE :: dataSrcAll(:)

    PRINT *, 'ENTERING interp_singleLevel_d...'
    IF (this%isActiveProc()) THEN
      IF (this%isBaseProc()) THEN
        ALLOCATE (dataSrcAll(SUM(this%ncount_group)))
        FORALL (i=1:SUM(this%ncount_group)) dataSrcAll(i) = datasrc(this%idxArray(i))! 数组重排
      END IF

      ! PRINT *, this%myrank, "this%ncount_group: ", this%ncount_group
      ! PRINT *, this%myrank, "this%disp_group: ", this%disp_group

      ! 分发dataSrcAll
      CALL MPI_SCATTERV(dataSrcAll, this%ncount_group, this%disp_group, MPI_DOUBLE_PRECISION, &
                        dataSrcEach, this%ncount_group(this%myrank), MPI_DOUBLE_PRECISION, 0, &
                        this%comm, this%ierr)

      CALL bilinear_interp(dataSrcEach, datatgt)

      IF (ALLOCATED(dataSrcAll)) DEALLOCATE (dataSrcAll)
      CALL this%barrier
    END IF
  END SUBROUTINE

  SUBROUTINE interp_singleLevel_s(this, datasrc, datatgt)
    CLASS(InterpHPNew_t) :: this
    REAL(r_single), INTENT(IN) :: datasrc(:)
    REAL(r_kind), INTENT(INOUT) :: datatgt(this%nTgtEach)
    REAL(r_single), ALLOCATABLE :: datatgtEach(:)
    INTEGER(i_kind) :: i, disp_group(this%nproc)
    REAL(r_single) :: dataSrcEach(this%ncount_group(this%myrank))
    REAL(r_single), ALLOCATABLE :: dataSrcAll(:)
    REAL(r_single) :: datatgt_s(this%nTgtEach)

    IF (this%isActiveProc()) THEN
      PRINT *, 'ENTERING interp_singleLevel_s...'

      IF (this%isBaseProc()) THEN
        ALLOCATE (dataSrcAll(SUM(this%ncount_group)))
        FORALL (i=1:SUM(this%ncount_group)) dataSrcAll(i) = datasrc(this%idxArray(i))! 数组重排
      END IF

      ! PRINT *, this%myrank, "this%ncount_group: ", this%ncount_group
      ! PRINT *, this%myrank, "this%disp_group: ", this%disp_group

      ! 分发dataSrcAll
      CALL MPI_SCATTERV(dataSrcAll, this%ncount_group, this%disp_group, MPI_REAL4, &
                        dataSrcEach, this%ncount_group(this%myrank), MPI_REAL4, 0, &
                        this%comm, this%ierr)

      CALL bilinear_interp(dataSrcEach, datatgt_s)

      datatgt = datatgt_s
      IF (ALLOCATED(dataSrcAll)) DEALLOCATE (dataSrcAll)
      CALL this%barrier
    END IF
  END SUBROUTINE

  SUBROUTINE interp_singleLevel_nearest_d(this, datasrc, datatgt)
    CLASS(InterpHPNew_t) :: this
    REAL(r_kind), INTENT(IN) :: datasrc(:)
    REAL(r_kind), INTENT(INOUT) :: datatgt(this%nTgtEach)
    REAL(r_kind), ALLOCATABLE :: datatgtEach(:)
    INTEGER(i_kind) :: i, disp_group(this%nproc)
    REAL(r_kind) :: dataSrcEach(this%ncount_group(this%myrank))
    REAL(r_kind), ALLOCATABLE :: dataSrcAll(:)

    IF (this%isActiveProc()) THEN
      IF (this%isBaseProc()) THEN
        ALLOCATE (dataSrcAll(SUM(this%ncount_group)))
        FORALL (i=1:SUM(this%ncount_group)) dataSrcAll(i) = datasrc(this%idxArray(i))! 数组重排
      END IF

      ! PRINT *, this%myrank, "this%ncount_group: ", this%ncount_group
      ! PRINT *, this%myrank, "this%disp_group: ", this%disp_group

      ! 分发dataSrcAll
      CALL MPI_SCATTERV(dataSrcAll, this%ncount_group, this%disp_group, MPI_DOUBLE_PRECISION, &
                        dataSrcEach, this%ncount_group(this%myrank), MPI_DOUBLE_PRECISION, 0, &
                        this%comm, this%ierr)

      FORALL (i=1:this%nTgtEach) datatgt(i) = dataSrcEach(tgt_grid%nn(1, i))

      IF (ALLOCATED(dataSrcAll)) DEALLOCATE (dataSrcAll)
      CALL this%barrier
    END IF
  END SUBROUTINE

  SUBROUTINE interp_singleLevel_nearest_s(this, datasrc, datatgt)
    CLASS(InterpHPNew_t) :: this
    REAL(r_single), INTENT(IN) :: datasrc(:)
    REAL(r_kind), INTENT(INOUT) :: datatgt(this%nTgtEach)
    REAL(r_single), ALLOCATABLE :: datatgtEach(:)
    INTEGER(i_kind) :: i, disp_group(this%nproc)
    REAL(r_single) :: dataSrcEach(this%ncount_group(this%myrank))
    REAL(r_single), ALLOCATABLE :: dataSrcAll(:)
    REAL(r_single) :: datatgt_s(this%nTgtEach)

    IF (this%isActiveProc()) THEN
      IF (this%isBaseProc()) THEN
        ALLOCATE (dataSrcAll(SUM(this%ncount_group)))
        FORALL (i=1:SUM(this%ncount_group)) dataSrcAll(i) = datasrc(this%idxArray(i))! 数组重排
      END IF

      ! PRINT *, this%myrank, "this%ncount_group: ", this%ncount_group
      ! PRINT *, this%myrank, "this%disp_group: ", this%disp_group

      ! 分发dataSrcAll
      CALL MPI_SCATTERV(dataSrcAll, this%ncount_group, this%disp_group, MPI_REAL4, &
                        dataSrcEach, this%ncount_group(this%myrank), MPI_REAL4, 0, &
                        this%comm, this%ierr)

      FORALL (i=1:this%nTgtEach) datatgt_s(i) = dataSrcEach(tgt_grid%nn(1, i))

      datatgt = datatgt_s
      IF (ALLOCATED(dataSrcAll)) DEALLOCATE (dataSrcAll)
      CALL this%barrier
    END IF
  END SUBROUTINE

  SUBROUTINE destructor(this)
    TYPE(InterpHPNew_t) :: this
    IF (ALLOCATED(this%ncount_group)) DEALLOCATE (this%ncount_group)
    IF (ALLOCATED(this%disp_group)) DEALLOCATE (this%disp_group)
    IF (ALLOCATED(this%idxArray)) DEALLOCATE (this%idxArray)

    PRINT *, "Finish the destructor of InterpHPNew_t"
  END SUBROUTINE destructor

END MODULE InterpHPNew_m
