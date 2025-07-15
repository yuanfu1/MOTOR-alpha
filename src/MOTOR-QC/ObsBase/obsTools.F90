! A library for observation tools
! Created by TingShu @2023-01-31
! Subroutines:
! griddedObsMap
! ObsMinusState
MODULE obsTools_m
  USE ObsBase_m
  USE kinds_m, ONLY: i_kind
  USE ObsConvention_m2, ONLY: ObsConvention_t2
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE YAMLRead_m
  USE mpObs_m, ONLY: mpObs_t

  IMPLICIT NONE
CONTAINS

  !> @brief
  !================================================================
  !! This routine checks the obsSet Y to fill in background in area
  !! where no observation associated with the background field.
  SUBROUTINE bkgdFillToObs(S, Y, Z, iRange_ramp)
    IMPLICIT NONE
    TYPE(ObsSet_t), INTENT(IN) :: Y
    TYPE(ObsSet_t), INTENT(INOUT) :: Z
    TYPE(state_t) :: S !< The state should be at model space, rather than control space

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, iu, iv, ix, nin, nv, nh, nt, imax, imin, jmax, jmin, kmax, kmin, nfill, idd
    REAL(i_kind) :: a
    TYPE(State_t) :: X

    INTEGER(i_kind), ALLOCATABLE :: idc(:, :)
    REAL(r_kind), ALLOCATABLE :: O(:)
    INTEGER(i_kind) :: iRange_ramp(:, :)

    X = S%zeroCopy()

    ! Check the rampRange:
    WRITE (*, 1) X%sg%gLevel, X%sg%vLevel, iRange_ramp(1, :), X%sg%mpddInfo_sg%myrank
    WRITE (*, 2) X%sg%gLevel, X%sg%num_cell, iRange_ramp(2, :), X%sg%mpddInfo_sg%myrank
    WRITE (*, 3) X%sg%gLevel, X%sg%tSlots, iRange_ramp(3, :), X%sg%mpddInfo_sg%myrank
1   FORMAT('BkgcFillin - G', I2, ' Vertical  levels:', I6, ' ramp: ', 20I6)
2   FORMAT('BkgcFillin - G', I2, ' Horizon num_cell:', I6, ' ramp: ', 20I6)
3   FORMAT('BkgcFillin - G', I2, ' Time Framework slots:', I6, ' ramp: ', 20I6)

    IF (MINVAL(iRange_ramp(1, :)) .GT. X%sg%vLevel .AND. &
        MINVAL(iRange_ramp(2, :)) .GT. X%sg%num_cell .AND. &
        MINVAL(iRange_ramp(3, :)) .GT. X%sg%tSlots) THEN
      Z = Y
      WRITE (*, 4) X%sg%gLevel, X%sg%mpddInfo_sg%myrank
4     FORMAT('No background filled in Glevel:', I3, ' pc', I2)
      RETURN ! No fill in
    END IF

    ! distributing Y obs to X:
    DO i = 1, UBOUND(Y%ObsFields, 1)
      ! Associated obs:
      IF (TRIM(Y%ObsFields(i)%Get_Name()) .EQ. 'cdir' .OR. &
          TRIM(Y%ObsFields(i)%Get_Name()) .EQ. 'sdir' .OR. &
          TRIM(Y%ObsFields(i)%Get_Name()) .EQ. 'wspd' .OR. &
          TRIM(Y%ObsFields(i)%Get_Name()) .EQ. 'vel') THEN
        iu = X%getVarIdx('uwnd')
        iv = X%getVarIdx('vwnd')
        DO j = 1, UBOUND(Y%ObsFields(i)%idx, 1)
          X%fields(iu)%DATA(Y%ObsFields(i)%idx(j)%vIdx, &
                            Y%ObsFields(i)%idx(j)%hIdx, &
                            Y%ObsFields(i)%idx(j)%tIdx) = 1.0D0
          X%fields(iv)%DATA(Y%ObsFields(i)%idx(j)%vIdx, &
                            Y%ObsFields(i)%idx(j)%hIdx, &
                            Y%ObsFields(i)%idx(j)%tIdx) = 1.0D0
        END DO
      ELSE
        ix = X%getVarIdx(TRIM(Y%ObsFields(i)%Get_Name()))
        IF (ix .GT. 0) THEN ! Found model state var
          DO j = 1, UBOUND(Y%ObsFields(i)%idx, 1)
            X%fields(ix)%DATA(Y%ObsFields(i)%idx(j)%vIdx, &
                              Y%ObsFields(i)%idx(j)%hIdx, &
                              Y%ObsFields(i)%idx(j)%tIdx) = 1.0D0
          END DO
        END IF
      END IF
    END DO

    ! Exchange the halo values to make sure each process see all possible obs:
    CALL X%exHalo()

    ! Spread the markers to fill the ramp regions of the horizontal grid:
    DO iv = 1, UBOUND(X%fields, 1)
      DO j = 1, iRange_ramp(2, iv)
        ! These markers values may over-count because of the neighbors overlaps
        CALL X%fields(iv)%Sum_Neighbors()
        ! After each ramp spread of the obs, the halo values need updated:
        ! CALL X%exhalo()
        CALL X%sg%ExchangeMatOnHaloForFieldGrid(X%sg%tSlots, X%sg%vLevel, &
                                                X%Fields(iv)%DATA)
      END DO
    END DO

    ! Fill in:
    nin = UBOUND(Y%obsFields, 1)
    Z = Y   ! set Z to Y as default and add background to the existing states

    ! Allocate O, and idc: setting them to the maximum possible fillins:
    ALLOCATE (O(UBOUND(X%fields(1)%DATA, 1) * UBOUND(X%fields(1)%DATA, 2) * UBOUND(X%fields(1)%DATA, 3)), &
              idc(3, UBOUND(X%fields(1)%DATA, 1) * UBOUND(X%fields(1)%DATA, 2) * UBOUND(X%fields(1)%DATA, 3)))

    idd = 0
    DO iv = 1, UBOUND(X%fields, 1)
      nfill = 0
      O = 0.0D0
      idc = 0
      nv = UBOUND(X%fields(iv)%DATA, 1)
      nh = X%sg%num_icell ! Only on the interior points to consider fillin
      nt = UBOUND(X%fields(iv)%DATA, 3)
      DO i = 1, 1 ! nv
        imin = 1 ! MAX(1, i - iRange_ramp(1, iv))
        imax = 1 !MIN(nv, i + iRange_ramp(1, iv))
        DO j = 1, nh
          jmin = MAX(1, j - iRange_ramp(2, iv))
          jmax = MIN(nh, j + iRange_ramp(2, iv))
          DO k = 1, nt
            kmin = MAX(1, k - iRange_ramp(3, iv))
            kmax = MIN(nt, k + iRange_ramp(3, iv))
            A = SUM(X%fields(iv)%DATA(imin:imax, j, kmin:kmax))

#ifdef TRACK_Fillins
            ! Debugging 2022-9-03:
            ! j values are at the gridpoint on different processor runs: e.g. 1 or 4 processors
            IF (iv .EQ. 1 .AND. j .EQ. 17 * 16 + 15 .OR. j .EQ. 34 * 16 + 32) THEN
              WRITE (*, 12) i, j, k, a, X%sg%cell_cntr(:, j), X%sg%mpddInfo_sg%myrank
12            FORMAT('Check the A values: ', 3I4, ' A: ', D10.2, ' LL', 2D15.8, ' pc:', I2)
            END IF
#endif

            IF (A .LT. 1.0D0) THEN
              ! No obs found in this ramp ranges:
              nfill = nfill + 1
              O(nfill) = S%fields(iv)%DATA(i, j, k)
              idc(1, nfill) = i
              idc(2, nfill) = j
              idc(3, nfill) = k
            END IF
          END DO
        END DO
      END DO

      ! Add or create Z ObsField:
      nv = 0
      iu = Y%getObsIdx(TRIM(X%fields(iv)%Get_Name()))
      IF (iu .GT. 0) THEN
        ! Y contains direct obs:
        nv = UBOUND(Y%ObsFields(iu)%values, 1)

        ! Release the data of Z and allocate addition memory for backround fillin:
        IF (ALLOCATED(Z%ObsFields(iu)%values)) DEALLOCATE (Z%ObsFields(iu)%values)
        IF (ALLOCATED(Z%ObsFields(iu)%idx)) DEALLOCATE (Z%ObsFields(iu)%idx)
        ALLOCATE (Z%ObsFields(iu)%values(nv + nfill), Z%ObsFields(iu)%idx(nv + nfill))
        Z%ObsFields(iu)%values(1:nv) = Y%ObsFields(iu)%values(1:nv)
        Z%ObsFields(iu)%values(nv + 1:nv + nfill) = O(1:nfill)
        Z%ObsFields(iu)%idx(1:nv) = Y%ObsFields(iu)%idx(1:nv)
        Z%ObsFields(iu)%idx(nv + 1:nv + nfill)%vIdx = idc(1, 1:nfill)
        Z%ObsFields(iu)%idx(nv + 1:nv + nfill)%hIdx = idc(2, 1:nfill)
        Z%ObsFields(iu)%idx(nv + 1:nv + nfill)%tIdx = idc(3, 1:nfill)

        CALL Z%ObsFields(iu)%Set_Name(TRIM(Y%ObsFields(iu)%Get_Name()))
        CALL Z%ObsFields(iu)%Set_ObsType(TRIM(Y%ObsFields(iu)%Get_ObsType()))
        Z%ObsFields(iu)%locObs = Y%ObsFields(iu)%locObs
      ELSE
        ! idd = idd + 1 ! One more var is additional to obs:
        ! ALLOCATE (Z%ObsFields(nin + idd)%values(nfill), Z%ObsFields(nin + idd)%idx(nfill))
        ! Z%ObsFields(nin + idd)%values(1:nfill) = O(1:nfill)
        ! Z%ObsFields(nin + idd)%idx(1:nfill)%vIdx = idc(1, 1:nfill)
        ! Z%ObsFields(nin + idd)%idx(1:nfill)%hIdx = idc(2, 1:nfill)
        ! Z%ObsFields(nin + idd)%idx(1:nfill)%tIdx = idc(3, 1:nfill)

        ! CALL Z%ObsFields(nin + idd)%Set_Name(TRIM(X%fields(iv)%Get_Name()))
        ! CALL Z%ObsFields(nin + idd)%Set_ObsType('bkgd')
        ! Z%ObsFields(nin + idd)%locObs(1:2) = X%sg%cell_cntr(1:2, 1)
        ! Z%ObsFields(nin + idd)%locObs(3) = 0.0D0 ! Temporarily set
        WRITE (*, 5) TRIM(X%fields(iv)%Get_Name()), X%mpddSub%myrank
5       FORMAT('No observation for this state variable: ', A, 1X, ' No need to fillin, pc: ', I4)
      END IF
    END DO

    ! Deallocate local variables:
    DEALLOCATE (O, idc)

  END SUBROUTINE bkgdFillToObs
  !> @brief
  !==================================================================
  !  This routine maps a gridded obs to a coarser grid.
  !   Input:
  !     sgCurrent: a singleGrid at the current resolution and procs.
  !     sgCoarse:  a singleGrid at a coarser resolution.
  !     yCurrent:  a gridded obs at the current resolution.
  !   Output:
  !     yCoarse:   a mapped gridded obs at a coarser resolution.
  !
  !   Created by Yuanfu Xie. 2022-06-21
  !
  SUBROUTINE griddedObsMap(my_obsbase, sgCurrent, sgCoarse, yCurrent, yCoarse, mpObs)
    IMPLICIT NONE

    CLASS(ObsBase_t) :: my_obsbase
    TYPE(SingleGrid_t), INTENT(IN) :: sgCurrent, sgCoarse
    TYPE(mpObs_t), INTENT(IN) :: mpObs
    TYPE(ObsSet_t), INTENT(IN) :: yCurrent
    TYPE(ObsSet_t), INTENT(OUT) :: yCoarse

    ! Local variables:
    INTEGER(i_kind) :: iv, ih, it, ic, inc(2), istatus, nc, in
    INTEGER(i_kind) :: myMask(sgCurrent%num_cell), nh, nt, numCoarseObs
    INTEGER(i_kind) :: numICells(sgCoarse%mpddGlob%nProc), disps(sgCoarse%mpddGlob%nProc)
    REAL(r_kind), ALLOCATABLE :: coarseGrid(:, :, :), coarseErrs(:, :, :), currentErrs(:), currentGrid(:)
    REAL(r_kind) :: t1, t2

    ! For more general processor setting. Inactive now.
    ! INCLUDE 'mpif.h'

    ALLOCATE (coarseGrid(sgCoarse%vLevel, sgCoarse%num_icell, sgCoarse%tSlots), &
              coarseErrs(sgCoarse%vLevel, sgCoarse%num_icell, sgCoarse%tSlots), &
              currentErrs(sgCurrent%num_icell), &
              currentGrid(sgCurrent%num_icell))

    ! Determine if vertical levels are the same:
    IF (sgCurrent%vLevel .EQ. sgCoarse%vLevel) THEN
      inc(1) = 1
    ELSE IF ((sgCurrent%vLevel - 1) / 2 + 1 .EQ. sgCoarse%vLevel) THEN
      inc(1) = 2
    ELSE
      PRINT *, ''
      WRITE (*, 1) sgCurrent%vLevel, sgCoarse%vLevel
      PRINT *, ''
1     FORMAT('griddedObsMap -- error: vertical levels are not suitable for a multigrid: Abort! ', 2I3)
      STOP
    END IF
    ! Determine if time frames are the same:
    IF (sgCurrent%tSlots .EQ. sgCoarse%tSlots) THEN
      inc(2) = 1
    ELSE IF ((sgCurrent%tSlots - 1) / 2 + 1 .EQ. sgCoarse%tSlots) THEN
      inc(2) = 2
    ELSE
      PRINT *, ''
      WRITE (*, 2) sgCurrent%vLevel, sgCoarse%vLevel
      PRINT *, ''
2     FORMAT('griddedObsMap -- error: time frames are not suitable for a multigrid: Abort! ', 2I3)
      STOP
    END IF

    ! Generating horizontal mask at finer grid for the coarse grid interpolation:
    ! The interpolation happens only if the mask is all 1.
    CALL CPU_TIME(t1)

    ! Gathering the num_icell: only for more general multi processor situation:
    ! Currently inactive:
    !disps = 0
    !CALL MPI_ALLGATHER(sgCoarse%num_icell,1,MPI_INTEGER, &
    !  numICells,1,MPI_INTEGER,sgCoarse%mpddInfo_sg%comm,istatus)
    !DO iv=sgCoarse%mpddGlob%nProc,1,-1
    !  disps(iv) = SUM(numICells(1:iv-1))+1
    !END DO
    !WRITE(*,111) disps,numICells
111 FORMAT('Displacement/NumIcells: ', 20I5)

    ! Initialize yCoarse:
    yCoarse = ObsSet_t(my_obsbase%configfile, mpObs)

    ALLOCATE (yCoarse%ObsFields(UBOUND(yCurrent%ObsFields, 1)))

    ! For all variables:
    DO iv = 1, UBOUND(yCurrent%ObsFields, 1)

      ! Initialize obsFields:
      yCoarse%ObsFields(iv) = ObsField_t(my_obsbase%configFile, mpObs)

      CALL yCoarse%ObsFields(iv)%set_name(yCurrent%ObsFields(iv)%get_name())
      CALL yCoarse%ObsFields(iv)%Set_obsType(yCurrent%ObsFields(iv)%Get_obsType())
      ! yCoarse%ObsFields(iv)%attr = yCurrent%ObsFields(iv)%attr
      yCoarse%ObsFields(iv)%ObsAttrSat = yCurrent%ObsFields(iv)%ObsAttrSat
      yCoarse%ObsFields(iv)%locObs = yCurrent%ObsFields(iv)%locObs

      ! Count the number of gridded obs at a coarse grid:
      numCoarseObs = 0
      coarseGrid = missing
      coarseErrs = missing
      DO ih = 1, sgCoarse%vLevel
        nh = COUNT(yCurrent%ObsFields(iv)%idx%vIdx .EQ. inc(1) * (ih - 1) + 1)
        IF (nh .LE. 0) CYCLE ! No matching vertical level
        DO it = 1, sgCoarse%tSlots
          nt = COUNT(yCurrent%ObsFields(iv)%idx%tIdx .EQ. inc(2) * (it - 1) + 1)
          IF (nt .LE. 0) CYCLE ! No matching temporal frame

          myMask = 0
          ! For all fine obs:
          DO ic = 1, UBOUND(yCurrent%ObsFields(iv)%idx, 1)
            IF (yCurrent%ObsFields(iv)%idx(ic)%vIdx .EQ. inc(1) * (ih - 1) + 1 .AND. &
                yCurrent%ObsFields(iv)%idx(ic)%tIdx .EQ. inc(2) * (it - 1) + 1) THEN
              ! Horizontal mask for all fine obs matches vertical and temporal grids:
              myMask(yCurrent%ObsFields(iv)%idx(ic)%hIdx) = 1

              ! Save the obs to a grid:
              currentGrid(yCurrent%ObsFields(iv)%idx(ic)%hIdx) = &
                yCurrent%ObsFields(iv)%values(ic)
              currentErrs(yCurrent%ObsFields(iv)%idx(ic)%hIdx) = &
                yCurrent%ObsFields(iv)%Errors(ic)
            END IF
          END DO

          ! Coarse grid obs:
          DO ic = 1, sgCoarse%num_icell
            ! For all child grid points:
            nc = 0

            ! Assuming the coarse grid point located at the center of a rectangle:
            DO in = 1, UBOUND(sgCoarse%c_tc_idx(:, ic), 1)
              ! Count only interior cells, c_tc_idx > 0:
              IF (sgCoarse%c_tc_idx(in, ic) .LE. 0) CYCLE
              IF (myMask(sgCoarse%c_tc_idx(in, ic)) .LE. 0) CYCLE

              ! Count the fine grid obs:
              nc = nc + myMask(sgCoarse%c_tc_idx(in, ic))

              ! There is a coarse grid obs and initializing it:
              IF (nc .GT. 0 .AND. coarseGrid(ih, ic, it) .EQ. missing) THEN
                coarseGrid(ih, ic, it) = 0.0D0
                coarseErrs(ih, ic, it) = 0.0D0
              END IF

              ! Accumulate the fine grid obs onto the coarse grid obs:
              coarseGrid(ih, ic, it) = coarseGrid(ih, ic, it) + currentGrid(sgCoarse%c_tc_idx(in, ic))
              coarseErrs(ih, ic, it) = coarseErrs(ih, ic, it) + currentErrs(sgCoarse%c_tc_idx(in, ic))
            END DO

            ! Take the average values:
            IF (nc .GT. 0) THEN
              ! Found a coarse grid point with a child/children on the
              ! finer grid points with gridded obs:
              numCoarseObs = numCoarseObs + 1

              ! Taking the average:
              coarseGrid(ih, ic, it) = coarseGrid(ih, ic, it) / nc
              coarseErrs(ih, ic, it) = coarseErrs(ih, ic, it) / nc

            END IF
          END DO

        END DO ! end of temporal
      END DO  ! end of vertical

      ! Construct the gridded obs at a coarse grid:
      ALLOCATE (yCoarse%ObsFields(iv)%idx(numCoarseObs), &
                yCoarse%ObsFields(iv)%values(numCoarseObs), &
                yCoarse%ObsFields(iv)%errors(numCoarseObs))
      numCoarseObs = 0
      DO ih = 1, sgCoarse%vLevel
        nh = COUNT(yCurrent%ObsFields(iv)%idx%vIdx .EQ. inc(1) * (ih - 1) + 1)
        IF (nh .LE. 0) CYCLE ! Only the vertical level matching with the fine grid

        DO it = 1, sgCoarse%tSlots
          nt = COUNT(yCurrent%ObsFields(iv)%idx%tIdx .EQ. inc(2) * (it - 1) + 1)

          ! Only if height and time match the fine grid:
          IF (nt .LE. 0) CYCLE ! Only the temporal frame matching with the fine grid

          DO ic = 1, sgCoarse%num_icell
            IF (coarseGrid(ih, ic, it) .NE. missing) THEN
              numCoarseObs = numCoarseObs + 1
              yCoarse%ObsFields(iv)%idx(numCoarseObs)%hIdx = ic
              yCoarse%ObsFields(iv)%idx(numCoarseObs)%vIdx = ih
              yCoarse%ObsFields(iv)%idx(numCoarseObs)%tIdx = it
              yCoarse%ObsFields(iv)%values(numCoarseObs) = coarseGrid(ih, ic, it)
              yCoarse%ObsFields(iv)%errors(numCoarseObs) = coarseErrs(ih, ic, it)
            END IF
          END DO
        END DO
      END DO
    END DO
    WRITE (*, 17) sgCurrent%num_icell_global, sgCurrent%num_icell_toCoarser, sgCurrent%num_icell_toFiner, sgCurrent%mpddInfo_sg%myrank
17  FORMAT('GLOBAL ICELL: ', 3I8, ' proc: ', I2)

    DEALLOCATE (coarseGrid, coarseErrs, currentErrs, currentGrid)
    CALL CPU_TIME(t2)
    WRITE (*, 10) t2 - t1
10  FORMAT('Time spent on griddedObs Mapping: ', D12.4)
  END SUBROUTINE griddedObsMap

  !> @brief
  !==================================================================
  !  This routine calculates the O-A for all analysis variables for a
  !   verification purpose.
  !   Input:
  !     state:          state_t with geometry and background fields
  !     verifyLatlon:   latlon region for verification
  !   Output:
  !     numVars:  Number of variables available
  !     varNames: Names of the variables
  !     OmX:      Obs minus analysis for all variables
  !     imx:      The indices where the O-A max location
  !     Xmm:      Max and min values of the analysis
  !
  !   Created by Yuanfu Xie. 2022-06-020
  !
  SUBROUTINE ObsMinusState(my_obsbase, state, verifyLatlon, numVars, OmX, imx, varNames, Xmm)
    IMPLICIT NONE

    CLASS(ObsBase_t) :: my_obsbase
    TYPE(State_t), INTENT(IN)     :: state
    REAL(r_kind), INTENT(IN)      :: verifyLatlon(4) ! 1-2: latMinMax; 3-4: lonMinMax
    INTEGER(i_kind), INTENT(OUT)   :: numVars
    CHARACTER(LEN=20), INTENT(OUT) :: varNames(:)
    INTEGER(i_kind), INTENT(OUT)  :: imx(:)       ! indices of the max O-X of variables:
    REAL(r_kind), INTENT(OUT)     :: OmX(:)       ! max ABS O-X
    REAL(r_kind), INTENT(OUT)     :: Xmm(:, :)

    ! Local variables:
    CHARACTER*20, PARAMETER :: header = 'ObsBase>ObsMinusState: '
    INTEGER(i_kind), PARAMETER :: header_len = 23

    ! Local variables:
    INTEGER(i_kind) :: numObs(my_obsbase%numVars)
    INTEGER(i_kind) :: nhstencil
    TYPE(domainCheck_t) :: domain

    ! Loop variables:
    INTEGER(i_kind) :: io, iv, i, j, k, istatus
    REAL(r_kind) :: bkgdAtObs, t1, t2, t11, t22, t111, t222, bkgd

    ! Consider to save these variables in ObsBase:
    INTEGER(i_kind) :: numValided, maskValided(my_obsbase%numobs)
    INTEGER(i_kind), ALLOCATABLE :: &
      idxGrdValided(:, :), idxHgtValided(:, :, :), idxTimValided(:, :)

    REAL(r_kind) :: domainLatlon(4)

    REAL(r_kind), ALLOCATABLE :: &
      coeGrdValided(:, :), coeHgtValided(:, :, :), coeTimValided(:, :), &
      forward(:, :, :) ! vertical, horizontal and time

    ! All parameters of valided obs:
    REAL(r_kind), ALLOCATABLE :: valided_obsData(:, :), valided_obsErrs(:, :), &
                                 valided_olatlon(:, :), valided_height(:), valided_tempor(:)

    REAL(r_kind) :: sum0 = 0.0D0, allSum, allNumObs
    ! Validating the obs and saving their interpolation coefficients:
    IF (TRIM(my_obsbase%obsType) == 'SOUND') THEN
      CALL domain%validation(state%sg, my_obsbase%numobs, my_obsbase%olatlon, my_obsbase%obsData(:, my_obsbase%presIdx), &
                             my_obsbase%obsTime, numValided, maskValided, &
                             idxGrdValided, idxHgtValided, idxTimValided, &
                             coeGrdValided, coeHgtValided, coeTimValided, &
                             my_obsbase%interpolation, nhstencil, my_obsbase%obsType)
    ELSE
      CALL domain%validation(state%sg, my_obsbase%numobs, my_obsbase%olatlon, my_obsbase%obsHght, &
                             my_obsbase%obsTime, numValided, maskValided, &
                             idxGrdValided, idxHgtValided, idxTimValided, &
                             coeGrdValided, coeHgtValided, coeTimValided, &
                             my_obsbase%interpolation, nhstencil, my_obsbase%obsType)
    END IF

    ! Allocate valided parameter arrays:
    ALLOCATE (valided_obsData(numValided, my_obsbase%numVars), &
              valided_obsErrs(numValided, my_obsbase%numVars), &
              valided_olatlon(2, numValided), &
              valided_height(numValided), &
              valided_tempor(numValided))

    ALLOCATE (forward(2, nhstencil, 2))  ! vertical, horizontal, temporal

    CALL CPU_TIME(t1)

    ! Pass the original arrays to the valided arrays:
    i = 0
    DO io = 1, my_obsbase%numObs
      IF (maskValided(io) .GT. 0) THEN
        i = i + 1
        valided_obsData(i, :) = my_obsbase%obsData(io, 1:my_obsbase%numVars)
        valided_obsErrs(i, :) = my_obsbase%obsErrs(io, 1:my_obsbase%numVars)
        valided_olatlon(:, i) = my_obsbase%olatlon(:, io)
        valided_height(i) = my_obsbase%obsHght(io)
        valided_tempor(i) = my_obsbase%obsTime(io)
      END IF
    END DO

    CALL CPU_TIME(t11)

    ! Thinning valid data only:
    OmX = 0.0D0
    imx = 0
    numVars = 0
    domainLatlon(1) = MINVAL(state%sg%cell_cntr(1, :))
    domainLatlon(2) = MAXVAL(state%sg%cell_cntr(1, :))
    domainLatlon(3) = MINVAL(state%sg%cell_cntr(2, :))
    domainLatlon(4) = MAXVAL(state%sg%cell_cntr(2, :))
    DO iv = 1, my_obsbase%numVars
      IF (.NOT. my_obsbase%Field2ObsIsExisted(state, TRIM(my_obsbase%varNames(iv)))) CYCLE

      ! Accepted variable:
      numVars = numVars + 1
      varnames(numVars) = TRIM(my_obsbase%varNames(iv))

      sum0 = 0.0D0
      DO io = 1, numValided

        ! Valided obs value required:
        IF (ABS(valided_obsData(io, iv)) .GE. ABS(missing)) CYCLE
        ! Only the verification area:
        IF (valided_olatlon(1, io) / degree2radian .GE. verifyLatlon(1) .AND. &
            valided_olatlon(1, io) / degree2radian .LE. verifyLatlon(2) .AND. &
            valided_olatlon(2, io) / degree2radian .GE. verifyLatlon(3) .AND. &
            valided_olatlon(2, io) / degree2radian .LE. verifyLatlon(4) .AND. &
            valided_olatlon(1, io) .GE. domainLatlon(1) .AND. &
            valided_olatlon(1, io) .LE. domainLatlon(2) .AND. &
            valided_olatlon(2, io) .GE. domainLatlon(3) .AND. &
            valided_olatlon(2, io) .LE. domainLatlon(4)) THEN

          ! Background value at obs location if required:
          bkgdAtObs = 0.0D0
          DO i = 1, 2 ! 2 time frames of interpolation
            DO j = 1, nhstencil ! 3 horizontal interpolation points
              DO k = 1, 2 ! 2 vertical levels
                forward(k, j, i) = my_obsbase%GetForwardValue(state, varnames(numVars), &
                                                              idxHgtValided(k, j, io), idxGrdValided(j, io), idxTimValided(i, io), iv)
                ! Interpolate the background at observation:
                bkgdAtObs = bkgdAtObs + forward(k, j, i) * &
                            coeHgtValided(k, j, io) * coeGrdValided(j, io) * coeTimValided(i, io)
              END DO
            END DO
          END DO

          sum0 = sum0 + ABS(valided_obsData(io, iv) - bkgdAtObs)
          ! Max error at iv and io:
          IF (OmX(numVars) .LT. ABS(valided_obsData(io, iv) - bkgdAtObs)) THEN
            OmX(numVars) = ABS(valided_obsData(io, iv) - bkgdAtObs)
            imx(numVars) = io
          END IF

        END IF ! Area check

      END DO

      PRINT *, 'Error average (A-O): ', sum0 / numValided, varnames(numVars)
    END DO

    ! Calculate the max/min values of X:
    ! For the variable:
    Xmm(1, 1:UBOUND(Xmm, 2)) = -1.0D10
    Xmm(2, 1:UBOUND(Xmm, 2)) = 1.0D10
    DO i = 1, state%sg%num_cell

      ! Only the verification area:
      IF (state%sg%cell_cntr(1, i) / degree2radian .GE. verifyLatlon(1) .AND. &
          state%sg%cell_cntr(1, i) / degree2radian .LE. verifyLatlon(2) .AND. &
          state%sg%cell_cntr(2, i) / degree2radian .GE. verifyLatlon(3) .AND. &
          state%sg%cell_cntr(2, i) / degree2radian .LE. verifyLatlon(4) .AND. &
          state%sg%cell_cntr(1, i) .GE. domainLatlon(1) .AND. &
          state%sg%cell_cntr(1, i) .LE. domainLatlon(2) .AND. &
          state%sg%cell_cntr(2, i) .GE. domainLatlon(3) .AND. &
          state%sg%cell_cntr(2, i) .LE. domainLatlon(4)) THEN

        ! Calculate the X max and min values over the verification area:
        DO iv = 1, numVars

          ! State value at obs location if required:
          DO j = 1, state%sg%vLevel
            DO k = 1, state%sg%tSlots
              bkgdAtObs = my_obsbase%GetForwardValue(state, varnames(iv), j, i, k, iv)

              IF (bkgdAtObs .GT. Xmm(1, iv)) Xmm(1, iv) = bkgdAtObs
              IF (bkgdAtObs .LT. Xmm(2, iv)) Xmm(2, iv) = bkgdAtObs
            END DO
          END DO

        END DO

      END IF

    END DO

    CALL CPU_TIME(t22)
    WRITE (*, 1) t22 - t11, state%sg%mpddInfo_sg%myrank
1   FORMAT('Time spent on ObsMinusState: ', D12.4, ' proc: ', I2)

  END SUBROUTINE ObsMinusState

  SUBROUTINE ObsControl(configFile, thinObs, X)
    IMPLICIT NONE

    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(ObsSet_t), INTENT(INOUT) :: thinObs
    TYPE(State_t), INTENT(IN) :: X

    ! may change to parameters in the furture:
    TYPE(mpObs_t) :: mpObs

    ! local variables
    TYPE(ObsConvention_t2) :: my_conv
    TYPE(ObsSet_t) :: local_thinObs(6)
    CHARACTER(LEN=10) :: typeName(6)
    INTEGER(i_kind) :: i, j, status, num_thin(6), j_idx
    LOGICAL :: typeExist(6)

    typeName(1) = "SYNOP"
    typeName(2) = "SHIP"
    typeName(3) = "BUOY"
    typeName(4) = "METAR"
    typeName(5) = "TEMP"
    typeName(6) = "PROFL"

    num_thin = 0
    DO i = 1, 6
      status = yaml_get_var(configFile, "obsList", typeName(i), typeExist(i))
      IF (typeExist(i)) THEN
        CALL my_conv%obsInitial(configFile, typeName(i))
        CALL my_conv%obsIngest(X)
        IF (UBOUND(my_conv%obsData, 1) .EQ. 0) THEN
          PRINT *, "Warnning!!! No data is loaded for ", my_conv%obsType
          CYCLE
        END IF

        CALL my_conv%obsQC(X)
        IF (UBOUND(my_conv%obsData, 1) .EQ. 0) THEN
          PRINT *, "Warnning!!! All data have been removed after QC for ", my_conv%obsType
          CYCLE
        END IF

        CALL my_conv%obsThinning(X, local_thinObs(i), mpObs, .FALSE.)
        num_thin(i) = SIZE(local_thinObs(i)%ObsFields)
        CALL my_conv%obsDeallocate()
      END IF
    END DO

    ! Concatenate thinned observation data
    ALLOCATE (thinObs%ObsFields(SUM(num_thin)))
    j_idx = 0
    DO i = 1, 6
      IF (typeExist(i)) THEN
        DO j = 1, num_thin(i)
          j_idx = j_idx + 1
          thinObs%ObsFields(j_idx) = local_thinObs(i)%ObsFields(j)
        END DO
      END IF
    END DO

    WRITE (*, *) 'Number of thinned observation variables: ', SIZE(thinObs%ObsFields)
  END SUBROUTINE ObsControl


  !> @brief
  !================================================================
  !! This routine checks the obsSet Y to fill in background in area
  !! where no observation associated with the background field.
  SUBROUTINE ObsSelection(Y, ObsSelectionConfig)
    IMPLICIT NONE
    TYPE(ObsSet_t), INTENT(INOUT) :: Y
    CHARACTER(LEN=*), INTENT(IN) :: ObsSelectionConfig

    

  END SUBROUTINE ObsSelection

END MODULE obsTools_m
