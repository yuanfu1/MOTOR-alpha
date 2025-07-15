!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com.com), 2021/10/27, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This is a module for a unit test of MOTOR-DA surface analysis
!!   it provides an analytic multiscale function to test and testing capability routine
!
MODULE unitTest_sfc_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE ObsSurface_m, ONLY: ObsSurface_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE parameters_m, ONLY: degree2radian

  TYPE unitTest_sfc_t
  CONTAINS
    PROCEDURE analytic
    PROCEDURE checking
    PROCEDURE, NOPASS:: randomObs
  END TYPE unitTest_sfc_t

CONTAINS
  ! Analytic function is defined:
  SUBROUTINE analytic(this, npts, xyt, func)

    IMPLICIT NONE

    CLASS(unitTest_sfc_t) :: this
    INTEGER(i_kind), INTENT(IN) :: npts
    REAL(r_kind), INTENT(IN) :: xyt(3, npts)
    REAL(r_kind), INTENT(OUT) :: func(npts)

    ! Local variables:
    INTEGER(i_kind) :: i
    REAL(r_kind) :: speed, pi, st, rg, fr, x, y, r, r2

    !speed = 10.0  ! About 74mph  if st=-14 rg=28 fr=20
    !speed = 7.6   ! About 90km/h if st=-14 rg=28 fr=20
    speed = 4.0   ! About 90km/h
    ! speed = 8.0   ! About 90km/h

    pi = 4.0 * ATAN(1.0)

    st = -3.5
    rg = 7.0
    fr = 5.0
    ! fr = 10.0

    ! For all observation locations, replace the obs values with the analytic:
    DO i = 1, npts

      x = st + rg * xyt(1, i)
      y = st + rg * xyt(2, i)
      r = ABS(x - y + fr + speed * xyt(3, i))

      r2 = SQRT((x + st)**2 + (y - st)**2)
      r = r - 0.0 - (0.0 + 0.8 * r2) * SQRT(2.0)

      IF (ABS(r) .GT. 1.0E-8) THEN
        func(i) = 0.8*(1.0+2.0*(1.0+TANH(r))*(1.0+1.0/r*SIN(r)**2))
        ! func(i) = 0.8 * (1.0 + 2.0 * (1.0 + TANH(r)) * (1.0 + 1.0 / r * SIN(10 * r)))
      ELSE
        func(i) = 0.8
      END IF
      ! func(i) = SIN(10*r)+1

    END DO

  END SUBROUTINE analytic

  ! Checking: choosing a set of random locations to check the analysis:
  SUBROUTINE checking(this, npts, xyt, truth, soltn)
    CLASS(unitTest_sfc_t) :: this
    INTEGER(i_kind), INTENT(IN) :: npts
    REAL(r_kind), INTENT(IN) :: xyt(3, npts)
    REAL(r_kind), INTENT(IN) :: truth(npts), soltn(npts)

    ! Local variables:
    INTEGER(i_kind) :: ipts, imx
    REAL(r_kind) :: amx

    amx = 0.0D0
    DO ipts = 1, npts
      IF (amx .LT. ABS(truth(ipts) - soltn(ipts))) THEN
        amx = ABS(truth(ipts) - soltn(ipts))

        imx = ipts
      END IF
    END DO

    ! Write out messages:
    WRITE (*, 1) amx, imx
1   FORMAT('Maximum analysis error: ', D12.4, ' at cell: ', I8)
  END SUBROUTINE checking

  SUBROUTINE randomObs(this, numObs, sg)
    IMPLICIT NONE

    TYPE(ObsSurface_t), INTENT(INOUT) :: this
    INTEGER(i_kind), INTENT(IN) :: numObs
    TYPE(SingleGrid_t), INTENT(IN) :: sg

    ! Local variables:
    INTEGER(i_kind) :: i, io, it
    REAL(r_kind) :: rans, aminll(2), amaxll(2)
    REAL(r_kind) :: minGlobLat, maxGlobLat, minGlobLon, maxGlobLon
    REAL(r_kind), ALLOCATABLE :: xyt(:, :)

    TYPE(unitTest_sfc_t) :: unitTest

    this%numObs = numObs * sg%tSlots  ! Assuming every time frame has the same obs

    ALLOCATE (this%obsData(this%numObs, this%numVars))
    ALLOCATE (this%obsErrs(this%numObs, this%numVars))
    ALLOCATE (this%olatlon(2, this%numObs))
    ALLOCATE (this%obsTime(this%numObs))
    ALLOCATE (this%obsHght(this%numObs))
    ALLOCATE (this%land(this%numObs))
    ALLOCATE (this%StNames(this%numObs))

    ALLOCATE (xyt(3, this%numObs))

    BLOCK
      REAL(r_kind) :: swap

      swap = MINVAL(sg%cell_cntr(1, 1:sg%num_icell))
      CALL sg%mpddInfo_sg%AllReduceMinReal(swap, minGlobLat)

      swap = MAXVAL(sg%cell_cntr(1, :))
      CALL sg%mpddInfo_sg%AllReduceMaxReal(swap, maxGlobLat)

      swap = MINVAL(sg%cell_cntr(2, :))
      CALL sg%mpddInfo_sg%AllReduceMinReal(swap, minGlobLon)

      swap = MAXVAL(sg%cell_cntr(2, :))
      CALL sg%mpddInfo_sg%AllReduceMaxReal(swap, maxGlobLon)
    END BLOCK

    ! Domain:
    aminll(1) = MINVAL(sg%cell_cntr(1, :))
    aminll(2) = MINVAL(sg%cell_cntr(2, :))
    amaxll(1) = MAXVAL(sg%cell_cntr(1, :))
    amaxll(2) = MAXVAL(sg%cell_cntr(2, :))
    WRITE (*, 4) aminll(:) / degree2radian, amaxll(:) / degree2radian
4   FORMAT('ranges of LL: ', 4D12.4)
    WRITE (*, 5) minGlobLat / degree2radian, maxGlobLat / degree2radian, &
      minGlobLon / degree2radian, maxGlobLon / degree2radian
5   FORMAT('ranges of Global LL: ', 4D12.4)

    ! For all time frames:
    io = 0
    DO it = 1, sg%tSlots
      DO i = 1, numObs
        io = io + 1
        CALL RANDOM_NUMBER(xyt(1, io))
        CALL RANDOM_NUMBER(xyt(2, io))
        this%olatlon(1, io) = aminll(1) + (amaxll(1) - aminll(1)) * xyt(1, io)
        this%olatlon(2, io) = aminll(2) + (amaxll(2) - aminll(2)) * xyt(2, io)
        xyt(1, io) = (this%olatlon(1, io) - minGlobLat) / (maxGlobLat - minGlobLat)
        xyt(2, io) = (this%olatlon(2, io) - minGlobLon) / (maxGlobLon - minGlobLon)

        xyt(3, io) = (sg%tt(it) - sg%tt(1)) / (sg%tt(sg%tSlots) - sg%tt(1))
      END DO
    END DO
    WRITE (*, 6) MINVAL(xyt(1, :)), MAXVAL(xyt(1, :)), &
      MINVAL(xyt(2, :)), MAXVAL(xyt(2, :)), &
      MINVAL(xyt(3, :)), MAXVAL(xyt(3, :)), sg%mpddInfo_sg%myRank
6   FORMAT('XYT ranges: ', 6D12.4, ' at Procs: ', I1)
    CALL unitTest%analytic(this%numObs, xyt, this%obsData(:, 1))

    WRITE (*, 8) this%obsData(1:3, 1), sg%mpddInfo_sg%myRank
8   FORMAT('OBS values: ', 3D12.4, ' at procs', I1)
    WRITE (*, 7) MINVAL(this%obsData(:, 1)), MAXVAL(this%obsData(:, 1)), sg%mpddInfo_sg%myRank
7   FORMAT('Min/Max OBS: ', 2D12.4, ' at Procs: ', I1)
!stop

    DO io = 1, this%numObs
      this%obsTime(io) = xyt(3, io) * (sg%tt(sg%tSlots) - sg%tt(1)) + sg%tt(1)

      this%obsErrs(io, :) = 1.0D0
      this%obsHght(io) = 0.0D0
      this%land(io) = 1.0D0
!       WRITE(this%StNames(io),1) io
! 1     FORMAT('Random',I6.6)
    END DO

    ! Convert obs latlon to degree:
    ! this%olatlon = this%olatlon/degree2radian

    !write(*,2) this%obsData(1,1),this%olatlon(:,1)/degree2radian,this%obsTime(1),this%obsTime(this%numObs)
2   FORMAT('Test obs: ', D12.4, ' latlon: ', 2D12.4, ' at: ', 2I11)
    !stop

    ! Deallocate local memory:
    DEALLOCATE (xyt)

  END SUBROUTINE randomObs
END MODULE unitTest_sfc_m
