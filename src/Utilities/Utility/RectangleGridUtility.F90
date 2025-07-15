!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-Repository
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2022-04-10   Created by Yuanfu Xie
!
!   Created by Yuanfu Xie (xieyf@gbamwf.com), 2022/04/10, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides basic utility functions designed for more efficiently handle regular grid.

MODULE RectangleGridUtility_m

  USE kinds_m, ONLY: i_kind, r_kind, r_double

  TYPE RectangleGridUtility_t
  CONTAINS
    PROCEDURE :: RectangleHorizontalIntp
  END TYPE RectangleGridUtility_t

CONTAINS

  !> @brief
  !------------------------------------------------------------------
  !  This routine interpolates a set of points assuming the underlying
  !  grid is regular, rectangular including latlon grid inside a
  !  singleGrid.
  !
  !  Assuming the srcll arranged as by rows, i+(numsrc(1))*(j-1) and
  !  the grid is parallel to latlon from west to east and south to
  !  north.
  !
  !  NOTE: The grid contains either fictitious points or halo points.
  SUBROUTINE rectangleHorizontalIntp(this, srcll, numsrc, tgtll, numtgt, &
                                     idxtgt, coetgt, iproc)
    USE parameters_m, ONLY: machineEps, degree2radian

    IMPLICIT NONE

    CLASS(RectangleGridUtility_t) :: this
    INTEGER(i_kind), INTENT(IN) :: numsrc, numtgt, iproc
    REAL(r_kind), INTENT(IN) :: srcll(numsrc, 2), tgtll(numtgt, 2)
    INTEGER(i_kind), INTENT(OUT) :: idxtgt(4, numtgt)
    REAL(r_kind), INTENT(OUT) :: coetgt(4, numtgt)

    ! Local variables:
    INTEGER(i_kind) :: ng(2), ix, iy, i, j, interior(2), interiorSure, interval(2)
    LOGICAL :: inBox
    REAL(r_kind) :: dlat, dlon, x, y, cornerLL(2, 2), southwestCorner(2)

    INTEGER(i_kind) :: iobs
    REAL(r_kind) :: obsll(2)

!#define TRACE_CELL
    ! Debugging:
    INTEGER :: imx
    REAL(r_kind) :: amx

#ifdef TRACE_CELL
    ! Debugging:
    REAL(r_kind) :: given_corner(2)
    given_corner(1) = -0.96494797378D-01
    given_corner(2) = 0.20833190282D+00
#endif

    ! Set latlon corners: bottom-left and top-right:
    cornerLL(1, 1) = MINVAL(srcll(:, 1))
    cornerLL(2, 1) = MINVAL(srcll(:, 2))
    cornerLL(1, 2) = MAXVAL(srcll(:, 1))
    cornerLL(2, 2) = MAXVAL(srcll(:, 2))

    ! Search the number of gridpoint along longitude: this should only end at interior points
    interior(1) = 1
    DO i = 2, numsrc
      interior(1) = interior(1) + 1   ! Count x direction points
      IF (ABS(srcll(interior(1), 1) - srcll(1, 1)) .GT. 1.0D5 * machineEps) EXIT ! Found different latitude, end first row
    END DO
    interior(1) = interior(1) - 1

#ifdef TRACE_CELL
    ! Debugging
    IF (iproc .EQ. 1) WRITE (*, 441) numsrc, interior(1)
441 FORMAT('Interior--: ', 3I6)
#endif

    ! Check to see if the number of interior points in x direction is acceptable:
    IF (interior(1) .GE. numsrc .OR. interior(1) .EQ. 1) THEN
      WRITE (*, 1) interior(1), numsrc
      STOP
    END IF

    ! Left and right halos:
    ng(1) = interior(1)
    IF (cornerLL(2, 1) .LT. srcll(1, 2)) ng(1) = ng(1) + 1  ! Left
    IF (cornerLL(2, 2) .GT. srcll(interior(1), 2)) ng(1) = ng(1) + 1  ! Right

    IF (MOD(numsrc, ng(1)) .NE. 0) THEN
      WRITE (*, 1) ng(1), numsrc, MOD(numsrc, ng(1))
      STOP
    END IF
1   FORMAT('This grid does not seem a regular parallel to latlon, check and rerun: ', I5, ' numSrc: ', 2I5)

    ! Find a regular latlon grid's number of gridpoints:
    ng(2) = numsrc / ng(1)    ! ng(2) counting all halo points in y direction

    ! Bottom halo: we cannot determine if there is top halo
    interior(2) = ng(2)
    IF (cornerLL(1, 1) .LT. srcll(1, 1)) interior(2) = interior(2) - 1  ! bottom

    ! We are not sure if there is a halo on top of the domain,
    ! but sure 1: interior(1)*(interior(2)-1) are all interiors:
    interiorSure = interior(1) * (interior(2) - 1)

#ifdef TRACE_CELL
    ! Debugging
    IF (iproc .EQ. 1) WRITE (*, 444) numsrc, interior, ng
444 FORMAT('Interior++: ', 3I6, ' NG: ', 2I4)
#endif

    dlon = srcll(2, 2) - srcll(1, 2)                ! two adjacent points longitude as Dlon
    dlat = srcll(interior(1) + 1, 1) - srcll(1, 1)    ! Use second row lat substract first row as Dlat
    IF (dlat .LE. 0.0D0 .OR. dlon .LE. 0.0D0) THEN
      WRITE (*, 2) dlat, dlon
2     FORMAT('RectangleHorizontalIntp: ', 2D12.4, ' currently only Lat-Lon grid is implemented! ')
    END IF

#ifdef TRACE_CELL
    ! Debugging:
    WRITE (*, 10) iproc, ng, MINVAL(srcll(:, 1)), MAXVAL(srcll(:, 1)), MINVAL(srcll(:, 2)), MAXVAL(srcll(:, 2))
10  FORMAT('Regular - proc: ', I1, ' Number of grids: ', 2I6, ' minMAXLL: ', 4D10.3)
#endif

    ! Tail of the grid: begin and end: there is possible 1 halo on the top in y direction:
    interval(1) = interior(1) * (interior(2) - 2) + 1 ! interior(2)-2 is to count the top interior row
    interval(2) = numsrc
    interval(1) = MAX(1, interval(1))

    ! Debugging
    IF (iproc .EQ. 1) WRITE (*, 411) numsrc, interval, interior
411 FORMAT('Interval: ', 3I6, ' interior: ', 2I4)

    ! As an initial try, we implement lat-lon grid here:
    idxtgt = 0
    coetgt = 0.0D0

#ifdef TRACE_CELL
    !debugging: 2022-08-26:
    iobs = 1 !92 ! The i-th obs to check
    obsll(1) = 0.225461921691895D+02    ! in degree
    obsll(2) = 0.113818359375000D+03    ! in degree
    obsll(1) = 0.21450D+02    ! in degree
    obsll(2) = 0.10918D+03    ! in degree

    !if (iproc .eq. 1) &
    WRITE (*, 115) tgtll(iobs, :), numtgt, machineEps, iproc
115 FORMAT('TRACE_CELL i-obsLL: ', 2D12.5, ' numtgt: ', I6, ' machineEps:', D15.8, ' pc:', I2)
#endif

    DO i = 1, numtgt

      IF (tgtll(i, 1) .LT. cornerLL(1, 1) .OR. tgtll(i, 1) .GT. cornerLL(1, 2) .OR. &
          tgtll(i, 2) .LT. cornerLL(2, 1) .OR. tgtll(i, 2) .GT. cornerLL(2, 2)) THEN
        CYCLE ! Outside of the domain
      ELSE
        ! Use of all interior points to interpolate:
        IF (tgtll(i, 1) .GE. srcll(1, 1) - machineEps .AND. &
            tgtll(i, 2) .GE. srcll(1, 2) - machineEps .AND. &
            tgtll(i, 1) .LE. srcll(interiorSure, 1) + machineEps .AND. &    ! Maybe halo on the top
            tgtll(i, 2) .LE. srcll(interior(1), 2) + machineEps) THEN
          ! The most halo rows in y direction is 2 so (ng(2)-2) is used:

          ! The initial indices of the array are interior points only:
          ix = INT8((tgtll(i, 2) - srcll(1, 2)) / dlon) + 1    ! Count the interior gridpoints in x
          iy = INT8((tgtll(i, 1) - srcll(1, 1)) / dlat) + 1    ! Count the interior gridpoints in y

          x = MOD((tgtll(i, 2) - srcll(1, 2)), dlon) / dlon
          y = MOD((tgtll(i, 1) - srcll(1, 1)), dlat) / dlat

          idxtgt(1, i) = ix + interior(1) * (iy - 1)
          idxtgt(2, i) = ix + 1 + interior(1) * (iy - 1)
          idxtgt(3, i) = ix + interior(1) * (iy)
          idxtgt(4, i) = ix + 1 + interior(1) * (iy)
          coetgt(1, i) = (1.0D0 - x) * (1.0D0 - y)
          coetgt(2, i) = x * (1.0D0 - y)
          coetgt(3, i) = (1.0D0 - x) * y
          coetgt(4, i) = x * y

#ifdef TRACE_CELL
          !debugging: 2022-08-26:
          IF (ABS(tgtll(i, 1) / degree2radian - obsll(1)) .LE. 1.0D-4 .AND. &
              ABS(tgtll(i, 2) / degree2radian - obsll(2)) .LE. 1.0D-4) THEN
            WRITE (*, 111) ix, iy, x, y, idxtgt(:, i), coetgt(:, i), iproc
111         FORMAT('TRACE_CELL interior - ixy:', 2I4, ' xy:', 2D22.15, ' idx', 4I6, ' coe', 4D12.5, ' pc:', I2)
          END IF
#endif

        ELSE
          ! Point is outside of the interior domain, search the halo region:
          inBox = .FALSE.
          x = MOD((tgtll(i, 2) - cornerLL(2, 1)), dlon) / dlon
          y = MOD((tgtll(i, 1) - cornerLL(1, 1)), dlat) / dlat
          southwestCorner(1) = INT((tgtll(i, 1) - cornerLL(1, 1)) / dlat) * dlat + cornerLL(1, 1)
          southwestCorner(2) = INT((tgtll(i, 2) - cornerLL(2, 1)) / dlon) * dlon + cornerLL(2, 1)

#ifdef TRACE_CELL
          ! debugging:
          IF (iproc .EQ. 3 .AND. ABS(x - given_corner(1)) .LE. 1.0D-8 .AND. ABS(y - given_corner(2)) .LE. 1.0D-8) &
            !IF (iproc .EQ. 1) &
            WRITE (*, 555) tgtll(i, 1:2), cornerLL(1:2, 1), cornerLL(1:2, 2), &
            southwestCorner, srcll(interiorSure, 1), dlat, dlon, (tgtll(i, 1) .LE. srcll(interiorSure, 1))
555       FORMAT('BottomEdge: tgt: ', 2D22.14, ' corn:', 4D16.8, ' SW: ', 2D16.8, ' SureLat: ', 1D22.14, ' dll: ', 2D16.8, ' LE?', L)

          !debugging: 2022-08-26:
          IF (ABS(tgtll(i, 1) / degree2radian - obsll(1)) .LE. 1.0D-4 .AND. &
              ABS(tgtll(i, 2) / degree2radian - obsll(2)) .LE. 1.0D-4) THEN
            WRITE (*, 222) ix, iy, x, y, tgtll(i, 1:2), srcll(1, 1:2), iproc
222         FORMAT('TRACE_CELL halo - ixy:', 2I4, ' xy:', 2D22.15, ' tgt', 2D12.4, ' src', 2D12.4, ' pc:', I2)
          END IF
#endif

          ! Search the possible halo and its adjacent interior points:
          !
          ! Idxtgt order:
          !   4   3
          !   1   2
          !
          ! 4 sides:
          ! 1. Bottom halo:
          IF (tgtll(i, 1) .LT. srcll(1, 1)) THEN
            ! It appears at the bottom halo:
            IF (tgtll(i, 2) .LT. srcll(1, 2)) THEN
              ! It appears at the bottom left corner:
              idxtgt(3, i) = 1 ! cell 1 is the top right corner of the box: 3
              inBox = .TRUE.
            ELSE IF (tgtll(i, 2) .GE. srcll(1, 2) .AND. tgtll(i, 2) .LE. srcll(interior(1), 2)) THEN
              ! It appears at the bottom halo but not at corners:
              idxtgt(3, i) = (tgtll(i, 2) - srcll(1, 2)) / dlon + 2
              idxtgt(4, i) = (tgtll(i, 2) - srcll(1, 2)) / dlon + 1 ! count grid points from the interior
              inBox = .TRUE.
            ELSE
              ! It appears at the bottom right corner:
              idxtgt(4, i) = interior(1) ! srcll(1,*) is the top left corner of the cell
              inBox = .TRUE.
            END IF
            ! 2. Top side:
          ELSE IF (tgtll(i, 1) .GT. srcll(interiorSure, 1)) THEN
            ! It appears on the top possible halo, &
            ! "possible" because ng(2)-2:numsrc may contains internal points
            ! Since we are not sure where the halos are, we simply search all 4 points:
            inBox = .TRUE.
            ! 3. Left side:
          ELSE IF (tgtll(i, 2) .LE. srcll(1, 2)) THEN
            ! It appears on the left side halo region:
            ! We know the two corners on the right wide, but not the left side of the cell:
            idxtgt(2, i) = 1 + interior(1) * INT((tgtll(i, 1) - srcll(1, 1)) / dlat)
            idxtgt(3, i) = 1 + interior(1) * INT((tgtll(i, 1) - srcll(1, 1)) / dlat + 1)
            inBox = .TRUE.
            ! 4. Right side:
          ELSE IF (tgtll(i, 2) .GE. srcll(interior(1), 2)) THEN
            idxtgt(1, i) = interior(1) + interior(1) * INT((tgtll(i, 1) - srcll(1, 1)) / dlat)
            idxtgt(4, i) = interior(1) + interior(1) * INT((tgtll(i, 1) - srcll(1, 1)) / dlat + 1)
            inBox = .TRUE.
          END IF

          ! Search the grid points and calcualte the interpolation coefficients:
          IF (inBox) THEN
            CALL cellSearch(4, srcll, x, y, dlon, dlat, interval, southwestCorner, idxtgt(:, i), coetgt(:, i), iproc)
          ELSE
            ! Debugging:
            ! write(*,13) i,tgtll(i,1:2),srcll(1,1:2),cornerLL(1:2)
13          FORMAT('outBox!: ', I2, ' obsLL: ', 2D14.7, ' grdLL: ', 2D14.7, ' Corner: ', 2D14.7)
          END IF

#ifdef TRACE_CELL
          !debugging: 2022-08-26:
          IF (ABS(tgtll(i, 1) / degree2radian - obsll(1)) .LE. 1.0D-4 .AND. &
              ABS(tgtll(i, 2) / degree2radian - obsll(2)) .LE. 1.0D-4) THEN
            WRITE (*, 113) x, y, idxtgt(:, i), coetgt(:, i), iproc
113         FORMAT('TRACE_CELL in halo: xy: ', 2D22.15, ' idx', 4I4, ' coe', 4D12.5, ' pc:', I2)
          END IF
#endif

        END IF
      END IF
    END DO

#ifdef TRACE_CELL
    !debugging: 2022-08-26:
    IF (iproc .EQ. 1) WRITE (*, 112) idxtgt(:, iobs), coetgt(:, iobs), iproc
112 FORMAT('TRACE_CELL idx: ', 4I6, ' coe', 4D12.5, ' pc:', I2)

    WRITE (*, 123) MAXVAL(coetgt), MINVAL(coetgt), iproc
123 FORMAT('MaxMin coetgt: ', 2E14.4, ' pc: ', I2)

    imx = 0
    amx = 0.0D0
    DO i = 1, numtgt
      IF (amx .LT. ABS(coetgt(1, i))) THEN
        imx = i
        amx = ABS(coetgt(1, i))
      END IF
    END DO
    PRINT *, 'XXX: ', imx, amx, numtgt, UBOUND(coetgt, 2)
#endif

  END SUBROUTINE rectangleHorizontalIntp

  SUBROUTINE cellSearch(np, srcll, x, y, dlon, dlat, interval, southwestCorner, idxtgt, coetgt, iproc)
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN) :: np, interval(2), iproc
    REAL(r_kind), INTENT(IN) :: srcll(interval(2), 2), x, y, dlon, dlat, southwestCorner(2)
    INTEGER(i_kind), INTENT(OUT) :: idxtgt(np)
    REAL(r_kind), INTENT(OUT) :: coetgt(np)

    ! Local variables:
    INTEGER(i_kind) :: j
    REAL(r_kind), PARAMETER :: threshold = 1.0D-6 ! 1 meter resolution, assuming 1 degree < 1000 km, 1km = 1000meter

    ! For given idxtgt values:
    IF (idxtgt(1) .NE. 0) coetgt(1) = (1.0D0 - x) * (1.0D0 - y)
    IF (idxtgt(2) .NE. 0) coetgt(2) = (x) * (1.0D0 - y)
    IF (idxtgt(3) .NE. 0) coetgt(3) = (x) * (y)
    IF (idxtgt(4) .NE. 0) coetgt(4) = (1.0D0 - x) * (y)

    DO j = interval(1), interval(2)
      IF (idxtgt(1) .EQ. 0 .AND. &
          ABS(srcll(j, 1) - southwestCorner(1)) .LT. threshold .AND. &
          ABS(srcll(j, 2) - southwestCorner(2)) .LT. threshold) THEN
        idxtgt(1) = j
        coetgt(1) = (1.0D0 - x) * (1.0D0 - y)
        IF (x .LE. threshold) THEN
          ! Fall on to this edge of y:
          idxtgt(2:3) = j - (/1, 2/) ! Fake, it will be used with zero interpolation coefficient
          coetgt(2:3) = 0.0D0
        END IF
        IF (y .LE. threshold) THEN
          ! Fall on to this edge of x:
          idxtgt(3:4) = j - (/1, 2/) ! Fake, it will be used with zero interpolation coefficient
          coetgt(3:4) = 0.0D0
        END IF
      END IF
      IF (idxtgt(2) .EQ. 0 .AND. &
          ABS(srcll(j, 1) - southwestCorner(1)) .LT. threshold .AND. &
          ABS(srcll(j, 2) - southwestCorner(2) - dlon) .LT. threshold) THEN
        idxtgt(2) = j
        coetgt(2) = (x) * (1.0D0 - y)
      END IF
      IF (idxtgt(3) .EQ. 0 .AND. &
          ABS(srcll(j, 1) - southwestCorner(1) - dlat) .LT. threshold .AND. &
          ABS(srcll(j, 2) - southwestCorner(2) - dlon) .LT. threshold) THEN
        idxtgt(3) = j
        coetgt(3) = (x) * (y)
      END IF
      IF (idxtgt(4) .EQ. 0 .AND. &
          ABS(srcll(j, 1) - southwestCorner(1) - dlat) .LT. threshold .AND. &
          ABS(srcll(j, 2) - southwestCorner(2)) .LT. threshold) THEN
        idxtgt(4) = j
        coetgt(4) = (1.0D0 - x) * (y)
      END IF
    END DO

#ifdef TRACE_CELL
    ! Check idxtgt:
    IF (MINVAL(idxtgt) .LE. 0) THEN
      PRINT *, ''
      PRINT *, '**************** Irregular exits ****************'
      WRITE (*, 1) interval, iproc, idxtgt, &
        MINVAL(srcll(interval(1):interval(2), 1)), &
        MAXVAL(srcll(interval(1):interval(2) - 4, 1)), &
        MINVAL(srcll(interval(1):interval(2), 2)), &
        MAXVAL(srcll(interval(1):interval(2) - 4, 2)), southwestCorner, x, y
1     FORMAT('cellSearch:', 2I6, ' proc: ', I1, ' invalid index: ', 4I6, ' src: ', 4D12.4 &
             ' XYcorner: ', 2D18.11, ' alpha beta: ', 2D18.11)
      PRINT *, ''
      STOP
    END IF
#endif

  END SUBROUTINE cellSearch

END MODULE RectangleGridUtility_m
