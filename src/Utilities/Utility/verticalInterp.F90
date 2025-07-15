!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-PS.TimeIntegration.TimeIntegrationAB3
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2024-09-12   Created by Yuanfu Xie
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides a general vertical interpolation module for 1D, 2D, 3D, and 4D grids.
MODULE verticalInterp_m
  USE, INTRINSIC :: iso_c_binding, ONLY: c_double, c_int

  TYPE verticalInterp_t
    CONTAINS
      GENERIC, PUBLIC :: verticalInterp =>  verticalInterp1D, &
                                            verticalInterp2D, &
                                            verticalInterp3D
      PROCEDURE, PUBLIC :: verticalInterp1D
      PROCEDURE, PUBLIC :: verticalInterp2D
      PROCEDURE, PUBLIC :: verticalInterp3D
  END TYPE verticalInterp_t

  CONTAINS
  SUBROUTINE verticalInterp1D(this, fin, vin, fout, vout, istatus, logOpt)
    CLASS(verticalInterp_t), INTENT(IN) :: this
    REAL(c_double), DIMENSION(:), INTENT(IN) :: fin
    REAL(c_double), DIMENSION(:), INTENT(IN) :: vin
    REAL(c_double), DIMENSION(:), INTENT(IN) :: vout
    REAL(c_double), DIMENSION(:), INTENT(OUT) :: fout
    INTEGER(c_int), INTENT(OUT) :: istatus
    INTEGER(c_int),INTENT(IN), OPTIONAL :: logOpt

    ! Local variables:
    INTEGER(c_int) :: levels_in, levels_out
    INTEGER(c_int) :: i, j, im1, ip1
    REAL(c_double), ALLOCATABLE :: fint(:), vint(:)

    levels_in = SIZE(fin)
    levels_out = SIZE(fout)
    istatus = 0
    IF (SIZE(vin) /= levels_in .OR. SIZE(vout) /= levels_out) THEN
      istatus = -1
      PRINT *, "Error: Input and output vertical levels do not match."
      RETURN
    END IF

    ! Implement 1D vertical interpolation logic here
    ALLOCATE(fint(levels_in))
    IF (PRESENT(logOpt)) THEN
      fint = DLOG(fin)
    ELSE
      fint = fin
    END IF

    ! Note: only the input coordinate and function need to be checked. No requirement for the output coordinate.
    ! Check if the coodinate decreases or increases:
    IF (vin(1) < vin(levels_in)) THEN
      ! Increasing vertical coordinate: height type of coordinate
      vint = vin(1:levels_in)
      fint = fint(1:levels_in)
    ELSE
      ! Decreasing vertical coordinate: pressure type of coordinate
      ! Reverse the input arrays for interpolation
      vint = vin(levels_in:1:-1)
      fint = fint(levels_in)
    END IF

    ! Only need to consider a coordinate increasing order only:
    fout = 0.0D0
    PRINT*,'verticalInterp1D: levels_in, levels_out:', levels_in, levels_out,UBOUND(fout, 1)
    DO i = 1, levels_out
      PRINT*,'v int 1:', i
      ! Simple linear interpolation for demonstration
      IF (vout(i) < vin(1)) THEN
        im1 = 1; ip1 = 2
      ELSE IF (vout(i) >= vin(levels_in)) THEN
        im1 = levels_in - 1; ip1 = levels_in
      ELSE
        DO j = 1, levels_in - 1
          IF (vout(i) >= vin(j) .AND. vout(i) < vin(j + 1)) THEN
            ! Linear interpolation
            im1 = j; ip1 = j + 1
            EXIT
          END IF
        END DO
      END IF
      PRINT*,'v int 2:', i, im1, ip1, vout(i), vint(im1), vint(ip1), fint(im1), fint(ip1)
      ! Extrapolate in second order:
      fout(i) = fint(im1)*(vint(ip1) - vout(i)) / (vint(ip1) - vint(im1)) + &
                fint(ip1)*(vint(im1) - vout(i)) / (vint(im1) - vint(ip1))
    END DO
    ! Clean up
    DEALLOCATE(fint)
  END SUBROUTINE verticalInterp1D

  SUBROUTINE verticalInterp2D(this, fin, vin, fout, vout, istatus, logOpt)
    CLASS(verticalInterp_t), INTENT(IN) :: this
    REAL(c_double), DIMENSION(:,:), INTENT(IN) :: fin
    REAL(c_double), DIMENSION(:,:), INTENT(IN) :: vin
    REAL(c_double), DIMENSION(:,:), INTENT(IN) :: vout
    REAL(c_double), DIMENSION(:,:), INTENT(OUT) :: fout
    INTEGER(c_int), INTENT(OUT) :: istatus
    INTEGER(c_int), INTENT(IN), OPTIONAL :: logOpt

    ! Local variables:
    INTEGER(c_int) :: levels_in, levels_out
    INTEGER(c_int) :: i, j
    REAL(c_double), ALLOCATABLE :: fint(:,:)

    levels_in = SIZE(fin, 1)
    levels_out = SIZE(fout, 1)
    istatus = 0
    IF (SIZE(vin, 1) /= levels_in .OR. SIZE(vout, 1) /= levels_out) THEN
      istatus = -1
      PRINT *, "Error: Input and output vertical levels do not match."
      RETURN
    END IF

    ! Implement 2D vertical interpolation logic here
    ALLOCATE(fint(levels_in, SIZE(fin, 2)))
    IF (PRESENT(logOpt)) THEN
      fint = DLOG(fin)
    ELSE
      fint = fin
    END IF
    fout = fint(1:levels_in, :)
    DEALLOCATE(fint)
  END SUBROUTINE verticalInterp2D

  SUBROUTINE verticalInterp3D(this, fin, vin, fout, vout, istatus, logOpt)
    CLASS(verticalInterp_t), INTENT(IN) :: this
    REAL(c_double), DIMENSION(:,:,:), INTENT(IN) :: fin
    REAL(c_double), DIMENSION(:,:,:), INTENT(IN) :: vin
    REAL(c_double), DIMENSION(:,:,:), INTENT(IN) :: vout
    REAL(c_double), DIMENSION(:,:,:), INTENT(OUT) :: fout
    INTEGER(c_int), INTENT(OUT) :: istatus
    INTEGER(c_int), INTENT(IN), OPTIONAL :: logOpt

    ! Local variables:
    INTEGER(c_int) :: levels_in, levels_out
    INTEGER(c_int) :: i, j, k
    REAL(c_double), ALLOCATABLE :: fint(:,:,:)

    levels_in = SIZE(fin, 1)
    levels_out = SIZE(fout, 1)
    istatus = 0
    IF (SIZE(vin, 1) /= levels_in .OR. SIZE(vout, 1) /= levels_out) THEN
      istatus = -1
      PRINT *, "Error: Input and output vertical levels do not match."
      RETURN
    END IF

    ! Implement 3D vertical interpolation logic here
    ALLOCATE(fint(levels_in, SIZE(fin, 2), SIZE(fin, 3)))
    IF (PRESENT(logOpt)) THEN
      fint = DLOG(fin)
    ELSE
      fint = fin
    END IF
    fout = fint(1:levels_in, :, :)
    DEALLOCATE(fint)
  END SUBROUTINE verticalInterp3D
END MODULE verticalInterp_m