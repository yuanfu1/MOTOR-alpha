!!-----------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.mgGen.edgeSetting_m
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring
!                   : Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/12/24, @GBA-MWF, Lafayette, CO
!!-----------------------------------------------------------------------------------------------------

!> @brief
!! This routine defines a mgGenLatLon edge's interpolation scheme of function, normal and tangential 
!! derivatives. It gives the interpolation stencil, coefficients.
!> @note We define the normal directions only for those corner cells of a latlon grid, leave tangential
!!  derivative interpolation coefficients to zeros. 
!! Input:
!!    ij:               2 elements of the i and j indices of a latlon grid
!!    nm:               2 elements of the maximum cell number in lat and lon directions
!!    numStencial:      Number of stencil cells for these interpolation schemes
!!    info:             edge orientation:
!!                                        1 - bottom
!!                                        2 - right
!!                                        3 - top
!!                                        4 - left
!! Older definition
!!                                        11 - left   edge of ij cell in interior
!!                                        12 - right  edge of ij cell in interior
!!                                        13 - left   edge of ij cell at left boundary
!!                                        14 - right  edge of ij cell at right boundary
!!                                        21 - bottom edge of ij cell in interior
!!                                        22 - top    edge of ij cell in interior
!!                                        23 - bottom edge of ij cell at bottom boundary
!!                                        24 - top    edge of ij cell at top boundary
!! Output:
!!    idxStencil:       indices of these stencils
!!    coeff_func:       interpolation coefficients of the function value at the edge center
!!    coeff_norm:       interpolation coefficients of the normal derivative at the edge center
!!    coeff_tang:       interpolation coefficients of the tangential derivative at the edge center
!> @warning
!> @attention
SUBROUTINE edgeSetting_s(ij,nm,ll,dll,numStencil,info,idxStencil, &
  coeff_func,coeff_norm,coeff_tang,istatus)
  USE kinds_m, ONLY: i_kind, r_kind
  USE parameters_m, ONLY: EarthRadius

  IMPLICIT NONE

  INTEGER(i_kind), INTENT(IN) :: ij(2),nm(2)    ! latlon cell index
  INTEGER(i_kind), INTENT(IN) :: numStencil     ! Number of stencils, currently is set to 7 for boundary interpolation
  INTEGER(i_kind), INTENT(IN) :: info
  REAL(r_kind), INTENT(IN) :: ll(2),dll(2)      ! cell LL and dlatlon at the equitor
  INTEGER(i_kind), INTENT(OUT) :: idxStencil(numStencil)
  INTEGER(i_kind), INTENT(OUT) :: istatus
  REAL(r_kind), INTENT(OUT) :: coeff_func(numStencil),coeff_norm(numStencil),coeff_tang(numStencil)

  ! Local variables:
  INTEGER(i_kind) :: inc(3),idn(2),maxIdx,iorder(2)
  REAL(r_kind) :: extrapolation(3,0:1),oneSideIntpln(3) ! Assume ds is positive
  DATA extrapolation/1.5D0, -0.5, 0.0D0, -2.0D0, 3.0D0, -1.0D0/
  DATA oneSideIntpln/-1.5D0, 2.0D0, -0.5D0/

  ! Check if the number of stencil matches the scheme in this routine:
  istatus = 0     ! Default as successful run
  IF (numStencil .NE. 7) THEN
    WRITE(*,1) numStencil
1   FORMAT('Current latlon edge interpolation with boundaries requires 7 stencil but we got ',I2)
    STOP
  END IF
  maxIdx = nm(1)*nm(2)

  iorder = ij/nm*2

  ! Edge along latitude line (vertical):
  IF (info .EQ. 2 .OR. info .EQ. 4) THEN
    idn(1) = 1
    idn(2) = 2
    inc = 0           ! inc(2) is used for edges on the corner
    SELECT CASE (info)
    CASE (2)      ! right edge
      inc(1) =  1
      IF (ij(1) .EQ. nm(1)) inc(1) = -1 ! right boundary
    CASE (4)     ! left edge
      inc(1) = -1
      IF (ij(1) .EQ. 1) inc(1) = 1      ! left boundary
    CASE DEFAULT
      WRITE(*,2) info
2     FORMAT('edgeSetting_s: Invalid latitude edge orientation: ',I2)
    END SELECT
    ! Hit the bottom or the top boundaries, shift the stencil:
    IF (ij(2) .EQ. 1) THEN
      inc(2) = 1; inc(3) = 0
    END IF
    IF (ij(2) .EQ. nm(2)) THEN
      inc(2) = 0; inc(3) = -1
    END IF
    ! stencil:
    ! +---+---+---+         +---+---+---+
    ! | 3 | 6 |   |         |   | 6 | 3 |
    ! +---+---+---+         +---+---+---+
    ! | 2 | 5 | 7 |   or    | 7 | 5 | 2 |
    ! +---+---+---+         +---+---+---+
    ! | 1 | 4 |   |         |   | 4 | 1 |
    ! +---+---+---+         +---+---+---+
    idxStencil(1) = ij(1) + nm(1)*(ij(2)-2 + inc(2) + inc(3))
    idxStencil(2) = ij(1) + nm(1)*(ij(2)-1 + inc(2) + inc(3))     ! current cell
    idxStencil(3) = ij(1) + nm(1)*(ij(2)   + inc(2) + inc(3))
    idxStencil(4) = ij(1) + inc(1) + nm(1)*(ij(2)-2 + inc(2) + inc(3))
    idxStencil(5) = ij(1) + inc(1) + nm(1)*(ij(2)-1 + inc(2) + inc(3))
    idxStencil(6) = ij(1) + inc(1) + nm(1)*(ij(2)   + inc(2) + inc(3))
    ! Default index 7 is unused except the boundary edge only:
    idxStencil(7) = ij(1) + nm(1)*(ij(2)-1)     ! Default to use the current cell

    ! Initialize the coefficients:
    coeff_func = 0.0D0; coeff_norm = 0.0D0; coeff_tang = 0.0D0
    IF ((ij(1) .EQ. 1 .AND. info .EQ. 4) .OR. &
        (ij(1) .EQ. nm(1) .AND. info .EQ. 2)) THEN
      ! extrapolation needs the 7 cell:
      idxStencil(7) = ij(1) + 2*inc(1) + nm(1)*(ij(2)-1)
      ! function: use extrapolation 0 index array
      coeff_func(2-inc(2)-inc(3)) = extrapolation(1,0)
      coeff_func(5-inc(2)-inc(3)) = extrapolation(2,0)
      ! normal derivative:use extrapolation 1 index array
      coeff_norm(2-inc(2)-inc(3)) = -DBLE((1-iorder(1))*(info-3))*extrapolation(1,1)/dll(idn(2))
      coeff_norm(5-inc(2)-inc(3)) = -DBLE((1-iorder(1))*(info-3))*extrapolation(2,1)/dll(idn(2))
      coeff_norm(7) =               -DBLE((1-iorder(1))*(info-3))*extrapolation(3,1)/dll(idn(2))
    ELSE
      coeff_func(5-inc(2)-inc(3)) =  0.5D0
      coeff_func(2-inc(2)-inc(3)) =  0.5D0
      coeff_norm(5-inc(2)-inc(3)) =  1.0D0/dll(idn(2))
      coeff_norm(2-inc(2)-inc(3)) = -1.0D0/dll(idn(2))
    END IF
    ! tangential derivative: except the corner cells where no tangential derivatives are defined
    ! at corner cells, two normal derivatives are perpendicular, representing two independent vectors
    IF (ij(2) .GT. 1 .AND. ij(2) .LT. nm(2)) THEN
      ! Interior edge:
      coeff_tang(1) =  0.25D0*DBLE(info-3)/dll(idn(1))
      coeff_tang(3) = -0.25D0*DBLE(info-3)/dll(idn(1))
      coeff_tang(4) =  0.25D0*DBLE(info-3)/dll(idn(1))
      coeff_tang(6) = -0.25D0*DBLE(info-3)/dll(idn(1))
    ELSE IF ((ij(1) .NE. 1 .OR. info .EQ. 2) .AND. &
             (ij(1) .NE. nm(1) .OR. info .EQ. 4)) THEN
      ! Boundary: f'(1) = 2f'(2)-f'(3)
      coeff_tang(1+iorder(2):3-iorder(2):1-iorder(2)) = 0.5D0*(1.0D0-DBLE(iorder(2)))*DBLE(inc(1))*oneSideIntpln/dll(idn(1))
      coeff_tang(4+iorder(2):6-iorder(2):1-iorder(2)) = 0.5D0*(1.0D0-DBLE(iorder(2)))*DBLE(inc(1))*oneSideIntpln/dll(idn(1))
    END IF
  ELSE IF (info .EQ. 1 .OR. info .EQ. 3) THEN
    idn(1) = 1
    idn(2) = 2
    inc = 0
    SELECT CASE (info)
    CASE (1)     ! bottom edge
      inc(1) = -1
      IF (ij(2) .EQ. 1) inc(1) = 1; 
    CASE (3)     ! top edge
      inc(1) =  1
      IF (ij(2) .EQ. nm(2)) inc(1) =  -1
    CASE DEFAULT
      WRITE(*,3) info
3     FORMAT('edgeSetting_s: Invalid longitude edge orientation: ',I2)
    END SELECT
    ! Setting the corners to an arbitrary interior with zero coefficient:
    IF (ij(1) .EQ. 1) THEN
      inc(2) = 1; inc(3) = 0
    END IF
    IF (ij(1) .EQ. nm(1)) THEN
      inc(2) = 0; inc(3) = -1
    END IF
    ! stencil:
    ! +---+---+---+         +---+---+---+
    ! |   | 7 |   |         | 1 | 2 | 3 |
    ! +---+---+---+         +---+---+---+
    ! | 4 | 5 | 6 |   or    | 4 | 5 | 6 |
    ! +---+---+---+         +---+---+---+
    ! | 1 | 2 | 3 |         |   | 7 |   |
    ! +---+---+---+         +---+---+---+
    
    idxStencil(1) = ij(1) - 1 + inc(2) + inc(3) + nm(1)*(ij(2)-1)
    idxStencil(2) = ij(1)     + inc(2) + inc(3) + nm(1)*(ij(2)-1)                   ! Current cell
    idxStencil(3) = ij(1) + 1 + inc(2) + inc(3) + nm(1)*(ij(2)-1)
    idxStencil(4) = ij(1) - 1 + inc(2) + inc(3) + nm(1)*(ij(2)-1 + inc(1))
    idxStencil(5) = ij(1)     + inc(2) + inc(3) + nm(1)*(ij(2)-1 + inc(1))
    idxStencil(6) = ij(1) + 1 + inc(2) + inc(3) + nm(1)*(ij(2)-1 + inc(1))
    ! Default index 7 is unused except the boundary edge only:
    idxStencil(7) = ij(1) + nm(1)*(ij(2)-1)     ! Default to use the current cell

    ! Initialize the coefficients:
    coeff_func = 0.0D0; coeff_norm = 0.0D0; coeff_tang = 0.0D0
    IF ((info .EQ. 1 .AND. ij(2) .EQ. 1) .OR. &
        (info .EQ. 3 .AND. ij(2) .EQ. nm(2))) THEN
      ! extrapolation needs the 7 cell:
      idxStencil(7) = ij(1) + nm(1)*(ij(2)-1 + 2*inc(1))
      ! function: use extrapolation 0 index array
      coeff_func(2-inc(2)-inc(3)) = extrapolation(1,0)
      coeff_func(5-inc(2)-inc(3)) = extrapolation(2,0)
      ! normal derivative: use extrapolation 1 index array
      coeff_norm(2-inc(2)-inc(3)) = (1.0D0-iorder(2))*DBLE(info-2)*extrapolation(1,1)/dll(idn(1))
      coeff_norm(5-inc(2)-inc(3)) = (1.0D0-iorder(2))*DBLE(info-2)*extrapolation(2,1)/dll(idn(1))
      coeff_norm(7) =               (1.0D0-iorder(2))*DBLE(info-2)*extrapolation(3,1)/dll(idn(1))
    ELSE
      coeff_func(5-inc(2)-inc(3)) =  0.5D0
      coeff_func(2-inc(2)-inc(3)) =  0.5D0
      coeff_norm(5-inc(2)-inc(3)) =  1.0D0/dll(idn(1))
      coeff_norm(2-inc(2)-inc(3)) = -1.0D0/dll(idn(1))
    END IF
    ! tangential derivative: except the corner cells where no tangential derivatives are defined
    ! at corner cells, two normal derivatives are perpendicular, representing two independent vectors
    IF (ij(1) .GT. 1 .AND. ij(1) .LT. nm(1)) THEN
      IF ((ij(2) .NE. 1 .AND. ij(2) .NE. nm(2)) .OR. &
        (info .EQ. 1 .AND. ij(2) .EQ. nm(2)) .OR. &
        (info .EQ. 3 .AND. ij(2) .EQ. 1)) THEN
        coeff_tang(1) =  0.25D0*DBLE(info-2)/dll(idn(2))
        coeff_tang(3) = -0.25D0*DBLE(info-2)/dll(idn(2))
        coeff_tang(4) =  0.25D0*DBLE(info-2)/dll(idn(2))
        coeff_tang(6) = -0.25D0*DBLE(info-2)/dll(idn(2))
      ELSE
        coeff_tang(1) =  0.75D0*DBLE(info-2)/dll(idn(2))
        coeff_tang(3) = -0.75D0*DBLE(info-2)/dll(idn(2))
        coeff_tang(4) = -0.25D0*DBLE(info-2)/dll(idn(2))
        coeff_tang(6) =  0.25D0*DBLE(info-2)/dll(idn(2))
      END IF
    ELSE IF ((ij(2) .NE. 1 .AND. info .NE. 3) .OR. &
             (ij(2) .NE. nm(2) .AND. info .NE. 1)) THEN
      coeff_tang(1+iorder(1):3-iorder(1):1-iorder(1)) = -0.5D0*(1.0D0-DBLE(iorder(1)))*DBLE(info-2)*oneSideIntpln(:)/dll(idn(2))
      coeff_tang(4+iorder(1):6-iorder(1):1-iorder(1)) = -0.5D0*(1.0D0-DBLE(iorder(1)))*DBLE(info-2)*oneSideIntpln(:)/dll(idn(2))
    ELSE IF ((ij(2) .EQ. 1 .AND. info .EQ. 1) .OR. &
             (ij(2) .EQ. nm(2) .AND. info .EQ. 3)) THEN
      coeff_tang(1+iorder(1):3-iorder(1):1-iorder(1)) = -1.5D0*(1.0D0-DBLE(iorder(1)))*DBLE(info-2)*oneSideIntpln(:)/dll(idn(2))
      coeff_tang(4+iorder(1):6-iorder(1):1-iorder(1)) =  0.5D0*(1.0D0-DBLE(iorder(1)))*DBLE(info-2)*oneSideIntpln(:)/dll(idn(2))
    END IF
  END IF
END SUBROUTINE edgeSetting_s