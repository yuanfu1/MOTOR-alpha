!!--------------------------------------------------------------------------------------------------
! PROJECT           : Grid generation
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : Beta 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2020/06/11, @GBA-MWF, Shenzhen
!   Modified by Zilong Qin (zilong.qin@gmail.com), 2020/11/23, @GBA-MWF, Shenzhen
!   Modified by Zilong Qin (zilong.qin@gmail.com), 2021/4/19, @GBA-MWF, Shenzhen, update the cell_stcl for boundary cells
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2025/02/25 for revising cell_adnb to match edge setting and adding it to boundary cells
!!--------------------------------------------------------------------------------------------------

!!===================================================================
!> @brief
!! # Latlon grid generation routine collection
!!
!!  *This file defines a set of routines implementing latlon grid.*
!!
!!  ## Important note: Latlon grid index is stored as lat first and
!!   lon second. However the grid cell is ordered in longitude and
!!   latitude, which meets the convention of domain north pointing
!!   up and east pointing to right*
!!
!! ## Routine set:
!!
!! ### getDomain_latlon
!!  *The routine reads domain parameters*
!!
!! @author Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
!!===================================================================
SUBROUTINE getDomain_latlon(configFile, nv, nc, ng, dm)
  USE kinds_m, ONLY: i_kind, r_kind
  USE YAMLRead_m

  IMPLICIT NONE

  CHARACTER(LEN=1024), INTENT(IN) :: configFile
  INTEGER(i_kind)                       :: rc

  INTEGER(i_kind), INTENT(OUT) :: nv, nc, ng(2)   ! Numbers of vertices, cells and grids
  REAL(r_kind), INTENT(OUT) :: dm(2, 2)            ! Latlon ranges in degrees

  ! Domain namelist variables:
  CHARACTER(LEN=256) :: path2static
  INTEGER(i_kind), ALLOCATABLE :: num_grid(:)         ! Numbers of grid points in latlon directions
  REAL(r_kind), ALLOCATABLE :: domain_latlon(:)         ! The domain range in degrees from the origin

  INTEGER(i_kind) :: istatus

  !NAMELIST /latlon_grid/ num_grid, domain_latlon

  ! Check whether file exists.
  INQUIRE (file=TRIM(configFile), iostat=rc)

  ! Read the namelist file
  !IF (rc .eq. 0) THEN
  !  OPEN (10, FILE=TRIM(configFile))
  !  READ (10, NML=latlon_grid)
  !  CLOSE (10)
  !ELSE
  !  ERROR STOP "Namelist file of *getDomain_latlon* does not exist!"
  !END IF
  IF (rc .EQ. 0) THEN
    istatus = yaml_get_var(TRIM(configFile), 'latlon_grid', 'num_grid', num_grid)
    istatus = yaml_get_var(TRIM(configFile), 'latlon_grid', 'domain_latlon', domain_latlon)
    PRINT *, "domain_latlon: ", domain_latlon
  ELSE
    ERROR STOP "Namelist file of *getDomain_latlon* does not exist!"
  END IF

  ! Namelist gives the actual grid points in each direction
  ! We need to add fictitious points on each direction:
  ! Number of vertices: (num_grid(1)+2)x(num_grid(2)+2)
  ! Number of cells: (num_grid(1)+1)x(num_grid(2)+1)
  nv = (num_grid(1) + 2) * (num_grid(2) + 2); nc = (num_grid(1) + 1) * (num_grid(2) + 1)
  ng = num_grid + 2
  dm = RESHAPE(domain_latlon, (/2, 2/))

END SUBROUTINE getDomain_latlon

!>
!!===================================================================
!!  This routine provides basic grid geometry location only, normally
!!  this part runs on a serial code.
!!
!!  Output:
!!    num_vertices: number of cell's vertices,  optional now
!!    num_gridcell: number of cells             optional now
!!    grid_latlons: the latitides and longitudes of these vertices
!!    cell_neighbr: the cell's cell and vertex neighbor information
!!
!!    cell_neighbr data structure: PRIVATE
!!      TYPE Cell
!!        INTEGER :: vextex
!!        INTEGER :: adjnbs
!!      END TYPE Cell
!!
!!  Assumption: Latlon grid is always symmetric to (0,0),dm(2) gives
!!              the ranges of latitude and longitude
!!
!!  Modified by Yuanfu Xie on 2025-02-25 for revising the cell_adnb
!!  matching with the edge order set in the latlon_precal where the
!!  edge stencil is set from bottom-right-top-left. In addition, he
!!  added the cell_adnb to the boundary cells as well.
!!===================================================================
!
SUBROUTINE latlon_geoGrid(ng, dm, para, geoq)
  USE kinds_m, ONLY: i_kind, r_kind
  ! USE cell_m, ONLY : cellLatlon_t
  USE gridStruct_m, ONLY: gridParams_t, gridGeoQty_t
  USE parameters_m, ONLY: degree2radian

  IMPLICIT NONE

  INTEGER(i_kind), INTENT(IN) :: ng(2)! Numbers of vertices and cells
  REAL(r_kind), INTENT(IN) :: dm(2, 2)         ! Latlon ranges in degrees

  TYPE(gridParams_t), INTENT(IN) :: para
  ! Structed type cannot be intented as OUT, e.g., geoq;
  ! otherwise its dimension will be unspecified.
  TYPE(gridGeoQty_t) :: geoq

  ! Local variables:

  !REAL(r_kind), ALLOCATABLE :: vertex_ll(:,:)
  !TYPE(cellLatlon_t), ALLOCATABLE :: grid_mesh(:)

  CHARACTER(LEN=256) :: path2static
  INTEGER(i_kind) :: i, j, nv, nc, im1,ip1,jm1,jp1
  REAL(r_kind) :: xy(2), dxy(2)

  ! Vertex and center latlon:
  ! Note: we order the vertices and centers in longitude first
  ! the latlon is ordered lat first and long second:
  !
  ! Domain extended with dxy for those fictitious points:
  nv = 0
  dxy = (dm(2, :) - dm(1, :)) / DBLE(ng - 3) ! actual interior cells are ng-3
  ! PRINT *,'dxy: ', dxy/degree2radian, dm(2, :)/degree2radian, dm(1, :)/degree2radian, DBLE(ng - 3)
  DO i = 1, ng(1)
    xy(1) = dm(1, 1) + DBLE(i - 2) * dxy(1)
    DO j = 1, ng(2)
      xy(2) = dm(1, 2) + DBLE(j - 2) * dxy(2)
      nv = nv + 1
      geoq%vrtx_lalo(1:2, nv) = xy
    END DO
  END DO
  IF (nv .NE. para%num_vrtx) THEN
    PRINT *, 'Throw this error out later: num_vrtx is not right!', &
      nv, para%num_vrtx
    STOP
  END IF

  ! Cell's neighbor information:
  nc = 0
  DO i = 1, ng(1) - 1
    DO j = 1, ng(2) - 1
      ! Cells are between vertices:
      nc = nc + 1
      ! Vertex information: right hand rule counter-clock
      ! left-bottom, right-bottom, right-top, right-left
      ! Note: the order of indices 3 and 4:
      geoq%cell_vrtx(1, nc) = j + ng(2) * (i - 1)
      geoq%cell_vrtx(2, nc) = j + 1 + ng(2) * (i - 1)
      geoq%cell_vrtx(4, nc) = j + ng(2) * (i)
      geoq%cell_vrtx(3, nc) = j + 1 + ng(2) * (i)

      ! Wrap the increments of i and j:
      im1 = MAX(i-1,0); ip1 = MOD(i+1,ng(1))
      jm1 = MAX(j-1,0); jp1 = MOD(j+1,ng(2))

      ! ! Cell information: Yuanfu Xie re-activated the following to define cell_adnb on interior only 2025-03-31
      geoq%cell_adnb(:, nc) = -1 ! No neighbors
      ! ! right hand rule counter-clock: left-bottom-right-top
      IF (i .GT. 1 .AND. i .LT. ng(1) - 1 .AND. &
          j .GT. 1 .AND. j .LT. ng(2) - 1) THEN
        geoq%cell_adnb(1, nc) = j - 1 + (ng(2) - 1) * (i - 1)
        geoq%cell_adnb(2, nc) = j + (ng(2) - 1) * (i - 2)
        geoq%cell_adnb(3, nc) = j + 1 + (ng(2) - 1) * (i - 1)
        geoq%cell_adnb(4, nc) = j + (ng(2) - 1) * (i)
      END IF
    END DO
  END DO
  IF (nc .NE. para%num_cell) THEN
    PRINT *, 'Throw this error out later: num_cell is not right!', &
      nc, para%num_cell
    STOP
  END IF

END SUBROUTINE latlon_geoGrid

!>
!!===================================================================
!!  This routine pre-calculates all necessary quantities of a finite
!!  volume scheme. Based on the analysis:
!!  LatLonModel.docs under
!!  projects/RegionalModel
!!  the centroid center is approximated by the middle latitude center
!!  in a second order accuracy. Therefore, we implemented the finite
!!  volume scheme of Poisson equation using middle latitude center,
!!  just like the regular rectangle grid.
!!
!!  @Author Yuanfu Xie
!!   History:
!!  created by Yuanfu Xie on 06-25-2020
!!  Modified by Yuanfu Xie on 07-23-2020, by considering the fictitious
!!  grid cells outside of the physical domain, so that the grid info
!!  will relate to the interior points.
!!
!! Below are modified by Zilong when coding MOTOR-DA
!! The edge_stcl, coef_func, coef_norm, coef_tgnt values on the boundary
!! are set to zero.
!!===================================================================
!
SUBROUTINE latlon_precals(nmgd, rsll, para, geoq)
  USE omp_lib
  USE kinds_m, ONLY: i_kind, r_kind
  ! USE cell_m, ONLY : cellLatlon_t
  USE gridStruct_m, ONLY: gridParams_t, gridGeoQty_t
  USE parameters_m, ONLY: degree2radian, EarthRadius

  IMPLICIT NONE

  INTEGER(i_kind) :: nmgd(2), ldim, rdim
  REAL(r_kind), INTENT(IN) :: rsll(2)
  TYPE(gridParams_t) :: para
  TYPE(gridGeoQty_t) :: geoq

  ! Local variables:
  CHARACTER(LEN=256) :: path2static
  INTEGER(i_kind) :: i, j, nc, nbeg(2), nend(2), ne, ie, nm(2) ! Yuanfu Xie added nm for calling edgeSetting putting lon in first 2025-01-08
  REAL(r_kind) :: rsln(2)
  INTEGER(i_kind) :: cell_stcl_width, m, n

  ! Use of a new routine edgeSetting:
  INTEGER(i_kind) :: ij(2),info,istatus

  ! Yuanfu Xie added the Gauss-Legendre quadrature weights:
  geoq%coef_gl(1:para%numQuadPerEdge) = 1.0D0

  ! Cell center latlon:
  !$omp PARALLEL DO
  DO i = 1, para%num_cell
    ! Assuming the cell vertices are arranged from
    ! left-bottom,right-bottom, right-top, left-top:
    geoq%cell_cntr(1, i) = 0.5D0 * ( &
                           geoq%vrtx_lalo(1, geoq%cell_vrtx(1, i)) + &
                           geoq%vrtx_lalo(1, geoq%cell_vrtx(4, i)))
    geoq%cell_cntr(2, i) = 0.5D0 * ( &
                           geoq%vrtx_lalo(2, geoq%cell_vrtx(1, i)) + &
                           geoq%vrtx_lalo(2, geoq%cell_vrtx(2, i)))

    geoq%cell_area(i) = EarthRadius**2 * &
                        ABS(geoq%vrtx_lalo(2, geoq%cell_vrtx(2, i)) - &
                            geoq%vrtx_lalo(2, geoq%cell_vrtx(1, i))) * &
                        ABS(DSIN(geoq%vrtx_lalo(1, geoq%cell_vrtx(4, i))) - &
                            DSIN(geoq%vrtx_lalo(1, geoq%cell_vrtx(1, i))))

    ! edges ordered as bottom, right, top and left:
    geoq%edge_lnth(1, i) = ABS(geoq%vrtx_lalo(2, geoq%cell_vrtx(2, i)) - &
                               geoq%vrtx_lalo(2, geoq%cell_vrtx(1, i))) * &
                           ABS(DCOS(geoq%vrtx_lalo(1, geoq%cell_vrtx(1, i)))) * EarthRadius
    geoq%edge_lnth(2, i) = ABS(geoq%vrtx_lalo(1, geoq%cell_vrtx(3, i)) - &
                               geoq%vrtx_lalo(1, geoq%cell_vrtx(2, i))) * EarthRadius
    geoq%edge_lnth(3, i) = ABS(geoq%vrtx_lalo(2, geoq%cell_vrtx(4, i)) - &
                               geoq%vrtx_lalo(2, geoq%cell_vrtx(3, i))) * &
                           ABS(DCOS(geoq%vrtx_lalo(1, geoq%cell_vrtx(3, i)))) * EarthRadius
    geoq%edge_lnth(4, i) = geoq%edge_lnth(2, i)
  END DO
  !$omp END PARALLEL DO

  geoq%cell_side = 0 ! Unspecified
  geoq%cell_stcl = 0 ! Default as boundary

  nbeg = 1
  nend = nmgd

  ! Do not count the fictitious points:
  ldim = nmgd(2) - 3
  rdim = nmgd(1) - 3
  nc = 0
  !$OMP PARALLEL DO PRIVATE(i,j,nc)
  DO i = nbeg(1), nend(1) - 1
    DO j = nbeg(2), nend(2) - 1 ! Longitude first
      !nc = nc + 1
      nc = j + (i - 1) * (nend(2) - nbeg(2))

      ! Multigrid info:
      ! By adding 0.2 to shift the parent cells to the right by 1
      geoq%prnt_indx(nc) = NINT((j - 1) / 2.0D0 + 0.2) + 1 + (ldim / 2 + 2) * NINT((i - 1) / 2.0D0 + 0.2)

      ! See the doc 2020-7-10 latlon model:
      ! geoq%chld_indx(1, nc) = (j - 1)*2 + ((nmgd(2) - 3)*2 + 2)*((i - 1)*2 - 1)
      ! IF ((j - 1)*2 .GT. (nmgd(2) - 3)*2 + 2 .OR. (j - 1)*2 .EQ. 0) geoq%chld_indx(1, nc) = 0
      ! IF ((i - 1)*2 - 1 .GE. (nmgd(1) - 3)*2 + 2 .OR. (i - 1)*2 - 1 .LT. 0) geoq%chld_indx(1, nc) = 0
      ! geoq%chld_indx(2, nc) = (j - 1)*2 + 1 + ((nmgd(2) - 3)*2 + 2)*((i - 1)*2 - 1)
      ! IF ((j - 1)*2 + 1 .GT. (nmgd(2) - 3)*2 + 2 .OR. (j - 1)*2 + 1 .EQ. 0) geoq%chld_indx(2, nc) = 0
      ! IF ((i - 1)*2 - 1 .GE. (nmgd(1) - 3)*2 + 2 .OR. (i - 1)*2 - 1 .LT. 0) geoq%chld_indx(2, nc) = 0
      ! geoq%chld_indx(3, nc) = (j - 1)*2 + ((nmgd(2) - 3)*2 + 2)*((i - 1)*2)
      ! IF ((j - 1)*2 .GT. (nmgd(2) - 3)*2 + 2 .OR. (j - 1)*2 .EQ. 0) geoq%chld_indx(3, nc) = 0
      ! IF ((i - 1)*2 .GE. (nmgd(1) - 3)*2 + 2 .OR. (i - 1)*2 .LT. 0) geoq%chld_indx(3, nc) = 0
      ! geoq%chld_indx(4, nc) = (j - 1)*2 + 1 + ((nmgd(2) - 3)*2 + 2)*((i - 1)*2)
      ! IF ((j - 1)*2 + 1 .GT. (nmgd(2) - 3)*2 + 2 .OR. (j - 1)*2 + 1 .EQ. 0) geoq%chld_indx(4, nc) = 0
      ! IF ((i - 1)*2 .GE. (nmgd(1) - 3)*2 + 2 .OR. (i - 1)*2 .LT. 0) geoq%chld_indx(4, nc) = 0

      ! New version by Qin Lilan from SZHPC
      !geoq%chld_indx(1, nc) = (j - 1)*2 + ((nmgd(2) - 3)*2 + 2)*((i - 1)*2 - 1)
      geoq%chld_indx(1, nc) = (j - 1) * 2 + (ldim * 2 + 2) * ((i - 1) * 2 - 1)
      !IF ((j - 1)*2 .GT. (nmgd(2) - 3)*2 + 2 .OR. (j - 1)*2 .EQ. 0) geoq%chld_indx(1, nc) = 0
      !IF ((i - 1)*2 - 1 .GE. (nmgd(1) - 3)*2 + 2 .OR. (i - 1)*2 - 1 .LT. 0) geoq%chld_indx(1, nc) = 0
      IF ((j - 1) * 2 .GT. ldim * 2 + 2 .OR. (j - 1) * 2 .EQ. 0) geoq%chld_indx(1, nc) = 0
      IF ((i - 1) * 2 - 1 .GE. rdim * 2 + 2 .OR. (i - 1) * 2 - 1 .LT. 0) geoq%chld_indx(1, nc) = 0
      !geoq%chld_indx(2, nc) = (j - 1)*2 + 1 + ((nmgd(2) - 3)*2 + 2)*((i - 1)*2 - 1)
      geoq%chld_indx(2, nc) = (j - 1) * 2 + 1 + (ldim * 2 + 2) * ((i - 1) * 2 - 1)
      !IF ((j - 1)*2 + 1 .GT. (nmgd(2) - 3)*2 + 2 .OR. (j - 1)*2 + 1 .EQ. 0) geoq%chld_indx(2, nc) = 0
      !IF ((i - 1)*2 - 1 .GE. (nmgd(1) - 3)*2 + 2 .OR. (i - 1)*2 - 1 .LT. 0) geoq%chld_indx(2, nc) = 0
      IF ((j - 1) * 2 + 1 .GT. ldim * 2 + 2) geoq%chld_indx(2, nc) = 0
      IF ((i - 1) * 2 - 1 .GE. rdim * 2 + 2 .OR. (i - 1) * 2 - 1 .LT. 0) geoq%chld_indx(2, nc) = 0
      !geoq%chld_indx(3, nc) = (j - 1)*2 + ((nmgd(2) - 3)*2 + 2)*((i - 1)*2)
      geoq%chld_indx(3, nc) = (j - 1) * 2 + (ldim * 2 + 2) * ((i - 1) * 2)
      !IF ((j - 1)*2 .GT. (nmgd(2) - 3)*2 + 2 .OR. (j - 1)*2 .EQ. 0) geoq%chld_indx(3, nc) = 0
      !IF ((i - 1)*2 .GE. (nmgd(1) - 3)*2 + 2 .OR. (i - 1)*2 .LT. 0) geoq%chld_indx(3, nc) = 0
      IF ((j - 1) * 2 .GT. ldim * 2 + 2 .OR. (j - 1) * 2 .EQ. 0) geoq%chld_indx(3, nc) = 0
      IF ((i - 1) * 2 .GE. rdim * 2 + 2) geoq%chld_indx(3, nc) = 0
      !geoq%chld_indx(4, nc) = (j - 1)*2 + 1 + ((nmgd(2) - 3)*2 + 2)*((i - 1)*2)
      geoq%chld_indx(4, nc) = (j - 1) * 2 + 1 + (ldim * 2 + 2) * ((i - 1) * 2)
      !IF ((j - 1)*2 + 1 .GT. (nmgd(2) - 3)*2 + 2 .OR. (j - 1)*2 + 1 .EQ. 0) geoq%chld_indx(4, nc) = 0
      !IF ((i - 1)*2 .GE. (nmgd(1) - 3)*2 + 2 .OR. (i - 1)*2 .LT. 0) geoq%chld_indx(4, nc) = 0
      IF ((j - 1) * 2 + 1 .GT. ldim * 2 + 2) geoq%chld_indx(4, nc) = 0
      IF ((i - 1) * 2 .GE. rdim * 2 + 2) geoq%chld_indx(4, nc) = 0

      !WRITE(*,21) nc,j,i,geoq%chld_indx(:,nc)
21    FORMAT('Children at NC: ', I8, 2I5, ' Indices: ', 4I8)

      ! Check where the cell is located in its coarser grid:
      ! Break the 4 coarser grid of prolongation into 16 fine
      ! grid cells, mark the fine grid location 1-16:
      !
      !      +----+----+----+----+----+----+
      !      | 31 | 32 | 33 | 34 | 35 | 36 |      <- Fictitious points
      !      +----+----+----+----+----+----+
      !      | 25 | 26 | 27 | 28 | 29 | 30 |
      !      +----+----+----+----+----+----+
      !      | 19 | 20 | 21 | 22 | 23 | 24 |
      !      +----+----+----+----+----+----+
      !      | 13 | 14 | 15 | 16 | 17 | 18 |
      !      +----+----+----+----+----+----|
      !      | 07 | 08 | 09 | 10 | 11 | 12 |
      !      +----+----+----+----+----+----+
      !      | 01 | 02 | 03 | 04 | 05 | 06 |      <- Fictitious points
      !      +----+----+----+----+----+----+
      !
      ! +---------+---------+---------+---------+
      ! |         |         |         |         |
      ! |    13   |    14   |    15   |    16   | <- Fictitious points
      ! |         |         |         |         |
      ! +---------+---------+---------+---------+
      ! |         |         |         |         |
      ! |    09   |    10   |    11   |    12   |
      ! |         |         |         |         |
      ! +---------+---------+---------+---------+
      ! |         |         |         |         |
      ! |    05   |    06   |    07   |    08   |
      ! |         |         |         |         |
      ! +---------+---------+---------+---------+
      ! |         |         |         |         |
      ! |    01   |    02   |    03   |    04   | <- Fictitious points
      ! |         |         |         |         |
      ! +---------+---------+---------+---------+
      !
      ! Find the location at its parent's cell:
      ! left-bottom,right-bottom,left-top,right-top
      ! then find the 4 parents' indices and coefficient to interpolate
      IF (MOD(i, 2) .EQ. 0) THEN   ! located at its parent bottom
        IF (MOD(j, 2) .EQ. 0) THEN ! located at its parent left-bottom
          geoq%cTof_stcl(1, nc) = geoq%prnt_indx(nc) - 1
          geoq%cTof_stcl(2, nc) = geoq%prnt_indx(nc)
          geoq%cTof_stcl(3, nc) = geoq%prnt_indx(nc) - 1 - (ldim / 2 + 2)
          geoq%cTof_stcl(4, nc) = geoq%prnt_indx(nc) - (ldim / 2 + 2)
          geoq%cTof_coef(1, nc) = 3.0D0 / 16.0D0
          geoq%cTof_coef(2, nc) = 9.0D0 / 16.0D0
          geoq%cTof_coef(3, nc) = 1.0D0 / 16.0D0
          geoq%cTof_coef(4, nc) = 3.0D0 / 16.0D0

          ! geoq%cTof_coef(1, nc) = 0.1323D0 ! IDW coefficients by QZL
          ! geoq%cTof_coef(2, nc) = 0.6618D0
          ! geoq%cTof_coef(3, nc) = 0.0735D0
          ! geoq%cTof_coef(4, nc) = 0.1323D0
        ELSE                      ! located at its parent right-bottom
          geoq%cTof_stcl(1, nc) = geoq%prnt_indx(nc)
          geoq%cTof_stcl(2, nc) = geoq%prnt_indx(nc) + 1
          geoq%cTof_stcl(3, nc) = geoq%prnt_indx(nc) - (ldim / 2 + 2)
          geoq%cTof_stcl(4, nc) = geoq%prnt_indx(nc) + 1 - (ldim / 2 + 2)
          geoq%cTof_coef(1, nc) = 9.0D0 / 16.0D0
          geoq%cTof_coef(2, nc) = 3.0D0 / 16.0D0
          geoq%cTof_coef(3, nc) = 3.0D0 / 16.0D0
          geoq%cTof_coef(4, nc) = 1.0D0 / 16.0D0

          ! geoq%cTof_coef(1, nc) = 0.6618D0
          ! geoq%cTof_coef(2, nc) = 0.1323D0
          ! geoq%cTof_coef(3, nc) = 0.1323D0
          ! geoq%cTof_coef(4, nc) = 0.0735D0
        END IF
      ELSE                        ! located at its parent top
        IF (MOD(j, 2) .EQ. 0) THEN ! located at its parent left-top
          geoq%cTof_stcl(1, nc) = geoq%prnt_indx(nc) - 1
          geoq%cTof_stcl(2, nc) = geoq%prnt_indx(nc)
          geoq%cTof_stcl(3, nc) = geoq%prnt_indx(nc) - 1 + (ldim / 2 + 2)
          geoq%cTof_stcl(4, nc) = geoq%prnt_indx(nc) + (ldim / 2 + 2)
          geoq%cTof_coef(1, nc) = 3.0D0 / 16.0D0
          geoq%cTof_coef(2, nc) = 9.0D0 / 16.0D0
          geoq%cTof_coef(3, nc) = 1.0D0 / 16.0D0
          geoq%cTof_coef(4, nc) = 3.0D0 / 16.0D0

          ! geoq%cTof_coef(1, nc) = 0.1323D0
          ! geoq%cTof_coef(2, nc) = 0.6618D0
          ! geoq%cTof_coef(3, nc) = 0.0735D0
          ! geoq%cTof_coef(4, nc) = 0.1323D0
        ELSE                      ! located at its parent right-top
          geoq%cTof_stcl(1, nc) = geoq%prnt_indx(nc)
          geoq%cTof_stcl(2, nc) = geoq%prnt_indx(nc) + 1
          geoq%cTof_stcl(3, nc) = geoq%prnt_indx(nc) + (ldim / 2 + 2)
          geoq%cTof_stcl(4, nc) = geoq%prnt_indx(nc) + 1 + (ldim / 2 + 2)
          geoq%cTof_coef(1, nc) = 9.0D0 / 16.0D0
          geoq%cTof_coef(2, nc) = 3.0D0 / 16.0D0
          geoq%cTof_coef(3, nc) = 3.0D0 / 16.0D0
          geoq%cTof_coef(4, nc) = 1.0D0 / 16.0D0

          ! geoq%cTof_coef(1, nc) = 0.6618D0
          ! geoq%cTof_coef(2, nc) = 0.1323D0
          ! geoq%cTof_coef(3, nc) = 0.1323D0
          ! geoq%cTof_coef(4, nc) = 0.0735D0
          !WRITE(*,2) nc,j,i,geoq%cTof_stcl(:,nc)
2         FORMAT('At ', I8, 2I5, ' STCL: ', 4I5)
          !WRITE(*,1) geoq%cTof_coef(:,nc)
1         FORMAT('with coef: ', 4F14.6)
        END IF
      END IF
    END DO
  END DO
  !$OMP END PARALLEL DO

  ! Leading dimension:
  ldim = nmgd(2) - 1

  nc = 0
  !$OMP PARALLEL DO PRIVATE(i,j,nc,rsln,ij,nm)

  DO i = nbeg(1), nend(1) - 1
    DO j = nbeg(2), nend(2) - 1 ! Longitude first
      ! nc = nc + 1
      ! nc = j + (i - 1)*(nend(2) - nbeg(2))
      nc = j + ldim * (i - 1)

      ! Yuanfu Xie 2025-01-26: moved these two statements inside the OMP loop otherwise, their values got changed
      nm(1) = nend(2)-1
      nm(2) = nend(1)-1 ! Yuanfu Xie added nm to put lon first 2025-01-08

      ! Resolution:
      rsln(1) = rsll(1) * EarthRadius
      rsln(2) = rsll(2) * EarthRadius

      ! Presume the cell is a boundary cell, correct it later
      geoq%cell_stcl(:, nc) = 0
      geoq%edge_stcl(:, :, nc) = 0
      geoq%edge_stcl(7,:,nc) = nc     ! Default value of the stencil 7 Yuanfu Xie 2025-01-8
      geoq%coef_func(:, :, :, nc) = 0.0D0
      geoq%coef_norm(:, :, :, nc) = 0.0D0
      geoq%coef_tgnt(:, :, :, nc) = 0.0D0
      geoq%tgnt_vct2(:, :, nc) = 0.0D0
      geoq%norm_vctr(:, :, nc) = 0.0D0

      ! Rotation vectors for latlon direction are universal for all grid cells:
      geoq%norm_vctr(1, 1, nc) = 0.0D0
      geoq%norm_vctr(2, 1, nc) = -1.0D0
      geoq%norm_vct2(1, 1, nc) = 0.0D0
      geoq%norm_vct2(2, 1, nc) = -1.0D0
      geoq%tgnt_vctr(1, 1, nc) = 1.0D0
      geoq%tgnt_vctr(2, 1, nc) = 0.0D0
      geoq%tgnt_vct2(1, 1, nc) = 1.0D0
      geoq%tgnt_vct2(2, 1, nc) = 0.0D0

      geoq%norm_vctr(1, 2, nc) = 1.0D0
      geoq%norm_vctr(2, 2, nc) = 0.0D0
      geoq%norm_vct2(1, 2, nc) = 1.0D0
      geoq%norm_vct2(2, 2, nc) = 0.0D0
      geoq%tgnt_vctr(1, 2, nc) = 0.0D0
      geoq%tgnt_vctr(2, 2, nc) = 1.0D0
      geoq%tgnt_vct2(1, 2, nc) = 0.0D0
      geoq%tgnt_vct2(2, 2, nc) = 1.0D0

      geoq%norm_vctr(1, 3, nc) = 0.0D0
      geoq%norm_vctr(2, 3, nc) = 1.0D0
      geoq%norm_vct2(1, 3, nc) = 0.0D0
      geoq%norm_vct2(2, 3, nc) = 1.0D0
      geoq%tgnt_vctr(1, 3, nc) = -1.0D0
      geoq%tgnt_vctr(2, 3, nc) = 0.0D0
      geoq%tgnt_vct2(1, 3, nc) = -1.0D0
      geoq%tgnt_vct2(2, 3, nc) = 0.0D0

      geoq%norm_vctr(1, 4, nc) = -1.0D0
      geoq%norm_vctr(2, 4, nc) = 0.0D0
      geoq%norm_vct2(1, 4, nc) = -1.0D0
      geoq%norm_vct2(2, 4, nc) = 0.0D0
      geoq%tgnt_vctr(1, 4, nc) = 0.0D0
      geoq%tgnt_vctr(2, 4, nc) = -1.0D0
      geoq%tgnt_vct2(1, 4, nc) = 0.0D0
      geoq%tgnt_vct2(2, 4, nc) = -1.0D0

      ! Edge stencil runs through 6 stencils/edge/cell
      ! The following is designed to allow Poisson solver to load the coef matrix:
      ! ABS(edge_stcl) is the cell index where this cell is used for that Poisson eqn
      IF (i .EQ. 1) THEN
        ! Bottom boundary edge:
        ! See the development note: 20200709_LatLonModel.docx or later version under
        ! /Users/xieyuanfu/projects/RegionalModel/coupledBC4LatlonGrid
        !
        ! Yuanfu Xie made new comments on 2024-12-20 on the coupled solver.
        ! He also provided a modification of this grid generation with extrapolations
        ! on the edges of the boundary cells.
        IF (j .EQ. 1) THEN
          ! Following Zilong's corner setting, we add an additional edge information
          ! in edge 3!!!
          ! see page 6 on his revision of 20200729_LatLonModel.docx
          !
          ! Cell stencil:
          ! geoq%cell_stcl(1, nc) = j - 1 + ldim*(i - 2)   ! left bottom
          ! geoq%cell_stcl(2, nc) = j + ldim*(i - 2)       ! bottom
          ! geoq%cell_stcl(3, nc) = j + 1 + ldim*(i - 2)   ! right bottom
          ! geoq%cell_stcl(4, nc) = j - 1 + ldim*(i - 1)   ! left
          geoq%cell_stcl(5, nc) = j + ldim * (i - 1)       ! self
          geoq%cell_stcl(6, nc) = j + 1 + ldim * (i - 1)   ! right
          ! geoq%cell_stcl(7, nc) = j - 1 + ldim*(i)       ! left top
          geoq%cell_stcl(8, nc) = j + ldim * (i)           ! top
          geoq%cell_stcl(9, nc) = j + 1 + ldim * (i)       ! right top

          ! Use of the new routine of edgeSetting:
          ij(1) = j; ij(2) = i
          DO ie = 1, 4
            rsln(2) = rsll(2) * EarthRadius * &
                      DCOS(geoq%cell_cntr(1, nc) + 0.5D0*DBLE(MOD(ie-2,2))*rsll(1))
            CALL edgeSetting_s(ij,nm,geoq%cell_cntr(1:2, nc),rsln, &
              7,ie,geoq%edge_stcl(:, ie, nc),geoq%coef_func(:, 1, ie, nc), &
              geoq%coef_norm(:, 1, ie, nc),geoq%coef_tgnt(:, 1, ie, nc),istatus)
          END DO

          ! WRITE(*,112) geoq%edge_stcl(1:7, 1, nc),nm
    112   FORMAT('cell 1 1 stencil bottom edge New: ',9I6)
        ELSE IF (j .EQ. nmgd(2) - 1) THEN    ! Use nmgd for boundaries, nend may not
          ! Stream function stencil:

          ! geoq%cell_stcl(1, nc) = j - 1 + ldim*(i - 2)   ! left bottom
          ! geoq%cell_stcl(2, nc) = j + ldim*(i - 2)       ! bottom
          ! geoq%cell_stcl(3, nc) = j + 1 + ldim*(i - 2)   ! right bottom
          geoq%cell_stcl(4, nc) = j - 1 + ldim * (i - 1)   ! left
          geoq%cell_stcl(5, nc) = j + ldim * (i - 1)       ! self
          ! geoq%cell_stcl(6, nc) = j + 1 + ldim*(i - 1)   ! right
          geoq%cell_stcl(7, nc) = j - 1 + ldim * (i)       ! left top
          geoq%cell_stcl(8, nc) = j + ldim * (i)           ! top
          ! geoq%cell_stcl(9, nc) = j + 1 + ldim*(i)       ! right top

          ! Use the edgeSetting routine: Yuanfu Xie 2025-01-08
          ij(1) = j; ij(2) = i
          DO ie = 1, 4
            rsln(2) = rsll(2) * EarthRadius * &
                      DCOS(geoq%cell_cntr(1, nc) + 0.5D0*DBLE(MOD(ie-2,2))*rsll(1))
            CALL edgeSetting_s(ij,nm,geoq%cell_cntr(1, nc),rsln, &
                7,ie,geoq%edge_stcl(:, ie, nc),geoq%coef_func(:, 1, ie, nc), &
                geoq%coef_norm(:, 1, ie, nc),geoq%coef_tgnt(:, 1, ie, nc),istatus)
          END DO
        ELSE
          ! Bottom boundary: boundery cells share upper edge with internal cells

          ! geoq%cell_stcl(1, nc) = j - 1 + ldim*(i - 2)   ! left bottom
          ! geoq%cell_stcl(2, nc) = j + ldim*(i - 2)       ! bottom
          ! geoq%cell_stcl(3, nc) = j + 1 + ldim*(i - 2)   ! right bottom
          geoq%cell_stcl(4, nc) = j - 1 + ldim * (i - 1)   ! left
          geoq%cell_stcl(5, nc) = j + ldim * (i - 1)       ! self
          geoq%cell_stcl(6, nc) = j + 1 + ldim * (i - 1)   ! right
          geoq%cell_stcl(7, nc) = j - 1 + ldim * (i)       ! left top
          geoq%cell_stcl(8, nc) = j + ldim * (i)           ! top
          geoq%cell_stcl(9, nc) = j + 1 + ldim * (i)       ! right top

          ! Use the edgeSetting routine: Yuanfu Xie 2025-01-08
          ij(1) = j; ij(2) = i
          DO ie = 1, 4
            rsln(2) = rsll(2) * EarthRadius * &
                      DCOS(geoq%cell_cntr(1, nc) + 0.5D0*DBLE(MOD(ie-2,2))*rsll(1))
            CALL edgeSetting_s(ij,nm,geoq%cell_cntr(1, nc),rsln, &
                7,ie,geoq%edge_stcl(:, ie, nc),geoq%coef_func(:, 1, ie, nc), &
                geoq%coef_norm(:, 1, ie, nc),geoq%coef_tgnt(:, 1, ie, nc),istatus)
          END DO
        END IF
      ELSE IF (i .EQ. nmgd(1) - 1) THEN
        ! Top boundary:
        IF (j .EQ. 1) THEN
          ! Left top corner:

          ! geoq%cell_stcl(1, nc) = j - 1 + ldim*(i - 2)   ! left bottom
          geoq%cell_stcl(2, nc) = j + ldim * (i - 2)       ! bottom
          geoq%cell_stcl(3, nc) = j + 1 + ldim * (i - 2)   ! right bottom
          ! geoq%cell_stcl(4, nc) = j - 1 + ldim*(i - 1)   ! left
          geoq%cell_stcl(5, nc) = j + ldim * (i - 1)       ! self
          geoq%cell_stcl(6, nc) = j + 1 + ldim * (i - 1)   ! right
          ! geoq%cell_stcl(7, nc) = j - 1 + ldim*(i)       ! left top
          ! geoq%cell_stcl(8, nc) = j + ldim*(i)           ! top
          ! geoq%cell_stcl(9, nc) = j + 1 + ldim*(i)       ! right top

          ! Use the edgeSetting routine: Yuanfu Xie 2025-01-08
          ij(1) = j; ij(2) = i
          DO ie = 1, 4
            rsln(2) = rsll(2) * EarthRadius * &
                      DCOS(geoq%cell_cntr(1, nc) + 0.5D0*DBLE(MOD(ie-2,2))*rsll(1))
            CALL edgeSetting_s(ij,nm,geoq%cell_cntr(1, nc),rsln, &
                7,ie,geoq%edge_stcl(:, ie, nc),geoq%coef_func(:, 1, ie, nc), &
                geoq%coef_norm(:, 1, ie, nc),geoq%coef_tgnt(:, 1, ie, nc),istatus)
          END DO
        ELSE IF (j .EQ. nmgd(2) - 1) THEN
          ! Right top corner:

          geoq%cell_stcl(1, nc) = j - 1 + ldim * (i - 2)   ! left bottom
          geoq%cell_stcl(2, nc) = j + ldim * (i - 2)       ! bottom
          ! geoq%cell_stcl(3, nc) = j + 1 + ldim*(i - 2)   ! right bottom
          geoq%cell_stcl(4, nc) = j - 1 + ldim * (i - 1)   ! left
          geoq%cell_stcl(5, nc) = j + ldim * (i - 1)       ! self
          ! geoq%cell_stcl(6, nc) = j + 1 + ldim*(i - 1)   ! right
          ! geoq%cell_stcl(7, nc) = j - 1 + ldim*(i)       ! left top
          ! geoq%cell_stcl(8, nc) = j + ldim*(i)           ! top
          ! geoq%cell_stcl(9, nc) = j + 1 + ldim*(i)       ! right top

          ! Use the edgeSetting routine: Yuanfu Xie 2025-01-08
          ij(1) = j; ij(2) = i
          DO ie = 1, 4
            rsln(2) = rsll(2) * EarthRadius * &
                      DCOS(geoq%cell_cntr(1, nc) + 0.5D0*DBLE(MOD(ie-2,2))*rsll(1))
            CALL edgeSetting_s(ij,nm,geoq%cell_cntr(1, nc),rsln, &
                7,ie,geoq%edge_stcl(:, ie, nc),geoq%coef_func(:, 1, ie, nc), &
                geoq%coef_norm(:, 1, ie, nc),geoq%coef_tgnt(:, 1, ie, nc),istatus)
          END DO
        ELSE
          ! Top boundary: boundary cells share lower edge with internal cells

          geoq%cell_stcl(1, nc) = j - 1 + ldim * (i - 2)   ! left bottom
          geoq%cell_stcl(2, nc) = j + ldim * (i - 2)       ! bottom
          geoq%cell_stcl(3, nc) = j + 1 + ldim * (i - 2)   ! right bottom
          geoq%cell_stcl(4, nc) = j - 1 + ldim * (i - 1)   ! left
          geoq%cell_stcl(5, nc) = j + ldim * (i - 1)       ! self
          geoq%cell_stcl(6, nc) = j + 1 + ldim * (i - 1)   ! right
          ! geoq%cell_stcl(7, nc) = j - 1 + ldim*(i)       ! left top
          ! geoq%cell_stcl(8, nc) = j + ldim*(i)           ! top
          ! geoq%cell_stcl(9, nc) = j + 1 + ldim*(i)       ! right top

          ! Use the edgeSetting routine: Yuanfu Xie 2025-01-08
          ij(1) = j; ij(2) = i
          DO ie = 1, 4
            rsln(2) = rsll(2) * EarthRadius * &
                      DCOS(geoq%cell_cntr(1, nc) + 0.5D0*DBLE(MOD(ie-2,2))*rsll(1))
            CALL edgeSetting_s(ij,nm,geoq%cell_cntr(1, nc),rsln, &
                7,ie,geoq%edge_stcl(:, ie, nc),geoq%coef_func(:, 1, ie, nc), &
                geoq%coef_norm(:, 1, ie, nc),geoq%coef_tgnt(:, 1, ie, nc),istatus)
          END DO
        END IF
      ELSE ! Here 1 < i < nmgd(1)-1
        IF (j .EQ. 1) THEN
          ! Left boundary: boundary cells share right edge with internal cells
          ! geoq%cell_stcl(1, nc) = j - 1 + ldim*(i - 2)   ! left bottom
          geoq%cell_stcl(2, nc) = j + ldim * (i - 2)       ! bottom
          geoq%cell_stcl(3, nc) = j + 1 + ldim * (i - 2)   ! right bottom
          ! geoq%cell_stcl(4, nc) = j - 1 + ldim*(i - 1)   ! left
          geoq%cell_stcl(5, nc) = j + ldim * (i - 1)       ! self
          geoq%cell_stcl(6, nc) = j + 1 + ldim * (i - 1)   ! right
          ! geoq%cell_stcl(7, nc) = j - 1 + ldim*(i)       ! left top
          geoq%cell_stcl(8, nc) = j + ldim * (i)           ! top
          geoq%cell_stcl(9, nc) = j + 1 + ldim * (i)       ! right top

          ! Use the edgeSetting routine: Yuanfu Xie 2025-01-08
          ij(1) = j; ij(2) = i
          DO ie = 1, 4
            rsln(2) = rsll(2) * EarthRadius * &
                      DCOS(geoq%cell_cntr(1, nc) + 0.5D0*DBLE(MOD(ie-2,2))*rsll(1))
            CALL edgeSetting_s(ij,nm,geoq%cell_cntr(1, nc),rsln, &
                7,ie,geoq%edge_stcl(:, ie, nc),geoq%coef_func(:, 1, ie, nc), &
                geoq%coef_norm(:, 1, ie, nc),geoq%coef_tgnt(:, 1, ie, nc),istatus)
          END DO
        ELSE IF (j .EQ. nmgd(2) - 1) THEN
          ! Right boundary: boundary cells share left edge with internal cells
          geoq%cell_stcl(1, nc) = j - 1 + ldim * (i - 2)   ! left bottom
          geoq%cell_stcl(2, nc) = j + ldim * (i - 2)       ! bottom
          ! geoq%cell_stcl(3, nc) = j + 1 + ldim*(i - 2)   ! right bottom
          geoq%cell_stcl(4, nc) = j - 1 + ldim * (i - 1)   ! left
          geoq%cell_stcl(5, nc) = j + ldim * (i - 1)       ! self
          ! geoq%cell_stcl(6, nc) = j + 1 + ldim*(i - 1)   ! right
          geoq%cell_stcl(7, nc) = j - 1 + ldim * (i)       ! left top
          geoq%cell_stcl(8, nc) = j + ldim * (i)           ! top
          ! geoq%cell_stcl(9, nc) = j + 1 + ldim*(i)       ! right top

          ! Use the edgeSetting routine: Yuanfu Xie 2025-01-08
          ij(1) = j; ij(2) = i
          DO ie = 1, 4
            rsln(2) = rsll(2) * EarthRadius * &
                      DCOS(geoq%cell_cntr(1, nc) + 0.5D0*DBLE(MOD(ie-2,2))*rsll(1))
            CALL edgeSetting_s(ij,nm,geoq%cell_cntr(1, nc),rsln, &
                7,ie,geoq%edge_stcl(:, ie, nc),geoq%coef_func(:, 1, ie, nc), &
                geoq%coef_norm(:, 1, ie, nc),geoq%coef_tgnt(:, 1, ie, nc),istatus)
          END DO
        ELSE
          ! Interior points:
          ! Cell stencil:
          geoq%cell_stcl(1, nc) = j - 1 + ldim * (i - 2)   ! left bottom
          geoq%cell_stcl(2, nc) = j + ldim * (i - 2)       ! bottom
          geoq%cell_stcl(3, nc) = j + 1 + ldim * (i - 2)   ! right bottom
          geoq%cell_stcl(4, nc) = j - 1 + ldim * (i - 1)   ! left
          geoq%cell_stcl(5, nc) = j + ldim * (i - 1)       ! self
          geoq%cell_stcl(6, nc) = j + 1 + ldim * (i - 1)   ! right
          geoq%cell_stcl(7, nc) = j - 1 + ldim * (i)       ! left top
          geoq%cell_stcl(8, nc) = j + ldim * (i)           ! top
          geoq%cell_stcl(9, nc) = j + 1 + ldim * (i)       ! right top

          ! Edge 1: at cell bottom
          rsln(2) = rsll(2) * EarthRadius * &
                    DCOS(0.5D0 * (geoq%cell_cntr(1, j + ldim * (i - 1)) + &
                                  geoq%cell_cntr(1, j + ldim * (i - 1 - 1))))
          geoq%edge_stcl(1, 1, nc) = j - 1 + ldim * (i - 2) ! left bottom
          geoq%edge_stcl(2, 1, nc) = j - 1 + ldim * (i - 1) ! left top
          geoq%edge_stcl(3, 1, nc) = j + ldim * (i - 2) ! mid bottom
          geoq%edge_stcl(4, 1, nc) = j + ldim * (i - 1) ! mid top
          geoq%edge_stcl(5, 1, nc) = j + 1 + ldim * (i - 2) ! right bottom
          geoq%edge_stcl(6, 1, nc) = j + 1 + ldim * (i - 1) ! right top
          geoq%coef_func(3, 1, 1, nc) = 0.5D0         ! mid bottm
          geoq%coef_func(4, 1, 1, nc) = 0.5D0         ! mid top
          geoq%coef_norm(3, 1, 1, nc) = 1.0D0 / rsln(1) ! mid bottom
          geoq%coef_norm(4, 1, 1, nc) = -1.0D0 / rsln(1) ! mid top
          geoq%coef_tgnt(1, 1, 1, nc) = -0.25D0 / rsln(2) ! left bottom
          geoq%coef_tgnt(2, 1, 1, nc) = -0.25D0 / rsln(2) ! left top
          geoq%coef_tgnt(5, 1, 1, nc) = 0.25D0 / rsln(2) ! right bottom
          geoq%coef_tgnt(6, 1, 1, nc) = 0.25D0 / rsln(2) ! right top
          ! geoq%norm_vctr(1, 1, nc) = 0.0D0
          ! geoq%norm_vctr(2, 1, nc) = -1.0D0
          ! geoq%norm_vct2(1, 1, nc) = 0.0D0
          ! geoq%norm_vct2(2, 1, nc) = -1.0D0
          ! geoq%tgnt_vctr(1, 1, nc) = 1.0D0
          ! geoq%tgnt_vctr(2, 1, nc) = 0.0D0
          ! geoq%tgnt_vct2(1, 1, nc) = 1.0D0
          ! geoq%tgnt_vct2(2, 1, nc) = 0.0D0

          ! Edge 2: at cell right
          rsln(2) = rsll(2) * EarthRadius * &
                    DCOS(geoq%cell_cntr(1, j + ldim * (i - 1)))
          geoq%edge_stcl(1, 2, nc) = j + ldim * (i - 2) ! left bottom
          geoq%edge_stcl(2, 2, nc) = j + 1 + ldim * (i - 2) ! right bottom
          geoq%edge_stcl(3, 2, nc) = j + ldim * (i - 1) ! left mid
          geoq%edge_stcl(4, 2, nc) = j + 1 + ldim * (i - 1) ! right mid
          geoq%edge_stcl(5, 2, nc) = j + ldim * (i) ! left top
          geoq%edge_stcl(6, 2, nc) = j + 1 + ldim * (i) ! right top
          geoq%coef_func(3, 1, 2, nc) = 0.5D0         ! left mid
          geoq%coef_func(4, 1, 2, nc) = 0.5D0         ! right mid
          geoq%coef_norm(3, 1, 2, nc) = -1.0D0 / rsln(2) ! left mid
          geoq%coef_norm(4, 1, 2, nc) = 1.0D0 / rsln(2) ! right mid
          geoq%coef_tgnt(1, 1, 2, nc) = -0.25D0 / rsln(1) ! left bottom
          geoq%coef_tgnt(2, 1, 2, nc) = -0.25D0 / rsln(1) ! right bottom
          geoq%coef_tgnt(5, 1, 2, nc) = 0.25D0 / rsln(1) ! left top
          geoq%coef_tgnt(6, 1, 2, nc) = 0.25D0 / rsln(1) ! right top
          ! geoq%norm_vctr(1, 2, nc) = 1.0D0
          ! geoq%norm_vctr(2, 2, nc) = 0.0D0
          ! geoq%norm_vct2(1, 2, nc) = 1.0D0
          ! geoq%norm_vct2(2, 2, nc) = 0.0D0
          ! geoq%tgnt_vctr(1, 2, nc) = 0.0D0
          ! geoq%tgnt_vctr(2, 2, nc) = 1.0D0
          ! geoq%tgnt_vct2(1, 2, nc) = 0.0D0
          ! geoq%tgnt_vct2(2, 2, nc) = 1.0D0

          ! Edge 3: at cell top
          rsln(2) = rsll(2) * EarthRadius * &
                    DCOS(0.5D0 * (geoq%cell_cntr(1, j + ldim * (i - 1)) + &
                                  geoq%cell_cntr(1, j + ldim * (i))))
          geoq%edge_stcl(1, 3, nc) = j + 1 + ldim * (i - 1) ! right bottom
          geoq%edge_stcl(2, 3, nc) = j + 1 + ldim * (i) ! right top
          geoq%edge_stcl(3, 3, nc) = j + ldim * (i - 1) ! mid bottom
          geoq%edge_stcl(4, 3, nc) = j + ldim * (i) ! mid top
          geoq%edge_stcl(5, 3, nc) = j - 1 + ldim * (i - 1) ! left bottom
          geoq%edge_stcl(6, 3, nc) = j - 1 + ldim * (i) ! left top
          geoq%coef_func(3, 1, 3, nc) = 0.5D0         ! mid bottom
          geoq%coef_func(4, 1, 3, nc) = 0.5D0         ! mid top
          geoq%coef_norm(3, 1, 3, nc) = -1.0D0 / rsln(1) ! mid bottom
          geoq%coef_norm(4, 1, 3, nc) = 1.0D0 / rsln(1) ! mid top
          geoq%coef_tgnt(1, 1, 3, nc) = -0.25D0 / rsln(2) ! right bottom
          geoq%coef_tgnt(2, 1, 3, nc) = -0.25D0 / rsln(2) ! right top
          geoq%coef_tgnt(5, 1, 3, nc) = 0.25D0 / rsln(2) ! left bottom
          geoq%coef_tgnt(6, 1, 3, nc) = 0.25D0 / rsln(2) ! left top
          ! geoq%norm_vctr(1, 3, nc) = 0.0D0
          ! geoq%norm_vctr(2, 3, nc) = 1.0D0
          ! geoq%norm_vct2(1, 3, nc) = 0.0D0
          ! geoq%norm_vct2(2, 3, nc) = 1.0D0
          ! geoq%tgnt_vctr(1, 3, nc) = -1.0D0
          ! geoq%tgnt_vctr(2, 3, nc) = 0.0D0
          ! geoq%tgnt_vct2(1, 3, nc) = -1.0D0
          ! geoq%tgnt_vct2(2, 3, nc) = 0.0D0

          ! Edge 4: at cell left
          rsln(2) = rsll(2) * EarthRadius * &
                    DCOS(geoq%cell_cntr(1, j + ldim * (i - 1)))
          geoq%edge_stcl(1, 4, nc) = j - 1 + ldim * (i) ! left top
          geoq%edge_stcl(2, 4, nc) = j + ldim * (i) ! right top
          geoq%edge_stcl(3, 4, nc) = j - 1 + ldim * (i - 1) ! left mid
          geoq%edge_stcl(4, 4, nc) = j + ldim * (i - 1) ! right mid
          geoq%edge_stcl(5, 4, nc) = j - 1 + ldim * (i - 2) ! left bottom
          geoq%edge_stcl(6, 4, nc) = j + ldim * (i - 2) ! right bottom
          geoq%coef_func(3, 1, 4, nc) = 0.5D0         ! mid left
          geoq%coef_func(4, 1, 4, nc) = 0.5D0         ! mid right
          geoq%coef_norm(3, 1, 4, nc) = 1.0D0 / rsln(2) ! mid left
          geoq%coef_norm(4, 1, 4, nc) = -1.0D0 / rsln(2) ! mid right
          geoq%coef_tgnt(1, 1, 4, nc) = -0.25D0 / rsln(1) ! left bottom
          geoq%coef_tgnt(2, 1, 4, nc) = -0.25D0 / rsln(1) ! left top
          geoq%coef_tgnt(5, 1, 4, nc) = 0.25D0 / rsln(1) ! left bottom
          geoq%coef_tgnt(6, 1, 4, nc) = 0.25D0 / rsln(1) ! right bottom
          ! geoq%norm_vctr(1, 4, nc) = -1.0D0
          ! geoq%norm_vctr(2, 4, nc) = 0.0D0
          ! geoq%norm_vct2(1, 4, nc) = -1.0D0
          ! geoq%norm_vct2(2, 4, nc) = 0.0D0
          ! geoq%tgnt_vctr(1, 4, nc) = 0.0D0
          ! geoq%tgnt_vctr(2, 4, nc) = -1.0D0
          ! geoq%tgnt_vct2(1, 4, nc) = 0.0D0
          ! geoq%tgnt_vct2(2, 4, nc) = -1.0D0
        END IF
      END IF

      !! IF halo_width > 1
      ! IF halo width = 1, Cell 05 current cell, 01 - 29 except 05 are neighber cells
      !      +----+----+----+
      !      | 07 | 08 | 09 |
      !      +----+----+----+
      !      | 04 | 05 | 06 |
      !      +----+----+----+
      !      | 01 | 02 | 03 |
      !      +----+----+----+
      !
      ! IF halo width = 2, Cell 05 is still current cell, 01 - 25 except 05 are neighber cells
      ! The first layer idx is from 1 to (1 + 2*layerIndex) ** 2, [1 to 9]
      ! The second layer idx is from (1 + 2*(layerIndex - 1)) ** 2 + 1 to (1 + 2*layerIndex) ** 2, [10, 25]
      ! ...
      !
      !      +----+----+----+----+----+
      !      | 21 | 22 | 23 | 24 | 25 |
      !      +----+----+----+----+----+
      !      | 19 | 07 | 08 | 09 | 20 |
      !      +----+----+----+----+----+
      !      | 17 | 04 | 05 | 06 | 18 |
      !      +----+----+----+----+----+
      !      | 15 | 01 | 02 | 03 | 16 |
      !      +----+----+----+----+----+
      !      | 10 | 11 | 12 | 13 | 14 |
      !      +----+----+----+----+----+

      ! IF halo width = 3
      !      +----+----+----+----+----+----+----+
      !      | 43 | 44 | 45 | 46 | 47 | 48 | 49 |
      !      +----+----+----+----+----+----+----+
      !      | 41 | 21 | 22 | 23 | 24 | 25 | 42 |
      !      +----+----+----+----+----+----+----+
      !      | 39 | 19 | 07 | 08 | 09 | 20 | 40 |
      !      +----+----+----+----+----+----+----+
      !      | 37 | 17 | 04 | 05 | 06 | 18 | 38 |
      !      +----+----+----+----+----+----+----+
      !      | 35 | 15 | 01 | 02 | 03 | 16 | 36 |
      !      +----+----+----+----+----+----+----+
      !      | 33 | 10 | 11 | 12 | 13 | 14 | 34 |
      !      +----+----+----+----+----+----+----+
      !      | 26 | 27 | 28 | 29 | 30 | 31 | 32 |
      !      +----+----+----+----+----+----+----+
      !
      ! FOR any input x, y, layer = MAXVAL(|x-i|, |y-j|)
      ! FOR left cell  x - i == -layer
      ! FOR right cell x - i == layer
      ! FOR top cell y - j == layer
      ! FOR bottom cell y -j == -layer
      ! The first cell x0, y0 = i - layer, j - layer

      cell_stcl_width = NINT((SQRT(DBLE(para%numStclPerCell)) - 1D0) / 2D0)

      ! PRINT*, 'HYJ+++ cell_stcl_width', cell_stcl_width, ldim, para%num_cell, ldim, rdim, nend

      ! i => nend(1) from bottom up, j => nend(2) from left right
      ! cell [i, j] is inner cell
      DO m = 2, cell_stcl_width

        !! bottom
        !! cond :       i - m >= nbeg(1) & (nbeg(2)<=j-m+n-1<=nend(2)-1)
        !!
        DO n = 1, 1 + 2 * m

          IF ((i - m) .GE. nbeg(1) .AND. ((j - m + n - 1) .LE. (nend(2) - 1)) .AND. ((j - m + n - 1) .GE. (nbeg(2)))) THEN
            geoq%cell_stcl((1 + 2 * (m - 1))**2 + n, nc) = (j - m) + ldim * (i - 1 - m) + n - 1
          ELSE
            geoq%cell_stcl((1 + 2 * (m - 1))**2 + n, nc) = 0
          END IF

        END DO

        !! LEFT
        !! cond:  j - m >= nbeg(2) &&  nbeg(1) <= i - m + n -<= nend(1) - 1
        DO n = 1, 1 + 2 * (m - 1)

          IF (((j - m) .GE. nbeg(2)) .AND. ((i - m + n - 1) .GE. nbeg(1)) .AND. ((i - m + n) .LE. (nend(1) - 1))) THEN
            geoq%cell_stcl((1 + 2 * (m - 1))**2 + 1 + 2 * m + 2 * n - 1, nc) = (j - m) + ldim * (i - 1 - m + n)
          ELSE
            geoq%cell_stcl((1 + 2 * (m - 1))**2 + 1 + 2 * m + 2 * n - 1, nc) = 0
          END IF

        END DO

        !! RIGHT
        !! cond:  j + m <= nbeg(2) - 1 &&  nbeg(1) <= i - m + n  <= nend(1) - 1
        DO n = 1, 1 + 2 * (m - 1)

          IF (((j + m) .LE. (nend(2) - 1)) .AND. ((i - m + n - 1) .GE. nbeg(1)) .AND. ((i - m + n) .LE. (nend(1) - 1))) THEN
            geoq%cell_stcl((1 + 2 * (m - 1))**2 + 1 + 2 * m + 2 * n, nc) = (j + m) + ldim * (i - 1 - m + n)
          ELSE
            geoq%cell_stcl((1 + 2 * (m - 1))**2 + 1 + 2 * m + 2 * n, nc) = 0
          END IF

          ! IF ( ((j + m ) .GT. (nend(2)- 1)) .OR. ((i - m ) .LT. (nbeg(1))) ) THEN
          !   geoq%cell_stcl((1+2*(m-1))**2 + 1 + 2*m + 2 * n, nc) = 0
          ! ELSE
          !   geoq%cell_stcl((1+2*(m-1))**2 + 1 + 2*m + 2 * n, nc) = (j + m) + ldim * (i - 1 - m + n)
          ! END IF
        END DO

        !! TOP
        !! cond i + m <= nend(1) - 1 && (nbeg(2)<=j-m+n-1<=nend(2)-1)
        DO n = 1, 1 + 2 * m

          IF ((i + m) .LE. (nend(1) - 1) .AND. ((j - m + n - 1) .LE. (nend(2) - 1)) .AND. ((j - m + n - 1) .GE. (nbeg(2)))) THEN
            geoq%cell_stcl((1 + 2 * m)**2 - (1 + 2 * m) + n, nc) = (j - m) + ldim * (i - 1 + m) + n - 1
          ELSE
            geoq%cell_stcl((1 + 2 * (m - 1))**2 + 1 + 2 * m + 2 * n, nc) = 0
          END IF

        END DO

      END DO

    END DO
  END DO
  !$OMP END PARALLEL DO

END SUBROUTINE latlon_precals
