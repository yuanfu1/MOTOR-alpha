!!--------------------------------------------------------------------------------------------------
! PROJECT           : Multi-Grid generation
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : Beta 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2020/05/29, @GBA-MWF, Shenzhen
!   Modified by Zilong Qin (zilong.qin@gmail.com), 2020/12/15, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! # Multi-Grid Generation
!!  *This module defines an abstraction of general grid generation for
!!  various grid types, including latlon, icosahedral etc. It assumes
!!  a finite volume scheme is used, where control cells, cell vertices
!!  , edge and area information are required*
!!
!! ## Latlon grid:
!!  *The latlon grid is assumed the grid is symmetric at the origin
!!  of the latlon center (0,0). Any latlon domain will be transformed
!!  by an Euler rotation and interpolated from the symmetric analysis
!!  and forecast domain by a slint package.*
!!
!! <center>![mgGen_multigrid](../../image/mgGen_multigrid.png)</center>
!! The index of the grid cell is listed as shown above, for Finer grid G2 and Coarser G1,
!! the inner cell region remains unchanged, while the boundary cells always surround the
!! region.
!!
!! @author Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
MODULE mgGen_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE gridStruct_m, ONLY: gridParams_t, gridGeoQty_t, multiLevel_t, &
                          lengthFile, file_params, file_geoqty

  IMPLICIT NONE

  PUBLIC

  PUBLIC :: mgGen_t, mgGenLatLon_t, mgGenIcosTr_t

  TYPE, ABSTRACT :: mgGen_t
    INTEGER(i_kind), PUBLIC :: mglvls, mgstts, mgends
    REAL(r_kind) :: resoln(2)
    TYPE(multiLevel_t) :: mlgrid
  CONTAINS
    PROCEDURE(geoGrid), DEFERRED :: mgGen_geoGrid   ! grid geometry locations
    PROCEDURE(preCals), DEFERRED :: mgGen_preCals   ! pre-calculate all needed information
    PROCEDURE, PASS :: destroy        !< Dealloc the geometry instance
  END TYPE mgGen_t

  TYPE, EXTENDS(mgGen_t) :: mgGenLatLon_t
    INTEGER(i_kind), PUBLIC :: numgrd(2) ! Number of grid points in lat and lon
    REAL(r_kind) :: domain(2, 2)  ! Range of lat and lon in degrees
    INTEGER(i_kind) :: halo_width
  CONTAINS
    PROCEDURE :: mgGen_geoGrid => mgGenLatLon_geoGrid
    PROCEDURE :: mgGen_preCals => mgGenLatLon_preCals
  END TYPE mgGenLatLon_t

  TYPE, EXTENDS(mgGen_t) :: mgGenIcosTr_t
    INTEGER(i_kind) :: odr_accu ! Order of accuracy
  CONTAINS
    PROCEDURE :: mgGen_geoGrid => mgGenIcosTr_geoGrid
    PROCEDURE :: mgGen_preCals => mgGenIcosTr_preCals
  END TYPE mgGenIcosTr_t

  ABSTRACT INTERFACE
    SUBROUTINE geoGrid(this, configFile)
      IMPORT :: mgGen_t
      CLASS(mgGen_t) :: this
      CHARACTER(LEN=1024), INTENT(IN) :: configFile
    END SUBROUTINE geoGrid

    SUBROUTINE preCals(this)
      IMPORT :: mgGen_t
      CLASS(mgGen_t) :: this
    END SUBROUTINE preCals
  END INTERFACE

CONTAINS
  SUBROUTINE destroy(this)
    CLASS(mgGen_t) :: this

    CALL this%mlgrid%delt_geoqty
    CALL this%mlgrid%delt_arrays
  END SUBROUTINE

  SUBROUTINE mgGenLatLon_geoGrid(this, configFile)
    USE parameters_m, ONLY: degree2Radian

    CLASS(mgGenLatLon_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile

    ! Local variables:
    LOGICAL :: divided  ! If a level can still be divided
    INTEGER(i_kind) :: m, ncount(2), nv, nc

    ! Get domain information:
    ! Clean up arrays for initialization:
    WRITE (*, 12)
12  FORMAT('Calling getDomain_latlon...')
    CALL getDomain_latlon(configFile, nv, nc, this%numgrd, this%domain)

    ! Determine the grid resolution:
    this%domain = this%domain * degree2radian
    ! numgrid records the total grid points including the fictitious points:
    ! In fact the actual grid points in the domain is numgrd-2, or numgrd-3 cell
    ! in each direction:
    this%resoln = (this%domain(2, :) - this%domain(1, :)) / DBLE(this%numgrd - 3)
    WRITE (*, 11) this%resoln, this%numgrd
11  FORMAT('Grid resolution: ', 2F14.6, /, &
           'Number of gridpoints (fictitious): ', 2I5)

    ! Determine how many multigrid levels:
    this%mglvls = 1
    this%mgstts = 1
    divided = .TRUE.
    ncount = this%numgrd - 3
    DO WHILE (.TRUE.)
      ! ncount counts the interior points only;
      ! thus <= 2 means further division will leave no interior points
      IF (MOD(ncount(1), 2) .NE. 0 .OR. ncount(1) .LE. 2 .OR. &
          MOD(ncount(2), 2) .NE. 0 .OR. ncount(2) .LE. 2) THEN
        EXIT
      END IF
      ncount = ncount / 2
      this%mglvls = this%mglvls + 1
    END DO
    this%mgends = this%mglvls
    WRITE (*, 10) this%mglvls, ncount + 1
10  FORMAT('Total MG-levels: ', i2, ' coarsest grid interior points: ', 2I6)

    ! Allocate all mg level grids:
    CALL this%mlgrid%allc_arrays(this%mglvls, this%mgstts, this%mgends)

    ! Latlon grid parameters:
    ! Yuanfu Xie changed numStclPerEdge to 7 for implementing boundary 
    ! interpolation on 2024-12-20.
    ! The 7th cell of a boundary cell is marked as following:
    ! +---+---+---+
    ! |   | 7 |   |
    ! +---+---+---+
    ! | 4 | 5 | 6 |
    ! +---+---+---+
    ! | 1 | 2 | 3 |
    ! +---+-x-+---+
    ! For all interior point, 
    ! +---+---+---+
    ! |   |   |   |
    ! +---+---+---+
    ! | 4 |5/7| 6 |
    ! +---+-x-+---+
    ! | 1 | 2 | 3 |
    ! +---+---+---+
    ! the 7th is set to side on the 5th cell with coefficient of zero.
    this%mlgrid%params(:)%numStclPerEdge = 7    ! tgnt requires 6, while norm and func need 2
    this%mlgrid%params(:)%numQuadPerEdge = 1    ! Number of quadrature points per edge
    this%mlgrid%params(:)%numStclPerCell = 9    ! each edge tgnt uses 6, cell involves 9
    this%mlgrid%params(:)%numStclPerVrtx = 4
    this%mlgrid%params(:)%numVrtxPerCell = 4

    this%mlgrid%params(:)%num_chld = 4

    ! Vertices and cells count those fictitious points:
    ncount = this%numgrd - 2
    DO m = this%mglvls, 1, -1
      ! Numbers of vertices and cells:
      this%mlgrid%params(m)%num_vrtx = (ncount(1) + 2) * (ncount(2) + 2)
      this%mlgrid%params(m)%num_cell = (ncount(1) + 1) * (ncount(2) + 1)
      ncount = (ncount - 1) / 2 + 1
    END DO

    WRITE (*, 9) this%mglvls
9   FORMAT('Calling allc_geoqty with number ', I2, ' of levels')
    CALL this%mlgrid%allc_geoqty

    ncount = this%numgrd
    DO m = this%mglvls, 1, -1
      WRITE (*, 8) m, ncount
8     FORMAT('Calling latlon_geoGrid at multigrid level: ', I2, &
             ' with gridpoints (Fic) ', 2I5)
      CALL latlon_geoGrid(ncount, this%domain, &
                          this%mlgrid%params(m), this%mlgrid%GeoQty(m))
      ! Count interior points first and add fictitious:
      ncount = (ncount - 3) / 2 + 3
    END DO

    ! Output grid structure only: using .TRUE.:
    CALL this%mlgrid%writ_mlgrid(.TRUE.)

    CALL this%mlgrid%delt_geoqty
    CALL this%mlgrid%delt_arrays

  END SUBROUTINE mgGenLatLon_geoGrid

  SUBROUTINE mgGenLatLon_precals(this)
    CLASS(mgGenLatLon_t) :: this

    ! Local variables:
    INTEGER(i_kind) :: m, ncount(2), npower

    ! Read in multigrid data file:
    CALL this%mlgrid%read_mlgrid(.TRUE.)

    ncount = this%numgrd
    npower = 0
    DO m = this%mlgrid%lvl_ends, this%mlgrid%lvl_stts, -1
      CALL latlon_precals(ncount, this%resoln * 2.0D0**npower, &
                          this%mlgrid%params(m), this%mlgrid%GeoQty(m))
      ! Count reductions are on interior cells only:
      ncount = (ncount - 3) / 2 + 3
      npower = npower + 1
    END DO

    ! Output multigrid using .FALSE.:
    CALL this%mlgrid%writ_mlgrid(.FALSE.)

    CALL this%mlgrid%delt_geoqty
    CALL this%mlgrid%delt_arrays

  END SUBROUTINE mgGenLatLon_precals

  SUBROUTINE mgGenIcosTr_geoGrid(this, configFile)
    CLASS(mgGenIcosTr_t) :: this

    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    ! Local variables:
    ! icosahedral level: G6 triangle grid is about 80km
    ! G12 is about 1.25km and G13 0.625km G14: 0.3125km
    ! G15: 0.15625 - 156 meters
    INTEGER(i_kind), PARAMETER :: ng_lvl = 8

    INTEGER(i_kind) :: num_vrtx(ng_lvl), num_cell(ng_lvl)
    DATA num_vrtx/0, 162, 642, 2562, 10242, 40962, 163842, 655362/, &
      num_cell/0, 320, 1280, 5120, 20480, 81920, 327680, 1310720/

    INTEGER(i_kind) :: i

    PRINT *, 'Calling getDomain_icosTr...'
    CALL getDomain_icosTr(configFile, this%mglvls, this%mgstts, this%mgends)
    PRINT *, 'Icos info: ', this%mgstts, &
      this%mgends, this%mglvls

    ! Allocate the parameter array:
    CALL this%mlgrid%allc_arrays(this%mglvls, this%mgstts, this%mgends)

    ! Icos triangle grid parameters:
    this%mlgrid%params(:)%numVrtxPerCell = 3
#ifdef MID_EDGE
    this%mlgrid%params(:)%numStclPerCell = 13    ! based on 10 stencil per edge
    this%mlgrid%params(:)%numStclPerEdge = 10    ! 10 stencil for a first order accuracy
    this%mlgrid%params(:)%numQuadPerEdge = 1     ! Number of quadrature points per edge
    this%mlgrid%params(:)%numStclPerVrtx = 12    ! For icos not sure yet???
#else
    this%mlgrid%params(:)%numStclPerCell = 24    ! based on 14 stencil per edge
    this%mlgrid%params(:)%numStclPerEdge = 14    ! 10 stencil for a second order accuracy
    this%mlgrid%params(:)%numQuadPerEdge = 2     ! Number of quadrature points per edge
    this%mlgrid%params(:)%numStclPerVrtx = 24    ! For icos not sure yet???
#endif

    this%mlgrid%params(:)%num_chld = 4

    ! Check if higher G level information is needed:
    IF (this%mgends .GT. ng_lvl) THEN
      PRINT *, 'mgGenIcosTr_geoGrid: need higher G level grid vertex and cell info!', &
        this%mgends, ng_lvl
      STOP
    END IF

    DO i = this%mgstts, this%mgends
      PRINT *, 'I::', i, num_vrtx(i), num_cell(i)
      this%mlgrid%Params(i)%num_vrtx = num_vrtx(i)
      this%mlgrid%Params(i)%num_cell = num_cell(i)
    END DO

    PRINT *, 'Calling allc_geoqty...', this%mglvls
    CALL this%mlgrid%allc_geoqty

    ! Output grid structure only: using .TRUE.:
    CALL this%mlgrid%writ_mlgrid(.TRUE.)

    CALL this%mlgrid%delt_geoqty
    CALL this%mlgrid%delt_arrays

  END SUBROUTINE mgGenIcosTr_geoGrid

  SUBROUTINE mgGenIcosTr_precals(this)
    CLASS(mgGenIcosTr_t) :: this

    CALL this%mlgrid%allc_arrays(1, 2, 2)
  END SUBROUTINE mgGenIcosTr_precals

END MODULE mgGen_m
