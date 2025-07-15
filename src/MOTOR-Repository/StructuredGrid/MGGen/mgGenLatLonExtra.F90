!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.mgGen.mgGenLatLonExtra_m
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2020/12/28, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module extends the mgGenLatLon_t, which add the function that could create the multi grid without creating files.
! @note
! @warning
! @attention
MODULE mgGenLatLonExtra_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE gridStruct_m, ONLY: gridParams_t, gridGeoQty_t, multiLevel_t, &
                          lengthFile, file_params, file_geoqty

  USE mgGen_m, ONLY: mgGenLatLon_t
  USE YAMLRead_m

  TYPE, EXTENDS(mgGenLatLon_t) :: mgGenLatLonExtra_t
  CONTAINS
    PROCEDURE, PUBLIC :: genGridAndPreCal

  END TYPE mgGenLatLonExtra_t
CONTAINS

  SUBROUTINE genGridAndPreCal(this, configFile, mlgrid)
    USE parameters_m, ONLY: degree2Radian

    CLASS(mgGenLatLonExtra_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(multiLevel_t), INTENT(INOUT) :: mlgrid

    ! Local variables:
    LOGICAL :: divided  ! If a level can still be divided
    INTEGER(i_kind) :: m, ncount(2), nv, nc, npower

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

    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'haloWidth', this%halo_width)

    ! PRINT *, 'HYJ+++ genGridAndPreCal halo width is set to ', this%halo_width

    IF (istatus .NE. 0) THEN
      this%halo_width = 1
    END IF

    ! Determine how many multigrid levels:
    this%mglvls = 1
    this%mgstts = 1
    divided = .TRUE.
    ncount = this%numgrd - 3
    DO WHILE (.TRUE.)
      ! ncount counts the interior points only;
      ! thus <= 2 means further division will leave no interior points
      ! IF (MOD(ncount(1), 2) .NE. 0 .OR. ncount(1) .LE. 2 .OR. &
      !     MOD(ncount(2), 2) .NE. 0 .OR. ncount(2) .LE. 2) THEN
      IF (MOD(ncount(1), 2) .NE. 0 .OR. ncount(1) .LE. 1 + this%halo_width .OR. &
          MOD(ncount(2), 2) .NE. 0 .OR. ncount(2) .LE. 1 + this%halo_width) THEN
        EXIT
      END IF
      ncount = ncount / 2
      this%mglvls = this%mglvls + 1
    END DO
    this%mgends = this%mglvls
    WRITE (*, 10) this%mglvls, ncount + 1
10  FORMAT('Total MG-levels: ', i2, ' coarsest grid interior points: ', 2I6)

    ! Allocate all mg level grids:
    CALL mlgrid%allc_arrays(this%mglvls, this%mgstts, this%mgends)

    ! Latlon grid parameters:
    mlgrid%params(:)%numStclPerEdge = 7    ! tgnt requires 6, while norm and func need 2
    mlgrid%params(:)%numQuadPerEdge = 1    ! Number of quadrature points per edge

    mlgrid%params(:)%numStclPerCell = (1 + 2 * this%halo_width)**2

    mlgrid%params(:)%numStclPerVrtx = 4
    mlgrid%params(:)%numVrtxPerCell = 4

    mlgrid%params(:)%num_chld = 4

    ! Vertices and cells count those fictitious points:
    ncount = this%numgrd - 2
    DO m = this%mglvls, 1, -1
      ! Numbers of vertices and cells:
      mlgrid%params(m)%num_vrtx = (ncount(1) + 2) * (ncount(2) + 2)
      mlgrid%params(m)%num_cell = (ncount(1) + 1) * (ncount(2) + 1)
      ncount = (ncount - 1) / 2 + 1
    END DO

    WRITE (*, 9) this%mglvls
9   FORMAT('Calling allc_geoqty with number ', I2, ' of levels')
    CALL mlgrid%allc_geoqty

    ncount = this%numgrd
    DO m = this%mglvls, 1, -1
      WRITE (*, 8) m, ncount
8     FORMAT('Calling latlon_geoGrid at multigrid level: ', I2, &
             ' with gridpoints (Fic) ', 2I5)
      CALL latlon_geoGrid(ncount, this%domain, &
                          mlgrid%params(m), mlgrid%GeoQty(m))
      ! Count interior points first and add fictitious:
      ncount = (ncount - 3) / 2 + 3
    END DO

    ncount = this%numgrd
    npower = 0
    DO m = mlgrid%lvl_ends, mlgrid%lvl_stts, -1
      CALL latlon_precals(ncount, this%resoln * 2.0D0**npower, &
                          mlgrid%params(m), mlgrid%GeoQty(m))
      ! Count reductions are on interior cells only:
      ncount = (ncount - 3) / 2 + 3
      npower = npower + 1
    END DO

  END SUBROUTINE

END MODULE mgGenLatLonExtra_m
