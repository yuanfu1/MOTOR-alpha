!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.geometry
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie, Ning Wang, Jinrong Jiang, Zilong Qin
! VERSION           : 0.1
! HISTORY           :
! This module is reforged from the original grid/icosTriangleGeom.F90
!!    2019-1 created by Yuanfu Xie
!!    2019-5 modified by Yuanfu Xie for reading the triangle grid
!!                                      in multigrid data format
!  Modified by Zilong Qin (zilong.qin@gmail.com), 2020/11/19, @GBA-MWF, Shenzhen, for encapsulation and parallelization
!  Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/06/09, for adding safeguard checking the
!     consistency of a multigrid and proc_layout as currently we require each processor handles the
!     same number of grid points. This can be relaxed with an improvement of MOTOR-DA multigrid.
!  Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2023/03/07, for adding Gauss-Legendre quadrature weights
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This is a MPDD based Lat-Lon geometry struct based on Xie's mgGen-LatLon grid.
!! @see
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
! @warning
! @attention
MODULE geometry_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE gridStruct_m, ONLY: gridParams_t, gridGeoQty_t, multiLevel_t, &
                          lengthFile, file_params, file_geoqty
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE multiGrid_m, ONLY: multiGrid_t
  USE mgGenLatLonExtra_m, ONLY: mgGenLatLonExtra_t
  USE namelist_gm_m
  USE YAMLRead_m

  PRIVATE
  PUBLIC :: geometry_t

!> @class geometry_t
!> @brief
!! The class which contains all the 2D geometry infos and parallel operators for model and DA application.
! @see
! @note
! @warning
! @attention
  TYPE geometry_t
    TYPE(multiGrid_t), ALLOCATABLE :: mg !< Multi-grid structure of the goemetry
    TYPE(mpddGlob_t), POINTER :: mpdd        !< Pointer of global multiprocessing domain decomposition class

  CONTAINS
    FINAL :: destructor
    PROCEDURE, PUBLIC :: initialize
    PROCEDURE, PUBLIC :: destroy
  END TYPE geometry_t

  ! INTERFACE geometry_t
  !   PROCEDURE :: constructor
  ! END interface geometry_t

CONTAINS

  SUBROUTINE boundingMG(mg, mg_pointer)
    TYPE(multiGrid_t), TARGET :: mg
    CLASS(*), POINTER :: mg_pointer

    mg_pointer => mg
  END SUBROUTINE boundingMG

!> @brief
!! Initialize the geometry object. This function should be call before using any other functions of this class.
! @see
! @note
! @warning
! @attention
  SUBROUTINE initialize(this, configFile, mpdd, opts)
    IMPLICIT NONE
    CLASS(geometry_t) this
    TYPE(mpddGlob_t), TARGET, INTENT(IN) :: mpdd        !< Global multiprocessing domain decomposition class pointer
    CHARACTER(LEN=1024), INTENT(IN) :: configFile   !< File name and path of configuration of geometry
    INTEGER(i_kind), OPTIONAL :: opts

    TYPE(mgGenLatLonExtra_t) :: mggen

    ! Temporary variables
    INTEGER(i_kind) :: i, num_children, num_icell_temp, mg_coarsest = 0, mg_finest = 0, istatus
    INTEGER(i_kind) ::  mgGenCoarsest = 0, mgGenFinest = 0
    INTEGER(i_kind) :: vLevel
    INTEGER(i_kind) :: haloWidth

    ! Auxiliary variables
    INTEGER(i_kind), ALLOCATABLE       :: proc_layout(:, :)         !< Layout of MPI/grid spliting
    INTEGER(i_kind), ALLOCATABLE       :: t_steps_mg(:)         !< Layout of MPI/grid spliting
    INTEGER(i_kind) :: dimCell_global(2)        !< Dimension of grid cells in Lat-lon grid
    REAL(r_kind) :: start_time, end_time
    this%mpdd => mpdd
    ! ALLOCATE multiGrit_t

    PRINT *, 'Start initialization of Geometry.'
    ! Read mgGen grid: .FALSE. indicates reading both grid structure
    ! and precal values: nread_mlgrid allocates memory for multi-level
    ! IF (mpdd%isBaseProc()) CALL mggen%read_mlgrid(.FALSE.)
    IF (mpdd%isBaseProc()) CALL mggen%genGridAndPreCal(configFile, mggen%mlgrid)

    IF (mpdd%isBaseProc()) mgGenCoarsest = mggen%mlgrid%lvl_stts
    IF (mpdd%isBaseProc()) mgGenFinest = mggen%mlgrid%lvl_ends
    CALL mpdd%bCastInt(mgGenCoarsest)
    CALL mpdd%bCastInt(mgGenFinest)
    ! Yuanfu Xie moved this to here as it is meaningless at its original position above: 2025-03-28
    PRINT *, 'Geometry setting: ', mgGenCoarsest, mgGenFinest 

    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', mg_coarsest)
    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', mg_finest)

    IF (mg_coarsest < mgGenCoarsest) mg_coarsest = mgGenCoarsest
    IF (mg_finest > mgGenFinest) mg_finest = mgGenFinest

    ! Passing the multi-level grid and precal to multigrid_t variable:
    !-------------------Begin multigrid assignment--------------------
    IF (PRESENT(opts)) THEN
      ! IF (opts == 1) mg_coarsest = mg_finest
    END IF

    ! Boardcast the mg_coarsest and mg_finest
    WRITE (*, 40) mg_coarsest, mg_finest
40  FORMAT('Geometry - Based on num grid, the coarsest level: ', I3, ' Finest level: ', I3)

    IF (.NOT. mpdd%isBaseProc()) THEN
      ALLOCATE (mggen%mlgrid%Params(mg_coarsest:mg_finest), &
                mggen%mlgrid%GeoQty(mg_coarsest:mg_finest))
    END IF
    ! Allocate the single grid in the multi grid
    ALLOCATE (proc_layout(mg_coarsest:mg_finest, 2))
    ALLOCATE (t_steps_mg(mg_coarsest:mg_finest))

    ! Read the parallel layout from namelist
    ! Get the dimension of the grid cells
    CALL namelist_gm(configFile, proc_layout, t_steps_mg, start_time, end_time, mg_coarsest, mg_finest, dimCell_global, vLevel)
    CALL mpdd%barrier

    IF (proc_layout(mg_finest, 1) * proc_layout(mg_finest, 2) /= mpdd%nProc) THEN
      IF (mpdd%nProc == 1) THEN
        PRINT *, 'Warning: The number of processors is not consistent with the parallel layout, using 1 X 1 layout.'
        proc_layout = 1
      ELSE
        PRINT *, "mpdd%nProc", mpdd%nProc
        PRINT *, "proc_layout(mg_finest, 1) * proc_layout(mg_finest, 2)", proc_layout(mg_finest, 1) * proc_layout(mg_finest, 2)
        PRINT *, 'Error: The number of processors is not consistent with the parallel layout.'
        STOP
      END IF
    END IF

    IF (PRESENT(opts)) THEN
      IF (opts == 1) proc_layout = 1
    END IF
    ! CALL mpdd%isNProcMatch(proc_layout(this%mg%mg_finest, 1)*proc_layout(this%mg%mg_finest, 2))

    DO i = mg_coarsest, mg_finest - 1
      IF (MOD((dimCell_global(1) - 2), proc_layout(i + 1, 1) * 2**(mgGenFinest - i)) .NE. 0) THEN
        proc_layout(i + 1, 1) = proc_layout(i + 1, 1) / 2
        PRINT *, 'Layout on g ', i, ' changed to ', proc_layout(i, :)
        ! EXIT  ! Acceptable coarsest level and continue to analysis
      END IF
      IF (MOD((dimCell_global(2) - 2), proc_layout(i + 1, 2) * 2**(mgGenFinest - i)) .NE. 0) THEN
        proc_layout(i + 1, 2) = proc_layout(i + 1, 2) / 2
        PRINT *, 'Layout on g ', i, ' changed to ', proc_layout(i, :)
        ! EXIT  ! Acceptable coarsest level and continue to analysis
      END IF
    END DO

    ! Check the HALO WIDTH setting is valid for multigrid
    ! The HALO WIDTH should less than NX and NY of inner cell in sub threads
    ! Note: check mg_coarsest only, if mg_coarsest is valid, the finner grids are also valid

    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'haloWidth', haloWidth)

    IF (istatus .NE. 0) THEN
      haloWidth = 1
    END IF

    PRINT *, 'Geometry: halo_width is set to value ', haloWidth, (dimCell_global(1) - 2), proc_layout(mg_coarsest, 1), (dimCell_global(2) - 2), proc_layout(mg_coarsest, 2)

    ! Check the consistency of multigrid and proc_layout by Yuanfu Xie 2022-06-09
    ! Note: dimCell_global adds 1 grid point to count the fictitious points on the right
    IF (mpdd%isBaseProc()) THEN
      DO i = mg_coarsest, mg_finest
        IF (MOD((dimCell_global(1) - 2), proc_layout(i, 1) * 2**(mg_finest - i)) .EQ. 0 .AND. &
            MOD((dimCell_global(2) - 2), proc_layout(i, 2) * 2**(mg_finest - i)) .EQ. 0) THEN
          EXIT  ! Acceptable coarsest level and continue to analysis
        END IF
      END DO

      ! Check the which level fits the layout:
      IF (i .NE. mg_coarsest) THEN
        WRITE (*, 10)
10      FORMAT('#==============================================================================#')
        WRITE (*, 12) i - 1
12      FORMAT('# Geometry - !!!! Aborting the analysis at Glevel: ', I2, ' ......!!!!               #')
        WRITE (*, 11) (dimCell_global - 1) / 2**(mg_finest - i + 1) + 1, proc_layout(i, :) ! i is the acceptable level
11      FORMAT('# Geometry - Multigrid@ ', 2I6, ' cannot fit the multiprocess layout:', 2I3, ' #'/, &
               '# Check your num_grid and layout setting and rerun. May set fewer mpi_layout.  #')
        WRITE (*, 10)
        ! CALL FLUSH
        STOP
      END IF
    END IF
    CALL mpdd%barrier
    ! End of consistency check by Yuanfu Xie

    IF (mpdd%isBaseProc()) PRINT *, 'Geometry - Grid dim at finest level with an extra fictitious is: ', dimCell_global

    ! ALLOCATE the mpi_info in mpdd
    ALLOCATE (this%mg)

    CALL this%mg%initialize(mpdd, mpdd%group, mpdd%rank, proc_layout, t_steps_mg, start_time, end_time, mg_coarsest, mg_finest, vLevel, configFile)

    IF (mpdd%isBaseProc()) THEN
      IF ((dimCell_global(1) * dimCell_global(2)) .NE. mggen%mlgrid%Params(mggen%mlgrid%lvl_ends)%num_cell) THEN
        PRINT *, 'Geometry - Finest grid does not match! Check and rerun'
        STOP
      END IF
    END IF

    CALL mpdd%barrier
    ! Allocate and assign each multigrid level:
    DO i = this%mg%mg_finest, this%mg%mg_coarsest, -1
      ASSOCIATE (sg => this%mg%sg(i))

        IF (mpdd%isBaseProc()) sg%num_edge = mggen%mlgrid%Params(i)%numVrtxPerCell
        IF (mpdd%isBaseProc()) sg%num_vrtx = mggen%mlgrid%Params(i)%num_vrtx
        IF (mpdd%isBaseProc()) sg%num_stcl = mggen%mlgrid%Params(i)%numStclPerEdge
        IF (mpdd%isBaseProc()) sg%numQuadPerEdge = mggen%mlgrid%Params(i)%numQuadPerEdge
        IF (mpdd%isBaseProc()) sg%dimCell_global = (dimCell_global - 2) / (2**(mgGenFinest - i)) + 2
        IF (mpdd%isBaseProc()) sg%num_icell_global = sg%dimCell_global(1) * sg%dimCell_global(2)

        CALL mpdd%bCast(sg%num_edge)
        CALL mpdd%bCast(sg%num_vrtx)
        CALL mpdd%bCast(sg%num_stcl)
        CALL mpdd%bCast(sg%numQuadPerEdge)
        CALL mpdd%bCast(sg%dimCell_global)
        CALL mpdd%bCast(sg%num_icell_global)

        PRINT *, 'Start the precal on grid: ', i, 'on proc: ', mpdd%myrank

        CALL sg%mpddInfo_sg%calCellSize(sg%dimCell_global, sg%num_iCell)
        CALL mpdd%barrier

        IF (sg%mpddInfo_sg%isActiveProc()) THEN
          ALLOCATE (sg%edge_stcl(sg%num_stcl, sg%num_edge, sg%num_iCell))
          ALLOCATE (sg%sp_t_g_idx(sg%num_icell))
          ALLOCATE (sg%g_t_sp_idx(sg%num_icell_global))
          ALLOCATE (sg%seq_store_idx(sg%num_icell_global))
          ALLOCATE (sg%cell_stcl(sg%num_cell_stcl, sg%num_icell))
          sg%edge_stcl = 0
          sg%sp_t_g_idx = 0
          sg%g_t_sp_idx = 0
          sg%seq_store_idx = 0
        END IF

        PRINT *, 'End the calCellSize on grid: ', i, 'on proc: ', mpdd%myrank
        CALL mpdd%barrier

        ! Genearte map index in mpddInfo
        CALL sg%mpddInfo_sg%genMapIdx(sg%sp_t_g_idx, sg%g_t_sp_idx, sg%seq_store_idx, sg%dimCell_global, sg%dimCell)

        PRINT *, 'End the genMapIdx on grid: ', i, 'on proc: ', mpdd%myrank
        ! CALL mpdd%barrier

        ! Distribute cell stencil
        CALL sg%distGridInt(mggen%mlgrid%geoqty(i)%cell_stcl, sg%cell_stcl, [sg%num_cell_stcl, sg%num_icell_global])

        CALL sg%distGridInt(mggen%mlgrid%geoqty(i)%edge_stcl, &
                            sg%edge_stcl, [sg%num_stcl, sg%num_edge, sg%num_icell_global])

        PRINT *, 'End the distGridInt on grid: ', i, 'on proc: ', mpdd%myrank
        ! CALL mpdd%barrier
        ! CALL FLUSH()

        !! generate halo info before mapping to subprocess idx
        ! Genearte halo info in mpddInfo
        CALL sg%mpddInfo_sg%generate_halo_info(sg%cell_stcl, sg%num_iCell, sg%num_cell_stcl, sg%dimCell_global)

        ! Remap cell stencil to local sub proc index
        CALL sg%mpddInfo_sg%RemapAllGridIdxToSP2D(sg%cell_stcl, sg%g_t_sp_idx, sg%dimCell_global)

        ! Set the value of cell sizes: inner cell size and halo size
        sg%num_cell = sg%num_iCell + sg%mpddInfo_sg%halo_size
        sg%num_hcell = sg%mpddInfo_sg%halo_size

        ! Yuanfu Xie modified this print for more info for debugging: 2022-06-07:
        WRITE (*, 1) mpdd%myrank, i, sg%num_cell, sg%mpddInfo_sg%halo_size, sg%isactiveproc()
1       FORMAT('Geometry -- num_cell on proc: ', I2, ' at G level: ', I3, ' is: ', I8, ' with Halo size: ', I4, ' active: ', L)

        ! Allocate the variables in single grid.
        ALLOCATE ( &
          sg%cell_area(sg%num_cell), &
          sg%edge_leng(sg%num_edge, sg%num_cell), &
          sg%edgeNorm2(2, sg%num_edge, sg%num_cell), &
          sg%edgeTang2(2, sg%num_edge, sg%num_cell), &
          sg%coef_func(sg%num_stcl, sg%numQuadPerEdge, &
                       sg%num_edge, &
                       sg%num_cell), &
          sg%coef_norm(sg%num_stcl, sg%numQuadPerEdge, &
                       sg%num_edge, &
                       sg%num_cell), &
          sg%coef_tang(sg%num_stcl, sg%numQuadPerEdge, &
                       sg%num_edge, &
                       sg%num_cell), &
          sg%cell_cntr(2, sg%num_cell), &
          sg%cell_type(sg%num_cell), &
          sg%nonZeroCellStcl(sg%num_icell), &
          sg%cell_adnb(sg%num_edge, sg%num_icell), &
          sg%coef_gl(sg%numQuadPerEdge))  ! Gauss-Legendre quadrature weights

        ! Yuanfu Xie added coef_gl to distribute over sg grid for GL gquadrature weights
        IF (mpdd%isBaseProc()) sg%coef_gl = mggen%mlgrid%geoqty(i)%coef_gl
        CALL mpdd%bCast(sg%coef_gl)

        ! Yuanfu Xie added the following command to mimic cell_stcl to distribute cell_adnb:
        CALL sg%distGridInt(mggen%mlgrid%geoqty(i)%cell_adnb, sg%cell_adnb, [sg%num_edge, sg%num_icell_global])
        CALL this%mpdd%barrier

        ! Yuanfu Xie added the following to mimic cell_stcl in remapping cell_adnb:
        CALL sg%mpddInfo_sg%RemapAllGridIdxToSP2D(sg%cell_adnb, sg%g_t_sp_idx, sg%dimCell_global)

        ! Generate the cell type according to edge stencil
        CALL sg%mpddInfo_sg%genCellType(sg%cell_type, sg%edge_stcl, sg%num_icell, sg%sp_t_g_idx, sg%dimCell_global)

        ! Update halo values to 4 if it sits on domain boundaries: Yuanfu Xie added 2022/11/06
        CALL sg%boundaryHalo()

        ! Remap the stencil index to single proc index
        ! After adding fictitious point interpolation indices to the edge_stcl, we need to remap the indices, RemapAllGridIdxToSP2D, for
        ! all cell indices: Yuanfu Xie 2025-02-13
        ! CALL sg%mpddInfo_sg%RemapInnerGridIdxToSP3D(sg%edge_stcl, sg%g_t_sp_idx, sg%cell_type, sg%dimCell_global)
        CALL sg%mpddInfo_sg%RemapAllGridIdxToSP3D(sg%edge_stcl, sg%g_t_sp_idx, sg%dimCell_global)

        IF (sg%mpddInfo_sg%isBaseProc()) THEN
          ALLOCATE (sg%lat1DAtBase(sg%dimCell_global(1)), sg%lon1DAtBase(sg%dimCell_global(2)))
          CALL sg%get_lat_lon_vec1D(mggen%mlgrid%geoqty(i)%cell_cntr, sg%lat1DAtBase, sg%lon1DAtBase, &
                                    sg%dimCell_global(1), sg%dimCell_global(2))

          sg%maxLatGlob = MAXVAL(mggen%mlgrid%geoqty(mg_coarsest)%cell_cntr(1, :))
          sg%minLatGlob = MINVAL(mggen%mlgrid%geoqty(mg_coarsest)%cell_cntr(1, :))
          sg%maxLonGlob = MAXVAL(mggen%mlgrid%geoqty(mg_coarsest)%cell_cntr(2, :))
          sg%minLonGlob = MINVAL(mggen%mlgrid%geoqty(mg_coarsest)%cell_cntr(2, :))
        END IF

        CALL mpdd%bCast(sg%maxLatGlob)
        CALL mpdd%bCast(sg%minLatGlob)
        CALL mpdd%bCast(sg%maxLonGlob)
        CALL mpdd%bCast(sg%minLonGlob)

        CALL sg%DistGridRealWithHaloEx(mggen%mlgrid%geoqty(i)%cell_area, sg%cell_area, [sg%num_iCell_global])
        CALL sg%DistGridRealWithHaloEx(mggen%mlgrid%geoqty(i)%cell_cntr, sg%cell_cntr, [2, sg%num_iCell_global])
        CALL sg%DistGridRealWithHaloEx(mggen%mlgrid%geoqty(i)%edge_lnth, sg%edge_leng, [sg%num_edge, sg%num_iCell_global])
        CALL sg%DistGridRealWithHaloEx(mggen%mlgrid%geoqty(i)%norm_vct2, sg%edgeNorm2, [2, sg%num_edge, sg%num_iCell_global])
        CALL sg%DistGridRealWithHaloEx(mggen%mlgrid%geoqty(i)%tgnt_vct2, sg%edgeTang2, [2, sg%num_edge, sg%num_iCell_global])
        CALL sg%DistGridRealWithHaloEx(mggen%mlgrid%geoqty(i)%coef_func, sg%coef_func, &
                                       [sg%num_stcl, sg%numQuadPerEdge, sg%num_edge, sg%num_iCell_global])
        CALL sg%DistGridRealWithHaloEx(mggen%mlgrid%geoqty(i)%coef_norm, sg%coef_norm, &
                                       [sg%num_stcl, sg%numQuadPerEdge, sg%num_edge, sg%num_iCell_global])
        CALL sg%DistGridRealWithHaloEx(mggen%mlgrid%geoqty(i)%coef_tgnt, sg%coef_tang, &
                                       [sg%num_stcl, sg%numQuadPerEdge, sg%num_edge, sg%num_iCell_global])

        ! Yuanfu Xie moved this call from above to here as it adds sigma_x and sigma_y that needs edge_stcl: 2022-08-17
        PRINT *, 'Calling m_gen_cart_v_coord to initiate vertical coordinate: ', i
        CALL sg%m_gen_cart_v_coord

        ! Calculate the dlat and dlon
        CALL sg%m_cal_dlatlon

        BLOCK
          INTEGER(i_kind) :: j, k
          DO j = 1, sg%num_cell
            DO k = 1, UBOUND(sg%edgeNorm2, 2)
            IF (ISNAN(sg%edgeNorm2(1, k, j)) .OR. ISNAN(sg%edgeNorm2(2, k, j)) .OR. &
                ISNAN(sg%edgeTang2(1, k, j)) .OR. ISNAN(sg%edgeTang2(2, k, j))) THEN
              WRITE (*, 727) j, k, sg%edgeNorm2(1, k, j), sg%edgeNorm2(2, k, j), &
                sg%edgeTang2(1, k, j), sg%edgeTang2(2, k, j), sg%mpddGlob%myrank
727           FORMAT('Geometry - ', I6, I2' NT nan: ', 4D12.4, ' proc: ', I2)
            END IF
            END DO
          END DO
        END BLOCK

        PRINT *, 'Finish the precal on grid: ', i, 'on proc: ', mpdd%myrank
        ! IF (sg%mpddInfo_sg%isBaseProc()) print *, mggen%mlgrid%geoqty(i)%coef_norm

        PRINT *, 'Cell allocate 1: ', i, ALLOCATED(sg%cell_cntr), mpdd%myrank
        CALL mpdd%barrier

      END ASSOCIATE
    END DO

    DO i = this%mg%mg_coarsest, this%mg%mg_finest
      num_children = 4

      IF (i < this%mg%mg_finest) THEN
        IF (this%mg%sg(i + 1)%mpddInfo_sg%isActiveProc()) THEN
          CALL this%mg%sg(i + 1)%mpddInfo_sg%calCellSize(this%mg%sg(i)%dimcell_global, this%mg%sg(i)%num_icell_toFiner)
          ALLOCATE (this%mg%sg(i)%c_tc_idx(num_children, this%mg%sg(i)%num_icell_toFiner))
          ALLOCATE (this%mg%sg(i)%sp_t_g_idx_toFiner(this%mg%sg(i)%num_icell_toFiner))
          ALLOCATE (this%mg%sg(i)%g_t_sp_idx_toFiner(this%mg%sg(i)%num_icell_global))
          ALLOCATE (this%mg%sg(i)%seq_store_idx_toFiner(this%mg%sg(i)%num_icell_global))

          CALL this%mg%sg(i + 1)%mpddInfo_sg%genMapIdx(this%mg%sg(i)%sp_t_g_idx_toFiner, &
                                                       this%mg%sg(i)%g_t_sp_idx_toFiner, &
                                                       this%mg%sg(i)%seq_store_idx_toFiner, &
                                                       this%mg%sg(i)%dimCell_global)

          CALL this%mg%sg(i + 1)%mpddInfo_sg%DistGridIntSeq(this%mg%sg(i)%num_icell_toFiner, &
                                                            mggen%mlgrid%geoqty(i)%chld_indx, &
                                                            this%mg%sg(i)%c_tc_idx, &
                                                            [num_children, this%mg%sg(i)%dimCell_global], &
                                                            this%mg%sg(i)%seq_store_idx_toFiner)

          CALL this%mg%sg(i + 1)%mpddInfo_sg%RemapAllGridIdxToSP2D(this%mg%sg(i)%c_tc_idx, &
                                                                   this%mg%sg(i + 1)%g_t_sp_idx, &
                                                                   this%mg%sg(i + 1)%dimCell_global)
!PRINT *, 'Pre-cal of restriction operators to coarser grid: ', i, 'on proc: ', mpdd%myrank, 'is done!'
        END IF
      END IF

      IF (i > this%mg%mg_coarsest) THEN
        IF (this%mg%sg(i - 1)%mpddInfo_sg%isActiveProc()) THEN
          CALL this%mg%sg(i - 1)%mpddInfo_sg%calCellSize(this%mg%sg(i)%dimcell_global, this%mg%sg(i)%num_icell_toCoarser)
          ALLOCATE (this%mg%sg(i)%c2f_intp_stcl(num_children, this%mg%sg(i)%num_icell_toCoarser))
          ALLOCATE (this%mg%sg(i)%c2f_intp_cf(num_children, this%mg%sg(i)%num_icell_toCoarser))
          ALLOCATE (this%mg%sg(i)%sp_t_g_idx_toCoarser(this%mg%sg(i)%num_icell_toCoarser))
          ALLOCATE (this%mg%sg(i)%g_t_sp_idx_toCoarser(this%mg%sg(i)%num_icell_global))
          ALLOCATE (this%mg%sg(i)%seq_store_idx_toCoarser(this%mg%sg(i)%num_icell_global))

          CALL this%mg%sg(i - 1)%mpddInfo_sg%genMapIdx(this%mg%sg(i)%sp_t_g_idx_toCoarser, &
                                                       this%mg%sg(i)%g_t_sp_idx_toCoarser, &
                                                       this%mg%sg(i)%seq_store_idx_toCoarser, &
                                                       this%mg%sg(i)%dimCell_global)

          CALL this%mg%sg(i - 1)%mpddInfo_sg%DistGridIntSeq(this%mg%sg(i)%num_icell_toCoarser, &
                                                            mggen%mlgrid%geoqty(i)%cTof_stcl, &
                                                            this%mg%sg(i)%c2f_intp_stcl, &
                                                            [num_children, this%mg%sg(i)%dimCell_global], &
                                                            this%mg%sg(i)%seq_store_idx_toCoarser)

          CALL this%mg%sg(i - 1)%mpddInfo_sg%DistGridRealSeq(this%mg%sg(i)%num_icell_toCoarser, &
                                                             mggen%mlgrid%geoqty(i)%cTof_coef, &
                                                             this%mg%sg(i)%c2f_intp_cf, &
                                                             [num_children, this%mg%sg(i)%dimCell_global], &
                                                             this%mg%sg(i)%seq_store_idx_toCoarser)

          CALL this%mg%sg(i - 1)%mpddInfo_sg%RemapAllGridIdxToSP2D(this%mg%sg(i)%c2f_intp_stcl, &
                                                                   this%mg%sg(i - 1)%g_t_sp_idx, &
                                                                   this%mg%sg(i - 1)%dimCell_global)

!          print *, 'Pre-cal of interpretation operators to finer grid: ', i, 'on proc: ', mpdd%myrank, 'is done!'
        END IF
      END IF

      CALL mpdd%barrier
    END DO

    ! Initial the pointer in single grid
    DO i = this%mg%mg_coarsest, this%mg%mg_finest
      CALL boundingMG(this%mg, this%mg%sg(i)%mgParent)
    END DO

    BLOCK
      USE parameters_m, ONLY: degree2Radian, ae
      CHARACTER(LEN=1024) :: output_dir
      REAL(r_kind) :: max_dlon_dis, min_dlon_dis, temp
      INTEGER(i_kind) :: iostat

      IF (yaml_get_var(TRIM(configFile), 'IO', 'output_dir', output_dir) == 0) THEN
        output_dir = TRIM(output_dir)//'/AnalysisSummary.txt'

        IF (mpdd%isBaseProc()) THEN

          OPEN (10, FILE=output_dir, STATUS='REPLACE', FORM='FORMATTED', IOSTAT=iostat)
          ! Check the file is opened
          IF (iostat /= 0) THEN
            PRINT *, ' Error: The file ', TRIM(output_dir), ' is not opened!'
            STOP
          END IF
          WRITE (10, 40) mg_coarsest, mg_finest
        END IF

        DO i = mg_coarsest, mg_finest
          IF (this%mg%sg(i)%isActiveProc()) THEN
            temp = MAXVAL(this%mg%sg(i)%cell_dist(6, 1:this%mg%sg(i)%num_icell))
            CALL this%mg%sg(i)%mpddInfo_sg%AllReduceMaxReal(temp, max_dlon_dis)
            temp = MINVAL(this%mg%sg(i)%cell_dist(6, 1:this%mg%sg(i)%num_icell), MASK=this%mg%sg(i)%cell_dist(6, 1:this%mg%sg(i)%num_icell) .NE. 0.0)
            CALL this%mg%sg(i)%mpddInfo_sg%AllReduceMinReal(temp, min_dlon_dis)

            IF (mpdd%isBaseProc()) THEN

              WRITE (10, 400) i, proc_layout(i, 1), proc_layout(i, 2), proc_layout(i, 1) * proc_layout(i, 2), &
                this%mg%sg(i)%dimCell_global(1), this%mg%sg(i)%dimCell_global(2), &
                this%mg%sg(i)%dimCell_global(1) * this%mg%sg(i)%dimCell_global(2), &
                this%mg%sg(i)%dimCell_global(1) - 2, this%mg%sg(i)%dimCell_global(2) - 2, &
                (this%mg%sg(i)%dimCell_global(1) - 2) * (this%mg%sg(i)%dimCell_global(2) - 2), &
                (this%mg%sg(i)%cell_cntr(1, this%mg%sg(i)%cell_stcl(8, 1)) - this%mg%sg(i)%cell_cntr(1, 1)) / degree2Radian, &
                (this%mg%sg(i)%cell_cntr(1, this%mg%sg(i)%cell_stcl(8, 1)) - this%mg%sg(i)%cell_cntr(1, 1)) * ae / 1000.0, &
                (this%mg%sg(i)%cell_cntr(2, this%mg%sg(i)%cell_stcl(6, 1)) - this%mg%sg(i)%cell_cntr(2, 1)) / degree2Radian, &
                min_dlon_dis / 1000.0, max_dlon_dis / 1000.0
400           FORMAT('At level G', I3, ' --------------------------------------------', /, &
                     'The MPI parallel layout is ', I3, '   * ', I3, '.   Total proc is ', I3, /, &
                     'The cell grid is ', I5, '   * ', I5, '.   Total cell is ', I9, /, &
                     'The inner cell grid is ', I5, '   * ', I5, '.   Total inner cell is ', I9, /, &
                     'Latitude resolution is ', F9.3, '   =   ', F9.3, ' km', /, &
                     'Longitude resolution is ', F9.3, '   =   ', F9.3, ' km~', F9.3, ' km' &
                     )
            END IF

          END IF
        END DO
        IF (mpdd%isBaseProc()) CLOSE (10)
      END IF
    END BLOCK

    CALL mggen%destroy

    CALL mpdd%barrier
    DEALLOCATE (proc_layout)
    DEALLOCATE (t_steps_mg)
    PRINT *, 'Proc rank:', mpdd%myrank, 'Pre-cal of grid is done!'

  END SUBROUTINE

!> @brief
!! Destory the instance of geometry class.
  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(geometry_t), INTENT(INOUT) :: this

    CALL this%destroy
  END SUBROUTINE destructor

  SUBROUTINE destroy(this)
    CLASS(geometry_t), INTENT(INOUT) :: this
    IF (ALLOCATED(this%mg)) THEN
      DEALLOCATE (this%mg)
    END IF
    PRINT *, 'The geometry is destroyed on proc: ', this%mpdd%myrank, '!'
  END SUBROUTINE

END MODULE geometry_m
