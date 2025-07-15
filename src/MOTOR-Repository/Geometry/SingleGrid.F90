!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.MultiGrid.SingleGrid
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.1
! HISTORY           :
!   Created by Yuanfu Xie, 2019/04, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie 2020-4 for adding boundary conditions
!   Modified by Zilong Qin (zilong.qin@gmail.com), 2020/11/24, @GBA-MWF, Shenzhen, for parallelization
!   Modified by Yuanfu Xie, 2022/08, @GBA-MWF, Shenzhen
!     For adding edge normal and tangential vectors, topo_x and topo_y, and their memory handling.
!   Modified by Yuanfu Xie, 2022/09/04, @GBA-MWF, Shenzhen
!     For adding adjacent neighbor's indices for accessing them in Background Fillin codes for now.
!   Modified by Yuanfu Xie, 2022/09/10, @GBA-MWF, Shenzhen
!     For adding vertical interpolation coefficients.
!   Modified by Yuanfu Xie, 2022/10/18, @GBA-MWF, Shenzhen
!     For modifying sigmax and sigmay as a function of height and horizon.
!   Modified by Yuanfu Xie, 2022/10/20, @GBA-MWF, Shenzhen
!     For adding a horizontal smoother and boundary interpolation/extrapolation.
!   Modified by Yuanfu Xie and Zilong Qin, 2022/11/06, @GBA-MWF, Shenzhen
!     For adding a routine updating halo's cell_type to 4 if it sits on domain boundaries as request
!     by the boundaryOI routine.
!   Modified by Yuanfu Xie, 2023/10/18, @GBA-MWF, Shenzhen
!     Adding a new static parameter weighing the Laplacian operator in the place of background cov.
!   Modified by Yuanfu Xie, 2025/06/10, @GBA-MWF, Shenzhen
!     Adding a new routine to update pressure coordinate from topo and sigma, for the z-grid model.
!!--------------------------------------------------------------------------------------------------

!> @brief
!! Data structure for single grid structures.
!! @author Zilong Qin, Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
!! @note The structure is forge from the original poisson solver.
! @warning
! @attention
MODULE SingleGrid_m
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE mpddInfo_sg_m, ONLY: mpddInfo_sg_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE SeedGen_m, ONLY: GenSeed
  ! USE Geometry_m, ONLY: Geometry_t
  ! USE MultiGrid_m, ONLY: MultiGrid_t
  USE YAMLRead_m

  ! USE mpi_f08
  INCLUDE "mpif.h"

!> @brief
!! The data struct of single grid inside the multigrid.
!! @see multiGrid_t
! @note
! @warning
! @attention
  TYPE GridIdx_t
    INTEGER(i_kind) :: tIdx
    INTEGER(i_kind) :: hIdx
    INTEGER(i_kind) :: vIdx
  END TYPE GridIdx_t

  TYPE ::  SingleGrid_t
    ! This mult-grid container is prepare for either lat-lon grid or
    INTEGER(i_kind)                 :: num_cell = 0, &                !< Total cell number on the self proc & grid level
                                       num_edge = 0, &                !< Number of stencils per cell
                                       num_stcl = 0, &                !< Number of stencils per edge
                                       num_vrtx = 0, &                !< Number of stencils per vertex
                                       num_icell = 0, &               !< Number of inner cells on the self proc & grid level
                                       num_hcell = 0, &               !< Number of halo cells on the self proc & grid level
                                       numQuadPerEdge = 0

    INTEGER(i_kind), ALLOCATABLE    :: edge_stcl(:, :, :)         !< Index of stencil of each edge. Index: stcl/edge/cell
    INTEGER(i_kind), ALLOCATABLE    :: cell_adnb(:, :)            !< adjacent neighbors indices
    REAL(r_kind), ALLOCATABLE       :: edgeNorm2(:, :, :)         !< edge normal vector in 2D, n_vec2 in
    REAL(r_kind), ALLOCATABLE       :: edgeTang2(:, :, :)         !< edge tangential vector in 2D, t_vec2 in
    REAL(r_kind), ALLOCATABLE       :: cell_area(:)               !< Area of each cell. Index: cell
    REAL(r_kind), ALLOCATABLE       :: cell_cntr(:, :)            !< Latitude and Longitude of each cell, index: latlon/cell, in radian
    REAL(r_kind), ALLOCATABLE       :: edge_leng(:, :)            !< Length of each edge of each cell, index: edge/cell
    REAL(r_kind), ALLOCATABLE       :: coef_func(:, :, :, :), &   !< Stencil coefficents of value of each egde of each cell, index: stcl/1/edge/cell
                                       coef_norm(:, :, :, :), &   !< Normal derivative coefficents at each edge of each cell, index: stcl/1/edge/cell
                                       coef_tang(:, :, :, :)      !< Tangential derivative coefficents at each edge of each cell, index: stcl/1/edge/cell
    REAL(r_kind), ALLOCATABLE       :: coef_vrtx(:, :, :)         !< Vertex interpolation coefficients, 12/7/numVrtx; 7: vertex+6 GL
    REAL(r_kind), ALLOCATABLE       :: coef_gl(:)                 !< Gauss-Legrendre coefficients

    ! For multigrid prolongation and projection (restriction):
    INTEGER(i_kind), ALLOCATABLE    :: p_tc_idx(:), &             !< Child to parent cell indices, finer to coarser level
                                       c_tc_idx(:, :), &          !< Parent to child cell indices, coarser to finer level
                                       c2f_intp_stcl(:, :)        !< Coarser to finer prolongation stencil
    REAL(r_kind), ALLOCATABLE       :: c2f_intp_cf(:, :)          !< Coarser to finer interpolation coefficents

    ! Identifer of cell type
    INTEGER(i_kind), ALLOCATABLE    :: cell_type(:)               !< Type sign of cells: 0, inner, 1, boundary, 2, boundary corner, 3, halo
    INTEGER(i_kind), ALLOCATABLE    :: cell_stcl(:, :)            !< Stencil index of one cell (9, num_cell)
    INTEGER(i_kind), ALLOCATABLE    :: nonZeroCellStcl(:)         !< Number of non zero distance cell stcl
    INTEGER(i_kind), ALLOCATABLE    :: bdy_type(:)
    INTEGER(i_kind)                 :: halo_width                 !< halo width
    INTEGER(i_kind)                 :: num_cell_stcl              !< number of cell stencils

    ! Indices of edge integral: storing edge or GL point around a vertex per edge per cell
    !   1-dim has two numbers: 1st: vertex number; 2nd: edge number attached the vertex.
    !   2-dim has 4 numbers: 1 and 4 are the end points of an edge, 2 and 3 are GL
    !   3-dim has 3 numbers: for 3 edges of a given cell
    !   4-dim is the cell index
    ! INTEGER(i_kind), ALLOCATABLE    :: idx_edge_itpl(:, :, :, :)     !< edge integral indices
    INTEGER(i_kind), ALLOCATABLE    :: ivtx_stcl(:, :)             !< indices of vertex value interpolation
    INTEGER(i_kind), ALLOCATABLE    :: num_tri(:)                 !< number of edges or cells attached to a vertex

    ! Cell information
    INTEGER(i_kind)                 :: dimCell(2) = [0, 0]         !< Dimensions (2D, [row, col]/[DimNS, DimEW], including boundary, without halo)
    !! of the cell region of **self proc & grid level**
    INTEGER(i_kind)                 :: dimCell_global(2) = [0, 0]  !< Dimensions (2D, [row, col]/[DimNS, DimEW], including boundary, without halo)
    !! of the entire cell region on all procs on this grid level
    INTEGER(i_kind)                 :: num_icell_global = 0

    ! Cell index mapping
    INTEGER(i_kind), ALLOCATABLE    :: sp_t_g_idx(:)              !< single proc to global grid index
    INTEGER(i_kind), ALLOCATABLE    :: g_t_sp_idx(:)              !< global proc to single grid index
    TYPE(mpddInfo_sg_t), ALLOCATABLE  :: mpddInfo_sg

    ! Mapping for Distribution
    INTEGER(i_kind), ALLOCATABLE    :: seq_store_idx(:)

    ! For prolongation and restriction
    INTEGER(i_kind) :: num_icell_toCoarser = 0
    INTEGER(i_kind), ALLOCATABLE    :: sp_t_g_idx_toCoarser(:)
    INTEGER(i_kind), ALLOCATABLE    :: g_t_sp_idx_toCoarser(:)
    INTEGER(i_kind), ALLOCATABLE    :: seq_store_idx_toCoarser(:)

    !
    INTEGER(i_kind) :: num_icell_toFiner = 0
    INTEGER(i_kind), ALLOCATABLE    :: sp_t_g_idx_toFiner(:)
    INTEGER(i_kind), ALLOCATABLE    :: g_t_sp_idx_toFiner(:)
    INTEGER(i_kind), ALLOCATABLE    :: seq_store_idx_toFiner(:)

    ! The grid level index in the multi-grid of this single grid.
    INTEGER(i_kind) :: gLevel = 0
    INTEGER(i_kind) :: vLevel = 0
    INTEGER(i_kind) :: tSlots = 0

    ! Active grid cells in 2D (Vertical x Horizontal) Cartesian coordinate:
    LOGICAL, ALLOCATABLE :: aboveTopo(:, :)

    ! Difference dlat and dlon
    REAL(r_kind) :: dlat_dis, dlon_dis
    REAL(r_kind) :: dlat, dlon
    REAL(r_kind), ALLOCATABLE :: cell_dist(:, :)          !< (9, num_cell) The distance from the base cell to adjacent cells, currently only the 2\4\6\8 elements has been filled.

    ! Vertical coordinate
    REAL(r_kind), ALLOCATABLE :: zHght(:, :)              !< 3D real height of model grid, in dimension of (vLevel, numCell)
    REAL(r_kind), ALLOCATABLE :: pres(:, :)               !< 3D background pressure from IO, in dimension of (vLevel, numCell)
    CHARACTER(len=1024) :: hybridFile                       !< Hybrid file name Yuanfu Xie 2025-06-10
    REAL(r_kind), ALLOCATABLE :: hybridCoeff(:, :)        !< Hybrid coefficient, in dimension of (vLevel, 2) for a 2D hybrid coordinate

    REAL(r_kind), ALLOCATABLE :: sigma(:)                 !< Sigma coordinate, in dimension of (vLevel)
    REAL(r_kind), ALLOCATABLE :: coef_sigma(:, :)         !< New array for a second order interpolation of uneven vertical grid
    REAL(r_kind), ALLOCATABLE :: tt(:)                    !< Posix time, in dimension of (tSlots)
    REAL(r_kind), ALLOCATABLE :: topo(:)                  !< topography of the model field, in dimension of (num_cell)
    REAL(r_kind), ALLOCATABLE :: sigmax(:, :), sigmay(:, :)  !< Terrain following derivatives, (vLevel, numCell)
    REAL(r_kind) :: ptop, ztop, vertcalBreak              !< ptop is the top pressure, ztop is the top height, and vertcalBreak is the vertical break point for vertical similarity

    REAL(r_kind), ALLOCATABLE :: tskin(:, :)          !< Skin temperature, in dimension of (num_cell)
    REAL(r_kind), ALLOCATABLE :: psfc(:, :)           !< Surface pressure, in dimension of (num_cell)
    REAL(r_kind), ALLOCATABLE :: u10m(:, :)           !< U wind at 10 m, in dimension of (num_cell)
    REAL(r_kind), ALLOCATABLE :: v10m(:, :)           !< V wind at 10 m, in dimension of (num_cell)
    REAL(r_kind), ALLOCATABLE :: landmask(:)       !< Sea/land (0/1) of the topography, in dimension of (num_cell)
    REAL(r_kind), ALLOCATABLE :: soil_type(:)      !< 16-category soil type, in dimension of (num_cell)
    REAL(r_kind), ALLOCATABLE :: veg_frac(:)       !< Monthly green fraction, in dimension of (num_cell)
    REAL(r_kind), ALLOCATABLE :: snowc(:, :)          !< snow cover, in dimension of (num_cell)
    REAL(r_kind), ALLOCATABLE :: xice(:, :)           !< xice cover , in dimension of (num_cell)
    ! Vertical similarity: parameter with 0 at the bottom 1 top, weighing horizonalSimilarity 100% at bottom, 0 top,
    ! with tanh center at vertcalBreak input from yaml, Yuanfu Xie 2023/10/21
    LOGICAL :: appliedSimilarity, flowDependentThinning ! Yuanfu Xie added landmask/topo similarity and flowDependent thinning options
    REAL(r_kind), ALLOCATABLE :: scales4Similarity(:)
    REAL(r_kind), ALLOCATABLE :: vertcalSimilarity(:)
    REAL(r_kind), ALLOCATABLE :: horizonSimilarity(:) !< Horizonal similarity: including landmask, top etc.
    !< If it is similar to its neighbor, it is 1; otherwise decays exponentially.

    ! Scaling factor with respected to log and height
    REAL(r_kind), ALLOCATABLE :: s1(:, :)
    REAL(r_kind), ALLOCATABLE :: SVapor(:, :), STemp(:, :), SVapor1D(:), STemp1D(:)
    REAL(r_kind), ALLOCATABLE :: SCloud(:, :), SIce(:, :), SRain(:, :), SSnow(:, :), SGraupel(:, :)
    REAL(r_kind), ALLOCATABLE :: FGPres(:, :, :)
    REAL(r_kind), ALLOCATABLE :: s_qvapor(:)

    TYPE(mpddGlob_t), POINTER :: mpddGlob        !< Pointer of global multiprocessing domain decomposition class
    INTEGER(i_kind) :: id

    CHARACTER(len=1024) :: configFile

    REAL(r_kind), ALLOCATABLE :: lat1DAtBase(:), lon1DAtBase(:)
    REAL(r_kind) :: maxLatGlob, maxLonGlob, minLatGlob, minLonGlob

    ! for z-grid model
    INTEGER(i_kind), ALLOCATABLE :: BufferHalo(:), & !<bdy buffer halo where Tendency is calculated>
                                    ks_interp(:, :)
    REAL(r_kind), ALLOCATABLE :: coef_fstdif(:, :, :), &  !< vertical coefs for first order derivative cal. in the full layer
                                 coef_scddif(:, :, :), &       !< vertical coefs for second order derivative in the full layer
                                 parz_parsigma(:, :), &        !< partial z partial sigma in full layer
                                 Hz(:, :), &        !< partial z partial sigma in full layer
                                 parHz_parsigma(:, :), &        !< partial z partial sigma in full layer
                                 Lz(:, :), & !< From here and below, the operators are the gzm operators
                                 F_z_z(:, :), & !<
                                 F_Hz_z(:, :), &
                                 F_invHz_z(:, :), &
                                 InterpCoef(:, :, :), &
                                 coef_fstdif_half(:, :, :), & !<>
                                 f(:, :), &
                                 BdyAlpha(:), BdyBeta(:) !<bdy tendency weight, alpha for bdy, beta for inner>
    REAL(r_kind) :: TMSL ! Mean sea level temperature
    CLASS(*), POINTER :: mgParent
    ! TYPE(MultiGrid_t), POINTER :: mg
  CONTAINS

    FINAL :: destructor
    PROCEDURE, PUBLIC :: initialize

!> @brief
!! Distribute the INTEGER(i_kind) buffer from the proc 0 to all procs on this grid.
    GENERIC, PUBLIC :: DistGridInt => DistGridInt3D, &
      DistGridInt2D, &
      DistGridInt1D, &
      DistGridInt4D
    PROCEDURE, PUBLIC :: DistGridInt4D
    PROCEDURE, PUBLIC :: DistGridInt3D
    PROCEDURE, PUBLIC :: DistGridInt2D
    PROCEDURE, PUBLIC :: DistGridInt1D

!> @brief
!! Distribute the REAL(r_kind) buffer from the proc 0 to all procs on this grid.
    GENERIC, PUBLIC :: DistGridReal => DistGridReal3D, &
      DistGridReal2D, &
      DistGridReal1D, &
      DistGridReal4D
    PROCEDURE, PUBLIC :: DistGridReal4D
    PROCEDURE, PUBLIC :: DistGridReal3D
    PROCEDURE, PUBLIC :: DistGridReal2D
    PROCEDURE, PUBLIC :: DistGridReal1D

!> @brief
!! Distribute the REAL(r_kind) buffer from the proc 0 to all procs on this grid, and update halos.
    GENERIC, PUBLIC :: DistGridRealWithHaloEx => DistGridRealWithHaloEx1D, &
      DistGridRealWithHaloEx2D, &
      DistGridRealWithHaloEx3D, &
      DistGridRealWithHaloEx4D
    PROCEDURE, PUBLIC :: DistGridRealWithHaloEx4D
    PROCEDURE, PUBLIC :: DistGridRealWithHaloEx3D
    PROCEDURE, PUBLIC :: DistGridRealWithHaloEx2D
    PROCEDURE, PUBLIC :: DistGridRealWithHaloEx1D

!> @brief
!! Fill the points on the boundary using linear interpolation.
    PROCEDURE, PUBLIC :: m_interp_points_on_bndy_linear

    GENERIC, PUBLIC :: aggrGridReal => aggrGridReal2D, aggrGridReal1D
    PROCEDURE, PUBLIC :: aggrGridReal2D, aggrGridReal1D
    ! PROCEDURE, PUBLIC :: aggrGridReal3D

    PROCEDURE, PUBLIC :: ExchangeMatOnHalo2D
    PROCEDURE, PUBLIC :: ExchangeMatOnHaloReverseSum2D

    ! These functions are designed specifically for manipulating state variables, for which the t index is set as the last dimension.
    PROCEDURE, PUBLIC :: DistGridRealForFieldGrid
    PROCEDURE, PUBLIC :: DistGridRealWithHaloExForFieldGrid
    PROCEDURE, PUBLIC :: DistGridRealWithHaloExForFieldGrid4D
    PROCEDURE, PUBLIC :: aggrGridRealForFieldGrid
    PROCEDURE, PUBLIC :: ExchangeMatOnHaloForFieldGrid
    PROCEDURE, PUBLIC :: ExchangeMatOnHaloReverseSumForFieldGrid

    PROCEDURE, PUBLIC :: isActiveProc
    PROCEDURE, PUBLIC :: isBaseProc
    PROCEDURE, PUBLIC :: allocateMat
    PROCEDURE, PUBLIC :: allocateMat3D          !added by Yu Xin
    PROCEDURE, PUBLIC :: allocateMatGlobalOnBase
    PROCEDURE, PUBLIC :: allocateMatGlobalOnBase3D   !added by Yu Xin
    PROCEDURE, PUBLIC :: m_cal_dlatlon
    PROCEDURE, PUBLIC :: m_gen_cart_v_coord
    PROCEDURE, PUBLIC, NOPASS :: get_lat_lon_vec1D
    PROCEDURE :: update_zHght_from_topo_and_sigma, update_press_from_topo_and_sigma !<-- Yuanfu Xie added press coordinate 2025-06-10
    PROCEDURE, PRIVATE, PASS(this) :: initialize_surface_parameter

    ! Horizontal smoother and boundary interpolation/extrapolation: Yuanfu Xie added 2022-10-20
    PROCEDURE, PUBLIC :: boundaryOI  ! Boundary cell value by an OI method
    PROCEDURE, PUBLIC :: boundaryHalo   ! Update cell_type = 4 for all halo sitting on domain boundaries
  END TYPE SingleGrid_t

CONTAINS

  SUBROUTINE get_lat_lon_vec1D(cell_cntr, lat, lon, latSize, lonSize)
    USE parameters_m, ONLY: degree2radian

    IMPLICIT NONE
    REAL(r_kind), INTENT(IN) :: cell_cntr(2, latSize * lonSize)
    REAL(r_kind), INTENT(OUT) :: lat(latSize), lon(lonSize)
    INTEGER(i_kind), INTENT(IN) :: latSize, lonSize
    INTEGER(i_kind) :: i, j, k

    DO i = 1, latSize
      lat(i) = cell_cntr(1, (i - 1) * lonSize + 1) / degree2radian
    END DO

    DO i = 1, lonSize
      lon(i) = cell_cntr(2, i) / degree2radian
    END DO
    ! print *, cell_cntr(:, 1), lat(1), lon(1)

  END SUBROUTINE get_lat_lon_vec1D

  SUBROUTINE m_gen_cart_v_coord(this)
    CLASS(SingleGrid_t) :: this
    INTEGER :: i

    ALLOCATE (this%sigma(this%vLevel), &
              this%vertcalSimilarity(this%vLevel), &
              this%zHght(this%vLevel, this%num_cell), &
              this%aboveTopo(this%vLevel, this%num_cell), &
              this%coef_sigma(3, this%vLevel), &  ! 3 coefs for a second order interpolation
              this%topo(this%num_cell), &
              this%pres(this%vLevel, this%num_cell), &
              this%sigmax(this%vLevel, this%num_cell), &
              this%sigmay(this%vLevel, this%num_cell), &
              this%s1(this%vLevel, this%num_cell), &
              this%SVapor(this%vLevel, this%num_cell), &
              this%SCloud(this%vLevel, this%num_cell), &
              this%SIce(this%vLevel, this%num_cell), &
              this%SRain(this%vLevel, this%num_cell), &
              this%SSnow(this%vLevel, this%num_cell), &
              this%SGraupel(this%vLevel, this%num_cell), &
              this%STemp(this%vLevel, this%num_cell), &
              this%SVapor1D(this%vLevel), &
              this%STemp1D(this%vLevel), &
              this%FGPres(this%vLevel, this%num_cell, this%tSlots), &
              this%s_qvapor(this%vLevel))

    IF (this%ztop < 10.0D0) this%ztop = 30000.00D0
    IF (this%vLevel > 1) THEN
      DO i = 1, this%vLevel
        this%sigma(i) = (i - 1) * 1.0D0 / (this%vLevel - 1) * this%ztop
        ! Yuanfu Xie added the vertcalSimilarity values:
        this%vertcalSimilarity(i) = (this%sigma(i) / this%ztop - this%vertcalBreak) / this%vertcalBreak * 6.0D0 ! Temporarily set curvature ratio to 6
        this%vertcalSimilarity(i) = (EXP(this%vertcalSimilarity(i)) - EXP(-this%vertcalSimilarity(i))) / &
                                    (EXP(this%vertcalSimilarity(i)) + EXP(-this%vertcalSimilarity(i)))
        this%vertcalSimilarity(i) = (1.0D0 + this%vertcalSimilarity(i)) / 2.0D0
      END DO
    ELSE
      this%sigma = 0.0D0; 
    END IF

    ! Yuanfu Xie added an initialization to zHght: 2022-08-22:
    this%zHght = 0.0D0
    this%coef_sigma = 0.0D0
    this%sigmax = 0.0D0 ! initialize
    this%sigmay = 0.0D0

    this%topo = 0.0D0
    this%s_qvapor = 1.0D0

    CALL this%update_zHght_from_topo_and_sigma()

    ! Initialize the pressure coordinate from topo and sigma:
    this%pres = 1.0D5 !< Default pressure is 1000 hPa
    CALL this%update_press_from_topo_and_sigma(this%pres(1,:)) !<-- Yuanfu Xie added press coordinate 2025-06-17

    CALL this%initialize_surface_parameter()
  END SUBROUTINE m_gen_cart_v_coord

  SUBROUTINE initialize_surface_parameter(this)
    CLASS(SingleGrid_t) :: this

    ALLOCATE (this%tskin(this%num_cell, this%tSlots), &
              this%psfc(this%num_cell, this%tSlots), &
              this%u10m(this%num_cell, this%tSlots), &
              this%v10m(this%num_cell, this%tSlots), &
              this%landmask(this%num_cell), &
              this%soil_type(this%num_cell), &
              this%veg_frac(this%num_cell), &
              this%snowc(this%num_cell, this%tSlots), &
              this%xice(this%num_cell, this%tSlots), &
              this%horizonSimilarity(this%num_cell) &
              )

    this%tskin = 0.0D0
    this%psfc = 0.0D0
    this%u10m = 0.0D0
    this%v10m = 0.0D0
    this%landMask = 1
    this%soil_type = 0.0D0

    this%veg_frac = 0.0D0
    this%snowc = 0.0D0
    this%xice = 0.0D0
    this%vertcalSimilarity = 1.0D0
    this%horizonSimilarity = 1.0D0  ! Default the current cell has identical similarity to its neighbors.
  END SUBROUTINE initialize_surface_parameter

  !>@brief
  !! Yuanfu Xie modified this routine for adding topo derivatives: 2022-08-17
  SUBROUTINE update_zHght_from_topo_and_sigma(this)
    IMPLICIT NONE
    CLASS(SingleGrid_t) :: this

    ! Local variables:
    INTEGER(i_kind) :: i, j, iedge, istcl
    REAL(r_kind) :: norm, tang

    INTEGER(i_kind) :: istatus
    CHARACTER(len=1024) :: VerticalCoord

    VerticalCoord = 'TerrainFollowing'
    istatus = yaml_get_var(TRIM(this%configFile), 'DASpace', 'VerticalCoord', VerticalCoord)

    PRINT *, 'MAX TOPO: ', MAXVAL(this%topo), this%sigma(this%vLevel), this%gLevel, this%vLevel

    IF (this%vLevel > 1) THEN
      this%ztop = this%sigma(this%vLevel)

      IF (VerticalCoord == 'TerrainFollowing') THEN
        DO i = 1, this%vLevel
          this%zHght(i, :) = this%sigma(i) * (this%ztop - this%topo) / this%ztop + this%topo
        END DO
      ELSE IF (VerticalCoord == 'Cartesian') THEN
        DO i = 1, this%vLevel
          this%zHght(i, :) = this%sigma(i) * (this%ztop - 0.0D0) / this%ztop + 0.0D0
        END DO
      END IF

      ASSOCIATE (h => this%sigma)
        DO i = 2, this%vLevel - 1
          this%coef_sigma(1, i) = (h(i + 1) - h(i)) / (h(i - 1) - h(i)) / (h(i + 1) - h(i - 1))
          this%coef_sigma(2, i) = (h(i + 1) + h(i - 1) - 2.0D0 * h(i)) / (h(i + 1) - h(i)) / (h(i) - h(i - 1))
          this%coef_sigma(3, i) = (h(i) - h(i - 1)) / (h(i + 1) - h(i)) / (h(i + 1) - h(i - 1))
        END DO
      END ASSOCIATE
    ELSE
      this%zHght(1, :) = this%topo
      this%sigma = 0.0D0
    END IF

    ! Calculate sigma_x and sigma_y: by Yuanfu Xie 2022-08-17
    WRITE (*, 2) UBOUND(this%topo, 1), MAXVAL(ABS(this%topo)), this%ztop, this%gLevel, this%mpddGlob%myrank
2   FORMAT('Update-zHght/topo/sigma: ', I8, ' Max topo: ', D12.4, ' Ztop: ', D12.4, ' at G: ', I2, ' proc: ', I2)
    this%sigmax = 0.0D0
    this%sigmay = 0.0D0
    DO i = 1, this%num_icell
      DO iedge = 1, UBOUND(this%edge_stcl, 2)
        norm = 0.0D0
        tang = 0.0D0
        DO istcl = 1, UBOUND(this%edge_stcl, 1)
          IF (this%edge_stcl(istcl, iedge, i) .GT. 0) THEN

            IF (this%edge_stcl(istcl, iedge, i) .GT. this%num_cell) THEN
              WRITE (*, 4545) this%edge_stcl(istcl, iedge, i), this%num_cell, &
                istcl, iedge, i, this%gLevel, this%mpddGlob%myrank
4545          FORMAT('Exceeding edge_stcl: ', 2I6, ' stcl:', I3, ' edge:', I2, ' cell:', I7, ' G', I1, ' pc: ', I1)
              STOP
            END IF

            norm = norm + &
                   this%coef_norm(istcl, 1, iedge, i) * this%topo(this%edge_stcl(istcl, iedge, i))
            tang = tang + &
                   this%coef_tang(istcl, 1, iedge, i) * this%topo(this%edge_stcl(istcl, iedge, i))
          END IF
        END DO
        IF (ISNAN(this%edgeNorm2(1, iedge, i)) .OR. ISNAN(this%edgeTang2(1, iedge, i)) .OR. &
            ISNAN(this%edgeNorm2(2, iedge, i)) .OR. ISNAN(this%edgeTang2(2, iedge, i)) .OR. &
            ISNAN(norm) .OR. ISNAN(tang)) THEN
          WRITE (*, 5) i, iedge, this%edgeNorm2(:, iedge, i), this%edgeTang2(:, iedge, i), &
            norm, tang, this%mpddGlob%myrank
5         FORMAT('update_zHght_from_topo_and_sigma: ', I6, I2, ' edgeNT: ', 4D12.4, ' nt: ', 2D12.4, ' proc ', I2)
        END IF
        this%sigmax(:, i) = this%sigmax(:, i) + &
                            this%edgeNorm2(1, iedge, i) * norm + this%edgeTang2(1, iedge, i) * tang
        this%sigmay(:, i) = this%sigmay(:, i) + &
                            this%edgeNorm2(2, iedge, i) * norm + this%edgeTang2(2, iedge, i) * tang
      END DO
      IF (this%ztop - this%topo(i) .LE. 0.0D0) THEN
        WRITE (*, 1) i, this%ztop - this%topo(i), this%mpddInfo_sg%myrank
1       FORMAT('update_zHght_from_topo_and_sigma - Ztop <= topo: ', D12.4, /, &
               'Check your ztop and topo and rerun! ')
        STOP
      END IF
      DO j = 1, this%vLevel
        this%sigmax(j, i) = this%ztop * (this%zHght(j, i) - this%ztop) * this%sigmax(j, i) / &
                            DBLE(UBOUND(this%edge_stcl, 2)) / (this%ztop - this%topo(i))**2
        this%sigmay(j, i) = this%ztop * (this%zHght(j, i) - this%ztop) * this%sigmay(j, i) / &
                            DBLE(UBOUND(this%edge_stcl, 2)) / (this%ztop - this%topo(i))**2

        IF (ISNAN(this%sigmax(j, i))) PRINT *, 'Sigmax: ', i, this%sigmax(j, i), this%mpddInfo_sg%myrank, this%gLevel
        IF (ISNAN(this%sigmay(j, i))) PRINT *, 'Sigmay: ', i, this%sigmay(j, i), this%mpddInfo_sg%myrank, this%gLevel
      END DO
    END DO

    this%s1 = 0.001 * EXP(-0.0001 * (this%zHght - 400.0))**(this%zHght / 35000.0 + 1)
    this%s_qvapor = 0.001 * EXP(-0.0001 * (this%sigma - 400.0)) ! sv2 TODO:
    ! this%s_qvapor = 0.001 ! sv2 TODO:
    this%sTemp = 1.0D0 ! TODO:
    this%SVapor = 10**(2.4 + (-20 / (LOG10(100000 * (1 - this%zHght / 44300)**8.256)))) ! wait for a smoothed terrain on coarse grids
    FORALL (i=1:this%num_cell) this%SCloud(:, i) = this%s_qvapor
    FORALL (i=1:this%num_cell) this%SIce(:, i) = this%s_qvapor
    FORALL (i=1:this%num_cell) this%SRain(:, i) = this%s_qvapor
    FORALL (i=1:this%num_cell) this%SSnow(:, i) = this%s_qvapor
    FORALL (i=1:this%num_cell) this%SGraupel(:, i) = this%s_qvapor
    ! this%SCloud = 0.0001D0
    ! this%SIce = 0.00001D0
    ! this%SCloud = 1.0D0
    ! ! this%SIce = 1.0D0
    ! this%SRain = 0.001D0
    ! this%SSnow = 0.001D0
    ! this%SGraupel = 0.001D0
    WHERE (this%SVapor < 1.0E-12) this%SVapor = 1.0E-12

    ! PRINT *, 'check sv: ', MAXVAL(this%SVapor), MINVAL(this%SVapor)
    ! SV: minval - 10^-6; maxval - 10^-2
  END SUBROUTINE update_zHght_from_topo_and_sigma

  !>@brief
  !! Yuanfu Xie modified this routine for adding topo derivatives: 2022-08-17
  SUBROUTINE update_press_from_topo_and_sigma(this, psfc)
    IMPLICIT NONE
    CLASS(SingleGrid_t) :: this
    REAL(r_kind), INTENT(IN) :: psfc(this%num_cell)

    INTEGER(i_kind) :: istatus, iv
    CHARACTER(len=1024) :: VerticalCoord

    VerticalCoord = 'TerrainFollowing'
    istatus = yaml_get_var(TRIM(this%configFile), 'DASpace', 'VerticalCoord', VerticalCoord)
    IF (istatus .NE. 0) THEN
      PRINT *, 'Error: VerticalCoord not set in the config file.'
      RETURN
    END IF

    ! Calculate the pressure and the sigmap:
    IF (this%ptop .GT. 0.0D0) THEN
      DO iv=1,this%vLevel
        this%pres(iv, :) = this%hybridCoeff(iv, 1) + this%hybridCoeff(iv, 2) * psfc(:)   !<-- Yuanfu Xie added press coordinate 2025-06-10
      END DO
      PRINT*,'Update hybrid pressure from hybrid coeffs: ', MAXVAL(this%pres), MINVAL(this%pres), this%gLevel, this%vLevel
    END IF
  END SUBROUTINE update_press_from_topo_and_sigma

  SUBROUTINE m_cal_dlatlon(this)
    USE parameters_m, ONLY: EarthRadius
    USE geoTools_m, ONLY: distance
    IMPLICIT NONE

    CLASS(singleGrid_t) :: this
    INTEGER :: i, j

    IF (.NOT. this%isActiveProc()) RETURN

    IF (ALLOCATED(this%cell_dist)) DEALLOCATE (this%cell_dist)
    ALLOCATE (this%cell_dist((this%halo_width * 2 + 1)**2, this%num_icell))

    this%cell_dist = 0.0D0
    this%dlat_dis = (this%cell_cntr(1, this%cell_stcl(8, 1)) - this%cell_cntr(1, 1)) * EarthRadius
    this%dlon_dis = ABS(this%cell_cntr(2, this%cell_stcl(6, 1)) - this%cell_cntr(2, 1)) * EarthRadius * DCOS(this%cell_cntr(1, 1))

    this%dlat = this%cell_cntr(1, this%cell_stcl(8, 1)) - this%cell_cntr(1, 1)
    this%dlon = this%cell_cntr(2, this%cell_stcl(6, 1)) - this%cell_cntr(2, 1)

    ! Yuanfu Xie modified this%cell_dist to all cells:
    DO i = 1, this%num_icell
      this%nonZeroCellStcl(i) = 0
      DO j = 1, UBOUND(this%cell_stcl, 1)
        IF (this%cell_stcl(j, i) .GT. 0) THEN
          this%nonZeroCellStcl(i) = this%nonZeroCellStcl(i) + 1
          this%cell_dist(j, i) = distance(this%cell_cntr(:, i), this%cell_cntr(:, this%cell_stcl(j, i)))
        END IF
      END DO
    END DO

  END SUBROUTINE m_cal_dlatlon

  SUBROUTINE m_interp_points_on_bndy_linear(this, vLevel, bufr)
    IMPLICIT NONE
    CLASS(SingleGrid_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vLevel
    REAL(r_kind), INTENT(INOUT) :: bufr(vLevel, this%num_cell)
    INTEGER(i_kind) :: i

    IF (this%num_iCell_global .EQ. 9) THEN
      PRINT *, 'No to interpolate the single cell grid.'
      RETURN
    END IF

    DO i = 1, this%num_icell
      IF ((this%cell_type(i) .EQ. 1)) THEN
        IF (this%cell_stcl(2, i) .EQ. 0) THEN
          bufr(:, i) = 2 * bufr(:, this%cell_stcl(8, i)) - bufr(:, this%cell_stcl(8, this%cell_stcl(8, i)))
          CYCLE
        ELSE IF (this%cell_stcl(6, i) .EQ. 0) THEN
          bufr(:, i) = 2 * bufr(:, this%cell_stcl(4, i)) - bufr(:, this%cell_stcl(4, this%cell_stcl(4, i)))
          CYCLE
        ELSE IF (this%cell_stcl(8, i) .EQ. 0) THEN
          bufr(:, i) = 2 * bufr(:, this%cell_stcl(2, i)) - bufr(:, this%cell_stcl(2, this%cell_stcl(2, i)))
          CYCLE
        ELSE IF (this%cell_stcl(4, i) .EQ. 0) THEN
          bufr(:, i) = 2 * bufr(:, this%cell_stcl(6, i)) - bufr(:, this%cell_stcl(6, this%cell_stcl(6, i)))
          CYCLE
        END IF
      ELSE IF (this%cell_type(i) .EQ. 2) THEN
        IF ((this%cell_stcl(7, i) .EQ. 0) .AND. (this%cell_stcl(1, i) .EQ. 0) .AND. (this%cell_stcl(3, i) .EQ. 0)) THEN
          bufr(:, i) = 2 * bufr(:, this%cell_stcl(9, i)) - bufr(:, this%cell_stcl(9, this%cell_stcl(9, i)))
          CYCLE
        ELSE IF ((this%cell_stcl(1, i) .EQ. 0) .AND. (this%cell_stcl(3, i) .EQ. 0) .AND. (this%cell_stcl(9, i) .EQ. 0)) THEN
          bufr(:, i) = 2 * bufr(:, this%cell_stcl(7, i)) - bufr(:, this%cell_stcl(7, this%cell_stcl(7, i)))
          CYCLE
        ELSE IF ((this%cell_stcl(3, i) .EQ. 0) .AND. (this%cell_stcl(9, i) .EQ. 0) .AND. (this%cell_stcl(7, i) .EQ. 0)) THEN
          bufr(:, i) = 2 * bufr(:, this%cell_stcl(1, i)) - bufr(:, this%cell_stcl(1, this%cell_stcl(1, i)))
          CYCLE
        ELSE IF ((this%cell_stcl(9, i) .EQ. 0) .AND. (this%cell_stcl(7, i) .EQ. 0) .AND. (this%cell_stcl(1, i) .EQ. 0)) THEN
          bufr(:, i) = 2 * bufr(:, this%cell_stcl(3, i)) - bufr(:, this%cell_stcl(3, this%cell_stcl(3, i)))
          CYCLE
        END IF
      END IF
    END DO

  END SUBROUTINE m_interp_points_on_bndy_linear

  SUBROUTINE initialize(this, mpddGlob, gLevel, vLevel, group, rank, proc_layout, tSlots, start_time, end_time, configFile)
    !USE NMLRead_m

    IMPLICIT NONE
    CLASS(SingleGrid_t) :: this
    TYPE(mpddGlob_t), TARGET, INTENT(IN) :: mpddGlob        !< Global multiprocessing domain decomposition class pointer
    INTEGER(i_kind), INTENT(IN) :: group, gLevel, vLevel
    INTEGER(i_kind), INTENT(IN) :: rank, proc_layout(2), tSlots
    REAL(r_kind), INTENT(IN) :: start_time, end_time
    INTEGER(i_kind) :: i, istatus
    CHARACTER(len=1024) configFile, staticDir

    REAL(r_kind), ALLOCATABLE :: coeff(:)

    PRINT *, 'Begin construction of singleGrid on G', gLevel

    this%mpddGlob => mpddGlob
    this%num_cell = 0
    this%num_edge = 0
    this%num_stcl = 0
    this%num_vrtx = 0
    this%num_icell = 0
    this%num_hcell = 0
    this%numQuadPerEdge = 0

    this%gLevel = gLevel
    this%vLevel = vLevel
    this%tSlots = tSlots

    ALLOCATE (this%tt(this%tSlots))
    IF (this%tSlots > 1) THEN
      FORALL (i=1:this%tSlots) this%tt(i) = start_time + (end_time - start_time) / (this%tSlots - 1) * (i - 1)
    ELSE
      this%tt(1) = start_time ! Single time slot analysis
    END IF

    ! PRINT *,'tt', start_time, end_time, this%tt
    this%id = GenSeed()

    ! IF (ALLOCATED(this%mpddInfo_sg)) DEALLOCATE (this%mpddInfo_sg)
    ALLOCATE (this%mpddInfo_sg)

    CALL this%mpddInfo_sg%initializemMPDDSub(this%gLevel, proc_layout, group, rank)

    this%configFile = configFile

    !PJM CALL namelist_read(configFile, 'ztop', this%ztop)
    this%ztop = 30000.00D0
    IF (this%vLevel .NE. 1) THEN
      istatus = yaml_get_var(TRIM(configFile), 'DASpace', 'ztop', this%ztop)

      ! Yuanfu Xie added a hybrid vertical coordinate option: 2025-06-17
      istatus = yaml_get_var(TRIM(configFile), 'DASpace', 'ptop', this%ptop)
      IF (istatus .EQ. 0) THEN
        WRITE(*,7) this%ptop, this%gLevel, this%mpddGlob%myrank
7       FORMAT('SingleGrid: ptop is set to value ', D12.4, ' at G', I2, ' proc: ', I2)
        istatus = yaml_get_var(TRIM(configFile), 'DASpace', 'hybridFile', this%hybridFile)
        IF (istatus .NE. 0) THEN
          WRITE (*, 9) istatus
9         FORMAT('SingleGrid: failed reading in hybridFile: ', I4)
        ELSE
          CALL GET_ENVIRONMENT_VARIABLE('STATIC_DIR',staticDir)
          PRINT*,'SingleGrid: staticDir is set to ', TRIM(staticDir)//'/'//TRIM(this%hybridFile)
          WRITE (*, 10) TRIM(this%hybridFile),TRIM(this%configFile)
10        FORMAT('SingleGrid: hybridFile is set to value ', 2A)
        END IF
        ! Read in hybrid coefficients:
        ALLOCATE (this%hybridCoeff(this%vLevel, 2))
        istatus = yaml_get_var(TRIM(staticDir)//'/'//TRIM(this%hybridFile), 'hybrid_levels', 'a', coeff)
        IF (UBOUND(coeff, 1) .NE. this%vLevel) THEN
          WRITE (*, 6) UBOUND(coeff, 1), this%vLevel
6         FORMAT('SingleGrid: hybridCoeff size mismatch: ', I4, ' vs. ', I4)
          ! Set istatus to 1 to indicate a problem
          istatus = 1
        ELSE
          this%hybridCoeff(:, 1) = coeff
          istatus = yaml_get_var(TRIM(staticDir)//'/'//TRIM(this%hybridFile), 'hybrid_levels', 'b', coeff)
          this%hybridCoeff(:, 2) = coeff(1:this%vLevel)
        END IF

        ! Not reading in hybrid coefficients, set default values:
        IF (istatus .NE. 0) THEN
          WRITE (*, 11) istatus
11        FORMAT('SingleGrid: failed reading in hybridCoeff: ', I4)
          ! Default hybrid coefficients:
          this%hybridCoeff = 0.0D0
          this%hybridCoeff(:, 1) = 1.0D0
          this%hybridCoeff(:, 1) = 0.0D0
        ELSE
          WRITE (*, 12) this%hybridCoeff(1,:),this%hybridCoeff(this%vLevel,:)
12        FORMAT('SingleGrid: hybridCoeff is set to value at 1: ', 2E12.4,' at vlevel: ', 2E12.4)
        END IF
      ELSE
        this%ptop = -1.0D0 !< Default ptop is -1.0D0, indicating no hybrid vertical coordinate
        WRITE (*, 8) istatus
8       FORMAT('SingleGrid: failed reading in ptop: ', I4)
      END IF
    END IF

    ! Yuanfu Xie added a parameter of vertcalBreak ingest from yaml for vertcalSimilarity setting:
    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'appliedSimilarity', this%appliedSimilarity)
    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'scales4Similarity', this%scales4Similarity)
    IF (istatus .NE. 0) THEN
      WRITE (*, 2) istatus
2       FORMAT('SingleGrid: failed reading in scales4Similarity: ', I4)
      IF (ALLOCATED(this%scales4Similarity)) DEALLOCATE (this%scales4Similarity)
      ALLOCATE (this%scales4Similarity(2))
      this%scales4Similarity = 0.0D0
      this%appliedSimilarity = .FALSE.
    ELSE
      WRITE (*, 3) this%scales4Similarity
3       FORMAT('SingleGrid: read in scales4Similarity: ', 2E12.4)
    END IF

    ! If appliedSimilarity is .FALSE., we need to set scales4Similarity to 0
    IF (.NOT. this%appliedSimilarity) this%scales4Similarity = 0.0D0

    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'vertcalBreak', this%vertcalBreak)
    istatus = yaml_get_var(TRIM(configFile), 'geometry', 'flowDependentThinning', this%flowDependentThinning)

    IF (yaml_get_var(TRIM(configFile), 'geometry', 'haloWidth', this%halo_width) .NE. 0) this%halo_width = 1
      
    WRITE (*, 4) this%halo_width
4   FORMAT('SingleGrid: halo_width is set to value ', I4)

    this%num_cell_stcl = this%mpddInfo_sg%cal_num_cell_stcl(this%halo_width)


    WRITE (*, 1) gLevel, this%mpddInfo_sg%myrank, this%flowDependentThinning, this%appliedSimilarity, this%scales4Similarity
1   FORMAT('End construction of singleGrid on G', 2I3, ' FlowDpdtThinning: ', L, ' Similarity: ', L, 8E14.6)

  END SUBROUTINE initialize

  SUBROUTINE aggrGridReal1D(this, bufSrc, bufDest, bufDestShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: bufDestShape
    REAL(r_kind), INTENT(in) :: bufSrc(:)
    REAL(r_kind), INTENT(out) :: bufDest(:)

    CALL this%mpddInfo_sg%aggrGridReal1D(this%num_icell, bufSrc, bufDest, bufDestShape, this%seq_store_idx)
  END SUBROUTINE

  SUBROUTINE aggrGridReal2D(this, bufSrc, bufDest, bufDestShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: bufDestShape(2)
    REAL(r_kind), INTENT(in) :: bufSrc(:, :)
    REAL(r_kind), INTENT(out) :: bufDest(:, :)

    CALL this%mpddInfo_sg%aggrGridReal2D(this%num_icell, bufSrc, bufDest, bufDestShape, this%seq_store_idx)
  END SUBROUTINE

  ! SUBROUTINE aggrGridReal3D(this, bufSrc, bufDest, bufDestShape)
  !   IMPLICIT NONE
  !   CLASS(SingleGrid_t) this
  !   INTEGER(i_kind), INTENT(in) :: bufDestShape(3)
  !   REAL(r_kind), INTENT(in) :: bufSrc(:, :, :)
  !   REAL(r_kind), INTENT(out) :: bufDest(:, :, :)

  !   CALL this%mpddInfo_sg%aggrGridReal3D(this%num_icell, bufSrc, bufDest, bufDestShape, this%seq_store_idx)
  ! END SUBROUTINE aggrGridReal3D

!> @brief
!! Allocate the grid mat which is compatiable with this geometry sg and mg
! @see
! @note
! @warning
! @attention
  SUBROUTINE allocateMat(this, vLevel, buf)
    CLASS(SingleGrid_t) this

    INTEGER(i_kind), INTENT(IN) :: vLevel    !< Vertical levels of grid
    REAL(r_kind), ALLOCATABLE :: buf(:, :)   !< Buffer address

    IF (ALLOCATED(buf)) DEALLOCATE (buf)
    ALLOCATE (buf(vLevel, this%num_cell))
  END SUBROUTINE allocateMat

  ! Allocate 3D matrix  added by Yu  Xin
  SUBROUTINE allocateMat3D(this, buf)
    CLASS(SingleGrid_t), INTENT(INOUT) :: this
    REAL(r_kind), ALLOCATABLE, INTENT(INOUT) :: buf(:, :, :)
    IF (ALLOCATED(buf)) DEALLOCATE (buf)
    ALLOCATE (buf(this%vLevel, this%num_cell, this%tSlots))
  END SUBROUTINE allocateMat3D

  SUBROUTINE allocateMatGlobalOnBase(this, vLevel, buf)
    CLASS(SingleGrid_t) this

    INTEGER(i_kind), INTENT(IN) :: vLevel
    REAL(r_kind), ALLOCATABLE :: buf(:, :)

    IF (.NOT. this%isBaseProc()) RETURN

    IF (ALLOCATED(buf)) DEALLOCATE (buf)
    ALLOCATE (buf(vLevel, this%num_iCell_global))
  END SUBROUTINE allocateMatGlobalOnBase

  ! Allocate 3D matrix on base processor added by Yu Xin
  SUBROUTINE allocateMatGlobalOnBase3D(this, buf)
    CLASS(SingleGrid_t), INTENT(INOUT) :: this
    REAL(r_kind), ALLOCATABLE, INTENT(INOUT) :: buf(:, :, :)

    IF (.NOT. this%isBaseProc()) RETURN

    IF (ALLOCATED(buf)) DEALLOCATE (buf)
    ALLOCATE (buf(this%vLevel, this%num_iCell_global, this%tSlots))
  END SUBROUTINE allocateMatGlobalOnBase3D

  LOGICAL FUNCTION isBaseProc(this)
    CLASS(SingleGrid_t) :: this
    isBaseProc = this%mpddInfo_sg%isBaseProc()
  END FUNCTION isBaseProc

!> @brief
!! Exchange and update the halos on the exBuf on each single proc on this single grid.
! @see
! @note
! @warning
! @attention
  SUBROUTINE ExchangeMatOnHalo2D(this, vlevel, exBuf)
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: vlevel                          !> Number of vertical level
    REAL(r_kind), INTENT(inout) :: exBuf(vlevel, this%num_cell)    !> Exhange buffer

    CALL this%mpddInfo_sg%ExchangeMatOnHalo(this%num_icell, vlevel, exBuf, this%g_t_sp_idx)

  END SUBROUTINE ExchangeMatOnHalo2D

  SUBROUTINE ExchangeMatOnHaloReverseSum2D(this, vlevel, exBuf)
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: vlevel                          !> Number of vertical level
    REAL(r_kind), INTENT(inout) :: exBuf(vlevel, this%num_cell)    !> Exhange buffer

    CALL this%mpddInfo_sg%ExchangeMatOnHaloReverseSum(this%num_icell, vlevel, exBuf, this%g_t_sp_idx)

  END SUBROUTINE ExchangeMatOnHaloReverseSum2D

  SUBROUTINE DistGridRealWithHaloEx1D(this, buf_src, buf_Dist, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) :: this
    INTEGER, INTENT(IN) :: bufShape(1)
    REAL(r_kind), ALLOCATABLE, INTENT(IN)    :: buf_src(:)
    REAL(r_kind), ALLOCATABLE, INTENT(INOUT) :: buf_Dist(:)
    INTEGER :: k1

    IF (.NOT. this%isActiveProc()) RETURN

    buf_Dist = 0
    CALL this%DistGridReal(buf_src, buf_Dist, bufShape)
    CALL this%mpddInfo_sg%ExchangeMatOnHalo(this%num_icell, 1, buf_Dist, this%g_t_sp_idx)
  END SUBROUTINE DistGridRealWithHaloEx1D

  SUBROUTINE DistGridRealWithHaloEx2D(this, buf_src, buf_Dist, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) :: this
    INTEGER, INTENT(IN) :: bufShape(2)
    REAL(r_kind), ALLOCATABLE, INTENT(IN)    :: buf_src(:, :)
    REAL(r_kind), INTENT(INOUT) :: buf_Dist(:, :)
    INTEGER :: k1

    IF (.NOT. this%isActiveProc()) RETURN

    CALL this%DistGridReal(buf_src, buf_Dist, bufShape)
    CALL this%mpddInfo_sg%ExchangeMatOnHalo(this%num_icell, bufShape(1), buf_Dist, this%g_t_sp_idx)
  END SUBROUTINE DistGridRealWithHaloEx2D

  SUBROUTINE DistGridRealWithHaloEx3D(this, buf_src, buf_Dist, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) :: this
    INTEGER, INTENT(IN) :: bufShape(3)
    REAL(r_kind), ALLOCATABLE, INTENT(IN)    :: buf_src(:, :, :)
    REAL(r_kind), ALLOCATABLE, INTENT(INOUT) :: buf_Dist(:, :, :)
    INTEGER :: i, j

    IF (.NOT. this%isActiveProc()) RETURN

    CALL this%DistGridReal(buf_src, buf_Dist, bufShape)
    DO i = 1, bufShape(1)
      CALL this%mpddInfo_sg%ExchangeMatOnHalo(this%num_icell, bufShape(2), buf_Dist(i, :, :), this%g_t_sp_idx)
    END DO
  END SUBROUTINE DistGridRealWithHaloEx3D

  SUBROUTINE DistGridRealWithHaloEx4D(this, buf_src, buf_Dist, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) :: this
    INTEGER, INTENT(IN) :: bufShape(4)
    REAL(r_kind), ALLOCATABLE, INTENT(IN)    :: buf_src(:, :, :, :)
    REAL(r_kind), ALLOCATABLE, INTENT(INOUT) :: buf_Dist(:, :, :, :)
    INTEGER :: i, j, k

    IF (.NOT. this%isActiveProc()) RETURN

    CALL this%DistGridReal(buf_src, buf_Dist, bufShape)
    DO i = 1, bufShape(1)
      DO j = 1, bufShape(2)
        CALL this%mpddInfo_sg%ExchangeMatOnHalo(this%num_icell, bufShape(3), buf_Dist(i, j, :, :), this%g_t_sp_idx)
      END DO
    END DO
  END SUBROUTINE DistGridRealWithHaloEx4D

  LOGICAL FUNCTION isActiveProc(this)
    CLASS(SingleGrid_t) :: this
    isActiveProc = this%mpddInfo_sg%isActiveProc()
  END FUNCTION isActiveProc

  SUBROUTINE DistGridInt4D(this, bufSrc, bufDest, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: bufSrc(:, :, :, :), bufShape(4)
    INTEGER(i_kind), INTENT(out) :: bufDest(:, :, :, :)

    CALL this%mpddInfo_sg%DistGridIntSeq(this%num_icell, bufSrc, bufDest, bufShape, this%seq_store_idx)
  END SUBROUTINE DistGridInt4D

  SUBROUTINE DistGridInt3D(this, bufSrc, bufDest, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: bufSrc(:, :, :), bufShape(3)
    INTEGER(i_kind), INTENT(out) :: bufDest(:, :, :)

    CALL this%mpddInfo_sg%DistGridIntSeq(this%num_icell, bufSrc, bufDest, bufShape, this%seq_store_idx)
  END SUBROUTINE DistGridInt3D

  SUBROUTINE DistGridInt2D(this, bufSrc, bufDest, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: bufSrc(:, :), bufShape(2)
    INTEGER(i_kind), INTENT(out) :: bufDest(:, :)

    CALL this%mpddInfo_sg%DistGridIntSeq(this%num_icell, bufSrc, bufDest, bufShape, this%seq_store_idx)
  END SUBROUTINE DistGridInt2D

!> @brief
!! Build the restriction matrix for the aggregation 1
! @see
! @note
! @warning
! @attention
  SUBROUTINE DistGridInt1D(this, bufSrc, bufDest, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: bufSrc(:), &                ! Buffer source
                                   bufShape(1)                 ! Shape of buffer source
    INTEGER(i_kind), INTENT(out) :: bufDest(:)                 ! Buffer destination

    CALL this%mpddInfo_sg%DistGridIntSeq(this%num_icell, bufSrc, bufDest, bufShape, this%seq_store_idx)
  END SUBROUTINE DistGridInt1D

  SUBROUTINE DistGridReal4D(this, bufSrc, bufDest, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: bufShape(4)
    REAL(r_kind), INTENT(in) :: bufSrc(:, :, :, :)
    REAL(r_kind), INTENT(out) :: bufDest(:, :, :, :)

    CALL this%mpddInfo_sg%DistGridRealSeq(this%num_icell, bufSrc, bufDest, bufShape, this%seq_store_idx)
  END SUBROUTINE DistGridReal4D

  SUBROUTINE DistGridReal3D(this, bufSrc, bufDest, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: bufShape(3)
    REAL(r_kind), INTENT(in) :: bufSrc(:, :, :)
    REAL(r_kind), INTENT(out) :: bufDest(:, :, :)

    CALL this%mpddInfo_sg%DistGridRealSeq(this%num_icell, bufSrc, bufDest, bufShape, this%seq_store_idx)
  END SUBROUTINE DistGridReal3D

  SUBROUTINE DistGridReal2D(this, bufSrc, bufDest, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: bufShape(2)
    REAL(r_kind), INTENT(in) :: bufSrc(:, :)
    REAL(r_kind), INTENT(out) :: bufDest(:, :)

    CALL this%mpddInfo_sg%DistGridRealSeq(this%num_icell, bufSrc, bufDest, bufShape, this%seq_store_idx)
  END SUBROUTINE DistGridReal2D

  SUBROUTINE DistGridReal1D(this, bufSrc, bufDest, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: bufShape(1)
    REAL(r_kind), INTENT(in) :: bufSrc(:)
    REAL(r_kind), INTENT(out) :: bufDest(:)

    CALL this%mpddInfo_sg%DistGridRealSeq(this%num_icell, bufSrc, bufDest, bufShape, this%seq_store_idx)
  END SUBROUTINE DistGridReal1D

  SUBROUTINE DistGridRealForFieldGrid(this, bufSrc, bufDest, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: bufShape(3)
    REAL(r_kind), INTENT(in) :: bufSrc(:, :, :)
    REAL(r_kind), INTENT(out) :: bufDest(:, :, :)
    INTEGER(r_kind) :: i

    IF (.NOT. this%isActiveProc()) RETURN

    DO i = 1, this%tSlots
      CALL this%mpddInfo_sg%DistGridRealSeq(this%num_icell, bufSrc(:, :, i), bufDest(:, :, i), bufShape(1:2), this%seq_store_idx)
    END DO
  END SUBROUTINE DistGridRealForFieldGrid

  SUBROUTINE DistGridRealWithHaloExForFieldGrid(this, buf_src, buf_Dist, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) :: this
    INTEGER, INTENT(IN) :: bufShape(3)
    ! real(r_kind), ALLOCATABLE, INTENT(IN)    :: buf_src(:, :, :)
    ! real(r_kind), ALLOCATABLE, INTENT(INOUT) :: buf_Dist(:, :, :)
    REAL(r_kind), INTENT(IN)    :: buf_src(:, :, :)
    REAL(r_kind), INTENT(INOUT) :: buf_Dist(:, :, :)
    INTEGER :: i, j

    IF (.NOT. this%isActiveProc()) RETURN

    CALL this%DistGridRealForFieldGrid(buf_src, buf_Dist, bufShape)
    DO i = 1, bufShape(3)
      CALL this%mpddInfo_sg%ExchangeMatOnHalo(this%num_icell, bufShape(1), buf_Dist(:, :, i), this%g_t_sp_idx)
    END DO

  END SUBROUTINE DistGridRealWithHaloExForFieldGrid

  SUBROUTINE DistGridRealWithHaloExForFieldGrid4D(this, buf_src, buf_Dist, bufShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) :: this
    INTEGER, INTENT(IN) :: bufShape(4)
    REAL(r_kind), ALLOCATABLE, INTENT(IN)    :: buf_src(:, :, :, :)
    REAL(r_kind), ALLOCATABLE, INTENT(INOUT) :: buf_Dist(:, :, :, :)
    INTEGER :: i, j

    IF (.NOT. this%isActiveProc()) RETURN

    DO j = 1, bufShape(4)
      CALL this%DistGridRealForFieldGrid(buf_src(:, :, :, j), buf_Dist(:, :, :, j), bufShape(1:3))

      DO i = 1, bufShape(3)
        CALL this%mpddInfo_sg%ExchangeMatOnHalo(this%num_icell, bufShape(1), buf_Dist(:, :, i, j), this%g_t_sp_idx)
      END DO
    END DO

  END SUBROUTINE DistGridRealWithHaloExForFieldGrid4D

  SUBROUTINE aggrGridRealForFieldGrid(this, bufSrc, bufDest, bufDestShape)
    IMPLICIT NONE
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: bufDestShape(3)
    REAL(r_kind), INTENT(in) :: bufSrc(:, :, :)
    REAL(r_kind), INTENT(out) :: bufDest(:, :, :)
    INTEGER :: i
    DO i = 1, bufDestShape(3)
      ! CALL this%mpddInfo_sg%aggrGridReal3D(this%num_icell, bufSrc, bufDest, bufDestShape, this%seq_store_idx)
      CALL this%mpddInfo_sg%aggrGridReal2D(this%num_icell, bufSrc(:, :, i), bufDest(:, :, i), bufDestShape(1:2), this%seq_store_idx)
    END DO
  END SUBROUTINE aggrGridRealForFieldGrid

  SUBROUTINE ExchangeMatOnHaloForFieldGrid(this, tSlots, vlevel, exBuf)
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: vlevel, tSlots                               !> Number of vertical level
    REAL(r_kind), INTENT(inout) :: exBuf(vlevel, this%num_cell, tSlots)    !> Exhange buffer
    INTEGER(r_kind) :: i

    IF (.NOT. this%isActiveProc()) RETURN

    DO i = 1, tSlots
      CALL this%mpddInfo_sg%ExchangeMatOnHalo(this%num_icell, vlevel, exBuf(:, :, i), this%g_t_sp_idx)
    END DO
  END SUBROUTINE ExchangeMatOnHaloForFieldGrid

  SUBROUTINE ExchangeMatOnHaloReverseSumForFieldGrid(this, tSlots, vlevel, exBuf)
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(in) :: vlevel, tSlots                               !> Number of vertical level
    REAL(r_kind), INTENT(inout) :: exBuf(vlevel, this%num_cell, tSlots)    !> Exhange buffer
    INTEGER(r_kind) :: i

    IF (.NOT. this%isActiveProc()) RETURN

    DO i = 1, tSlots
      CALL this%mpddInfo_sg%ExchangeMatOnHaloReverseSum(this%num_icell, vlevel, exBuf(:, :, i), this%g_t_sp_idx)
    END DO
  END SUBROUTINE

  ! Yuanfu Xie added boundaryItp: using closest neighbor cell's weighted average to extrapolate:
  SUBROUTINE boundaryOI(this, vlvl, tslot, bndy)
    USE parameters_m, ONLY: EarthRadius
    USE geoTools_m, ONLY: distance
    IMPLICIT NONE
    CLASS(SingleGrid_t) this
    INTEGER(i_kind), INTENT(IN) :: vlvl, tslot
    REAL(r_kind), INTENT(INOUT) :: bndy(vlvl, this%num_cell, tslot)

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, neighbor
    LOGICAL :: itpl
    REAL(r_kind) :: closest, total, sigma
    REAL(r_kind), ALLOCATABLE :: average(:, :), dist(:)

    ! Allocate average array with vlvl and tslot:
    ALLOCATE (average(vlvl, tslot), dist(UBOUND(this%cell_stcl, 1)))

    DO j = 1, this%num_icell
      ! Not interior or halo cell:
      IF (this%cell_type(j) .NE. 0) THEN
        i = 0
        closest = EarthRadius
        ! Find the closest interior cell
        DO k = 1, UBOUND(this%cell_stcl, 1)
          IF (this%cell_stcl(k, j) .GT. 0 .AND. &
              this%cell_stcl(k, j) .NE. j .AND. & ! Not j cell itself
              this%cell_dist(k, j) .LE. closest) THEN

            IF (this%cell_type(this%cell_stcl(k, j)) .EQ. 0 .OR. &
                this%cell_type(this%cell_stcl(k, j)) .EQ. 3) THEN
              IF (this%cell_dist(k, j) .LT. closest) THEN ! Use the closest interior or point only
                i = k
                closest = this%cell_dist(k, j)
              END IF
            END IF
          END IF
        END DO
        ! Change the stcl index to grid index:
        i = this%cell_stcl(i, j)

        ! Total distance weight at i:
        total = 0.0D0 ! count in the extrapolated point's weight
        average = 0.0D0
        neighbor = 0
        sigma = 0.0D0
        dist = 0.0D0
        DO k = 1, UBOUND(this%cell_stcl, 1)
          IF (this%cell_stcl(k, i) .GT. 0) THEN
            IF (this%cell_type(this%cell_stcl(k, i)) .NE. 1 .AND. &
                this%cell_type(this%cell_stcl(k, i)) .NE. 2) THEN ! Not use boundary points
              average = average + bndy(:, this%cell_stcl(k, i), :)
              neighbor = neighbor + 1
              dist(k) = distance(this%cell_cntr(:, this%cell_stcl(k, i)), this%cell_cntr(:, j))
              IF (sigma .LT. dist(k)) sigma = dist(k)
            END IF
          END IF
        END DO
        average = average / DBLE(neighbor)

        !sigma = sigma*2.0D0

        DO k = 1, UBOUND(this%cell_stcl, 1)
          IF (this%cell_stcl(k, i) .GT. 0) THEN
            IF (this%cell_type(this%cell_stcl(k, i)) .NE. 1 .AND. &
                this%cell_type(this%cell_stcl(k, i)) .NE. 2) THEN ! Not use boundary points
              total = total + EXP(-0.5D0 * (dist(k) / sigma)**2)
            END IF
          END IF
        END DO

        ! Interpolation/extrapolation:
        bndy(:, j, :) = 0.0D0
        DO k = 1, UBOUND(this%cell_stcl, 1)
          IF (this%cell_stcl(k, i) .GT. 0) THEN
            IF (this%cell_type(this%cell_stcl(k, i)) .NE. 1 .AND. &
                this%cell_type(this%cell_stcl(k, i)) .NE. 2) THEN ! Not use boundary or corner points
              bndy(:, j, :) = bndy(:, j, :) + &
                              (bndy(:, this%cell_stcl(k, i), :) - average) * EXP(-0.5D0 * (dist(k) / sigma)**2)
            END IF
          END IF
        END DO
        bndy(:, j, :) = bndy(:, j, :) / total + average
      END IF
    END DO
    !print*,'Final OI: ',bndy(1,jdebug,1)

    ! Deallocate memory:
    DEALLOCATE (average, dist)
  END SUBROUTINE boundaryOI

  ! Yuanfu Xie added boundaryHalo: update halo cell_type = 4 if it sits on domain bounderies
  SUBROUTINE boundaryHalo(this)
    IMPLICIT NONE
    CLASS(SingleGrid_t) this

    ! Local variables:
    INTEGER(i_kind) :: i, j
    REAL(r_kind), ALLOCATABLE :: TYPE(:)

    ALLOCATE (TYPE(this%num_cell))

    TYPE = this%cell_type

    CALL this%ExchangeMatOnHalo2D(1, TYPE)

    DO i = this%num_icell + 1, this%num_cell
      IF (TYPE(i) .GT. 0.0D0) this%cell_type(i) = 4
    END DO

    DEALLOCATE (TYPE)
  END SUBROUTINE boundaryHalo

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(SingleGrid_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%edge_stcl)) DEALLOCATE (this%edge_stcl)
    IF (ALLOCATED(this%cell_area)) DEALLOCATE (this%cell_area)
    IF (ALLOCATED(this%cell_cntr)) DEALLOCATE (this%cell_cntr)
    IF (ALLOCATED(this%edge_leng)) DEALLOCATE (this%edge_leng)
    IF (ALLOCATED(this%cell_adnb)) DEALLOCATE (this%cell_adnb)  ! Yuanfu Xie added on 09/04/2022
    IF (ALLOCATED(this%edgeNorm2)) DEALLOCATE (this%edgeNorm2)  ! Yuanfu Xie added on 08/24/2022
    IF (ALLOCATED(this%edgeTang2)) DEALLOCATE (this%edgeTang2)  ! Yuanfu Xie added on 08/24/2022
    IF (ALLOCATED(this%coef_func)) DEALLOCATE (this%coef_func)
    IF (ALLOCATED(this%coef_norm)) DEALLOCATE (this%coef_norm)
    IF (ALLOCATED(this%coef_tang)) DEALLOCATE (this%coef_tang)
    IF (ALLOCATED(this%coef_vrtx)) DEALLOCATE (this%coef_vrtx)
    IF (ALLOCATED(this%coef_gl)) DEALLOCATE (this%coef_gl)

    IF (ALLOCATED(this%p_tc_idx)) DEALLOCATE (this%p_tc_idx)
    IF (ALLOCATED(this%c_tc_idx)) DEALLOCATE (this%c_tc_idx)

    ! if (ALLOCATED(this%c_tp_idx)) DEALLOCATE (this%c_tp_idx)
    IF (ALLOCATED(this%c2f_intp_stcl)) DEALLOCATE (this%c2f_intp_stcl)
    IF (ALLOCATED(this%c2f_intp_cf)) DEALLOCATE (this%c2f_intp_cf)
    IF (ALLOCATED(this%cell_type)) DEALLOCATE (this%cell_type)

    IF (ALLOCATED(this%sp_t_g_idx)) DEALLOCATE (this%sp_t_g_idx)
    IF (ALLOCATED(this%g_t_sp_idx)) DEALLOCATE (this%g_t_sp_idx)
    IF (ALLOCATED(this%seq_store_idx)) DEALLOCATE (this%seq_store_idx)

    IF (ALLOCATED(this%sp_t_g_idx_toCoarser)) DEALLOCATE (this%sp_t_g_idx_toCoarser)
    IF (ALLOCATED(this%g_t_sp_idx_toCoarser)) DEALLOCATE (this%g_t_sp_idx_toCoarser)
    IF (ALLOCATED(this%seq_store_idx_toCoarser)) DEALLOCATE (this%seq_store_idx_toCoarser)

    IF (ALLOCATED(this%sp_t_g_idx_toFiner)) DEALLOCATE (this%sp_t_g_idx_toFiner)
    IF (ALLOCATED(this%g_t_sp_idx_toFiner)) DEALLOCATE (this%g_t_sp_idx_toFiner)
    IF (ALLOCATED(this%seq_store_idx_toFiner)) DEALLOCATE (this%seq_store_idx_toFiner)
    IF (ALLOCATED(this%cell_stcl)) DEALLOCATE (this%cell_stcl)
    IF (ALLOCATED(this%nonZeroCellStcl)) DEALLOCATE (this%nonZeroCellStcl)

    IF (ALLOCATED(this%mpddInfo_sg)) DEALLOCATE (this%mpddInfo_sg)
    IF (ALLOCATED(this%cell_dist)) DEALLOCATE (this%cell_dist)

    IF (ALLOCATED(this%zHght)) DEALLOCATE (this%zHght)
    IF (ALLOCATED(this%aboveTopo)) DEALLOCATE (this%aboveTopo)
    IF (ALLOCATED(this%pres)) DEALLOCATE (this%pres)

    IF (ALLOCATED(this%sigma)) DEALLOCATE (this%sigma)
    IF (ALLOCATED(this%coef_sigma)) DEALLOCATE (this%coef_sigma)
    IF (ALLOCATED(this%topo)) DEALLOCATE (this%topo)
    IF (ALLOCATED(this%sigmax)) DEALLOCATE (this%sigmax)
    IF (ALLOCATED(this%sigmay)) DEALLOCATE (this%sigmay)
    IF (ALLOCATED(this%tt)) DEALLOCATE (this%tt)

    IF (ALLOCATED(this%s1)) DEALLOCATE (this%s1)
    IF (ALLOCATED(this%SVapor)) DEALLOCATE (this%SVapor)
    IF (ALLOCATED(this%STemp)) DEALLOCATE (this%STemp)
    IF (ALLOCATED(this%SVapor1D)) DEALLOCATE (this%SVapor1D)
    IF (ALLOCATED(this%STemp1D)) DEALLOCATE (this%STemp1D)
    IF (ALLOCATED(this%FGPres)) DEALLOCATE (this%FGPres)
    IF (ALLOCATED(this%s_qvapor)) DEALLOCATE (this%s_qvapor)

    IF (ALLOCATED(this%tskin)) DEALLOCATE (this%tskin)
    IF (ALLOCATED(this%psfc)) DEALLOCATE (this%psfc)
    IF (ALLOCATED(this%u10m)) DEALLOCATE (this%u10m)
    IF (ALLOCATED(this%v10m)) DEALLOCATE (this%v10m)
    IF (ALLOCATED(this%landmask)) DEALLOCATE (this%landmask)
    IF (ALLOCATED(this%soil_type)) DEALLOCATE (this%soil_type)
    IF (ALLOCATED(this%veg_frac)) DEALLOCATE (this%veg_frac)
    IF (ALLOCATED(this%snowc)) DEALLOCATE (this%snowc)
    IF (ALLOCATED(this%xice)) DEALLOCATE (this%xice)
    IF (ALLOCATED(this%scales4Similarity)) DEALLOCATE (this%scales4Similarity)
    IF (ALLOCATED(this%vertcalSimilarity)) DEALLOCATE (this%vertcalSimilarity)
    IF (ALLOCATED(this%horizonSimilarity)) DEALLOCATE (this%horizonSimilarity)

    IF (ALLOCATED(this%lat1DAtBase)) DEALLOCATE (this%lat1DAtBase)
    IF (ALLOCATED(this%lon1DAtBase)) DEALLOCATE (this%lon1DAtBase)

    IF (ALLOCATED(this%BufferHalo)) DEALLOCATE (this%BufferHalo)
    IF (ALLOCATED(this%F_invHz_z)) DEALLOCATE (this%F_invHz_z)
    IF (ALLOCATED(this%F_Hz_z)) DEALLOCATE (this%F_Hz_z)
    IF (ALLOCATED(this%F_z_z)) DEALLOCATE (this%F_z_z)
    IF (ALLOCATED(this%Lz)) DEALLOCATE (this%Lz)
    IF (ALLOCATED(this%parHz_parsigma)) DEALLOCATE (this%parHz_parsigma)
    IF (ALLOCATED(this%Hz)) DEALLOCATE (this%Hz)
    IF (ALLOCATED(this%parz_parsigma)) DEALLOCATE (this%parz_parsigma)
    IF (ALLOCATED(this%coef_scddif)) DEALLOCATE (this%coef_scddif)
    IF (ALLOCATED(this%coef_fstdif)) DEALLOCATE (this%coef_fstdif)
    IF (ALLOCATED(this%InterpCoef)) DEALLOCATE (this%InterpCoef)
    IF (ALLOCATED(this%coef_fstdif_half)) DEALLOCATE (this%coef_fstdif_half)
    IF (ALLOCATED(this%f)) DEALLOCATE (this%f)
    IF (ALLOCATED(this%BdyBeta)) DEALLOCATE (this%BdyBeta)
    IF (ALLOCATED(this%BdyAlpha)) DEALLOCATE (this%BdyAlpha)
    IF (ALLOCATED(this%bdy_type)) DEALLOCATE (this%bdy_type)

  END SUBROUTINE destructor

END MODULE SingleGrid_m
