!!--------------------------------------------------------------------------------------------------
! PROJECT           : Grid generation
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : Beta 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2020/05/29, @GBA-MWF, Shenzhen
!  Modified by Zilong Qin (zilong.qin@gmail.com), 2020/11/23, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!!===================================================================
!> @brief
!! # Grid Structure Module
!!
!!  *This module defines data structures for grid specification*
!!
!! ## Grid Parameters:
!!  *It defines all necessary parameters for a given grid. It depends
!!   on a macro MID_EDGE, which specifies if a finite volume uses a
!!   middle point of an edge to discretize the normal and tangential
!!   gradients*
!!
!! ## Grid Geometry Quantities:
!!  * It provides all needed arrays of grid geometry, such as cell,
!!    edge, prolongation coefficients etc.
!!
!! @author Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
!!===================================================================
MODULE gridStruct_m
  USE kinds_m, ONLY: i_kind, r_kind

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: gridParams_t, gridGeoQty_t, multiLevel_t, &
            lengthFile, file_params, file_geoqty

  CHARACTER(LEN=256) :: path2static
  INTEGER(i_kind), PARAMETER :: lengthFile = 18
  CHARACTER(LEN=lengthFile), PARAMETER :: file_params = 'gridParameters.dat'
  CHARACTER(LEN=lengthFile), PARAMETER :: file_geoqty = 'gridGeoquantty.dat'

  TYPE :: gridParams_t
    INTEGER(i_kind) :: num_vrtx, num_cell, num_chld
    INTEGER(i_kind) :: numStclPerEdge, numQuadPerEdge
    INTEGER(i_kind) :: numStclPerCell, numStclPerVrtx, numVrtxPerCell
  END TYPE gridParams_t

  TYPE :: gridGeoQty_t
    ! INTEGER(i_kind) :: max_c_stcl
    ! This following equivalency is referred to the ITmeshGeoQnty by Ning Wang's icosahedron
    INTEGER(i_kind), ALLOCATABLE :: cell_stcl(:, :)    ! equivalent to cv_stcl
    INTEGER(i_kind), ALLOCATABLE :: cell_vrtx(:, :)    ! equivalent to cv_vrtx
    INTEGER(i_kind), ALLOCATABLE :: cell_opsd(:, :)    ! equivalent to nb
    INTEGER(i_kind), ALLOCATABLE :: cell_adnb(:, :)    ! equivalent to adjnb
    INTEGER(i_kind), ALLOCATABLE :: cell_side(:, :)    ! equivalent to adjnb_si
    INTEGER(i_kind), ALLOCATABLE :: edge_stcl(:, :, :)  ! equivalent to edge_stcl
    INTEGER(i_kind), ALLOCATABLE :: edge_vrtx(:, :, :)  ! equivalent to edge_vrtx
    INTEGER(i_kind), ALLOCATABLE :: vrtx_stcl(:, :, :)  ! equivalent to vrtx_stcl
    INTEGER(i_kind), ALLOCATABLE :: prnt_indx(:)      ! equivalent to p_tc_idx
    INTEGER(i_kind), ALLOCATABLE :: chld_indx(:, :)    ! equivalent to c_tc_idx
    INTEGER(i_kind), ALLOCATABLE :: cTof_stcl(:, :)    ! equivalent to c2f_intp_stcl
    ! Note: the above cTof_stcl is from the view of the fine grid, the current level
    ! It states which indices from the grid above are used to interpolate the this
    ! cell. The corresponding coefficients are in the REAL declaration below

    ! Gauss-Legrendre quadrature weights:
    REAL(r_kind), ALLOCATABLE :: coef_gl(:)

    REAL(r_kind), PUBLIC, ALLOCATABLE :: vrtx_lalo(:, :)    ! equivalent to icos_grid; dim: 2 x num_vrtx
    REAL(r_kind), ALLOCATABLE :: cell_area(:)      ! equivalent to area num_cell
    REAL(r_kind), ALLOCATABLE :: cell_cntr(:, :)    ! equivalent to center (lat/lon) 2 x num_cell
    REAL(r_kind), ALLOCATABLE :: edge_lnth(:, :)    ! equivalent to lnth numVrtxPerCell x num_cell
    ! Consider to save them in num_vrtx
    REAL(r_kind), ALLOCATABLE :: tgnt_vctr(:, :, :)  ! equivalent to t_vec
    REAL(r_kind), ALLOCATABLE :: norm_vctr(:, :, :)  ! equivalent to n_vec
    REAL(r_kind), ALLOCATABLE :: tgnt_vct2(:, :, :)  ! equivalent to t_vec2d
    REAL(r_kind), ALLOCATABLE :: norm_vct2(:, :, :)  ! equivalent to n_vec2d
    REAL(r_kind), ALLOCATABLE :: coef_norm(:, :, :, :)  ! equivalent to itpl_cf_nd
    REAL(r_kind), ALLOCATABLE :: coef_tgnt(:, :, :, :)  ! equivalent to itpl_cf_td
    REAL(r_kind), ALLOCATABLE :: coef_func(:, :, :, :)  ! equivalent to itpl_cf_fv
    REAL(r_kind), ALLOCATABLE :: coef_vrtx(:, :, :)  ! equivalent to itpl_v_cf_fv in V2 gridGen
    REAL(r_kind), ALLOCATABLE :: cTof_coef(:, :)    ! equivalent to c2f_intp_cf

  CONTAINS
    PROCEDURE, PUBLIC :: initial
    PROCEDURE, PUBLIC :: destroy
  END TYPE gridGeoQty_t

  TYPE :: multiLevel_t
    INTEGER(i_kind) :: num_lvls, lvl_stts, lvl_ends
    TYPE(gridParams_t), ALLOCATABLE :: Params(:)
    TYPE(gridGeoQty_t), ALLOCATABLE :: GeoQty(:)

  CONTAINS
    PROCEDURE, PUBLIC :: allc_arrays  ! Allocatae the arrays
    PROCEDURE, PUBLIC :: delt_arrays  ! Deallocate the arrays
    PROCEDURE, PUBLIC :: allc_geoqty  ! Allocatae the geoqty
    PROCEDURE, PUBLIC :: delt_geoqty  ! Deallocate the geoqty
    PROCEDURE, PUBLIC :: read_mlgrid
    PROCEDURE, PUBLIC :: writ_mlgrid

  END TYPE multiLevel_t

CONTAINS

  SUBROUTINE allc_arrays(this, mglvls, mgstts, mgends)
    CLASS(multiLevel_t) :: this
    INTEGER(i_kind), INTENT(IN) :: mglvls, mgstts, mgends

    ! Local variables:
    INTEGER(i_kind) :: m

    this%num_lvls = mglvls
    this%lvl_stts = mgstts
    this%lvl_ends = mgends

    ALLOCATE (this%Params(this%lvl_stts:this%lvl_ends), &
              this%Geoqty(this%lvl_stts:this%lvl_ends))

  END SUBROUTINE allc_arrays

  SUBROUTINE delt_arrays(this)
    CLASS(multiLevel_t) :: this

    IF (ALLOCATED(this%Params)) DEALLOCATE (this%Params)
    IF (ALLOCATED(this%Geoqty)) DEALLOCATE (this%Geoqty)

  END SUBROUTINE delt_arrays

  SUBROUTINE allc_geoqty(this)
    CLASS(multiLevel_t) :: this

    ! Local variables:
    INTEGER(i_kind) :: m

    ! Allocate memory for all elements:
    DO m = this%lvl_stts, this%lvl_ends
      CALL this%GeoQty(m)%initial(this%Params(m))
    END DO
  END SUBROUTINE allc_geoqty

  SUBROUTINE delt_geoqty(this)
    CLASS(multiLevel_t) :: this

    ! Local variables:
    INTEGER(i_kind) :: m

    ! Deallocate memory for all elements:
    IF (ALLOCATED(this%Geoqty)) THEN
      DO m = LBOUND(this%GeoQty, 1), UBOUND(this%GeoQty, 1)
        CALL this%GeoQty(m)%destroy
      END DO
    END IF
  END SUBROUTINE delt_geoqty

  SUBROUTINE read_mlgrid(this, grdOnly)
    CLASS(multiLevel_t) :: this
    LOGICAL, INTENT(IN) :: grdOnly

    ! Local variables:
    INTEGER(i_kind) :: m

    ! Open multigrid files to read:
    CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", path2static)
    OPEN (10, FILE=TRIM(path2static)//"/"//TRIM(file_params(1:lengthFile)), FORM='UNFORMATTED')
    READ (10) this%num_lvls, this%lvl_stts, this%lvl_ends
    ALLOCATE (this%Params(this%lvl_stts:this%lvl_ends), &
              this%GeoQty(this%lvl_stts:this%lvl_ends))
    DO m = this%lvl_stts, this%lvl_ends
      READ (10) this%Params(m)
    END DO
    CLOSE (10)
    DO m = this%lvl_stts, this%lvl_ends
      CALL this%GeoQty(m)%initial(this%Params(m))
    END DO
    CLOSE (10)

    OPEN (10, FILE=TRIM(path2static)//"/"//TRIM(file_geoqty(1:lengthFile)), FORM='UNFORMATTED')
    IF (grdOnly) THEN
      DO m = this%lvl_stts, this%lvl_ends
        READ (10) this%GeoQty(m)%vrtx_lalo(1, 1:this%Params(m)%num_vrtx), &
          this%GeoQty(m)%vrtx_lalo(2, 1:this%Params(m)%num_vrtx), &
          this%GeoQty(m)%cell_vrtx(1:this%Params(m)%numVrtxPerCell, &
                                   1:this%Params(m)%num_cell), &
          this%GeoQty(m)%cell_adnb(1:this%Params(m)%numVrtxPerCell, &
                                   1:this%Params(m)%num_cell)
      END DO
    ELSE
      DO m = this%lvl_stts, this%lvl_ends
        !print*,'vrtx_lalo:',m,sizeof(this%GeoQty(m)%vrtx_lalo)
        READ (10) this%GeoQty(m)%vrtx_lalo
        !print*,'area/cntr: ',sizeof(this%GeoQty(m)%cell_area), &
        !sizeof(this%GeoQty(m)%cell_cntr)
        READ (10) this%GeoQty(m)%cell_area, this%GeoQty(m)%cell_cntr
        !print*,'cell_stcl/vrtx: ',sizeof(this%GeoQty(m)%cell_stcl), &
        !sizeof(this%GeoQty(m)%cell_vrtx)
        READ (10) this%GeoQty(m)%cell_stcl, this%GeoQty(m)%cell_vrtx
        !print*,'cell_adnb/side: ',sizeof(this%GeoQty(m)%cell_adnb), &
        !sizeof(this%GeoQty(m)%cell_side)
        READ (10) this%GeoQty(m)%cell_adnb, this%GeoQty(m)%cell_side
        !print*,'edge_lnth/opsd:',sizeof(this%GeoQty(m)%edge_lnth), &
        !sizeof(this%GeoQty(m)%cell_opsd)
        READ (10) this%GeoQty(m)%edge_lnth, this%GeoQty(m)%cell_opsd
        !print*,'edge_stcl/vrtx: ',sizeof(this%GeoQty(m)%edge_stcl), &
        !sizeof(this%GeoQty(m)%edge_vrtx)
        READ (10) this%GeoQty(m)%edge_stcl, this%GeoQty(m)%edge_vrtx
        !print*,'vrtx_stcl:',sizeof(this%GeoQty(m)%vrtx_stcl)
        READ (10) this%GeoQty(m)%vrtx_stcl
        READ (10) this%GeoQty(m)%tgnt_vctr, this%GeoQty(m)%norm_vctr
        READ (10) this%GeoQty(m)%tgnt_vct2, this%GeoQty(m)%norm_vct2
        READ (10) this%GeoQty(m)%coef_norm, this%GeoQty(m)%coef_tgnt
        READ (10) this%GeoQty(m)%coef_func
#ifndef MID_EDGE
        READ (10) this%GeoQty(m)%coef_vrtx
#endif
        READ (10) this%GeoQty(m)%prnt_indx
        READ (10) this%GeoQty(m)%chld_indx
        READ (10) this%GeoQty(m)%cTof_stcl
        READ (10) this%GeoQty(m)%cTof_coef
      END DO
    END IF
    CLOSE (10)

  END SUBROUTINE read_mlgrid

  SUBROUTINE writ_mlgrid(this, grdOnly)
    CLASS(multiLevel_t) :: this
    LOGICAL, INTENT(IN) :: grdOnly

    ! Local variables:
    INTEGER(i_kind) :: m, i

    ! Open multigrid files to read:
    CALL GET_ENVIRONMENT_VARIABLE("GRID_DIR", path2static)
    OPEN (10, FILE=TRIM(path2static)//"/"//TRIM(file_params(1:lengthFile)), FORM='UNFORMATTED')
    WRITE (10) this%num_lvls, this%lvl_stts, this%lvl_ends
    DO m = this%lvl_stts, this%lvl_ends
      WRITE (10) this%Params(m)
    END DO
    CLOSE (10)

    OPEN (10, FILE=TRIM(path2static)//"/"//TRIM(file_geoqty(1:lengthFile)), FORM='UNFORMATTED')
    IF (grdOnly) THEN
      DO m = this%lvl_stts, this%lvl_ends
        WRITE (10) this%GeoQty(m)%vrtx_lalo(1, 1:this%Params(m)%num_vrtx), &
          this%GeoQty(m)%vrtx_lalo(2, 1:this%Params(m)%num_vrtx), &
          this%GeoQty(m)%cell_vrtx(1:this%Params(m)%numVrtxPerCell, &
                                   1:this%Params(m)%num_cell), &
          this%GeoQty(m)%cell_adnb(1:this%Params(m)%numVrtxPerCell, &
                                   1:this%Params(m)%num_cell)
      END DO
    ELSE
      DO m = this%lvl_stts, this%lvl_ends
        WRITE (10) this%GeoQty(m)%vrtx_lalo
        WRITE (10) this%GeoQty(m)%cell_area, this%GeoQty(m)%cell_cntr
        WRITE (10) this%GeoQty(m)%cell_stcl, this%GeoQty(m)%cell_vrtx
        WRITE (10) this%GeoQty(m)%cell_adnb, this%GeoQty(m)%cell_side
        WRITE (10) this%GeoQty(m)%edge_lnth, this%GeoQty(m)%cell_opsd
        WRITE (10) this%GeoQty(m)%edge_stcl, this%GeoQty(m)%edge_vrtx
        WRITE (10) this%GeoQty(m)%vrtx_stcl
        WRITE (10) this%GeoQty(m)%tgnt_vctr, this%GeoQty(m)%norm_vctr
        WRITE (10) this%GeoQty(m)%tgnt_vct2, this%GeoQty(m)%norm_vct2
        WRITE (10) this%GeoQty(m)%coef_norm, this%GeoQty(m)%coef_tgnt
        WRITE (10) this%GeoQty(m)%coef_func
#ifndef MID_EDGE
        WRITE (10) this%GeoQty(m)%coef_vrtx
#endif
        WRITE (10) this%GeoQty(m)%prnt_indx
        WRITE (10) this%GeoQty(m)%chld_indx
        WRITE (10) this%GeoQty(m)%cTof_stcl
        WRITE (10) this%GeoQty(m)%cTof_coef
      END DO
    END IF
    CLOSE (10)

  END SUBROUTINE writ_mlgrid

  SUBROUTINE initial(this, para)
    CLASS(gridGeoQty_t) :: this
    TYPE(gridParams_t) :: para

    ALLOCATE (this%cell_stcl(para%numStclPerCell, para%num_cell), &
              this%cell_vrtx(para%numVrtxPerCell, para%num_cell), &
              this%cell_opsd(para%numVrtxPerCell, para%num_cell), &
              this%cell_adnb(para%numVrtxPerCell, para%num_cell), &
              this%cell_side(para%numVrtxPerCell, para%num_cell), &
              this%edge_stcl(para%numStclPerEdge, para%numVrtxPerCell, para%num_cell), &
              this%edge_vrtx(2, para%numVrtxPerCell, para%num_cell), &
              this%vrtx_stcl(para%numStclPerVrtx, para%numVrtxPerCell, para%num_cell), &
              this%prnt_indx(para%num_cell), &
              this%chld_indx(para%num_chld, para%num_cell), &
              this%cTof_stcl(para%num_chld, para%num_cell))

    ALLOCATE (this%vrtx_lalo(2, para%num_vrtx), this%cell_area(para%num_cell), &
              this%cell_cntr(2, para%num_cell), &
              this%edge_lnth(para%numVrtxPerCell, para%num_cell), &
              this%tgnt_vctr(3, para%numVrtxPerCell, para%num_cell), & ! 3D vctr
              this%norm_vctr(3, para%numVrtxPerCell, para%num_cell), & ! 3D vctr
              this%tgnt_vct2(2, para%numVrtxPerCell, para%num_cell), &
              this%norm_vct2(2, para%numVrtxPerCell, para%num_cell), &
              ! If not using mid_edge for a finite volume, there are 2 points on an
              ! edge selected for function values, and normal, tangential derivatives.
              this%coef_norm(para%numStclPerEdge, para%numQuadPerEdge, &
                             para%numVrtxPerCell, para%num_cell), &
              this%coef_tgnt(para%numStclPerEdge, para%numQuadPerEdge, &
                             para%numVrtxPerCell, para%num_cell), &
              this%coef_func(para%numStclPerEdge, para%numQuadPerEdge, &
                             para%numVrtxPerCell, para%num_cell), &
              this%coef_vrtx(12, para%numVrtxPerCell, para%num_cell), &
              this%cTof_coef(para%num_chld, para%num_cell), &
              this%coef_gl(para%numQuadPerEdge))
  END SUBROUTINE initial

  SUBROUTINE destroy(This)
    CLASS(gridGeoQty_t) :: this

    IF (ALLOCATED(this%cell_stcl)) THEN
      DEALLOCATE (this%cell_stcl, this%cell_vrtx, this%cell_opsd, this%cell_adnb, &
                  this%cell_side, this%edge_stcl, this%edge_vrtx, this%vrtx_stcl, &
                  this%prnt_indx, this%chld_indx, this%cTof_stcl)
      DEALLOCATE (this%vrtx_lalo, this%cell_area, this%cell_cntr, this%edge_lnth, &
                  this%tgnt_vctr, this%norm_vctr, this%tgnt_vct2, this%norm_vct2, &
                  this%coef_norm, this%coef_tgnt, this%coef_func, this%coef_vrtx, &
                  this%cTof_coef, this%coef_gl)
    END IF
  END SUBROUTINE destroy
END MODULE gridStruct_m
