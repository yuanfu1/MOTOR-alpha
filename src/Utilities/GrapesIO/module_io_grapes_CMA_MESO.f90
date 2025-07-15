!!--------------------------------------------------------------------------------------------------
! PROJECT           : GRAPES IO
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Sanshan Tu
! VERSION           : Beta 0.0
! HISTORY           :
!   Created by Sanshan Tu (tss71618@163.com), 2020/12/31, @SZSC, Shenzhen
!   Modified by Zhao Liu (liuzhao@nsccsz.cn), 2021/3/18, @SZSC, Shenzhen
!   Modified by Sanshan Tu (tss71618@163.com), 2021/11/01, @SZSC, Shenzhen
!  Modified by Zilong Qin (zilong.qin@gmail.com), 2022/5/13, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!!===================================================================
!> @brief
!! # GRAPES IO Module
!!
!!  *This module defines the subroutines for GRAPES IO*
!! @author Sanshan Tu
!! @copyright (C) 2023 GBA-MWF, All rights reserved.
!!===================================================================
MODULE module_io_grapes_CMA_MESO_m
  USE module_domain
  USE module_configure
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE module_io_grapes_m, ONLY: module_io_grapes_t
  USE NameValueMap_m
  IMPLICIT NONE

  TYPE, EXTENDS(module_io_grapes_t) :: module_io_grapes_CMA_MESO_t

  CONTAINS
    PROCEDURE, PUBLIC, NOPASS :: init_grapes_map_modes
    PROCEDURE, PUBLIC, NOPASS :: input_grapes_sequential
    PROCEDURE, PUBLIC, NOPASS :: input_grapes_hydrometeors
    GENERIC, PUBLIC :: read_value => read_value_1d, read_value_2d, read_value_3d
    GENERIC, PUBLIC :: write_value => write_value_1d, write_value_2d, write_value_3d
    PROCEDURE, PUBLIC :: read_value_1d
    PROCEDURE, PUBLIC :: read_value_2d
    PROCEDURE, PUBLIC :: read_value_3d
    PROCEDURE, PUBLIC :: write_value_1d
    PROCEDURE, PUBLIC :: write_value_2d
    PROCEDURE, PUBLIC :: write_value_3d
    PROCEDURE, PUBLIC :: initialize
  END TYPE

CONTAINS
  SUBROUTINE initialize(this, hgrid, nlFileName, giFileName)
    CLASS(module_io_grapes_CMA_MESO_t)   :: this
    TYPE(domain_t), TARGET      :: hgrid
    CHARACTER(len=*)            :: giFileName, nlFileName

    CALL this%module_io_grapes_t%initialize(hgrid, nlFileName, giFileName)
    this%valName2FileDisp = init_grapes_map_modes(nlFileName, giFileName, hgrid, 'sequential')
  END SUBROUTINE

  SUBROUTINE input_grapes_hydrometeors(HydroFileName, grid)
    USE module_configure
    IMPLICIT NONE
    CHARACTER(*)            :: HydroFileName
    TYPE(domain_t)          :: grid

    INTEGER :: i, j, k

    OPEN(138,FILE=TRIM(HydroFileName),FORM='unformatted',ACCESS='sequential',   &
             convert='big_endian',STATUS='unknown')
    ! gcas only contains qv, qc, qr   
    DO k = grid%kds, grid%kde
      READ(138) ((grid%qc(i,k,j),i=grid%ids,grid%ide),j=grid%jds,grid%jde) ! qv
    END DO
    DO k = grid%kds, grid%kde
      READ(138) ((grid%qc(i,k,j),i=grid%ids,grid%ide),j=grid%jds,grid%jde)
    END DO
    DO k = grid%kds, grid%kde
      READ(138) ((grid%qr(i,k,j),i=grid%ids,grid%ide),j=grid%jds,grid%jde)
    END DO
    DO k = grid%kds, grid%kde
      READ(138) ((grid%qi(i,k,j),i=grid%ids,grid%ide),j=grid%jds,grid%jde)
    END DO
    DO k = grid%kds, grid%kde
      READ(138) ((grid%qs(i,k,j),i=grid%ids,grid%ide),j=grid%jds,grid%jde)
    END DO
    DO k = grid%kds, grid%kde
      READ(138) ((grid%qg(i,k,j),i=grid%ids,grid%ide),j=grid%jds,grid%jde)
    END DO
    print *, 'shape of grid%qc ', shape(grid%qc)

    CLOSE(138)

  END SUBROUTINE input_grapes_hydrometeors
  
  SUBROUTINE input_grapes_sequential(nlFileName, inputfn, grid)
    USE module_configure
    IMPLICIT NONE
    CHARACTER*(*)           :: inputfn, nlFileName
    INTEGER                 :: fh !file handle
    TYPE(domain_t)          :: grid
    INTEGER*4               :: mdate(5)
    INTEGER*8               :: disp
    INTEGER(i_kind)         :: ids, ide, jds, jde, kds, kde, i, j, k
    INTEGER(i_kind)         :: offset_3d, offset_1d_j, offset_1d_k, offset_1d_soil, offset_2d, offset_3d_soil

    ids = grid%ids
    ide = grid%ide
    jds = grid%jds
    jde = grid%jde
    kds = grid%kds
    kde = grid%kde

    CALL initial_config(TRIM(nlFileName), grid%config)
    CALL getLocalIndices(grid)
    CALL alloc_grid_memory(grid)

    disp = 1

    fh = 333 !file handle
    OPEN (fh, file=inputfn, FORM='unformatted', STATUS='unknown', ACCESS='STREAM', ACTION='read')

    ! Read the datetime
    READ (fh, pos=disp) mdate
    ! disp = disp + 20 + 2*4

    WRITE (*, *) "datetime from grapesinput = ", mdate

    grid%config%start_year = mdate(1)
    grid%config%start_month = mdate(2)
    grid%config%start_day = mdate(3)
    grid%config%start_hour = mdate(4)
    grid%config%start_minute = mdate(5)
    grid%config%start_second = 0

    INQUIRE (fh, pos=disp)

    offset_1d_j = (jde - jds + 1) * 4
    offset_1d_k = (kde - kds + 1) * 4
    offset_1d_soil = grid%config%num_soil_layers * 4
    offset_2d = (ide - ids + 1) * (jde - jds + 1) * 4
    offset_3d = (ide - ids + 1) * (jde - jds + 1) * (kde - kds + 1) * 4
    offset_3d_soil = (ide - ids + 1) * (jde - jds + 1) * grid%config%num_soil_layers * 4

    PRINT *, 'Disp is:', disp, SIZE(grid%zz), offset_3d
    ! Read the file data
    CALL skip_vars(disp, 1, offset_3d) ! dthdz
    CALL skip_vars(disp, 1, offset_3d) ! thref
    CALL skip_vars(disp, 1, offset_3d) ! piref
    CALL skip_vars(disp, 1, offset_3d) ! zsx
    CALL skip_vars(disp, 1, offset_3d) ! zsy

    READ (fh, pos=disp) grid%zz; CALL skip_vars(disp, 1, offset_3d); ! zz
    READ (fh, pos=disp) grid%pip; CALL skip_vars(disp, 1, offset_3d); ! pip
    READ (fh, pos=disp) grid%pi; CALL skip_vars(disp, 1, offset_3d); ! pi
    READ (fh, pos=disp) grid%u; CALL skip_vars(disp, 1, offset_3d); ! u
    READ (fh, pos=disp) grid%v; CALL skip_vars(disp, 1, offset_3d); ! v
    CALL skip_vars(disp, 1, offset_3d); ! read (fh, pos=disp) grid%wzet ! wzet
    READ (fh, pos=disp) grid%thp; CALL skip_vars(disp, 1, offset_3d); ! thp
    READ (fh, pos=disp) grid%th; CALL skip_vars(disp, 1, offset_3d); ! th
    READ (fh, pos=disp) grid%q; CALL skip_vars(disp, 1, offset_3d); ! mois_2
    CALL skip_vars(disp, 1, offset_3d); ! resv

    CALL skip_vars(disp, 1, offset_3d_soil); ! tslb
    CALL skip_vars(disp, 1, offset_3d_soil); ! smois

    CALL skip_vars(disp, 1, offset_2d); ! zs
    CALL skip_vars(disp, 1, offset_2d); ! zst
    CALL skip_vars(disp, 1, offset_2d); ! fisx
    CALL skip_vars(disp, 1, offset_2d); ! fisy

    READ (fh, pos=disp) grid%tsk; CALL skip_vars(disp, 1, offset_2d); ! tsk

    CALL skip_vars(disp, 1, offset_2d); ! ps
    CALL skip_vars(disp, 1, offset_2d); ! psea
    CALL skip_vars(disp, 1, offset_2d); ! u10
    CALL skip_vars(disp, 1, offset_2d); ! v10
    CALL skip_vars(disp, 1, offset_2d); ! t2
    CALL skip_vars(disp, 1, offset_2d); ! q2

    READ (fh, pos=disp) grid%ht; CALL skip_vars(disp, 1, offset_2d); ! ht
    CALL skip_vars(disp, 1, offset_2d); ! Tmn
    READ (fh, pos=disp) grid%xland; CALL skip_vars(disp, 1, offset_2d); !xland
    CALL skip_vars(disp, 1, offset_2d); ! lu_index

    READ (fh, pos=disp) grid%snowc; CALL skip_vars(disp, 1, offset_2d); !snowc
    READ (fh, pos=disp) grid%xice; CALL skip_vars(disp, 1, offset_2d); !xice
    CALL skip_vars(disp, 1, offset_2d); ! read (fh, pos=disp) grid%xlat !xlat
    CALL skip_vars(disp, 1, offset_2d); ! read (fh, pos=disp) grid%xlong !xlong
    READ (fh, pos=disp) grid%soil_type; CALL skip_vars(disp, 1, offset_2d); !soil_type
    READ (fh, pos=disp) grid%veg_fraction; CALL skip_vars(disp, 1, offset_2d); !veg_fraction

    CALL skip_vars(disp, 1, offset_2d); ! uslope
    CALL skip_vars(disp, 1, offset_2d); ! vslope
    CALL skip_vars(disp, 1, offset_2d); ! wslope
    CALL skip_vars(disp, 1, offset_2d); ! albedo
    CALL skip_vars(disp, 1, offset_2d); ! snoalb

    CALL skip_vars(disp, 8, offset_1d_j); 
    READ (fh, pos=disp) grid%zh !zh

    PRINT *, 'pos: ', disp, ' grid%zh: ', grid%zh

    CLOSE (fh)
  END SUBROUTINE input_grapes_sequential

  SUBROUTINE skip_vars(disp, nVars, sizeVars)
    IMPLICIT NONE
    INTEGER(8)                :: disp !file handle
    INTEGER(i_kind) :: nVars !< Skips
    INTEGER(i_kind) :: sizeVars !< Skips

    disp = disp + nVars * sizeVars
  END SUBROUTINE

  FUNCTION init_grapes_map_modes(nlFileName, inputfn, grid, mode)
    USE module_configure
    USE NameValueMap_m
    IMPLICIT NONE
    CHARACTER*(*)           :: inputfn, nlFileName, mode
    INTEGER                 :: fh !file handle
    TYPE(domain_t)         :: grid
    TYPE(NameValueMap_t)  :: init_grapes_map_modes
    INTEGER*8               :: disp, offset_3d, i, offset_2d, offset_3d_soil, offset_1d, &
                               head, offset_1d_j, offset_1d_k, offset_1d_soil
    INTEGER(i_kind)         :: ids, ide, jds, jde, kds, kde

    init_grapes_map_modes = create_map_list(100)
    ids = grid%ids
    ide = grid%ide
    jds = grid%jds
    jde = grid%jde
    kds = grid%kds
    kde = grid%kde

    offset_1d_j = (jde - jds + 1) * 4
    offset_1d_k = (kde - kds + 1) * 4
    offset_1d_soil = grid%config%num_soil_layers * 4
    offset_2d = (ide - ids + 1) * (jde - jds + 1) * 4
    offset_3d = (ide - ids + 1) * (jde - jds + 1) * (kde - kds + 1) * 4
    offset_3d_soil = (ide - ids + 1) * (jde - jds + 1) * grid%config%num_soil_layers * 4

    disp = 21

    CALL init_grapes_map_modes%add("dthdz", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("thref", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("piref", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("zsx", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("zsy", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("zz", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("pip", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("pi", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("u", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("v", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("wzet", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("thp", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("th", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("q", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("resv", disp); disp = disp + offset_3d; 
    CALL init_grapes_map_modes%add("tslb", disp); disp = disp + offset_3d_soil; 
    CALL init_grapes_map_modes%add("smois", disp); disp = disp + offset_3d_soil; 
    CALL init_grapes_map_modes%add("zs", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("zst", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("fisx", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("fisy", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("tsk", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("ps", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("psea", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("u10", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("v10", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("t2", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("q2", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("ht", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("tmn", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("xland", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("lu_index", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("snowc", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("xice", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("xlat", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("xlong", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("soil_type", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("vegfra", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("uslope", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("vslope", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("wslope", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("albedo", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("snoalb", disp); disp = disp + offset_2d; 
    CALL init_grapes_map_modes%add("dx", disp); disp = disp + offset_1d_j; 
    CALL init_grapes_map_modes%add("dy", disp); disp = disp + offset_1d_j; 
    CALL init_grapes_map_modes%add("fu", disp); disp = disp + offset_1d_j; 
    CALL init_grapes_map_modes%add("fv", disp); disp = disp + offset_1d_j; 
    CALL init_grapes_map_modes%add("cosv", disp); disp = disp + offset_1d_j; 
    CALL init_grapes_map_modes%add("cosfi", disp); disp = disp + offset_1d_j; 
    CALL init_grapes_map_modes%add("tanu", disp); disp = disp + offset_1d_j; 
    CALL init_grapes_map_modes%add("tanv", disp); disp = disp + offset_1d_j; 
    PRINT *, 'in Map, disp is: ', disp
    CALL init_grapes_map_modes%add("zh", disp); disp = disp + offset_1d_k; 
    CALL init_grapes_map_modes%add("dk", disp); disp = disp + offset_1d_k; 
    CALL init_grapes_map_modes%add("d2k", disp); disp = disp + offset_1d_k; 
    ! CALL init_grapes_map_modes%add("zzs", disp)
    ! CALL init_grapes_map_modes%add("dzs", disp)

  END FUNCTION init_grapes_map_modes

  SUBROUTINE read_value_3d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_CMA_MESO_t) :: this
    INTEGER                 :: fh !file handle
    INTEGER*8               :: disp
    INTEGER(i_kind)         :: ids, ide, jds, jde, kds, kde
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE :: res
    TYPE(domain_t)          :: hgrid

    ids = hgrid%ids
    ide = hgrid%ide
    jds = hgrid%jds
    jde = hgrid%jde
    kds = hgrid%kds
    kde = hgrid%kde

    fh = 333 !file handle
    OPEN (fh, file=this%giFileName, FORM='unformatted', STATUS='unknown', &
          ACCESS='STREAM', ACTION='read')

    IF (ALLOCATED(res)) DEALLOCATE (res)
    ALLOCATE (res(ids:ide, kds:kde, jds:jde))
    disp = this%valName2FileDisp%get_value_by_name(vname)
    READ (fh, pos=disp) res

    CLOSE (fh)

  END SUBROUTINE

  SUBROUTINE read_value_2d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_CMA_MESO_t) :: this
    INTEGER                 :: fh !file handle
    INTEGER*8               :: disp
    INTEGER(i_kind)         :: ids, ide, jds, jde, kds, kde
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:, :), ALLOCATABLE :: res
    TYPE(domain_t)          :: hgrid

    ids = hgrid%ids
    ide = hgrid%ide
    jds = hgrid%jds
    jde = hgrid%jde
    kds = hgrid%kds
    kde = hgrid%kde

    fh = 333 !file handle
    OPEN (fh, file=this%giFileName, FORM='unformatted', STATUS='unknown', &
          ACCESS='STREAM', ACTION='read')

    IF (ALLOCATED(res)) DEALLOCATE (res)
    ALLOCATE (res(ids:ide, jds:jde))
    disp = this%valName2FileDisp%get_value_by_name(vname)
    READ (fh, pos=disp) res
    CLOSE (fh)

  END SUBROUTINE

  SUBROUTINE read_value_1d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_CMA_MESO_t) :: this
    INTEGER                 :: fh !file handle
    INTEGER*8               :: disp
    INTEGER(i_kind)         :: ids, ide, jds, jde, kds, kde
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:), ALLOCATABLE :: res
    TYPE(domain_t)          :: hgrid

    ids = hgrid%ids
    ide = hgrid%ide
    jds = hgrid%jds
    jde = hgrid%jde
    kds = hgrid%kds
    kde = hgrid%kde

    fh = 333 !file handle
    OPEN (fh, file=this%giFileName, FORM='unformatted', STATUS='unknown', &
          ACCESS='STREAM', ACTION='read')

    IF (ALLOCATED(res)) DEALLOCATE (res)
    ALLOCATE (res(kds:kde))
    disp = this%valName2FileDisp%get_value_by_name(vname)
    READ (fh, pos=disp) res

    CLOSE (fh)
  END SUBROUTINE

  SUBROUTINE write_value_3d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_CMA_MESO_t) :: this
    INTEGER                 :: fh !file handle
    INTEGER*8               :: disp
    INTEGER(i_kind)         :: ids, ide, jds, jde, kds, kde
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:, :, :) :: res
    TYPE(domain_t)          :: hgrid

    ids = hgrid%ids
    ide = hgrid%ide
    jds = hgrid%jds
    jde = hgrid%jde
    kds = hgrid%kds
    kde = hgrid%kde

    fh = 333 !file handle
    OPEN (fh, file=this%giFileName, FORM='unformatted', STATUS='unknown', &
          ACCESS='STREAM', ACTION='write')

    disp = this%valName2FileDisp%get_value_by_name(vname)
    PRINT *, 'disp is: ', disp, TRIM(this%giFileName)
    WRITE (fh, pos=disp) res

    CLOSE (fh)

  END SUBROUTINE

  SUBROUTINE write_value_2d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_CMA_MESO_t) :: this
    INTEGER                 :: fh !file handle
    INTEGER*8               :: disp
    INTEGER(i_kind)         :: ids, ide, jds, jde, kds, kde
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:, :) :: res
    TYPE(domain_t)          :: hgrid

    ids = hgrid%ids
    ide = hgrid%ide
    jds = hgrid%jds
    jde = hgrid%jde
    kds = hgrid%kds
    kde = hgrid%kde

    fh = 333 !file handle
    OPEN (fh, file=this%giFileName, FORM='unformatted', STATUS='unknown', &
          ACCESS='STREAM', ACTION='write')

    disp = this%valName2FileDisp%get_value_by_name(vname)
    PRINT *, disp, "disp"
    WRITE (fh, pos=disp) res

    CLOSE (fh)
  END SUBROUTINE

  SUBROUTINE write_value_1d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_CMA_MESO_t) :: this
    INTEGER                 :: fh !file handle
    INTEGER*8               :: disp
    INTEGER(i_kind)         :: ids, ide, jds, jde, kds, kde
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:) :: res
    TYPE(domain_t)          :: hgrid

    ids = hgrid%ids
    ide = hgrid%ide
    jds = hgrid%jds
    jde = hgrid%jde
    kds = hgrid%kds
    kde = hgrid%kde

    fh = 333 !file handle
    OPEN (fh, file=this%giFileName, FORM='unformatted', STATUS='unknown', &
          ACCESS='STREAM', ACTION='write')

    disp = this%valName2FileDisp%get_value_by_name(vname)
    WRITE (fh, pos=disp) res

    CLOSE (fh)
  END SUBROUTINE

END MODULE
