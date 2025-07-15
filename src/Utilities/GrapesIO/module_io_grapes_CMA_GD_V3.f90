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
!  Modified by Zilong Qin (zilong.qin@gmail.com), 2023/6/16, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!!===================================================================
!> @brief
!! # GRAPES IO Module
!!
!!  *This module defines the subroutines for GRAPES IO*
!! @author Sanshan Tu
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
!!===================================================================
MODULE module_io_grapes_CMA_GD_V3_m
  USE module_domain
  USE module_configure
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE module_io_grapes_m, ONLY: module_io_grapes_t
  USE NameValueMap_m

  TYPE, EXTENDS(module_io_grapes_t) :: module_io_grapes_CMA_GD_V3_t

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
    CLASS(module_io_grapes_CMA_GD_V3_t)   :: this
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
    print *, 'grids are: ', grid%kds, grid%kde, grid%ids,grid%ide, grid%jds,grid%jde
    DO k = grid%kds+1, grid%kde
      READ(138) ((grid%qc(i,k,j),i=grid%ids,grid%ide),j=grid%jds,grid%jde) ! qv
    END DO
    DO k = grid%kds+1, grid%kde
      READ(138) ((grid%qc(i,k,j),i=grid%ids,grid%ide),j=grid%jds,grid%jde)
    END DO
    DO k = grid%kds+1, grid%kde
      READ(138) ((grid%qr(i,k,j),i=grid%ids,grid%ide),j=grid%jds,grid%jde)
    END DO
    DO k = grid%kds+1, grid%kde
      READ(138) ((grid%qi(i,k,j),i=grid%ids,grid%ide),j=grid%jds,grid%jde)
    END DO
    DO k = grid%kds+1, grid%kde
      READ(138) ((grid%qs(i,k,j),i=grid%ids,grid%ide),j=grid%jds,grid%jde)
    END DO
    DO k = grid%kds+1, grid%kde
      READ(138) ((grid%qg(i,k,j),i=grid%ids,grid%ide),j=grid%jds,grid%jde)
    END DO
    print *, 'shape of grid%qc ', shape(grid%qc)
    
    grid%qc(grid%ids:grid%ids+grid%config%spec_bdy_width-1,:,:) = 0.0D0
    grid%qc(grid%ide-grid%config%spec_bdy_width+1:grid%ide,:,:) = 0.0D0
    grid%qc(:,:,grid%jds:grid%jds+grid%config%spec_bdy_width-1) = 0.0D0
    grid%qc(:,:,grid%jde-grid%config%spec_bdy_width+1:grid%jde) = 0.0D0

    grid%qr(grid%ids:grid%ids+grid%config%spec_bdy_width-1,:,:) = 0.0D0
    grid%qr(grid%ide-grid%config%spec_bdy_width+1:grid%ide,:,:) = 0.0D0
    grid%qr(:,:,grid%jds:grid%jds+grid%config%spec_bdy_width-1) = 0.0D0
    grid%qr(:,:,grid%jde-grid%config%spec_bdy_width+1:grid%jde) = 0.0D0

    grid%qi(grid%ids:grid%ids+grid%config%spec_bdy_width-1,:,:) = 0.0D0
    grid%qi(grid%ide-grid%config%spec_bdy_width+1:grid%ide,:,:) = 0.0D0
    grid%qi(:,:,grid%jds:grid%jds+grid%config%spec_bdy_width-1) = 0.0D0
    grid%qi(:,:,grid%jde-grid%config%spec_bdy_width+1:grid%jde) = 0.0D0

    grid%qs(grid%ids:grid%ids+grid%config%spec_bdy_width-1,:,:) = 0.0D0
    grid%qs(grid%ide-grid%config%spec_bdy_width+1:grid%ide,:,:) = 0.0D0
    grid%qs(:,:,grid%jds:grid%jds+grid%config%spec_bdy_width-1) = 0.0D0
    grid%qs(:,:,grid%jde-grid%config%spec_bdy_width+1:grid%jde) = 0.0D0

    grid%qg(grid%ids:grid%ids+grid%config%spec_bdy_width-1,:,:) = 0.0D0
    grid%qg(grid%ide-grid%config%spec_bdy_width+1:grid%ide,:,:) = 0.0D0
    grid%qg(:,:,grid%jds:grid%jds+grid%config%spec_bdy_width-1) = 0.0D0
    grid%qg(:,:,grid%jde-grid%config%spec_bdy_width+1:grid%jde) = 0.0D0
    CLOSE(138)

  END SUBROUTINE input_grapes_hydrometeors

  SUBROUTINE input_grapes_sequential(nlFileName, inputfn, grid)
    USE module_configure
    IMPLICIT NONE
    CHARACTER*(*)           :: inputfn, nlFileName
    INTEGER                 :: fh !file handle
    TYPE(domain_t)          :: grid
    INTEGER*4               :: mdate(5)
    INTEGER*8              :: disp

    CALL initial_config(TRIM(nlFileName), grid%config)
    CALL getLocalIndices(grid)
    CALL alloc_grid_memory(grid)

    disp = 5

    fh = 333 !file handle
    OPEN (fh, file=inputfn, FORM='unformatted', STATUS='unknown', ACCESS='STREAM', ACTION='read', convert="big_endian")

    ! Read the datetime
    READ (fh, pos=disp) mdate
    disp = disp + 20 + 2 * 4

    WRITE (*, *) "datetime from grapesinput = ", mdate

    grid%config%start_year = mdate(1)
    grid%config%start_month = mdate(2)
    grid%config%start_day = mdate(3)
    grid%config%start_hour = mdate(4)
    grid%config%start_minute = mdate(5)
    grid%config%start_second = 0

    ! Read the file data
    disp = disp - 8
    CALL skip_records(fh, disp, 9); READ (fh, pos=disp + 8) grid%zz ! zz
    CALL skip_records(fh, disp); READ (fh, pos=disp + 8) grid%pip ! pip
    CALL skip_records(fh, disp); READ (fh, pos=disp + 8) grid%pi ! pi
    CALL skip_records(fh, disp); READ (fh, pos=disp + 8) grid%u ! u
    CALL skip_records(fh, disp); READ (fh, pos=disp + 8) grid%v ! v
    CALL skip_records(fh, disp); ! read (fh, pos=disp + 8) grid%wzet ! wzet
    CALL skip_records(fh, disp); READ (fh, pos=disp + 8) grid%thp ! thp
    CALL skip_records(fh, disp); READ (fh, pos=disp + 8) grid%th ! th
    CALL skip_records(fh, disp); READ (fh, pos=disp + 8) grid%q ! mois_2

    CALL skip_records(fh, disp, 2)
    CALL skip_records(fh, disp, 2)
    CALL skip_records(fh, disp, 5)

    CALL skip_records(fh, disp); READ (fh, pos=disp + 8) grid%tsk ! tsk
    CALL skip_records(fh, disp); ! read (fh, pos=disp + 8) grid%ps ! ps

    CALL skip_records(fh, disp)      !  psea
    CALL skip_records(fh, disp)      !  u10
    CALL skip_records(fh, disp)      !  v10
    CALL skip_records(fh, disp)      !  t2
    CALL skip_records(fh, disp)      !  q2
    CALL skip_records(fh, disp); READ (fh, pos=disp + 8) grid%ht ! ht

    CALL skip_records(fh, disp)     ! Tmn
    CALL skip_records(fh, disp); READ (fh, pos=disp + 8) grid%xland !xland
    CALL skip_records(fh, disp)
    CALL skip_records(fh, disp); READ (fh, pos=disp + 8) grid%snowc !snowc
    CALL skip_records(fh, disp); READ (fh, pos=disp + 8) grid%xice !xice
    CALL skip_records(fh, disp); ! read (fh, pos=disp + 8) grid%xlat !xlat
    CALL skip_records(fh, disp); ! read (fh, pos=disp + 8) grid%xlong !xlong
    CALL skip_records(fh, disp); READ (fh, pos=disp + 8) grid%soil_type !soil_type
    CALL skip_records(fh, disp); READ (fh, pos=disp + 8) grid%veg_fraction !veg_fraction

    CALL skip_records(fh, disp, 4)
    CALL skip_records(fh, disp, 13) !for CMA-GD v3.1
    ! CALL skip_records(fh, disp, 22) !for CMA-GD v3.2
    CALL skip_records(fh, disp, 10)
    CALL skip_records(fh, disp); 
    READ (fh, pos=disp + 8) grid%zh !zh
    ! PRINT *, 'in read, zh disp is: ', disp, grid%zh
    ! CALL skip_records(fh, disp);
    ! read (fh, pos=disp + 8) grid%zh !zh
    ! PRINT *, 'in read, zh disp is: ', disp, grid%zh
    ! CALL skip_records(fh, disp);
    ! read (fh, pos=disp + 8) grid%zh !zh
    ! PRINT *, 'in read, zh disp is: ', disp, grid%zh

    CLOSE (fh)
  END SUBROUTINE input_grapes_sequential

  SUBROUTINE skip_records(fh, disp, nVars)
    IMPLICIT NONE
    INTEGER                   :: fh !file handle
    INTEGER(8)                :: disp !file handle
    INTEGER(i_kind), OPTIONAL :: nVars !< Skips
    INTEGER :: nSkip, i
    INTEGER*4  :: trailMarker = 0, leadMarker

    IF (PRESENT(nVars)) THEN
      nSkip = nVars
    ELSE
      nSkip = 1
    END IF

    DO i = 1, nSkip
      READ (fh, pos=disp) trailMarker !xland
      disp = disp + 4
      READ (fh, pos=disp) leadMarker !xland
      PRINT *, 'trailMarker', trailMarker, 'leadMarker', leadMarker
      disp = disp + 4 + leadMarker
    END DO

  END SUBROUTINE

  FUNCTION init_grapes_map(nlFileName, inputfn, grid)
    USE module_configure
    USE NameValueMap_m
    IMPLICIT NONE
    CHARACTER*(*)           :: inputfn, nlFileName
    INTEGER                 :: fh !file handle
    TYPE(domain_t)         :: grid
    TYPE(NameValueMap_t)  :: init_grapes_map
    INTEGER*8               :: disp, offset_3d, i, offset_2d, offset_3d_soil, offset_1d
    INTEGER(i_kind)         :: ids, ide, jds, jde, kds, kde

    fh = 333 !file handle
    OPEN (fh, file=inputfn, FORM='unformatted', STATUS='unknown', ACCESS='STREAM', ACTION='read', convert="big_endian")

    init_grapes_map = create_map_list(100)

    disp = 5
    disp = disp + 20 + 2 * 4
    disp = disp - 8

    CALL init_grapes_map%add("dthdz", disp + 8)
    CALL skip_records(fh, disp, 3); CALL init_grapes_map%add("thref", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("piref", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("zsx", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("zsy", disp + 8)
    CALL skip_records(fh, disp, 3); CALL init_grapes_map%add("zz", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("pip", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("pi", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("u", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("v", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("wzet", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("thp", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("th", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("q", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("resv", disp + 8)
    CALL skip_records(fh, disp); 
    CALL skip_records(fh, disp, 2); CALL init_grapes_map%add("tslb", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("smois", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("zs", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("zst", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("fisx", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("fisy", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("tsk", disp + 8)

    CALL skip_records(fh, disp); CALL init_grapes_map%add("ps", disp + 8)
    CALL skip_records(fh, disp, 6); CALL init_grapes_map%add("ht", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("tmn", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("xland", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("lu_index", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("snowc", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("xice", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("xlat", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("xlong", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("soil_type", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("vegfra", disp + 8)

    CALL skip_records(fh, disp, 5)
    CALL skip_records(fh, disp, 13); CALL init_grapes_map%add("dx", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("dy", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("fu", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("fv", disp + 8)
    CALL skip_records(fh, disp, 3); CALL init_grapes_map%add("cosv", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("cosfi", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("tanu", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("tanv", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("zh", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("dk", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("d2k", disp + 8)
    CALL skip_records(fh, disp, 4); CALL init_grapes_map%add("zzs", disp + 8)
    CALL skip_records(fh, disp); CALL init_grapes_map%add("dzs", disp + 8)

    CLOSE (fh)
  END FUNCTION init_grapes_map

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

    ids = grid%ids
    ide = grid%ide
    jds = grid%jds
    jde = grid%jde
    kds = grid%kds
    kde = grid%kde

    IF (mode == 'sequential') head = 8
    disp = 5
    offset_1d_j = (jde - jds + 1) * 4 + head
    PRINT *, "offset_1d_j", offset_1d_j
    offset_1d_k = (kde - kds + 1) * 4 + head
    offset_1d_soil = grid%config%num_soil_layers * 4 + head
    offset_2d = (ide - ids + 1) * (jde - jds + 1) * 4 + head
    offset_3d = (ide - ids + 1) * (jde - jds + 1) * (kde - kds + 1) * 4 + head
    offset_3d_soil = (ide - ids + 1) * (jde - jds + 1) * grid%config%num_soil_layers * 4 + head

    fh = 333 !file handle
    OPEN (fh, file=inputfn, FORM='unformatted', STATUS='unknown', ACCESS='STREAM', ACTION='read', convert="big_endian")

    init_grapes_map_modes = create_map_list(100)

    disp = 5
    disp = disp + 20 + 2 * 4
    disp = disp - 8

    CALL init_grapes_map_modes%add("dthdz", disp + 8)

    disp = disp + offset_3d * 3; CALL init_grapes_map_modes%add("thref", disp + 8)
    disp = disp + offset_3d; CALL init_grapes_map_modes%add("piref", disp + 8)
    disp = disp + offset_3d; CALL init_grapes_map_modes%add("zsx", disp + 8)
    disp = disp + offset_3d; CALL init_grapes_map_modes%add("zsy", disp + 8)
    disp = disp + offset_3d * 3; CALL init_grapes_map_modes%add("zz", disp + 8)
    disp = disp + offset_3d; CALL init_grapes_map_modes%add("pip", disp + 8)
    disp = disp + offset_3d; CALL init_grapes_map_modes%add("pi", disp + 8)
    disp = disp + offset_3d; CALL init_grapes_map_modes%add("u", disp + 8)
    disp = disp + offset_3d; CALL init_grapes_map_modes%add("v", disp + 8)
    disp = disp + offset_3d; CALL init_grapes_map_modes%add("wzet", disp + 8)
    disp = disp + offset_3d; CALL init_grapes_map_modes%add("thp", disp + 8)
    disp = disp + offset_3d; CALL init_grapes_map_modes%add("th", disp + 8)
    disp = disp + offset_3d; CALL init_grapes_map_modes%add("q", disp + 8)
    disp = disp + offset_3d; CALL init_grapes_map_modes%add("resv", disp + 8)

    disp = disp + offset_3d
    disp = disp + offset_3d_soil * 2; CALL init_grapes_map_modes%add("tslb", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("smois", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("zs", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("zst", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("fisx", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("fisy", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("tsk", disp + 8)

    disp = disp + offset_2d; CALL init_grapes_map_modes%add("ps", disp + 8)
    disp = disp + offset_2d * 6; CALL init_grapes_map_modes%add("ht", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("tmn", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("xland", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("lu_index", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("snowc", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("xice", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("xlat", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("xlong", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("soil_type", disp + 8)
    disp = disp + offset_2d; CALL init_grapes_map_modes%add("vegfra", disp + 8)
    disp = disp + offset_2d * 5
    disp = disp + offset_2d * 13; CALL init_grapes_map_modes%add("dx", disp + 8)
    disp = disp + offset_1d_j; CALL init_grapes_map_modes%add("dy", disp + 8)
    disp = disp + offset_1d_j; CALL init_grapes_map_modes%add("fu", disp + 8)
    disp = disp + offset_1d_j; CALL init_grapes_map_modes%add("fv", disp + 8)
    disp = disp + offset_1d_j * 3; CALL init_grapes_map_modes%add("cosv", disp + 8)
    disp = disp + offset_1d_j; CALL init_grapes_map_modes%add("cosfi", disp + 8)
    disp = disp + offset_1d_j; CALL init_grapes_map_modes%add("tanu", disp + 8)
    disp = disp + offset_1d_j; CALL init_grapes_map_modes%add("tanv", disp + 8)
    disp = disp + offset_1d_j; 
    PRINT *, 'in Map, disp is: ', disp
    CALL init_grapes_map_modes%add("zh", disp + 8)
    disp = disp + offset_1d_k; CALL init_grapes_map_modes%add("dk", disp + 8)
    disp = disp + offset_1d_k; CALL init_grapes_map_modes%add("d2k", disp + 8)
    disp = disp + offset_1d_k + offset_1d_soil * 2 + 6 * 4 + head; CALL init_grapes_map_modes%add("zzs", disp + 8)
    disp = disp + 4 + head; CALL init_grapes_map_modes%add("dzs", disp + 8)

    CLOSE (fh)
  END FUNCTION init_grapes_map_modes

  SUBROUTINE read_value_3d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_CMA_GD_V3_t) :: this
    INTEGER                 :: fh !file handle
    INTEGER*8               :: disp
    INTEGER(i_kind)         :: ids, ide, jds, jde, kds, kde
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE :: res
    TYPE(domain_t)          :: hgrid

    ! Yuanfu Xie added an error check in opening the output file:
    INTEGER(i_kind) :: ierror

    ids = hgrid%ids
    ide = hgrid%ide
    jds = hgrid%jds
    jde = hgrid%jde
    kds = hgrid%kds
    kde = hgrid%kde

    fh = 333 !file handle
    ! Yuanfu Xie added an error check in opening the output file: 2023/10/27
    OPEN (fh, file=this%giFileName, FORM='unformatted', STATUS='unknown', &
          ACCESS='STREAM', ACTION='read', convert="big_endian", iostat=ierror)
    IF (ierror /= 0) THEN
      WRITE (*, 1) TRIM(this%giFileName), ierror
1     FORMAT('You need this file: ', A, ' for grid info, please copy the proper file to your output dir and rerun!', I5)
      STOP
    END IF

    IF (ALLOCATED(res)) DEALLOCATE (res)
    ALLOCATE (res(ids:ide, kds:kde, jds:jde))
    disp = this%valName2FileDisp%get_value_by_name(vname)
    READ (fh, pos=disp) res

    CLOSE (fh)

  END SUBROUTINE

  SUBROUTINE read_value_2d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_CMA_GD_V3_t) :: this
    INTEGER                 :: fh !file handle
    INTEGER*8               :: disp
    INTEGER(i_kind)         :: ids, ide, jds, jde, kds, kde
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:, :), ALLOCATABLE :: res
    TYPE(domain_t)          :: hgrid

    ! Yuanfu Xie added an error check in opening the output file:
    INTEGER(i_kind) :: ierror

    ids = hgrid%ids
    ide = hgrid%ide
    jds = hgrid%jds
    jde = hgrid%jde
    kds = hgrid%kds
    kde = hgrid%kde

    fh = 333 !file handle
    ! Yuanfu Xie added an error check in opening the output file: 2023/10/27
    OPEN (fh, file=this%giFileName, FORM='unformatted', STATUS='unknown', &
          ACCESS='STREAM', ACTION='read', convert="big_endian", iostat=ierror)
    IF (ierror /= 0) THEN
      WRITE (*, 1) TRIM(this%giFileName), ierror
1     FORMAT('You need this file: ', A, ' for grid info, please copy the proper file to your output dir and rerun!', I5)
      STOP
    END IF

    IF (ALLOCATED(res)) DEALLOCATE (res)
    ALLOCATE (res(ids:ide, jds:jde))
    disp = this%valName2FileDisp%get_value_by_name(vname)
    READ (fh, pos=disp) res
    CLOSE (fh)

  END SUBROUTINE

  SUBROUTINE read_value_1d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_CMA_GD_V3_t) :: this
    INTEGER                 :: fh !file handle
    INTEGER*8               :: disp
    INTEGER(i_kind)         :: ids, ide, jds, jde, kds, kde
    CHARACTER*(*)           :: vname
    REAL(r_single), DIMENSION(:), ALLOCATABLE :: res
    TYPE(domain_t)          :: hgrid

    ! Yuanfu Xie added an error check in opening the output file:
    INTEGER(i_kind) :: ierror

    ids = hgrid%ids
    ide = hgrid%ide
    jds = hgrid%jds
    jde = hgrid%jde
    kds = hgrid%kds
    kde = hgrid%kde

    fh = 333 !file handle
    ! Yuanfu Xie added an error check in opening the output file: 2023/10/27
    OPEN (fh, file=this%giFileName, FORM='unformatted', STATUS='unknown', &
          ACCESS='STREAM', ACTION='read', convert="big_endian", iostat=ierror)
    IF (ierror /= 0) THEN
      WRITE (*, 1) TRIM(this%giFileName), ierror
1     FORMAT('You need this file: ', A, ' for grid info, please copy the proper file to your output dir and rerun!', I5)
      STOP
    END IF

    IF (ALLOCATED(res)) DEALLOCATE (res)
    ALLOCATE (res(kds:kde))
    disp = this%valName2FileDisp%get_value_by_name(vname)
    READ (fh, pos=disp) res

    CLOSE (fh)
  END SUBROUTINE

  SUBROUTINE write_value_3d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_CMA_GD_V3_t) :: this
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
          ACCESS='STREAM', ACTION='write', convert="big_endian")

    disp = this%valName2FileDisp%get_value_by_name(vname)
    PRINT *, 'disp is: ', disp, TRIM(this%giFileName)
    WRITE (fh, pos=disp) res

    CLOSE (fh)

  END SUBROUTINE

  SUBROUTINE write_value_2d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_CMA_GD_V3_t) :: this
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
          ACCESS='STREAM', ACTION='write', convert="big_endian")

    disp = this%valName2FileDisp%get_value_by_name(vname)
    PRINT *, disp, "disp"
    WRITE (fh, pos=disp) res

    CLOSE (fh)
  END SUBROUTINE

  SUBROUTINE write_value_1d(this, vname, res, hgrid)
    IMPLICIT NONE
    CLASS(module_io_grapes_CMA_GD_V3_t) :: this
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
          ACCESS='STREAM', ACTION='write', convert="big_endian")

    disp = this%valName2FileDisp%get_value_by_name(vname)
    WRITE (fh, pos=disp) res

    CLOSE (fh)
  END SUBROUTINE

END MODULE
