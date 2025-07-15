!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 1.0
! HISTORY           :
!   Created by Jiongming Pang, 2023/8/8, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

MODULE GrapesPostvarIO_m
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE parameters_m, ONLY: spec_heat_const_pres, dry_air_gas_const, surface_ref_pres, degree2radian
  USE Interp1D_m, ONLY: interp1d_3D_idx2_single, interp1d_3D_idx2
  ! USE Filter_m, ONLY: smoothField, guidedfilter
  USE AdvanceTime_m

  IMPLICIT NONE

  TYPE GrapesPostvarIO_t
    CHARACTER(LEN=1024), ALLOCATABLE  ::  postvarFns(:)
    CHARACTER(LEN=1024)               ::  postvarCtlFn
    CHARACTER(LEN=1024)               ::  DEMFileName
    INTEGER(i_kind)                   ::  num_files   ! the numbers of background files
    REAL(r_kind), ALLOCATABLE       ::  lat3D(:, :, :), lon3D(:, :, :), lat2D(:, :), lon2D(:, :), lat1D(:), lon1D(:)
    REAL(r_single), ALLOCATABLE       ::  zRHght(:, :, :), &
                                         zs(:, :)              ! geo-height
    REAL(r_kind), ALLOCATABLE       ::  zRHght_s(:, :, :), topo(:, :)
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE ::  ps, &
                                                       psl, &
                                                       ts, &
                                                       rh2, &
                                                       t2m, &
                                                       u10m, &
                                                       v10m, &
                                                       q2m, &
                                                       psfc, &
                                                       rainc, &
                                                       tmp2d
    REAL(r_single), DIMENSION(:, :, :, :), ALLOCATABLE :: uwnd, vwnd, temp, H, pres, tmp3d, qvapor

    INTEGER(i_kind) :: idn, jdn, kdn, tdn
    ! INTEGER(i_kind) :: ids, ide, jds, jde, kds, kde
    INTEGER(i_kind), ALLOCATABLE :: time_unix(:), time_gmt(:, :)

  CONTAINS
    FINAL :: destructor
    PROCEDURE :: get_domain_from_DEM
    PROCEDURE, PRIVATE :: get_domain
    PROCEDURE, PRIVATE :: read_value
  END TYPE GrapesPostvarIO_t

  INTERFACE GrapesPostvarIO_t
    PROCEDURE :: constructor
  END INTERFACE

CONTAINS
  FUNCTION constructor(postvarFns, postvarCtlFn, DEMFileName, mdate) RESULT(this)
    TYPE(GrapesPostvarIO_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: postvarFns(:), postvarCtlFn, DEMFileName
    INTEGER(i_kind), INTENT(IN) :: mdate(:, :)

    this%postvarFns = postvarFns
    this%postvarCtlFn = TRIM(postvarCtlFn)
    this%DEMFileName = TRIM(DEMFileName)
    this%time_gmt = mdate

    ! get domain info from ctl files
    CALL this%get_domain()

    ! read vars' value from modelvar and HeightField files
    CALL this%read_value()

  END FUNCTION constructor

  SUBROUTINE get_domain_from_DEM(this, latTopo, lonTopo, valueTopo, NTopo)
    IMPLICIT NONE
    CLASS(GrapesPostvarIO_t)  :: this
    INTEGER(i_kind) :: n, NTopo
    ! REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE :: dem
    REAL(r_kind), ALLOCATABLE :: latTopo(:), lonTopo(:), valueTopo(:)

    INTEGER :: fid, iostat, i, j, k
    CHARACTER(len=256) :: line
    CHARACTER(len=256) :: tmp

    n = 3601

    PRINT *, "************************ get_domain_from_DEM ************************"

    ! 打开文件
    OPEN (fid, file=TRIM(this%DEMFileName), status='old', action='read', iostat=iostat)
    IF (iostat /= 0) THEN
      PRINT *, 'Error opening file'
      STOP
    END IF

    ! skip the first line
    READ (fid, '(A)', iostat=iostat) line
    PRINT *, "line:", line

    ALLOCATE (latTopo(n * n), lonTopo(n * n), valueTopo(n * n))

    !按行读取文件
    DO i = 1, n
      DO j = 1, n
        ! 读取新的行
        READ (fid, '(A)', iostat=iostat) line
        ! 解析每行的内容
        READ (line, *, iostat=iostat) lonTopo((i - 1) * n + j), latTopo((i - 1) * n + j), valueTopo((i - 1) * n + j)
        IF (iostat /= 0) THEN
          PRINT *, 'Error parsing line ', i, ' data point ', j
          STOP
        END IF
      END DO
    END DO

    ! 关闭文件
    CLOSE (fid)

    NTopo = n * n

    ! ! 打印一些数据以验证读取结果
    ! PRINT *, 'dem(1,1,1):', dem(1,1,1)
    ! PRINT *, 'dem(1,1,2):', dem(1,1,2)
    ! PRINT *, 'dem(1,2,1):', dem(1,2,1)
    ! PRINT *, 'dem(1, 3601, 3601):', dem(1, 3601, 3601)
    ! PRINT *, 'dem(2, 3601, 3601):', dem(2, 3601, 3601)
    ! PRINT *, 'dem(3, 3601, 3601):', dem(3, 3601, 3601)

  END SUBROUTINE get_domain_from_DEM

  SUBROUTINE get_domain(this)
    IMPLICIT NONE
    CLASS(GrapesPostvarIO_t) :: this
    INTEGER           :: fid  !file id
    INTEGER           :: i
    REAL(r_kind)      :: lon_s, lon_int, lat_s, lat_int
    CHARACTER(1024)   :: tmp

    this%num_files = UBOUND(this%postvarFns, 1)
    this%tdn = this%num_files

    ALLOCATE (this%time_unix(this%tdn))
    DO i = 1, this%tdn
      PRINT *, "GMT Time: ", this%time_gmt(:, i)
      CALL Time_GMT_to_Unix(this%time_gmt(:, i), this%time_unix(i))
      PRINT *, "Unix Time: ", this%time_unix(i)
    END DO

    fid = 120    !file id
    OPEN (fid, FILE=TRIM(this%postvarCtlFn), ACTION='read')

    DO i = 1, 3
      READ (fid, *)
    END DO

    ! initial
    this%idn = 0
    this%jdn = 0
    this%kdn = 0

    READ (fid, *) tmp, this%idn, tmp, lon_s, lon_int   ! get longitute info
    READ (fid, *) tmp, this%jdn, tmp, lat_s, lat_int   ! get latitute info
    PRINT *, 'idn, lon_s, lon_int:', this%idn, lon_s, lon_int
    PRINT *, 'jdn, lat_s, lat_int:', this%jdn, lat_s, lat_int

    ! get the amount of z-level
    READ (fid, *) tmp, this%kdn
    PRINT *, 'kdn:', this%kdn

    ! generate lon and lat,3D-array, 2D-array, and 1D-array
    ALLOCATE (this%lon3D(this%idn, this%jdn, this%kdn), this%lat3D(this%idn, this%jdn, this%kdn))
    ALLOCATE (this%lon2D(this%idn, this%jdn), this%lat2D(this%idn, this%jdn))
    ALLOCATE (this%lon1D(this%idn), this%lat1D(this%jdn))
    DO i = 1, this%idn
      this%lon2D(i, :) = lon_s + lon_int * (i - 1)
      this%lon1D(i) = lon_s + lon_int * (i - 1)
    END DO
    DO i = 1, this%jdn
      this%lat2D(:, i) = lat_s + lat_int * (i - 1)
      this%lat1D(i) = lat_s + lat_int * (i - 1)
    END DO
    DO i = 1, this%kdn
      this%lat3D(:, :, i) = this%lat2D
      this%lon3D(:, :, i) = this%lon2D
    END DO

    PRINT *, 'min and max of lat:', MINVAL(this%lat2D(:, :)), MAXVAL(this%lat2D(:, :))
    PRINT *, 'min and max of lon:', MINVAL(this%lon2D(:, :)), MAXVAL(this%lon2D(:, :))

    ! convert degree to radian

    this%lon3D = this%lon3D * degree2radian
    this%lat3D = this%lat3D * degree2radian
    this%lon2D = this%lon2D * degree2radian
    this%lat2D = this%lat2D * degree2radian
    this%lon1D = this%lon1D * degree2radian
    this%lat1D = this%lat1D * degree2radian

    CLOSE (fid)

    PRINT *, 'idn, jdn, kdn, tdn:', this%idn, this%jdn, this%kdn, this%tdn

  END SUBROUTINE get_domain

  SUBROUTINE read_value(this)
    IMPLICIT NONE
    CLASS(GrapesPostvarIO_t)  :: this
    INTEGER                   :: fid   ! file id
    INTEGER(i_kind)           :: idn, jdn, kdn, nt
    INTEGER(i_kind)           :: i, j, k, t, skip, skips, iprec

    idn = this%idn
    jdn = this%jdn
    kdn = this%kdn
    nt = this%tdn

    ! ALLOCATE (this%ps(idn, jdn, nt))
    ! ALLOCATE (this%psl(idn, jdn, nt))
    ! ALLOCATE (this%ts(idn, jdn, nt))
    ! ALLOCATE (this%rh2(idn, jdn, nt))
    ! ALLOCATE (this%t2m(idn, jdn, nt))
    ! ALLOCATE (this%u10m(idn, jdn, nt))
    ! ALLOCATE (this%v10m(idn, jdn, nt))
    ! ALLOCATE (this%psfc(idn, jdn, nt))
    ! ALLOCATE (this%tmp2d(idn, jdn, nt))
    ! ALLOCATE (this%zs(idn, jdn))
    ALLOCATE (this%uwnd(idn, jdn, kdn, nt))
    ALLOCATE (this%vwnd(idn, jdn, kdn, nt))
    ALLOCATE (this%temp(idn, jdn, kdn, nt))
    ALLOCATE (this%pres(idn, jdn, kdn, nt))
    ALLOCATE (this%tmp3d(idn, jdn, kdn, nt))
    ALLOCATE (this%H(idn, jdn, kdn, nt))
    ALLOCATE (this%qvapor(idn, jdn, kdn, nt))

    fid = 110   ! file id
    PRINT *, "idn, jdn,", idn, jdn
    DO t = 1, nt
      ! open modelvar data file
      OPEN (fid, file=TRIM(this%postvarFns(t)), FORM='unformatted', STATUS='old', ACCESS='stream', ACTION='read')

      iprec = 1
      !read tmp2d
      DO skip = 1, 17 
        READ (fid, pos=iprec) ((this%tmp3d(i, j, 1, t), i=1, idn), j=1, jdn)
        iprec = iprec + idn * jdn * 4 
        PRINT *, "iprec:", iprec
        PRINT *, 'min and max of tmp3d:', MINVAL(this%tmp3d(:, :, 1, t)), MAXVAL(this%tmp3d(:, :, 1, t))
        PRINT *, 'done tmp3d reading ...'
      END DO

      ! read uwnd
      READ (fid, pos=iprec) (((this%uwnd(i, j, k, t), i=1, idn), j=1, jdn), k=1, kdn)
      iprec = iprec + idn * jdn * kdn * 4
      PRINT *, "***********   **************", kdn
      PRINT *, "iprec:", iprec
      PRINT *, 'min and max of uwnd:', MINVAL(this%uwnd(:, :, :, t)), MAXVAL(this%uwnd(:, :, :, t))
      PRINT *, 'done uwnd reading ...'

      ! ! read vwnd
      READ (fid, pos=iprec) (((this%vwnd(i, j, k, t), i=1, idn), j=1, jdn), k=1, kdn)
      iprec = iprec + idn * jdn * kdn * 4
      PRINT *, "***********   **************", kdn
      PRINT *, "iprec:", iprec
      PRINT *, 'min and max of vwnd:', MINVAL(this%vwnd(:, :, :, t)), MAXVAL(this%vwnd(:, :, :, t))
      PRINT *, 'done vwnd reading ...'

      ! ! read temp
      READ (fid, pos=iprec) (((this%temp(i, j, k, t), i=1, idn), j=1, jdn), k=1, kdn)
      iprec = iprec + idn * jdn * kdn * 4
      PRINT *, "***********   **************", kdn
      PRINT *, "iprec:", iprec
      PRINT *, 'min and max of temp:', MINVAL(this%temp(:, :, :, t)), MAXVAL(this%temp(:, :, :, t))
      PRINT *, 'done temp reading ...'

      BLOCK
        REAL(r_kind) :: pres_list(17) = (/1000.00, 925.00, 850.00, 700.00, 600.00, 500.00, 400.00, 300.00, 250.00, 200.00, 150.00, 100.00, 70.00, 50.00, 30.00, 20.00, 10.00/)

        DO k = 1, kdn
          this%pres(:, :, k, t) = pres_list(k) * 100.00D0
        END DO
      END BLOCK

      ! ! read height
      READ (fid, pos=iprec) (((this%H(i, j, k, t), i=1, idn), j=1, jdn), k=1, kdn)
      iprec = iprec + idn * jdn * kdn * 4
      PRINT *, "***********   **************", kdn
      PRINT *, "iprec:", iprec
      PRINT *, 'min and max of H:', MINVAL(this%H(:, :, :, t)), MAXVAL(this%H(:, :, :, t))
      PRINT *, 'done height reading ...'

      ! ! read qvapor
      READ (fid, pos=iprec) (((this%qvapor(i, j, k, t), i=1, idn), j=1, jdn), k=1, kdn)
      iprec = iprec + idn * jdn * kdn * 4
      PRINT *, "***********   **************", kdn
      PRINT *, "iprec:", iprec
      PRINT *, 'min and max of qvapor:', MINVAL(this%qvapor(:, :, :, t)), MAXVAL(this%qvapor(:, :, :, t))
      PRINT *, 'done qvapor reading ...'
      ! close modelvar data file
      CLOSE (fid)
      PRINT *, 'done reading modelvar files:', t, '/', nt
    END DO
    ! STOP
  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(GrapesPostvarIO_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%tmp3d)) DEALLOCATE (this%tmp3d)
    IF (ALLOCATED(this%uwnd)) DEALLOCATE (this%uwnd)
    IF (ALLOCATED(this%vwnd)) DEALLOCATE (this%vwnd)
    IF (ALLOCATED(this%temp)) DEALLOCATE (this%temp)
    IF (ALLOCATED(this%pres)) DEALLOCATE (this%pres)
    IF (ALLOCATED(this%pres)) DEALLOCATE (this%qvapor)
    IF (ALLOCATED(this%H)) DEALLOCATE (this%H)

    PRINT *, 'destructor works'

  END SUBROUTINE destructor

END MODULE GrapesPostvarIO_m
