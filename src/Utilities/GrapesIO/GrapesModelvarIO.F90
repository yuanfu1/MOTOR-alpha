!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 1.0
! HISTORY           :
!   Created by Jiongming Pang, 2023/6/10, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

MODULE GrapesModelvarIO_m
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE parameters_m, ONLY: spec_heat_const_pres, dry_air_gas_const, surface_ref_pres, degree2radian
  USE Interp1D_m, ONLY: interp1d_3D_idx2_single, interp1d_3D_idx2
  USE AdvanceTime_m

  IMPLICIT NONE

  TYPE GrapesModelvarIO_t
    CHARACTER(LEN=1024), ALLOCATABLE  ::  modvarFns(:)
    CHARACTER(LEN=1024)               ::  modvarCtlFn, hfFn
    INTEGER(i_kind)                   ::  num_files
    REAL(r_kind), ALLOCATABLE         ::  lat2D(:, :), lon2D(:, :), lat1D(:), lon1D(:)
    REAL(r_single), ALLOCATABLE       ::  zRHght_u(:, :, :), zs(:, :)
    REAL(r_kind), ALLOCATABLE         ::  zRHght_s(:, :, :), topo(:, :)
    REAL(r_single), DIMENSION(:, :, :, :), ALLOCATABLE :: pip, &
                                                          pi, &
                                                          thp, &
                                                          th, &
                                                          u, &
                                                          v, &
                                                          w, &
                                                          qv
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE   :: temp, &
                                                          pres, &
                                                          uwnd, &
                                                          vwnd, &
                                                          qvapor
    REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE      :: t2m, &
                                                          u10m, &
                                                          v10m, &
                                                          psfc, &
                                                          q2m
    INTEGER(i_kind)              :: idn, jdn, kdn, tdn
    INTEGER(i_kind)              :: ids, ide, jds, jde, kds, kde
    INTEGER(i_kind), ALLOCATABLE :: time_unix(:), time_gmt(:, :)

  CONTAINS
    FINAL :: destructor
    PROCEDURE, PRIVATE :: get_domain
    PROCEDURE, PRIVATE :: read_value
    PROCEDURE, PRIVATE :: u_c2a
    PROCEDURE, PRIVATE :: v_c2a
    PROCEDURE, PRIVATE :: compute_t_p
    PROCEDURE, PRIVATE :: compute_zHght_s
  END TYPE GrapesModelvarIO_t

  INTERFACE GrapesModelvarIO_t
    PROCEDURE :: constructor
  END INTERFACE

CONTAINS
  FUNCTION constructor(modvarFns, modvarCtlFn, hfFn, mdate) RESULT(this)
    TYPE(GrapesModelvarIO_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: modvarFns(:), modvarCtlFn, hfFn
    INTEGER(i_kind), INTENT(IN)  :: mdate(:, :)
    INTEGER(i_kind)              :: count_start, count_end, count_rate
    REAL(r_kind)                 :: ts

    this%modvarFns = modvarFns
    this%modvarCtlFn = TRIM(modvarCtlFn)
    this%hfFn = TRIM(hfFn)
    this%time_gmt = mdate

    ! get domain info from ctl files
    CALL this%get_domain()

    ! read vars' value from modelvar and HeightField files
    CALL system_clock(count_start, count_rate)
    CALL this%read_value()
    CALL system_clock(count_end)
    ts = REAL(count_end - count_start) / REAL(count_rate)
    PRINT *, 'Time spent in reading modelvar and HeighField files: ', ts

    ! compute zHght_u
    CALL this%compute_zHght_s()

    ! compute temp and press
    CALL this%compute_t_p()

    ! convert uwnd and vwnd from c-grid to a-grid
    CALL this%u_c2a()
    CALL this%v_c2a()

  END FUNCTION constructor

  SUBROUTINE get_domain(this)
    IMPLICIT NONE
    CLASS(GrapesModelvarIO_t) :: this
    INTEGER                   :: fid
    INTEGER                   :: i
    REAL(r_kind)              :: lon_s, lon_int, lat_s, lat_int
    CHARACTER(1024)           :: tmp

    this%num_files = UBOUND(this%modvarFns, 1)
    this%tdn = this%num_files

    ALLOCATE (this%time_unix(this%tdn))
    DO i = 1, this%tdn
      PRINT *, "GMT Time: ", this%time_gmt(:, i)
      CALL Time_GMT_to_Unix(this%time_gmt(:, i), this%time_unix(i))
      PRINT *, "Unix Time: ", this%time_unix(i)
    END DO

    fid = 120
    OPEN (fid, FILE=TRIM(this%modvarCtlFn), ACTION='read')

    DO i = 1, 5
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

    ! generate lon and lat, 2D-array, and 1D-array
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

    this%lon2D = this%lon2D * degree2radian
    this%lat2D = this%lat2D * degree2radian
    this%lon1D = this%lon1D * degree2radian
    this%lat1D = this%lat1D * degree2radian

    ! get the amount of z-level
    READ (fid, *) tmp, this%kdn
    PRINT *, 'kdn:', this%kdn

    CLOSE (fid)

    this%ids = 1
    this%ide = this%idn
    this%jds = 1
    this%jde = this%jdn
    this%kds = 1
    this%kde = this%kdn
    PRINT *, 'ids, ide:', this%ids, this%ide
    PRINT *, 'jds, jde:', this%jds, this%jde
    PRINT *, 'kds, kde:', this%kds, this%kde
    PRINT *, 'idn, jdn, kdn, tdn:', this%idn, this%jdn, this%kdn, this%tdn

  END SUBROUTINE get_domain

  SUBROUTINE read_value(this)
    IMPLICIT NONE
    CLASS(GrapesModelvarIO_t) :: this
    INTEGER                   :: fid
    INTEGER(i_kind)           :: ids, ide, jds, jde, kds, kde, idn, jdn, kdn, nt
    INTEGER(i_kind)           :: i, j, k, t, iprec

    ids = this%ids
    ide = this%ide
    jds = this%jds
    jde = this%jde
    kds = this%kds
    kde = this%kde

    idn = this%idn
    jdn = this%jdn
    kdn = this%kdn
    nt = this%tdn

    ALLOCATE (this%pip(idn, jdn, kdn, nt))
    ALLOCATE (this%pi(idn, jdn, kdn, nt))
    ALLOCATE (this%thp(idn, jdn, kdn - 1, nt))
    ALLOCATE (this%th(idn, jdn, kdn - 1, nt))
    ALLOCATE (this%u(idn, jdn, kdn - 2, nt))
    ALLOCATE (this%v(idn, jdn, kdn - 2, nt))
    ALLOCATE (this%w(idn, jdn, kdn - 1, nt))
    ALLOCATE (this%qv(idn, jdn, kdn - 1, nt))

    fid = 110   ! file id

    CALL omp_set_num_threads(nt)

    !$OMP PARALLEL DO PRIVATE(t, fid, iprec, k, i, j) SHARED(nt, idn, jdn, kdn)
    DO t = 1, nt
      ! open modelvar data file
      OPEN (newunit=fid, file=TRIM(this%modvarFns(t)), FORM='unformatted', STATUS='old', ACCESS='stream', ACTION='read', CONVERT='LITTLE_ENDIAN')

      iprec = 1

      ! read pip
      DO k = 1, kdn
        READ (fid, pos=iprec) ((this%pip(i, j, k, t), i=1, idn), j=1, jdn)
        iprec = iprec + idn * jdn * 4
        ! PRINT *, 'min and max of pip in k-level:',k, MINVAL(this%pip(:,:,k,t)), MAXVAL(this%pip(:,:,k,t))
      END DO
      PRINT *, 'min and max of pip at tid:', t, MINVAL(this%pip(:,:,:,t)), MAXVAL(this%pip(:,:,:,t))
      PRINT *, 'done pip reading ...'

      ! read pi
      DO k = 1, kdn
        READ (fid, pos=iprec) ((this%pi(i, j, k, t), i=1, idn), j=1, jdn)
        iprec = iprec + idn * jdn * 4
      END DO
      PRINT *, 'min and max of pi at tid:', t, MINVAL(this%pi(:,:,:,t)), MAXVAL(this%pi(:,:,:,t))
      PRINT *, 'done pi reading ...'

      ! read thp
      DO k = 1, kdn - 1
        READ (fid, pos=iprec) ((this%thp(i, j, k, t), i=1, idn), j=1, jdn)
        iprec = iprec + idn * jdn * 4
      END DO
      PRINT *, 'min and max of thp at tid:', t, MINVAL(this%thp(:,:,:,t)), MAXVAL(this%thp(:,:,:,t))
      PRINT *, 'done thp reading ...'

      ! read th
      DO k = 1, kdn - 1
        READ (fid, pos=iprec) ((this%th(i, j, k, t), i=1, idn), j=1, jdn)
        iprec = iprec + idn * jdn * 4
      END DO
      PRINT *, 'min and max of th at tid:', t, MINVAL(this%th(:,:,:,t)), MAXVAL(this%th(:,:,:,t))
      PRINT *, 'done th reading ...'

      ! read u
      DO k = 1, kdn - 2
        READ (fid, pos=iprec) ((this%u(i, j, k, t), i=1, idn), j=1, jdn)
        iprec = iprec + idn * jdn * 4
      END DO
      PRINT *, 'min and max of u at tid:', t, MINVAL(this%u(:,:,:,t)), MAXVAL(this%u(:,:,:,t))
      PRINT *, 'done u reading ...'

      ! read v
      DO k = 1, kdn - 2
        READ (fid, pos=iprec) ((this%v(i, j, k, t), i=1, idn), j=1, jdn)
        iprec = iprec + idn * jdn * 4
      END DO
      PRINT *, 'min and max of v at tid:', t, MINVAL(this%v(:,:,:,t)), MAXVAL(this%v(:,:,:,t))
      PRINT *, 'done v reading ...'

      ! read w
      DO k = 1, kdn - 1
        READ (fid, pos=iprec) ((this%w(i, j, k, t), i=1, idn), j=1, jdn)
        iprec = iprec + idn * jdn * 4
      END DO
      PRINT *, 'min and max of w at tid:', t, MINVAL(this%w(:,:,:,t)), MAXVAL(this%w(:,:,:,t))
      PRINT *, 'done w reading ...'

      ! read qv-moist
      DO k = 1, kdn - 1
        READ (fid, pos=iprec) ((this%qv(i, j, k, t), i=1, idn), j=1, jdn)
        iprec = iprec + idn * jdn * 4
      END DO
      PRINT *, 'min and max of qv at tid:', t, MINVAL(this%qv(:,:,:,t)), MAXVAL(this%qv(:,:,:,t))
      this%qvapor = this%qv
      PRINT *, 'done qvapor reading ...'

      ! close modelvar data file
      CLOSE (fid)
      PRINT *, 'done reading modelvar files:', t, '/', nt
    END DO
    !$OMP END PARALLEL DO

    ! open HeightField data file
    OPEN (fid, file=TRIM(this%hfFn), FORM='unformatted', ACTION='read', STATUS='old', ACCESS='stream', CONVERT='BIG_ENDIAN')

    iprec = 1

    ! read zz, the real-height of domain in each z-level
    ALLOCATE (this%zRHght_u(ids:ide, jds:jde, kds:kde))
    DO k = kds, kde
      READ (fid, pos=iprec) ((this%zRHght_u(i, j, k), i=ids, ide), j=jds, jde)
      iprec = iprec + idn * jdn * 4 + 2 * 4
    END DO
    PRINT *, 'max and min of zRHght_u:', MAXVAL(this%zRHght_u), MINVAL(this%zRHght_u)
    PRINT *, 'done zz reading ...'

    ! assign topo from zz
    ALLOCATE (this%topo(ids:ide, jds:jde))
    this%topo = this%zRHght_u(:, :, kds + 1)
    PRINT *, 'done topo assigning ...'

    ! close HeightField data file
    CLOSE (fid)

  END SUBROUTINE

  SUBROUTINE compute_zHght_s(this)
    IMPLICIT NONE

    CLASS(GrapesModelvarIO_t) :: this
    INTEGER(i_kind)           :: ids, ide, jds, jde, kds, kde
    INTEGER(i_kind)           :: i, j, k

    ids = this%ids
    ide = this%ide
    jds = this%jds
    jde = this%jde
    kds = this%kds
    kde = this%kde

    ALLOCATE (this%zRHght_s(ids:ide, jds:jde, kds:kde - 1))
    this%zRHght_s = (this%zRHght_u(:, :, kds:kde - 1) + this%zRHght_u(:, :, kds + 1:kde)) / 2.0D0
    PRINT *, 'done zRHght_s computation ...'

  END SUBROUTINE compute_zHght_s

  SUBROUTINE compute_t_p(this)
    IMPLICIT NONE
    CLASS(GrapesModelvarIO_t) :: this
    INTEGER(i_kind)           :: ids, ide, jds, jde, kds, kde, tdn
    INTEGER(i_kind)           :: i, j, k
    REAL(r_kind)              :: r_d, cp

    ids = this%ids
    ide = this%ide
    jds = this%jds
    jde = this%jde
    kds = this%kds
    kde = this%kde
    tdn = this%tdn

    r_d = 287.0D0
    cp = 7.0D0 * r_d / 2.0D0

    ALLOCATE (this%temp(ids:ide, jds:jde, kds:kde - 1, 1:tdn), &
              this%pres(ids:ide, jds:jde, kds:kde - 1, 1:tdn))

    this%temp = this%th(:, :, kds:kde - 1, :) * ((this%pi(:, :, kds + 1:kde, :) + this%pi(:, :, kds:kde - 1, :)) * 0.5D0)
    this%pres = EXP(LOG(1000.0D0) - cp * LOG(this%th(:, :, kds:kde - 1, :) / this%temp(:, :, kds:kde - 1, :)) / r_d) * 100.0D0

    PRINT *, 'min and max of temp:', MINVAL(this%temp), MAXVAL(this%temp)
    PRINT *, 'min and max of pres:', MINVAL(this%pres), MAXVAL(this%pres)

    PRINT *, 'done temp and pres computation ...'

  END SUBROUTINE compute_t_p

  SUBROUTINE u_c2a(this)
    CLASS(GrapesModelvarIO_t) :: this
    INTEGER(i_kind)           :: ids, ide, jds, jde, kds, kde, tdn
    INTEGER(i_kind)           :: i, j, k
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: tmp

    ids = this%ids
    ide = this%ide
    jds = this%jds
    jde = this%jde
    kds = this%kds
    kde = this%kde
    tdn = this%tdn

    ALLOCATE (tmp(ids:ide, jds:jde, kds:kde - 1, 1:tdn))
    tmp(:, :, kds, :) = this%u(:, :, 1, :)
    tmp(:, :, kde - 1, :) = this%u(:, :, kde - 2, :)
    tmp(:, :, 2:kde - 2, :) = (this%u(:, :, 2:kde - 2, :) + this%u(:, :, 1:kde - 3, :)) / 2.0

    ALLOCATE (this%uwnd(ids:ide, jds:jde, kds:kde - 1, 1:tdn))
    this%uwnd(1, :, :, :) = tmp(1, :, :, :)
    this%uwnd(2:ide, :, :, :) = (tmp(1:ide - 1, :, :, :) + tmp(2:ide, :, :, :)) / 2.0

    DEALLOCATE (tmp)
    PRINT *, 'done u_c2a ...'
  END SUBROUTINE u_c2a

  SUBROUTINE v_c2a(this)
    CLASS(GrapesModelvarIO_t) :: this
    INTEGER(i_kind)           :: ids, ide, jds, jde, kds, kde, tdn
    INTEGER(i_kind)           :: i, j, k
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: tmp

    ids = this%ids
    ide = this%ide
    jds = this%jds
    jde = this%jde
    kds = this%kds
    kde = this%kde
    tdn = this%tdn

    ALLOCATE (tmp(ids:ide, jds:jde, kds:kde - 1, 1:tdn))
    tmp(:, :, kds, :) = this%v(:, :, 1, :)
    tmp(:, :, kde - 1, :) = this%v(:, :, kde - 2, :)
    tmp(:, :, 2:kde - 2, :) = (this%v(:, :, 2:kde - 2, :) + this%v(:, :, 1:kde - 3, :)) / 2.0

    ALLOCATE (this%vwnd(ids:ide, jds:jde, kds:kde - 1, 1:tdn))
    this%vwnd(:, 1, :, :) = tmp(:, 1, :, :)
    this%vwnd(:, 2:jde, :, :) = (tmp(:, 1:jde - 1, :, :) + tmp(:, 2:jde, :, :)) / 2.0

    DEALLOCATE (tmp)
    PRINT *, 'done v_c2a ...'
  END SUBROUTINE v_c2a

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(GrapesModelvarIO_t), INTENT(INOUT) :: this

    !DEALLOCATE (this%pip, this%pi, this%thp, this%th, this%u, this%v, this%w, this%qv)

    IF (ALLOCATED(this%pip)) DEALLOCATE(this%pip)
    IF (ALLOCATED(this%pi)) DEALLOCATE(this%pi)
    IF (ALLOCATED(this%thp)) DEALLOCATE(this%thp)
    IF (ALLOCATED(this%th)) DEALLOCATE(this%th)
    IF (ALLOCATED(this%u)) DEALLOCATE(this%u)
    IF (ALLOCATED(this%v)) DEALLOCATE(this%v)
    IF (ALLOCATED(this%w)) DEALLOCATE(this%w)
    IF (ALLOCATED(this%qv)) DEALLOCATE(this%qv)


    PRINT *, 'destructor works'

  END SUBROUTINE destructor

END MODULE GrapesModelvarIO_m
