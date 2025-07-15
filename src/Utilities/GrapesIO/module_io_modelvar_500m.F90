!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 1.0
! HISTORY           :
!   Created by Jiongming Pang, 2023/6/10, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

MODULE module_io_modelvar_500m_m
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE parameters_m, ONLY: spec_heat_const_pres, dry_air_gas_const, surface_ref_pres, degree2radian
  USE Interp1D_m, ONLY: interp1d_3D_idx2_single, interp1d_3D_idx2
  USE AdvanceTime_m

  IMPLICIT NONE

  TYPE module_io_modelvar_500m_t
    CHARACTER(LEN=1024), ALLOCATABLE  ::  modvarFns(:)
    CHARACTER(LEN=1024)               ::  modvarCtlFn
    INTEGER(i_kind)                   ::  num_files, vLevel
    REAL(r_kind), ALLOCATABLE         ::  lat2D(:, :), lon2D(:, :), lat1D(:), lon1D(:)
    REAL(r_kind), ALLOCATABLE       ::  zRHght_u(:, :, :), zs(:, :)
    REAL(r_kind), ALLOCATABLE         ::  zRHght_s(:, :, :), topo(:, :)
    REAL(r_single), DIMENSION(:, :, :, :), ALLOCATABLE :: pi, &
                                                          th, &
                                                          u, &
                                                          v, &
                                                          w, &
                                                          qv, &
                                                          zz
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE   :: temp, &
                                                          pres, &
                                                          uwnd, &
                                                          vwnd, &
                                                          wwnd, &
                                                          qvapor, &
                                                          height
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE    :: rainc_s, &
                                                          rainnc_s, &
                                                          t2m_s, &
                                                          q2m_s, &
                                                          u10m_s, &
                                                          v10m_s, &
                                                          vis_s, &
                                                          psl_s, &
                                                          psfc_s
    REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE      :: rainc, &
                                                          rainnc, &
                                                          t2m, &
                                                          q2m, &
                                                          u10m, &
                                                          v10m, &
                                                          vis, &
                                                          psl, &
                                                          psfc
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
    PROCEDURE, PUBLIC :: writeBackToGrapesInput
  END TYPE module_io_modelvar_500m_t

  INTERFACE module_io_modelvar_500m_t
    PROCEDURE :: constructor
  END INTERFACE

CONTAINS
  SUBROUTINE init_modelvar(modvarFns, modvarCtlFn, mdate, vLevel)
    TYPE(module_io_modelvar_500m_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: modvarFns(:), modvarCtlFn
    INTEGER(i_kind), INTENT(IN)  :: mdate(:, :)
    INTEGER(i_kind), INTENT(IN) :: vLevel

    this%modvarFns = modvarFns
    this%modvarCtlFn = TRIM(modvarCtlFn)
    this%time_gmt = mdate
    this%vLevel = vLevel
    PRINT *, "GMT Time: ", this%time_gmt(:, 1)
    ! get domain info from ctl files
    CALL this%get_domain()

    ! read vars' value from modelvar and HeightField files
    CALL this%read_value()
    IF (this%vLevel == 1) THEN
      PRINT *, "Read modelvar surface data"
    ELSE
      PRINT *, "Read modelvar level data"
      CALL this%compute_zHght_s()
      CALL this%compute_t_p()
      CALL this%u_c2a()
      CALL this%v_c2a()
    END IF
   

  END SUBROUTINE init_modelvar

  SUBROUTINE get_domain(this)
    IMPLICIT NONE
    CLASS(module_io_modelvar_500m_t) :: this
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
    CLASS(module_io_modelvar_500m_t) :: this
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

    ALLOCATE (this%pi(idn, jdn, kdn, nt))
    ALLOCATE (this%th(idn, jdn, kdn - 1, nt))
    ALLOCATE (this%u(idn, jdn, kdn - 2, nt))
    ALLOCATE (this%v(idn, jdn, kdn - 2, nt))
    ALLOCATE (this%w(idn, jdn, kdn - 1, nt))
    ALLOCATE (this%qv(idn, jdn, kdn - 1, nt))
    ALLOCATE (this%zz(idn, jdn, kdn - 1, nt))
    ALLOCATE (this%rainc_s(idn, jdn, nt))
    ALLOCATE (this%rainnc_s(idn, jdn, nt))
    ALLOCATE (this%t2m_s(idn, jdn, nt))
    ALLOCATE (this%q2m_s(idn, jdn, nt))
    ALLOCATE (this%u10m_s(idn, jdn, nt))
    ALLOCATE (this%v10m_s(idn, jdn, nt))
    ALLOCATE (this%vis_s(idn, jdn, nt))
    ALLOCATE (this%psl_s(idn, jdn, nt))

    fid = 110   ! file id
    this%zz = 0.0D0

    IF (this%vLevel == 1) THEN
      DO t = 1, nt

        ! open modelvar data file
        OPEN (fid, file=TRIM(this%modvarFns(t)), FORM='unformatted', STATUS='old', ACCESS='stream', ACTION='read', CONVERT='LITTLE_ENDIAN')
        iprec = 1 + idn * jdn * 4 * 67 + idn * jdn * 4 * 66 * 4 + idn * jdn * 4 * 65 * 2

        ! read surface vars
        !read rainc
        READ (fid, pos=iprec) this%rainc_s(:,:,t)
        this%rainc = this%rainc_s
        iprec = iprec + idn * jdn * 4
        PRINT *, 'min and max of rainc:', MINVAL(this%rainc), MAXVAL(this%rainc)
        PRINT *, 'done rainc reading ...'

        ! read rainnc
        READ (fid, pos=iprec) this%rainnc_s(:,:,t)
        this%rainnc = this%rainnc_s
        iprec = iprec + idn * jdn * 4
        PRINT *, 'min and max of rainnc:', MINVAL(this%rainnc), MAXVAL(this%rainnc)
        PRINT *, 'done rainnc reading ...'

        !read t2m
        READ (fid, pos=iprec) this%t2m_s(:,:,t)
        this%t2m = this%t2m_s
        iprec = iprec + idn * jdn * 4
        PRINT *, 'min and max of t2m:', MINVAL(this%t2m), MAXVAL(this%t2m)
        PRINT *, 'done t2m reading ...'

        ! read q2m
        READ (fid, pos=iprec) this%q2m_s(:, :, t)
        this%q2m = this%q2m_s
        iprec = iprec + idn * jdn * 4
        PRINT *, 'min and max of q2m:', MINVAL(this%q2m), MAXVAL(this%q2m)
        PRINT *, 'done q2m reading ...'


        ! read u10m
        READ (fid, pos=iprec) this%u10m_s(:,:,t)
        this%u10m = this%u10m_s
        iprec = iprec + idn * jdn * 4
        PRINT *, 'min and max of u10m:', MINVAL(this%u10m), MAXVAL(this%u10m)
        PRINT *, 'done u10m reading ...'

        ! read v10m
        READ (fid, pos=iprec) this%v10m_s(:,:,t)
        this%v10m = this%v10m_s
        iprec = iprec + idn * jdn * 4
        PRINT *, 'min and max of v10m:', MINVAL(this%v10m), MAXVAL(this%v10m)
        PRINT *, 'done v10m reading ...'

        ! read vis
        READ (fid, pos=iprec) this%vis_s(:,:,t)
        this%vis = this%vis_s 
        iprec = iprec + idn * jdn * 4
        PRINT *, 'min and max of vis:', MINVAL(this%vis), MAXVAL(this%vis)
        PRINT *, 'done vis reading ...'

        ! read psl
        READ (fid, pos=iprec) this%psl_s(:, :, t)
        this%psl = this%psl_s
        iprec = iprec + idn * jdn * 4
        PRINT *, 'min and max of psl:', MINVAL(this%psl), MAXVAL(this%psl)
        PRINT *, 'done psl reading ...'
      END DO
      CLOSE (fid)
      PRINT *, 'done reading modelvar files:', t, '/', nt
    ELSE 
      DO t = 1, nt
        ! open modelvar data file
        OPEN (fid, file=TRIM(this%modvarFns(t)), FORM='unformatted', STATUS='old', ACCESS='stream', ACTION='read', CONVERT='LITTLE_ENDIAN')
        iprec = 1

        ! read pi
        DO k = 1, kdn
          READ (fid, pos=iprec) ((this%pi(i, j, k, t), i=1, idn), j=1, jdn)
          iprec = iprec + idn * jdn * 4
        END DO
        
        PRINT *, 'min and max of pi:', MINVAL(this%pi), MAXVAL(this%pi)
        PRINT *, 'done pi reading ...'

        ! read th
        DO k = 1, kdn - 1
          READ (fid, pos=iprec) ((this%th(i, j, k, t), i=1, idn), j=1, jdn)
          iprec = iprec + idn * jdn * 4
        END DO
        PRINT *, 'min and max of th:', MINVAL(this%th), MAXVAL(this%th)
        PRINT *, 'done th reading ...'

        ! read u
        DO k = 1, kdn - 2
          READ (fid, pos=iprec) ((this%u(i, j, k, t), i=1, idn), j=1, jdn)
          iprec = iprec + idn * jdn * 4
        END DO
        PRINT *, 'min and max of u:', MINVAL(this%u), MAXVAL(this%u)
        PRINT *, 'done u reading ...'

        ! read v
        DO k = 1, kdn - 2
          READ (fid, pos=iprec) ((this%v(i, j, k, t), i=1, idn), j=1, jdn)
          iprec = iprec + idn * jdn * 4
        END DO
        PRINT *, 'min and max of v:', MINVAL(this%v), MAXVAL(this%v)
        PRINT *, 'done v reading ...'

        ! read w
        DO k = 1, kdn - 1
          READ (fid, pos=iprec) ((this%w(i, j, k, t), i=1, idn), j=1, jdn)
          iprec = iprec + idn * jdn * 4
        END DO
        PRINT *, 'min and max of w:', MINVAL(this%w), MAXVAL(this%w)
        this%wwnd = this%w
        PRINT *, 'done w reading ...'

        ! read qv-moist
        DO k = 1, kdn - 1
          READ (fid, pos=iprec) ((this%qv(i, j, k, t), i=1, idn), j=1, jdn)
          iprec = iprec + idn * jdn * 4
        END DO
        PRINT *, 'min and max of qv:', MINVAL(this%qv), MAXVAL(this%qv)
        this%qvapor = this%qv
        PRINT *, 'done qvapor reading ...'

        ! read height
        READ (fid, pos=iprec) this%zz(:,:,:,t)
        this%height = this%zz
        iprec = iprec + idn * jdn * (kdn - 1) * 4
        PRINT *, 'min and max of height:', MINVAL(this%height), MAXVAL(this%height)
        PRINT *, 'done height reading ...'
        ! close modelvar data file
      END DO

      ! assign topo from zz
      this%zRHght_u = this%height(:, :, :, 1)
      PRINT *, minval(this%zRHght_u), maxval(this%zRHght_u)
      PRINT *, 'done zRHght_u assigning ...'

      ALLOCATE (this%topo(ids:ide, jds:jde))
      this%topo = this%zRHght_u(:, :, kds)
      PRINT *, 'done topo assigning ...'

      CLOSE (fid)
    END IF

  END SUBROUTINE

  SUBROUTINE compute_zHght_s(this)
    IMPLICIT NONE

    CLASS(module_io_modelvar_500m_t) :: this
    INTEGER(i_kind)           :: ids, ide, jds, jde, kds, kde
    INTEGER(i_kind)           :: i, j, k

    ids = this%ids
    ide = this%ide
    jds = this%jds
    jde = this%jde
    kds = this%kds
    kde = this%kde

    ALLOCATE (this%zRHght_s(ids:ide, jds:jde, kds:kde))
    
    this%zRHght_s(: , :, kds+1:kde-1) = (this%zRHght_u(:, :, kds:kde - 2) + this%zRHght_u(:, :, kds + 1:kde-1)) / 2.0D0
    this%zRHght_s(: , :, kds) = this%zRHght_u(:, :, kds) - (this%zRHght_s(: , :, kds + 1) - this%zRHght_u(:, :, kds))
    this%zRHght_s(: , :, kde) = this%zRHght_u(:, :, kde-1) + (this%zRHght_u(:, :, kde-1) - this%zRHght_s(: , :, kde-1))

    PRINT*, this%zRHght_u(1,1,:)
    PRINT*, this%zRHght_s(1,1,:)
    PRINT *, 'done zRHght_s computation ...'

  END SUBROUTINE compute_zHght_s

  SUBROUTINE compute_t_p(this) 
    USE parameters_m, ONLY: spec_heat_const_pres, dry_air_gas_const, surface_ref_pres
        IMPLICIT NONE

    CLASS(module_io_modelvar_500m_t) :: this
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
    this%pres = surface_ref_pres * 100.0D0 * ((this%pi(:, :, kds + 1:kde, :) + this%pi(:, :, kds:kde - 1, :)) * 0.5D0) ** (spec_heat_const_pres / dry_air_gas_const)
    PRINT *, 'min and max of temp:', MINVAL(this%temp), MAXVAL(this%temp)
    PRINT *, 'min and max of pres:', MINVAL(this%pres), MAXVAL(this%pres)

    PRINT *, 'done temp and pres computation ...'
    ! STOP

  END SUBROUTINE compute_t_p

  SUBROUTINE u_c2a(this)
    CLASS(module_io_modelvar_500m_t) :: this
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

    ALLOCATE (tmp(ids:ide, jds:jde, kds:kde-1, 1:tdn))
    ALLOCATE (this%uwnd(ids:ide, jds:jde, kds:kde - 1, 1:tdn))
    tmp(:, :, kds+1:kde-2, :) = (this%u(:, :, kds:kde-3, :) + this%u(:, :, kds+1:kde-2, :)) / 2.0D0
    tmp(:, :, kds, :) = this%u(:, :, kds, :) - (tmp(:, :, kds + 1, :) - this%u(:, :, kds, :))
    tmp(:, :, kde-1, :) = this%u(:, :, kde-2, :) + (this%u(:, :, kde-2, :) - tmp(:, :, kde-2, :))

    this%uwnd(2:ide-1, :, :, :) = (tmp(1:ide-2, :, :, :) + tmp(2:ide-1, :, :, :)) / 2.0
    this%uwnd(1, :, :, :) = tmp(1, :, :, :) - (this%uwnd(2, :, :, :) - tmp(1, :, :, :))
    this%uwnd(ide, :, :, :) = tmp(ide-1, :, :, :) + (tmp(ide-1, :, :, :) - this%uwnd(ide-1, :, :, :))


    ! tmp(:, :, kds, :) = this%u(:, :, 1, :)
    ! tmp(:, :, kde - 1, :) = this%u(:, :, kde - 2, :)
    ! tmp(:, :, 2:kde - 2, :) = (this%u(:, :, 2:kde - 2, :) + this%u(:, :, 1:kde - 3, :)) / 2.0

    ! ALLOCATE (this%uwnd(ids:ide, jds:jde, kds:kde - 1, 1:tdn))
    ! this%uwnd(1, :, :, :) = tmp(1, :, :, :)
    ! this%uwnd(2:ide, :, :, :) = (tmp(1:ide - 1, :, :, :) + tmp(2:ide, :, :, :)) / 2.0

    DEALLOCATE (tmp)
    PRINT *, 'done u_c2a ...'
  END SUBROUTINE u_c2a

  SUBROUTINE v_c2a(this)
    CLASS(module_io_modelvar_500m_t) :: this
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

    ALLOCATE (tmp(ids:ide, jds:jde, kds:kde-1, 1:tdn))
    ALLOCATE (this%vwnd(ids:ide, jds:jde, kds:kde - 1, 1:tdn))
    tmp(:, :, kds+1:kde-2, :) = (this%v(:, :, kds:kde-3, :) + this%v(:, :, kds+1:kde-2, :)) / 2.0D0
    tmp(:, :, kds, :) = this%v(:, :, kds, :) - (tmp(:, :, kds + 1, :) - this%v(:, :, kds, :))
    tmp(:, :, kde-1, :) = this%v(:, :, kde-2, :) + (this%v(:, :, kde-2, :) - tmp(:, :, kde-2, :))

    this%vwnd(:, 2:jde-1, :, :) = (tmp(:, 1:jde-2, :, :) + tmp(:, 2:jde-1, :, :)) / 2.0
    this%vwnd(:, 1, :, :) = tmp(:, 1, :, :) - (this%vwnd(:, 2, :, :) - tmp(:, 1, :, :))
    this%vwnd(:, jde, :, :) = tmp(:, jde-1, :, :) + (tmp(:, jde-1, :, :) - this%vwnd(:, jde-1, :, :))

    ! ALLOCATE (tmp(ids:ide, jds:jde, kds:kde - 1, 1:tdn))
    ! tmp(:, :, kds, :) = this%v(:, :, 1, :)
    ! tmp(:, :, kde - 1, :) = this%v(:, :, kde - 2, :)
    ! tmp(:, :, 2:kde - 2, :) = (this%v(:, :, 2:kde - 2, :) + this%v(:, :, 1:kde - 3, :)) / 2.0

    ! ALLOCATE (this%vwnd(ids:ide, jds:jde, kds:kde - 1, 1:tdn))
    ! this%vwnd(:, 1, :, :) = tmp(:, 1, :, :)
    ! this%vwnd(:, 2:jde, :, :) = (tmp(:, 1:jde - 1, :, :) + tmp(:, 2:jde, :, :)) / 2.0

    DEALLOCATE (tmp)
    PRINT *, 'done v_c2a ...'
  END SUBROUTINE v_c2a

  SUBROUTINE writeBackToGrapesInput(this)
    ! USE Filter_m, ONLY: guidedfilter

    CLASS(module_io_modelvar_500m_t) :: this
    INTEGER(i_kind) :: i, j, k
    REAL(r_single) :: swap
    REAL(r_single), ALLOCATABLE :: oriSwap(:, :, :)
    REAL(r_kind), ALLOCATABLE :: cc(:, :)
    ALLOCATE (cc(this%idn, this%jdn))

    ! Update the filename
    this%module_io_grapes%giFileName = this%giFileName

    IF (ALLOCATED(this%hgrid%u)) THEN
      CALL this%ua2uc(this%hgrid%u)
      CALL this%module_io_grapes%read_value_3d('u', oriSwap, this%hgrid)
      PRINT *, ALLOCATED(oriSwap), SHAPE(oriSwap)
      this%hgrid%u = oriSwap + this%hgrid%u

      !!!!!!!!
      ! BLOCK
      !   DO i = this%kds, this%kde
      !     cc = this%hgrid%u(:, i, :)*1.0D0
      !     CALL guidedfilter(cc, cc, 1, 0.4D0, cc)
      !     this%hgrid%u(:, i, :) = cc
      !   END DO
      ! END BLOCK
      !!!!!!!!

      CALL this%module_io_grapes%write_value_3d('u', this%hgrid%u, this%hgrid)
      PRINT *, 'u is flushed to the file.'
    END IF

    IF (ALLOCATED(this%hgrid%v)) THEN
      CALL this%va2vc(this%hgrid%v)
      CALL this%module_io_grapes%read_value_3d('v', oriSwap, this%hgrid)
      this%hgrid%v = oriSwap + this%hgrid%v

      !!!!!!!!
      ! BLOCK
      !   DO i = this%kds, this%kde
      !     cc = this%hgrid%v(:, i, :)*1.0D0
      !     CALL guidedfilter(cc, cc, 1, 0.4D0, cc)
      !     this%hgrid%v(:, i, :) = cc
      !   END DO
      ! END BLOCK
      !!!!!!!!

      CALL this%module_io_grapes%write_value_3d('v', this%hgrid%v, this%hgrid)
      PRINT *, 'v is flushed to the file.'
    END IF

    IF (ALLOCATED(this%hgrid%q)) THEN
      CALL this%module_io_grapes%read_value_3d('q', oriSwap, this%hgrid)
      this%hgrid%q = oriSwap + this%hgrid%q
      this%hgrid%q(:, this%kds, :) = 0.0

      WHERE(this%hgrid%q <= 0.0D0) this%hgrid%q = 1E-7  !Add for passing CMA MESO forecast

      !!!!!!!!
      ! BLOCK
      !   DO i = this%kds, this%kde
      !     cc = this%hgrid%q(:, i, :)*1.0D0
      !     CALL guidedfilter(cc, cc, 1, 0.4D0, cc)
      !     this%hgrid%q(:, i, :) = cc
      !   END DO
      ! END BLOCK
      !!!!!!!!

      CALL this%module_io_grapes%write_value_3d('q', this%hgrid%q, this%hgrid)
      PRINT *, 'q is flushed to the file.'
    END IF

    ! IF (ALLOCATED(this%hgrid%ps)) THEN
    !   CALL this%write_value_2d('ps', this%hgrid%ps)
    !   PRINT *, 'ps is flushed to the file.'
    ! END IF

    IF (ALLOCATED(this%hgrid%th)) THEN
      this%hgrid%th(:, this%kds, :) = 0.0
      CALL this%module_io_grapes%read_value_3d('th', oriSwap, this%hgrid)
      this%hgrid%th = oriSwap + this%hgrid%th
      this%hgrid%th(:, this%kds, :) = 0.0

      !!!!!!!!
      ! BLOCK
      !   DO i = this%kds, this%kde
      !     cc = this%hgrid%th(:, i, :)*1.0D0
      !     CALL guidedfilter(cc, cc, 1, 0.4D0, cc)
      !     this%hgrid%th(:, i, :) = cc
      !   END DO
      ! END BLOCK
      !!!!!!!!

      CALL this%module_io_grapes%write_value_3d('th', this%hgrid%th, this%hgrid)
      PRINT *, 'th is flushed to the file.'

      CALL this%module_io_grapes%read_value_3d('thp', this%hgrid%thp, this%hgrid)
      this%hgrid%thp = this%hgrid%th - (oriSwap - this%hgrid%thp)
      this%hgrid%thp(:, this%kds, :) = 0.0
      CALL this%module_io_grapes%write_value_3d('thp', this%hgrid%thp, this%hgrid)
      PRINT *, 'thp is flushed to the file.'
    END IF

    IF (ALLOCATED(this%hgrid%pi)) THEN
      CALL this%module_io_grapes%read_value_3d('pi', oriSwap, this%hgrid)
      this%hgrid%pi = oriSwap + this%hgrid%pi

      !!!!!!!!
      ! BLOCK
      !   DO i = this%kds, this%kde
      !     cc = DLOG(this%hgrid%pi(:, i, :)*1.0D0)
      !     CALL guidedfilter(cc, cc, 1, 0.4D0, cc)
      !     this%hgrid%pi(:, i, :) = EXP(cc)
      !   END DO
      ! END BLOCK
      !!!!!!!!

      ! Use the hydrostatic to get the first and the end level of pi and pip
      this%hgrid%pi(:, this%kds_u, :) = this%hgrid%pi(:, this%kds_u + 1, :) - &
                                        (-1.0D0 * g / spec_heat_const_pres / (this%hgrid%th(:, this%kds_s, :) &
                                                                              * (1 + 0.608 * this%hgrid%q(:, this%kds_s, :)))) &
                                        * (this%zRHght_u(:, this%kds_u + 1, :) - this%zRHght_u(:, this%kds_u, :))

      this%hgrid%pi(:, this%kde_u, :) = this%hgrid%pi(:, this%kde_u - 1, :) + &
                                        (-1.0D0 * g / spec_heat_const_pres / (this%hgrid%th(:, this%kde_s, :) &
                                                                              * (1 + 0.608 * this%hgrid%q(:, this%kde_s, :)))) &
                                        * (this%zRHght_u(:, this%kde_u, :) - this%zRHght_u(:, this%kde_u - 1, :))

      !-------------------------------------------------

      CALL this%module_io_grapes%write_value_3d('pi', this%hgrid%pi, this%hgrid)
      PRINT *, 'pi is flushed to the file.'

      CALL this%module_io_grapes%read_value_3d('pip', this%hgrid%pip, this%hgrid)
      this%hgrid%pip = this%hgrid%pi - (oriSwap - this%hgrid%pip)
      CALL this%module_io_grapes%write_value_3d('pip', this%hgrid%pip, this%hgrid)
      PRINT *, 'pip is flushed to the file.'
    END IF

    IF (ALLOCATED(oriSwap)) DEALLOCATE (oriSwap)
    DEALLOCATE (cc)
  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(module_io_modelvar_500m_t), INTENT(INOUT) :: this
    IF (ALLOCATED(this%pi)) DEALLOCATE(this%pi)
    IF (ALLOCATED(this%th)) DEALLOCATE(this%th)
    IF (ALLOCATED(this%u)) DEALLOCATE(this%u)
    IF (ALLOCATED(this%v)) DEALLOCATE(this%v)
    IF (ALLOCATED(this%w)) DEALLOCATE(this%w)
    IF (ALLOCATED(this%qv)) DEALLOCATE(this%qv)
    IF (ALLOCATED(this%rainc_s)) DEALLOCATE(this%rainc_s)
    IF (ALLOCATED(this%rainnc_s)) DEALLOCATE(this%rainnc_s)
    IF (ALLOCATED(this%t2m_s)) DEALLOCATE(this%t2m_s)
    IF (allocated(this%q2m_s)) DEALLOCATE(this%q2m_s)
    IF (ALLOCATED(this%u10m_s)) DEALLOCATE(this%u10m_s)
    IF (ALLOCATED(this%v10m_s)) DEALLOCATE(this%v10m_s)
    IF (ALLOCATED(this%vis_s)) DEALLOCATE(this%vis_s)
    IF (ALLOCATED(this%psl_s)) DEALLOCATE(this%psl_s)

    PRINT *, 'destructor works'

  END SUBROUTINE destructor

END MODULE module_io_modelvar_500m_m
