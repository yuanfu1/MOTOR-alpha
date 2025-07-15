!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.BMatrix.BField
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2024/1/24, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE BFieldHybInc_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE AuxTypeSG_m, ONLY: AuxTypeSG_t
  USE BFieldBase_m, ONLY: BFieldBase_t
  USE BKErr_m, ONLY: BKErr_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE YAMLRead_m
  USE geoTools_m

  TYPE, EXTENDS(BFieldBase_t):: BFieldHybInc_t
    TYPE(ObsSet_t), POINTER   :: Y
    REAL(r_kind), ALLOCATABLE :: fmatz(:, :, :, :), fmat0z(:, :, :), znorm_new(:, :)
    REAL(r_kind), ALLOCATABLE :: fmatx(:, :, :, :, :), fmat0x(:, :, :, :), xnorm_new(:, :, :)
    REAL(r_kind), ALLOCATABLE :: fmaty(:, :, :, :), fmat0y(:, :, :), ynorm_new(:, :)
    REAL(r_kind)    :: s_ens_hv, s_ens_vv
    REAL(r_kind)    :: s_ens_h_gu_x, s_ens_h_gu_y, s_ens_h, s_ens_v, lozRadius
    REAL(r_kind)    :: beta_s, beta_e
    INTEGER(i_kind) :: nc   ! the number of control variables
    INTEGER(i_kind) :: nz   ! the number of non-zero eigenvalues
    INTEGER(i_kind) :: kl
    LOGICAL         :: s_ens_vlv
    LOGICAL         :: init_flag
  CONTAINS
    PROCEDURE, PRIVATE :: loadBMatFiles
    PROCEDURE, PRIVATE :: convert_km_to_grid_units
    PROCEDURE, PRIVATE :: init_rf_z
    PROCEDURE, PRIVATE :: init_rf_x
    PROCEDURE, PRIVATE :: init_rf_y
    PROCEDURE, PRIVATE :: normal_new_factorization_rf_z
    PROCEDURE, PRIVATE :: normal_new_factorization_rf_x
    PROCEDURE, PRIVATE :: normal_new_factorization_rf_y
    PROCEDURE, PRIVATE :: new_factorization_rf_z
    PROCEDURE, PRIVATE :: new_factorization_rf_x
    PROCEDURE, PRIVATE :: new_factorization_rf_y
    PROCEDURE, PRIVATE :: bkgcov_a_en_new_factorization

    PROCEDURE, PUBLIC  :: sqrt_inverse_multiply
    PROCEDURE, PUBLIC  :: sqrt_inverse_multiply_adjoint

    PROCEDURE :: sqrtInvMul => sqrt_inverse_multiply
    PROCEDURE :: sqrtInvMulAdj => sqrt_inverse_multiply_adjoint

    PROCEDURE, PUBLIC :: inverse_multiply

    PROCEDURE, PUBLIC :: initialize
    FINAL :: destructor
  END TYPE BFieldHybInc_t

CONTAINS

  SUBROUTINE initialize(this, configFile, sg, varName, Y)
    IMPLICIT NONE

    CLASS(BFieldHybInc_t) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg

    CHARACTER(LEN=1024), INTENT(IN)  :: configFile
    CHARACTER(LEN=10), INTENT(IN)  :: varName
    TYPE(ObsSet_t), TARGET, OPTIONAL :: Y

    INTEGER(i_kind) :: istatus, numCell

    CALL this%AuxTypeSG_t%aux_initialize(sg)
    this%name = varName

    this%scaleParaX = 1.0D0
    this%scaleParaY = 1.0D0
    this%scaleParaZ = 1.0D0
    this%scaleParaT = 1.0D0

    istatus = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaX', this%scaleParaX)
    istatus = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaY', this%scaleParaY)
    istatus = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaZ', this%scaleParaZ)
    istatus = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaT', this%scaleParaT)

    IF (PRESENT(Y)) THEN
      this%Y => Y
    END IF

    ! For hybens localization setup
    istatus = yaml_get_var(TRIM(configFile), 'BMat', 's_ens_h', this%s_ens_h)
    istatus = yaml_get_var(TRIM(configFile), 'BMat', 's_ens_v', this%s_ens_v)

    this%s_ens_hv = this%s_ens_h * 1.0D0
    this%s_ens_vv = this%s_ens_v * 1.0D0
    this%kl = 1

    istatus = yaml_get_var(TRIM(configFile), 'BMat', 's_ens_vlv', this%s_ens_vlv)

    ! Set up localization filters
    CALL this%convert_km_to_grid_units   ! convert s_ens_h from km to grid units.

    IF (this%s_ens_v > 0) THEN
      CALL this%init_rf_z
      CALL this%normal_new_factorization_rf_z
    END IF

    CALL this%init_rf_x
    CALL this%init_rf_y
    CALL this%normal_new_factorization_rf_x
    CALL this%normal_new_factorization_rf_y

    this%init_flag = .TRUE.
    CALL this%loadBMatFiles(varName)
    this%init_flag = .FALSE.

  END SUBROUTINE initialize

  SUBROUTINE loadBMatFiles(this, varName)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    CHARACTER(LEN=10), INTENT(IN) :: varName
    TYPE(BKErr_t)                 :: br
    CHARACTER(LEN=1024)           :: outputDir, datFileName
    CHARACTER(LEN=3)              :: gid
    INTEGER(i_kind)               :: status, nlvl, nhor, ntim, nz, nc, i, j, info, k
    REAL(r_kind)                  :: tmp(this%sg%vLevel, this%sg%num_cell, this%sg%tSlots)
    REAL(r_kind)                  :: tmpvalue

    ALLOCATE (this%sigma(this%sg%vLevel, this%sg%num_cell, this%sg%tSlots))
    this%sigma = 1.0D0 !0.000000001D0!1.0D0 !0.00001D0 !100.0D0 !10000.0D0
    tmp = this%sigma

    CALL this%bkgcov_a_en_new_factorization(tmp)
    tmpvalue = getMagnitude(MAXVAL(tmp))
    this%sigma = tmpvalue

  END SUBROUTINE loadBMatFiles

  SUBROUTINE inverse_multiply(this, field)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    TYPE(Field_t), INTENT(INOUT) :: field

  END SUBROUTINE inverse_multiply

  SUBROUTINE sqrt_inverse_multiply(this, field)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    TYPE(Field_t), INTENT(INOUT) :: field
    REAL(r_kind), ALLOCATABLE    :: delta_x(:, :, :), datatmp(:, :, :), w(:, :)

    ALLOCATE (delta_x(this%sg%vLevel, this%sg%num_cell, this%sg%tSlots))
    delta_x = field%DATA
    datatmp = delta_x

    CALL this%bkgcov_a_en_new_factorization(delta_x)
    field%DATA = delta_x / this%sigma

    DEALLOCATE (delta_x)

    ! PRINT *, 'DONE BField%sqrt_inverse_multiply ...'

  END SUBROUTINE sqrt_inverse_multiply

  SUBROUTINE sqrt_inverse_multiply_adjoint(this, field, currentField)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    TYPE(Field_t), OPTIONAL, INTENT(IN) :: currentField
    TYPE(Field_t), INTENT(INOUT)        :: field

    field%DATA = field%DATA

    ! PRINT *, 'DONE BField%sqrt_inverse_multiply_adjoint ...'

  END SUBROUTINE sqrt_inverse_multiply_adjoint

  SUBROUTINE convert_km_to_grid_units(this)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    REAL(r_kind), ALLOCATABLE  :: lat_2d(:, :), lon_2d(:, :)
    REAL(r_kind)               :: griddist_dx, griddist_dy
    INTEGER(i_kind)            :: nx, ny

    ny = this%sg%dimCell_global(1)
    nx = this%sg%dimCell_global(2)
    ALLOCATE (lat_2d(nx, ny), lon_2d(nx, ny))
    lat_2d = RESHAPE(this%sg%cell_cntr(1, :), (/nx, ny/))
    lon_2d = RESHAPE(this%sg%cell_cntr(2, :), (/nx, ny/))
    griddist_dy = distance_hav(lat_2d(nx / 2, ny / 2), lon_2d(nx / 2, ny / 2), lat_2d(nx / 2, ny / 2 + 1), lon_2d(nx / 2, ny / 2 + 1))
    griddist_dx = distance_hav(lat_2d(nx / 2, ny / 2), lon_2d(nx / 2, ny / 2), lat_2d(nx / 2 + 1, ny / 2), lon_2d(nx / 2 + 1, ny / 2))

    this%s_ens_h_gu_y = CEILING(this%s_ens_hv / griddist_dy)
    this%s_ens_h_gu_x = CEILING(this%s_ens_hv / griddist_dx)

    IF (ALLOCATED(lat_2d)) DEALLOCATE (lat_2d)
    IF (ALLOCATED(lon_2d)) DEALLOCATE (lon_2d)

  END SUBROUTINE convert_km_to_grid_units

  SUBROUTINE init_rf_z(this)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    INTEGER(i_kind)                               ::  k, nxy, i, ii, jj, j, l, nsig
    REAL(r_kind)                                  ::  dlnp, kap1, kapr, d1
    REAL(r_kind), ALLOCATABLE                     :: aspect(:)
    REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE :: fmatz_tmp
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE    :: fmat0z_tmp

    nxy = this%sg%dimCell_global(1) * this%sg%dimCell_global(2)
    nsig = this%sg%vLevel

    ! use new factorization:
    ALLOCATE (this%fmatz(nxy, 2, nsig, 2), this%fmat0z(nxy, nsig, 2))
    ALLOCATE (fmatz_tmp(2, nsig, 2), fmat0z_tmp(nsig, 2))
    ALLOCATE (aspect(nsig))

    IF (this%s_ens_v > 0.0D0) THEN

      ! s_ens_vv is in grid units
      IF (.NOT. this%s_ens_vlv) THEN
        DO k = 1, nsig
          aspect(k) = this%s_ens_vv**2
        END DO
      ELSE
        DO k = 1, INT(this%s_ens_vv)
          aspect(k) = (this%s_ens_vv - k + 1)**2
        END DO
        DO k = INT(this%s_ens_vv) + 1, nsig
          aspect(k) = 1.0D0
        END DO
      END IF


      DO i = 1, nxy

        CALL get_new_alpha_beta(aspect, nsig, fmatz_tmp, fmat0z_tmp)
        DO l = 1, 2
          DO k = 1, nsig
            DO j = 1, 2
              this%fmatz(i, j, k, l) = fmatz_tmp(j, k, l)
            END DO
          END DO
        END DO

        DO l = 1, 2
          DO k = 1, nsig
            this%fmat0z(i, k, l) = fmat0z_tmp(k, l)
          END DO
        END DO

      END DO

    END IF

    DEALLOCATE (fmatz_tmp, fmat0z_tmp)
    RETURN

  END SUBROUTINE init_rf_z

  SUBROUTINE init_rf_x(this)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    INTEGER(i_kind)           :: i, j, k, l, kk, nx, ny, kl
    REAL(r_kind), ALLOCATABLE :: aspect(:), fmatc(:, :, :), fmat0c(:, :)

    nx = this%sg%dimCell_global(2)
    ny = this%sg%dimCell_global(1)
    kl = this%kl

    ! use new factorization:
    IF (ALLOCATED(this%fmatx)) DEALLOCATE (this%fmatx)
    IF (ALLOCATED(this%fmat0x)) DEALLOCATE (this%fmat0x)
    ALLOCATE (this%fmatx(ny, 2, nx, 2, kl), this%fmat0x(ny, nx, 2, kl))
    ALLOCATE (fmatc(2, nx, 2), fmat0c(nx, 2), aspect(nx))
    DO k = 1, kl
      DO i = 1, ny

        DO j = 1, nx
          aspect(j) = this%s_ens_h_gu_x**2
        END DO

        CALL get_new_alpha_beta(aspect, nx, fmatc, fmat0c)

        DO kk = 1, 2
          DO j = 1, nx

            DO l = 1, 2
              this%fmatx(i, l, j, kk, k) = fmatc(l, j, kk)
            END DO
            this%fmat0x(i, j, kk, k) = fmat0c(j, kk)

          END DO
        END DO

      END DO
    END DO

    DEALLOCATE (fmatc, fmat0c, aspect)
    RETURN

  END SUBROUTINE init_rf_x

  SUBROUTINE init_rf_y(this)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    REAL(r_kind), ALLOCATABLE :: aspect(:)
    INTEGER(i_kind)           :: j, k, l, kk, ny, kl

    ny = this%sg%dimCell_global(1)
    kl = this%kl

    ! use new factorization:
    IF (ALLOCATED(this%fmaty)) DEALLOCATE (this%fmaty)
    IF (ALLOCATED(this%fmat0y)) DEALLOCATE (this%fmat0y)
    ALLOCATE (this%fmaty(2, ny, 2, kl), this%fmat0y(ny, 2, kl))
    ALLOCATE (aspect(ny))

    DO k = 1, kl

      DO j = 1, ny
        aspect(j) = this%s_ens_h_gu_y**2
      END DO

      CALL get_new_alpha_beta(aspect, ny, this%fmaty(1, 1, 1, k), this%fmat0y(1, 1, k))

    END DO

    DEALLOCATE (aspect)
    RETURN

  END SUBROUTINE init_rf_y

  SUBROUTINE normal_new_factorization_rf_z(this)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    INTEGER(i_kind)           :: k, iadvance, iback, nxy, nsig
    REAL(r_kind), ALLOCATABLE :: f(:, :), diag(:, :)

    nxy = this%sg%dimCell_global(1) * this%sg%dimCell_global(2)
    nsig = this%sg%vLevel

    IF (ALLOCATED(this%znorm_new)) DEALLOCATE (this%znorm_new)
    ALLOCATE (this%znorm_new(nxy, nsig))
    ALLOCATE (f(nsig, nxy), diag(nsig, nxy))

    this%znorm_new = 1.0D0

    DO k = 1, nsig
      f = 0.0D0
      ! f(:,k) = 1.0D0
      f(k, :) = 1.0D0

      iadvance = 1
      iback = 2
      CALL this%new_factorization_rf_z(f, iadvance, iback)
      iadvance = 2
      iback = 1
      CALL this%new_factorization_rf_z(f, iadvance, iback)

      diag(k, :) = SQRT(1.0D0 / f(k, :))
    END DO

    DO k = 1, nsig
      this%znorm_new(:, k) = diag(k, :)
    END DO

    DEALLOCATE (f, diag)
    RETURN

  END SUBROUTINE normal_new_factorization_rf_z

  SUBROUTINE normal_new_factorization_rf_x(this)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    INTEGER(i_kind)           :: i, j, k, iadvance, iback, kl, nx, ny
    REAL(r_kind), ALLOCATABLE :: f(:, :, :), diag(:, :, :)

    ny = this%sg%dimCell_global(1)
    nx = this%sg%dimCell_global(2)
    kl = this%kl

    IF (ALLOCATED(this%xnorm_new)) DEALLOCATE (this%xnorm_new)
    ALLOCATE (this%xnorm_new(ny, nx, kl))
    IF (ALLOCATED(diag)) DEALLOCATE (diag)
    ALLOCATE (diag(ny, nx, kl))
    IF (ALLOCATED(f)) DEALLOCATE (f)
    ALLOCATE (f(ny, nx, kl))

    this%xnorm_new = 1.0D0

    DO j = 1, nx
      f = 0.0D0
      DO k = 1, kl
        DO i = 1, ny
          f(i, j, k) = 1.0D0
        END DO
      END DO
      iadvance = 1
      iback = 2
      CALL this%new_factorization_rf_x(f, iadvance, iback, kl)
      iadvance = 2
      iback = 1
      CALL this%new_factorization_rf_x(f, iadvance, iback, kl)
      DO k = 1, kl
        DO i = 1, ny
          diag(i, j, k) = SQRT(1.0D0 / f(i, j, k))
        END DO
      END DO
    END DO

    DO k = 1, kl
      DO j = 1, nx
        DO i = 1, ny
          this%xnorm_new(i, j, k) = diag(i, j, k)
        END DO
      END DO
    END DO

    DEALLOCATE (diag, f)
    RETURN

  END SUBROUTINE normal_new_factorization_rf_x

  SUBROUTINE normal_new_factorization_rf_y(this)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    INTEGER(i_kind)           :: i, k, lend, lcount, iadvance, iback, kl, loop, ll, iend, ny, nx
    REAL(r_kind), ALLOCATABLE :: f(:, :, :), diag(:, :)

    ny = this%sg%dimCell_global(1)
    nx = this%sg%dimCell_global(2)
    kl = this%kl

    IF (ALLOCATED(this%ynorm_new)) DEALLOCATE (this%ynorm_new)
    ALLOCATE (this%ynorm_new(ny, kl))
    IF (ALLOCATED(diag)) DEALLOCATE (diag)
    ALLOCATE (diag(ny, kl))
    IF (ALLOCATED(f)) DEALLOCATE (f)
    ALLOCATE (f(ny, nx, kl))

    this%ynorm_new = 1.0D0

    IF (ny <= nx) THEN
      lend = 1
      iend = ny
    ELSE
      lend = ny / nx + 1
      iend = nx
    END IF

    DO loop = 1, lend
      ll = (loop - 1) * iend
      f = 0.0D0
      DO k = 1, kl
        DO i = 1, iend
          lcount = ll + i
          f(lcount, i, k) = 1.0D0
          IF (lcount == ny) EXIT
        END DO
      END DO

      iadvance = 1
      iback = 2
      CALL this%new_factorization_rf_y(f, iadvance, iback, kl)
      iadvance = 2
      iback = 1
      CALL this%new_factorization_rf_y(f, iadvance, iback, kl)

      DO k = 1, kl
        DO i = 1, iend
          lcount = ll + i
          diag(lcount, k) = SQRT(1.0D0 / f(lcount, i, k))
          this%ynorm_new(lcount, k) = diag(lcount, k)
          IF (lcount == ny) EXIT
        END DO
      END DO
    END DO

    DEALLOCATE (diag, f)
    RETURN

  END SUBROUTINE normal_new_factorization_rf_y

  SUBROUTINE new_factorization_rf_z(this, f, iadvance, iback)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    INTEGER(i_kind), INTENT(IN) :: iadvance, iback
    REAL(r_kind), INTENT(INOUT) :: f(:, :)
    ! REAL(r_kind)   , INTENT(INOUT) :: f(:,:)

    INTEGER(i_kind) :: i, k, l, nxy, nz

    nxy = this%sg%dimCell_global(1) * this%sg%dimCell_global(2)
    nz = this%sg%vLevel

    ! PRINT *, 'PJM-DEBUG --- in Sub new_factorization_rf_z: nz = ', nz

    IF (iadvance == 1) THEN
      DO k = 1, nz
        DO i = 1, nxy
          ! f(i,k) = this%znorm_new(i,k) * f(i,k)
          f(k, i) = this%znorm_new(i, k) * f(k, i)
        END DO
      END DO
    END IF
    DO k = 1, nz
      DO l = 1, MIN(2, k - 1)
        DO i = 1, nxy
          ! f(i,k) = f(i,k) - this%fmatz(i,l,k,iadvance)*f(i,k-l)
          f(k, i) = f(k, i) - this%fmatz(i, l, k, iadvance) * f(k - l, i)
        END DO
      END DO
      DO i = 1, nxy
        ! f(i,k) = this%fmat0z(i,k,iadvance) * f(i,k)
        f(k, i) = this%fmat0z(i, k, iadvance) * f(k, i)
      END DO
    END DO
    DO k = nz, 1, -1
      DO l = 1, MIN(2, nz - k)
        DO i = 1, nxy
          ! f(i,k) = f(i,k) - this%fmatz(i,l,k+l,iback)*f(i,k+l)
          f(k, i) = f(k, i) - this%fmatz(i, l, k + l, iback) * f(k + l, i)
        END DO
      END DO
      DO i = 1, nxy
        ! f(i,k) = this%fmat0z(i,k,iback) * f(i,k)
        f(k, i) = this%fmat0z(i, k, iback) * f(k, i)
      END DO
    END DO
    IF (iadvance == 2) THEN
      DO k = 1, nz
        DO i = 1, nxy
          ! f(i,k) = this%znorm_new(i,k) * f(i,k)
          f(k, i) = this%znorm_new(i, k) * f(k, i)
        END DO
      END DO
    END IF

    RETURN

  END SUBROUTINE new_factorization_rf_z

  SUBROUTINE new_factorization_rf_x(this, f, iadvance, iback, nlevs)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    INTEGER(i_kind), INTENT(IN) :: iadvance, iback, nlevs
    REAL(r_kind), INTENT(INOUT) :: f(:, :, :)

    INTEGER(i_kind) :: i, j, k, l, ny, nx, nz

    ny = this%sg%dimCell_global(1)
    nx = this%sg%dimCell_global(2)
    nz = nlevs

    DO k = 1, nz

      IF (iadvance == 1) THEN
        DO j = 1, nx
          DO i = 1, ny
            f(i, j, k) = this%xnorm_new(i, j, 1) * f(i, j, k)
          END DO
        END DO
      END IF

      DO j = 1, nx
        DO l = 1, MIN(2, j - 1)
          DO i = 1, ny
            f(i, j, k) = f(i, j, k) - this%fmatx(i, l, j, iadvance, 1) * f(i, j - l, k)
          END DO
        END DO
        DO i = 1, ny
          f(i, j, k) = this%fmat0x(i, j, iadvance, 1) * f(i, j, k)
        END DO
      END DO

      DO j = nx, 1, -1
        DO l = 1, MIN(2, nx - j)
          DO i = 1, ny
            f(i, j, k) = f(i, j, k) - this%fmatx(i, l, j + l, iback, 1) * f(i, j + l, k)
          END DO
        END DO
        DO i = 1, ny
          f(i, j, k) = this%fmat0x(i, j, iback, 1) * f(i, j, k)
        END DO
      END DO

      IF (iadvance == 2) THEN
        DO j = 1, nx
          DO i = 1, ny
            f(i, j, k) = this%xnorm_new(i, j, 1) * f(i, j, k)
          END DO
        END DO
      END IF

    END DO
    RETURN
  END SUBROUTINE new_factorization_rf_x

  SUBROUTINE new_factorization_rf_y(this, f, iadvance, iback, nlevs)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    INTEGER(i_kind), INTENT(IN) :: iadvance, iback, nlevs
    REAL(r_kind), INTENT(INOUT) :: f(:, :, :)

    INTEGER(i_kind) :: i, j, k, l, nx, ny, nz

    ny = this%sg%dimCell_global(1)
    nx = this%sg%dimCell_global(2)
    nz = nlevs

    DO k = 1, nz
      DO j = 1, nx

        IF (iadvance == 1) THEN
          DO i = 1, ny
            f(i, j, k) = this%ynorm_new(i, 1) * f(i, j, k)
          END DO
        END IF

        DO i = 1, ny
          DO l = 1, MIN(2, i - 1)
            f(i, j, k) = f(i, j, k) - this%fmaty(l, i, iadvance, 1) * f(i - l, j, k)
          END DO
          f(i, j, k) = this%fmat0y(i, iadvance, 1) * f(i, j, k)
        END DO

        DO i = ny, 1, -1
          DO l = 1, MIN(2, ny - i)
            f(i, j, k) = f(i, j, k) - this%fmaty(l, i + l, iback, 1) * f(i + l, j, k)
          END DO
          f(i, j, k) = this%fmat0y(i, iback, 1) * f(i, j, k)
        END DO

        IF (iadvance == 2) THEN
          DO i = 1, ny
            f(i, j, k) = this%ynorm_new(i, 1) * f(i, j, k)
          END DO
        END IF

      END DO
    END DO
    RETURN
  END SUBROUTINE new_factorization_rf_y

  SUBROUTINE get_new_alpha_beta(aspect, ng, fmat_out, fmat0_out)
    IMPLICIT NONE

    INTEGER(i_kind), INTENT(IN)  :: ng
    REAL(r_kind), INTENT(IN)     :: aspect(ng)
    REAL(r_kind), INTENT(OUT)    :: fmat_out(2, ng, 2), fmat0_out(ng, 2)
    INTEGER(i_kind)              :: i
    REAL(r_kind)                 :: sig(0:ng - 1), fmat(0:ng - 1, -2:0, 2)

    DO i = 1, ng
      sig(i - 1) = SQRT(aspect(i))
    END DO
    CALL stringop(ng - 1, sig, fmat)

    DO i = 1, ng
      fmat_out(2, i, 1) = fmat(i - 1, -2, 1)
      fmat_out(1, i, 1) = fmat(i - 1, -1, 1)
      fmat0_out(i, 1) = 1.0D0 / fmat(i - 1, 0, 1)
      fmat_out(2, i, 2) = fmat(i - 1, -2, 2)
      fmat_out(1, i, 2) = fmat(i - 1, -1, 2)
      fmat0_out(i, 2) = 1.0D0 / fmat(i - 1, 0, 2)
    END DO
    RETURN

  END SUBROUTINE get_new_alpha_beta

  SUBROUTINE bkgcov_a_en_new_factorization(this, a_en)
    IMPLICIT NONE
    CLASS(BFieldHybInc_t) :: this
    REAL(r_kind), INTENT(INOUT) :: a_en(:, :, :)
    INTEGER(i_kind)             :: ii, k, iflg, iadvance, iback, is, ie, ipnt, istatus, nx, ny, kl, nt, t, nxy, nz
    REAL(r_kind), ALLOCATABLE   :: a_en_tmp(:, :, :, :)

    ny = this%sg%dimCell_global(1)
    nx = this%sg%dimCell_global(2)
    nz = this%sg%vLevel
    kl = this%kl
    nt = this%sg%tSlots
    nxy = nx * ny

    ALLOCATE (a_en_tmp(ny, nx, nz, nt))

    ! Apply vertical smoother on each ensemble member
    IF (this%s_ens_v > 0) THEN
      iadvance = 1
      iback = 2
      DO t = 1, nt
        CALL this%new_factorization_rf_z(a_en(:, :, t), iadvance, iback)
      END DO
    END IF

    ! Convert from domain to full horizontal field
    DO t = 1, nt
      DO k = 1, nz
        a_en_tmp(:, :, k, t) = TRANSPOSE(RESHAPE(a_en(k, :, t), (/nx, ny/)))
      END DO
    END DO

    ! Apply horizontal smoother for number of horizontal scales
    DO t = 1, nt
      iadvance = 1
      iback = 2
      CALL this%new_factorization_rf_x(a_en_tmp(:, :, :, t), iadvance, iback, nz)
      CALL this%new_factorization_rf_y(a_en_tmp(:, :, :, t), iadvance, iback, nz)

      iadvance = 2
      iback = 1
      CALL this%new_factorization_rf_y(a_en_tmp(:, :, :, t), iadvance, iback, nz)
      CALL this%new_factorization_rf_x(a_en_tmp(:, :, :, t), iadvance, iback, nz)
    END DO

    ! Convert back from full horizontal field to domain
    DO t = 1, nt
      DO k = 1, nz
        a_en(k, :, t) = RESHAPE(TRANSPOSE(a_en_tmp(:, :, k, t)), (/nxy/))
      END DO
    END DO

    ! Retrieve ensemble components from long vector
    ! Apply vertical smoother on each ensemble member
    IF (this%s_ens_v > 0) THEN
      iadvance = 2
      iback = 1
      DO t = 1, nt
        CALL this%new_factorization_rf_z(a_en(:, :, t), iadvance, iback)
      END DO
    END IF

    DEALLOCATE (a_en_tmp)
    RETURN

  END SUBROUTINE bkgcov_a_en_new_factorization

  SUBROUTINE stringop(n, sig, fmat)
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN)                        :: n
    REAL(r_kind), DIMENSION(0:n), INTENT(IN)           :: sig
    REAL(r_kind), DIMENSION(0:n, -2:0, 2), INTENT(OUT) :: fmat
    COMPLEX(r_kind), DIMENSION(2)                      :: rts
    REAL(r_kind), DIMENSION(0:n, -1:1)                 :: dmat
    REAL(r_kind)                                       :: dmatt
    INTEGER(i_kind)                                    :: i, im
    DATA rts/ &
      (-3.4588884621354108D0, 1.7779487522437312D0), &
      (-0.54111153786458903D0, 5.0095518087248694D0)/
      

    dmat = 0 !.0D0
    DO i = 1, n
      im = i - 1
      dmatt = sig(im) * sig(i)
      dmat(im, 1) = -dmatt
      dmat(i, -1) = -dmatt
      dmat(im, 0) = dmat(im, 0) + dmatt
      dmat(i, 0) = dmat(i, 0) + dmatt
    END DO
    DO i = 1, 2
      CALL quadfil(n, dmat, rts(i), fmat(:, :, i))
    END DO
  END SUBROUTINE stringop

  SUBROUTINE quadfil(n, dmat, rts, fmat)
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN)                     :: n
    REAL(r_kind), DIMENSION(0:n, -1:1), INTENT(IN)  :: dmat
    COMPLEX(r_kind), INTENT(IN)                     :: rts
    REAL(r_kind), DIMENSION(0:n, -2:0), INTENT(OUT) :: fmat
    REAL(r_kind), DIMENSION(0:n, -1:1)              :: tmat
    REAL(r_kind), DIMENSION(0:n, -2:2)              :: umat
    REAL(r_kind)                                    :: a1, a2, rrtsi, qrtsi
    COMPLEX(r_kind)                                 :: rtsi
    INTEGER(i_kind)                                 :: np
    

    np = n + 1
    rtsi = 1.0D0 / rts
    rrtsi = REAL(rtsi)
    qrtsi = AIMAG(rtsi)
    a1 = -2.0D0 * rrtsi
    a2 = rrtsi**2 + qrtsi**2
    tmat(0:n, -1:1) = a2 * dmat
    tmat(0:n, 0) = tmat(0:n, 0) + a1
    CALL mulbb(dmat, tmat(0:n, -1:1), umat(0:n, -2:2), np, np, 1, 1, 1, 1, 2, 2)
    umat(0:n, 0) = umat(0:n, 0) + 1.0D0
    IF (n == 0) THEN
      fmat = 0.0D0
      fmat(0, 0) = umat(0, 0)
    ELSE
      CALL l1lb(umat, fmat, np, 2)
    END IF

  END SUBROUTINE quadfil

  SUBROUTINE mulbb(a, b, c, m1, m2, mah1, mah2, mbh1, mbh2, mch1, mch2)
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN)  :: m1, m2, mah1, mah2, mbh1, mbh2, mch1, mch2
    REAL(r_kind), INTENT(IN)     :: a(m1, -mah1:mah2), b(m2, -mbh1:mbh2)
    REAL(r_kind), INTENT(INOUT)  :: c(m1, -mch1:mch2)
    INTEGER(i_kind)              :: nch1, nch2, j, k, jpk, i1, i2

    c = 0.0D0
    nch1 = mah1 + mbh1
    nch2 = mah2 + mbh2
    IF (nch1 /= mch1 .OR. nch2 /= mch2) STOP 'In MULBB, dimensions inconsistent'
    DO j = -mah1, mah2
      DO k = -mbh1, mbh2
        jpk = j + k
        i1 = MAX(1, 1 - j)
        i2 = MIN(m1, m2 - j)
        c(i1:i2, jpk) = c(i1:i2, jpk) + a(i1:i2, j) * b(j + i1:j + i2, k)
      END DO
    END DO
  END SUBROUTINE mulbb

  SUBROUTINE l1lb(a, b, m, mah)   ! Cholesky LU decomposition of Banded.
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN)  :: m, mah
    REAL(r_kind), INTENT(IN)     :: a(m, -mah:mah)
    REAL(r_kind), INTENT(OUT)    :: b(m, -mah:0)
    INTEGER(i_kind)              :: i, j, jmi
    REAL(r_kind)                 :: s

    CALL clib(b, m, m, mah, 0)
    DO j = 1, m
      s = a(j, 0) - DOT_PRODUCT(b(j, -mah:-1), b(j, -mah:-1))
      IF (s <= 0.0D0) THEN
        PRINT '(" L1LB detects non-positivity at diagonal index",i5)', j
        STOP
      END IF
      s = SQRT(s)
      b(j, 0) = s
      s = 1.0D0 / s
      DO i = j + 1, MIN(m, j + mah)
        jmi = j - i
        b(i, jmi) = s * (a(i, jmi) - DOT_PRODUCT(b(i, -mah:jmi - 1), b(j, -mah - jmi:-1)))
      END DO
    END DO
  END SUBROUTINE l1lb

  SUBROUTINE clib(a, m1, m2, mah1, mah2) ! Clip the dead space of the band matrix, a
    IMPLICIT NONE
    INTEGER(i_kind), INTENT(IN)  :: m1, m2, mah1, mah2
    REAL(r_kind), INTENT(INOUT)  :: a(m1, -mah1:mah2)
    INTEGER(i_kind)              :: j

    IF (m2 - m1 + mah1 < 0) STOP 'In CLIB, form of band matrix implies redundant rows'
    DO j = 1, mah1
      a(1:MIN(m1, j), -j) = 0.0D0
    END DO
    DO j = m2 - m1 + 1, mah2
      a(MAX(1, m2 - j + 1):m1, j) = 0.0D0
    END DO
  END SUBROUTINE clib

  FUNCTION getMagnitude(x) RESULT(magnitude)
    REAL(r_kind), INTENT(IN) :: x
    REAL(r_kind) :: magnitude

    magnitude = 10.0D0**(INT(LOG10(REAL(x))))

  END FUNCTION getMagnitude

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(BFieldHybInc_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%sigma)) DEALLOCATE (this%sigma)
  END SUBROUTINE destructor

END MODULE BFieldHybInc_m
