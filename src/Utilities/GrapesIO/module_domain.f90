!!--------------------------------------------------------------------------------------------------
! PROJECT           : GRAPES IO
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Sanshan Tu
! VERSION           : Beta 0.0
! HISTORY           :
!   Created by Sanshan Tu (tss71618@163.com), 2020/12/31, @SZSC, Shenzhen
!   Modified by Sanshan Tu (tss71618@163.com), 2021/11/01, @SZSC, Shenzhen
!  Modified by Zilong Qin (zilong.qin@gmail.com), 2022/5/13, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!!===================================================================
!> @brief
!! # GRAPES IO Module
!!
!!  *This module defines data structures for GRAPES input*
!! @author Sanshan Tu
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
!!===================================================================
MODULE module_domain
  USE module_configure
  USE kinds_m
  IMPLICIT NONE
  TYPE domain_t

    type(grid_config_rec_type)::config
    REAL(r_single), DIMENSION(:), ALLOCATABLE   :: lon !! longitudes
    REAL(r_single), DIMENSION(:), ALLOCATABLE   :: lat !! latitudes

    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: u !! u-component of wind speed (m/s)
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: v !! v-component of wind speed (m/s)

    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: q  !! moisture, mois_2 in grapes
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: qc  !! moisture, mois_3 in grapes
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: qr  !! moisture, mois_4 in grapes
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: qi  !! moisture, mois_5 in grapes
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: qs  !! moisture, mois_6 in grapes
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: qg  !! moisture, mois_7 in grapes

    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: zz !! zz

    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: pip !! pip
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: pi !! pi^M

    ! REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE   :: T  !! Temperature in K from mois_2 and zz in grapes

    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: thp !! thp
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: th !! th

    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: wzet !! wzet

    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE   :: res   !! used to store the return of flexible reading
    ! REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE   :: rho   !! density of the air
    ! REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE   :: p   !! pressure of the air

    REAL(r_single), DIMENSION(:, :), ALLOCATABLE   :: tsk!! tst
    REAL(r_single), DIMENSION(:, :), ALLOCATABLE   :: xland !! xland

    REAL(r_single), DIMENSION(:, :), ALLOCATABLE   :: snowc !! snowc
    REAL(r_single), DIMENSION(:, :), ALLOCATABLE   :: xice !! xice

    REAL(r_single), DIMENSION(:, :), ALLOCATABLE   :: soil_type !! soil_type
    REAL(r_single), DIMENSION(:, :), ALLOCATABLE   :: veg_fraction !! veg_fraction

    REAL(r_single), DIMENSION(:, :), ALLOCATABLE   :: ht !! ht
    REAL(r_single), DIMENSION(:, :), ALLOCATABLE   :: xlat !! xlat
    REAL(r_single), DIMENSION(:, :), ALLOCATABLE   :: xlong !! xlong

    REAL(r_single), DIMENSION(:, :), ALLOCATABLE   :: u10 !! xlong
    REAL(r_single), DIMENSION(:, :), ALLOCATABLE   :: v10 !! xlong
    REAL(r_single), DIMENSION(:, :), ALLOCATABLE   :: t2 !! xlong
    REAL(r_single), DIMENSION(:, :), ALLOCATABLE   :: q2 !! xlong
    REAL(r_single), DIMENSION(:, :), ALLOCATABLE   :: ps !! atmospheric pressure on ground from zz in grapes

    REAL(r_single), DIMENSION(:), ALLOCATABLE   :: zh

    INTEGER(i_kind)    :: ids, ide, jds, jde, kds, kde

  CONTAINS
    FINAL :: domain_destructor
    PROCEDURE, PUBLIC, PASS :: getLocalIndices
    PROCEDURE, PUBLIC, PASS :: alloc_grid_memory
    PROCEDURE, PUBLIC, PASS :: dealloc_grid_memory

  END TYPE domain_t

  INTERFACE domain_t
    PROCEDURE :: domain_constructor
  END INTERFACE

CONTAINS
  FUNCTION domain_constructor(nlFileName) result(this)
    type(domain_t) :: this
    CHARACTER*(*), INTENT(IN) :: nlFileName

    CALL initial_config(TRIM(nlFileName), this%config)
    CALL this%getLocalIndices()

  END FUNCTION domain_constructor

  SUBROUTINE getLocalIndices(this)
    IMPLICIT NONE
    CLASS(domain_t) :: this

    this%ids = this%config%s_we
    this%ide = this%config%e_we
    this%jds = this%config%s_sn
    this%jde = this%config%e_sn
    this%kds = this%config%s_vert - 1
    this%kde = this%config%e_vert + 1

  END SUBROUTINE getLocalIndices

  SUBROUTINE alloc_grid_memory(this)
    IMPLICIT NONE
    CLASS(domain_t) :: this
    INTEGER(i_kind) i, j, ids, ide, kds, kde, jds, jde
    INTEGER ierr

    ids = this%ids
    ide = this%ide
    kds = this%kds
    kde = this%kde
    jds = this%jds
    jde = this%jde

    ALLOCATE (this%lon(ids:ide))
    ALLOCATE (this%lat(jds:jde))
    DO i = this%ids, this%ide
      this%lon(i) = this%config%xs_we + (i - 1)*this%config%xd
    END DO

    DO j = this%jds, this%jde
      this%lat(j) = this%config%ys_sn + (j - 1)*this%config%yd
    END DO

    ALLOCATE (this%u(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%v(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ! ALLOCATE (this%T(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%q(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%qc(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%qr(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%qi(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%qs(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%qg(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%zz(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%pip(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%pi(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%wzet(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%thp(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%th(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%tsk(ids:ide, jds:jde), STAT=ierr)
    ALLOCATE (this%xland(ids:ide, jds:jde), STAT=ierr)
    ALLOCATE (this%ht(ids:ide, jds:jde), STAT=ierr)
    ALLOCATE (this%zh(kds:kde), STAT=ierr)
    ALLOCATE (this%res(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ALLOCATE (this%xlat(ids:ide, jds:jde), STAT=ierr)
    ALLOCATE (this%xlong(ids:ide, jds:jde), STAT=ierr)

    ALLOCATE (this%u10(ids:ide, jds:jde), STAT=ierr)
    ALLOCATE (this%v10(ids:ide, jds:jde), STAT=ierr)
    ALLOCATE (this%t2(ids:ide, jds:jde), STAT=ierr)
    ALLOCATE (this%q2(ids:ide, jds:jde), STAT=ierr)
    ALLOCATE (this%ps(ids:ide, jds:jde), STAT=ierr)
    ! ALLOCATE (this%rho(ids:ide, kds:kde, jds:jde), STAT=ierr)
    ! ALLOCATE (this%p(ids:ide, kds:kde, jds:jde), STAT=ierr)

    ALLOCATE (this%snowc(ids:ide, jds:jde), STAT=ierr)
    ALLOCATE (this%xice(ids:ide, jds:jde), STAT=ierr)
    ALLOCATE (this%soil_type(ids:ide, jds:jde), STAT=ierr)
    ALLOCATE (this%veg_fraction(ids:ide, jds:jde), STAT=ierr)
  END SUBROUTINE alloc_grid_memory

  SUBROUTINE dealloc_grid_memory(this)
    IMPLICIT NONE
    CLASS(domain_t) :: this
    integer::ierr

    IF (ALLOCATED(this%u)) DEALLOCATE (this%u, STAT=ierr)
    IF (ALLOCATED(this%v)) DEALLOCATE (this%v, STAT=ierr)
    ! IF (ALLOCATED(this%T)) DEALLOCATE (this%T, STAT=ierr)
    IF (ALLOCATED(this%ps)) DEALLOCATE (this%ps, STAT=ierr)
    IF (ALLOCATED(this%q)) DEALLOCATE (this%q, STAT=ierr)
    IF (ALLOCATED(this%qc)) DEALLOCATE (this%qc, STAT=ierr)
    IF (ALLOCATED(this%qr)) DEALLOCATE (this%qr, STAT=ierr)
    IF (ALLOCATED(this%qi)) DEALLOCATE (this%qi, STAT=ierr)
    IF (ALLOCATED(this%qs)) DEALLOCATE (this%qs, STAT=ierr)
    IF (ALLOCATED(this%qg)) DEALLOCATE (this%qg, STAT=ierr)
    IF (ALLOCATED(this%zz)) DEALLOCATE (this%zz, STAT=ierr)
    IF (ALLOCATED(this%pip)) DEALLOCATE (this%pip, STAT=ierr)
    IF (ALLOCATED(this%pi)) DEALLOCATE (this%pi, STAT=ierr)
    IF (ALLOCATED(this%wzet)) DEALLOCATE (this%wzet, STAT=ierr)
    IF (ALLOCATED(this%wzet)) DEALLOCATE (this%res, STAT=ierr)
    IF (ALLOCATED(this%thp)) DEALLOCATE (this%thp, STAT=ierr)
    IF (ALLOCATED(this%th)) DEALLOCATE (this%th, STAT=ierr)
    IF (ALLOCATED(this%tsk)) DEALLOCATE (this%wzet, STAT=ierr)
    IF (ALLOCATED(this%xland)) DEALLOCATE (this%xland, STAT=ierr)
    IF (ALLOCATED(this%zh)) DEALLOCATE (this%zh, STAT=ierr)
    IF (ALLOCATED(this%ht)) DEALLOCATE (this%ht, STAT=ierr)
    IF (ALLOCATED(this%xlat)) DEALLOCATE (this%xlat, STAT=ierr)
    IF (ALLOCATED(this%xlong)) DEALLOCATE (this%xlong, STAT=ierr)

    IF (ALLOCATED(this%u10)) DEALLOCATE (this%u10, STAT=ierr)
    IF (ALLOCATED(this%v10)) DEALLOCATE (this%v10, STAT=ierr)
    IF (ALLOCATED(this%t2)) DEALLOCATE (this%t2, STAT=ierr)
    IF (ALLOCATED(this%q2)) DEALLOCATE (this%q2, STAT=ierr)
    IF (ALLOCATED(this%ps)) DEALLOCATE (this%ps, STAT=ierr)
    ! IF (ALLOCATED(this%rho)) DEALLOCATE (this%rho, STAT=ierr)
    ! IF (ALLOCATED(this%p)) DEALLOCATE (this%p, STAT=ierr)

    IF (ALLOCATED(this%snowc)) DEALLOCATE (this%snowc, STAT=ierr)
    IF (ALLOCATED(this%xice)) DEALLOCATE (this%xice, STAT=ierr)
    IF (ALLOCATED(this%soil_type)) DEALLOCATE (this%soil_type, STAT=ierr)
    IF (ALLOCATED(this%veg_fraction)) DEALLOCATE (this%veg_fraction, STAT=ierr)

  END SUBROUTINE dealloc_grid_memory

  IMPURE ELEMENTAL SUBROUTINE domain_destructor(this)
    IMPLICIT NONE
    TYPE(domain_t), INTENT(INOUT) :: this
    call this%dealloc_grid_memory()
  END SUBROUTINE domain_destructor

END MODULE module_domain
