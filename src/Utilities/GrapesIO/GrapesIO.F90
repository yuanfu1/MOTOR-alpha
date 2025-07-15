!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Sanshan Tu
! VERSION           : V 1.0
! HISTORY           :
!   Created by Sanshan Tu, 2021/3/22, @GBA-MWF, Shenzhen
! Modified by Zilong Qin (zilong.qin@gmail.com), 2021/3/22, @GBA-MWF, Shenzhen
! Modified by Zilong Qin (zilong.qin@gmail.com), 2021/10/22, @GBA-MWF, Shenzhen
! Modified by Sanshan Tu (tss71618@163.com), 2021/11/01, @SZSC, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
! @note
!! In version V1.0, the c-grid u and v are converted to a-grid u and v horizontally. In the vertical
!! direction, the input is remaining c-grid. Here are the list of staggered and unstaggered variables:
!! staggered:
!! Unstaggered:
! @warning
! @attention
MODULE GrapesIO_m
  USE module_domain, ONLY: domain_t
  USE module_configure
  USE NameValueMap_m
  USE module_io_grapes_m
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE parameters_m, ONLY: spec_heat_const_pres, dry_air_gas_const, surface_ref_pres, g
  USE Interp1D_m, ONLY: interp1d_3D_idx2_single, interp1d_3D_idx2
  ! USE Filter_m, ONLY: smoothField, guidedfilter
  USE module_io_grapes_m, ONLY: module_io_grapes_t
  USE module_io_grapes_CMA_GD_V3p2_m, ONLY: module_io_grapes_CMA_GD_V3p2_t
  USE module_io_grapes_CMA_GD_V3_m, ONLY: module_io_grapes_CMA_GD_V3_t
  USE module_io_grapes_CMA_MESO_m, ONLY: module_io_grapes_CMA_MESO_t

  IMPLICIT NONE

  TYPE GrapesIO_t

    TYPE(domain_t)            :: hgrid
    TYPE(NameValueMap_t)      :: valName2FileDisp
    CHARACTER(LEN=1024)       :: nlFileName, giFileName
    CLASS(module_io_grapes_t), ALLOCATABLE :: module_io_grapes

    REAL(r_single), ALLOCATABLE :: lat(:), lon(:), &
                                   zSigma_s(:), &
                                   zSigma_u(:), &
                                   zRHght_s(:, :, :), &            !< In dimension of (vLevel, numlatlon)
                                   zRHght_u(:, :, :)!               < 2D to 1D meshGrid
    INTEGER(i_kind) :: idn, jdn, kdn, kds_s, kds_u, kde_s, kde_u
    INTEGER(i_kind) :: ids, ide, jds, jde, kds, kde

  CONTAINS
    FINAL :: destructor
    PROCEDURE, PRIVATE :: uc2ua
    PROCEDURE, PRIVATE :: vc2va
    PROCEDURE, PRIVATE :: ua2uc
    PROCEDURE, PRIVATE :: va2vc
    PROCEDURE, PRIVATE :: cal_zHghtSigma
    PROCEDURE, PUBLIC :: releaseHGrid
    PROCEDURE, PUBLIC :: writeBackToQcqrInput
    PROCEDURE, PUBLIC :: writeBackToGrapesInput
    PROCEDURE, PRIVATE :: interpUnstaggerBottomTopGrid
  END TYPE GrapesIO_t

  INTERFACE GrapesIO_t
    PROCEDURE :: constructor
  END INTERFACE

CONTAINS
  FUNCTION constructor(nlFileName, giFileName, qcqrFileName, isTest, grapes_model_type) result(this)
    TYPE(GrapesIO_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: nlFileName, giFileName
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: qcqrFileName
    INTEGER*4  :: mdate(5), century_year, month, day, hour, i, j, k
    INTEGER                 :: fh !file handle
    INTEGER*8               :: disp
    REAL(r_kind), ALLOCATABLE :: lat2D(:, :), lon2D(:, :)
    LOGICAL, OPTIONAL :: isTest
    CHARACTER(LEN=1024) :: model_type
    CHARACTER(LEN=*), OPTIONAL :: grapes_model_type

    this%hgrid = domain_t(TRIM(nlFileName))
    ! this%valName2FileDisp = init_grapes_map_modes(nlFileName, giFileName, this%hgrid, 'sequential')

    this%giFileName = giFileName
    this%nlFileName = nlFileName

    IF (PRESENT(grapes_model_type)) THEN
      model_type = TRIM(grapes_model_type)
    ELSE
      model_type = 'CMA-GD'
    END IF

    IF (model_type == 'CMA-GD') THEN
      ALLOCATE (module_io_grapes_CMA_GD_V3p2_t:: this%module_io_grapes)
    ELSEIF (model_type == 'CMA-GD-V3') THEN
      ALLOCATE (module_io_grapes_CMA_GD_V3_t:: this%module_io_grapes)
    ELSEIF (model_type == 'CMA-MESO') THEN
      ALLOCATE (module_io_grapes_CMA_MESO_t:: this%module_io_grapes)
    ELSE
      WRITE (*, *) 'model_type is not supported!'
      STOP
    END IF

    CALL this%module_io_grapes%initialize(this%hgrid, nlFileName, giFileName)
    CALL this%module_io_grapes%input_grapes_sequential(TRIM(nlFileName), TRIM(giFileName), this%hgrid)

    this%idn = this%hgrid%ide - this%hgrid%ids + 1
    this%jdn = this%hgrid%jde - this%hgrid%jds + 1
    this%kdn = this%hgrid%kde - this%hgrid%kds + 1
    this%ids = this%hgrid%ids; this%ide = this%hgrid%ide
    this%jds = this%hgrid%jds; this%jde = this%hgrid%jde
    this%kds = this%hgrid%kds; this%kde = this%hgrid%kde

    ! Define the
    ALLOCATE (this%lat(this%idn * this%jdn), this%lon(this%idn * this%jdn))

    ! Define the latitude and longitude: 1D and 2D
    FORALL (i=this%hgrid%ids:this%hgrid%ide, j=this%hgrid%jds:this%hgrid%jde)
      this%lat((j - 1) * this%idn + i) = this%hgrid%config%ys_sn + (j - 1) * this%hgrid%config%yd
      this%lon((j - 1) * this%idn + i) = this%hgrid%config%xs_we + (i - 1) * this%hgrid%config%xd
    END FORALL

    ! Calculate the height of sigma coordinates
    CALL this%cal_zHghtSigma

    IF (.NOT. PRESENT(isTest)) THEN
      ! uc to ua, and vc to va
      CALL this%uc2ua(this%hgrid%u)
      CALL this%vc2va(this%hgrid%v)

      CALL this%interpUnstaggerBottomTopGrid(this%hgrid%u)
      CALL this%interpUnstaggerBottomTopGrid(this%hgrid%v)
      CALL this%interpUnstaggerBottomTopGrid(this%hgrid%pi)
      CALL this%interpUnstaggerBottomTopGrid(this%hgrid%pip)

      FORALL (i=this%hgrid%ids:this%hgrid%ide, j=this%hgrid%jds:this%hgrid%jde, &
              k=this%hgrid%kds:this%hgrid%kde, this%hgrid%q(i, k, j) <= 0) this%hgrid%q(i, k, j) = 1.0D-10
      PRINT *, 'MINVAL(this%hgrid%q): ', MINVAL(this%hgrid%q)
    END IF

    IF (PRESENT(qcqrFileName)) THEN
      CALL this%module_io_grapes%input_grapes_hydrometeors(TRIM(qcqrFileName), this%hgrid)
    END IF
    ! ! Concatenate the surface with non-surface levels.
    ! CALL this%concat_sfc_to_3dfield
  END FUNCTION constructor

  SUBROUTINE interpUnstaggerBottomTopGrid(this, buf3D)
    CLASS(GrapesIO_t) :: this
    REAL(r_single) :: buf3D(this%hgrid%ids:this%hgrid%ide, &
                            this%hgrid%kds:this%hgrid%kde, &
                            this%hgrid%jds:this%hgrid%jde)
    REAL(r_single) :: ipCoe

    ipCoe = (this%zSigma_u(this%kds_u) - this%zSigma_u(this%kds_u + 1)) &
            / (this%zSigma_u(this%kds_u + 1) - this%zSigma_u(this%kds_u + 2))

    buf3D(:, this%kds_u, :) = (ipCoe + 1) * buf3D(:, this%kds_u + 1, :) - &
                              ipCoe * buf3D(:, this%kds_u + 2, :)

    ipCoe = (this%zSigma_u(this%kde_u) - this%zSigma_u(this%kde_u - 1)) &
            / (this%zSigma_u(this%kde_u - 1) - this%zSigma_u(this%kde_u - 2))

    buf3D(:, this%kde_u, :) = (ipCoe + 1) * buf3D(:, this%kde_u - 1, :) - &
                              ipCoe * buf3D(:, this%kde_u - 2, :)
  END SUBROUTINE

  SUBROUTINE cal_zHghtSigma(this)
    IMPLICIT NONE

    CLASS(GrapesIO_t) :: this
    INTEGER(i_kind) :: i, j, k
    REAL(r_kind) :: ztop

    this%kds_s = this%hgrid%kds + 1
    this%kde_s = this%hgrid%kde

    this%kds_u = this%hgrid%kds
    this%kde_u = this%hgrid%kde

    PRINT *, "read zSigma_s begin"
    ! call this%read_value("zh", this%zSigma_s, this%hgrid)
    this%zSigma_s = this%hgrid%zh
    ! this%zSigma_s = (this%hgrid%zz(1, :, 1)-this%hgrid%zz(1, this%kds_s, 1))*this%hgrid%zz(1, this%kde_s, 1)/(this%hgrid%zz(1, this%kde_s, 1)- this%hgrid%zz(1, this%kds_s, 1))
    PRINT *, "this%zSigma_s: ", this%zSigma_s

    ztop = this%zSigma_s(this%kde_s)

    this%zSigma_u = this%zSigma_s
    this%zSigma_u(this%kds_u + 1:this%kde_u - 1) = (this%zSigma_s(this%kds_s + 1:this%kde_s) + this%zSigma_s(this%kds_s:this%kde_s - 1)) / 2.0D0
    ! this%zSigma_u(this%kds_u) = 0 !-1.0*this%zSigma_u(this%kds_u + 1)
    ! this%zSigma_u(this%kde_u) = this%zSigma_s(this%kde_s) ! 2.0*this%zSigma_s(this%kde_s) - this%zSigma_u(this%kde_u - 1)
    this%zSigma_u(this%kds_u) = -1.0 * this%zSigma_u(this%kds_u + 1)
    this%zSigma_u(this%kde_u) = 2.0 * this%zSigma_s(this%kde_s) - this%zSigma_u(this%kde_u - 1)

    ! PRINT *, 'this%zSigma_u', this%zSigma_u
    ! PRINT *, 'this%zSigma_s', this%zSigma_s

    PRINT *, "read zRHght_s begin"
    ALLOCATE (this%zRHght_s(this%ids:this%ide, this%kds:this%kde, this%jds:this%jde), &
              this%zRHght_u(this%ids:this%ide, this%kds:this%kde, this%jds:this%jde))

    this%zRHght_s = this%hgrid%zz

    PRINT *, 'ztop is: ', ztop

    ! PRINT*, this%zRHght_s(:,2,:)
    this%zRHght_u(:, this%kds_u + 1:this%kde_u - 1, :) = (this%zRHght_s(:, this%kds_s + 1:this%kde_s, :) + &
                                                          this%zRHght_s(:, this%kds_s:this%kde_s - 1, :)) / 2.0D0
    this%zRHght_u(:, this%kds_u, :) = this%zSigma_u(this%kds_u) / ztop * (ztop - this%hgrid%ht) + this%hgrid%ht
    this%zRHght_u(:, this%kde_u, :) = this%zSigma_u(this%kde_u) / ztop * (ztop - this%hgrid%ht) + this%hgrid%ht

    PRINT *, 'Done here!'
  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(GrapesIO_t), INTENT(INOUT) :: this

    PRINT *, 'destructor works'

    IF (ALLOCATED(this%lat)) DEALLOCATE (this%lat)
    IF (ALLOCATED(this%lon)) DEALLOCATE (this%lon)

    IF (ALLOCATED(this%zSigma_s)) DEALLOCATE (this%zSigma_s)
    IF (ALLOCATED(this%zSigma_u)) DEALLOCATE (this%zSigma_u)
    IF (ALLOCATED(this%zRHght_s)) DEALLOCATE (this%zRHght_s)
    IF (ALLOCATED(this%zRHght_u)) DEALLOCATE (this%zRHght_u)

    CALL this%hgrid%dealloc_grid_memory()
  END SUBROUTINE destructor

  SUBROUTINE releaseHGrid(this)
    CLASS(GrapesIO_t) :: this

    IF (ALLOCATED(this%lat)) DEALLOCATE (this%lat)
    IF (ALLOCATED(this%lon)) DEALLOCATE (this%lon)

    IF (ALLOCATED(this%zSigma_s)) DEALLOCATE (this%zSigma_s)
    ! IF (ALLOCATED(this%zSigma_u)) DEALLOCATE (this%zSigma_u)
    IF (ALLOCATED(this%zRHght_s)) DEALLOCATE (this%zRHght_s)
    ! IF (ALLOCATED(this%zRHght_u)) DEALLOCATE (this%zRHght_u)
    CALL this%hgrid%dealloc_grid_memory()
  END SUBROUTINE

  SUBROUTINE uc2ua(this, uwnd)
    CLASS(GrapesIO_t) :: this
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE :: uwnd

    uwnd(2:this%hgrid%ide, :, :) = (uwnd(1:this%hgrid%ide - 1, :, :) + uwnd(2:this%hgrid%ide, :, :)) / 2.0D0
    uwnd(1, :, :) = 2 * uwnd(1, :, :) - uwnd(2, :, :)
  END SUBROUTINE uc2ua

  SUBROUTINE vc2va(this, vwnd)
    CLASS(GrapesIO_t) :: this
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE :: vwnd

    vwnd(:, :, 2:this%hgrid%jde) = (vwnd(:, :, 1:this%hgrid%jde - 1) + vwnd(:, :, 2:this%hgrid%jde)) / 2.0D0
    vwnd(:, :, 1) = 2 * vwnd(:, :, 1) - vwnd(:, :, 2)
  END SUBROUTINE vc2va

  SUBROUTINE ua2uc(this, uwnd)
    CLASS(GrapesIO_t) :: this
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE :: uwnd

    uwnd(1:this%hgrid%ide - 1, :, :) = &
      (uwnd(1:this%hgrid%ide - 1, :, :) + uwnd(2:this%hgrid%ide, :, :)) / 2.0D0

    uwnd(this%hgrid%ide, :, :) = 2 * uwnd(this%hgrid%ide, :, :) - uwnd(this%hgrid%ide - 1, :, :)
  END SUBROUTINE ua2uc

  SUBROUTINE va2vc(this, vwnd)
    CLASS(GrapesIO_t) :: this
    real(r_single), DIMENSION(:, :, :), ALLOCATABLE :: vwnd

    vwnd(:, :, 1:this%hgrid%jde - 1) = &
      (vwnd(:, :, 1:this%hgrid%jde - 1) + vwnd(:, :, 2:this%hgrid%jde))/2.0D0

    vwnd(:, :, this%hgrid%jde) = 2*vwnd(:, :, this%hgrid%jde) - vwnd(:, :, this%hgrid%jde - 1)
  END SUBROUTINE va2vc

  SUBROUTINE writeBackToQcqrInput(this, inputQcqrName, OutputQcqrName)
    CLASS(GrapesIO_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: inputQcqrName, OutputQcqrName
    INTEGER(i_kind) :: i, j, k, nx, ny, nz
    REAL(r_single), ALLOCATABLE :: oriSwap(:, :, :)
    TYPE(domain_t) :: hgrid_tmp

    hgrid_tmp = this%hgrid
    IF (ALLOCATED(hgrid_tmp%qc) .AND. ALLOCATED(hgrid_tmp%qr)) THEN
           
      OPEN(137,FILE=TRIM(inputQcqrName),FORM='unformatted',ACCESS='sequential',   &
             convert='big_endian',STATUS='unknown')

      OPEN(138,FILE=TRIM(OutputQcqrName),FORM='unformatted',ACCESS='sequential',   &
             convert='big_endian',STATUS='unknown')

      ALLOCATE(oriSwap(this%ids:this%ide, this%kds:this%kde, this%jds:this%jde))
      PRINT *, 'out grids: ', this%ids, this%ide, this%kds, this%kde, this%jds, this%jde
      ! gcas only contains qv, qc, qr   
      DO k = this%kds+1, this%kde
        READ(137) ((oriSwap(i,k,j),i=this%ids,this%ide),j=this%jds,this%jde) ! qv
      END DO
      hgrid_tmp%q = oriSwap + hgrid_tmp%q
      hgrid_tmp%q(:, this%kds, :) = 0.0
      WHERE(hgrid_tmp%q <=0.0D0) hgrid_tmp%q = 1E-7
      DO k = this%kds+1, this%kde
        WRITE(138) ((hgrid_tmp%q(i,k,j),i=this%ids,this%ide),j=this%jds,this%jde)
      END DO
      PRINT *, 'QV is re-written in qcqr-MOTORDA'

      DO k = this%kds+1, this%kde
        READ(137) ((oriSwap(i,k,j),i=this%ids,this%ide),j=this%jds,this%jde) ! qc
      END DO
      hgrid_tmp%qc = oriSwap + hgrid_tmp%qc
      hgrid_tmp%qc(:, this%kds, :) = 0.0
      WHERE(hgrid_tmp%qc <=0.0D0) hgrid_tmp%qc = 1E-7
      DO k = this%kds+1, this%kde
        WRITE(138) ((hgrid_tmp%qc(i,k,j),i=this%ids,this%ide),j=this%jds,this%jde)
      END DO
      PRINT *, 'QC is re-written in qcqr-MOTORDA'

      DO k = this%kds+1, this%kde
        READ(137) ((oriSwap(i,k,j),i=this%ids,this%ide),j=this%jds,this%jde) ! qr
      END DO
      hgrid_tmp%qr = oriSwap + hgrid_tmp%qr
      hgrid_tmp%qr(:, this%kds, :) = 0.0
      WHERE(hgrid_tmp%qr <=0.0D0) hgrid_tmp%qr = 1E-7
      DO k = this%kds+1, this%kde
        WRITE(138) ((hgrid_tmp%qr(i,k,j),i=this%ids,this%ide),j=this%jds,this%jde)
      END DO
      PRINT *, 'QR is re-written in qcqr-MOTORDA'

      DO k = this%kds+1, this%kde
        READ(137) ((oriSwap(i,k,j),i=this%ids,this%ide),j=this%jds,this%jde) ! qi
      END DO
      hgrid_tmp%qi = oriSwap + hgrid_tmp%qi
      hgrid_tmp%qi(:, this%kds, :) = 0.0
      WHERE(hgrid_tmp%qi <=0.0D0) hgrid_tmp%qi = 1E-7
      DO k = this%kds+1, this%kde
        WRITE(138) ((hgrid_tmp%qi(i,k,j),i=this%ids,this%ide),j=this%jds,this%jde)
      END DO
      PRINT *, 'QI is re-written in qcqr-MOTORDA'

      DO k = this%kds+1, this%kde
        READ(137) ((oriSwap(i,k,j),i=this%ids,this%ide),j=this%jds,this%jde) ! qs
      END DO
      hgrid_tmp%qs = oriSwap + hgrid_tmp%qs
      hgrid_tmp%qs(:, this%kds, :) = 0.0
      WHERE(hgrid_tmp%qs <=0.0D0) hgrid_tmp%qs = 1E-7
      DO k = this%kds+1, this%kde
        WRITE(138) ((hgrid_tmp%qs(i,k,j),i=this%ids,this%ide),j=this%jds,this%jde)
      END DO
      PRINT *, 'QS is re-written in qcqr-MOTORDA'

      DO k = this%kds+1, this%kde
        READ(137) ((oriSwap(i,k,j),i=this%ids,this%ide),j=this%jds,this%jde) ! qg
      END DO
      hgrid_tmp%qg = oriSwap + hgrid_tmp%qg
      hgrid_tmp%qg(:, this%kds, :) = 0.0
      WHERE(hgrid_tmp%qg <=0.0D0) hgrid_tmp%qg = 1E-7
      DO k = this%kds+1, this%kde
        WRITE(138) ((hgrid_tmp%qg(i,k,j),i=this%ids,this%ide),j=this%jds,this%jde)
      END DO
      PRINT *, 'QG is re-written in qcqr-MOTORDA'

      DEALLOCATE(oriSwap)
      CLOSE(137)
      CLOSE(138)

      PRINT *, 'qcqr is flushed to the file.'
    END IF

    IF (ALLOCATED(oriSwap)) DEALLOCATE (oriSwap)
  END SUBROUTINE writeBackToQcqrInput

  SUBROUTINE writeBackToGrapesInput(this)
    ! USE Filter_m, ONLY: guidedfilter

    CLASS(GrapesIO_t) :: this
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

END MODULE GrapesIO_m
