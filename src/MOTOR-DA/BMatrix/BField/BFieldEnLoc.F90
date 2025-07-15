!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.BMatrix.BField
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2023/3/1, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE BFieldEnLoc_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE AuxTypeSG_m, ONLY: AuxTypeSG_t
  USE BFieldBase_m, ONLY: BFieldBase_t
  USE BKErr_m, ONLY: BKErr_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE YAMLRead_m
  !USE onemkl_blas_omp_offload_lp64

  TYPE, EXTENDS(BFieldBase_t) :: BFieldEnLoc_t
    TYPE(ObsSet_t), POINTER    :: Y
    TYPE(Field_t), ALLOCATABLE :: EnLoc_ADJ(:)
    INTEGER(i_kind) :: nc          ! the number of control variables
    INTEGER(i_kind) :: nz          ! the number of non-zero eigenvalues
    INTEGER(i_kind) :: numsgicell  ! vlel*num_icell*tslots in each sg   
  CONTAINS
    PROCEDURE, PRIVATE :: loadBMatFiles

    PROCEDURE, PUBLIC  :: sqrt_inverse_multiply
    PROCEDURE, PUBLIC  :: sqrt_inverse_multiply_adjoint

    PROCEDURE :: sqrtInvMul => sqrt_inverse_multiply
    PROCEDURE :: sqrtInvMulAdj => sqrt_inverse_multiply_adjoint

    PROCEDURE, PUBLIC :: inverse_multiply

    PROCEDURE, PUBLIC :: initialize
    FINAL :: destructor
  END TYPE BFieldEnLoc_t

CONTAINS

  SUBROUTINE initialize(this, configFile, sg, varName, Y)
    IMPLICIT NONE

    CLASS(BFieldEnLoc_t) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    CHARACTER(LEN=1024), INTENT(IN)        :: configFile
    CHARACTER(LEN=10), INTENT(IN)          :: varName
    TYPE(ObsSet_t), TARGET, OPTIONAL       :: Y
    INTEGER(i_kind)                        :: ifile, numCell

    CALL this%AuxTypeSG_t%aux_initialize(sg)
    this%name = varName

    this%scaleParaX = 1.0D0
    this%scaleParaY = 1.0D0
    this%scaleParaZ = 1.0D0
    this%scaleParaT = 1.0D0

    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaX', this%scaleParaX)
    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaY', this%scaleParaY)
    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaZ', this%scaleParaZ)
    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaT', this%scaleParaT)

    IF (this%sg%isBaseProc()) THEN
       this%ts_DataLoad = 0.0D0
       this%ts_DataAggr = 0.0D0
       this%ts_DataDist = 0.0D0
       this%ts_DataCalc = 0.0D0
       this%ts_DataCalc_FWD = 0.0D0
       this%ts_DataCalc_ADJ = 0.0D0
    END IF


    !CALL this%loadBMatFiles(TRIM(varName))
    CALL this%loadBMatFiles(varName)

    IF (PRESENT(Y)) THEN
      this%Y => Y
    END IF

  END SUBROUTINE

  SUBROUTINE loadBMatFiles(this,  varName)
    IMPLICIT NONE
    CLASS(BFieldEnLoc_t) :: this
    CHARACTER(LEN=10), INTENT(IN)   :: varName
    TYPE(BKErr_t)                   :: br
    CHARACTER(LEN=1024)             :: outputDir, datFileName
    CHARACTER(LEN=3)                :: gid
    INTEGER(i_kind)                 :: status, nlvl, nhor, ntim, nz, nc, i, j, info, ifile
    INTEGER(i_kind)                 :: count_start, count_end, count_rate
    REAL(r_kind)                    :: ts
    REAL(r_kind), ALLOCATABLE       :: EnLoc_ADJ_tmp(:,:,:,:)

    ALLOCATE (this%sigma(this%sg%vLevel, this%sg%num_cell, this%sg%tSlots))
    this%sigma = 1.0D0 !0.000000001D0!1.0D0 !0.00001D0 !100.0D0 !10000.0D0

    ! read ADJ matrix file
    CALL GET_ENVIRONMENT_VARIABLE("OUTPUT_DIR", outputDir)
    WRITE (gid, '(A1,I2.2)') "G", this%sg%gLevel
    datFileName = TRIM(outputDir)//'/BEC_EnLoc_ADJ_'//TRIM(this%name)//'_'//gid
    CALL br%b_adjoint_read_dims(TRIM(datFileName))

    nlvl = this%sg%vLevel
    nhor = this%sg%num_icell_global
    ntim = this%sg%tSlots
    nz   = br%nz
    nc   = br%nc

    this%nc = nc
    this%nz = nz

    this%numsgicell = this%sg%vlevel*this%sg%num_icell*this%sg%tSlots

    IF (.NOT. ALLOCATED(this%EnLoc_ADJ)) ALLOCATE (this%EnLoc_ADJ(nz))

    IF (this%sg%isBaseProc()) THEN
      CALL system_clock(count_start, count_rate)
      CALL br%b_adjoint_read_data(TRIM(datFileName))
      CALL system_clock(count_end)

      ALLOCATE(EnLoc_ADJ_tmp(nz, nlvl, nhor, ntim))
      DO i = 1, nz
         EnLoc_ADJ_tmp(i,:,:,:) = RESHAPE(br%QRinvTI(:,i), (/nlvl, nhor, ntim/))
      END DO

      ts = REAL(count_end - count_start) / REAL(count_rate)       
      this%ts_DataLoad = this%ts_DataLoad + ts
    END IF

    DO j = 1, nz
      ALLOCATE(this%EnLoc_ADJ(j)%data(this%sg%vlevel, this%sg%num_cell, this%sg%tSlots))
      CALL this%sg%DistGridRealWithHaloExForFieldGrid(EnLoc_ADJ_tmp(j,:,:,:), this%EnLoc_ADJ(j)%data, &
              [this%sg%vLevel, this%sg%num_icell_global, this%sg%tSlots])
    END DO

    IF (this%sg%isBaseProc()) THEN
      IF (ALLOCATED(EnLoc_ADJ_tmp)) DEALLOCATE(EnLoc_ADJ_tmp)
    END IF

    PRINT *, 'DONE BField%loadBMatFiles ...'

  END SUBROUTINE loadBMatFiles

  SUBROUTINE inverse_multiply(this, field)
    IMPLICIT NONE
    CLASS(BFieldEnLoc_t) :: this
    TYPE(Field_t), INTENT(INOUT) :: field

  END SUBROUTINE inverse_multiply

  SUBROUTINE sqrt_inverse_multiply(this, field)
    IMPLICIT NONE
    CLASS(BFieldEnLoc_t) :: this
    TYPE(Field_t), INTENT(INOUT) :: field
    REAL(r_kind), ALLOCATABLE    :: datatmp(:, :, :)
    INTEGER(i_kind)              :: n, i, j, jj
    INTEGER(i_kind)              :: count_start, count_end, count_rate
    REAL(r_kind)                 :: ts
    REAL(r_kind), ALLOCATABLE    :: tmp_SUM(:)
    REAL(r_kind)                 :: tmp
    REAL(r_kind), ALLOCATABLE    :: w_tmp(:)

    !field%data = field%data * this%sigma

    IF (this%sg%isBaseProc()) THEN
      ALLOCATE(datatmp(this%sg%vLevel, this%sg%num_icell_global, this%sg%tSlots))
      datatmp = 0.0D0
    END IF

    ALLOCATE(w_tmp(this%numsgicell), tmp_SUM(this%nz))

    w_tmp = RESHAPE(field%data(:,1:this%sg%num_icell,:), (/this%numsgicell/))
    field%data = 0.0D0
    CALL this%sg%ExchangeMatOnHaloForFieldGrid(this%sg%tSlots, this%sg%vLevel, field%data)

    CALL system_clock(count_start, count_rate)
    DO i = 1, this%nz
      tmp = DOT_PRODUCT( RESHAPE(this%EnLoc_ADJ(i)%data(:,1:this%sg%num_icell,:),   &
                                  (/this%numsgicell/)),                             &
                         w_tmp )
      CALL this%sg%mpddInfo_sg%AllReduceSumReal(tmp, tmp_SUM(i))
    END DO
    CALL system_clock(count_end)

    IF (this%sg%isBaseProc()) THEN
      datatmp(1:this%nz,1,1) = tmp_SUM(:)

      ts = REAL(count_end - count_start) / REAL(count_rate)
      this%ts_DataCalc     = this%ts_DataCalc + ts
      this%ts_DataCalc_FWD = this%ts_DataCalc_FWD + ts
    END IF

    CALL system_clock(count_start, count_rate)
    CALL this%sg%DistGridRealWithHaloExForFieldGrid(datatmp, field%data, [this%sg%vLevel, this%sg%num_icell_global, this%sg%tSlots])
    CALL system_clock(count_end)

    IF (ALLOCATED(w_tmp))   DEALLOCATE(w_tmp)
    IF (ALLOCATED(tmp_SUM)) DEALLOCATE(tmp_SUM)

    IF (this%sg%isBaseProc()) THEN
      ts = REAL(count_end - count_start) / REAL(count_rate)
      this%ts_DataDist = this%ts_DataDist + ts

      IF (ALLOCATED(datatmp)) DEALLOCATE(datatmp)
    END IF

  END SUBROUTINE sqrt_inverse_multiply

  SUBROUTINE sqrt_inverse_multiply_adjoint(this, field, currentField)
    IMPLICIT NONE
    CLASS(BFieldEnLoc_t) :: this
    TYPE(Field_t), OPTIONAL, INTENT(IN) :: currentField
    TYPE(Field_t), INTENT(INOUT)        :: field
    REAL(r_kind), ALLOCATABLE           :: dataSwap(:, :, :)
    REAL(r_kind), ALLOCATABLE           :: datatmp(:), delta_x(:)
    INTEGER(i_kind)                     :: count_start, count_end, count_rate, i, j
    REAL(r_kind)                        :: ts
    REAL(r_kind), ALLOCATABLE           :: FWD_tmp(:, :)

    !field%data = field%data * this%sigma

    IF (this%sg%isBaseProc()) THEN
      ALLOCATE (dataSwap(this%sg%vLevel, this%sg%num_icell_global, this%sg%tSlots))
    END IF
    ALLOCATE (delta_x(this%nz))

    CALL system_clock(count_start, count_rate)
    CALL this%sg%aggrGridRealForFieldGrid(field%data, dataSwap, [this%sg%vLevel, this%sg%num_icell_global, this%sg%tSlots])
    CALL system_clock(count_end)
    IF (this%sg%isBaseProc()) THEN
      ts = REAL(count_end - count_start) / REAL(count_rate)
      this%ts_DataAggr = this%ts_DataAggr + ts
      delta_x(:) = dataSwap(1:this%nz, 1, 1)
    END IF

    CALL system_clock(count_start, count_rate)
    CALL this%sg%mpddInfo_sg%bCastReal1D(delta_x)

    ALLOCATE (FWD_tmp(this%numsgicell, this%nz))
    FWD_tmp = 0.0D0
    DO i = 1, this%nz
      FWD_tmp(:,i) = RESHAPE(this%EnLoc_ADJ(i)%data(:,1:this%sg%num_icell,:), (/this%numsgicell/))
    END DO
   
    ALLOCATE (datatmp(this%numsgicell))

    CALL DGEMV('N', this%numsgicell, this%nz,   &
         1.0d0, FWD_tmp, this%numsgicell,       &    
         delta_x, 1,                         &                     
         0.0d0, datatmp, 1)

    field%data(:,1:this%sg%num_icell,:) = RESHAPE(datatmp, (/this%sg%vLevel, this%sg%num_icell, this%sg%tSlots/))
    CALL this%sg%ExchangeMatOnHaloForFieldGrid(this%sg%tSlots, this%sg%vLevel, field%data)
    CALL system_clock(count_end)

    IF (this%sg%isBaseProc()) THEN
      ts = REAL(count_end - count_start) / REAL(count_rate)
      this%ts_DataCalc = this%ts_DataCalc + ts
      this%ts_DataCalc_ADJ = this%ts_DataCalc_ADJ + ts

      IF (ALLOCATED(dataSwap)) DEALLOCATE(dataSwap)
    END IF
    IF (ALLOCATED(delta_x)) DEALLOCATE(delta_x)
    IF (ALLOCATED(datatmp)) DEALLOCATE(datatmp)
    IF (ALLOCATED(FWD_tmp)) DEALLOCATE(FWD_tmp)


  END SUBROUTINE sqrt_inverse_multiply_adjoint

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(BFieldEnLoc_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%sigma)) DEALLOCATE (this%sigma)
    IF (ALLOCATED(this%EnLoc_ADJ)) DEALLOCATE (this%EnLoc_ADJ)
  END SUBROUTINE destructor

END MODULE BFieldEnLoc_m
