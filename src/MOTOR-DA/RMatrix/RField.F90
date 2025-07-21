!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.BMatrix.RField
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2022/3/10, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE RField_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE ObsField_m, ONLY: ObsField_t
  USE AuxTypeObs_m, ONLY: AuxTypeObs_t
  USE MPObs_m, ONLY: MPObs_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE MultiGrid_m, ONLY: MultiGrid_t
  USE YAMLRead_m

  TYPE, EXTENDS(AuxTypeObs_t):: RField_t
    REAL(r_kind), ALLOCATABLE :: sigma(:)
    REAL(r_kind) :: scaleObs = 1.0, ScaleTBB = 1.0D0
    CHARACTER(len=20) :: name
    CHARACTER(LEN=1024) :: configFile
    LOGICAL :: isEye
    INTEGER :: numObs
    REAL(r_kind), ALLOCATABLE :: z_varied_sigma(:)

  CONTAINS
    PROCEDURE, PUBLIC :: initialize
    PROCEDURE, PRIVATE :: loadBMatFiles
    PROCEDURE, PUBLIC :: sqrt_inverse_multiply
    PROCEDURE, PUBLIC :: sqrt_inverse_multiply_adjoint
    PROCEDURE, PUBLIC :: Get_Name
    FINAL :: destructor
  END TYPE RField_t

CONTAINS

  SUBROUTINE initialize(this, configFile, mpObs, Y, sg, Nocoarest)
    IMPLICIT NONE
    CLASS(RField_t) :: this
    TYPE(ObsField_t), INTENT(IN) :: Y
    TYPE(MPObs_t), TARGET, INTENT(IN) :: mpObs
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    INTEGER(i_kind), OPTIONAL, INTENT(IN) :: Nocoarest

    INTEGER :: vidx, hidx, iobs, Ratio, istatus
    REAL :: RatioErr

    this%AuxTypeObs_t = AuxTypeObs_t(mpObs)
    this%name = Y%Get_Id_Name()
    this%configFile = configFile

    ALLOCATE (this%z_varied_sigma(SIZE(Y%values)))
    this%z_varied_sigma = 0.001D0

    istatus = yaml_get_var(TRIM(this%configFile), 'BMat', 'ScaleTBB', this%ScaleTBB)

    FORALL (iobs=1:SIZE(Y%values)) this%z_varied_sigma(iobs) = sg%SVapor(Y%idx(iobs)%vidx, Y%idx(iobs)%hidx)
    ! DO iobs = 1, size(Y%values)
    !   this%z_varied_sigma(iobs) = max(0.1D-5, this%z_varied_sigma(iobs))
    ! END DO
    ! DO iobs = 1, size(Y%values)
    !   IF (TRIM(this%name) .NE. 'SYNOP_qvapor' .AND. TRIM(this%name) .NE. 'SOUND_qvapor') CYCLE
    !     ! IF (Y%values(iobs) > 0.0D0) THEN))
    !     !   Ratio = nint(log10(Y%values(iobs)))
    !     !   RatioErr = 10.0D0**(Ratio)
    !     !   this%z_varied_sigma(iobs) = RatioErr
    !     !   this%z_varied_sigma(iobs) = min(0.001D0,this%z_varied_sigma(iobs))
    !     !   this%z_varied_sigma(iobs) = max(0.1D-3, this%z_varied_sigma(iobs))
    !     !   IF (Y%idx(iobs)%vidx >10) PRINT *, 'check z_varied_sigma: ', this%name, Y%idx(iobs)%vidx, Y%values(iobs), this%z_varied_sigma(iobs)
    !     ! END IF

    !   DO vidx = 1, sg%vLevel
    !     IF (Y%idx(iobs)%vidx .EQ. vidx .AND. Y%values(iobs) > 0.0D0) THEN
    !       this%z_varied_sigma(iobs) = sg%SVapor1D(vidx)
    !       IF (Y%idx(iobs)%vidx >10) PRINT *, 'check z_varied_sigma: ', this%name, Y%idx(iobs)%vidx, Y%values(iobs), this%z_varied_sigma(iobs)
    !       CYCLE
    !     END IF
    !   END DO
    ! END DO

    ! Here is the code of monk for simulating the fileinput of B mat files
    ALLOCATE (this%sigma(SIZE(Y%values))); 
    ! For normalization
    ! IF (PRESENT(Nocoarest)) THEN
    !   IF (Nocoarest .GT. 0) THEN
    !     PRINT *, 'Check normalization: ', ' this%name = ', TRIM(this%name), ' No = ', this%numObs
    !     this%sigma = this%sigma * SQRT(DBLE(this%numObs) / DBLE(Nocoarest))
    !   END IF
    ! END IF

    IF (PRESENT(Nocoarest)) THEN
      IF (Nocoarest .GT. 0) THEN
        CALL Y%mpObs%AllReduceSumInt(SIZE(Y%values), this%numObs)
        PRINT *, 'Check normalization: ', ' this%name = ', TRIM(this%name), ' No = ', this%numObs
        this%sigma = this%sigma * SQRT(DBLE(this%numObs) / DBLE(Nocoarest))
      END IF
    END IF
    ! this%sigma = this%sigma * (DBLE(this%numObs) / sg%num_icell /sg%tSlots/ sg%vLevel)

    ! IF (TRIM(this%name) .EQ. 'SOUND_qvapor' .AND. sg%gLevel >5 ) this%sigma = this%sigma * 0.1D0
    ! For scaling
    CALL Y%mpObs%AllReduceSumInt(SIZE(Y%values), this%numObs)

    ! ---------------------------------------- For testing ----------------------------------------
    this%scaleObs = DBLE(sg%num_icell_global * sg%vLevel * sg%tSlots) / DBLE(this%numObs)
    IF (this%numObs > 1) THEN
      SELECT CASE (TRIM(this%name(1:5)))
      CASE ('SYNOP')
        this%scaleObs = DBLE(sg%num_icell_global * sg%vLevel * sg%tSlots) / DBLE(this%numObs)
        PRINT *, 'NUM SYNOP: ', this%numObs, this%scaleObs
      CASE ('SOUND')
        this%scaleObs = DBLE(sg%num_icell_global * sg%vLevel * sg%tSlots) / DBLE(this%numObs) / 5.0D0
      CASE ('RAD_Z')
        this%scaleObs = DBLE(sg%num_icell_global * sg%vLevel * sg%tSlots) / DBLE(this%numObs) / 700.0D0
        FORALL (iobs=1:SIZE(Y%values))
          this%z_varied_sigma(iobs) = 2.0 - 1.0D0 / ((1.0D0 + EXP(-0.0012D0 * (sg%sigma(Y%idx(iobs)%vidx) - 3500.0D0)))) * &
                                      (1.0D0 - (1.0D0 / ((1.0D0 + EXP(-0.0012D0 * (sg%sigma(Y%idx(iobs)%vidx) - 17000.0D0))))))
        END FORALL
      CASE ('fy4_1')
        this%scaleObs = DBLE(sg%num_icell_global * sg%vLevel * sg%tSlots) / DBLE(this%numObs) / this%ScaleTBB ! 200.0D0
        PRINT *, 'NUM fy4_1: ', this%numObs, this%scaleObs
      CASE ('fy4_2')
        this%scaleObs = DBLE(sg%num_icell_global * sg%vLevel * sg%tSlots) / DBLE(this%numObs) / this%ScaleTBB ! 200.0D0
        ! PRINT*, this%scaleObs
      CASE ('VWPW_')
        this%scaleObs = DBLE(sg%num_icell_global * sg%vLevel * sg%tSlots) / DBLE(this%numObs) / 500.0D0
      CASE ('GNSSR')
        this%scaleObs = DBLE(sg%num_icell_global * sg%vLevel * sg%tSlots) / DBLE(this%numObs) / 30.0D0
        ! 5.0 => scaleobs = ~10
        ! STOP
      END SELECT

      SELECT TYPE (mg => sg%mgParent)
      TYPE is (MultiGrid_t)
        PRINT *, 'MultiGrid_t', mg%mg_coarsest, mg%mg_finest
        SELECT CASE (TRIM(this%name))
        CASE ('SOUND_temp')
          this%scaleObs = this%scaleObs / (1.7**(sg%gLevel - mg%mg_coarsest))
        CASE ('VWPW_uwnd')
          this%scaleObs = this%scaleObs * (2.3**(sg%gLevel - mg%mg_coarsest))
        CASE ('VWPW_vwnd')
          this%scaleObs = this%scaleObs * (2.3**(sg%gLevel - mg%mg_coarsest))
        CASE ('SOUND_uwnd')
          FORALL (iobs=1:SIZE(Y%values))
            ! ff = -2./((1+exp(-0.001*(h-2000)))).*(1-(1./((1+exp(-0.001*(h-24000))))))+4;
            this%z_varied_sigma(iobs) = 4.0 - 2.0D0 / ((1.0D0 + EXP(-0.0012D0 * (sg%sigma(Y%idx(iobs)%vidx) - 3500.0D0)))) * &
                                        (1.0D0 - (1.0D0 / ((1.0D0 + EXP(-0.0012D0 * (sg%sigma(Y%idx(iobs)%vidx) - 17000.0D0))))))
          END FORALL
          PRINT *, MAXVAL(this%z_varied_sigma), MINVAL(this%z_varied_sigma)
        CASE ('SOUND_vwnd')
          FORALL (iobs=1:SIZE(Y%values))
            this%z_varied_sigma(iobs) = 4.0 - 2.0D0 / ((1.0D0 + EXP(-0.0012D0 * (sg%sigma(Y%idx(iobs)%vidx) - 3500.0D0)))) * &
                                        (1.0D0 - (1.0D0 / ((1.0D0 + EXP(-0.0012D0 * (sg%sigma(Y%idx(iobs)%vidx) - 17000.0D0))))))
          END FORALL
        END SELECT

      END SELECT
    END IF

    CALL this%loadBMatFiles(Y)

    ! PRINT*, this%scaleObs
    ! STOP
  END SUBROUTINE

  SUBROUTINE loadBMatFiles(this, Y)
    CLASS(RField_t) :: this
    TYPE(ObsField_t), INTENT(IN) :: Y

    ! Set all the sigma to 1.0D0
    this%sigma = 1.0D0
    this%isEye = .FALSE.
    ! this%sigma = 0.01D0
    ! this%isEye = .FALSE.
    ! print *, 'AA:', this%name

    IF (TRIM(this%name) .EQ. 'RADAR_ref') THEN
      this%sigma = 1.0D0
      this%isEye = .FALSE.
    END IF

    IF (INDEX(TRIM(this%name), 'tbb') .NE. 0) THEN
      this%isEye = .FALSE.
      ! this%sigma = 1.0D0
      this%sigma = Y%errors !TODO: temporarily comment out
      PRINT *, 'check tbb sigma ', MAXVAL(this%sigma), MINVAL(this%sigma)
      ! this%sigma = this%sigma * 14.0
    END IF

    IF (this%name(1:4) == 'RAD_' .AND. this%name(11:14) == 'rwnd') THEN
      this%isEye = .FALSE.
      this%sigma = this%z_varied_sigma
    END IF

    SELECT CASE (TRIM(this%name))
    CASE ('SYNOP_pres')
      this%isEye = .FALSE.
      this%sigma = 100.0D0
    CASE ('SYNOP_psl')
      this%isEye = .FALSE.
      this%sigma = 100.0D0
    CASE ('SYNOP_qvapor')
      this%isEye = .FALSE.
      this%sigma = this%z_varied_sigma * 1.0D0 !0.5D0, 0.8D0, 0.7D0, 0.5D0, 0.6D0, 0.4D0
      ! lower, weaker
      ! this%sigma = 0.001D0
      ! CASE ('RADAR_vel')
      !   this%isEye = .FALSE.
      !   this%sigma = 50.0D0
    CASE ('SOUND_temp')
      this%isEye = .FALSE.
      this%sigma = 1.00D0
    CASE ('SOUND_qvapor')
      this%isEye = .FALSE.
      ! this%sigma = 0.001D0
      this%sigma = this%z_varied_sigma * 0.25D0! * 2.0D0 !0.5D0; 1.5D0; 1.0D0; 0.8D0; 0.6D0; 0.5D0; 0.4D0
      ! higher, stronger
      !print *, 'check sigma of SOUND_qvapor', this%sigma
    CASE ('SOUND_uwnd')
      this%isEye = .FALSE.
      this%sigma = this%z_varied_sigma !4.0D0
    CASE ('SOUND_vwnd')
      this%isEye = .FALSE.
      this%sigma = this%z_varied_sigma !4.0D0
    CASE ('VWPW_uwnd')
      this%isEye = .FALSE.
      this%sigma = 2.0D0
    CASE ('VWPW_vwnd')
      this%isEye = .FALSE.
      this%sigma = 2.0D0
    CASE ('SYNOP_uwnd')
      this%isEye = .FALSE.
      this%sigma = 1.0D0
    CASE ('SYNOP_vwnd')
      this%isEye = .FALSE.
      this%sigma = 1.0D0
    CASE ('SYNOP_lnp')
      this%isEye = .FALSE.
      this%sigma = 1.00D0
    CASE ('SYNOP_temp')
      this%isEye = .FALSE.
      this%sigma = 1.0D0
      ! CASE ('SYNOP_qvapor')
      !   this%isEye = .FALSE.
      !   this%sigma = 2.0D0
    CASE ('GNSSRO_refractivity')
      this%isEye = .FALSE.
      this%sigma = this%z_varied_sigma * 0.25D0

    CASE ('RAD_Sweep_102372_rwn')
      this%isEye = .FALSE.
      this%sigma = 0.1
    END SELECT

    ! IF (TRIM(this%name) .EQ. 'RADAR_ref') THEN
    !   this%sigma = 1.0D0
    !   this%isEye = .FALSE.
    ! END IF
  END SUBROUTINE loadBMatFiles

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(RField_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%sigma)) DEALLOCATE (this%sigma)
    IF (ALLOCATED(this%z_varied_sigma)) DEALLOCATE (this%z_varied_sigma)
  END SUBROUTINE destructor

  SUBROUTINE sqrt_inverse_multiply(this, Y)
    IMPLICIT NONE
    CLASS(RField_t) :: this
    TYPE(ObsField_t), INTENT(INOUT) :: Y

    IF (this%isEye) RETURN
    !Y%values = Y%values/(this%sigma*this%z_varied_sigma)
    Y%values = Y%values / this%sigma * this%scaleObs

  END SUBROUTINE sqrt_inverse_multiply

  SUBROUTINE sqrt_inverse_multiply_adjoint(this, Y)
    IMPLICIT NONE
    CLASS(RField_t) :: this
    TYPE(ObsField_t), INTENT(INOUT) :: Y

    IF (this%isEye) RETURN
    !Y%values = Y%values/(this%sigma*this%z_varied_sigma)
    Y%values = Y%values / this%sigma * this%scaleObs
  END SUBROUTINE sqrt_inverse_multiply_adjoint

  !> @brief
!! Return the type of this field.
  FUNCTION Get_Name(this) RESULT(name)
    CLASS(RField_t) :: this
    CHARACTER(LEN=20) :: name

    name = this%name
  END FUNCTION
END MODULE RField_m
