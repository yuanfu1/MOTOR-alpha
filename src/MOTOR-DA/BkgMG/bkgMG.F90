!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA/bkgMG
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2024-01-25   Created by Yuanfu Xie
!
!   Created by Yuanfu Xie (xieyf@gbamwf.com), 2024/01/25, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides an interface of accessing background fields for multigrid.
MODULE bkgMG_m
  USE kinds_m
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State_m, ONLY: State_t
  USE YAMLRead_m
  USE IOModel_m, ONLY: IOModel_t
  USE IOWRF_m, ONLY: IOWRF_t
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE State2NC_m

  USE MGOpts_m
  USE Ctl2State_m, ONLY: Ctl2State_t

  TYPE :: bkgMG_t
    CHARACTER(LEN=8) :: bkgdModel
    TYPE(geometry_t), POINTER :: geo
    CLASS(IOModel_t), ALLOCATABLE :: IOModel
    TYPE(Ctl2State_t) :: Ctl2State

    ! Backgrounds:
    TYPE(State_t), ALLOCATABLE :: XbMG(:)
  CONTAINS
    PROCEDURE :: initialize_s
    PROCEDURE :: writeBkgd_s
    PROCEDURE :: destroy_s
    PROCEDURE :: getBackgrd_s

  END TYPE bkgMG_t

  ABSTRACT INTERFACE
    SUBROUTINE model(this, yamlFile)
      IMPORT :: bkgMG_t, state_t, geometry_t
      CLASS(bkgMG_t) :: this
      CHARACTER(LEN=1024), INTENT(IN) :: yamlFile
    END SUBROUTINE model
  END INTERFACE

CONTAINS
  !> @brief
  !! This routine initializes background module:
  SUBROUTINE initialize_s(this, geometry, yamlFile)
    IMPLICIT NONE
    CLASS(bkgMG_t) :: this
    TYPE(geometry_t), TARGET, INTENT(IN) :: geometry
    CHARACTER(LEN=1024), INTENT(IN) :: yamlFile

    ! Local variables:
    INTEGER(i_kind) :: istatus

    this%geo => geometry

    ! Select a background model:
    istatus = yaml_get_var(TRIM(yamlFile), 'IO', 'bk_model', this%bkgdModel)

    ! Initialize the geometry
    ALLOCATE (this%XbMG(this%geo%mg%mg_coarsest:this%geo%mg%mg_finest))

    ! Initialize the Scaling implementation
    CALL this%Ctl2State%initialize(yamlFile)
  END SUBROUTINE initialize_s

  !> @brief
        !! This routine writes out background fields to NC files:
  SUBROUTINE writeBkgd_s(this, yamlFile)
    IMPLICIT NONE
    CLASS(bkgMG_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: yamlFile

    ! Local variables:
    CHARACTER(LEN=1024) :: ncOutputFile
    CHARACTER(LEN=20) :: task
    INTEGER(i_kind) :: i, istatus

    ! Find out the directory to write:
    istatus = yaml_get_var(TRIM(yamlFile), 'IO', 'output_dir', ncOutputFile)
    istatus = yaml_get_var(TRIM(yamlFile), 'RunMode', 'Task', task)

    ! Write XbMG to NC files:
    DO i = this%geo%mg%mg_coarsest, this%geo%mg%mg_finest
      CALL Output_NC_State_AV(this%XbMG(i), ncOutputFile, &
                              TRIM(task)//"_bkg", .TRUE., .TRUE.)
    END DO
  END SUBROUTINE writeBkgd_s

  SUBROUTINE destroy_s(this)
    CLASS(bkgMG_t) :: this

    DEALLOCATE (this%XbMG)
  END SUBROUTINE destroy_s

  SUBROUTINE getBackgrd_s(this, yamlFile)
    IMPLICIT NONE
    CLASS(bkgMG_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: yamlFile

    ! Local variables:
    INTEGER(i_kind) :: mgStart, mgEnd, i
    TYPE(IOWRF_t), TARGET :: ioWRF
    TYPE(IOGrapes_t), TARGET :: ioGrapes

    mgStart = this%geo%mg%mg_coarsest ! Output array starting position, corresponding to this%geo%mg%mg_coarsest
    mgEnd = this%geo%mg%mg_finest

    ASSOCIATE (sg => this%geo%mg%sg(this%geo%mg%mg_finest))
      CALL this%XbMG(mgEnd)%initialize(yamlFile, sg)
      CALL this%XbMG(mgEnd)%ShowInfo

      BLOCK
        IF (TRIM(this%bkgdModel) == 'WRF') THEN
          ALLOCATE (IOWRF_t:: this%IOModel)

          ASSOCIATE (IOModel => this%IOModel)
            SELECT TYPE (IOModel)
            TYPE IS (IOWRF_t)
              CALL IOModel%initialize(yamlFile, this%geo)
            END SELECT
          END ASSOCIATE
          CALL this%IOModel%m_read_bcg_into_Xm_Ens(this%XbMG(mgEnd), sg)
          PRINT *, TRIM(this%bkgdModel), ' reading is done!'

        ELSE IF (TRIM(this%bkgdModel) == 'GRAPES') THEN
          ALLOCATE (IOGrapes_t:: this%IOModel)

          ASSOCIATE (IOModel => this%IOModel)
            SELECT TYPE (IOModel)
            TYPE IS (IOGrapes_t)
              CALL IOModel%initialize(yamlFile, this%geo)
            END SELECT
          END ASSOCIATE

          PRINT *, 'READING finest grid level from ', TRIM(this%bkgdModel)
          CALL this%IOModel%m_read_bcg_into_Xm(this%XbMG(mgEnd), sg)
          PRINT *, TRIM(this%bkgdModel), ' reading is done!'
        END IF
      END BLOCK
    END ASSOCIATE

    ! Restrict the finer resolution background to coarser grid:
    DO i = mgEnd - 1, mgStart, -1
      ASSOCIATE (sgFiner => this%geo%mg%sg(i + 1), sgCoarser => this%geo%mg%sg(i))
        CALL this%XbMG(i)%initialize(yamlFile, sgCoarser)
        CALL restrictionMG(this%XbMG(i), this%XbMG(i + 1), this%geo%mg)
        CALL this%geo%mg%restrictionOfStatics(sgFiner, sgCoarser)
        CALL update_restrictionOfStatics(yamlFile, this%XbMG(i), sgCoarser)
      END ASSOCIATE
    END DO

    ! Temporarily follow the operation run setting:
    ! DO i = mgEnd, mgStart, -1
    !   CALL this%Ctl2State%transBackward(this%XbMG(i))
    ! END DO
  END SUBROUTINE getBackgrd_s
END MODULE bkgMG_m
