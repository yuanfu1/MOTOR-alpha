MODULE GenContainers_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE GenSgDiffCoef_m, ONLY: GenSgDiffCoef_t
  USE InterpValue_m, ONLY: InterpValue_t
  USE YAMLRead_m
  USE parameters_m, ONLY: Omega
  USE GenBdy_m, ONLY: GenBdy_t

  TYPE :: GenContainers_t
    TYPE(mpddGlob_t) :: mpddGlob
    TYPE(GenBdy_t) :: GenBdy
    CHARACTER(LEN=1024) :: configFile
    REAL(r_kind) :: verifyLatlon

  CONTAINS
    PROCEDURE :: GenGeometry
    PROCEDURE :: GenStates
    PROCEDURE :: updateSg ! update single layer sg
    PROCEDURE :: updateGeometryAllSg ! update all the sg in the geometry
    PROCEDURE :: destroy
  END TYPE
  INTERFACE GenContainers_t
    PROCEDURE :: constructor
  END INTERFACE GenContainers_t

CONTAINS

  FUNCTION constructor(configFile) RESULT(this)
    IMPLICIT NONE
    TYPE(GenContainers_t) :: this
    CHARACTER(*), INTENT(IN) :: configFile
    this%configFile = TRIM(configFile)

    ! Yuanfu Xie temporarily turned this off as the main test program Test_PS.F90 has already initialized a MPDD on 2025-02-13
    !CALL this%mpddGlob%initialize()                                   ! Initialize the mpdd

    this%GenBdy = GenBdy_t(this%configFile)
    PRINT *, 'Finish initialization of mpddGlob: ', this%mpddGlob%myrank

    ! Reading the verification area latlon: Yuanfu Xie added on 2022-06-11
    IF (this%mpddGlob%isBaseProc()) THEN
      WRITE (*, 19)
19    FORMAT('#================================================================================#')
      WRITE (*, 18)
18    FORMAT('# If you have not set Verify parameters in your yaml file, this may crash!!! #', /, &
             '# If so, please set it and rerun!                                                #', /, &
             '# otherwise, you can ignore this message.                                        #')
      WRITE (*, 19)
    END IF

  END FUNCTION

  SUBROUTINE GenGeometry(this, geometry, sigma)
    IMPLICIT NONE
    CLASS(GenContainers_t) :: this
    TYPE(geometry_t), INTENT(INOUT) :: geometry
    REAL(r_kind), INTENT(IN), OPTIONAL :: sigma(:)
    REAL(r_kind), ALLOCATABLE :: sigma2(:)
    INTEGER(i_kind) :: i, ik, istatus

    ! Initialize geometry
    ! CALL geometry%initialize(this%configFile, this%mpddGlob)

    PRINT *, 'geometry is set!!!', this%mpddGlob%myrank,this%mpddGlob%rank

    DO i = geometry%mg%mg_coarsest, geometry%mg%mg_finest
      ASSOCIATE (sg => geometry%mg%sg(i))
        IF (PRESENT(sigma)) THEN
          DEALLOCATE (sg%sigma)
          ALLOCATE (sg%sigma(SIZE(sigma)))
          sg%sigma = sigma
          sg%vLevel = SIZE(sigma)
        ELSE
          istatus = yaml_get_var(TRIM(this%configFile), 'TimeIntegral', 'sigma', sigma2)
          IF (istatus .NE. 0) THEN
            PRINT *, 'warning: sigma is not set, use defualt'
            DEALLOCATE (sigma2)
            ALLOCATE (sigma2(sg%vLevel))
            DO ik = 1, sg%vLevel
              sigma2(ik) = (ik - 1.0D0) / (sg%vLevel - 1.0D0) * sg%ztop
            END DO
          END IF
          ! PRINT *, 'sigma2 is ', sigma2
          sg%sigma = sigma2
          DEALLOCATE (sigma2)
        END IF
      END ASSOCIATE
    END DO

    CALL this%updateGeometryAllSg(geometry, geometry)

  END SUBROUTINE

  SUBROUTINE updateGeometryAllSg(this, geometry_this, geometry_other)
    IMPLICIT NONE

    CLASS(GenContainers_t) :: this
    TYPE(geometry_t), INTENT(INOUT) :: geometry_this
    TYPE(geometry_t), INTENT(IN) :: geometry_other
    TYPE(GenSgDiffCoef_t) :: GenSgDiffCoef
    TYPE(InterpValue_t) :: InterpValue
    INTEGER(i_kind) :: i, k

    DO i = geometry_this%mg%mg_coarsest, geometry_this%mg%mg_finest

      CALL this%updateSg(geometry_this%mg%sg(i), geometry_other%mg%sg(i))

    END DO

  END SUBROUTINE

  SUBROUTINE updateSg(this, sg_this, sg_other)
    IMPLICIT NONE

    CLASS(GenContainers_t) :: this
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg_this
    TYPE(SingleGrid_t), INTENT(IN) :: sg_other
    TYPE(GenSgDiffCoef_t) :: GenSgDiffCoef
    TYPE(InterpValue_t) :: InterpValue
    INTEGER(i_kind) :: i, k

    CALL this%GenBdy%GenBdyIdx(sg_this)
    CALL this%GenBdy%GenBdyBufferWeight(sg_this)
    CALL GenSgDiffCoef%GenSgVerDiffCoefs(sg_other, sg_this)
    CALL InterpValue%UpdateSgInterpCoef(sg_other, sg_this)

    ALLOCATE (sg_this%f(sg_this%vLevel, sg_this%num_cell))

    DO k = 1, sg_this%vLevel
      sg_this%f(k, :) = 2.0D0 * Omega * DSIN(sg_this%cell_cntr(1, :))
    END DO

  END SUBROUTINE

  SUBROUTINE GenStates(this, sg, X)
    IMPLICIT NONE
    CLASS(GenContainers_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    TYPE(SingleGrid_t), INTENT(IN) :: sg

    CALL X%initialize(this%configFile, sg)

  END SUBROUTINE

  SUBROUTINE destroy(this)
    IMPLICIT NONE
    CLASS(GenContainers_t) :: this

    CALL this%mpddGlob%barrier

    ! Finalize
    CALL this%mpddGlob%finalize

  END SUBROUTINE destroy

END MODULE
