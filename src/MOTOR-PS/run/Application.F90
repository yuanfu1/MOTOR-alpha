MODULE Application_m
  USE rk4_m, ONLY: rk4_t
  USE GenContainers_m, ONLY: GenContainers_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State_m, ONLY: State_t
  USE IOGrapes_m, ONLY: IOGrapes_t
  USE YAMLRead_m
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind

  TYPE Application_t
    TYPE(rk4_t) :: rk4
    TYPE(GenContainers_t) :: GenContainers
    TYPE(mpddGlob_t) :: mpddGlob
    TYPE(geometry_t) :: GeometryHalf, GeometryFull
    TYPE(State_t) :: XHalf, XFull
    CHARACTER(LEN=1024) :: configFile
    INTEGER(i_kind) :: gLevel
    REAL(r_kind) :: dt
  CONTAINS
    PROCEDURE :: initial
    !   PROCEDURE :: PreData
    PROCEDURE :: Integration
    PROCEDURE :: destroy
  END TYPE

CONTAINS

  SUBROUTINE initial(this, configFile)
    IMPLICIT NONE
    CLASS(Application_t) :: this
    CHARACTER(LEN=*), INTENT(IN) :: configFile
    REAL(r_kind), ALLOCATABLE :: sigmaHalf(:)
    INTEGER(i_kind) :: istatus

    this%configFile = configFile
    this%GenContainers = GenContainers_t(TRIM(this%configFile))

    istatus = yaml_get_var(TRIM(this%configFile), 'TimeIntegral', 'gLevel', this%gLevel)
    IF (istatus .NE. 0) THEN
      PRINT *, 'ERROR: gLevel for the model simulation is not set'
      STOP
    END IF
    istatus = yaml_get_var(TRIM(this%configFile), 'TimeIntegral', 'time_step', this%dt)
    IF (istatus .NE. 0) THEN
      PRINT *, 'ERROR: time_step for the model simulation is not set'
      STOP
    END IF

    CALL this%GenContainers%GenGeometry(this%GeometryFull)
    ASSOCIATE (sigma => this%GeometryFull%mg%sg(this%GeometryFull%mg%mg_coarsest)%sigma, &
               vLevel => this%GeometryFull%mg%sg(this%GeometryFull%mg%mg_coarsest)%vLevel)
      ALLOCATE (sigmaHalf(vLevel - 1))
      sigmaHalf = (sigma(2:vLevel) + sigma(1:vLevel - 1)) / 2.0D0
      CALL this%GenContainers%GenGeometry(this%GeometryHalf, sigmaHalf)

         !! the sg%ztopo has not been given, thus the following two functions have problems for terrain following

      CALL this%GenContainers%updateGeometryAllSg(this%geometryFull, this%geometryHalf)
      CALL this%GenContainers%updateGeometryAllSg(this%geometryHalf, this%GeometryFull)

      CALL this%GenContainers%GenStates(this%geometryFull%mg%sg(this%gLevel), this%XFull)
      CALL this%GenContainers%GenStates(this%GeometryHalf%mg%sg(this%gLevel), this%XHalf)

    END ASSOCIATE

  END SUBROUTINE

!    SUBROUTINE PreData()
!    END SUBROUTINE

  SUBROUTINE Integration(this)
    IMPLICIT NONE
    CLASS(Application_t) :: this

    this%rk4 = rk4_t(this%configfile, this%glevel, this%GeometryHalf, this%GeometryFull)

    CALL this%rk4%timeIntegrals(this%XHalf, this%XFull)

  END SUBROUTINE

  SUBROUTINE destroy(this)
    IMPLICIT NONE
    CLASS(Application_t) :: this

    CALL this%GeometryFull%destroy
    CALL this%GeometryHalf%destroy
    CALL this%mpddGlob%barrier
    CALL this%GenContainers%destroy

    ! Finalize
    CALL this%mpddGlob%finalize

  END SUBROUTINE destroy

END MODULE
