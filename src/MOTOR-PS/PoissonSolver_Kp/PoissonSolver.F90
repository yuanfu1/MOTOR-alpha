!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.poissonSolver_Kp
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie, Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module is the container of 2D poisson solver with multigrid method. \n
!! Created by Yuanfu Xie, in grids/gzGrid.F90 \n
!! Created by Yuanfu Xie, 2019/04, @GBA-MWF, Shenzhen \n
!! Modified by Yuanfu Xie 2020-4 for adding boundary conditions \n
!! Reforged by Zilong Qin (zilong.qin@gmail.com), 2020/12/17, @GBA-MWF, Shenzhen, add parallelization \n
!! Modified by Zilong Qin, add adjoint option 2024/10/27 \n
!! @copyright (C) 2024 Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!! Shenzhen Institute of Meteorological Innovation (SIMI), All rights reserved.
!! @author Yuanfu Xie, Zilong Qin
! @note
! @warning
! @attention
MODULE poissonSolver_m
  USE geometry_m, ONLY: geometry_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE psMatrix_m, ONLY: psMatrix_t
  USE psSolutn_m, ONLY: psSolutn_t
  USE namelist_ps_m
  USE parameters_m, ONLY: EarthRadius

!> @brief
!! The module of 2D poisson solver for data structure of geometry. The boundary conditions are dirichlet type.
! @see
! @note
! @warning
! @attention
  TYPE :: poissonSolver_t
    TYPE(geometry_t), POINTER :: geometry
    TYPE(psMatrix_t), ALLOCATABLE :: psMat(:)
    CHARACTER(LEN=3) :: solver
    INTEGER(i_kind) :: vLevel
    REAL(r_kind), ALLOCATABLE :: omegas(:)
    INTEGER(i_kind) :: &
      !  ioOpts, &
      !  imethd, &
      nRelax, &
      nCycle, &
      nIterPre, &
      nIterPost, &
      max_band

  CONTAINS

    FINAL :: destructor
    PROCEDURE, PUBLIC :: PoissonSol
    PROCEDURE, PUBLIC :: PoissonSol_sphere
    PROCEDURE, PRIVATE :: fullMcycle
    PROCEDURE, PRIVATE :: singleVcyc
    PROCEDURE, PRIVATE :: iteration

  END TYPE poissonSolver_t
  INTERFACE poissonSolver_t
    PROCEDURE :: constructor
  END INTERFACE poissonSolver_t

CONTAINS

  FUNCTION constructor(configFile, geometry) RESULT(this)
    IMPLICIT NONE
    TYPE(poissonSolver_t) :: this
    TYPE(geometry_t), TARGET, INTENT(IN) :: geometry !< Geometry data structure
    CHARACTER(LEN=1024), INTENT(IN) :: configFile    !< Configuration file
    INTEGER(i_kind) :: i, j, k

    ! transfer the pointer of global geometry to this class
    PRINT *, 'Constructing a Poisson solver... inside constructor'
    this%geometry => geometry
    ALLOCATE (this%psMat(this%geometry%mg%mg_coarsest:this%geometry%mg%mg_finest))

    ! ALLOCATE (this%omegas(10))
    ! Read the configuration from the namelist files
    PRINT *, 'read in poisson namelist... ', TRIM(configFile)
    CALL namelist_ps(configFile, this%solver, this%nCycle, this%nIterPre, &
                     this%nIterPost, this%nRelax, this%omegas, this%max_band)

    ! this%omegas(1) = 1.3895D0
    ! this%omegas(2) = 0.5617D0

    ! this%omegas(1) = 1.951627D0
    ! this%omegas(2) = 0.569688D0

    ! Initialize the matrix
    IF (this%geometry%mpdd%isBaseProc()) WRITE(*,1) this%max_band, this%geometry%mg%mg_coarsest, this%geometry%mg%mg_finest
1   FORMAT('Inside Poisson constructor - max_band: ',I2,' mg_coarsest: ',I2,' mg_finest: ',I2) 
    DO i = this%geometry%mg%mg_coarsest, this%geometry%mg%mg_finest
      this%psMat(i) = psMatrix_t(this%max_band, this%geometry%mg%sg(i))
    END DO
    IF (this%geometry%mpdd%isBaseProc()) PRINT *, 'Poisson initialization finished...'

  END FUNCTION constructor

  SUBROUTINE PoissonSol(this, gLevel, vLevel, rights, solutn, operatorType)
    CLASS(poissonSolver_t) :: this
    INTEGER(i_kind), INTENT(IN) :: gLevel !< Grid level
    INTEGER(i_kind), INTENT(IN) :: vLevel !< Number of vertical levels
    ! Yuanfu Xie has changed the INTENT of rights to IN only replacing original INOUT 2024/12/03
    REAL(r_kind), INTENT(IN) :: rights(vlevel, this%geometry%mg%sg(gLevel)%num_cell) !< Right hand side, with dirichlet boundary conditions at which cell_type /= 0
    REAL(r_kind), INTENT(OUT)   :: solutn(vlevel, this%geometry%mg%sg(gLevel)%num_cell) !< Solution
    CHARACTER(LEN=*), INTENT(IN) :: operatorType !< "forward" or "adjoint"

    ! Local variables:
    INTEGER(i_kind) :: icycle, i, oprType
    TYPE(psSolutn_t), ALLOCATABLE :: psSolutn(:)

    ! Set the operator type
    IF (operatorType .EQ. 'forward') THEN
      oprType = 1
    ELSE IF (operatorType .EQ. 'adjoint') THEN
      oprType = 2
    ELSE
      STOP
      WRITE (*, *) 'No operator type is defined'
    END IF

    ! Initialize the temp solution container
    ! Yuanfu Xie changed this from this%geometry%mg%mg_coarses:gLevel as
    ! it causes problems when this%geometry%mg%mg_coarses != 1
    ALLOCATE (psSolutn(gLevel))
    DO i = this%geometry%mg%mg_coarsest, gLevel
      psSolutn(i) = psSolutn_t(vlevel, this%geometry%mg%sg(i)%num_cell)
    END DO

    psSolutn(gLevel)%rights = rights
    psSolutn(gLevel)%solutn = solutn

    ! Set the Dirichlet boundary conditions to the solution array
    DO i = 1, this%geometry%mg%sg(gLevel)%num_icell
      IF (this%geometry%mg%sg(gLevel)%cell_type(i) .EQ. 1) THEN         ! On the boundary cells, rhs is the stream function
        psSolutn(gLevel)%solutn(:, i) = psSolutn(gLevel)%rights(:, i)
      ELSE IF (this%geometry%mg%sg(gLevel)%cell_type(i) .EQ. 2) THEN         ! On the boundary cells, rhs is the stream function
        psSolutn(gLevel)%solutn(:, i) = psSolutn(gLevel)%rights(:, i)
      END IF
    END DO

    ! Start the solver
    DO icycle = 1, this%nCycle
      IF (this%solver .EQ. 'FMV') THEN
        PRINT *, 'Start gzGrid Poisson FMV cycle...', icycle
        CALL this%fullMcycle(gLevel, vlevel, psSolutn, oprType)

      ELSE IF (this%solver .EQ. 'VVV') THEN
        PRINT *, 'Start gzGrid Poisson VVV cycle...', icycle
        CALL this%singleVcyc(gLevel, vlevel, psSolutn, oprType)
      ELSE
        WRITE (*, *) 'No solver is defned'
      END IF
    END DO

    solutn = psSolutn(gLevel)%solutn

    ! Release the temp.
    DEALLOCATE (psSolutn)
  END SUBROUTINE PoissonSol

  SUBROUTINE PoissonSol_sphere(this, gLevel, vLevel, rights, solutn)
    CLASS(poissonSolver_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vLevel, gLevel
    ! Yuanfu Xie has changed the INTENT of rights to IN only replacing original INOUT 2024/12/03
    REAL(r_kind), INTENT(IN) :: rights(vlevel, this%geometry%mg%sg(gLevel)%num_cell)
    REAL(r_kind), INTENT(OUT) :: solutn(vlevel, this%geometry%mg%sg(gLevel)%num_cell)

    ! Local variables:
    INTEGER(i_kind) :: k
    REAL(r_kind) :: ratio

    ASSOCIATE (sg => this%geometry%mg%sg(gLevel))
      CALL this%PoissonSol(gLevel, vLevel, rights, solutn, "forward")
      DO k = 1, sg%vLevel
        ratio = (EarthRadius + sg%sigma(k)) / EarthRadius
        solutn(k, :) = solutn(k, :) * ratio**4.0D0
      END DO
    END ASSOCIATE

  END SUBROUTINE PoissonSol_sphere

  RECURSIVE SUBROUTINE fullMcycle(this, gLevel, vlevel, psSolutn, oprType)
    CLASS(poissonSolver_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vlevel, gLevel, &
                                   oprType !< 1: forward, 2: adjoint
    TYPE(psSolutn_t), INTENT(INOUT) :: psSolutn(:)

    REAL(r_kind), ALLOCATABLE :: tempState(:, :), worksp(:, :)

    ASSOCIATE (sg => this%geometry%mg%sg(gLevel))

      ! IF (sg%isBaseProc()) THEN
      !   PRINT *, 'Start fullMcycle in gLevel: ', gLevel
      ! END IF

      ALLOCATE (tempState(vLevel, sg%num_cell), &
                worksp(vLevel, sg%num_cell))

      IF (gLevel .LE. this%geometry%mg%mg_coarsest) THEN
        CALL this%singleVcyc(gLevel, vlevel, psSolutn, oprType)
      ELSE
        ! Calculate the residual:
        worksp = 0.0D0
        CALL this%psMat(gLevel)%CalResidus(vlevel, &
                                           psSolutn(gLevel)%rights, &
                                           psSolutn(gLevel)%solutn, &
                                           worksp, &
                                           sg, oprType)
        ! To coarser grid:
        psSolutn(gLevel - 1)%rights = 0.0D0
        CALL this%geometry%mg%restrictionAtGLevel(gLevel, vLevel, worksp, psSolutn(gLevel - 1)%rights, .FALSE.)

        ! Call coarser M cycle: recursive
        psSolutn(gLevel - 1)%solutn = 0.0D0
        CALL this%fullMcycle(gLevel - 1, vLevel, psSolutn, oprType)

        ! Interpolate to coarser solution to current multigrid
        worksp = 0.0D0

        CALL this%geometry%mg%prolongationAtGLevel(gLevel - 1, vLevel, worksp, psSolutn(gLevel - 1)%solutn, .FALSE.)

        ! and correct the solution:
        psSolutn(gLevel)%solutn = psSolutn(gLevel)%solutn + worksp

        CALL this%singleVcyc(gLevel, vLevel, psSolutn, oprType)
      END IF

      DEALLOCATE (tempState, worksp)
    END ASSOCIATE

  END SUBROUTINE fullMcycle

  SUBROUTINE iteration(this, gLevel, vLevel, psSolutn, nIter, oprType)
    CLASS(poissonSolver_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vlevel, gLevel, nIter, &
                                   oprType !< 1: forward, 2: adjoint

    TYPE(psSolutn_t), INTENT(INOUT) :: psSolutn(:)

    REAL(r_kind), ALLOCATABLE :: tempState(:, :)
    INTEGER(i_kind) :: i

    ASSOCIATE (sg => this%geometry%mg%sg(gLevel))
      ALLOCATE (tempState(vLevel, sg%num_cell))

      DO iover = 1, nIter
        tempState = psSolutn(gLevel)%solutn

        ! DO i = 1, sg%num_icell
        !   IF ((sg%cell_type(i) .EQ. 1) .OR. (sg%cell_type(i) .EQ. 2)) THEN         ! On the boundary cells, rhs is the stream function
        !     tempState(:, i) = psSolutn(gLevel)%rights(:, i)
        !     ! psSolutn(gLevel)%solutn(:, i) = psSolutn(gLevel)%rights(:, i)
        !   END IF
        ! END DO
        ! ! CALL sg%ExchangeMatOnHalo2D(vLevel, psSolutn(gLevel)%Solutn)
        ! CALL sg%ExchangeMatOnHalo2D(vLevel, tempState)

        DO itr = 1, this%nRelax
          IF (.TRUE.) THEN
            CALL Jacobi(this%psMat(gLevel)%num_eqns, &
                        sg%num_cell, &
                        vlevel, &
                        this%psMat(gLevel)%max_band, &
                        this%psMat(gLevel)%num_elem, &
                        this%psMat(gLevel)%idx_cols, &
                        this%psMat(gLevel)%row_elem, &
                        this%omegas(itr), &
                        psSolutn(gLevel)%rights, &
                        tempState, &
                        psSolutn(gLevel)%Solutn, oprType)
            IF (oprType == 2) &
              CALL sg%ExchangeMatOnHaloReverseSum2D(vLevel, psSolutn(gLevel)%Solutn)
            CALL Relaxation(psSolutn(gLevel)%Solutn, &
                            tempState, &
                            this%omegas(itr), &
                            sg%num_cell, &
                            this%psMat(gLevel)%num_eqns, &
                            vlevel)
            CALL sg%ExchangeMatOnHalo2D(vLevel, psSolutn(gLevel)%Solutn)
          ELSE
            CALL GaussSeidel(this%psMat(gLevel)%num_eqns, &
                             vlevel, &
                             this%psMat(gLevel)%max_band, &
                             this%psMat(gLevel)%num_elem, &
                             this%psMat(gLevel)%idx_cols, &
                             this%psMat(gLevel)%row_elem, &
                             this%omegas(itr), &
                             psSolutn(gLevel)%rights, &
                             tempState, &
                             psSolutn(gLevel)%Solutn, oprType)

            CALL Relaxation(psSolutn(gLevel)%Solutn, &
                            tempState, &
                            this%omegas(itr), &
                            sg%num_cell, &
                            this%psMat(gLevel)%num_eqns, &
                            vlevel)
          END IF
          !
          IF (itr .LT. this%nRelax) tempState = psSolutn(gLevel)%solutn
        END DO
      END DO

      DEALLOCATE (tempState)
    END ASSOCIATE

  END SUBROUTINE iteration

  RECURSIVE SUBROUTINE singleVcyc(this, gLevel, vlevel, psSolutn, oprType)
    CLASS(poissonSolver_t) :: this
    INTEGER(i_kind), INTENT(IN) :: vlevel, gLevel, &
                                   oprType !< 1: forward, 2: adjoint
    TYPE(psSolutn_t), INTENT(INOUT) :: psSolutn(:)

    REAL(r_kind), ALLOCATABLE :: worksp(:, :)

    ! CALL this%geometry%mpdd%barrier
    ASSOCIATE (sg => this%geometry%mg%sg(gLevel))
      ! ALLOCATE (worksp(vLevel, sg%num_cell))
      CALL sg%allocateMat(vLevel, worksp)

      ! IF (sg%isBaseProc()) THEN
      !   PRINT *, 'Start singleVcyc in gLevel: ', gLevel
      ! END IF

      ! Relaxation 1:
      CALL this%iteration(gLevel, vlevel, psSolutn, this%nIterPre, oprType)

      ! Multigrid start:
      IF (gLevel .GT. this%geometry%mg%mg_coarsest) THEN

        ! Calculate the residual:
        worksp = 0.0D0

        CALL this%psMat(gLevel)%CalResidus(vlevel, &
                                           psSolutn(gLevel)%rights, &
                                           psSolutn(gLevel)%solutn, &
                                           worksp, &
                                           sg, oprType)

        ! To coarser grid:
        psSolutn(gLevel - 1)%rights = 0.0D0
        CALL this%geometry%mg%restrictionAtGLevel(gLevel, vLevel, worksp, psSolutn(gLevel - 1)%rights, .FALSE.)

        ! Call coarser V cycle: recursive
        psSolutn(gLevel - 1)%solutn = 0.0D0
        CALL this%singleVcyc(gLevel - 1, vLevel, psSolutn, oprType)

        ! Interpolate to coarser solution to current multigrid
        worksp = 0.0D0
        CALL this%geometry%mg%prolongationAtGLevel(gLevel - 1, vLevel, worksp, psSolutn(gLevel - 1)%solutn, .FALSE.)

        ! and correct the solution:
        psSolutn(gLevel)%solutn = psSolutn(gLevel)%solutn + worksp
      END IF

      ! Relaxation 2:
      IF (gLevel .GT. this%geometry%mg%mg_coarsest) THEN
        CALL this%iteration(gLevel, vlevel, psSolutn, this%nIterPost, oprType)
      ELSE
        CALL this%iteration(gLevel, vlevel, psSolutn, this%nIterPost, oprType)
      END IF

      ! worksp = 0.0D0
      ! CALL this%psMat(gLevel)%CalResidus(vlevel, &
      !                                    psSolutn(gLevel)%rights, &
      !                                    psSolutn(gLevel)%solutn, &
      !                                    worksp, &
      !                                    sg, oprType)

      ! IF (sg%isBaseProc()) THEN
      !   Print *, 'Max residul is: ', maxval(abs(worksp(1, :))), ' at gLevel: ', gLevel
      !   Print *, 'End cycle in gLevel: ', gLevel
      ! END IF
      ! IF (sg%isBaseProc()) THEN
      !   PRINT *, 'Finish singleVcyc in gLevel: ', gLevel
      ! END IF
      DEALLOCATE (worksp)
    END ASSOCIATE
  END SUBROUTINE singleVcyc

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(poissonSolver_t), INTENT(INOUT) :: this
    IF (ALLOCATED(this%psMat)) DEALLOCATE (this%psMat)
    IF (ALLOCATED(this%omegas)) DEALLOCATE (this%omegas)
  END SUBROUTINE destructor
END MODULE
