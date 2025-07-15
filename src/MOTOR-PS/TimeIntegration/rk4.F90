!!---------------------------------------------------------------------------------------
! PROJECT           : MOTOR-PS
! AFFILIATION       : Self-employed
! AUTOHR(S)         : Jilong CHEN, Yanfu Xie
! VERSION           : V1
! HISTORY           :
!   Created  by Jilong CHEN
!!---------------------------------------------------------------------------------------

MODULE rk4_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE nwp_m, ONLY: nwp_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE YAMLRead_m
  USE CalPres_m, ONLY: CalPres_t
  USE Parameters_m, ONLY: EarthRadius
  USE State2NC_m

  TYPE :: rk4_t
    TYPE(State_t), ALLOCATABLE :: stepsHalf(:), stepsFull(:)
    TYPE(geometry_t), POINTER :: geometryHalf, geometryFull
    ! a Poisson solver:
    TYPE(poissonSolver_t) :: poisson
    TYPE(nwp_t) :: nwp
    INTEGER(i_kind) :: gLevel, numSteps, kt
    REAL(r_kind) :: dt
    INTEGER(i_kind), ALLOCATABLE :: bdy_idx(:)
    CHARACTER(LEN=1024) :: configFile, outputDir, filenameFull, filenameHalf
  CONTAINS
    FINAL :: destructor
    PROCEDURE :: timeIntegrals => timeIntegrals_rk4
  END TYPE rk4_t

  INTERFACE rk4_t
    PROCEDURE :: constructor
  END INTERFACE rk4_t

CONTAINS

  FUNCTION constructor(configFile, gLevel, geometryHalf, geometryFull) RESULT(this)
    IMPLICIT NONE

    TYPE(rk4_t) :: this
    TYPE(geometry_t), TARGET :: geometryHalf, geometryFull
    INTEGER(i_kind), INTENT(IN) :: gLevel
    INTEGER(i_kind) :: istatus, i, ifile
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    REAL(r_kind) :: size_tt_in
    REAL(r_kind), ALLOCATABLE :: temp(:)

    ALLOCATE (this%stepsFull(4), this%stepsHalf(4))

    ifile = yaml_get_var(TRIM(configFile), 'IO', 'output_dir', this%outputDir)

    this%geometryFull => geometryFull
    this%geometryHalf => geometryHalf
    this%configFile = TRIM(configFile)
    this%gLevel = gLevel
    this%poisson = poissonSolver_t(configFile, geometryHalf)
    istatus = 0
    istatus = yaml_get_var(TRIM(this%configFile), 'TimeIntegral', 'time_step', this%dt)
    IF (istatus .NE. 0) THEN
      PRINT *, 'Error, time_step is not set'
      STOP
    END IF
    ASSOCIATE (tt_in => this%geometryFull%mg%sg(gLevel)%tt)
      size_tt_in = SIZE(tt_in)
      this%numSteps = (tt_in(SIZE(tt_in)) - tt_in(1)) / this%dt + 1
      this%kt = (tt_in(2) - tt_in(1)) / this%dt
      PRINT *, 'this%dt is: ', this%dt, this%numSteps, this%kt
    END ASSOCIATE

    this%nwp = nwp_t(this%geometryFull%mg%sg(gLevel), this%geometryHalf%mg%sg(gLevel))

    ASSOCIATE (sg => geometryHalf%mg%sg(gLevel))
      ! temp = PACK(sg%cell_type, MASK=sg%cell_type > 0)
      temp = PACK([(i, i=1, SIZE(sg%cell_type))], sg%cell_type > 0)
      IF (SIZE(temp) .GE. 1) THEN
        ALLOCATE (this%bdy_idx(SIZE(temp)))
        this%bdy_idx = temp
        DEALLOCATE (temp)
      END IF
    END ASSOCIATE

  END FUNCTION

  SUBROUTINE timeIntegrals_rk4(this, XHalf, XFull)
    IMPLICIT NONE

    CLASS(rk4_t) :: this
    TYPE(State_t), INTENT(INOUT) :: XHalf, XFull
    TYPE(State_t) :: XtHalf, XtFull

    ! Local variables:
    INTEGER(i_kind) :: it, i_x, k, i, iv, ik
    REAL(r_kind) :: dk(3), t0
    REAL :: start_time, end_time, start_time_rk4, end_time_rk4, &
            start_time_write, end_time_write, &
            start_time_nwp, end_time_nwp
    CHARACTER(20) :: out_i1, out_i2
    CHARACTER(7) :: varNamesHalf(7) = ['vor    ', 'div    ', 'rho    ', &
                                       'qvapor ', 'pres   ', 'psi    ', 'chi    '], &
                    varNamesFull(2) = ['theta  ', 'w      ']
    CHARACTER(10) :: TenVarNamesHalf(5) = ['tenvor    ', 'tendiv    ', 'tenrho    ', 'tenqvapor ', 'tenpres   '], &
                     TenVarNamesFull(2) = ['tentheta  ', 'tenw      ']
    INTEGER(i_kind) :: div_idx, vor_idx, psi_idx, chi_idx, chi_idxt, pres_idx

    REAL(r_kind), ALLOCATABLE :: vortemp(:, :), divtemp(:, :), ratio(:)

    ! RK4: y := y + (k1 + 2 k2 + 2 k3 + k4) * (dt/6)
    dk = (/0.5D0, 0.5D0, 1.0D0/) * this%dt

    ALLOCATE (ratio(XHalf%sg%vLevel))
    DO i = 1, XHalf%sg%vLevel
      ratio(i) = ((EarthRadius + XHalf%sg%sigma(i)) / EarthRadius)**2.0D0
    END DO

    XtHalf = XHalf
    XtFull = XFull

    DO i = 1, 4
      this%stepsFull(i) = XFull
      this%stepsHalf(i) = XHalf
    END DO

    vor_idx = XtHalf%getVarIdx('vor')
    div_idx = XtHalf%getVarIdx('div')
    pres_idx = XtHalf%getVarIdx('pres')
    psi_idx = XtHalf%getVarIdx('psi')
    chi_idx = XtHalf%getVarIdx('chi')

    ALLOCATE (vortemp(XHalf%sg%vLevel, XHalf%sg%num_cell), &
              divtemp(XHalf%sg%vLevel, XHalf%sg%num_cell))

    ASSOCIATE (num_icell => XFull%sg%num_icell)

      i_x = 2

      DO it = 2, this%numSteps
        PRINT *, 'it is: ', it

        CALL CPU_TIME(start_time_rk4)

        this%stepsHalf(1) = XtHalf
        this%stepsFull(1) = XtFull

        CALL this%nwp%righthands(this%stepsFull(1), this%stepsHalf(1))
        CALL this%nwp%UpdateBdy(this%stepsFull(1), XFull, i_x - 1, 'FULL')
        CALL this%nwp%UpdateBdy(this%stepsHalf(1), XHalf, i_x - 1, 'HALF')

        WRITE (out_i1, '(I0)') it - 1
        WRITE (out_i2, '(I0)') 1

        this%filenameFull = TRIM('Test_warm_bubble_Full_'//TRIM(out_i1)//'_'//TRIM(out_i1))
        this%filenameHalf = TRIM('Test_warm_bubble_Half_'//TRIM(out_i1)//'_'//TRIM(out_i2))

        ! CALL Output_NC_State_AV(this%stepsFull(1), this%outputDir, this%filenameFull)
        ! CALL Output_NC_State_AV(this%stepsHalf(1), this%outputDir, this%filenameHalf)

        ! RK4 for K2-K4
        DO k = 2, 4

          DO iv = 1, 4
            this%stepsHalf(k)%Fields(this%stepsHalf(k)%getVarIdx(varNamesHalf(iv)))%DATA(:, 1:num_icell, 1) = &
              this%stepsHalf(1)%Fields(this%stepsHalf(1)%getVarIdx(varNamesHalf(iv)))%DATA(:, 1:num_icell, 1) + &
              dk(k - 1) * &
              this%stepsHalf(k - 1)%Fields(this%stepsHalf(k - 1)%getVarIdx(TenVarNamesHalf(iv)))%DATA(:, 1:num_icell, 1)

            CALL this%stepsHalf(k)%sg%ExchangeMatOnHalo2D(this%stepsHalf(k)%sg%vLevel, &
                                                          this%stepsHalf(k)%Fields(this%stepsHalf(k)%getVarIdx(varNamesHalf(iv)))%DATA(:, :, 1))
          END DO
          DO iv = 1, 2
            this%stepsFull(k)%Fields(this%stepsFull(k)%getVarIdx(varNamesFull(iv)))%DATA(:, 1:num_icell, 1) = &
              this%stepsFull(1)%Fields(this%stepsFull(1)%getVarIdx(varNamesFull(iv)))%DATA(:, 1:num_icell, 1) + &
              dk(k - 1) * &
              this%stepsFull(k - 1)%Fields(this%stepsFull(k - 1)%getVarIdx(TenVarNamesFull(iv)))%DATA(:, 1:num_icell, 1)

            CALL this%stepsFull(k)%sg%ExchangeMatOnHalo2D(this%stepsFull(k)%sg%vLevel, &
                                                          this%stepsFull(k)%Fields(this%stepsFull(k)%getVarIdx(varNamesFull(iv)))%DATA(:, :, 1))
          END DO

          CALL this%nwp%CalPres%CalPres(this%stepsHalf(k)%Fields(this%stepsHalf(k)%getVarIdx('rho'))%DATA(:, :, 1), &
                                        this%stepsFull(k)%Fields(this%stepsFull(k)%getVarIdx('theta'))%DATA(:, :, 1), &
                                        this%stepsHalf(k)%Fields(this%stepsHalf(k)%getVarIdx('pres'))%DATA(:, :, 1))

          CALL this%stepsHalf(k)%sg%ExchangeMatOnHalo2D(this%stepsHalf(k)%sg%vLevel, &
                                                        this%stepsHalf(k)%Fields(this%stepsHalf(k)%getVarIdx('pres'))%DATA(:, :, 1))

          WRITE (out_i1, '(I0)') it - 1
          WRITE (out_i2, '(I0)') k

          this%filenameFull = TRIM('Test_warm_bubble_Full_'//TRIM(out_i1)//'_'//TRIM(out_i2)//'_before_Poisson')
          this%filenameHalf = TRIM('Test_warm_bubble_Half_'//TRIM(out_i1)//'_'//TRIM(out_i2)//'_before_Poisson')

          ! CALL Output_NC_State_AV(this%stepsFull(k), this%outputDir, this%filenameFull)
          ! CALL Output_NC_State_AV(this%stepsHalf(k), this%outputDir, this%filenameHalf)

          vortemp = this%stepsHalf(k)%Fields(vor_idx)%DATA(:, :, 1)
          divtemp = this%stepsHalf(k)%Fields(div_idx)%DATA(:, :, 1)

          DO ik = 1, XHalf%sg%vLevel
            vortemp(ik, :) = vortemp(ik, :) * ratio(ik)
            divtemp(ik, :) = divtemp(ik, :) * ratio(ik)
          END DO

          IF (ALLOCATED(this%bdy_idx)) THEN
            vortemp(:, this%bdy_idx) = this%stepsHalf(k)%Fields(psi_idx)%DATA(:, this%bdy_idx, 1)
            divtemp(:, this%bdy_idx) = this%stepsHalf(k)%Fields(chi_idx)%DATA(:, this%bdy_idx, 1)
          END IF

          CALL this%stepsHalf(k)%sg%ExchangeMatOnHalo2D(this%stepsHalf(k)%sg%vLevel, vortemp)
          CALL this%stepsHalf(k)%sg%ExchangeMatOnHalo2D(this%stepsHalf(k)%sg%vLevel, divtemp)

          this%stepsHalf(k)%Fields(chi_idx)%DATA(:, :, 1) = this%stepsHalf(k - 1)%Fields(chi_idx)%DATA(:, :, 1)
          this%stepsHalf(k)%Fields(psi_idx)%DATA(:, :, 1) = this%stepsHalf(k - 1)%Fields(psi_idx)%DATA(:, :, 1)

          CALL CPU_TIME(start_time)
          CALL this%poisson%PoissonSol(this%gLevel, this%stepsHalf(k)%sg%vLevel, &
                                       divtemp, &
                                       this%stepsHalf(k)%Fields(chi_idx)%DATA(:, :, 1))
          CALL this%poisson%PoissonSol(this%gLevel, this%stepsHalf(k)%sg%vLevel, &
                                       vortemp, &
                                       this%stepsHalf(k)%Fields(psi_idx)%DATA(:, :, 1))
          CALL CPU_TIME(end_time)
          PRINT *, 'poisson solver takes time', end_time - start_time, 'S'
          CALL this%stepsHalf(k)%sg%ExchangeMatOnHalo2D(this%stepsHalf(k)%sg%vLevel, &
                                                        this%stepsHalf(k)%Fields(chi_idx)%DATA(:, :, 1))
          CALL this%stepsHalf(k)%sg%ExchangeMatOnHalo2D(this%stepsHalf(k)%sg%vLevel, &
                                                        this%stepsHalf(k)%Fields(psi_idx)%DATA(:, :, 1))
          CALL CPU_TIME(start_time_nwp)
          CALL this%nwp%righthands(this%stepsFull(k), this%stepsHalf(k))
          CALL CPU_TIME(end_time_nwp)
          PRINT *, 'nwp takes time', end_time_nwp - start_time_nwp, 'S'
          CALL this%nwp%UpdateBdy(this%stepsFull(k), XFull, i_x - 1, 'FULL')
          CALL this%nwp%UpdateBdy(this%stepsHalf(k), XHalf, i_x - 1, 'HALF')

          this%filenameFull = TRIM('Test_warm_bubble_Full_'//TRIM(out_i1)//'_'//TRIM(out_i2)//'_after_Poisson')
          this%filenameHalf = TRIM('Test_warm_bubble_Half_'//TRIM(out_i1)//'_'//TRIM(out_i2)//'_after_Poisson')

          ! CALL Output_NC_State_AV(this%stepsFull(k), this%outputDir, this%filenameFull)
          ! CALL Output_NC_State_AV(this%stepsHalf(k), this%outputDir, this%filenameHalf)

        END DO

        DO iv = 1, 4
          XtHalf%Fields(XtHalf%getVarIdx(varNamesHalf(iv)))%DATA(:, 1:num_icell, 1) = &
            this%stepsHalf(1)%Fields(this%stepsHalf(1)%getVarIdx(varNamesHalf(iv)))%DATA(:, 1:num_icell, 1) + &
            (this%stepsHalf(1)%Fields(this%stepsHalf(1)%getVarIdx(TenVarNamesHalf(iv)))%DATA(:, 1:num_icell, 1) + &
             this%stepsHalf(2)%Fields(this%stepsHalf(2)%getVarIdx(TenVarNamesHalf(iv)))%DATA(:, 1:num_icell, 1) * 2.0D0 + &
             this%stepsHalf(3)%Fields(this%stepsHalf(3)%getVarIdx(TenVarNamesHalf(iv)))%DATA(:, 1:num_icell, 1) * 2.0D0 + &
             this%stepsHalf(4)%Fields(this%stepsHalf(4)%getVarIdx(TenVarNamesHalf(iv)))%DATA(:, 1:num_icell, 1)) &
            * this%dt / 6.0D0

          CALL XtHalf%sg%ExchangeMatOnHalo2D(XtHalf%sg%vLevel, &
                                             XtHalf%Fields(XtHalf%getVarIdx(varNamesHalf(iv)))%DATA(:, :, 1))
        END DO
        DO iv = 1, 2
          XtFull%Fields(XtFull%getVarIdx(varNamesFull(iv)))%DATA(:, 1:num_icell, 1) = &
            this%stepsFull(1)%Fields(this%stepsFull(1)%getVarIdx(varNamesFull(iv)))%DATA(:, 1:num_icell, 1) + &
            (this%stepsFull(1)%Fields(this%stepsFull(1)%getVarIdx(TenVarNamesFull(iv)))%DATA(:, 1:num_icell, 1) + &
             this%stepsFull(2)%Fields(this%stepsFull(2)%getVarIdx(TenVarNamesFull(iv)))%DATA(:, 1:num_icell, 1) * 2.0D0 + &
             this%stepsFull(3)%Fields(this%stepsFull(3)%getVarIdx(TenVarNamesFull(iv)))%DATA(:, 1:num_icell, 1) * 2.0D0 + &
             this%stepsFull(4)%Fields(this%stepsFull(4)%getVarIdx(TenVarNamesFull(iv)))%DATA(:, 1:num_icell, 1)) &
            * this%dt / 6.0D0

          CALL XtFull%sg%ExchangeMatOnHalo2D(XtFull%sg%vLevel, &
                                             XtFull%Fields(XtFull%getVarIdx(varNamesFull(iv)))%DATA(:, :, 1))
        END DO

        CALL this%nwp%CalPres%CalPres(XtHalf%Fields(XtHalf%getVarIdx('rho'))%DATA(:, :, 1), &
                                      XtFull%Fields(XtFull%getVarIdx('theta'))%DATA(:, :, 1), &
                                      XtHalf%Fields(XtHalf%getVarIdx('pres'))%DATA(:, :, 1))

        CALL XtHalf%sg%ExchangeMatOnHalo2D(XtHalf%sg%vLevel, XtHalf%Fields(XtHalf%getVarIdx('pres'))%DATA(:, :, 1))

        WRITE (out_i1, '(I0)') it
        ! this%filenameFull = TRIM('Test_warm_bubble_Full_Xt_'//TRIM(out_i1)//'_before_Poisson')
        ! this%filenameHalf = TRIM('Test_warm_bubble_Half_Xt_'//TRIM(out_i1)//'_before_Poisson')

        ! CALL Output_NC_State_AV(XtFull, this%outputDir, this%filenameFull)
        ! CALL Output_NC_State_AV(XtHalf, this%outputDir, this%filenameHalf)

        vortemp = XtHalf%Fields(vor_idx)%DATA(:, :, 1)
        divtemp = XtHalf%Fields(div_idx)%DATA(:, :, 1)

        DO ik = 1, XHalf%sg%vLevel
          vortemp(ik, :) = vortemp(ik, :) * ratio(ik)
          divtemp(ik, :) = divtemp(ik, :) * ratio(ik)
        END DO

        IF (ALLOCATED(this%bdy_idx)) THEN
          vortemp(:, this%bdy_idx) = this%stepsHalf(1)%Fields(psi_idx)%DATA(:, this%bdy_idx, 1)
          divtemp(:, this%bdy_idx) = this%stepsHalf(1)%Fields(chi_idx)%DATA(:, this%bdy_idx, 1)
        END IF
        CALL CPU_TIME(start_time)
        CALL this%poisson%PoissonSol(this%gLevel, XtHalf%sg%vLevel, &
                                     divtemp, &
                                     XtHalf%Fields(chi_idx)%DATA(:, :, 1))
        CALL this%poisson%PoissonSol(this%gLevel, XtHalf%sg%vLevel, &
                                     vortemp, &
                                     XtHalf%Fields(psi_idx)%DATA(:, :, 1))
        CALL CPU_TIME(end_time)
        PRINT *, 'poisson solver takes time', end_time - start_time, 'S'
        CALL XtHalf%sg%ExchangeMatOnHalo2D(XtHalf%sg%vLevel, XtHalf%Fields(psi_idx)%DATA(:, :, 1))
        CALL XtHalf%sg%ExchangeMatOnHalo2D(XtHalf%sg%vLevel, XtHalf%Fields(chi_idx)%DATA(:, :, 1))

        PRINT *, "max of w in ", it, "step is: ", MAXVAL(dabs(XtFull%Fields(XtFull%getVarIdx(varNamesFull(2)))%DATA(:, 1:num_icell, 1)))

        this%filenameFull = TRIM('Test_warm_bubble_Full_Xt_'//TRIM(out_i1)//'_after_Poisson')
        this%filenameHalf = TRIM('Test_warm_bubble_Half_Xt_'//TRIM(out_i1)//'_after_Poisson')
        IF (MOD(it, 50) .EQ. 0) THEN
          CALL CPU_TIME(start_time_write)
          CALL Output_NC_State_AV(XtFull, this%outputDir, this%filenameFull)
          CALL Output_NC_State_AV(XtHalf, this%outputDir, this%filenameHalf)
          CALL CPU_TIME(end_time_write)
          PRINT *, 'write file takes time', end_time_write - start_time_write, 'S'
        END IF

        CALL CPU_TIME(end_time_rk4)
        PRINT *, 'rk4 takes time', end_time_rk4 - start_time_rk4, 'S'

        ! Fill the integrated vars to X Fields
        IF (MOD(it, this%kt) .EQ. 1) THEN
          DO iv = 1, SIZE(varNamesHalf)
            XHalf%Fields(XHalf%getVarIdx(varNamesHalf(iv)))%DATA(:, :, i_x) = XtHalf%Fields(XHalf%getVarIdx(varNamesHalf(iv)))%DATA(:, :, 1)
          END DO
          DO iv = 1, SIZE(varNamesFull)
            XFull%Fields(XFull%getVarIdx(varNamesFull(iv)))%DATA(:, :, i_x) = XtFull%Fields(XFull%getVarIdx(varNamesFull(iv)))%DATA(:, :, 1)
          END DO

          i_x = i_x + 1
        END IF

      END DO

    END ASSOCIATE

  END SUBROUTINE timeIntegrals_rk4

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(rk4_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%stepsFull)) DEALLOCATE (this%stepsFull)
    IF (ALLOCATED(this%stepsHalf)) DEALLOCATE (this%stepsHalf)
    IF (ALLOCATED(this%bdy_idx)) DEALLOCATE (this%bdy_idx)
    IF (ASSOCIATED(this%geometryHalf)) NULLIFY (this%geometryHalf)
    IF (ASSOCIATED(this%geometryFull)) NULLIFY (this%geometryFull)

  END SUBROUTINE destructor

END MODULE rk4_m
