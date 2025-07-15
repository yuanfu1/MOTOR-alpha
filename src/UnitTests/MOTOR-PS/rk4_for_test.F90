!!---------------------------------------------------------------------------------------
! PROJECT           : MOTOR-PS
! AFFILIATION       : Self-employed
! AUTOHR(S)         : Jilong CHEN, Yanfu Xie
! VERSION           : V 1
! HISTORY           :
!   Created  by Jilong CHEN
!!---------------------------------------------------------------------------------------

MODULE rk4test_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE parameters_m, ONLY: k_d, r_d
  USE nwp_m, ONLY: nwp_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE poissonSolver_m, ONLY: poissonSolver_t
  USE YAMLRead_m
  USE CalPres_m, ONLY: CalPres_t
  USE UnitTestDyn_m, ONLY: UnitTestDyn_t
  USE State2NC_m
  USE parameters_m, ONLY: EarthRadius

  TYPE :: rk4_t
    TYPE(State_t), ALLOCATABLE :: stepsHalf(:), stepsFull(:)
    TYPE(geometry_t), POINTER :: geometryHalf, geometryFull
    ! a Poisson solver:
    TYPE(poissonSolver_t) :: poisson
    TYPE(nwp_t) :: nwp
    TYPE(UnitTestDyn_t) :: RossbyHalf, RossbyFull
    INTEGER(i_kind) :: gLevel, numSteps, kt
    REAL(r_kind) :: dt
    INTEGER(i_kind), ALLOCATABLE :: bdy_idx(:)
    CHARACTER(LEN=1024) :: configFile, outputDir, filenameFull, filenameHalf

  CONTAINS
    FINAL :: destructor
    PROCEDURE :: timeIntegrals => timeIntegrals_rk4
    PROCEDURE :: tendency_ana_minus_simu
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
    PRINT *, 'Debug 4-0'
    this%poisson = poissonSolver_t(configFile, geometryHalf)
    PRINT *, 'Debug 4-01'
    istatus = 0
    istatus = yaml_get_var(TRIM(this%configFile), 'TimeIntegral', 'time_step', this%dt)
    IF (istatus .NE. 0) THEN
      PRINT *, 'time_step is not set'
      STOP
    END IF

    ASSOCIATE (tt_in => this%geometryFull%mg%sg(gLevel)%tt)
      size_tt_in = SIZE(tt_in)
      this%numSteps = (tt_in(SIZE(tt_in)) - tt_in(1)) / this%dt + 1
      this%kt = (tt_in(2) - tt_in(1)) / this%dt
    END ASSOCIATE

    this%dt = this%dt / 1000.D0

    this%nwp = nwp_t(geometryFull%mg%sg(gLevel), geometryHalf%mg%sg(gLevel))

    ASSOCIATE (sg => geometryHalf%mg%sg(gLevel))
      CALL this%RossbyHalf%initial(geometryHalf%mg%sg(gLevel))
      CALL this%RossbyHalf%Hght_func(sg%cell_cntr, this%RossbyHalf%zhght, this%RossbyHalf%z_s, this%RossbyHalf%Hz)
      PRINT *, "max of zhght and Hz in RossbyHalf are: ", MAXVAL(this%RossbyHalf%zhght(1, :)), MAXVAL(this%RossbyHalf%Hz(1, :))
      PRINT *, "max of zhght and Hz in geometryHalf are: ", MAXVAL(sg%zhght(1, :)), MAXVAL(sg%Hz(1, :))

    END ASSOCIATE

    ASSOCIATE (sg => geometryFull%mg%sg(gLevel))
      CALL this%RossbyFull%initial(geometryFull%mg%sg(gLevel))
      CALL this%RossbyFull%Hght_func(sg%cell_cntr, this%RossbyFull%zhght, this%RossbyFull%z_s, this%RossbyFull%Hz)
      PRINT *, "max of zhght and Hz in RossbyFull are: ", MAXVAL(this%RossbyFull%zhght(1, :)), MAXVAL(this%RossbyFull%Hz(1, :))
      PRINT *, "max of zhght and Hz in geometryFull are: ", MAXVAL(sg%zhght(1, :)), MAXVAL(sg%Hz(1, :))
    END ASSOCIATE

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

  SUBROUTINE tendency_ana_minus_simu(this, XFull, XHalf)
    IMPLICIT NONE

    CLASS(rk4_t) :: this
    TYPE(State_t), INTENT(INOUT) :: XHalf, XFull

    ! Local variables:
    INTEGER(i_kind) :: it, i_x, k, i
    CHARACTER(10) :: TenVarNamesHalf(5) = ['tenvor    ', 'tendiv    ', 'tenrho    ', 'tenqvapor ', 'tenpres   '], &
                     TenVarNamesHalft(5) = ['tenvort   ', 'tendivt   ', 'tenrhot   ', 'tenqvaport', 'tenprest  '], &
                     TenVarNamesFull(2) = ['tentheta  ', 'tenw      '], &
                     TenVarNamesFullt(2) = ['tenthetat ', 'tenwt     ']
    CHARACTER(7) :: mVarNamesHalf(5) = ['mvor   ', 'mdiv   ', 'mrho   ', 'mqvapor', 'mpres  '], &
                    mVarNamesFull(2) = ['mtheta ', 'mw     ']
    IF (XHalf%getVarIdx(mVarNamesHalf(1)) .EQ. 0) THEN
      DO i = 1, SIZE(mVarNamesHalf, 1)
        CALL XHalf%addVar(mVarNamesHalf(i))
      END DO
      DO i = 1, SIZE(mVarNamesFull, 1)
        CALL XFull%addVar(mVarNamesFull(i))
      END DO
    END IF

    DO i = 1, SIZE(TenVarNamesHalf, 1)
      XHalf%Fields(XHalf%getVarIdx(mVarNamesHalf(i)))%DATA(:, :, 1) = &
        XHalf%Fields(XHalf%getVarIdx(TenVarNamesHalft(i)))%DATA(:, :, 1) - &
        XHalf%Fields(XHalf%getVarIdx(TenVarNamesHalf(i)))%DATA(:, :, 1)
    END DO

    DO i = 1, SIZE(TenVarNamesFull, 1)
      XFull%Fields(XFull%getVarIdx(mVarNamesFull(i)))%DATA(:, :, 1) = &
        XFull%Fields(XFull%getVarIdx(TenVarNamesFullt(i)))%DATA(:, :, 1) - &
        XFull%Fields(XFull%getVarIdx(TenVarNamesFull(i)))%DATA(:, :, 1)
    END DO

  END SUBROUTINE

  CHARACTER(LEN=20) FUNCTION INT_TO_STRING(num)
    INTEGER, INTENT(IN) :: num
    WRITE (INT_TO_STRING, '(I0)') num
  END FUNCTION INT_TO_STRING

  SUBROUTINE UpdateBdy_t(Xin, XBdy, layertype)
    TYPE(State_t), INTENT(IN) :: XBdy
    TYPE(State_t), INTENT(INOUT) :: Xin
    ! INTEGER(i_kind), INTENT(IN) :: k
    CHARACTER(4), INTENT(IN) :: layertype
    INTEGER(i_kind) :: i, j, st, sHalf, sFull
    REAL(r_kind) :: dt
    REAL(r_kind), ALLOCATABLE :: BufferRegion
    CHARACTER(10) :: TenVarNamesHalf(4) = ['tenvort   ', 'tendivt   ', 'tenrhot   ', 'tenqvaport'], &
                     TenVarNamesFull(2) = ['tenwt     ', 'tenthetat ']
    CHARACTER(9) :: varNamesHalf(4) = ['tenvor   ', 'tendiv   ', 'tenrho   ', 'tenqvapor'], &
                    varNamesFull(2) = ['tenw     ', 'tentheta ']

    sHalf = SIZE(varNamesHalf)
    sFull = SIZE(varNamesFull)

    ASSOCIATE (sg => Xin%sg)

      IF ((.NOT. ALLOCATED(sg%bdy_type)) .OR. (.NOT. ALLOCATED(sg%BdyAlpha))) THEN
        PRINT *, '#==============================='
        PRINT *, 'ERROR, sg%bdy_type should be set and BufferWeight should be calculated before UpdateBdy'
        PRINT *, 'Use GenBdYIdx in GenBdy module to set sg%bdy_type'
        PRINT *, 'Use GenBdyBufferWeight in GenBdy module to set sg%BdyAlpha and sg%BdyBeta'
        PRINT *, '#==============================='
        STOP
      END IF

      PRINT *, 'sFull is: ', sFull

      IF (.NOT. ALLOCATED(sg%BufferHalo)) RETURN

      PRINT *, 'Debug 100'
      IF (layertype .EQ. 'HALF') THEN
        DO i = 1, sHalf
          PRINT *, 'max of bdy values are: ', TRIM(TenVarNamesHalf(i)), MAXVAL(XBdy%Fields(XBdy%getVarIdx(TRIM(TenVarNamesHalf(i))))%DATA(:, sg%BufferHalo, 1))
          DO j = 1, sg%vLevel
            Xin%Fields(Xin%getVarIdx(TRIM(varNamesHalf(i))))%DATA(j, sg%BufferHalo, 1) = &
              Xin%Fields(Xin%getVarIdx(TRIM(varNamesHalf(i))))%DATA(j, sg%BufferHalo, 1) * sg%BdyBeta(sg%BufferHalo) &
              + XBdy%Fields(XBdy%getVarIdx(TRIM(TenVarNamesHalf(i))))%DATA(j, sg%BufferHalo, 1) * sg%BdyAlpha(sg%BufferHalo)
          END DO
        END DO
      ELSEIF (layertype .EQ. 'FULL') THEN
        DO i = 1, sFull
          PRINT *, 'max of bdy values are: ', TRIM(TenVarNamesHalf(i)), MAXVAL(XBdy%Fields(XBdy%getVarIdx(TRIM(TenVarNamesFull(i))))%DATA(:, sg%BufferHalo, 1))
          DO j = 1, sg%vLevel
            Xin%Fields(Xin%getVarIdx(TRIM(varNamesFull(i))))%DATA(j, sg%BufferHalo, 1) = &
              Xin%Fields(Xin%getVarIdx(TRIM(varNamesFull(i))))%DATA(j, sg%BufferHalo, 1) * sg%BdyBeta(sg%BufferHalo) &
              + XBdy%Fields(XBdy%getVarIdx(TRIM(TenVarNamesFull(i))))%DATA(j, sg%BufferHalo, 1) * sg%BdyAlpha(sg%BufferHalo)
          END DO
        END DO
      ELSE
        PRINT *, 'ERROR: Layer type should be FULL or HALF'
        STOP
      END IF
    END ASSOCIATE

  END SUBROUTINE

  SUBROUTINE timeIntegrals_rk4(this, XHalf, XFull)
    IMPLICIT NONE

    CLASS(rk4_t) :: this
    TYPE(State_t), INTENT(INOUT) :: XHalf, XFull
    TYPE(State_t) :: XtHalf, XtFull
    ! Total time integration steps and output intervals:
    ! REAL(r_kind), INTENT(IN) :: dt

    ! Local variables:
    INTEGER(i_kind) :: it, i_x, k, i, iv, ik
    REAL(r_kind) :: dk(3), t0, error, error1
    REAL :: start_time, end_time
    CHARACTER(7) :: varNamesHalf(7) = ['vor    ', 'div    ', 'rho    ', &
                                       'qvapor ', 'pres   ', 'psi    ', 'chi    '], &
                    varNamesHalft(7) = ['vort   ', 'divt   ', 'rhot   ', &
                                        'qvaport', 'prest  ', 'psit   ', 'chit   '], &
                    varNamesFull(2) = ['theta  ', 'w      '], &
                    varNamesFullt(2) = ['thetat ', 'wt     '], &
                    mVarNamesHalf(5) = ['mvor   ', 'mdiv   ', 'mrho   ', 'mqvapor', 'mpres  '], &
                    mVarNamesFull(2) = ['mtheta ', 'mw     ']
    CHARACTER(10) :: TenVarNamesHalf(5) = ['tenvor    ', 'tendiv    ', 'tenrho    ', 'tenqvapor ', 'tenpres   '], &
                     TenVarNamesHalft(5) = ['tenvort   ', 'tendivt   ', 'tenrhot   ', 'tenqvaport', 'tenprest  '], &
                     TenVarNamesFull(2) = ['tentheta  ', 'tenw      '], &
                     TenVarNamesFullt(2) = ['tenthetat ', 'tenwt     '], &
                     ParVarNamesHalft(5) = ['parvort   ', 'pardivt   ', 'parrhot   ', 'parqvaport', 'parprest  '], &
                     ParVarNamesFullt(2) = ['parthetat ', 'parwt     ']

    INTEGER(i_kind) :: div_idx, vor_idx, psi_idx, chi_idx, psi_idxt, chi_idxt, pres_idx!, q_idx, rho_idx, &
    !                    pres_idx, w_idx, theta_idx, &
    !                    tdiv_idx, tvor_idx, tq_idx, trho_idx, &
    !                    tw_idx, ttheta_idx, tpres_idx, &
    !                    div_idxt, vor_idxt, psi_idxt, chi_idxt, q_idxt, rho_idxt, &
    !                    pres_idxt, w_idxt, theta_idxt, &
    !                    tdiv_idxt, tvor_idxt, tq_idxt, trho_idxt, &
    !                    tw_idxt, ttheta_idxt, &
    !                    pdiv_idxt, pvor_idxt, pq_idxt, prho_idxt, &
    !                    pw_idxt, ptheta_idxt, ppres_idxt
    REAL(r_kind), ALLOCATABLE :: vortemp(:, :), divtemp(:, :), ratio(:)

    ! RK4: y := y + (k1 + 2 k2 + 2 k3 + k4) * (dt/6)
    dk = (/0.5D0, 0.5D0, 1.0D0/)

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
    ! rho_idx = XtHalf%getVarIdx('rho')
    ! q_idx = XtHalf%getVarIdx('qvapor')
    pres_idx = XtHalf%getVarIdx('pres')
    psi_idx = XtHalf%getVarIdx('psi')
    chi_idx = XtHalf%getVarIdx('chi')
    ! w_idx = XtFull%getVarIdx('w')
    ! theta_idx = XtFull%getVarIdx('theta')

    ! tvor_idx = XtHalf%getVarIdx('tenvor')
    ! tdiv_idx = XtHalf%getVarIdx('tendiv')
    ! trho_idx = XtHalf%getVarIdx('tenrho')
    ! tq_idx = XtHalf%getVarIdx('tenqvapor')
    ! tw_idx = XtFull%getVarIdx('tenw')
    ! ttheta_idx = XtFull%getVarIdx('tentheta')
    ! tpres_idx = XtHalf%getVarIdx('tenpres')

    psi_idxt = XtHalf%getVarIdx('psit')
    chi_idxt = XtHalf%getVarIdx('chit')
    ! vor_idxt = XtHalf%getVarIdx('vort')
    ! div_idxt = XtHalf%getVarIdx('divt')
    ! rho_idxt = XtHalf%getVarIdx('rhot')
    ! q_idxt = XtHalf%getVarIdx('qvaport')
    ! pres_idxt = XtHalf%getVarIdx('prest')
    ! w_idxt = XtFull%getVarIdx('wt')
    ! theta_idxt = XtFull%getVarIdx('thetat')

    ! tvor_idxt = XtHalf%getVarIdx('tenvort')
    ! tdiv_idxt = XtHalf%getVarIdx('tendivt')
    ! trho_idxt = XtHalf%getVarIdx('tenrhot')
    ! tq_idxt = XtHalf%getVarIdx('tenqvaport')
    ! tw_idxt = XtFull%getVarIdx('tenwt')
    ! ttheta_idxt = XtFull%getVarIdx('tenthetat')

    ! pvor_idxt = XtHalf%getVarIdx('parvort')
    ! pdiv_idxt = XtHalf%getVarIdx('pardivt')
    ! prho_idxt = XtHalf%getVarIdx('parrhot')
    ! pq_idxt = XtHalf%getVarIdx('parqvaport')
    ! ptheta_idxt = XtFull%getVarIdx('parthetat')
    ! ppres_idxt = XtFull%getVarIdx('parprest')

    ! PRINT *, 'div_idx is: ', chi_idxt, vor_idx

    t0 = 0.0D0
    ALLOCATE (vortemp(XHalf%sg%vLevel, XHalf%sg%num_cell), &
              divtemp(XHalf%sg%vLevel, XHalf%sg%num_cell))

    ASSOCIATE (num_icell => XFull%sg%num_icell)

      i_x = 2

      DO it = 2, this%numSteps
        ! K1 = f(t_n,y_n):
        PRINT *, 'i_x is: ', i_x

        this%stepsHalf(1) = XtHalf
        this%stepsFull(1) = XtFull

        CALL this%nwp%righthands(this%stepsFull(1), this%stepsHalf(1))

        PRINT *, 'dt in time step ', it - 1, ' rk step 1 ', 'is ', t0 + (it - 2) * this%dt

        ! CALL this%RossbyHalf%righthandsHalf(t0 + (it - 2)*this%dt, this%stepsHalf(1))
        ! CALL this%RossbyFull%righthandsFull(t0 + (it - 2)*this%dt, this%stepsFull(1))

        CALL this%RossbyHalf%righthandsHalf(t0, this%stepsHalf(1))
        CALL this%RossbyFull%righthandsFull(t0, this%stepsFull(1))

        this%filenameFull = TRIM('Test_rk4_Full_d1_'//TRIM(INT_TO_STRING(it - 1))//'_'//TRIM(INT_TO_STRING(1))//'_before_udtbdy')
        this%filenameHalf = TRIM('Test_rk4_Half_d1_'//TRIM(INT_TO_STRING(it - 1))//'_'//TRIM(INT_TO_STRING(1))//'_before_udtbdy')

        ! CALL Output_NC_State_AV(this%stepsFull(1), this%outputDir, this%filenameFull)
        ! CALL Output_NC_State_AV(this%stepsHalf(1), this%outputDir, this%filenameHalf)

        ! CALL this%nwp%UpdateBdy(this%stepsFull(1), XFull, i_x - 1, 'FULL')
        ! CALL this%nwp%UpdateBdy(this%stepsHalf(1), XHalf, i_x - 1, 'HALF')

        PRINT *, 'Debug udpatebdy_t ready'

        CALL UpdateBdy_t(this%stepsFull(1), this%stepsFull(1), 'FULL')
        CALL UpdateBdy_t(this%stepsHalf(1), this%stepsHalf(1), 'HALF')

        PRINT *, 'Debug udpatebdy_t finish'

        this%stepsHalf(1)%Fields(this%stepsHalf(1)%getVarIdx(TenVarNamesHalf(5)))%DATA(:, 1:num_icell, 1) = &
          this%stepsHalf(1)%Fields(this%stepsHalf(1)%getVarIdx(ParVarNamesHalft(5)))%DATA(:, 1:num_icell, 1)

        CALL this%tendency_ana_minus_simu(this%stepsFull(1), this%stepsHalf(1))

        ! PRINT *, 'Test in_to_string: ', INT_TO_STRING(it - 1)

        this%filenameFull = TRIM('Test_rk4_Full_d1_'//TRIM(INT_TO_STRING(it - 1))//'_'//TRIM(INT_TO_STRING(1)))
        this%filenameHalf = TRIM('Test_rk4_Half_d1_'//TRIM(INT_TO_STRING(it - 1))//'_'//TRIM(INT_TO_STRING(1)))

        ! CALL Output_NC_State_AV(this%stepsFull(1), this%outputDir, this%filenameFull)
        ! CALL Output_NC_State_AV(this%stepsHalf(1), this%outputDir, this%filenameHalf)

        ! RK4 for K2-K4
        DO k = 2, 4

          DO iv = 1, 4
            this%stepsHalf(k)%Fields(this%stepsHalf(k)%getVarIdx(varNamesHalf(iv)))%DATA(:, 1:num_icell, 1) = &
              this%stepsHalf(1)%Fields(this%stepsHalf(1)%getVarIdx(varNamesHalf(iv)))%DATA(:, 1:num_icell, 1) + &
              dk(k - 1) * &
              this%stepsHalf(k - 1)%Fields(this%stepsHalf(k - 1)%getVarIdx(mvarNamesHalf(iv)))%DATA(:, 1:num_icell, 1)

            CALL this%stepsHalf(k)%sg%ExchangeMatOnHalo2D(this%stepsHalf(k)%sg%vLevel, &
                                                          this%stepsHalf(k)%Fields(this%stepsHalf(k)%getVarIdx(varNamesHalf(iv)))%DATA(:, :, 1))
          END DO
          DO iv = 1, 1
            this%stepsFull(k)%Fields(this%stepsFull(k)%getVarIdx(varNamesFull(iv)))%DATA(:, 1:num_icell, 1) = &
              this%stepsFull(1)%Fields(this%stepsFull(1)%getVarIdx(varNamesFull(iv)))%DATA(:, 1:num_icell, 1) + &
              dk(k - 1) * &
              this%stepsFull(k - 1)%Fields(this%stepsFull(k - 1)%getVarIdx(mVarNamesFull(iv)))%DATA(:, 1:num_icell, 1)

            CALL this%stepsFull(k)%sg%ExchangeMatOnHalo2D(this%stepsFull(k)%sg%vLevel, &
                                                          this%stepsFull(k)%Fields(this%stepsFull(k)%getVarIdx(varNamesFull(iv)))%DATA(:, :, 1))
          END DO
          this%stepsFull(k)%Fields(this%stepsFull(k)%getVarIdx(varNamesFull(2)))%DATA(:, :, 1) = 0.0D0

          ! this%stepsFull(k)%Fields(w_idx)%data(:, 1:num_icell, 1) = &
          !    this%stepsFull(k - 1)%Fields(w_idx)%data(:, 1:num_icell, 1) &
          !    + dk(k - 1)* &
          !    (this%stepsFull(k - 1)%Fields(tw_idx)%data(:, 1:num_icell, 1) &
          !     - this%stepsFull(k - 1)%Fields(tw_idxt)%data(:, 1:num_icell, 1))

          ! this%stepsFull(k)%Fields(this%stepsFull(k)%getVarIdx(varNamesFull(2)))%data = 0.0D0

          this%stepsHalf(k)%Fields(this%stepsHalf(k)%getVarIdx(TenVarNamesHalf(5)))%DATA(:, 1:num_icell, 1) = &
            this%stepsHalf(k)%Fields(this%stepsHalf(k)%getVarIdx(ParVarNamesHalft(5)))%DATA(:, 1:num_icell, 1)

          this%stepsHalf(k)%Fields(this%stepsHalf(k)%getVarIdx(varNamesHalf(5)))%DATA(:, 1:num_icell, 1) = &
            this%stepsHalf(1)%Fields(this%stepsHalf(1)%getVarIdx(varNamesHalf(5)))%DATA(:, 1:num_icell, 1) &
            + dk(k - 1) * &
            this%stepsHalf(k - 1)%Fields(this%stepsHalf(k - 1)%getVarIdx(TenvarNamesHalf(5)))%DATA(:, 1:num_icell, 1) ! &
          ! + this%stepsHalf(k)%Fields(pres_idx)%data(:, 1:num_icell, 1) &
          ! - 10.0D0**(-5.0D0*k_d/(1.0D0 - k_d)) &
          ! *(this%stepsHalf(k)%Fields(rho_idxt)%data(:, 1:num_icell, 1) &
          !   *r_d**this%stepsHalf(k)%Fields(theta_idxt)%data(:, 1:num_icell, 1)) &
          ! **(1/(1.0D0 - k_d))

          CALL this%stepsHalf(k)%sg%ExchangeMatOnHalo2D(this%stepsHalf(k)%sg%vLevel, &
                                                        this%stepsHalf(k)%Fields(this%stepsHalf(k)%getVarIdx(varNamesHalf(5)))%DATA(:, :, 1))
          CALL this%stepsHalf(k)%sg%ExchangeMatOnHalo2D(this%stepsHalf(k)%sg%vLevel, &
                                                        this%stepsHalf(k)%Fields(this%stepsHalf(k)%getVarIdx(TenVarNamesHalf(5)))%DATA(:, :, 1))

          ! CALL this%RossbyHalf%righthandsHalf(t0 + (it - 2)*this%dt + dk(k - 1)*this%dt, this%stepsHalf(k))
          ! CALL this%RossbyFull%righthandsFull(t0 + (it - 2)*this%dt + dk(k - 1)*this%dt, this%stepsFull(k))
          CALL this%RossbyHalf%righthandsHalf(t0, this%stepsHalf(k))
          CALL this%RossbyFull%righthandsFull(t0, this%stepsFull(k))

          PRINT *, 'dt in time step ', it - 1, ' rk step ', k, ' is ', t0 + (it - 2) * this%dt + dk(k - 1) * this%dt

          this%filenameFull = TRIM('Test_rk4_Full_d1_'//TRIM(INT_TO_STRING(it - 1))//'_'//TRIM(INT_TO_STRING(k))//'_before_Poisson')
          this%filenameHalf = TRIM('Test_rk4_Half_d1_'//TRIM(INT_TO_STRING(it - 1))//'_'//TRIM(INT_TO_STRING(k))//'_before_Poisson')

          ! CALL Output_NC_State_AV(this%stepsFull(k), this%outputDir, this%filenameFull)
          ! CALL Output_NC_State_AV(this%stepsHalf(k), this%outputDir, this%filenameHalf)
          CALL CPU_TIME(start_time)

          vortemp = this%stepsHalf(k)%Fields(vor_idx)%DATA(:, :, 1)
          divtemp = this%stepsHalf(k)%Fields(div_idx)%DATA(:, :, 1)

          DO ik = 1, XHalf%sg%vLevel
            vortemp(ik, :) = vortemp(ik, :) * ratio(ik)
            divtemp(ik, :) = divtemp(ik, :) * ratio(ik)
          END DO

          IF (ALLOCATED(this%bdy_idx)) THEN
            vortemp(:, this%bdy_idx) = this%stepsHalf(k)%Fields(psi_idxt)%DATA(:, this%bdy_idx, 1)
            divtemp(:, this%bdy_idx) = this%stepsHalf(k)%Fields(chi_idxt)%DATA(:, this%bdy_idx, 1)
          END IF

          CALL this%stepsHalf(k)%addVar('vortemp')
          this%stepsHalf(k)%fields(this%stepsHalf(k)%getVarIdx('vortemp'))%DATA(:, :, 1) = vortemp

          CALL this%stepsHalf(k)%sg%ExchangeMatOnHalo2D(this%stepsHalf(k)%sg%vLevel, vortemp)
          CALL this%stepsHalf(k)%sg%ExchangeMatOnHalo2D(this%stepsHalf(k)%sg%vLevel, divtemp)

          this%stepsHalf(k)%Fields(chi_idx)%DATA(:, :, 1) = this%stepsHalf(k - 1)%Fields(chi_idx)%DATA(:, :, 1)
          this%stepsHalf(k)%Fields(psi_idx)%DATA(:, :, 1) = this%stepsHalf(k - 1)%Fields(psi_idx)%DATA(:, :, 1)

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

          CALL this%nwp%righthands(this%stepsFull(k), this%stepsHalf(k))

          CALL UpdateBdy_t(this%stepsFull(k), this%stepsFull(k), 'FULL')
          CALL UpdateBdy_t(this%stepsHalf(k), this%stepsHalf(k), 'HALF')
          CALL this%tendency_ana_minus_simu(this%stepsFull(k), this%stepsHalf(k))
          ! CALL this%nwp%UpdateBdy(this%stepsFull(k), XFull, i_x - 1, 'FULL')
          ! CALL this%nwp%UpdateBdy(this%stepsHalf(k), XHalf, i_x - 1, 'HALF')

          this%filenameFull = TRIM('Test_rk4_Full_d1_'//TRIM(INT_TO_STRING(it - 1))//'_'//TRIM(INT_TO_STRING(k))//'_after_Poisson')
          this%filenameHalf = TRIM('Test_rk4_Half_d1_'//TRIM(INT_TO_STRING(it - 1))//'_'//TRIM(INT_TO_STRING(k))//'_after_Poisson')

          ! CALL Output_NC_State_AV(this%stepsFull(k), this%outputDir, this%filenameFull)
          ! CALL Output_NC_State_AV(this%stepsHalf(k), this%outputDir, this%filenameHalf)

        END DO

        DO iv = 1, 4
          XtHalf%Fields(XtHalf%getVarIdx(varNamesHalf(iv)))%DATA(:, 1:num_icell, 1) = &
            this%stepsHalf(1)%Fields(this%stepsHalf(1)%getVarIdx(varNamesHalf(iv)))%DATA(:, 1:num_icell, 1) + &
            (this%stepsHalf(1)%Fields(this%stepsHalf(1)%getVarIdx(mvarNamesHalf(iv)))%DATA(:, 1:num_icell, 1) + &
             this%stepsHalf(2)%Fields(this%stepsHalf(2)%getVarIdx(mvarNamesHalf(iv)))%DATA(:, 1:num_icell, 1) * 2.0D0 + &
             this%stepsHalf(3)%Fields(this%stepsHalf(3)%getVarIdx(mvarNamesHalf(iv)))%DATA(:, 1:num_icell, 1) * 2.0D0 + &
             this%stepsHalf(4)%Fields(this%stepsHalf(4)%getVarIdx(mvarNamesHalf(iv)))%DATA(:, 1:num_icell, 1)) &
            * this%dt / 6.0D0

          CALL XtHalf%sg%ExchangeMatOnHalo2D(XtHalf%sg%vLevel, &
                                             XtHalf%Fields(XtHalf%getVarIdx(varNamesHalf(iv)))%DATA(:, :, 1))
        END DO
        DO iv = 1, 2
          XtFull%Fields(XtFull%getVarIdx(varNamesFull(iv)))%DATA(:, 1:num_icell, 1) = &
            this%stepsFull(1)%Fields(this%stepsFull(1)%getVarIdx(varNamesFull(iv)))%DATA(:, 1:num_icell, 1) + &
            (this%stepsFull(1)%Fields(this%stepsFull(1)%getVarIdx(mvarNamesFull(iv)))%DATA(:, 1:num_icell, 1) + &
             this%stepsFull(2)%Fields(this%stepsFull(2)%getVarIdx(mvarNamesFull(iv)))%DATA(:, 1:num_icell, 1) * 2.0D0 + &
             this%stepsFull(3)%Fields(this%stepsFull(3)%getVarIdx(mvarNamesFull(iv)))%DATA(:, 1:num_icell, 1) * 2.0D0 + &
             this%stepsFull(4)%Fields(this%stepsFull(4)%getVarIdx(mvarNamesFull(iv)))%DATA(:, 1:num_icell, 1)) &
            * this%dt / 6.0D0

          CALL XtFull%sg%ExchangeMatOnHalo2D(XtFull%sg%vLevel, &
                                             XtFull%Fields(XtFull%getVarIdx(varNamesFull(iv)))%DATA(:, :, 1))
        END DO

        ! XtFull%Fields(w_idx)%data(:, 1:num_icell, 1) = this%stepsFull(1)%Fields(w_idx)%data(:, 1:num_icell, 1) + &
        !                                                (this%stepsFull(1)%Fields(tw_idx)%data(:, 1:num_icell, 1) &
        !                                                 + this%stepsFull(2)%Fields(tw_idx)%data(:, 1:num_icell, 1)*2.0D0 &
        !                                                 + this%stepsFull(3)%Fields(tw_idx)%data(:, 1:num_icell, 1)*2.0D0 &
        !                                                + this%stepsFull(4)%Fields(tw_idx)%data(:, 1:num_icell, 1))*this%dt/6.0D0
        XtFull%Fields(XtFull%getVarIdx(varNamesFull(2)))%DATA = 0.0D0

        XtHalf%Fields(XtHalf%getVarIdx(varNamesHalf(5)))%DATA(:, 1:num_icell, 1) = &
          this%stepsHalf(1)%Fields(XtHalf%getVarIdx(varNamesHalf(5)))%DATA(:, 1:num_icell, 1) + &
          (this%stepsHalf(1)%Fields(XtHalf%getVarIdx(TenvarNamesHalf(5)))%DATA(:, 1:num_icell, 1) &
           + this%stepsHalf(2)%Fields(XtHalf%getVarIdx(TenvarNamesHalf(5)))%DATA(:, 1:num_icell, 1) * 2.0D0 &
           + this%stepsHalf(3)%Fields(XtHalf%getVarIdx(TenvarNamesHalf(5)))%DATA(:, 1:num_icell, 1) * 2.0D0 &
           + this%stepsHalf(4)%Fields(XtHalf%getVarIdx(TenvarNamesHalf(5)))%DATA(:, 1:num_icell, 1)) * this%dt / 6.0D0

        CALL XtHalf%sg%ExchangeMatOnHalo2D(XtHalf%sg%vLevel, XtHalf%Fields(pres_idx)%DATA(:, :, 1))

        this%filenameFull = TRIM('Test_rk4_Full_d1_Xt_'//TRIM(INT_TO_STRING(it))//'_before_Poisson')
        this%filenameHalf = TRIM('Test_rk4_Half_d1_Xt_'//TRIM(INT_TO_STRING(it))//'_before_Poisson')

        ! CALL Output_NC_State_AV(XtFull, this%outputDir, this%filenameFull)
        ! CALL Output_NC_State_AV(XtHalf, this%outputDir, this%filenameHalf)

        CALL CPU_TIME(start_time)

        vortemp = XtHalf%Fields(vor_idx)%DATA(:, :, 1)
        divtemp = XtHalf%Fields(div_idx)%DATA(:, :, 1)

        DO ik = 1, XHalf%sg%vLevel
          vortemp(ik, :) = vortemp(ik, :) * ratio(ik)
          divtemp(ik, :) = divtemp(ik, :) * ratio(ik)
        END DO

        IF (ALLOCATED(this%bdy_idx)) THEN
          vortemp(:, this%bdy_idx) = this%stepsHalf(k - 1)%Fields(psi_idxt)%DATA(:, this%bdy_idx, 1)
          divtemp(:, this%bdy_idx) = this%stepsHalf(k - 1)%Fields(chi_idxt)%DATA(:, this%bdy_idx, 1)
        END IF

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

        this%filenameFull = TRIM('Test_rk4_Full_d1_Xt_'//TRIM(INT_TO_STRING(it))//'_after_Poisson')
        this%filenameHalf = TRIM('Test_rk4_Half_d1_Xt_'//TRIM(INT_TO_STRING(it))//'_after_Poisson')

        ! CALL Output_NC_State_AV(XtFull, this%outputDir, this%filenameFull)
        ! CALL Output_NC_State_AV(XtHalf, this%outputDir, this%filenameHalf)

        ! CALL this%RossbyHalf%righthandsHalf(t0 + (it - 1)*this%dt, XtHalf)
        ! CALL this%RossbyFull%righthandsFull(t0 + (it - 1)*this%dt, XtFull)
        CALL this%RossbyHalf%righthandsHalf(t0, XtHalf)
        CALL this%RossbyFull%righthandsFull(t0, XtFull)
        PRINT *, 'begin write error'
        DO iv = 1, 4
          error = 0.0D0
          error1 = 0.0D0

          DO k = 1, XtHalf%sg%vLevel
            DO i = 1, XtHalf%sg%num_icell
              IF (XtHalf%sg%bdy_type(i) .LT. -1) THEN
                error1 = DABS(XtHalf%Fields(XtHalf%getVarIdx(varNamesHalf(iv)))%DATA(k, i, 1) - &
                              XtHalf%Fields(XtHalf%getVarIdx(varNamesHalft(iv)))%DATA(k, i, 1))
                IF (error .LE. error1) THEN
                  error = error1
                END IF
              END IF
            END DO
          END DO
          OPEN (unit=10, file=TRIM('error_Half_'//TRIM(INT_TO_STRING(iv))//'.txt'), &
                status='unknown', action='write', access='append', form='formatted')
          WRITE (10, *) error
          CLOSE (10)
        END DO

        DO iv = 1, 1
          error = 0.0D0
          error1 = 0.0D0

          DO k = 1, XtFull%sg%vLevel
            DO i = 1, XtFull%sg%num_icell
              IF (XtFull%sg%bdy_type(i) .LT. -1) THEN
                error1 = DABS(XtFull%Fields(XtFull%getVarIdx(varNamesFull(iv)))%DATA(k, i, 1) - &
                              XtFull%Fields(XtFull%getVarIdx(varNamesFullt(iv)))%DATA(k, i, 1))
                IF (error .LE. error1) THEN
                  error = error1
                END IF
              END IF
            END DO
          END DO
          OPEN (unit=10, file=TRIM('error_Full_'//TRIM(INT_TO_STRING(iv))//'.txt'), &
                status='unknown', action='write', access='append', form='formatted')
          WRITE (10, *) error
          CLOSE (10)
        END DO
        PRINT *, 'end write error'

        PRINT *, 'dt in Xt_', it, ' is ', t0 + (it - 1) * this%dt
        ! CALL UpdateBdy_t(XtFull, XFull, i_x - 1, 'FULL')
        ! CALL UpdateBdy_t(XtHalf, XHalf, i_x - 1, 'HALF')

        ! Fill the integrated vars to X Fields
        IF (MOD(it, this%kt) .EQ. 1) THEN
          DO iv = 1, SIZE(varNamesHalf)
            XHalf%Fields(XHalf%getVarIdx(varNamesHalf(iv)))%DATA(:, :, i_x) = XtHalf%Fields(XHalf%getVarIdx(varNamesHalf(iv)))%DATA(:, :, 1)
            XHalf%Fields(XHalf%getVarIdx(varNamesHalft(iv)))%DATA(:, :, i_x) = XtHalf%Fields(XHalf%getVarIdx(varNamesHalft(iv)))%DATA(:, :, 1)
          END DO
          DO iv = 1, SIZE(varNamesFull)
            XFull%Fields(XFull%getVarIdx(varNamesFull(iv)))%DATA(:, :, i_x) = XtFull%Fields(XFull%getVarIdx(varNamesFull(iv)))%DATA(:, :, 1)
            XFull%Fields(XFull%getVarIdx(varNamesFullt(iv)))%DATA(:, :, i_x) = XtFull%Fields(XFull%getVarIdx(varNamesFullt(iv)))%DATA(:, :, 1)
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
    IF (ASSOCIATED(this%geometryHalf)) NULLIFY (this%geometryHalf)
    IF (ASSOCIATED(this%geometryFull)) NULLIFY (this%geometryFull)

  END SUBROUTINE destructor
END MODULE rk4test_m
