MODULE ObsConvention_m2
  USE kinds_m, ONLY: i_kind, r_kind
  USE ObsBase_m2, ONLY: ObsBase_t2
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE parameters_m, ONLY: missing, invalid, degree2radian, dry_air_gas_const
  USE conversions_m, ONLY: wind_to_uv, Td_to_qvapor, hydroHeightUp, hydroHeightDown
  USE YAMLRead_m
  ! Ingest subroutine needs:
  USE geoTools_m, ONLY: GeoBox_t
  ! QC subroutine needs:
  USE domainCheck_m, ONLY: domainCheck_t
  ! Thinning subroutine needs:
  USE mpObs_m, ONLY: mpObs_t
  USE ObsField_m, ONLY: ObsField_t
  IMPLICIT NONE

  TYPE, EXTENDS(ObsBase_t2) :: ObsConvention_t2
    INTEGER(i_kind) :: numSt, numLevel ! number of stations, variables and levels
    CHARACTER(LEN=20), ALLOCATABLE :: nameSt(:) ! an array for station names (or ID)
    REAL(r_kind), ALLOCATABLE, DIMENSION(:, :) :: obs4d ! an array for station time and location (time, lat, lon, station height)
    REAL(r_kind), ALLOCATABLE, DIMENSION(:, :) :: obsHgt ! observation height for temp and sound observations
    REAL(r_kind), ALLOCATABLE, DIMENSION(:, :, :) :: obsData, obsErrs
    REAL(r_kind), ALLOCATABLE, DIMENSION(:) :: qcThreshold

  CONTAINS
    PROCEDURE, PUBLIC :: obsIngest => ingest_conv
    PROCEDURE, PUBLIC :: obsQC_conv => qc_conv
    PROCEDURE, PUBLIC :: obsThinning => thinning_conv
    PROCEDURE, PUBLIC :: obsDeallocate => deallocate_conv
  END TYPE ObsConvention_t2

  PRIVATE :: ingest_conv, ingest_conv_each, qc_conv, qc_synop, qc_ship, qc_buoy, qc_metar, qc_temp, qc_profl

CONTAINS
  SUBROUTINE ingest_conv(this, X)
    IMPLICIT NONE

    CLASS(ObsConvention_t2) :: this
    TYPE(State_t) :: X

    ! temporal data
    CHARACTER(LEN=20), ALLOCATABLE :: t_nameSt(:)
    REAL(r_kind), ALLOCATABLE :: t_obs4d(:, :), t_obsData(:, :, :), t_obsErrs(:, :, :), t_obsHgt(:, :)

    INTEGER(i_kind) :: i, j, status, i_start, i_end, info_num
    INTEGER(i_kind) :: num_st_vec(this%numFiles + 1), num_level_vec(this%numFiles)
    TYPE(ObsConvention_t2) :: local_conv(this%numFiles)

    REAL(r_kind) :: missingValue = -9999.0D0
    REAL(r_kind), PARAMETER :: PI = 3.1415926535897932
    INTEGER(i_kind) :: wd_idx, ws_idx, dp_idx
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: wd_value, ws_value, wd_err, ws_err

    ! local varaibles for checking obs inside processor's domain
    TYPE(GeoBox_t) :: geoBox
    ! check if keeping this observation or not
    LOGICAL, ALLOCATABLE :: index_unique_judge(:)
    INTEGER(i_kind), ALLOCATABLE :: index_unique_tmp(:), index_unique(:)
    INTEGER(i_kind) :: num_unique

    ! maximum height
    ! REAL(r_kind) :: ztop
    ! black list related variables
    CHARACTER(LEN=1024) :: blackList
    CHARACTER(LEN=5), ALLOCATABLE :: blackListNames(:)
    LOGICAL :: blackFlag, blackListExist
    INTEGER(i_kind) :: blackLen

    num_st_vec(1) = 0
    ! read data from nc file
    DO i = 1, this%numFiles
      CALL ingest_conv_each(this, local_conv(i), i)
      num_st_vec(i + 1) = local_conv(i)%numSt
      num_level_vec(i) = local_conv(i)%numLevel
    END DO

    !put conventional observation data into one object
    this%numSt = SUM(num_st_vec)
    this%numLevel = MAXVAL(num_level_vec)

    ! get information number, 4 for single level observation; 5 for mutiple levels observations
    info_num = UBOUND(local_conv(1)%obs4d, 2)

    ! allocate variables
    ALLOCATE (t_nameSt(this%numSt), STAT=status)
    ALLOCATE (t_obs4d(this%numSt, info_num), STAT=status)
    ALLOCATE (t_obsData(this%numSt, this%numVars, this%numLevel), STAT=status)
    ALLOCATE (t_obsErrs(this%numSt, this%numVars, this%numLevel), STAT=status)

    ! get observation height for temp and sound observation
    IF (info_num .EQ. 5) ALLOCATE (t_obsHgt(this%numSt, this%numLevel), STAT=status)

    DO i = 1, this%numFiles
      i_start = SUM(num_st_vec(1:i)) + 1
      i_end = SUM(num_st_vec(1:i + 1))
      t_nameSt(i_start:i_end) = local_conv(i)%nameSt
      t_obs4d(i_start:i_end, :) = local_conv(i)%obs4d
      t_obsData(i_start:i_end, :, 1:num_level_vec(i)) = local_conv(i)%obsData
      t_obsErrs(i_start:i_end, :, 1:num_level_vec(i)) = local_conv(i)%obsErrs
      ! get observation height
      IF (info_num .EQ. 5) t_obsHgt(i_start:i_end, 1:num_level_vec(i)) = local_conv(i)%obsHgt
    END DO

    ! convert wind direction and speed to u, w wind
    ALLOCATE (wd_value(this%numSt, this%numLevel), STAT=status)
    ALLOCATE (ws_value(this%numSt, this%numLevel), STAT=status)
    ALLOCATE (wd_err(this%numSt, this%numLevel), STAT=status)
    ALLOCATE (ws_err(this%numSt, this%numLevel), STAT=status)
    wd_idx = FINDLOC(this%varNames, 'wd', 1)
    ws_idx = FINDLOC(this%varNames, 'ws', 1)
    IF (wd_idx .NE. 0 .AND. ws_idx .NE. 0) THEN
      wd_value = t_obsData(:, wd_idx, :)
      ws_value = t_obsData(:, ws_idx, :)
      CALL wind_to_uv(ws_value, wd_value, t_obsData(:, wd_idx, :), t_obsData(:, ws_idx, :))
      wd_err = t_obsErrs(:, wd_idx, :)
      ws_err = t_obsErrs(:, ws_idx, :)
      CALL wind_to_uv(ws_err, wd_err, t_obsErrs(:, wd_idx, :), t_obsErrs(:, ws_idx, :))
      ! change variable name
      this%varNames(wd_idx) = 'uwnd'
      this%varNames(ws_idx) = 'vwnd'
    END IF

    ! Make invalid data
    wd_idx = FINDLOC(this%varNames, 'uwnd', 1)
    ws_idx = FINDLOC(this%varNames, 'vwnd', 1)
    WHERE (ABS(t_obsData(:, wd_idx, :)) .GT. ABS(missing) .OR. &
           ABS(t_obsData(:, ws_idx, :)) .GT. ABS(missing) .OR. &
           t_obsData(:, ws_idx, :) .LT. 0.0D0 .OR. &
           t_obsData(:, ws_idx, :) .GT. 500.0D0) t_obsData(:, wd_idx, :) = invalid

    WHERE (ABS(t_obsData(:, wd_idx, :)) .GT. ABS(missing) .OR. &
           ABS(t_obsData(:, ws_idx, :)) .GT. ABS(missing) .OR. &
           t_obsData(:, ws_idx, :) .LT. 0.0D0 .OR. &
           t_obsData(:, ws_idx, :) .GT. 500.0D0) t_obsData(:, ws_idx, :) = invalid

    ! check missing data
    WHERE (ABS(t_obsData - missingValue) .LT. 0.01D0) t_obsData = missing

    ! convert degree to radian for lat and lon
    t_obs4d(:, 2:3) = t_obs4d(:, 2:3) * degree2radian
    geoBox = GeoBox_t(X%sg%maxLatGlob, X%sg%minLatGlob, X%sg%maxLonGlob, X%sg%minLonGlob)

    ! allocate check array
    ALLOCATE (index_unique_judge(this%numSt), STAT=status)
    ALLOCATE (index_unique_tmp(this%numSt), STAT=status)
    index_unique_judge = .TRUE.
    index_unique_tmp = 0

    WHERE (.NOT. geoBox%inBoxElem(t_obs4d(:, 2), t_obs4d(:, 3), X%sg%maxLatGlob, X%sg%minLatGlob, &
                                  X%sg%maxLonGlob, X%sg%minLonGlob)) index_unique_judge = .FALSE.

    ! remove duplicated observations
    DO i = 1, this%numSt
      IF (index_unique_judge(i)) THEN
        DO j = i + 1, this%numSt
          ! do not check false observation
          IF (.NOT. index_unique_judge(j)) CYCLE

          ! compare ID and time
          IF (t_nameSt(i) .EQ. t_nameSt(j) .AND. t_obs4d(i, 1) .EQ. t_obs4d(j, 1)) index_unique_judge(j) = .FALSE.
        END DO
      END IF
    END DO

    ! get black list names from yaml file
    status = yaml_get_var(TRIM(this%configFile), TRIM(this%obsType)//"_info", 'blackExist', blackFlag)
    IF (blackFlag) status = yaml_get_var(TRIM(this%configFile), 'Verify', TRIM(this%obsType)//'_BlackListAdd', blackList)
    INQUIRE (FILE=TRIM(blackList), EXIST=blackListExist)
    IF (blackListExist) status = yaml_get_var(TRIM(blackList), 'sites', 'name', blackListNames)
    blackLen = UBOUND(blackListNames, 1)  ! length of black list
    IF (blackListExist .AND. blackLen .GT. 0) THEN
      DO i = 1, this%numSt
        ! jump 'false' observation
        IF (.NOT. index_unique_judge(i)) CYCLE

        ! remove black station
        DO j = 1, blackLen
          IF (TRIM(t_nameSt(i)) .EQ. TRIM(blackListNames(j))) index_unique_judge(i) = .FALSE.
        END DO
      END DO
    END IF

    ! assigne index_unique_tmp elements to be i for those 'true' observations
    FORALL (i=1:this%numSt, index_unique_judge(i)) index_unique_tmp(i) = i

    ! remove false observations
    num_unique = COUNT(index_unique_tmp .NE. 0)
    ALLOCATE (index_unique(num_unique), STAT=status)
    index_unique = PACK(index_unique_tmp, index_unique_tmp .NE. 0)

    ALLOCATE (this%nameSt(num_unique), STAT=status)
    ALLOCATE (this%obs4d(num_unique, info_num), STAT=status)
    ALLOCATE (this%obsData(num_unique, this%numVars, this%numLevel), STAT=status)
    ALLOCATE (this%obsErrs(num_unique, this%numVars, this%numLevel), STAT=status)
    IF (info_num .EQ. 5) THEN
      ALLOCATE (this%obsHgt(num_unique, this%numLevel), STAT=status)
      this%obsHgt = t_obsHgt(index_unique, :)
    END IF
    this%numSt = num_unique

    this%nameSt = t_nameSt(index_unique)
    this%obs4d = t_obs4d(index_unique, :)
    this%obsData = t_obsData(index_unique, :, :)
    this%obsErrs = t_obsErrs(index_unique, :, :)

    IF (ALLOCATED(t_nameSt)) DEALLOCATE (t_nameSt)
    IF (ALLOCATED(t_obs4d)) DEALLOCATE (t_obs4d)
    IF (ALLOCATED(t_obsData)) DEALLOCATE (t_obsData)
    IF (ALLOCATED(t_obsErrs)) DEALLOCATE (t_obsErrs)
    IF (ALLOCATED(t_obsHgt)) DEALLOCATE (t_obsHgt)
    IF (ALLOCATED(wd_value)) DEALLOCATE (wd_value)
    IF (ALLOCATED(ws_value)) DEALLOCATE (ws_value)
    IF (ALLOCATED(wd_err)) DEALLOCATE (wd_err)
    IF (ALLOCATED(ws_err)) DEALLOCATE (ws_err)
    IF (ALLOCATED(index_unique_judge)) DEALLOCATE (index_unique_judge)
    IF (ALLOCATED(index_unique_tmp)) DEALLOCATE (index_unique_tmp)
    IF (ALLOCATED(index_unique)) DEALLOCATE (index_unique)
    IF (ALLOCATED(blackListNames)) DEALLOCATE (blackListNames)
  END SUBROUTINE ingest_conv

  SUBROUTINE qc_conv(this, state)
    IMPLICIT NONE
    CLASS(ObsConvention_t2) :: this
    TYPE(State_t), INTENT(IN) :: state

    INTEGER(i_kind) :: i, dp_idx
    IF (.NOT. ALLOCATED(this%qcThreshold)) ALLOCATE (this%qcThreshold(this%numVars))
    DO i = 1, this%numVars
      SELECT CASE (TRIM(this%varNames(i)))
      CASE ('t')
        this%qcThreshold(i) = 20.0D0
      CASE ('qv')
        this%qcThreshold(i) = 0.02D0
      CASE ('uwnd')
        this%qcThreshold(i) = 20.0D0
      CASE ('vwnd')
        this%qcThreshold(i) = 20.0D0
      CASE ('p')
        this%qcThreshold(i) = 6000.0D0
      CASE DEFAULT
        this%qcThreshold(i) = 500.0D0
      END SELECT
    END DO

    SELECT CASE (this%obsType)
    CASE ("SYNOP")
      CALL qc_synop(this)
    CASE ("SHIP")
      CALL qc_ship(this)
    CASE ("BUOY")
      CALL qc_buoy(this)
    CASE ("METAR")
      CALL qc_metar(this)
    CASE ("TEMP")
      CALL qc_temp(this, state)
    CASE ("PROFL")
      CALL qc_profl(this)
    END SELECT

    ! Convert dp to qv
    dp_idx = FINDLOC(this%varNames, 'dp', 1)
    IF (dp_idx .NE. 0) THEN
      this%obsData(:, dp_idx, :) = Td_to_qvapor(this%obsData(:, dp_idx, :), this%obsData(:, FINDLOC(this%varNames, 'p', 1), :))
      this%varNames(dp_idx) = 'qv'
    END IF
  END SUBROUTINE qc_conv

  SUBROUTINE thinning_conv(this, state, thinObs, mpObs, useBkg)
    IMPLICIT NONE
    CLASS(ObsConvention_t2) :: this
    TYPE(State_t), INTENT(IN) :: state
    TYPE(ObsSet_t), INTENT(INOUT) :: thinObs
    TYPE(mpObs_t), INTENT(IN) :: mpObs
    LOGICAL, INTENT(IN) :: useBkg

    ! local variables
    INTEGER(i_kind) :: i, j, k, io, iv, status, num_obs, i_start, i_end, t_id, p_id, q_id
    INTEGER(i_kind) :: var_state_id, h_idx, g_idx, t_idx, num_thin_obs, num_state_var, i_state
    REAL(r_kind) :: bkg_at_obs, gaussian
    CHARACTER(LEN=20) :: varname
    REAL(r_kind), PARAMETER :: sigma2 = 4.0D0

    ! Thinning weights, observation weights, weighted obs, weighted err
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: weights, wghtobs, wghterr

    ! Transpose and variables
    REAL(r_kind), ALLOCATABLE, DIMENSION(:, :) :: trans_data, trans_err, trans_pos, valided_obsData, valided_obsErrs, valided_pos
    INTEGER(i_kind), ALLOCATABLE, DIMENSION(:) :: trans_time, valided_time
    CHARACTER(LEN=20), ALLOCATABLE, DIMENSION(:) :: trans_ID, valided_ID, state_var_list

    ! valid related variables
    INTEGER(i_kind) :: num_valid, nhstencil
    INTEGER(i_kind), ALLOCATABLE :: all_idx(:), valided_idx(:), mask_valid(:), topo(:)
    INTEGER(i_kind), ALLOCATABLE ::  idx_grd_valid(:, :), idx_hgt_valid(:, :, :), idx_t_valid(:, :)
    REAL(r_kind), ALLOCATABLE ::  coe_grd_valid(:, :), coe_hgt_valid(:, :, :), coe_t_valid(:, :)
    REAL(r_kind), ALLOCATABLE :: forward(:, :, :)

    ! domain check variables
    CHARACTER(LEN=20) :: d_obsType
    TYPE(domainCheck_t) :: domain
    ! An array for saving if this variable exist in the state or not
    LOGICAL :: var_exist(this%numVars)

    IF (.NOT. state%sg%isActiveProc()) RETURN

    ! Output observation data information (minimal and maximal values)
    DO i = 1, this%numVars
      varname = this%varNames(i)
      var_exist(i) = state%getVarIdx(TRIM(varname)) .NE. 0
      ! Check observation variable exist in state variables or not
      IF (.NOT. var_exist(i)) CYCLE

      WRITE (*, 1) TRIM(varname), MAXVAL(this%obsData(:, i, :), this%obsData(:, i, :) .LT. 1.0D8), &
        MINVAL(this%obsData(:, i, :), this%obsData(:, i, :) .GT. -1.0D8), &
        state%sg%mpddInfo_sg%myrank, TRIM(this%obsType), state%sg%gLevel
1     FORMAT('MAX/MIN obs: ', A, ' before thinning: ', 2D12.4, ' at proc: ', I1, ' Type: ', A, ' Glevel: ', I2)
    END DO

    num_state_var = COUNT(var_exist)
    ALLOCATE (state_var_list(num_state_var))
    i = 0
    DO iv = 1, this%numVars
      IF (.NOT. var_exist(iv)) CYCLE
      i = i + 1
      state_var_list(i) = TRIM(this%varNames(i))
    END DO

    ! Intialialize thinObs
    thinObs = ObsSet_t(this%configFile, mpObs)

    ! Allocate memory for thinning weights
    ALLOCATE ( &
      weights(state%sg%vLevel, state%sg%num_cell, state%sg%tSlots, num_state_var), &
      wghtobs(state%sg%vLevel, state%sg%num_cell, state%sg%tSlots, num_state_var), &
      wghterr(state%sg%vLevel, state%sg%num_cell, state%sg%tSlots, num_state_var), &
      STAT=status)
    IF (status .NE. 0) THEN
      WRITE (*, *) 'Observation thinning: cannot allocate memory for weights'
      STOP
    END IF

    ! Initializing values:
    weights = 0.0D0
    wghtobs = 0.0D0
    wghterr = 0.0D0

    ! Convert observation data from numSt*numVars*numLevel to numObs*numVars
    IF (UBOUND(this%obs4d, 2) .EQ. 4) THEN
      num_obs = this%numSt
      ALLOCATE (trans_time(num_obs), trans_ID(num_obs), trans_pos(3, num_obs), &
                trans_data(num_obs, this%numVars), trans_err(num_obs, this%numVars), STAT=status)
      trans_time = INT(this%obs4d(:, 1))
      trans_ID = this%nameSt
      trans_pos = TRANSPOSE(this%obs4d(:, 2:4))
      trans_data = this%obsData(:, :, 1)
      trans_err = this%obsErrs(:, :, 1)
    ELSE ! for sound and temp observation data
      num_obs = INT(SUM(this%obs4d(:, 5)))
      ALLOCATE (trans_time(num_obs), trans_ID(num_obs), trans_pos(3, num_obs), &
                trans_data(num_obs, this%numVars), trans_err(num_obs, this%numVars), STAT=status)
      i_start = 1
      i_end = INT(this%obs4d(1, 5))
      DO i = 1, this%numSt
        IF (i > 1) i_start = 1 + INT(SUM(this%obs4d(1:i - 1, 5)))
        i_end = INT(SUM(this%obs4d(1:i, 5)))
        trans_time(i_start:i_end) = INT(this%obs4d(i, 1))
        trans_ID(i_start:i_end) = this%nameSt(i)
        trans_pos(1, i_start:i_end) = this%obs4d(i, 2)
        trans_pos(2, i_start:i_end) = this%obs4d(i, 3)
        trans_pos(3, i_start:i_end) = this%obsHgt(i, 1:INT(this%obs4d(i, 5)))
        trans_data(i_start:i_end, :) = TRANSPOSE(this%obsData(i, :, 1:INT(this%obs4d(i, 5))))
        trans_err(i_start:i_end, :) = TRANSPOSE(this%obsErrs(i, :, 1:INT(this%obs4d(i, 5))))
      END DO
    END IF

    ! Domain check
    ! Add d for sound observation
    d_obsType = TRIM(this%obsType)
    IF (d_obsType .EQ. 'TEMP') d_obsType = 'SOUND_d'
    ALLOCATE (mask_valid(num_obs), STAT=status)
    CALL domain%validation(state%sg, num_obs, trans_pos(1:2, :) * degree2radian, trans_pos(3, :), &
                           trans_time, num_valid, mask_valid, idx_grd_valid, &
                           idx_hgt_valid, idx_t_valid, &
                           coe_grd_valid, coe_hgt_valid, &
                           coe_t_valid, 1, nhstencil, d_obsType)

    ! Allocate valided variables:
    ALLOCATE (valided_time(num_valid), valided_ID(num_valid), &
              valided_pos(3, num_valid), valided_obsData(num_valid, this%numVars), &
              valided_obsErrs(num_valid, this%numVars), &
              STAT=status)
    ALLOCATE (forward(2, nhstencil, 2), STAT=status)

    ! Pass values from observation data to valided data
    ALLOCATE (all_idx(num_obs), valided_idx(num_valid), STAT=status)
    all_idx = [(i, i=1, num_obs)]
    WHERE (mask_valid .LE. 0) all_idx = 0
    valided_idx = PACK(mask_valid, mask_valid .NE. 0)
    valided_time = trans_time(valided_idx)
    valided_ID = trans_ID(valided_idx)
    valided_pos = trans_pos(:, valided_idx)
    valided_obsData = trans_data(valided_idx, :)
    valided_obsErrs = trans_err(valided_idx, :)

    ! Correct the pressure and temperature of the surface observations
    IF (TRIM(this%obsType) .EQ. 'SYNOP' .AND. state%sg%vLevel > 1) THEN
      ALLOCATE (topo(num_valid), STAT=status)
      topo = 0.0D0

      ! Get topography height
      DO i = 1, num_valid
        DO j = 1, nhstencil
          topo(i) = topo(i) + state%sg%topo(idx_grd_valid(j, i)) * coe_grd_valid(j, i)
        END DO
      END DO

      t_id = FINDLOC(this%varNames, 't', 1)
      p_id = FINDLOC(this%varNames, 'p', 1)
      q_id = FINDLOC(this%varNames, 'qv', 1)
      ! Correct temperature
      IF (t_id .NE. 0) valided_obsData(:, t_id) = valided_obsData(:, t_id) - 0.65 * (topo - valided_pos(3, :)) / 100.0D0
      ! Correct pressure
      IF (t_id .NE. 0 .AND. p_id .NE. 0 .AND. q_id .NE. 0) THEN
        valided_obsData(:, p_id) = valided_obsData(:, p_id) * EXP(-(topo - valided_pos(3, :)) * 9.80655D0 / &
                                                                  (dry_air_gas_const * ((valided_obsData(:, t_id) + 0.65 * (topo - valided_pos(3, :)) / 200.0D0) * &
                                                                                        (1.0D0 + 0.608 * valided_obsData(:, q_id)))))
        FORALL (i=1:num_valid, topo(i) - valided_pos(3, i) .GE. 2000.0D0)
          valided_obsData(i, t_id) = 1E9
          valided_obsData(i, p_id) = 1E9
        END FORALL
      END IF
    END IF

    ! Just thinning valid data and variable
    i_state = 0
    DO iv = 1, this%numVars
      ! State does not have this variable
      var_state_id = state%getVarIdx(TRIM(this%varNames(iv)))
      IF (.NOT. var_exist(iv)) CYCLE

      i_state = i_state + 1
      ! For state variable, thinning one by one
      DO io = 1, num_valid
        ! Jump missing variable
        IF (ABS(valided_obsData(io, iv)) .GE. ABS(missing)) CYCLE

        bkg_at_obs = 0.0D0
        IF (useBkg) THEN
          DO i = 1, 2
            DO j = 1, nhstencil
              DO k = 1, 2
                forward(k, j, i) = state%Fields(var_state_id)%DATA( &
                                   idx_hgt_valid(k, j, io), idx_grd_valid(j, io), idx_t_valid(i, io))
                bkg_at_obs = bkg_at_obs + forward(k, j, i) * &
                             coe_hgt_valid(k, j, io) * coe_grd_valid(j, io) * coe_t_valid(i, io)
              END DO
            END DO
          END DO
        END IF

        IF (ABS(bkg_at_obs - valided_obsData(io, iv)) .LE. this%qcThreshold(iv)) THEN
          DO i = 1, 2
            DO j = 1, nhstencil
              DO k = 1, 2
                gaussian = EXP(-(1.0D0 - coe_hgt_valid(k, j, io))**2 / sigma2 &
                               - (1.0D0 - coe_grd_valid(j, io))**2 / sigma2 &
                               - (1.0D0 - coe_t_valid(i, io))**2 / sigma2)
                IF (gaussian .GE. this%correlation_threshold) THEN
                  h_idx = idx_hgt_valid(k, j, io)
                  g_idx = idx_grd_valid(j, io)
                  t_idx = idx_t_valid(i, io)

                  weights(h_idx, g_idx, t_idx, i_state) = &
                    weights(h_idx, g_idx, t_idx, i_state) + gaussian
                  wghtobs(h_idx, g_idx, t_idx, i_state) = &
                    wghtobs(h_idx, g_idx, t_idx, i_state) + gaussian * (valided_obsData(io, iv) - bkg_at_obs)
                  wghterr(h_idx, g_idx, t_idx, i_state) = &
                    wghterr(h_idx, g_idx, t_idx, i_state) + gaussian / valided_obsErrs(io, iv)
                END IF
              END DO
            END DO
          END DO
        END IF
      END DO
    END DO

    ALLOCATE (thinObs%ObsFields(num_state_var))
    DO i = 1, num_state_var
      num_thin_obs = COUNT(weights(:, 1:state%sg%num_icell, :, i) .GT. 0.0)
      thinObs%ObsFields(i) = ObsField_t(this%configFile, mpObs)
      ALLOCATE (thinObs%ObsFields(i)%values(num_thin_obs), &
                thinObs%ObsFields(i)%errors(num_thin_obs), &
                thinObs%ObsFields(i)%idx(num_thin_obs))
      CALL thinObs%ObsFields(i)%Set_Name(state_var_list(i))
      CALL thinObs%ObsFields(i)%Set_ObsType(TRIM(this%obsType))
    END DO

    DO iv = 1, num_state_var
      io = 0
      DO i = 1, state%sg%tSlots
        DO j = 1, state%sg%num_icell
          DO k = 1, state%sg%vLevel
            IF (weights(k, j, i, iv) .GT. 0.0) THEN
              io = io + 1
              ! Assign index values:
              thinObs%ObsFields(iv)%idx(io)%hIdx = j
              thinObs%obsFields(iv)%idx(io)%tIdx = i
              thinObs%ObsFields(iv)%idx(io)%vIdx = k
              bkg_at_obs = 0.0D0
              IF (useBkg) bkg_at_obs = state%Fields(state%getVarIdx(state_var_list(iv)))%DATA(k, j, i)
              thinObs%ObsFields(iv)%values(io) = bkg_at_obs + wghtobs(k, j, i, iv) / weights(k, j, i, iv)
              ! Thinned observation errors:
              thinObs%ObsFields(iv)%errors(io) = 1.0 / weights(k, j, i, iv) / MAX(wghterr(k, j, i, iv), 1.0E-6)
            END IF
          END DO
        END DO
      END DO
    END DO

    WRITE (*, *) "Finished thinning observation: "//TRIM(this%obsType)
  END SUBROUTINE thinning_conv

  SUBROUTINE deallocate_conv(this)
    IMPLICIT NONE
    CLASS(ObsConvention_t2) :: this

    ! local variables
    INTEGER(i_kind) :: status

    ! deallocate variables
    DEALLOCATE (this%nameSt, STAT=status)
    DEALLOCATE (this%obs4d, STAT=status)
    DEALLOCATE (this%obsHgt, STAT=status)
    DEALLOCATE (this%obsData, STAT=status)
    DEALLOCATE (this%obsErrs, STAT=status)
    DEALLOCATE (this%qcThreshold, STAT=status)
  END SUBROUTINE deallocate_conv

  SUBROUTINE ingest_conv_each(this, one_conv, fileID)
    USE mo_netcdf, ONLY: NcDataset, NcDimension, NcVariable
    IMPLICIT NONE

    CLASS(ObsConvention_t2) :: this
    TYPE(ObsConvention_t2), INTENT(INOUT) :: one_conv
    INTEGER(i_kind), INTENT(IN) :: fileID

    ! local variables
    TYPE(NcDataset) :: nc_dataset
    TYPE(NcDimension) :: nc_dim
    TYPE(NcVariable) :: nc_var
    CHARACTER(LEN=20), ALLOCATABLE :: st_in_name(:)
    INTEGER(i_kind) :: i, status, info_num, data_len(2)
    REAL(r_kind), ALLOCATABLE :: data1d(:), data2d(:, :)
    REAL(r_kind) :: data0d

    ! open nc dataset from local file
    nc_dataset = NcDataset(TRIM(this%fileNames(fileID)), "r")
    ! get station number (dimension)
    nc_dim = nc_dataset%getDimension("snum")
    one_conv%numSt = nc_dim%getLength()
    ! get level number for multiple level data (dimension)
    IF (this%obsType .EQ. 'TEMP' .OR. this%obsType .EQ. 'PROFL') THEN
      nc_dim = nc_dataset%getDimension("nlevel")
      one_conv%numLevel = nc_dim%getLength()
      ! for multiple levels, each site has different level numbers
      info_num = 5

      ! for multiple levels, exist observation height
      nc_var = nc_dataset%getVariable("obsh")
      CALL nc_var%getData(one_conv%obsHgt)
      IF (UBOUND(one_conv%obsHgt, 1) /= one_conv%numSt) one_conv%obsHgt = TRANSPOSE(one_conv%obsHgt)
    ELSE
      one_conv%numLevel = 1
      info_num = 4
    END IF

    ! get station id/name
    nc_var = nc_dataset%getVariable("sid")
    CALL nc_var%getData(one_conv%nameSt)

    ALLOCATE (one_conv%obs4d(one_conv%numSt, info_num), STAT=status)
    ALLOCATE (one_conv%obsData(one_conv%numSt, this%numVars, one_conv%numLevel), STAT=status)
    ALLOCATE (one_conv%obsErrs(one_conv%numSt, this%numVars, one_conv%numLevel), STAT=status)

    ! get time, latitude, longitude, height,
    ! and level number for each site when TEMP and PROFL
    ALLOCATE (st_in_name(info_num), STAT=status)
    st_in_name(1) = 'time'
    st_in_name(2) = 'lat'
    st_in_name(3) = 'lon'
    st_in_name(4) = 'height'
    IF (info_num .EQ. 5) st_in_name(5) = 'level'
    DO i = 1, info_num
      ! BUOY station height is set to 0 meter
      IF (this%obsType .EQ. 'BUOY' .AND. i .EQ. 4) THEN
        one_conv%obs4d(:, i) = 0
        CYCLE
      END IF
      nc_var = nc_dataset%getVariable(TRIM(st_in_name(i)))
      IF (one_conv%numSt .EQ. 1) THEN
        CALL nc_var%getData(data0d)
        one_conv%obs4d(:, i) = data0d
      ELSE
        CALL nc_var%getData(data1d)
        one_conv%obs4d(:, i) = data1d
      END IF
    END DO

    ! get variables one by one
    DO i = 1, this%numVars
      nc_var = nc_dataset%getVariable(TRIM(this%varNames(i)))
      IF (one_conv%numLevel == 1) THEN
        CALL nc_var%getData(data1d)
      ELSE
        CALL nc_var%getData(data2d)
        ! 2D dimension mismatch
        data_len = SHAPE(data2d)
        IF (data_len(1) /= one_conv%numSt) data2d = TRANSPOSE(data2d)
        one_conv%obsData(:, i, :) = data2d
      END IF
    END DO
            !!!!!! If get errors, should change the code!!!!!
    one_conv%obsErrs = 1
    CALL nc_dataset%CLOSE()
  END SUBROUTINE ingest_conv_each

  SUBROUTINE qc_synop(this)
    IMPLICIT NONE
    CLASS(ObsConvention_t2) :: this
    WRITE (*, *) "00: Finish ", TRIM(this%obsType), " quality control."
  END SUBROUTINE qc_synop

  SUBROUTINE qc_ship(this)
    IMPLICIT NONE
    CLASS(ObsConvention_t2) :: this
    WRITE (*, *) "01: Finish ", TRIM(this%obsType), " quality control."
  END SUBROUTINE qc_ship

  SUBROUTINE qc_buoy(this)
    IMPLICIT NONE
    CLASS(ObsConvention_t2) :: this
    WRITE (*, *) "02: Finish ", TRIM(this%obsType), " quality control."
  END SUBROUTINE qc_buoy

  SUBROUTINE qc_metar(this)
    IMPLICIT NONE
    CLASS(ObsConvention_t2) :: this
    WRITE (*, *) "03: Finish ", TRIM(this%obsType), " quality control."
  END SUBROUTINE qc_metar

  SUBROUTINE qc_temp(this, X)
    IMPLICIT NONE
    CLASS(ObsConvention_t2) :: this
    TYPE(State_t) :: X

    ! Local variables
    INTEGER(i_kind), PARAMETER :: nPolynomial = 6, nx = nPolynomial + 1
    INTEGER(i_kind) :: n_obs, rank, n_use, i, j, k, idx, numFit, nfill, lowHeight, uppHeight, lowBound
    INTEGER(i_kind), ALLOCATABLE :: t_obs(:), mask(:)
    INTEGER(i_kind) :: st_i, p_id, var_i, other_i

    REAL(r_kind), PARAMETER :: dryHydro = 50000.0D0
    REAL(r_kind), ALLOCATABLE :: r_obs(:, :), h_obs(:), llobs(:, :), error(:, :), lands(:), qvapor(:)
    REAL(r_kind) :: coef(nx), fit_range(2), xxx, temp(5)
    CHARACTER(LEN=20) :: site

    ! QR decomposition:
    INTEGER :: job, info, jpvt(nx)
    DOUBLE PRECISION, ALLOCATABLE :: cmatrix(:, :), fitting(:, :), qy(:), qty(:), xb(:), rsd(:)
    DOUBLE PRECISION :: xmin, xmax
    DOUBLE PRECISION :: aux(nx), work(nx), sol(nx)

    rank = X%mpddGlob%myrank

    DO st_i = 1, this%numSt
      n_obs = INT(this%obs4d(st_i, 5))
      ! Allocate data
      ALLOCATE (t_obs(n_obs), mask(n_obs), r_obs(n_obs, 5), h_obs(n_obs), llobs(1:2, n_obs), &
                error(n_obs, 1:5), lands(n_obs), fitting(n_obs, 2), qy(n_obs), qty(n_obs), &
                xb(n_obs), rsd(n_obs), qvapor(n_obs))
      site = TRIM(this%nameSt(st_i))
      ! Assign data
      t_obs = INT(this%obs4d(st_i, 1))
      ! r_obs = this%obsData(st_i, :, 1:n_obs)
      h_obs = this%obsHgt(st_i, 1:n_obs)
      llobs(1, :) = this%obs4d(st_i, 2)
      llobs(2, :) = this%obs4d(st_i, 3)
      ! error = this%obsErrs(st_i, :, 1:n_obs)
      lands = 1.0D0
      other_i = 4
      DO var_i = 1, this%numVars
        SELECT CASE (TRIM(this%varNames(var_i)))
        CASE ('t')
          r_obs(:, 1) = this%obsData(st_i, var_i, 1:n_obs)
          error(:, 1) = this%obsErrs(st_i, var_i, 1:n_obs)
        CASE ('dp')
          r_obs(:, 2) = this%obsData(st_i, var_i, 1:n_obs)
          error(:, 2) = this%obsErrs(st_i, var_i, 1:n_obs)
        CASE ('p')
          p_id = var_i
          r_obs(:, 3) = this%obsData(st_i, var_i, 1:n_obs)
          error(:, 3) = this%obsErrs(st_i, var_i, 1:n_obs)
        CASE DEFAULT
          r_obs(:, other_i) = this%obsData(st_i, var_i, 1:n_obs)
          error(:, other_i) = this%obsErrs(st_i, var_i, 1:n_obs)
          other_i = other_i + 1
        END SELECT
      END DO

      ! Merge obs with identical pressure levels:
      DO j = 1, n_obs
        DO i = j + 1, n_obs
        IF (ABS(r_obs(j, 3) - r_obs(i, 3)) .LT. 1.0D-8 .AND. r_obs(j, 3) .LT. missing) THEN
          ! For identical pressure levels:
          DO k = 1, 5  ! Merge the obs fields
          IF (k .NE. 3) THEN ! Pressure values do not change
            IF (r_obs(i, k) .LT. missing .AND. r_obs(j, k) .LT. missing) THEN
            IF (h_obs(j) .LT. missing) THEN ! Trust the obs with height information
              r_obs(i, k) = r_obs(j, k)
            ELSE
              r_obs(i, k) = 0.5D0 * (r_obs(i, k) + r_obs(j, k))
            END IF
            ELSE IF (r_obs(j, k) .LT. missing) THEN
            r_obs(i, k) = r_obs(j, k)
            END IF
          END IF
          IF (error(i, k) .LT. missing .AND. error(j, k) .LT. missing) THEN
            error(i, k) = 0.5D0 * (error(i, k) + error(j, k))
          ELSE IF (error(j, k) .LT. missing) THEN
            error(i, k) = error(j, k)
          END IF
          END DO
          h_obs(i) = MIN(h_obs(j), h_obs(i))
          t_obs(i) = MIN(t_obs(j), t_obs(i))
          llobs(1, i) = MIN(llobs(1, j), llobs(1, i))
          llobs(2, i) = MIN(llobs(2, j), llobs(2, i))
          lands(i) = MIN(lands(j), lands(i))

          r_obs(j, 1:5) = missing ! void the duplicated obs
          h_obs(j) = missing ! Void the height obs as well
        END IF
        END DO
      END DO

      ! Remove unused levels:
      n_use = 0
      DO j = 1, n_obs
        ! Valid pressure:
        IF (r_obs(j, 3) .LT. missing) THEN
          n_use = n_use + 1
          r_obs(n_use, 1:5) = r_obs(j, 1:5)
          h_obs(n_use) = h_obs(j)
          t_obs(n_use) = t_obs(j)
          llobs(1:2, n_use) = llobs(1:2, j)
          error(n_use, 1:5) = error(j, 1:5)
          lands(n_use) = lands(j)
        ELSE
#ifdef DEBUG_SOUNDING
          IF (h_obs(j) .LT. missing) &
            WRITE (*, 11) j, n_obs, r_obs(j, 3), h_obs(j), rank, TRIM(site)
11        FORMAT('obs height - missing pressure: ', I3, I4, ' HghtPres: ', 2D12.4, ' pc:', I2, ' stn: ', A)
#endif
        END IF
      END DO

      ! Possible sorting according to the pressure levels: from ground to the top
      DO j = 1, n_use
        DO i = 1, n_use - j
        IF (r_obs(i, 3) .LT. r_obs(i + 1, 3)) THEN
          temp(1:5) = r_obs(i, 1:5)
          r_obs(i, 1:5) = r_obs(i + 1, 1:5)
          r_obs(i + 1, 1:5) = temp(1:5)

          temp(1) = h_obs(i)
          h_obs(i) = h_obs(i + 1)
          h_obs(i + 1) = temp(1)

          temp(1) = t_obs(i)
          t_obs(i) = t_obs(i + 1)
          t_obs(i + 1) = temp(1)

          temp(1:2) = llobs(1:2, i)
          llobs(1:2, i) = llobs(1:2, i + 1)
          llobs(1:2, i + 1) = temp(1:2)

          temp(1:5) = error(i, 1:5)
          error(i, 1:5) = error(i + 1, 1:5)
          error(i + 1, 1:5) = temp(1:5)

          temp(1) = lands(i)
          lands(i) = lands(i + 1)
          lands(i + 1) = temp(1)
        END IF
        END DO
      END DO

      ! Fill in temperature and dew point:
      DO idx = 1, 2
        numFit = 0
        nfill = 0
        fitting = 0.0D0
        fit_range(1) = 1.0D10 ! Range of the fitting in vertical
        fit_range(2) = 0.0D0
        DO j = 1, n_use
        IF (r_obs(j, idx) .LT. missing) THEN
          ! Count fitting points:
          numFit = numFit + 1
          fitting(numFit, 1) = LOG(r_obs(j, 3))
          fitting(numFit, 2) = r_obs(j, idx)
          IF (r_obs(j, 3) .GT. fit_range(2)) fit_range(2) = r_obs(j, 3)
          IF (r_obs(j, 3) .LT. fit_range(1)) fit_range(1) = r_obs(j, 3)
        ELSE
          nfill = nfill + 1
        END IF
        END DO
        ! Cannot fit with too little obs:
        IF (numFit .LT. nx) THEN
          WRITE (*, 1) TRIM(site), idx, numFit, nPolynomial, n_use, n_obs, r_obs(1, idx), rank
1         FORMAT('FillingMissing: ', A, ' var: ', I1, ' - too little data: ', 2I2, &
                 ' in use: ', I3, ' total: ', I3, ' values: ', D10.2, ' pc', I2)
          CYCLE
        END IF

        ! Fitting the profile:
        ALLOCATE (cmatrix(numFit, nx))
        xmin = LOG(fit_range(1)) ! MINVAL(fitting(1:numFit,1))
        xmax = LOG(fit_range(2)) ! MAXVAL(fitting(1:numFit,1))
        DO j = 1, numFit
        DO i = 0, nPolynomial
          cmatrix(j, i + 1) = ((fitting(j, 1) - xmin) / (xmax - xmin))**i
        END DO
        END DO
        ! Solve the coefficients:
        job = 1 ! Pivoting
        jpvt = 0
        CALL dqrdc(cmatrix, numFit, numFit, nx, aux, jpvt, work, job)
        job = 110 ! Just compute the xbb only
        rsd = 0.0D0
        CALL dqrsl(cmatrix, numFit, numFit, nx, aux, fitting(:, 2), qy, qty, sol, rsd, xb, job, info)

        ! Pivot back:
        DO i = 1, nx
          coef(jpvt(i)) = sol(i)
        END DO

#ifdef DEBUG_SOUNDING
        WRITE (*, 2) nPolynomial, info, sol(1:nx)
        WRITE (*, 3) nPolynomial, info, coef(1:nx)
2       FORMAT('Least square solution: ', I2, ' info: ', I3, ' sol: ', 10D10.2)
3       FORMAT('Least square pivoted : ', I2, ' info: ', I3, ' sol: ', 10D10.2)
        WRITE (*, 4) idx, MAXVAL(rsd), MAXVAL(fitting(:, 2)), x%mpddGlob%myrank
4       FORMAT('Max residual of ', I1, ' is ', D12.4, ' Max values: ', D12.4, ' pc: ', I2)

        IF (x%mpddGlob%myrank .EQ. 2) THEN
        IF (idx .EQ. 2) THEN
          DO i = 1, numFit
            xb(i) = 0.0D0
            DO j = 1, nx
              xb(i) = xb(i) + coef(j) * ((fitting(i, 1) - xmin) / (xmax - xmin))**(j - 1)
            END DO
            !WRITE(11,5) i,(fitting(i,1)-xmin)/(xmax-xmin), &
            !  fitting(i,2)/MAXVAL(fitting(:,2)),xb(i)/MAXVAL(fitting(:,2))
            WRITE (11, 5) i, fitting(i, 1), LOG(fitting(i, 1)), &
              xb(i), fitting(i, 2)
5           FORMAT(I4, 4F12.4)
          END DO

          ! Plotting the whole used obs:

          nfill = 0
          DO i = 1, n_use
          IF (r_obs(i, idx) .GT. missing - 1.0D0) THEN
            nfill = nfill + 1
            PRINT *, 'Missing temDew: ', i, n_use, idx, r_obs(i, 3), r_obs(i, idx)
            r_obs(i, idx) = 190.0D0
          END IF
          END DO
          PRINT *, 'Total of filling: ', nfill
          !xmin = MINVAL(LOG(r_obs(1:n_use,3)))
          !xmax = MAXVAL(LOG(r_obs(1:n_use,3)))
          !print*,'MMX2: ',xmin,xmax,rank,idx,MAXVAL(r_obs(1:n_use,3)),MINVAL(r_obs(1:n_use,3))
          DO i = 1, n_use
            xb(i) = 0.0D0
            DO j = 1, nx
              xxx = LOG(r_obs(i, 3))
              xb(i) = xb(i) + coef(j) * ((xxx - xmin) / (xmax - xmin))**(j - 1)
            END DO

            ! No extrapolation of the polynomial:
            IF (r_obs(i, 3) .LT. fit_range(1) .OR. r_obs(i, 3) .GT. fit_range(2)) CYCLE

            IF (r_obs(i, idx) .GT. missing - 1.0D0) THEN
              WRITE (12, 6) i, r_obs(i, 3), LOG(r_obs(i, 3)), xb(i)
            ELSE
              WRITE (12, 6) i, r_obs(i, 3), LOG(r_obs(i, 3)), xb(i), r_obs(i, idx)
            END IF
6           FORMAT(I4, 4F12.4)
          END DO
        END IF
        END IF
#endif

        ! Filling in temperature and dewpoint:
        DO i = 1, n_use
          ! No extrapolation of the polynomial:
          IF (r_obs(i, 3) .LT. fit_range(1) .OR. r_obs(i, 3) .GT. fit_range(2)) CYCLE

          ! No replacement for dewpoint if the obs is valid:
          IF (idx .EQ. 2 .AND. r_obs(i, idx) .LT. missing) CYCLE

          xb(i) = 0.0D0
          DO j = 1, nx
            xxx = LOG(r_obs(i, 3))
            xb(i) = xb(i) + coef(j) * ((xxx - xmin) / (xmax - xmin))**(j - 1)
          END DO
          ! Replace all temperature or dewpoint obs with the polynomial fitting:
          r_obs(i, idx) = xb(i)
        END DO

        DEALLOCATE (cmatrix)
      END DO

      ! Heights fill in:

      ! Masking the vertical levels with no fill: 0; need fill: 1; and with height: 2
      mask = 0
      qvapor = 0.0D0
      DO i = 1, n_use
        IF (r_obs(i, 3) .LT. dryHydro) THEN
          ! No qvapor is needed for hydrostatic height:
          IF (r_obs(i, 1) .LT. missing) THEN
            qvapor(i) = 0.0D0
            IF (h_obs(i) .GT. missing - 1.0D0) THEN
              mask(i) = 1
            ELSE
              mask(i) = 2
            END IF
          END IF
        ELSE
          ! Qvapor is needed for a hydrostatic height:
          IF (r_obs(i, 1) .LT. missing .AND. r_obs(i, 2) .LT. missing) THEN
            qvapor(i) = Td_to_qvapor(r_obs(i, 2), r_obs(i, 3))
            IF (h_obs(i) .GT. missing - 1.0D0) THEN
              mask(i) = 1
            ELSE
              mask(i) = 2
            END IF
          END IF
        END IF
      END DO

      ! Filling in heights: assuming the pressure levels have been sorted above
      lowHeight = 0
      uppHeight = -1
      lowBound = 0
      nfill = 0
      DO i = 1, n_use
        ! Find the consecutive fillin interval:
        IF (mask(i) .EQ. 1) THEN
        IF (uppHeight .LT. lowHeight) THEN
          lowHeight = i; uppHeight = i  ! First encounter of fillin section
        END IF
        ! Hitting the top:
        IF (i .EQ. n_use) THEN
          IF (lowBound .GT. 0) THEN
            CALL hydroHeightUp(i - lowBound + 1, r_obs(lowBound:i, 3), &
                               r_obs(lowBound:i, 1), qvapor(lowBound:i), h_obs(lowBound:i))
            nfill = nfill + n_use - lowBound
          END IF
        ELSE
          uppHeight = i
        END IF
        ELSE IF (mask(i) .EQ. 0) THEN ! No fill level
        IF (lowBound .GT. 0) THEN
          CALL hydroHeightUp(i - lowBound, r_obs(lowBound:i - 1, 3), &
                             r_obs(lowBound:i - 1, 1), qvapor(lowBound:i - 1), h_obs(lowBound:i - 1)) ! no fill in current level
          nfill = nfill + uppHeight - lowHeight + 1
        END IF
        ! Reset:
        lowBound = 0
        lowHeight = 0
        uppHeight = -1
        ELSE ! mask = 2
        ! Filling the section:
        IF (uppHeight .GE. lowHeight) THEN
          CALL hydroHeightDown(i - lowHeight + 1, r_obs(lowHeight:i, 3), &
                               r_obs(lowHeight:i, 1), qvapor(lowHeight:i), h_obs(lowHeight:i))
          nfill = nfill + uppHeight - lowHeight + 1
        END IF
        lowHeight = i + 1; uppHeight = i; lowBound = i
        END IF
      END DO
#ifdef DEBUG
      WRITE (*, 7) nfill, x%mpddGlob%myrank, x%sg%gLevel
7     FORMAT('fillinSoundHeights - total fillin: ', I4, ' pc: ', I2, ' Glvl: ', I2)
#endif
      ! Assign modified values to obsData at p_id
      this%obsData(st_i, p_id, :) = missing
      this%obsData(st_i, p_id, 1:n_use) = r_obs(1:n_use, 3)
      this%obs4d(st_i, 5) = n_use
      DEALLOCATE (t_obs, mask, r_obs, h_obs, llobs, error, &
                  lands, fitting, qy, qty, xb, rsd, qvapor)
    END DO

    WRITE (*, *) "04: Finish ", TRIM(this%obsType), " quality control."
  END SUBROUTINE qc_temp

  SUBROUTINE qc_profl(this)
    IMPLICIT NONE
    CLASS(ObsConvention_t2) :: this
    WRITE (*, *) "05: Finish ", TRIM(this%obsType), " quality control."
  END SUBROUTINE qc_profl
END MODULE ObsConvention_m2
