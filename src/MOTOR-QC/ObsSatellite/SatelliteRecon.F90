!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.SatelliteReconstruction
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Ting Shu
! FUNCTION          : Wavelet-based satellite reconstruction
! VERSION           : V 0.0
! HISTORY           :
!   Created by Ting Shu (shuting@gbamwf.com), 2023/08/09, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements Wavelet-based satellite reconstruction.
MODULE satelliteRecon_m
    USE kinds_m, ONLY: r_kind, i_kind
    USE wavelet_m, ONLY: wavelet_t
    USE satelliteW_m, ONLY: get_latlon, grid_idx, decompress_bi_array
    USE slint, ONLY: slint_init, bilinear_interp, tgt_grid
    USE parameters_m, ONLY: degree2radian, Missing
    USE netcdf
    USE YAMLRead_m
    USE State_m, ONLY: State_t
    USE ObsSet_m, ONLY: ObsSet_t
    USe MPObs_m, ONLY: MPObs_t
    USE ObsField_m, ONLY: obsField_t
    IMPLICIT NONE

    TYPE :: satCoeff
        CHARACTER(LEN=:), ALLOCATABLE :: cfileName, ifileName
        TYPE(wavelet_t), ALLOCATABLE :: dataCoeffs(:)
        REAL(r_kind), ALLOCATABLE :: fineInfo(:, :, :), dataRange(:, :)
        INTEGER(i_kind), ALLOCATABLE :: dataDecomFlag(:), missValPos(:, :, :)
        INTEGER(i_kind) :: numDeData  ! Number of decomposed data
    END TYPE satCoeff
    
    TYPE :: satelliteRecon_t
        CHARACTER(LEN=1024) :: configFile
        CHARACTER(LEN=10) :: obsType
        CHARACTER(LEN=20), ALLOCATABLE :: varNames(:), chanList(:)
        INTEGER(i_kind) :: mgStart, mgEnd, numFiles, numVars, numInfo
        INTEGER(i_kind), DIMENSION(:), ALLOCATABLE :: latlon_times
        REAL(r_kind), DIMENSION(:), ALLOCATABLE :: fineLat, fineLon, time
        ! Size of allCoeffs equal to numFiles
        TYPE(satCoeff), ALLOCATABLE :: allCoeffs(:)
        ! Flag of processing error via wavelet
        LOGICAL :: processErr

        CONTAINS
        PROCEDURE, PUBLIC :: ingest => sat_recon_ingest
        PROCEDURE, PUBLIC :: recon => sat_recon_data
        PROCEDURE, PUBLIC :: destroy => sat_deallocate
    END TYPE satelliteRecon_t
    PRIVATE :: sat_recon_ingest, sat_recon_data, load_nc, interpolate_data
    
    CONTAINS
    SUBROUTINE sat_recon_ingest(this, configFile, tSlots, sat_platform_name, sat_inst_name)
        IMPLICIT NONE
        CLASS(satelliteRecon_t), INTENT(INOUT) :: this
        INTEGER(i_kind), INTENT(IN) :: tSlots
        CHARACTER(LEN=*), INTENT(IN) :: configFile, sat_platform_name, sat_inst_name
        ! Local variables
        CHARACTER(LEN=1024) :: obsFileDir, temp_fname
        CHARACTER(LEN=1024), ALLOCATABLE :: obsFileNames(:)
        REAL(r_kind), DIMENSION(:), ALLOCATABLE :: latlon_domain
        INTEGER(i_kind) :: status, i, f_len, latlon_len(2), temp_i
        REAL(r_kind) :: temp_d, dxdy(2), new_latlon_domain(4)
        CHARACTER(LEN=1) :: tStr
        LOGICAL :: fileExist
        REAL(r_kind) :: t1, t2

        CALL CPU_TIME(t1)
        ! Assign configFile
        this%configFile = configFile
        ! Get wavelet coefficient nc file names
        status = yaml_get_var(configFile, "IO", 'output_dir', obsFileDir)
        status = yaml_get_var(configFile, "SatelliteInfo", 'varNames', this%varNames)
        status = yaml_get_var(configFile, "SatelliteInfo", 'obsType', this%obsType)
        ! Get channel number list
        status = yaml_get_var(configFile, "SatelliteInfo", 'chanList', this%chanList)

        ALLOCATE(obsFileNames(tSlots), STAT=status)
        temp_i = 0
        DO i = 1, tSlots
            WRITE(tStr, "(I1)") i
            temp_fname = TRIM(obsFileDir)//'/'//TRIM(sat_platform_name)//&
                         '_'//TRIM(sat_inst_name)//'_after_QC_BC_'//tStr//"_wcoeff.nc"
            INQUIRE(file=TRIM(temp_fname), exist=fileExist)
            IF (fileExist) THEN
                temp_i = temp_i + 1
                obsFileNames(temp_i) = TRIM(temp_fname)
            END IF
        END DO
        this%numFiles = temp_i
        ALLOCATE(this%allCoeffs(this%numFiles), this%time(this%numFiles), STAT=status)
        DO i = 1, this%numFiles
            this%allCoeffs(i)%cfileName = TRIM(obsFileNames(i))
            f_len = LEN(TRIM(this%allCoeffs(i)%cfileName)) - 10
            this%allCoeffs(i)%ifileName = this%allCoeffs(i)%cfileName(1 : f_len)//".nc"
        END DO

        ! Get start and end grid id
        status = yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', this%mgStart)
        status = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', this%mgEnd)
        ! Get the flag wheather needs to process error via wavelet
        status = yaml_get_var(TRIM(configFile), 'Wavelet', 'processErr', this%processErr)
        ! Get latitude and longitude domain and times
        status = yaml_get_var(TRIM(configFile), 'latlon_grid', 'domain_latlon', latlon_domain)
        status = yaml_get_var(TRIM(configFile), 'Wavelet', 'latlon_times', this%latlon_times)
        CALL CPU_TIME(t2)
        PRINT *, "************Ingest: running time of reading yaml file:", t2-t1
        
        ! Compute the latitude and longitude of the finest grid
        CALL get_latlon(this%mgStart, this%mgEnd, this%latlon_times, latlon_domain, this%fineLat, this%fineLon)
        CALL CPU_TIME(t1)
        PRINT *, "************Ingest: running time of calculating finest latitude and longitude:", t1-t2

        ! Load nc files to assign data allCoeffs
        CALL load_nc(this)
        CALL CPU_TIME(t2)
        PRINT *, "************Ingest: running time of loading and converting wavelet coefficients:", t2-t1

        ! Deallocate variables
        IF (ALLOCATED(obsFileNames)) DEALLOCATE(obsFileNames, STAT=status)
        IF (ALLOCATED(latlon_domain)) DEALLOCATE(latlon_domain, STAT=status)
        CALL CPU_TIME(t1)
        PRINT *, "************Ingest: deallocating variables:", t1-t2
    END SUBROUTINE sat_recon_ingest

    SUBROUTINE sat_recon_data(this, state, thinObs, mpObs)
        IMPLICIT NONE
        CLASS(satelliteRecon_t), INTENT(INOUT) :: this
        TYPE(State_t), INTENT(IN) :: state
        TYPE(ObsSet_t), INTENT(INOUT) :: thinObs
        TYPE(MPObs_t), INTENT(IN) :: mpObs

        ! Local variables
        INTEGER(i_kind) :: i, j, k, l, numC, dim1, dim2, reLevel, status, temp_idx, temp_i, mgNo
        TYPE(wavelet_t) :: coeff
        REAL(r_kind), ALLOCATABLE :: temp_x1(:, :), temp_data(:, :), temp_error(:, :), temp_info(:, :, :), obs_lat(:), obs_lon(:)
        INTEGER(i_kind), ALLOCATABLE :: latIdx(:), lonIdx(:), state_timeIdx(:), obs_timeIdx(:), miss_array(:, :)
        REAL(r_kind) :: min, max
        INTEGER(i_kind) :: all_stateIdx(state%sg%tSlots), all_obsIdx(state%sg%tSlots), num_slots, num_obs, &
                           len_lat, len_lon, temp_count, temp_hIdx, numDecomData, numOneWave
        REAL(r_kind) :: t1, t2
        
        CALL CPU_TIME(t1)
        mgNo = state%sg%gLevel
        ! Observation file avaliable index
        all_obsIdx = 0
        DO i = 1, state%sg%tSlots
            CALL get_coeff_id(this, state%sg%tt(i), temp_idx)
            all_stateIdx(i) = i
            all_obsIdx(i) = temp_idx
        END DO

        num_slots = COUNT(all_obsIdx .NE. 0)
        IF (num_slots .EQ. 0) THEN
            PRINT *, "All the loading satellite data are outside the analysis time period."
            RETURN
        END IF
        CALL CPU_TIME(t2)
        PRINT *, " ************Recon: running time of checking avariable time slots:", t2-t1

        PRINT *, "==== Start wavelet-based stallite thinning. ===="
        ALLOCATE(state_timeIdx(num_slots), obs_timeIdx(num_slots), STAT=status)
        ! State avaliable time index
        state_timeIdx = PACK(all_stateIdx, all_obsIdx .NE. 0)
        obs_timeIdx = all_obsIdx(state_timeIdx)
        ! Reconstruct satellite data according to state
        reLevel = mgNo - this%mgStart
        CALL grid_idx(this%mgStart, this%mgEnd, this%latlon_times, mgNo, latIdx, lonIdx)
        len_lat = SIZE(latIdx)
        len_lon = SIZE(lonIdx)
        ! Assign observation latitude and longitude in radian
        ALLOCATE(obs_lat(len_lat), obs_lon(len_lon), STAT=status)
        obs_lat = this%fineLat(latIdx) * degree2radian
        obs_lon = this%fineLon(lonIdx) * degree2radian
        ! ALLOCATE(thinObs%obsFields(this%numVars), STAT=status)
        CALL CPU_TIME(t1)
        PRINT *, " ************Recon: running time of getting corresponding index for the current grid:", t1-t2
        
        temp_i = 0
        DO i = 1, this%numVars
            PRINT *, "Start wavelet-based stallite thinning for varibale: ", this%varNames(i)
            ! Calculate the number of observation for the current variable
            num_obs = 0
            DO j = 1, num_slots
                IF (this%allCoeffs(obs_timeIdx(j))%dataDecomFlag(i) .NE. 0) &
                    num_obs = num_obs + COUNT(this%allCoeffs(j)%missValPos(i, latIdx, lonIdx) .EQ. 0)
            END DO
            thinObs%obsFields(i) = obsField_t(this%configFile, mpObs)
            ! Allocate obsSet data structure
            ALLOCATE(thinObs%obsFields(i)%idx(num_obs), thinObs%obsFields(i)%values(num_obs), &
                     thinObs%obsFields(i)%errors(num_obs), STAT=status)
            ! Allocate satellite related angles
            ALLOCATE(thinObs%obsFields(i)%ObsAttrSat%azangle(num_obs), thinObs%obsFields(i)%ObsAttrSat%zenangle(num_obs), &
                     thinObs%obsFields(i)%ObsAttrSat%sunazangle(num_obs), thinObs%obsFields(i)%ObsAttrSat%sunzenangle(num_obs), STAT=status)
            ! CALL thinObs%obsFields(i)%Set_Name(TRIM(this%varNames(i)))
            CALL thinObs%obsFields(i)%Set_Name(TRIM(this%varNames(i)(6:8)))
            CALL thinObs%obsFields(i)%Set_ObsType(TRIM(this%obsType)//'_'//TRIM(this%varNames(i)(1:4)))
            ! Set variable type, 2023-08-18, just set SP for satellite
            ! CALL thinObs%obsFields(i)%Set_ValType("SP")
            IF (num_obs .NE. 0) THEN
                temp_i = temp_i + 1
                temp_count = 0
                DO j = 1, num_slots
                    temp_idx = obs_timeIdx(j)
                    ALLOCATE(temp_data(len_lat, len_lon), temp_error(len_lat, len_lon), &
                    temp_info(this%numInfo, len_lat, len_lon), STAT=status)
                    temp_info = this%allCoeffs(temp_idx)%fineInfo(:, latIdx, lonIdx)
                    ! Missing value position information
                    ALLOCATE(miss_array(len_lat, len_lon), STAT=status)
                    miss_array = this%allCoeffs(temp_idx)%missValPos(i, latIdx, lonIdx)
                    ! Get number of decompressed data
                    numDecomData = this%allCoeffs(temp_idx)%numDeData
                    IF (this%processErr) THEN
                        numOneWave = 2
                    ELSE  ! Error is in temp_info
                        numOneWave = 1
                        temp_error = temp_info(4+temp_i, :, :)
                    END IF
                    DO k = 1, numOneWave
                        coeff = this%allCoeffs(temp_idx)%dataCoeffs((k-1)*numDecomData+temp_i)
                        CALL CPU_TIME(t1)
                        CALL coeff%recon(reLevel, temp_x1)
                        CALL CPU_TIME(t2)
                        PRINT *, "************Recon: running time of reconstructing data:", t2-t1, i, j, k
                        min = this%allCoeffs(temp_idx)%dataRange((k-1)*this%numVars+i, 1)
                        max = this%allCoeffs(temp_idx)%dataRange((k-1)*this%numVars+i, 2)
                        WHERE(temp_x1 .LT. min) temp_x1=min
                        WHERE(temp_x1 .GT. max) temp_x1=max
                        IF (k .EQ. 1) THEN
                            temp_data = temp_x1(latIdx, lonIdx)
                        ELSE
                            temp_error = temp_x1(latIdx, lonIdx)
                        END IF
                        IF (ALLOCATED(temp_x1)) DEALLOCATE(temp_x1, STAT=status)
                    END DO
                    temp_hIdx = 0
                    CALL CPU_TIME(t1)
                    DO k = 1, len_lat
                        DO l = 1, len_lon
                            temp_hIdx = temp_hIdx + 1
                            ! Do not save missing value
                            IF (miss_array(k, l) .EQ. 0) THEN
                                ! Check wheather at the same grid or not, and the latitude and longitude values just keep five decimal points
                                IF (ABS(state%sg%cell_cntr(1, temp_hIdx)-obs_lat(k)) .LT. 1E-5 .AND. &
                                ABS(state%sg%cell_cntr(2, temp_hIdx)-obs_lon(l)) .LT. 1E-5) THEN
                                    temp_count = temp_count + 1
                                    ! Assign state time, horizontal, and vectical index
                                    thinObs%obsFields(i)%idx(temp_count)%tIdx = state_timeIdx(j)
                                    thinObs%obsFields(i)%idx(temp_count)%hIdx = temp_hIdx
                                    ! Satellite data just has one level along the vertical direction
                                    thinObs%obsFields(i)%idx(temp_count)%vIdx = 1

                                    ! Assign values
                                    thinObs%obsFields(i)%values(temp_count) = temp_data(k, l)

                                    ! Assign errors
                                    thinObs%obsFields(i)%errors(temp_count) = temp_error(k, l)

                                    ! Assign four angles
                                    thinObs%obsFields(i)%ObsAttrSat%azangle(temp_count) = temp_info(1, k, l) * degree2radian
                                    thinObs%obsFields(i)%ObsAttrSat%zenangle(temp_count) = temp_info(2, k, l) * degree2radian
                                    thinObs%obsFields(i)%ObsAttrSat%sunazangle(temp_count) = temp_info(3, k, l) * degree2radian
                                    thinObs%obsFields(i)%ObsAttrSat%sunzenangle(temp_count) = temp_info(4, k, l) * degree2radian
                                END IF
                            END IF
                        END DO
                    END DO
                    CALL CPU_TIME(t2)
                    PRINT *, "************Recon: running time of putting data into thinObs:", t2-t1, i, j
                    IF (ALLOCATED(temp_data)) DEALLOCATE(temp_data, STAT=status)
                    IF (ALLOCATED(temp_error)) DEALLOCATE(temp_error, STAT=status)
                    IF (ALLOCATED(temp_info)) DEALLOCATE(temp_info, STAT=status)
                    IF (ALLOCATED(miss_array)) DEALLOCATE(miss_array, STAT=status)
                    CALL coeff%release()
                    CALL CPU_TIME(t1)
                    PRINT *, "************Recon: running time of deallocating do-loop variables:", t1-t2, i, j
                END DO
            END IF
        END DO

        CALL CPU_TIME(t1)
        IF (ALLOCATED(obs_lat)) DEALLOCATE(obs_lat, STAT=status)
        IF (ALLOCATED(obs_lon)) DEALLOCATE(obs_lon, STAT=status)
        IF (ALLOCATED(latIdx)) DEALLOCATE(latIdx, STAT=status)
        IF (ALLOCATED(lonIdx)) DEALLOCATE(lonIdx, STAT=status)
        IF (ALLOCATED(state_timeIdx)) DEALLOCATE(state_timeIdx, STAT=status)
        IF (ALLOCATED(obs_timeIdx)) DEALLOCATE(obs_timeIdx, STAT=status)
        CALL CPU_TIME(t2)
        PRINT *, "************Recon: running time of deallocating subroutine variables:", t2-t1
        PRINT *, "Wavelet-based satellite thinning is done for grid ", state%sg%gLevel
    END SUBROUTINE sat_recon_data

    ! Load and rebuild coefficients from nc file
    ! Need call rebuild_coeff and interpolate_data subroutines
    SUBROUTINE load_nc(this)
        IMPLICIT NONE
        CLASS(satelliteRecon_t), INTENT(INOUT) :: this
        ! Local variables
        INTEGER(i_kind) :: i, j, status, nc_id, dim_id, var_id, len1, &
                           posNLen, wdLen, dLen, temp_i, temp_j
        ! Interpolation related variables
        REAL(r_kind), ALLOCATABLE :: temp_info(:, :), temp_var(:)
        REAL(r_kind) :: obsTime(this%numFiles)
        ! Information name list
        CHARACTER(LEN=10), PARAMETER :: info_name(6) = (/"lat       ", "lon       ",&
                                        "sat_azi   ", "sat_zenith", "sol_azi   ", "sol_zenith"/)
        ! Compressed missing position information related variables
        INTEGER(i_kind), ALLOCATABLE :: temp_com_miss(:), temp_miss(:)
        CHARACTER(LEN=2) :: idStr
        INTEGER(i_kind) :: com_len, fLatS, fLonS
        REAL(r_kind) :: t1, t2

        DO i = 1, this%numFiles
            ! Open satellite wavelet coefficient files one by one
            status = NF90_OPEN(TRIM(this%allCoeffs(i)%cfileName), NF90_NOWRITE, nc_id)
            ! Get data coefficient
            status = NF90_INQ_DIMID(nc_id, "posN", dim_id)
            status = NF90_INQUIRE_DIMENSION(nc_id, dim_id, LEN=posNLen)
            status = NF90_INQ_DIMID(nc_id, "dlen2times", dim_id)
            status = NF90_INQUIRE_DIMENSION(nc_id, dim_id, LEN=wdLen)
            ! number of variables
            IF (i .EQ. 1) THEN
                IF (this%processErr) THEN
                    dLen = wdLen/2
                ELSE
                    dLen = wdLen
                END IF
                this%numVars = dLen
            END IF
            ! Allocate data coefficient and dataDecomFlag
            ALLOCATE(this%allCoeffs(i)%dataDecomFlag(wdLen), STAT=status)
            status = NF90_INQ_VARID(nc_id, "dataDecomFlag", var_id)
            status = NF90_GET_VAR(nc_id, var_id, this%allCoeffs(i)%dataDecomFlag)
            IF (this%processErr) THEN  ! Include error
                this%allCoeffs(i)%numDeData = COUNT(this%allCoeffs(i)%dataDecomFlag .NE. 0)/2
                ALLOCATE(this%allCoeffs(i)%dataCoeffs(this%allCoeffs(i)%numDeData*2), STAT=status)
            ELSE  ! No error
                this%allCoeffs(i)%numDeData = COUNT(this%allCoeffs(i)%dataDecomFlag .NE. 0)
                ALLOCATE(this%allCoeffs(i)%dataCoeffs(this%allCoeffs(i)%numDeData), STAT=status)
            END IF
            ! Get data range (minimum and maximum of original satellite data)
            ALLOCATE(this%allCoeffs(i)%dataRange(wdLen, posNLen), STAT=status)
            status = NF90_INQ_VARID(nc_id, "dataRange", var_id)
            status = NF90_GET_VAR(nc_id, var_id, this%allCoeffs(i)%dataRange)
            ! Allocate missing position array
            fLatS = SIZE(this%fineLat)
            fLonS = SIZE(this%fineLon)
            ALLOCATE(this%allCoeffs(i)%missValPos(this%numVars, fLatS, fLonS), STAT=status)
            ! If the i_th variable is missing, then missValPos(i, :, :) = 1
            this%allCoeffs(i)%missValPos = 1
            temp_i = 0
            DO j = 1, wdLen
                IF (this%allCoeffs(i)%dataDecomFlag(j) .NE. 0) THEN
                    temp_i = temp_i + 1
                    CALL CPU_TIME(t1)
                    CALL rebuild_coeff(this, nc_id, posNLen, i, j, temp_i)
                    CALL CPU_TIME(t2)
                    PRINT *, "************LoadingNC: running time of rebuilding coefficient:", t2-t1, i, j
                    ! Get the missing position compressed array
                    IF (j .LE. dLen) THEN
                        ! Convert data id from number to string
                        IF (j < 10) THEN
                            WRITE(idStr, "(A1I1)") "0", j
                        ELSE
                            WRITE(idStr, "(I2)") j
                        END IF
                        status = NF90_INQ_DIMID(nc_id, "comMissN"//TRIM(idStr), dim_id)
                        status = NF90_INQUIRE_DIMENSION(nc_id, dim_id, LEN=com_len)
                        status = NF90_INQ_VARID(nc_id, "comMissArray"//TRIM(idStr), var_id)
                        ALLOCATE(temp_com_miss(com_len), STAT=status)
                        status = NF90_GET_VAR(nc_id, var_id, temp_com_miss)
                        ! Decompress missing position array
                        CALL CPU_TIME(t1)
                        CALL decompress_bi_array(temp_com_miss, temp_miss)
                        CALL CPU_TIME(t2)
                        PRINT *, "************LoadingNC: running time of rebuilding mask:", t2-t1, i, j
                        IF (SIZE(temp_miss) .NE. fLatS*fLonS) THEN
                            PRINT *, "Error: we did not get the position array of missing values."
                            PRINT *, "Eomething wrong with "//"comMissArray"//TRIM(idStr), " in ",  TRIM(this%allCoeffs(i)%cfileName)
                            PRINT *, "But if you donot use qc for satellite data, please ignore this message."
                            this%allCoeffs(i)%missValPos(j, :, :) = 0
                        ELSE
                            ! Only variable has no missing values, missValPos is assign
                            this%allCoeffs(i)%missValPos(j, :, :) = RESHAPE(temp_miss, (/fLatS, fLonS/))
                        END IF
                        IF (ALLOCATED(temp_com_miss)) DEALLOCATE(temp_com_miss, STAT=status)
                        IF (ALLOCATED(temp_miss)) DEALLOCATE(temp_miss, STAT=status)
                    END IF
                END IF
            END DO
            ! Load satellite information related variables, lat, lon, four angels, errors
            status = NF90_GET_ATT(nc_id, NF90_GLOBAL, "observationTime", this%time(i))
            status = NF90_CLOSE(nc_id)

            ! Load information nc file
            status = NF90_OPEN(TRIM(this%allCoeffs(i)%ifileName), NF90_NOWRITE, nc_id)
            ! get dimension length
            status = NF90_INQ_DIMID(nc_id, "NOBS", dim_id)
            status = NF90_INQUIRE_DIMENSION(nc_id, dim_id, LEN=len1)
            ! get values for satellite information data
            IF (this%processErr) THEN
                ! latitude, longitude, satellite azi, satellite zenith, sun azi, sun zenith
                ALLOCATE(temp_info(len1, SIZE(info_name)), temp_var(len1), STAT=status)
            ELSE
                ! latitude, longitude, satellite azi, satellite zenith, sun azi, sun zenith, error
                ALLOCATE(temp_info(len1, SIZE(info_name)+this%allCoeffs(i)%numDeData), temp_var(len1), STAT=status)
                ! Load error data from nc file, ERROR
                temp_j = 0
                DO j = 1, SIZE(this%chanList)
                    IF (this%allCoeffs(i)%dataDecomFlag(j) .NE. 0) THEN
                        temp_j = temp_j + 1
                        status = NF90_INQ_VARID(nc_id, "ERROR"//TRIM(this%chanList(j)), var_id)
                        status = NF90_GET_VAR(nc_id, var_id, temp_var)
                        temp_info(:, SIZE(info_name)+temp_j) = temp_var
                    END IF
                END DO
            END IF
            ! latitude, longitude, satellite azi, satellite zenith, sun azi, sun zenith
            DO j = 1, SIZE(info_name)
                status = NF90_INQ_VARID(nc_id, TRIM(info_name(j)), var_id)
                status = NF90_GET_VAR(nc_id, var_id, temp_var)
                temp_info(:, j) = temp_var
            END DO
            ! Close nc file
            status = NF90_CLOSE(nc_id)
            CALL CPU_TIME(t1)
            CALL interpolate_data(this, temp_info, this%allCoeffs(i)%fineInfo)
            CALL CPU_TIME(t2)
            PRINT *, "************LoadingNC: running time of interpolating latlon and angles:", t2-t1

            ! Deallocate temp data
            IF (ALLOCATED(temp_info)) DEALLOCATE(temp_info, STAT=status)
            IF (ALLOCATED(temp_var)) DEALLOCATE(temp_var, STAT=status)
        END DO
    END SUBROUTINE load_nc

    SUBROUTINE interpolate_data(this, infoData, interpData)
        IMPLICIT NONE
        CLASS(satelliteRecon_t), INTENT(INOUT) :: this
        REAL(r_kind), ALLOCATABLE, INTENT(IN) :: infoData(:, :)
        REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: interpData(:, :, :)
        ! Interpolation related variables
        REAL(r_kind), ALLOCATABLE :: tLatLon(:, :), sLatlon(:, :), tempData(:)
        INTEGER(i_kind) :: tLatS, tLonS, tSize, sSize, i, j, k, status, slintIdx(3)
        REAL(r_kind) :: slintCoeff(3), slintValue, t1, t2

        ! interpolation for information related data
        ! Get length of target latitude and longitude from the finest grid
        tLatS = SIZE(this%fineLat)
        tLonS = SIZE(this%fineLon)
        tSize = tLatS * tLonS
        ! Allocate target latitude and longitude array
        ALLOCATE(tLatlon(tSize, 2), STAT=status)
        ! Assign this array, latitude, longitude
        tLatlon(:, 1) = RESHAPE((SPREAD(this%fineLat, 2, tLonS)), (/tSize/))
        tLatlon(:, 2) = RESHAPE((SPREAD(this%fineLon, 1, tLatS)), (/tSize/))
        ! Slint function is based on radian
        tLatlon = tLatlon * degree2radian

        ! Get source latitude and longitude from observation
        sSize = UBOUND(infoData, 1)
        ! Allocate and assign source laitute and longitude array
        ALLOCATE(sLatlon(sSize, 2), STAT=status)
        sLatlon = infoData(:, 1:2)
        ! Convert degree to radian
        sLatlon = sLatlon * degree2radian

        ! Use slint function to interpolate source to target
        CALL CPU_TIME(t1)
        CALL slint_init(sLatlon, sSize, tLatlon, tSize)
        CALL CPU_TIME(t2)
        PRINT *, "************LoadingNC: running time of slint in interpolate_data:", t2-t1
        ! Interpolate data
        this%numInfo = UBOUND(infoData, 2) - 2
        ALLOCATE(tempData(tSize), interpData(this%numInfo, tLatS, tLonS), STAT=status)
        DO i = 1, this%numInfo
            DO j = 1, tSize
                slintCoeff = tgt_grid%coeffs(:, j)
                slintIdx = tgt_grid%nn(:, j)
                slintValue = 0.0D0
                Do k = 1, 3
                    IF(slintIdx(k) .GT. 0 .AND. slintIdx(k) .LE. sSize .AND. slintCoeff(k) .GT. 0) &
                            slintValue = slintValue + slintCoeff(k)*infoData(slintIdx(k), i+2)
                END DO

                IF (slintValue .LE. 0) THEN
                    tempData(j) = Missing
                ELSE
                    tempData(j) = slintValue
                END IF
            END DO
            interpData(i, :, :) = RESHAPE(tempData, (/tLatS, tLonS/))
        END DO
        IF (ALLOCATED(tLatLon)) DEALLOCATE(tLatLon, STAT=status)
        IF (ALLOCATED(sLatlon)) DEALLOCATE(sLatlon, STAT=status)
        IF (ALLOCATED(tempData)) DEALLOCATE(tempData, STAT=status)
    END SUBROUTINE interpolate_data

    ! Build coefficient matrix from nc file
    SUBROUTINE rebuild_coeff(this, nc_id, posNLen, fileId, dataId, coeffId)
        IMPLICIT NONE
        CLASS(satelliteRecon_t), INTENT(INOUT) :: this
        INTEGER(i_kind), INTENT(IN) :: nc_id, posNLen, fileId, dataId, coeffId
        ! Local variables
        CHARACTER(LEN=2) :: idStr
        INTEGER(i_kind) :: status, dim_id, var_dim, &
                           dim0, dim1, i, var_id
        REAL(r_kind), ALLOCATABLE :: comValue(:), coeffMatrix(:, :)
        INTEGER(i_kind), ALLOCATABLE :: comIdx(:, :), levelIdx(:, :)

        ! Convert data id from number to string
        IF (dataId < 10) THEN
            WRITE(idStr, "(A1I1)") "0", dataId
        ELSE
            WRITE(idStr, "(I2)") dataId
        END IF

        ! Get dimensions
        status = NF90_INQ_DIMID(nc_id, "levelP1", dim_id)
        status = NF90_INQUIRE_DIMENSION(nc_id, dim_id, LEN=dim0)
        status = NF90_INQ_DIMID(nc_id, "compressionN"//TRIM(idStr), dim_id)
        status = NF90_INQUIRE_DIMENSION(nc_id, dim_id, LEN=dim1)
        ! Get variables
        ALLOCATE(comValue(dim1), comIdx(dim1, posNLen), &
                 levelIdx(dim0, posNLen), STAT=status)
        status = NF90_INQ_VARID(nc_id, "comValue"//TRIM(idStr), var_id)
        status = NF90_GET_VAR(nc_id, var_id, comValue)
        status = NF90_INQ_VARID(nc_id, "comIdx"//TRIM(idStr), var_id)
        status = NF90_GET_VAR(nc_id, var_id, comIdx)
        status = NF90_INQ_VARID(nc_id, "levelIdx"//TRIM(idStr), var_id)
        status = NF90_GET_VAR(nc_id, var_id, levelIdx)

        ! Reconstruct satellite data
        ! Allocate coefficient matrix
        ALLOCATE(coeffMatrix(SUM(levelIdx(:, 1)), SUM(levelIdx(:, 2))), STAT=status)
        coeffMatrix = 0.0D0
        ! Put compressed values to coefficient matrix
        DO i = 1, dim1
            coeffMatrix(comIdx(i, 1), comIdx(i, 2)) = comValue(i)
        END DO
        ! Convert coefficient matrix to list
        CALL this%allCoeffs(fileId)%dataCoeffs(coeffId)%matrix2coeff(coeffMatrix, levelIdx)

        ! Deallocate variables
        IF (ALLOCATED(comValue)) DEALLOCATE(comValue, STAT=status)
        IF (ALLOCATED(comIdx)) DEALLOCATE(comIdx, STAT=status)
        IF (ALLOCATED(levelIdx)) DEALLOCATE(levelIdx, STAT=status)
        IF (ALLOCATED(coeffMatrix)) DEALLOCATE(coeffMatrix, STAT=status)
    END SUBROUTINE rebuild_coeff

    SUBROUTINE get_coeff_id(this, curTime, coeffIdx)
        IMPLICIT NONE
        CLASS(satelliteRecon_t), INTENT(IN) :: this
        REAL(r_kind), INTENT(IN) :: curTime
        INTEGER(i_kind), INTENT(OUT) :: coeffIdx
        ! Local variables
        INTEGER(i_kind) :: i

        ! If the observation time does not equal any state%sg%tt, then assign 0
        coeffIdx = 0
        DO i = 1, this%numFiles
            IF (this%time(i) .EQ. curTime) THEN
                coeffIdx = i
                EXIT
            END IF
        END DO
    END SUBROUTINE get_coeff_id

    SUBROUTINE sat_deallocate(this)
        IMPLICIT NONE
        CLASS(satelliteRecon_t), INTENT(INOUT) :: this
        INTEGER(i_kind) :: i, j
        REAL(r_kind) :: t1, t2

        CALL CPU_TIME(t1)
        IF (ALLOCATED(this%varNames)) DEALLOCATE(this%varNames)
        IF (ALLOCATED(this%chanList)) DEALLOCATE(this%chanList)
        IF (ALLOCATED(this%latlon_times)) DEALLOCATE(this%latlon_times)
        IF (ALLOCATED(this%fineLat)) DEALLOCATE(this%fineLat)
        IF (ALLOCATED(this%fineLon)) DEALLOCATE(this%fineLon)
        IF (ALLOCATED(this%time)) DEALLOCATE(this%time)

        IF (ALLOCATED(this%allCoeffs)) THEN
            DO i = 1, SIZE(this%allCoeffs)
                IF (ALLOCATED(this%allCoeffs(i)%cfileName)) DEALLOCATE(this%allCoeffs(i)%cfileName)
                IF (ALLOCATED(this%allCoeffs(i)%ifileName)) DEALLOCATE(this%allCoeffs(i)%ifileName)
                IF (ALLOCATED(this%allCoeffs(i)%fineInfo)) DEALLOCATE(this%allCoeffs(i)%fineInfo)
                IF (ALLOCATED(this%allCoeffs(i)%dataRange)) DEALLOCATE(this%allCoeffs(i)%dataRange)
                IF (ALLOCATED(this%allCoeffs(i)%dataDecomFlag)) DEALLOCATE(this%allCoeffs(i)%dataDecomFlag)
                IF (ALLOCATED(this%allCoeffs(i)%missValPos)) DEALLOCATE(this%allCoeffs(i)%missValPos)
                IF (ALLOCATED(this%allCoeffs(i)%dataCoeffs)) THEN
                    DO j = 1, SIZE(this%allCoeffs(i)%dataCoeffs)
                        CALL this%allCoeffs(i)%dataCoeffs(j)%release()
                    END DO
                    DEALLOCATE(this%allCoeffs(i)%dataCoeffs)
                END IF
            END DO
            DEALLOCATE(this%allCoeffs)
        END IF
        CALL CPU_TIME(t2)
        PRINT *, " ************Destroy: running time:", t2-t1
    END SUBROUTINE sat_deallocate
END MODULE satelliteRecon_m