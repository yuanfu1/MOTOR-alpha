!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.SatelliteDecomposition
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Ting Shu
! FUNCTION          : Wavelet-based satellite decomposition
! VERSION           : V 0.0
! HISTORY           :
!   Created by Ting Shu (shuting@gbamwf.com), 2023/08/08, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements Wavelet-based satellite decomposition.
MODULE satelliteDecom_m
    USE kinds_m, ONLY: r_kind, i_kind
    USE YAMLRead_m
    USE netcdf
    USE parameters_m, ONLY: degree2radian, Missing
    USE slint, ONLY: slint_init, bilinear_interp, tgt_grid
    USE wavelet_m, ONLY: wavelet_t
    USE satelliteW_m, ONLY: get_latlon, insert_sort, compress_bi_array, quick_sort
    IMPLICIT NONE
    
    TYPE :: satelliteDecom_t
        CHARACTER(LEN=1024) :: configFile
        INTEGER(i_kind) :: mgStart, mgEnd, decomLevel
        REAL(r_kind) :: keepRate
        REAL(r_kind), DIMENSION(:), ALLOCATABLE :: fineLat, fineLon
        LOGICAL :: processErr

        CONTAINS
        PROCEDURE, PUBLIC :: initialize => sat_decom_ini
        PROCEDURE, PUBLIC :: decom => sat_decom
        PROCEDURE, PUBLIC :: destroy => sat_deallocate
    END TYPE satelliteDecom_t
    PRIVATE :: sat_decom_ini, sat_decom, sat_deallocate

    CONTAINS

    SUBROUTINE sat_decom_ini(this, configFile)
        IMPLICIT NONE
        CLASS(satelliteDecom_t), INTENT(INOUT) :: this
        CHARACTER(LEN=1024), INTENT(IN) :: configFile
        ! Local variables
        INTEGER(i_kind), DIMENSION(:), ALLOCATABLE :: latlon_times
        REAL(r_kind), DIMENSION(:), ALLOCATABLE :: latlon_domain
        INTEGER(i_kind) :: status, i, latlon_len(2)
        REAL(r_kind) :: temp_d, dxdy(2), new_latlon_domain(4)

        ! Assign configFile
        this%configFile = configFile
        ! Get wavelet compression keep rate from yaml file
        status = yaml_get_var(TRIM(configFile), 'Wavelet', 'keepRate', this%keepRate)
        ! Get the flag wheather needs to process error via wavelet
        status = yaml_get_var(TRIM(configFile), 'Wavelet', 'processErr', this%processErr)
        ! Get start and end grid id
        status = yaml_get_var(TRIM(configFile), 'geometry', 'mgStart', this%mgStart)
        status = yaml_get_var(TRIM(configFile), 'geometry', 'mgEnd', this%mgEnd)
        ! Wavelet decomposition level
        this%decomLevel = this%mgEnd - this%mgStart + 1
        ! Get latitude and longitude domain and times
        status = yaml_get_var(TRIM(configFile), 'latlon_grid', 'domain_latlon', latlon_domain)
        status = yaml_get_var(TRIM(configFile), 'Wavelet', 'latlon_times', latlon_times)
        ! Compute the latitude and longitude of the finest grid
        CALL get_latlon(this%mgStart, this%mgEnd, latlon_times, latlon_domain, this%fineLat, this%fineLon)
    END SUBROUTINE sat_decom_ini

    SUBROUTINE sat_decom(this, sat_latlon, sat_time, sat_data, sat_errors, saveFileName, t0)
        IMPLICIT NONE
        CLASS(satelliteDecom_t), INTENT(INOUT) :: this
        REAL(r_kind), INTENT(IN) :: sat_latlon(:, :), sat_data(:, :), sat_errors(:, :), sat_time, t0
        CHARACTER(LEN=*), INTENT(IN) :: saveFileName
        ! Local variables
        INTEGER(i_kind) :: tLatS, tLonS, tSize, sSize, dLen, n, wdLen
        INTEGER(i_kind) :: i, j, k, l, status
        REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: tLatlon, sLatlon, data_errors, sat2d
        REAL(r_kind), DIMENSION(:), ALLOCATABLE :: tempData
        ! Coefficient related variables
        TYPE(wavelet_t) :: coeff
        REAL(r_kind), ALLOCATABLE :: coeffMatrix(:, :), coeff1d(:), compressValue(:), dataMinMax(:, :)
        INTEGER(i_kind), ALLOCATABLE :: coeffIdx(:, :), corder(:), compressIdx(:, :)
        REAL(r_kind) :: threshold, cur_var(3), slintCoeff(3)
        ! Nc file related variables
        CHARACTER(LEN=2) :: var_name
        INTEGER(i_kind) :: nc_id, dim_id0, dim_id1, dim_id2, var_id, miss_i, miss_count
        ! Missing data flag array, missing value positions array
        INTEGER(i_kind), ALLOCATABLE :: data_decom_flag(:), missValPos(:), com_miss(:)
        REAL(r_kind) :: t1, t2

        PRINT *, "************Start calculate running time of each step************************"
        CALL CPU_TIME(t1)
        PRINT *, "*********Running time of declaring variables: ", t1-t0
        ! Assign satellite data
        sSize = UBOUND(sat_latlon, 2)
        dLen = UBOUND(sat_data, 2)

        IF (this%processErr) THEN
            wdLen = dLen*2
            ! Combine satellite data and errors in one array for decomposition one by one
            ALLOCATE(data_errors(sSize, wdLen), STAT=status)
            data_errors(:, 1:dLen) = sat_data
            data_errors(:, dLen+1:wdLen) = sat_errors
        ELSE
            wdLen = dLen
            ! Just process data, no error
            ALLOCATE(data_errors(sSize, wdLen), STAT=status)
            data_errors = sat_data
        END IF

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
        ! Allocate and assign source laitute and longitude array
        ALLOCATE(sLatlon(sSize, 2), STAT=status)
        sLatlon(:, 1) = sat_latlon(1, :)
        sLatlon(:, 2) = sat_latlon(2, :)
        ! Convert degree to radian
        ! sLatlon = sLatlon * degree2radian

        ! Use slint function to interpolate source to target
        CALL CPU_TIME(t1)
        CALL slint_init(sLatlon, sSize, tLatlon, tSize)
        CALL CPU_TIME(t2)
        PRINT *, "*********Running time of calling slint: ", sSize, tSize, t2-t1

        ! Create a nc file
        status = NF90_CREATE(TRIM(saveFileName), NF90_NETCDF4, nc_id)
        status = NF90_DEF_DIM(nc_id, "levelP1", this%decomLevel+1, dim_id0)
        status = NF90_DEF_DIM(nc_id, "posN", 2, dim_id2)
        ! Allocate interpolation data
        ALLOCATE(tempData(tSize), sat2d(tLatS, tLonS), STAT=status)
        ! Original data minimal and maximal value
        ALLOCATE(dataMinMax(wdLen, 2), data_decom_flag(wdLen), STAT=status)
        ! Missing value positions array, 1 for missing, 0 for presenting
        ALLOCATE(missValPos(tSize), STAT=status)
        dataMinMax = 0.0D0
        data_decom_flag = 0
        PRINT *, "========= Start decomposing satellite data one by one ======"
        DO i = 1, wdLen
            missValPos = 1
            ! Print minimum and maximum of 1D array after remove missing value
            ! Determine the data is missing or not
            IF (i .LE. dLen) THEN
                IF (COUNT(sat_data(:, i) .EQ. Missing) .NE. sSize) THEN
                    data_decom_flag(i) = i ! for data
                    IF(this%processErr) data_decom_flag(i+dLen) = i+dLen  ! for error
                END IF
            END IF
            ! IF the data is missing, then do not decompose data and error
            IF (data_decom_flag(i) .NE. 0) THEN
                CALL CPU_TIME(t1)
                ! Interpolate data
                DO j = 1, tSize
                    miss_count = 0
                    slintCoeff = tgt_grid%coeffs(:, j)
                    DO miss_i = 1, 3
                        ! Sometimes, slint function does not assign tgt_grid%nn appropriately
                        IF (tgt_grid%nn(miss_i, j) .GT. sSize .OR. &
                        tgt_grid%nn(miss_i, j) .LE. 0 .OR. &
                        slintCoeff(miss_i) .LE. 0) THEN
                            slintCoeff(miss_i) = 0.0D0
                            cur_var(miss_i) = 0.0D0
                            miss_count = miss_count + 1
                        ELSE
                            cur_var(miss_i) = data_errors(tgt_grid%nn(miss_i, j), i)
                            IF (ABS(ABS(cur_var(miss_i)) - Missing) .LE. 1) THEN
                                cur_var(miss_i) = 0.0D0
                                miss_count = miss_count + 1
                            END IF
                        END IF
                    END DO
                    ! The variable is missing
                    IF (i .LE. dLen .AND. miss_count .EQ. 3) THEN
                        missValPos(j) = 1
                        tempData(j) = 0.0D0
                    ELSE
                        missValPos(j) = 0
                        tempData(j) = SUM(slintCoeff * cur_var)
                    END IF
                END DO
                CALL CPU_TIME(t2)
                PRINT *, "*******Running time of interpolation: ", t2-t1, i
                sat2d = RESHAPE(tempData, (/tLatS, tLonS/))
                ! Get minimum and maximum of original data
                dataMinMax(i, 1) = MINVAL(sat2d)
                dataMinMax(i, 2) = MAXVAL(sat2d)
                ! Decompose 2d satellite data
                CALL coeff%decom(sat2d, this%decomLevel)
                CALL CPU_TIME(t1)
                PRINT *, "*******Running time of wavelet decomposition: ", t1-t2, i
                ! Convert coefficient to matrix
                CALL coeff%coeff2matrix(coeffMatrix, coeffIdx)
                CALL CPU_TIME(t2)
                PRINT *, "*******Running time of converting coefficient to matrix: ", t2-t1, i
                ! Sort coefficient matrix
                n = SIZE(coeffMatrix)
                ALLOCATE(coeff1d(n), corder(n), STAT=status)
                coeff1d = RESHAPE(coeffMatrix, (/n/))
                coeff1d = ABS(coeff1d)
                corder = (/(k, k=1, n)/)
                ! CALL insert_sort(coeff1d, corder)
                CALL quick_sort(coeff1d, corder, 1, n)
                CALL CPU_TIME(t1)
                PRINT *, "*******Running time of coefficient matrix decompression1, sorting: ", t1-t2, i
                ! Calculate the threshold for the keep rate
                threshold = coeff1d(FLOOR((1 - this%keepRate)*n))
                ! Make the coefficients tobe zeros where smaller than and equal to threshold
                WHERE(ABS(coeffMatrix) .LE. threshold) coeffMatrix = 0

                ! Comprese coefficient
                n = COUNT(coeffMatrix .NE. 0)
                ALLOCATE(compressValue(n), compressIdx(n, 2), STAT=status)
                status = 0
                DO k = 1, UBOUND(coeffMatrix, 1)
                    DO l = 1, UBOUND(coeffMatrix, 2)
                        IF (coeffMatrix(k, l) .NE. 0) THEN
                            status = status + 1
                            compressValue(status) = coeffMatrix(k, l)
                            compressIdx(status, :) = (/k, l/)
                        END IF
                    END DO
                END DO
                CALL CPU_TIME(t2)
                PRINT *, "*******Running time of coefficient matrix decompression2, converting: ", t2-t1, i
                ! Save coefficient into nc file
                if (i < 10) THEN
                    WRITE(var_name, "(A1I1)") "0", i
                ELSE
                    WRITE(var_name, "(I2)") i
                END IF
                status = NF90_DEF_VAR(nc_id, "levelIdx"//TRIM(var_name), NF90_INT, (/dim_id0, dim_id2/), var_id)
                status = NF90_PUT_VAR(nc_id, var_id, coeffIdx)

                status = NF90_DEF_DIM(nc_id, "compressionN"//TRIM(var_name), n, dim_id1)
                status = NF90_DEF_VAR(nc_id, "comValue"//TRIM(var_name), NF90_FLOAT, (/dim_id1/), var_id)
                status = NF90_PUT_VAR(nc_id, var_id, compressValue)
                status = NF90_DEF_VAR(nc_id, "comIdx"//TRIM(var_name), NF90_INT, (/dim_id1, dim_id2/), var_id)
                status = NF90_PUT_VAR(nc_id, var_id, compressIdx)

                ! Save missing value flag vector
                IF (i .LE. dlen) THEN
                    CALL compress_bi_array(missValPos, com_miss)
                    status = NF90_DEF_DIM(nc_id, "comMissN"//TRIM(var_name), SIZE(com_miss), dim_id1)
                    status = NF90_DEF_VAR(nc_id, "comMissArray"//TRIM(var_name), NF90_SHORT, (/dim_id1/), var_id)
                    status = NF90_PUT_VAR(nc_id, var_id, com_miss)
                    IF (ALLOCATED(com_miss)) DEALLOCATE(com_miss, STAT=status)
                END IF
                CALL CPU_TIME(t1)
                PRINT *, "*******Running time of saving data: ", t1-t2, i
                ! Deallocate variables
                CALL coeff%release()
                IF (ALLOCATED(coeffMatrix)) DEALLOCATE(coeffMatrix, STAT=status)
                IF (ALLOCATED(coeffIdx)) DEALLOCATE(coeffIdx, STAT=status)
                IF (ALLOCATED(coeff1d)) DEALLOCATE(coeff1d, STAT=status)
                IF (ALLOCATED(corder)) DEALLOCATE(corder, STAT=status)
                IF (ALLOCATED(compressValue)) DEALLOCATE(compressValue, STAT=status)
                IF (ALLOCATED(compressIdx)) DEALLOCATE(compressIdx, STAT=status)
                CALL CPU_TIME(t2)
                PRINT *, "*******Running time of deallocating do-loop variables: ", t2-t1, i
            END IF
        END DO
        IF (this%processErr) THEN
            PRINT "(I3AI3A)", COUNT(data_decom_flag .NE. 0)/2, " of ", dLen, " varibales have been decomposed."
        ELSE
            PRINT "(I3AI3A)", COUNT(data_decom_flag .NE. 0), " of ", dLen, " varibales have been decomposed."
        END IF
        PRINT *, "========= Finish decomposing satellite data ======"
        
        CALL CPU_TIME(t1)
        ! Save minimum and maximum of original data
        status = NF90_DEF_DIM(nc_id, "dlen2times", wdLen, dim_id1)
        status = NF90_DEF_VAR(nc_id, "dataRange", NF90_FLOAT, (/dim_id1, dim_id2/), var_id)
        status = NF90_PUT_VAR(nc_id, var_id, dataMinMax)
        ! Save decomposition data flag
        status = NF90_DEF_VAR(nc_id, "dataDecomFlag", NF90_INT, (/dim_id1/), var_id)
        status = NF90_PUT_VAR(nc_id, var_id, data_decom_flag)

        ! Add time information to the nc file using attribute
        status = NF90_PUT_ATT(nc_id, NF90_GLOBAL, "observationTime", sat_time)
        status = NF90_CLOSE(nc_id)
        CALL CPU_TIME(t2)
        PRINT *, "*******Running time of saving variables into ncfile: ", t2-t1

        IF (ALLOCATED(dataMinMax)) DEALLOCATE(dataMinMax, STAT=status)
        IF (ALLOCATED(missValPos)) DEALLOCATE(missValPos, STAT=status)
    END SUBROUTINE sat_decom

    SUBROUTINE sat_deallocate(this)
        IMPLICIT NONE
        CLASS(satelliteDecom_t), INTENT(INOUT) :: this
        
        IF(ALLOCATED(this%fineLat)) DEALLOCATE(this%fineLat)
        IF(ALLOCATED(this%fineLon)) DEALLOCATE(this%fineLon)
    END SUBROUTINE sat_deallocate
END MODULE satelliteDecom_m