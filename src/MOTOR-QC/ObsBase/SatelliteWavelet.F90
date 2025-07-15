!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.SatelliteWavelet
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Ting Shu
! FUNCTION          : Wavelet-based satellite processing tool module
! VERSION           : V 0.0
! HISTORY           :
!   Created by Ting Shu (shuting@gbamwf.com), 2023/08/09, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements Wavelet-based satellite processing tools.
MODULE satelliteW_m
    USE kinds_m, ONLY: r_kind, i_kind
    USE parameters_m, ONLY: degree2radian
    USE wavelet_m, ONLY: wavelet_t
    IMPLICIT NONE

CONTAINS
    SUBROUTINE get_latlon(mgStart, mgEnd, latlon_times, latlon_domain, fineLat, fineLon)
        IMPLICIT NONE
        INTEGER(i_kind), INTENT(IN) :: mgStart, mgEnd, latlon_times(:)
        REAL(r_kind), INTENT(IN) :: latlon_domain(:)
        REAL(r_kind), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: fineLat, fineLon
        ! Local variables
        INTEGER(i_kind) :: i, status, latlon_len(2)
        REAL(r_kind), DIMENSION(4) :: new_latlon_domain
        
        REAL(r_kind) :: temp_d, dxdy(2)

        ! Calculate the dx and dy for mgStart grid
        DO i = 1, 2
            temp_d = (latlon_domain(i*2) - latlon_domain(i*2-1)) / (2 ** mgStart * latlon_times(i))
            ! Latitude and longitude of new domain
            new_latlon_domain(i*2-1) = latlon_domain(i*2-1) - temp_d
            new_latlon_domain(i*2) = latlon_domain(i*2) + temp_d
            ! Dx and dy for high fine grid (final grid)
            dxdy(i) = (latlon_domain(i*2-1) - new_latlon_domain(i*2-1)) / (2 ** (mgEnd - mgStart))
            ! Length of final latitude and longitude
            latlon_len(i) = 2 ** (mgEnd - mgStart + 1) + 2 ** mgEnd * latlon_times(i) + 1
        END DO

        ! Allocate fineLat and fineLon
        ALLOCATE(fineLat(latlon_len(1)), fineLon(latlon_len(2)), STAT=status)
        fineLat = new_latlon_domain(1) + (/(i, i=0, latlon_len(1)-1)/) * dxdy(1)
        fineLon = new_latlon_domain(3) + (/(i, i=0, latlon_len(2)-1)/) * dxdy(2)
    END SUBROUTINE get_latlon

    SUBROUTINE grid_idx(mgStart, mgEnd, latlon_times, curMgId, latIdx, lonIdx)
        IMPLICIT NONE
        INTEGER(i_kind), INTENT(IN) :: mgStart, mgEnd, curMgId, latlon_times(:)
        INTEGER(i_kind), DIMENSION(:), ALLOCATABLE, INTENT(OUT) :: latIdx, lonIdx
        ! Local variables
        INTEGER(i_kind) :: i, status, startIdx, mgTimes, latLen, lonLen

        ! Grid start index
        startIdx = 2 ** (mgEnd - mgStart) + 1
        ! Current grid times of the lowest-resolution grid
        mgTimes = 2 ** (mgEnd - curMgId)
        latLen = 2 ** (curMgId - 1) * latlon_times(1) + 2
        lonLen = 2 ** (curMgId - 1) * latlon_times(2) + 2
        ! Allocate latitude and longitude index array
        ALLOCATE(latIdx(latLen), lonIdx(lonLen), STAT=status)
        latIdx = startIdx - mgTimes + (/(i, i=0, latLen-1)/) * 2 * mgTimes
        lonIdx = startIdx - mgTimes + (/(i, i=0, lonLen-1)/) * 2 * mgTimes
    END SUBROUTINE grid_idx

    SUBROUTINE insert_sort(array, idex)
        IMPLICIT NONE
        REAL(r_kind), INTENT(INOUT) :: array(:)
        INTEGER(i_kind), INTENT(INOUT) :: idex(:)
        ! Local variables
        REAL(r_kind) :: temp
        INTEGER(i_kind) :: i, j, d_len, temp_i

        ! Get length of array
        d_len = SIZE(array)
        idex = (/(i, i=1, d_len)/)
        ! two loops for comparing array elements
        DO i = 2, d_len
            ! Get current value and index
            temp = array(i)
            temp_i = idex(i)
            ! From i-1 to 1
            DO j = i-1, 1, -1
                ! If array(i) is the biggest,
                ! compared with its previous elements, done
                IF(array(j) .LE. temp) EXIT
                ! If not, then move, 
                ! insert the temp to its position
                array(j+1) = array(j)
                idex(j+1) = idex(j)
            END DO
            ! j = 0
            array(j+1) = temp
            idex(j+1) = temp_i
        END DO
    END SUBROUTINE insert_sort

    RECURSIVE SUBROUTINE quick_sort(array, idx, first, last)
        IMPLICIT NONE
        REAL(r_kind), DIMENSION(:), INTENT(INOUT) :: array
        INTEGER(i_kind), DIMENSION(:), INTENT(INOUT) :: idx
        INTEGER(i_kind), INTENT(IN) :: first, last
        REAL(r_kind) :: x, t
        INTEGER(i_kind) :: i, j, tmp_idx

        IF(SIZE(array) .NE. SIZE(idx)) THEN
            PRINT *, "Wrong with calling quick_sort, please input the same size of array and index."
            RETURN
        END IF

        x = array((first+last)/2)
        i = first
        j = last
        DO
            ! IF the middle x is larger than left varaible
            ! Low pointer moves to right
            DO WHILE(array(i) < x)
                i = i + 1
            END DO
            ! IF the middle x is smaller than right varaible
            ! High pointer moves to left
            DO WHILE(x < array(j))
                j = j - 1
            END DO
            ! If low pointer moves to the right of the high pointer
            ! sorting is done
            IF (i >= j) EXIT
            ! Else, swap low and high's array and index
            t = array(i)
            array(i) = array(j)
            array(j) = t
            tmp_idx = idx(i)
            idx(i) = idx(j)
            idx(j) = tmp_idx

            i = i + 1
            j = j - 1
        END DO
        IF (first < i - 1) CALL quick_sort(array, idx, first, i-1)
        IF (j+1 < last) CALL quick_sort(array, idx, j+1, last)
    END SUBROUTINE quick_sort

    ! Compress binary array, for example
    ! 10011101 -> 112311
    SUBROUTINE compress_bi_array(array, com_array)
        IMPLICIT NONE
        INTEGER(i_kind), INTENT(IN), DIMENSION(:) :: array
        INTEGER(i_kind), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: com_array
        ! Local variables
        INTEGER(i_kind) :: a_len, com_len, status, i, temp_value, temp_count, temp_idx
        INTEGER(i_kind), ALLOCATABLE :: temp_array(:)

        ! Get length of the input array
        a_len = SIZE(array)
        ! Allocate temporal array
        ALLOCATE(temp_array(a_len), STAT=status)

        ! Compress array one by one
        ! Assign the first value
        temp_array(1) = array(1)
        temp_value = array(1)
        temp_count = 0
        temp_idx = 1
        DO i = 1, a_len
            IF (array(i) .NE. 0 .AND. array(i) .NE. 1) THEN
                PRINT *, "Error: compress_bi_array subroutine just processes binary array!"
                RETURN
            END IF
            IF (array(i) .EQ. temp_value) THEN
                temp_count = temp_count + 1
            ELSE
                ! The value of array changes
                temp_idx = temp_idx + 1
                temp_array(temp_idx) = temp_count
                ! 1 -> 0, 0 -> 1
                temp_value = 1 - temp_value
                ! Restart counting
                temp_count = 1
            END IF
        END DO
        ! Save the number of the last value
        IF (array(a_len) .EQ. temp_value) THEN
            temp_idx = temp_idx + 1
            temp_array(temp_idx) = temp_count
        END IF

        ! Save the compressed array
        ALLOCATE(com_array(temp_idx), STAT=status)
        com_array = temp_array(1 : temp_idx)
        
        IF(ALLOCATED(temp_array)) DEALLOCATE(temp_array, STAT=status)
    END SUBROUTINE compress_bi_array

    ! Decompress binary array, for example
    ! 112311 -> 10011101
    SUBROUTINE decompress_bi_array(com_array, array)
        IMPLICIT NONE
        INTEGER(i_kind), INTENT(IN), DIMENSION(:) :: com_array
        INTEGER(i_kind), INTENT(OUT), ALLOCATABLE, DIMENSION(:) :: array
        ! Local variables
        INTEGER(i_kind) :: com_len, a_len, status, i, temp_value, start_idx

        ! Get length of the input compressed array
        com_len = SIZE(com_array)
        ! Calcualte the length of the decompressed array
        a_len = SUM(com_array(2 : com_len))
        ! Allocate decompressed array
        ALLOCATE(array(a_len), STAT=status)

        ! Decompress array one by one
        ! Get the first value
        temp_value = com_array(1)
        start_idx = 1
        DO i = 2, com_len
            array(start_idx : start_idx + com_array(i) - 1) = temp_value
            temp_value = 1 - temp_value
            start_idx = start_idx + com_array(i)
        END DO
    END SUBROUTINE decompress_bi_array
END MODULE satelliteW_m