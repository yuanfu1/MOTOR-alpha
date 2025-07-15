!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : WRFDA
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jiongming Pang (pang.j.m@hotmail.com), 2021-10-8, @GBA-MWF, Shenzhen
!----------------------------------------------------------------------------------------

MODULE InterpolateTime_m
  USE kinds_m
  IMPLICIT NONE

  TYPE InterpolateData_Time_t
    REAL(r_kind), ALLOCATABLE :: data_intperpolate(:, :, :) ! v x h x t, as same as state
    INTEGER(i_kind), ALLOCATABLE :: time_slots_in(:)     ! unix time in sec
    INTEGER(i_kind), ALLOCATABLE :: time_slots_target(:) ! unix time in sec

  CONTAINS
    PROCEDURE, PUBLIC :: IntpInitialize => initialize_sub
    PROCEDURE, PUBLIC :: IntpLinearTime => linearTime_sub
  END TYPE InterpolateData_Time_t

  PRIVATE :: initialize_sub, linearTime_sub

CONTAINS

  SUBROUTINE initialize_sub(this, inData, inTimeSlots, targetTimeSlots)
    CLASS(InterpolateData_Time_t) :: this
    REAL(r_kind), INTENT(IN)  :: inData(:, :, :)
    INTEGER(i_kind), INTENT(IN)  :: inTimeSlots(:), targetTimeSlots(:)
    INTEGER(i_kind) :: shape_indata(3)

    shape_indata = SHAPE(inData)
    IF (shape_indata(3) /= SIZE(inTimeSlots)) THEN
      PRINT *, "time scale of inData and inTimeSlots:", shape_indata(3), SIZE(inTimeSlots)
      STOP "The dimension of the input time slots doesn't match the input array!"
    END IF

    ALLOCATE (this%data_intperpolate(shape_indata(1), shape_indata(2), SIZE(targetTimeSlots)))

    ALLOCATE (this%time_slots_in(SIZE(inTimeSlots)), this%time_slots_target(SIZE(targetTimeSlots)))

    this%time_slots_in = inTimeSlots
    this%time_slots_target = targetTimeSlots

  END SUBROUTINE initialize_sub

  SUBROUTINE linearTime_sub(this, inData)
    CLASS(InterpolateData_Time_t), INTENT(INOUT) :: this
    REAL(r_kind), INTENT(IN)  :: inData(:, :, :)
    REAL(r_kind), ALLOCATABLE :: delta_inTimeSlots(:)
    REAL(r_kind)    :: delta_left, delta_right
    REAL(r_kind)    :: ratio_left, ratio_right
    INTEGER(i_kind) :: index_left, index_right, index_center
    INTEGER(i_kind) :: shape_indata(3)
    INTEGER(i_kind) :: i, j

    shape_indata = SHAPE(inData)

    DO i = 1, SIZE(this%time_slots_target)

      ALLOCATE (delta_inTimeSlots(SIZE(this%time_slots_in)))
      delta_inTimeSlots(:) = this%time_slots_in(:) - this%time_slots_target(i)

      ! GET THE ADJOINT TIME OF THE TARGET TIME, THEN INTERPOLATE TO THE TARGET TIME
      ! two conditions, firstly, the target time is the same as the input time,
      ! look for the index (index_center), then assignt the value of INPUT to TARGET;
      ! secondly, get the indexes of the two adjoint time, then interpolate.
      index_center = 0
      DO j = 1, SIZE(this%time_slots_in)
        IF (delta_inTimeSlots(j) .EQ. 0.0D0) THEN
          index_center = j
          EXIT
        END IF
      END DO

      IF (index_center .NE. 0) THEN
        this%data_intperpolate(:, :, i) = inData(:, :, index_center)
      ELSE
        ! get the value and index of the left adjoint time
        delta_left = MAXVAL(delta_inTimeSlots, MASK=(delta_inTimeSlots < 0), DIM=1)
        index_left = MAXLOC(delta_inTimeSlots, MASK=(delta_inTimeSlots < 0), DIM=1)

        ! get the value and index of the right adjoint time
        delta_right = MINVAL(delta_inTimeSlots, MASK=(delta_inTimeSlots > 0), DIM=1)
        index_right = MINLOC(delta_inTimeSlots, MASK=(delta_inTimeSlots > 0), DIM=1)

        ! calculate the ratio of the two adjoint time based on time-distance
        ASSOCIATE (dist_left => ABS(delta_left), dist_right => ABS(delta_right))
          ratio_left = dist_right / (dist_left + dist_right)
          ratio_right = dist_left / (dist_left + dist_right)
        END ASSOCIATE

        ! interpolate to the target time
        this%data_intperpolate(:, :, i) = ratio_left * inData(:, :, index_left) + ratio_right * inData(:, :, index_right)
      END IF

      DEALLOCATE (delta_inTimeSlots)

    END DO

  END SUBROUTINE linearTime_sub

  ! SUBROUTINE IntpLinearIndex(ind_origin, ind_target, intp_ind_origin, intp_ind_target)
  SUBROUTINE IntpLinearIndex(num_origin, num_target, intp_ind_origin, intp_ind_target)
    IMPLICIT NONE

    INTEGER(i_kind), INTENT(IN) :: num_origin, num_target
    INTEGER(i_kind), INTENT(INOUT) :: intp_ind_origin(:)
    INTEGER(i_kind), INTENT(INOUT) :: intp_ind_target(:)
    INTEGER(i_kind) :: num_final, int_length, i

    INTEGER(i_kind), ALLOCATABLE :: ind_origin(:), ind_target(:)

    ALLOCATE (ind_origin(num_origin), ind_target(num_target))
    DO i = 1, num_origin
      ind_origin(i) = i
    END DO
    DO i = 1, num_target
      ind_target(i) = i
    END DO
    ! num_origin = SIZE(ind_origin)
    ! num_target = SIZE(ind_target)
    num_final = MAXVAL((/num_origin, num_target/))

    IF (num_target >= num_origin) THEN
      int_length = (num_final - 1) / (num_origin - 1); 
    ELSE
      int_length = (num_final - 1) / (num_target - 1); 
    END IF

    IF (num_target >= num_final) THEN
      intp_ind_target = ind_target
      DO i = 1, num_origin
        intp_ind_origin(i) = ind_origin(1) + int_length * (i - 1)
      END DO
    ELSE
      intp_ind_origin = ind_origin
      DO i = 1, num_target
        intp_ind_target(i) = ind_target(1) + int_length * (i - 1)
      END DO
    END IF

  END SUBROUTINE IntpLinearIndex

  SUBROUTINE IntpLinearTimeSlots(inTimeSlots, targetTimeSlots)
    INTEGER(i_kind), INTENT(IN)    :: inTimeSlots(:)
    INTEGER(i_kind), INTENT(INOUT) :: targetTimeSlots(:)
    INTEGER(i_kind) :: numTimeSlots, deltaTime, i

    numTimeSlots = SIZE(targetTimeSlots)
    deltaTime = (MAXVAL(inTimeSlots, DIM=1) - MINVAL(inTimeSlots, DIM=1)) / (numTimeSlots - 1)

    targetTimeSlots(1) = MINVAL(inTimeSlots, DIM=1)
    targetTimeSlots(numTimeSlots) = MAXVAL(inTimeSlots, DIM=1)
    DO i = 2, numTimeSlots - 1
      targetTimeSlots(i) = targetTimeSlots(1) + deltaTime * (i - 1)
    END DO

  END SUBROUTINE IntpLinearTimeSlots

END MODULE InterpolateTime_m
