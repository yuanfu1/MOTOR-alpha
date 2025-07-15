MODULE kinds_m
  IMPLICIT NONE
  PRIVATE

! Integer type definitions below

! Integer types
  INTEGER, PARAMETER, PUBLIC  :: i_byte = SELECTED_INT_KIND(1)      ! byte  integer
  INTEGER, PARAMETER, PUBLIC  :: i_short = SELECTED_INT_KIND(4)      ! short integer
  INTEGER, PARAMETER, PUBLIC  :: i_long = SELECTED_INT_KIND(8)      ! long  integer
  INTEGER, PARAMETER, PRIVATE :: llong_t = SELECTED_INT_KIND(16)     ! llong integer
  INTEGER, PARAMETER, PUBLIC  :: i_llong = MAX(llong_t, i_long)

! Expected 8-bit byte sizes of the integer kinds
  INTEGER, PARAMETER, PUBLIC :: num_bytes_for_i_byte = 1
  INTEGER, PARAMETER, PUBLIC :: num_bytes_for_i_short = 2
  INTEGER, PARAMETER, PUBLIC :: num_bytes_for_i_long = 4
  INTEGER, PARAMETER, PUBLIC :: num_bytes_for_i_llong = 8

! Define arrays for default definition
  INTEGER, PARAMETER, PRIVATE :: num_i_kinds = 4
  INTEGER, PARAMETER, DIMENSION(num_i_kinds), PRIVATE :: integer_types = (/ &
                                                         i_byte, i_short, i_long, i_llong/)
  INTEGER, PARAMETER, DIMENSION(num_i_kinds), PRIVATE :: integer_byte_sizes = (/ &
                                                         num_bytes_for_i_byte, num_bytes_for_i_short, &
                                                         num_bytes_for_i_long, num_bytes_for_i_llong/)

! Default values
! **** CHANGE THE FOLLOWING TO CHANGE THE DEFAULT INTEGER TYPE KIND ***
  INTEGER, PARAMETER, PRIVATE :: default_integer = 3  ! 1=byte,
  ! 2=short,
  ! 3=long,
  ! 4=llong
  INTEGER, PARAMETER, PUBLIC  :: i_kind = integer_types(default_integer)
  INTEGER, PARAMETER, PUBLIC  :: num_bytes_for_i_kind = &
                                 integer_byte_sizes(default_integer)

! Real definitions below

! Real types
  INTEGER, PARAMETER, PUBLIC  :: r_single = SELECTED_REAL_KIND(6)  ! single precision
!  INTEGER, PARAMETER, PUBLIC  :: r_double = SELECTED_REAL_KIND(15) ! double precision
  INTEGER, PARAMETER, PUBLIC  :: r_double = SELECTED_REAL_KIND(15) ! double precision
  INTEGER, PARAMETER, PRIVATE :: quad_t = SELECTED_REAL_KIND(20) ! quad precision
  INTEGER, PARAMETER, PUBLIC  :: r_quad = MAX(quad_t, r_double)

! Expected 8-bit byte sizes of the real kinds
  INTEGER, PARAMETER, PUBLIC :: num_bytes_for_r_single = 4
  INTEGER, PARAMETER, PUBLIC :: num_bytes_for_r_double = 8
  INTEGER, PARAMETER, PUBLIC :: num_bytes_for_r_quad = 16

! Define arrays for default definition
  INTEGER, PARAMETER, PRIVATE :: num_r_kinds = 3
  INTEGER, PARAMETER, DIMENSION(num_r_kinds), PRIVATE :: REAL_KINDS = (/ &
                                                         r_single, r_double, r_quad/)
  INTEGER, PARAMETER, DIMENSION(num_r_kinds), PRIVATE :: real_byte_sizes = (/ &
                                                         num_bytes_for_r_single, num_bytes_for_r_double, &
                                                         num_bytes_for_r_quad/)

! Default values
! **** CHANGE THE FOLLOWING TO CHANGE THE DEFAULT REAL TYPE KIND ***
#ifdef _REAL4_
  INTEGER, PARAMETER, PRIVATE :: default_real = 1  ! 1=single,
#endif
#ifdef _REAL8_
  INTEGER, PARAMETER, PRIVATE :: default_real = 2  ! 2=double,
!  INTEGER, PARAMETER, PRIVATE :: default_real = 1  ! 2=double,
#endif
#ifdef _REAL16_
  INTEGER, PARAMETER, PRIVATE :: default_real = 3  ! 3=quad
#endif
  INTEGER, PARAMETER, PUBLIC  :: r_kind = REAL_KINDS(default_real)
  INTEGER, PARAMETER, PUBLIC  :: num_bytes_for_r_kind = &
                                 real_byte_sizes(default_real)

END MODULE kinds_m
