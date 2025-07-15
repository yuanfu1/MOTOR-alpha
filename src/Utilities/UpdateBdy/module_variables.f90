MODULE module_variables
  IMPLICIT NONE

  INTEGER, PARAMETER :: P_XSB = 1
  INTEGER, PARAMETER :: P_XEB = 2
  INTEGER, PARAMETER :: P_YSB = 3
  INTEGER, PARAMETER :: P_YEB = 4

  REAL, DIMENSION(:, :, :), ALLOCATABLE   :: pi, th, u, v, q
  REAL, DIMENSION(:, :, :, :), ALLOCATABLE   :: pi_b, pi_b_old, pi_b_1, pi_b_2, pi_bt, pi_bt_old
  REAL, DIMENSION(:, :, :, :), ALLOCATABLE   :: th_b, th_b_old, th_b_1, th_b_2, th_bt, th_bt_old
  REAL, DIMENSION(:, :, :, :), ALLOCATABLE   :: u_b, u_b_old, u_b_1, u_b_2, u_bt, u_bt_old
  REAL, DIMENSION(:, :, :, :), ALLOCATABLE   :: v_b, v_b_old, v_b_1, v_b_2, v_bt, v_bt_old
  REAL, DIMENSION(:, :, :, :), ALLOCATABLE   :: q_b, q_b_old, q_b_1, q_b_2, q_bt, q_bt_old

END MODULE
