PROGRAM Test_ECMIO
  USE ECMIO_m
  USE kinds_m, ONLY: i_kind
  IMPLICIT NONE

  ! 定义测试数据
  CHARACTER(LEN=40), DIMENSION(2) :: ECMFileNames = [ &
                                     "/mnt/d/Data/hagongdaAnalysis/20230103.nc", "/mnt/d/Data/hagongdaAnalysis/20230103.nc"]
  CHARACTER(LEN=20), DIMENSION(10) :: ECMVarNames = [ &
                                      "t2m                 ", "q2m                 ", "precipitation       ", &
                                      "u10m                ", "v10m                ", "u200m               ", &
                                      "v200m               ", "u100m               ", "v100m               ", &
                                      "pressure            "]
  INTEGER(i_kind), DIMENSION(6, 2) :: mdate = RESHAPE([ &
                                                      2023, 1, 1, 0, 0, 0, &
                                                      2023, 1, 2, 0, 0, 0], [6, 2])

  ! 创建 ECMIO_t 实例
  TYPE(ECMIO_t) :: ecmio

  ! 构造函数
  ecmio = constructor(ECMFileNames, ECMVarNames, mdate)

  ! 验证数据
  CALL verify_data(ecmio)

  PRINT *, "Test completed successfully."

CONTAINS

  SUBROUTINE verify_data(ecmio)
    TYPE(ECMIO_t), INTENT(IN) :: ecmio
    INTEGER(i_kind) :: i, j, k

    ! 验证时间戳
    DO i = 1, ecmio%nt
      PRINT *, "Time Unix: ", ecmio%time_unix(i)
    END DO

    ! 验证部分数据

    PRINT *, "t2m(", 1, ",", 1, ",", 1, "): ", ecmio%t2m(1, 1, 1)
    PRINT *, "t2m(", 1, ",", 1, ",", 2, "): ", ecmio%t2m(1, 1, 2)
    PRINT *, "q2m(", 1, ",", 1, ",", 1, "): ", ecmio%q2m(1, 1, 1)
    PRINT *, "q2m(", 1, ",", 1, ",", 2, "): ", ecmio%q2m(1, 1, 2)
    PRINT *, "precipitation(", 1, ",", 1, ",", 1, "): ", ecmio%precipitation(1, 1, 1)
    PRINT *, "precipitation(", 1, ",", 1, ",", 2, "): ", ecmio%precipitation(1, 1, 2)
    PRINT *, "u10m(", 1, ",", 1, ",", 1, "): ", ecmio%u10m(1, 1, 1)
    PRINT *, "u10m(", 1, ",", 1, ",", 2, "): ", ecmio%u10m(1, 1, 2)
    PRINT *, "v10m(", 1, ",", 1, ",", 1, "): ", ecmio%v10m(1, 1, 1)
    PRINT *, "v10m(", 1, ",", 1, ",", 2, "): ", ecmio%v10m(1, 1, 2)
    PRINT *, "u200m(", 1, ",", 1, ",", 1, "): ", ecmio%u200m(1, 1, 1)
    PRINT *, "u200m(", 1, ",", 1, ",", 2, "): ", ecmio%u200m(1, 1, 2)
    PRINT *, "v200m(", 1, ",", 1, ",", 1, "): ", ecmio%v200m(1, 1, 1)
    PRINT *, "v200m(", 1, ",", 1, ",", 2, "): ", ecmio%v200m(1, 1, 2)
    PRINT *, "u100m(", 1, ",", 1, ",", 1, "): ", ecmio%u100m(1, 1, 1)
    PRINT *, "u100m(", 1, ",", 1, ",", 2, "): ", ecmio%u100m(1, 1, 2)
    PRINT *, "v100m(", 1, ",", 1, ",", 1, "): ", ecmio%v100m(1, 1, 1)
    PRINT *, "v100m(", 1, ",", 1, ",", 2, "): ", ecmio%v100m(1, 1, 2)
    PRINT *, "pressure(", 1, ",", 1, ",", 1, "): ", ecmio%pressure(1, 1, 1)
    PRINT *, "pressure(", 1, ",", 1, ",", 2, "): ", ecmio%pressure(1, 1, 2)

  END SUBROUTINE verify_data

END PROGRAM Test_ECMIO
