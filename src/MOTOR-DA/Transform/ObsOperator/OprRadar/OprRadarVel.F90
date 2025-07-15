!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE OprRadarVel_m
  USE State_m, ONLY: State_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE ObsSet_m, ONLY: ObsSet_t
  USE M2OBase_m, ONLY: M2OBase_t
  IMPLICIT NONE

  TYPE attrRadar_t
    REAL(r_kind), ALLOCATABLE :: factor(:, :)
  END TYPE attrRadar_t

  TYPE, EXTENDS(M2OBase_t) :: OprRadarVel_t
    TYPE(ObsSet_t), POINTER :: Y
    TYPE(State_t), POINTER :: X
    REAL(r_kind) :: a1, a2
    TYPE(attrRadar_t), ALLOCATABLE :: attrRadars(:)

  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: transBackward

    PROCEDURE, PUBLIC, PASS(this) :: transFwdNonLinear_opr
    PROCEDURE, PUBLIC, PASS(this) :: transFwdTanLinear_opr
    PROCEDURE, PUBLIC, PASS(this) :: transAdjMultiply_opr

    PROCEDURE, PUBLIC, PASS(this) :: transFwdNonLinear
    PROCEDURE, PUBLIC, PASS(this) :: transFwdTanLinear
    PROCEDURE, PUBLIC, PASS(this) :: transAdjMultiply

    PROCEDURE :: fwdNL_opr => transFwdNonLinear_opr
    PROCEDURE :: fwdTL_opr => transFwdTanLinear_opr
    PROCEDURE :: adjMul_opr => transAdjMultiply_opr

    PROCEDURE :: fwdNL => transFwdNonLinear
    PROCEDURE :: fwdTL => transFwdTanLinear
    PROCEDURE :: adjMul => transAdjMultiply

    PROCEDURE, PUBLIC, NOPASS :: transElemForwardForDA
    PROCEDURE, PUBLIC, NOPASS :: transElemForwardForThinning

    PROCEDURE, PRIVATE, NOPASS :: transElemAdjointMultiply
    PROCEDURE, PUBLIC, NOPASS :: calRadarVelFactorElem
    PROCEDURE, PRIVATE :: calRadarVelFactor
    ! GENERIC :: OPERATOR(.yaml.) => transFwdNonLinear_opr
    ! GENERIC :: OPERATOR(.ADJ.) => transAdjMultiply_opr
    FINAL :: destructor
  END TYPE OprRadarVel_t

  INTERFACE OprRadarVel_t
    PROCEDURE :: constructor
  END INTERFACE OprRadarVel_t

CONTAINS

  FUNCTION constructor(configFile, X, Y) RESULT(this)
    IMPLICIT NONE
    TYPE(OprRadarVel_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X
    TYPE(ObsSet_t), OPTIONAL, TARGET, INTENT(IN) :: Y

    this%X => X
    IF (PRESENT(Y)) THEN
      this%Y => Y
    END IF

    CALL this%calRadarVelFactor()
  END FUNCTION constructor

  SUBROUTINE calRadarVelFactor(this)
    CLASS(OprRadarVel_t) :: this
    INTEGER(i_kind) :: i, j, k, l, m, n

    ALLOCATE (this%attrRadars(SIZE(this%Y%ObsFields)))
    DO i = 1, SIZE(this%Y%ObsFields)
      IF (TRIM(this%Y%ObsFields(i)%Get_Name()) .EQ. 'rwnd') THEN
        ASSOCIATE (rwnd => this%Y%ObsFields(i))
          ALLOCATE (this%attrRadars(i)%factor(3, SIZE(rwnd%values)))
          DO j = 1, SIZE(rwnd%values)
            this%attrRadars(i)%factor(:, j) = this%calRadarVelFactorElem( &
                                              rwnd%locObs, &
                                              (/this%X%sg%cell_cntr(1, rwnd%idx(j)%hIdx), this%X%sg%cell_cntr(2, rwnd%idx(j)%hIdx), this%X%sg%zHght(rwnd%idx(j)%vIdx, rwnd%idx(j)%hIdx)/))
          END DO
        END ASSOCIATE
      END IF
    END DO

    ! Multi radar QC
    BLOCK
      USE parameters_m, ONLY: degree2radian
      REAL(r_kind) theta_thres
      INTEGER(i_kind) :: temp(this%X%sg%vLevel, this%X%sg%num_icell, this%X%sg%tSlots)
      INTEGER(i_kind), ALLOCATABLE :: radarIndices(:, :, :)
      INTEGER(i_kind) :: maxCount, idx, countValid, l, m, n, countColine
      LOGICAL :: maskTemp(this%X%sg%vLevel, this%X%sg%num_icell, this%X%sg%tSlots)

      temp = 0
      theta_thres = DCOS(70.0 * degree2radian)

      ! 统计每个网格点的雷达观测数量
      DO l = 1, SIZE(this%Y%ObsFields)
        IF (TRIM(this%Y%ObsFields(l)%Get_Name()) .EQ. 'rwnd') THEN
          DO m = 1, SIZE(this%Y%ObsFields(l)%values)
            temp(this%Y%ObsFields(l)%idx(m)%vidx, this%Y%ObsFields(l)%idx(m)%hidx, this%Y%ObsFields(l)%idx(m)%tidx) = &
              temp(this%Y%ObsFields(l)%idx(m)%vidx, this%Y%ObsFields(l)%idx(m)%hidx, this%Y%ObsFields(l)%idx(m)%tidx) + 1
          END DO
        END IF
      END DO

      ! 找出最大观测数量以分配数组
      maxCount = MAXVAL(temp)
      IF (maxCount <= 1) RETURN ! 如果没有多重观测点，直接返回

      ! 创建掩码以识别多重观测点
      maskTemp = temp > 1

      ! 为每个多重观测点分配索引存储
      ALLOCATE (radarIndices(maxCount, 2, COUNT(maskTemp)))

      ! 遍历多重观测点
      idx = 0
      DO i = 1, this%X%sg%vLevel
        DO j = 1, this%X%sg%num_icell
          DO k = 1, this%X%sg%tSlots
            IF (.NOT. maskTemp(i, j, k)) CYCLE

            idx = idx + 1
            countValid = 0

            ! 收集该点的所有观测
            DO l = 1, SIZE(this%Y%ObsFields)
              IF (TRIM(this%Y%ObsFields(l)%Get_Name()) .NE. 'rwnd') CYCLE

              DO m = 1, SIZE(this%Y%ObsFields(l)%values)
                IF (this%Y%ObsFields(l)%idx(m)%vidx == i .AND. &
                    this%Y%ObsFields(l)%idx(m)%hidx == j .AND. &
                    this%Y%ObsFields(l)%idx(m)%tidx == k) THEN

                  countValid = countValid + 1
                  radarIndices(countValid, 1, idx) = l
                  radarIndices(countValid, 2, idx) = m

                  IF (countValid == temp(i, j, k)) EXIT
                END IF
              END DO

              IF (countValid == temp(i, j, k)) EXIT
            END DO

            ! 筛选相似方向的观测
            countColine = 1
            DO l = 2, countValid
              ! 检查是否与之前的观测相似
              DO m = 1, l - 1
                ASSOCIATE (vec1 => this%attrRadars(radarIndices(l, 1, idx))%factor(:, radarIndices(l, 2, idx)), &
                           vec2 => this%attrRadars(radarIndices(m, 1, idx))%factor(:, radarIndices(m, 2, idx)))
                  IF (vec2(1) .EQ. 0 .OR. vec1(1) .EQ. 0) CYCLE
                  IF (ABS(SUM(vec1 * vec2)) > theta_thres) THEN
                    vec1 = 0.0D0
                    this%Y%ObsFields(radarIndices(l, 1, idx))%values(radarIndices(l, 2, idx)) = 0.0D0
                    EXIT ! 一旦找到相似性，可以跳出内循环
                  END IF
                END ASSOCIATE
              END DO

              ! 如果方向不相似，增加计数
              IF (this%attrRadars(radarIndices(l, 1, idx))%factor(1, radarIndices(l, 2, idx)) .NE. 0) THEN
                countColine = countColine + 1
              END IF

              ! 超过3个不同方向的观测，移除其余的
              IF (countColine == 3) THEN
                DO n = l + 1, countValid
                  this%attrRadars(radarIndices(n, 1, idx))%factor(:, radarIndices(n, 2, idx)) = 0
                  this%Y%ObsFields(radarIndices(n, 1, idx))%values(radarIndices(n, 2, idx)) = 0.0D0
                END DO
                EXIT
              END IF
            END DO

          END DO
        END DO
      END DO

      DEALLOCATE (radarIndices)
    END BLOCK

  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(OprRadarVel_t), INTENT(INOUT) :: this
    INTEGER(i_kind) :: i

    IF (ALLOCATED(this%attrRadars)) THEN
      DO i = 1, SIZE(this%attrRadars)
        IF (ALLOCATED(this%attrRadars(i)%factor)) DEALLOCATE (this%attrRadars(i)%factor)
      END DO
      DEALLOCATE (this%attrRadars)
    END IF
  END SUBROUTINE destructor

  SUBROUTINE transBackward(this, X)
    IMPLICIT NONE
    CLASS(OprRadarVel_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X
    INTEGER(i_kind) :: i, j, k

  END SUBROUTINE

  SUBROUTINE transFwdNonLinear(this, X, Y)
    IMPLICIT NONE
    CLASS(OprRadarVel_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(ObsSet_t), INTENT(INOUT) :: Y

    INTEGER(i_kind) :: i, j, k
    LOGICAL :: hasWwnd = .TRUE.

    IF ((X%getVarIdx(TRIM('uwnd')) .EQ. 0) .OR. &
        (X%getVarIdx(TRIM('vwnd')) .EQ. 0) .OR. &
        (Y%getObsIdx(TRIM('rwnd')) .EQ. 0)) THEN
      PRINT *, 'Error use the OprRadarVel_t in transFwdNonLinear, STOP! '
      STOP
    END IF

    IF (X%getVarIdx(TRIM('wwnd')) .EQ. 0) THEN
      hasWwnd = .FALSE.
      CALL X%addVar('wwnd')
    END IF

    ASSOCIATE (uwnd => X%Fields(X%getVarIdx(TRIM('uwnd'))), &
               vwnd => X%Fields(X%getVarIdx(TRIM('vwnd'))), &
               wwnd => X%Fields(X%getVarIdx(TRIM('wwnd'))))
      DO i = 1, SIZE(Y%ObsFields)
        IF (TRIM(this%Y%ObsFields(i)%Get_Name()) .EQ. 'rwnd') THEN
          ASSOCIATE (rwnd => Y%ObsFields(i))
            DO k = LBOUND(rwnd%idx, 1), UBOUND(rwnd%idx, 1)
              rwnd%values(k) = transElemForwardForDA(uwnd%Get_Value(rwnd%idx(k)), &
                                                     vwnd%Get_Value(rwnd%idx(k)), &
                                                     wwnd%Get_Value(rwnd%idx(k)), &
                                                     this%attrRadars(i)%factor(:, k))
            END DO
          END ASSOCIATE
        END IF
      END DO
    END ASSOCIATE

    IF (.NOT. hasWwnd) CALL X%rmVar('wwnd')
  END SUBROUTINE

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(Y)
    IMPLICIT NONE
    CLASS(OprRadarVel_t) :: this
    TYPE(State_t), INTENT(in) :: X
    TYPE(ObsSet_t) :: Y

    Y = this%Y%zeroCopy()
    CALL this%transFwdNonLinear(X, Y)
  END FUNCTION

  SUBROUTINE transAdjMultiply(this, D, dX, X)
    IMPLICIT NONE
    CLASS(OprRadarVel_t) :: this
    TYPE(ObsSet_t), INTENT(IN) :: D
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    LOGICAL :: hasWwnd = .TRUE.

    INTEGER(i_kind) :: i, j, k

    IF ((dX%getVarIdx(TRIM('uwnd')) .EQ. 0) .OR. &
        (dX%getVarIdx(TRIM('vwnd')) .EQ. 0) .OR. &
        (D%getObsIdx(TRIM('rwnd')) .EQ. 0)) THEN
      PRINT *, 'Error use the OprRadarVel_t in transAdjMultiply, STOP! '
      STOP
    END IF

    CALL dX%clearHalo()
    ! Add a wwnd to dX if not have
    IF (dX%getVarIdx(TRIM('wwnd')) .EQ. 0) THEN
      hasWwnd = .FALSE.
      CALL dX%addVar('wwnd')
    END IF

    ASSOCIATE (uwnd => dX%Fields(dX%getVarIdx(TRIM('uwnd'))), &
               vwnd => dX%Fields(dX%getVarIdx(TRIM('vwnd'))), &
               wwnd => dX%Fields(dX%getVarIdx(TRIM('wwnd'))))

      DO i = 1, SIZE(D%ObsFields)
        IF (TRIM(D%ObsFields(i)%Get_Name()) .EQ. 'rwnd') THEN
          ASSOCIATE (rwnd => D%ObsFields(i))
            DO k = LBOUND(rwnd%idx, 1), UBOUND(rwnd%idx, 1)
              BLOCK
                REAL(r_kind) :: out(3)
                out = this%transElemAdjointMultiply(rwnd%values(k), this%attrRadars(i)%factor(:, k))
                ! PRINT *, 'out: ', out
                CALL uwnd%Add_Value(rwnd%idx(k), out(1))
                CALL vwnd%Add_Value(rwnd%idx(k), out(2))
                CALL wwnd%Add_Value(rwnd%idx(k), out(3))
                ! PRINT*, 'here'
                ! STOP
              END BLOCK
            END DO
          END ASSOCIATE
        END IF
      END DO
      wwnd%DATA(1, :, :) = 0.0D0
    END ASSOCIATE

    ! Remove wwnd from dX if not have
    IF (.NOT. hasWwnd) CALL dX%rmVar('wwnd')

    CALL dX%exHaloRevSum()
    CALL dX%exHalo()

  END SUBROUTINE transAdjMultiply

  FUNCTION transAdjMultiply_opr(this, D, X) RESULT(dX)
    IMPLICIT NONE
    CLASS(OprRadarVel_t) :: this
    TYPE(ObsSet_t), INTENT(IN) :: D
    TYPE(State_t) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    dX = X%zeroCopy()
    CALL this%transAdjMultiply(D, dX, X)
  END FUNCTION

  FUNCTION calRadarVelFactorElem(locRadar, locGrid) RESULT(factors)
    USE GeoTools_m, ONLY: distance
    USE parameters_m, ONLY: EarthRadius, degree2radian

    REAL(r_kind), INTENT(IN) :: &
      locRadar(3), & ! lat, lon, alt
      locGrid(3)
    REAL(r_kind) :: factors(3)  ! at u, v, w direction
    REAL(r_kind) :: AngXOE, lOV, l
    REAL(r_kind) :: AE

    ! AE = EarthRadius
    ! ! AE = 8500000.00D0

    ! AngXOE = distance(locRadar(1:2), locGrid(1:2))/AE
    ! lOV = (AE + locGrid(3))*DCOS(AngXOE)
    ! ! PRINT *, 'LOV: ', lOV, 'EarthRadius: ', EarthRadius, 'EarthRadius + locGrid(3): ', EarthRadius + locGrid(3)

    ! factors(3) = lOV - AE - locRadar(3)
    ! factors(1) = lOV*DTAN(locGrid(2) - locRadar(2))
    ! factors(2) = lOV*DTAN(locGrid(1) - locRadar(1))

    ! l = sum(factors*factors)**(0.5)
    ! ! PRINT *, 'vec: ', l, factors
    ! factors = factors/l
    REAL(r_kind) :: azimuth, slant_range, elev
    CALL latlon_to_radar(locGrid(1) / degree2radian, locGrid(2) / degree2radian, locGrid(3), &
                         azimuth, slant_range, elev, &
                         locRadar(1) / degree2radian, locRadar(2) / degree2radian, locRadar(3))

    factors(3) = sind(elev)
    factors(2) = cosd(elev) * cosd(azimuth)
    factors(1) = cosd(elev) * sind(azimuth)

  END FUNCTION

  SUBROUTINE latlon_to_radar(lat_grid, lon_grid, height_grid, azimuth, slant_range, elev, rlat_radar, rlon_radar, rheight_radar)
    REAL(r_kind) :: lat_grid, lon_grid, height_grid, azimuth, slant_range, elev, rlat_radar, rlon_radar, rheight_radar
    REAL(r_kind) :: rpd, mpd, radius_earth, radius_earth_8_thirds, diff, delta_lat, delta_lon
    REAL(r_kind) :: cos_factor, delta_y, delta_x, height_factor, hor_dist, delta_z, curvature, height_diff

    IF (rlat_radar .EQ. 0.0) THEN
      WRITE (6, *) ' Warning, Radar Coords NOT Initialized'
    END IF

    rpd = 3.141592653589 / 180.
    mpd = 111194.
    radius_earth = 6371.E3
    radius_earth_8_thirds = 6371.E3 * 2.6666666
    diff = 0.

    delta_lat = lat_grid - rlat_radar
    delta_lon = lon_grid - rlon_radar
    cos_factor = cosd(0.5 * (rlat_radar + lat_grid))
    delta_y = delta_lat * mpd
    delta_x = delta_lon * mpd * cos_factor
    height_factor = (radius_earth + 0.5 * (rheight_radar + height_grid)) / radius_earth
    hor_dist = SQRT(delta_x**2 + delta_y**2) * height_factor
    delta_z = height_grid - rheight_radar
    slant_range = SQRT(hor_dist**2 + delta_z**2)

    IF (hor_dist .GT. 0.0) THEN
      azimuth = ATAN2(delta_x, delta_y) / rpd
    ELSE
      azimuth = 0.0
    END IF

    curvature = hor_dist**2 / radius_earth_8_thirds
    height_diff = delta_z - curvature

    IF (slant_range .GT. 0.0) THEN
      elev = ATAN2(height_diff, hor_dist) / rpd
    ELSE
      elev = 0.0
    END IF
  END

  FUNCTION transElemForwardForDA(uwnd, vwnd, wwnd, factor) RESULT(rWnd)
    REAL(r_kind), INTENT(IN) :: uwnd, vwnd, wwnd, factor(3)
    REAL(r_kind) :: rWnd  ! at u, v, w direction

    rwnd = factor(1) * uwnd + factor(2) * vwnd + factor(3) * wwnd
  END FUNCTION

  FUNCTION transElemForwardForThinning(uwnd, vwnd, wwnd, locRadar, locGrid) RESULT(rWnd)
    REAL(r_kind), INTENT(IN) :: uwnd, vwnd, wwnd, &
                                locRadar(3), & ! lat, lon, alt
                                locGrid(3)
    REAL(r_kind):: rWnd
    REAL(r_kind) :: vec(3)  ! at u, v, w direction

    vec = calRadarVelFactorElem(locRadar, locGrid)
    rwnd = vec(1) * uwnd + vec(2) * vwnd + vec(3) * wwnd
  END FUNCTION

  FUNCTION transElemAdjointMultiply(in, factor) RESULT(out)
    REAL(r_kind) :: in, factor(3)
    REAL(r_kind) :: out(3)

    out = factor * in
  END FUNCTION

  SUBROUTINE transFwdTanLinear(this, dX, dY, X)
    CLASS(OprRadarVel_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(ObsSet_t), INTENT(INOUT) :: dY
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    CALL this%transFwdNonLinear(dX, dY)
  END SUBROUTINE

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dY)
    CLASS(OprRadarVel_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(ObsSet_t) :: dY

    dY = this%transFwdNonLinear_opr(dX)
  END FUNCTION

END MODULE OprRadarVel_m
