!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-Repository
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           : 2022-03-19   Created by Yuanfu Xie
!
!   Created by Yuanfu Xie (xieyf@gbamwf.com), 2022/03/19, @GBA-MWF, Shenzhen
!               Adapted from Zilong Qin's GeoBox_t for general grid and 4D.
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module provides a domain validation check, interpolation coefficients, indices in 4D domain
!! for a given a set of 4D points.

MODULE domainCheck_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE parameters_m
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE slint, ONLY: slint_init, tgt_grid
  USE RectangleGridUtility_m
  IMPLICIT NONE
  TYPE domainCheck_t
  CONTAINS
    PROCEDURE :: horizontalIntp
    PROCEDURE :: validation
  END TYPE domainCheck_t

  INTERFACE domainCheck_t
    PROCEDURE :: constructor
  END INTERFACE domainCheck_t

CONTAINS

  FUNCTION constructor() RESULT(this)
    TYPE(domainCheck_t) :: this
  END FUNCTION

  !> @brief
  !------------------------------------------------------------------
  !  This routine validates a given set of data along with it latlon,
  !  heights and times and
  SUBROUTINE horizontalIntp(this, sg, numPoints, latlon, &
                            idx, coe, interpolation, nstencil)

    IMPLICIT NONE

    CLASS(domainCheck_t) :: this
    TYPE(SingleGrid_t) :: sg
    INTEGER(i_kind), INTENT(IN) :: numPoints, interpolation
    REAL(r_kind), INTENT(IN) :: latlon(2, numPoints)
    INTEGER(i_kind), INTENT(OUT) :: nstencil
    INTEGER(i_kind), ALLOCATABLE, INTENT(OUT) :: idx(:, :)
    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: coe(:, :)

    ! Local variables:
    INTEGER(i_kind) :: i
    REAL(r_kind), ALLOCATABLE :: srcLL(:, :), tgtLL(:, :)

    ! Regular grid utility:
    TYPE(RectangleGridUtility_t) :: regular

    ! Allocate memory for interpolation:
    SELECT CASE (interpolation)
    CASE (1)
      nstencil = 3  ! Slint stencil
    CASE (2)
      nstencil = 4  ! Regular stencil
    CASE DEFAULT
      WRITE (*, 1) interpolation
1     FORMAT('No such interpolation scheme is implemented! check the namelist and rerun! ', I2)
      STOP
    END SELECT

    ! Memory for output arrays:
    ALLOCATE (idx(nstencil, numPoints))
    ALLOCATE (srcLL(sg%num_cell, 2), tgtLL(numPoints, 2), coe(nstencil, numPoints))
    srcLL(:, 1) = sg%cell_cntr(1, 1:sg%num_cell)
    srcLL(:, 2) = sg%cell_cntr(2, 1:sg%num_cell)
    PRINT *, 'srcLL is set', sg%mpddInfo_sg%myrank

    tgtLL(:, 1) = latlon(1, :)
    tgtLL(:, 2) = latlon(2, :)

    ! Interpolation:
    SELECT CASE (interpolation)
    CASE (1)
      ! Slint interpolation:
      CALL slint_init(srcLL, sg%num_cell, tgtLL, numPoints)
      WRITE (*, 123) numPoints, sg%mpddInfo_sg%myrank, sg%gLevel
123   FORMAT('Horizontal interpolating using slint - numPoints: ', I8, ' proc ', I1, ' G: ', I2)
      idx = tgt_grid%nn
      coe = tgt_grid%coeffs
      WRITE (*, 124) idx(:, 1), coe(:, 1), tgtll(1, :) / degree2radian, MINVAL(srcll(:, 1)) / degree2radian, &
        MAXVAL(srcll(:, 1)) / degree2radian, MINVAL(srcll(:, 2)) / degree2radian, MAXVAL(srcll(:, 2)) / degree2radian
124   FORMAT('IDX: ', 3I4, ' COEF: ', 3D12.4, ' tgtLL1:', 2D12.4, ' DomainLL: ', 4D12.4)
    CASE (2)
      CALL regular%RectangleHorizontalIntp(srcll, sg%num_cell, tgtll, numPoints, &
                                           idx, coe, sg%mpddInfo_sg%myrank)
    END SELECT
    WRITE (*, 125) interpolation
125 FORMAT('HorizontalIntp finished using interpolation scheme: ', I1)

  END SUBROUTINE horizontalIntp

  !> @brief
  !------------------------------------------------------------------
  !  This routine validates a given set of data along with it latlon,
  !  heights and times and
  SUBROUTINE validation(this, sg, &
                        numToValid, toValidLatlon, toValidVCoord, toValidTempor, &
                        numValided, maskValided, &
                        idxGrdValided, idxHgtValided, idxTimValided, &
                        coeGrdValided, coeHgtValided, coeTimValided, interpolation, nstencil, obsType, &
                        minLat_in, maxLat_in, minLon_in, maxLon_in)

    IMPLICIT NONE

    CLASS(domainCheck_t) :: this
    TYPE(SingleGrid_t) :: sg
    INTEGER(i_kind), INTENT(IN) :: numToValid, &
                                   toValidTempor(numToValid), interpolation
    REAL(r_kind), INTENT(IN) :: &
      toValidLatlon(2, numToValid), &
      toValidVCoord(numToValid)
    INTEGER(i_kind), INTENT(OUT) :: numValided, maskValided(numToValid), nstencil
    INTEGER(i_kind), ALLOCATABLE, INTENT(OUT) :: &
      idxGrdValided(:, :), idxHgtValided(:, :, :), idxTimValided(:, :) ! horizontal and height in one array
    REAL(r_kind), ALLOCATABLE, INTENT(OUT) :: &
      coeGrdValided(:, :), coeHgtValided(:, :, :), coeTimValided(:, :) ! horizontal and height in one array
    REAL(r_kind), INTENT(IN), OPTIONAL :: minLat_in, maxLat_in, minLon_in, maxLon_in ! if you want to specify a domain by yourself

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, preFilteredCount, maxCount
    INTEGER(i_kind), ALLOCATABLE :: idxValid(:), tempTempor(:), idxtgt(:, :)
    REAL(r_kind) :: dheight, dt
    ! Used to first round elimination of obs outside of domain on this processor
    REAL(r_kind) :: minLat, maxLat, minLon, maxLon
    REAL(r_kind), ALLOCATABLE :: tempVCoord(:), tempLatlon(:, :)
    REAL(r_kind), ALLOCATABLE :: srcLL(:, :), tgtLL(:, :), coetgt(:, :)

    REAL(r_kind) :: t1, t2, t_start
    CHARACTER(len=*) :: obsType
    LOGICAL :: isSfcObs = .FALSE.
    LOGICAL :: isSoundObs = .FALSE.
    LOGICAL, ALLOCATABLE :: timeValid(:), domainValid(:)
        ! 垂直域检查
        LOGICAL :: withinDom
    ! Regular grid utility:
    TYPE(RectangleGridUtility_t) :: regular

    CALL CPU_TIME(t_start) ! 记录总时间

    ! 设置观测类型标志 - 只做一次
    isSfcObs = (TRIM(obsType) == 'SYNOP')
    isSoundObs = (TRIM(obsType) == 'SOUND')

    ! 获取域的边界 - 只计算一次
    IF (PRESENT(minLat_in) .AND. PRESENT(maxLat_in) .AND. PRESENT(minLon_in) .AND. PRESENT(maxLon_in)) THEN
      minLat = minLat_in
      maxLat = maxLat_in
      minLon = minLon_in
      maxLon = maxLon_in
    ELSE
      minLat = MINVAL(sg%cell_cntr(1, :))
      maxLat = MAXVAL(sg%cell_cntr(1, :))
      minLon = MINVAL(sg%cell_cntr(2, :))
      maxLon = MAXVAL(sg%cell_cntr(2, :))
    END IF

    ! 设置插值点数 - 仅需设置一次
    SELECT CASE (interpolation)
    CASE (1)
      nstencil = 3  ! Slint stencil
    CASE (2)
      nstencil = 4  ! Regular stencil
    CASE DEFAULT
      WRITE (*, 1) interpolation
1     FORMAT('No such interpolation scheme is implemented! check the namelist and rerun! ', I2)
      STOP
    END SELECT

    ! 使用逻辑数组进行预过滤，比单独计数更高效
    ALLOCATE(timeValid(numToValid), domainValid(numToValid))
    
    ! 预过滤1: 时间窗口检查
    timeValid = (toValidTempor .GE. sg%tt(1)) .AND. (toValidTempor .LE. sg%tt(sg%tSlots))
    
    ! 预过滤2: 域边界和高度检查
    domainValid = timeValid .AND. &
                 (toValidLatlon(1, :) .GE. minLat) .AND. (toValidLatlon(1, :) .LE. maxLat) .AND. &
                 (toValidLatlon(2, :) .GE. minLon) .AND. (toValidLatlon(2, :) .LE. maxLon) .AND. &
                 ((sg%vLevel .EQ. 1) .OR. (toValidVCoord .LE. sg%zHght(sg%vLevel, 1)) .OR. isSoundObs)
    
    ! 计算有效点数并设置掩码
    maskValided = 0
    WHERE (domainValid) maskValided = [(i, i=1,numToValid)]
    preFilteredCount = COUNT(domainValid)
    
    ! 输出过滤统计信息
    WRITE (*, 100) sg%mpddInfo_sg%myrank, COUNT(.NOT. timeValid), numToValid, TRIM(obsType)
100 FORMAT('domainCheck - Data off time window at pc: ', I1, ' is : ', I6, ' of total: ', I10, ' ObsType:', A)
    WRITE (*, 300) sg%mpddInfo_sg%myrank, COUNT(timeValid .AND. .NOT. domainValid), numToValid, sg%zHght(UBOUND(sg%zHght, 1), 1), TRIM(obsType)
300 FORMAT('domainCheck - Data over bounds at pc: ', I1, ' is : ', I6, ' of total: ', I10, ' DomainTop', D12.4, ' ObsType:', A)

    ! 如果没有有效点，提前返回
    IF (preFilteredCount == 0) THEN
      numValided = 0
      ALLOCATE(idxGrdValided(nstencil, 0), idxHgtValided(2, nstencil, 0), idxTimValided(2, 0))
      ALLOCATE(coeGrdValided(nstencil, 0), coeHgtValided(2, nstencil, 0), coeTimValided(2, 0))
      DEALLOCATE(timeValid, domainValid)
      CALL CPU_TIME(t2)
      WRITE (*, 111) t2-t_start, sg%mpddInfo_sg%myrank, numToValid, 0, sg%num_cell
111   FORMAT('Time validation: ', D12.4, ' proc ', I1, ' NtoValid: ', I6, ' NValided: ', I6, ' Ncell: ', I8)
      RETURN
    END IF

    ! 为预过滤后的数据分配内存
    ALLOCATE (idxValid(preFilteredCount), idxtgt(nstencil, preFilteredCount))
    ALLOCATE (tempLatlon(2, preFilteredCount), tempTempor(preFilteredCount), tempVCoord(preFilteredCount))
    ALLOCATE (srcLL(sg%num_cell, 2), tgtLL(preFilteredCount, 2), coetgt(nstencil, preFilteredCount))

    ! 提取预过滤后的有效数据
    idxValid = PACK(maskValided, maskValided .NE. 0)
    tempVCoord = toValidVCoord(idxValid)
    tempTempor = toValidTempor(idxValid)
    tempLatlon(1:2, :) = toValidLatlon(1:2, idxValid)

    ! 设置插值源和目标点
    srcLL(:, 1) = sg%cell_cntr(1, 1:sg%num_cell)
    srcLL(:, 2) = sg%cell_cntr(2, 1:sg%num_cell)
    tgtLL(:, 1) = tempLatlon(1, :)
    tgtLL(:, 2) = tempLatlon(2, :)

    ! 执行水平插值
    SELECT CASE (interpolation)
    CASE (1)
      ! Slint插值
      CALL CPU_TIME(t1)
      CALL slint_init(srcLL, sg%num_cell, tgtLL, preFilteredCount)
      CALL CPU_TIME(t2)
      WRITE (*, 123) preFilteredCount, numToValid, sg%mpddInfo_sg%myrank, sg%gLevel, sg%vLevel, t2 - t1
123   FORMAT('Validating using slint - NValid: ', I8, ' NAll: ', I8, ' proc ', I1, ' G: ', I2, ' vLevels: ', I3, ' Slintime: ', D12.4)
      idxtgt = tgt_grid%nn
      coetgt = tgt_grid%coeffs
    CASE (2)
      CALL regular%RectangleHorizontalIntp(srcll, sg%num_cell, tgtll, preFilteredCount, &
                                         idxtgt, coetgt, sg%mpddInfo_sg%myrank)
    END SELECT

    CALL CPU_TIME(t1)

    ! 进一步过滤和计数有效点
    maxCount = preFilteredCount  ! 预分配最大可能数量
    
    ! 创建临时数组以存储最终有效数据的索引
    ALLOCATE(idxGrdValided(nstencil, maxCount), idxHgtValided(2, nstencil, maxCount), &
             idxTimValided(2, maxCount), coeGrdValided(nstencil, maxCount), &
             coeHgtValided(2, nstencil, maxCount), coeTimValided(2, maxCount))

    numValided = 0
    DO i = 1, preFilteredCount
      ! 检查是否在水平网格内
      IF (MAXVAL(idxtgt(:, i)) .LE. sg%num_cell .AND. &
          MINVAL(idxtgt(:, i)) .GT. 0 .AND. &
          COUNT(idxtgt(:, i) .EQ. MINVAL(idxtgt(:, i))) .EQ. 1 .AND. &
          COUNT(idxtgt(:, i) .EQ. MAXVAL(idxtgt(:, i))) .EQ. 1) THEN
          

        
        IF (isSoundObs) THEN
          withinDom = sg%vLevel .EQ. 1 .OR. &
                      (tempVCoord(i) .LE. MAXVAL(sg%pres(1, idxtgt(:, i))) .AND. &
                       tempVCoord(i) .GE. MINVAL(sg%pres(sg%vLevel, idxtgt(:, i))))
        ELSE
          withinDom = sg%vLevel .EQ. 1 .OR. &
                      (tempVCoord(i) .GE. MINVAL(sg%zHght(1, idxtgt(:, i))) .AND. &
                       tempVCoord(i) .LE. MAXVAL(sg%zHght(sg%vLevel, idxtgt(:, i)))) &
                      .OR. tempVCoord(i) .LE. -1.0E7 .OR. isSfcObs
        END IF
        
        IF (withinDom) THEN
          numValided = numValided + 1
          
          ! 保存水平插值信息
          idxGrdValided(:, numValided) = idxtgt(:, i)
          coeGrdValided(:, numValided) = coetgt(:, i)
          
          ! 计算和保存高度插值系数
          DO j = 1, nstencil
            idxHgtValided(1, j, numValided) = 1
            
            IF (sg%vLevel .EQ. 1 .OR. tempVCoord(i) .LE. -1.0E7 .OR. isSfcObs) THEN
              idxHgtValided(1:2, j, numValided) = 1
              coeHgtValided(1, j, numValided) = 1.0D0
              coeHgtValided(2, j, numValided) = 0.0D0
            ELSE
              IF (isSoundObs) THEN
                ! 查找包含该点的压力层
                DO k = 1, sg%vLevel
                  IF (tempVCoord(i) .GT. sg%pres(k, idxtgt(j, i))) THEN
                    idxHgtValided(1, j, numValided) = MAX(1, k - 1)
                    idxHgtValided(2, j, numValided) = k
                    
                    ! 计算压力插值系数
                    IF (k - 1 .LT. 1) THEN
                      coeHgtValided(1:2, j, numValided) = 0.5D0
                    ELSE
                      dheight = sg%pres(k, idxtgt(j, i)) - sg%pres(k - 1, idxtgt(j, i))
                      coeHgtValided(1, j, numValided) = (sg%pres(k, idxtgt(j, i)) - tempVCoord(i)) / dheight
                      coeHgtValided(2, j, numValided) = (tempVCoord(i) - sg%pres(k - 1, idxtgt(j, i))) / dheight
                    END IF
                    
                    EXIT
                  ELSE IF (tempVCoord(i) .GE. sg%pres(k, idxtgt(j, i)) .AND. k .EQ. sg%vLevel) THEN
                    idxHgtValided(1:2, j, numValided) = k
                    coeHgtValided(1:2, j, numValided) = 0.5D0
                  END IF
                END DO
              ELSE
                ! 查找包含该点的高度层
                DO k = 1, sg%vLevel
                  IF (tempVCoord(i) .LT. sg%zHght(k, idxtgt(j, i))) THEN
                    idxHgtValided(1, j, numValided) = MAX(1, k - 1)
                    idxHgtValided(2, j, numValided) = k
                    
                    ! 计算高度插值系数
                    IF (k - 1 .LT. 1) THEN
                      coeHgtValided(1:2, j, numValided) = 0.5D0
                    ELSE
                      dheight = sg%zHght(k, idxtgt(j, i)) - sg%zHght(k - 1, idxtgt(j, i))
                      coeHgtValided(1, j, numValided) = (sg%zHght(k, idxtgt(j, i)) - tempVCoord(i)) / dheight
                      coeHgtValided(2, j, numValided) = (tempVCoord(i) - sg%zHght(k - 1, idxtgt(j, i))) / dheight
                    END IF
                    
                    EXIT
                  ELSE IF (ABS(tempVCoord(i) - sg%zHght(k, idxtgt(j, i))) < 1E-6 .AND. k .EQ. sg%vLevel) THEN
                    idxHgtValided(1:2, j, numValided) = k
                    coeHgtValided(1:2, j, numValided) = 0.5D0
                  END IF
                END DO
                
                ! 检查错误
                IF (MINVAL(idxHgtValided(:, j, numValided)) .LE. 0) THEN
                  WRITE (*, 5) idxHgtValided(1:2, j, numValided), tempVCoord(i), &
                    MAXVAL(sg%zHght(sg%vLevel, idxtgt(:, i)))
5                 FORMAT('domainCheck error: ', 2I6, ' obsHght: ', D17.6, ' domainHght: ', D17.6)
                  PRINT *, 'For test', tempVCoord(i) .LE. sg%zHght(sg%vLevel, idxtgt(j, i)), isSfcObs, &
                           tempVCoord(i), sg%zHght(sg%vLevel, idxtgt(j, i))
                  STOP
                END IF
              END IF
            END IF
          END DO
          
          ! 计算时间插值系数 - 使用二分搜索来提高效率
          idxTimValided(1:2, numValided) = 1
          coeTimValided(1:2, numValided) = 0.0D0
          
          ! 简化的时间搜索 - 这里我们知道时间是递增排序的
          DO k = sg%tSlots, 1, -1
            IF (tempTempor(i) .GE. sg%tt(k)) THEN
              idxTimValided(1, numValided) = k
              idxTimValided(2, numValided) = MIN(k + 1, sg%tSlots)
              
              dt = sg%tt(idxTimValided(2, numValided)) - sg%tt(k)
              
              IF (dt .LE. machineEps) THEN
                coeTimValided(1, numValided) = 0.0D0
                coeTimValided(2, numValided) = 1.0D0
              ELSE
                coeTimValided(1, numValided) = (sg%tt(idxTimValided(2, numValided)) - tempTempor(i)) / dt
                coeTimValided(2, numValided) = (tempTempor(i) - sg%tt(k)) / dt
              END IF
              
              EXIT
            END IF
          END DO
          
          ! 更新掩码表示此点完全有效
          maskValided(idxValid(i)) = idxValid(i)
        ELSE
          maskValided(idxValid(i)) = 0
        END IF
      ELSE
        maskValided(idxValid(i)) = 0
      END IF
    END DO
    
    ! 如果没有有效点，释放所有临时数组
    IF (numValided == 0) THEN
      DEALLOCATE(idxGrdValided, idxHgtValided, idxTimValided, coeGrdValided, coeHgtValided, coeTimValided)
      ALLOCATE(idxGrdValided(nstencil, 0), idxHgtValided(2, nstencil, 0), idxTimValided(2, 0))
      ALLOCATE(coeGrdValided(nstencil, 0), coeHgtValided(2, nstencil, 0), coeTimValided(2, 0))
    ELSE IF (numValided < maxCount) THEN
      ! 重新分配为实际大小，避免浪费内存
      BLOCK
INTEGER(i_kind), ALLOCATABLE :: tempIdx1(:,:), tempIdx2(:,:,:), tempIdx3(:,:)
      REAL(r_kind), ALLOCATABLE :: tempCoe1(:,:), tempCoe2(:,:,:), tempCoe3(:,:)
      
      ! 保存当前数据
      ALLOCATE(tempIdx1(nstencil, numValided), tempIdx2(2, nstencil, numValided), tempIdx3(2, numValided))
      ALLOCATE(tempCoe1(nstencil, numValided), tempCoe2(2, nstencil, numValided), tempCoe3(2, numValided))
      
      tempIdx1 = idxGrdValided(:,1:numValided)
      tempIdx2 = idxHgtValided(:,:,1:numValided)
      tempIdx3 = idxTimValided(:,1:numValided)
      tempCoe1 = coeGrdValided(:,1:numValided)
      tempCoe2 = coeHgtValided(:,:,1:numValided) 
      tempCoe3 = coeTimValided(:,1:numValided)
      
      ! 重新分配到正确大小
      DEALLOCATE(idxGrdValided, idxHgtValided, idxTimValided)
      DEALLOCATE(coeGrdValided, coeHgtValided, coeTimValided)
      
      ALLOCATE(idxGrdValided(nstencil, numValided), idxHgtValided(2, nstencil, numValided), idxTimValided(2, numValided))
      ALLOCATE(coeGrdValided(nstencil, numValided), coeHgtValided(2, nstencil, numValided), coeTimValided(2, numValided))
      
      ! 复制数据
      idxGrdValided = tempIdx1
      idxHgtValided = tempIdx2
      idxTimValided = tempIdx3
      coeGrdValided = tempCoe1
      coeHgtValided = tempCoe2
      coeTimValided = tempCoe3
      
      ! 清理临时数组
      DEALLOCATE(tempIdx1, tempIdx2, tempIdx3, tempCoe1, tempCoe2, tempCoe3)
      END BLOCK
    END IF
    
    ! 释放临时数组
    DEALLOCATE(idxValid, idxtgt, tempLatlon, tempTempor, tempVCoord, srcLL, tgtLL, coetgt, timeValid, domainValid)
    
    CALL CPU_TIME(t2)
    WRITE (*, 110) t2-t_start, sg%mpddInfo_sg%myrank, numToValid, numValided, sg%num_cell
110 FORMAT('Time validation: ', D12.4, ' proc ', I1, ' NtoValid: ', I6, ' NValided: ', I6, ' Ncell: ', I8)
    
  END SUBROUTINE validation
END MODULE domainCheck_m
