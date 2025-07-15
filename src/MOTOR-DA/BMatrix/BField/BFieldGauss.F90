!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.BMatrix.BField
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/4/23, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/10/19, @GBA-MWF, Shenzhen
!     for adding an option currentField preparing to add a conversion of pressure obs to log(pressure).
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/04/18, @GBA-MWF, Shenzhen
!     Changing the Laplacian operators applied to control variables only, i.e., maskVertical,
!     maskHorizontal, maskTemporal equal to 1.
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE BFieldLaplace_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE MultiGrid_m, ONLY: MultiGrid_t
  USE AuxTypeSG_m, ONLY: AuxTypeSG_t
  USE BFieldBase_m, ONLY: BFieldBase_t
  USE ObsSet_m, ONLY: ObsSet_t
  IMPLICIT NONE

  TYPE, EXTENDS(BFieldBase_t):: BFieldLaplace_t
    TYPE(ObsSet_t), POINTER :: Y
    LOGICAL :: useNeumCondOnBndy = .TRUE.
    LOGICAL :: disableBmatForUV = .FALSE.
    REAL(r_kind) ::  c1 = 1.0D0, c2 = -2.0, c3 = 1.0D0

  CONTAINS
    PROCEDURE, PRIVATE :: loadBMatFiles

    PROCEDURE, PUBLIC :: sqrt_inverse_multiply
    PROCEDURE, PUBLIC :: sqrt_inverse_multiply_adjoint

    PROCEDURE :: sqrtInvMul => sqrt_inverse_multiply
    PROCEDURE :: sqrtInvMulAdj => sqrt_inverse_multiply_adjoint

    PROCEDURE, PUBLIC :: inverse_multiply

    PROCEDURE, PUBLIC :: Initialize
    FINAL :: destructor
  END TYPE BFieldLaplace_t

CONTAINS

  SUBROUTINE Initialize(this, configFile, sg, varName, Y)
    !USE NMLRead_m
    USE YAMLRead_m
    IMPLICIT NONE

    CLASS(BFieldLaplace_t) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg

    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    CHARACTER(LEN=10), INTENT(IN) :: varName
    REAL(r_kind) :: RelativeWeightJb2Jo

    INTEGER(i_kind) :: ifile, numCell, i
    TYPE(ObsSet_t), TARGET, OPTIONAL :: Y

    !this%AuxTypeSG_t = AuxTypeSG_t(sg)
    CALL this%AuxTypeSG_t%aux_initialize(sg)
    this%name = varName

    this%scaleParaX = 1.0D0
    this%scaleParaY = 1.0D0
    this%scaleParaZ = 1.0D0
    this%scaleParaT = 1.0D0

    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaX', this%scaleParaX)
    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaY', this%scaleParaY)
    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaZ', this%scaleParaZ)
    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaT', this%scaleParaT)

    IF (this%scaleParaX >= 500000) this%scaleParaX = this%scaleParaX / sg%dlon_dis / sg%dlon_dis
    IF (this%scaleParaY >= 500000) this%scaleParaY = this%scaleParaY / sg%dlat_dis / sg%dlat_dis
    IF (this%ScaleParaT >= 500000) THEN
      IF (sg%tSlots > 1) this%ScaleParaT = this%ScaleParaT / (this%sg%tt(2) - this%sg%tt(1)) / (this%sg%tt(2) - this%sg%tt(1))
    END IF
    CALL this%loadBMatFiles(varName)

    PRINT *, 'loadBMatFiles of ', TRIM(varName), ' min Sigma:', MINVAL(this%sigma)
    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'RelativeWeightJb2Jo', RelativeWeightJb2Jo)
    IF (ifile /= 0) THEN
      PRINT *, 'RelativeWeightJb2Jo is not found in the yaml file, set to 200.0'
      RelativeWeightJb2Jo = 200.0
    END IF

    ! Normlize the parameters
    RelativeWeightJb2Jo = 1.0 / (this%scaleParaX + this%scaleParaY + this%scaleParaZ + this%scaleParaT) * RelativeWeightJb2Jo
    IF(this%scaleParaX + this%scaleParaY + this%scaleParaZ + this%scaleParaT == 0) RelativeWeightJb2Jo = 0.0


    ! Disable the Bmat for uwnd and vwnd
    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'disableBmatForUV', this%disableBmatForUV)
    IF (ifile /= 0) THEN
      PRINT *, 'disableBmatForUV is not found in the yaml file, set to .FALSE.'
      this%disableBmatForUV = .FALSE.
    END IF

    IF (this%disableBmatForUV) THEN
      IF (TRIM(varName) == 'uwnd') THEN
        RelativeWeightJb2Jo = 0.0
      END IF
      IF (TRIM(varName) == 'vwnd') THEN
        RelativeWeightJb2Jo = 0.0
      END IF
    END IF
    PRINT *, 'BFieldLaplace_t: ', TRIM(varName), ' this%scaleParaX, this%scaleParaY, this%scaleParaZ, this%scaleParaT', this%scaleParaX, this%scaleParaY, this%scaleParaZ, this%scaleParaT

    this%scaleParaX = this%scaleParaX * RelativeWeightJb2Jo
    this%scaleParaY = this%scaleParaY * RelativeWeightJb2Jo
    this%scaleParaZ = this%scaleParaZ * RelativeWeightJb2Jo
    this%scaleParaT = this%scaleParaT * RelativeWeightJb2Jo

    ! BLOCK
    !   REAL (r_kind) :: mean_x, mean_y, totalWeight
    !   CALL sg%mpddInfo_sg%AllReduceSumReal(SUM(sg%cell_dist(4, 1:sg%num_icell)), mean_x)
    !   mean_x = mean_x / sg%dimCell_global(1)/sg%dimCell_global(2)
    !   CALL sg%mpddInfo_sg%AllReduceSumReal(SUM(sg%cell_dist(2, 1:sg%num_icell)), mean_y)
    !   mean_y = mean_y / sg%dimCell_global(1)/sg%dimCell_global(2)

    !   totalWeight = this%scaleParaX+this%scaleParaY
    !   mean_x =  1/mean_x/mean_x
    !   mean_y =  1/mean_y/mean_y

    !   this%scaleParaX = (mean_x)/(mean_x + mean_y)*totalWeight
    !   this%scaleParaY = (mean_y)/(mean_x + mean_y)*totalWeight
    !   PRINT *, 'mean_x, mean_y, totalWeight, this%scaleParaX, this%scaleParaY', mean_x, mean_y, totalWeight, this%scaleParaX, this%scaleParaY
    !   ! STOP
    ! END BLOCK

    SELECT TYPE (mg => sg%mgParent)
    TYPE is (MultiGrid_t)
      PRINT *, 'MultiGrid_t', mg%mg_coarsest, mg%mg_finest
      IF (sg%gLevel < mg%mg_finest) THEN
        DO i = mg%mg_finest, sg%gLevel + 1, -1
          this%scaleParaX = this%scaleParaX / 4.0
          this%scaleParaY = this%scaleParaY / 4.0
          IF (sg%gLevel > 1) this%scaleParaZ = this%scaleParaZ / (((mg%sg(i)%vLevel - 1) / (mg%sg(i - 1)%vLevel - 1))**2.0D0)
          this%scaleParaT = this%scaleParaT / (((mg%sg(i)%tSlots - 1) / (mg%sg(i - 1)%tSlots - 1))**2.0D0)
        END DO
      END IF
    END SELECT

    ! Yuanfu Xie modified the print statement for more info: : 2023/10/20
    WRITE (*, 1) sg%dlon_dis, sg%dlat_dis, sg%gLevel, sg%mpddInfo_sg%myrank
1   FORMAT('Laplace B constructor- gridSizes in X & Y:', 2E12.4, ' Glvl: ', I2, ' pc:', I3)

    ! ---------------------------------------- For testing ----------------------------------------
    ! IF (TRIM(varname) == 'wwnd') this%scaleParaX = 0
    ! IF (TRIM(varname) == 'wwnd') this%scaleParaY = 0
    ! IF (TRIM(varname) == 'wwnd') this%scaleParaZ = 0
    ! IF (TRIM(varname) == 'wwnd') this%scaleParaT = 0
    ! ---------------------------------------- For testing ----------------------------------------

    IF (PRESENT(Y)) THEN
      this%Y => Y
    END IF

  END SUBROUTINE

  SUBROUTINE loadBMatFiles(this, varName)
    CLASS(BFieldLaplace_t) :: this
    CHARACTER(LEN=10), INTENT(IN) :: varName
    INTEGER(i_kind) :: i, k

    ! Here is the code of monk for simulating the fileinput of B mat files
    ALLOCATE (this%sigma(this%sg%vLevel, this%sg%num_cell, this%sg%tSlots))
    ! Set all the sigma to 1.0D0
    this%sigma = 1.0D0
    SELECT CASE (TRIM(varName))
    CASE ('pres')
      this%sigma = 100.0D0
    CASE ('psl')
      this%sigma = 100.0D0
    CASE ('qvapor')
      ! this%sigma = 0.001D0
      DO k = 1, this%sg%tSlots
        this%sigma(:, :, k) = this%sg%SVapor
      END DO
    CASE ('temp')
      DO k = 1, this%sg%tSlots
        this%sigma(:, :, k) = this%sg%STemp
      END DO
    END SELECT

  END SUBROUTINE loadBMatFiles

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(BFieldLaplace_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%sigma)) DEALLOCATE (this%sigma)
  END SUBROUTINE destructor

  SUBROUTINE inverse_multiply(this, field)
    IMPLICIT NONE
    CLASS(BFieldLaplace_t) :: this
    TYPE(Field_t), INTENT(INOUT) :: field

  END SUBROUTINE inverse_multiply

  ! Adding an option currentField preparing to add a conversion of pressure obs to log(pressure) obs
  SUBROUTINE sqrt_inverse_multiply_adjoint(this, field, currentField)
    IMPLICIT NONE
    CLASS(BFieldLaplace_t) :: this
    TYPE(Field_t), INTENT(INOUT) :: field
    TYPE(Field_t), OPTIONAL, INTENT(IN) :: currentField
    INTEGER :: i, ti, ni, k, it, j, t, kk, jj, tt
    REAL(r_kind), ALLOCATABLE :: swap(:, :, :)
    REAL(r_kind) :: h1, h2, h3, h4, a1, a2, a3, a4, as(4, this%sg%vLevel)

    ! IF(this%mpddSub%isBaseProc()) PRINT *, 'this%scalePara', this%scaleParaX, this%scaleParaY, this%scaleParaZ

    ALLOCATE (swap(this%sg%vLevel, this%sg%num_cell, this%sg%tSlots))
    ASSOCIATE (sg => this%sg)

      ! 1/sigma*X
      swap = field%DATA / this%sigma
      CALL sg%ExchangeMatOnHaloForFieldGrid(sg%tSlots, sg%vLevel, swap)

      field%DATA = 0.0D0
            BLOCK
        REAL(r_kind) :: gaussian(3,9,3) ! vLevel x num_cell x tSlots  
        REAL(r_kind) :: sigma_x = 300e3, sigma_y = 300e3, sigma_h = 2e3, sigma_t = 600, dx, dy
        sigma_x = sigma_x * sigma_x
        sigma_y = sigma_y * sigma_y
        sigma_h = sigma_h * sigma_h
        sigma_t = sigma_t * sigma_t

      DO k = 1, sg%vLevel
        DO j = 1, sg%num_icell
          DO t = 1, sg%tSlots
            
            gaussian = 1.0D0

            IF(k > 1 .AND. k < sg%vLevel) THEN
              gaussian(1, :, :) = gaussian(1, :, :) * exp((sg%sigma(k) - sg%sigma(k-1))**2/ sigma_h)
              gaussian(2, :, :) = gaussian(2, :, :) * exp((sg%sigma(k) - sg%sigma(k))**2/ sigma_h)
              gaussian(3, :, :) = gaussian(3, :, :) * exp((sg%sigma(k) - sg%sigma(k+1))**2/ sigma_h)
            ELSE IF(k == 1)THEN
              gaussian(1, :, :) = gaussian(1, :, :) * exp((sg%sigma(k) - sg%sigma(k+1))**2/ sigma_h)
              gaussian(2, :, :) = gaussian(2, :, :) * exp((sg%sigma(k) - sg%sigma(k))**2/ sigma_h)
              gaussian(3, :, :) = gaussian(3, :, :) * exp((sg%sigma(k) - sg%sigma(k+1))**2/ sigma_h)
            ELSE IF(k == sg%vLevel)THEN
              gaussian(1, :, :) = gaussian(1, :, :) * exp((sg%sigma(k) - sg%sigma(k-1))**2/ sigma_h)
              gaussian(2, :, :) = gaussian(2, :, :) * exp((sg%sigma(k) - sg%sigma(k))**2/ sigma_h)
              gaussian(3, :, :) = gaussian(3, :, :) * exp((sg%sigma(k) - sg%sigma(k-1))**2/ sigma_h)
            END IF
            
            IF(t > 1 .AND. t < sg%tSlots) THEN
              gaussian(:, :, 1) = gaussian(:, :, 1) * exp((sg%tt(t) - sg%tt(t-1))**2/ sigma_t)
              gaussian(:, :, 2) = gaussian(:, :, 2) * exp((sg%tt(t) - sg%tt(t))**2/ sigma_t)
              gaussian(:, :, 3) = gaussian(:, :, 3) * exp((sg%tt(t) - sg%tt(t+1))**2/ sigma_t)
            ELSE IF(t == 1)THEN
              gaussian(:, :, 1) = gaussian(:, :, 1) * exp((sg%tt(t) - sg%tt(t+1))**2/ sigma_t)
              gaussian(:, :, 2) = gaussian(:, :, 2) * exp((sg%tt(t) - sg%tt(t))**2/ sigma_t)
              gaussian(:, :, 3) = gaussian(:, :, 3) * exp((sg%tt(t) - sg%tt(t+1))**2/ sigma_t)
            ELSE IF(t == sg%tSlots)THEN
              gaussian(:, :, 1) = gaussian(:, :, 1) * exp((sg%tt(t) - sg%tt(t-1))**2/ sigma_t)
              gaussian(:, :, 2) = gaussian(:, :, 2) * exp((sg%tt(t) - sg%tt(t))**2/ sigma_t)
              gaussian(:, :, 3) = gaussian(:, :, 3) * exp((sg%tt(t) - sg%tt(t-1))**2/ sigma_t)
            END IF

            IF(sg%cell_dist(2, j) == 0) THEN 
              dy = sg%cell_dist(8, j)
            ELSE
              dy = sg%cell_dist(2, j)
            ENDIF

            IF(sg%cell_dist(4, j) == 0) THEN 
              dx = sg%cell_dist(6, j)
            ELSE
              dx = sg%cell_dist(4, j)
            ENDIF

            gaussian(:, 1, :) = gaussian(:, 1, :) * exp( dx ** 2 / sigma_x + dy **2/sigma_y)
            gaussian(:, 2, :) = gaussian(:, 2, :) * exp( dy ** 2 / sigma_y)
            gaussian(:, 3, :) = gaussian(:, 3, :) * exp( dx ** 2 / sigma_x + dy **2/sigma_y)
            gaussian(:, 4, :) = gaussian(:, 4, :) * exp( dx ** 2 / sigma_x)
            gaussian(:, 6, :) = gaussian(:, 6, :) * exp( dx ** 2 / sigma_x)
            gaussian(:, 7, :) = gaussian(:, 7, :) * exp( dx ** 2 / sigma_x + dy **2/sigma_y)
            gaussian(:, 8, :) = gaussian(:, 8, :) * exp( dy ** 2 / sigma_y)
            gaussian(:, 9, :) = gaussian(:, 9, :) * exp( dx ** 2 / sigma_x + dy **2/sigma_y)

            gaussian = gaussian/SUM(gaussian)
            gaussian (2, 5, 2) = gaussian(2, 5, 2) - 1.0

            DO kk = k-1, k+1
              DO jj = 1, 9
                DO tt = t-1, t+1
                  IF (kk > 0 .AND. kk <= sg%vLevel .AND. sg%cell_stcl(jj, j) /= 0 .AND. tt > 0 .AND. tt <= sg%tSlots) THEN
                    field%DATA(kk, sg%cell_stcl(jj, j), tt) = field%DATA(kk, sg%cell_stcl(jj, j), tt) + swap(k, j, t) * gaussian(kk-k+2, jj, tt-t+2)
                  END IF
                END DO
              END DO
            END DO

          END DO
        END DO
      END DO
      END BLOCK
     
      CALL sg%ExchangeMatOnHaloReverseSumForFieldGrid(sg%tSlots, sg%vLevel, field%DATA)! Exchange hale
      CALL sg%ExchangeMatOnHaloForFieldGrid(sg%tSlots, sg%vLevel, field%DATA)! Exchange hale

    END ASSOCIATE
    DEALLOCATE (swap)

  END SUBROUTINE

  SUBROUTINE sqrt_inverse_multiply(this, field)
    IMPLICIT NONE
    CLASS(BFieldLaplace_t) :: this
    TYPE(Field_t), INTENT(INOUT) :: field
    INTEGER :: i, ti, ni, k, j, t, kk, jj, tt
    REAL(r_kind), ALLOCATABLE :: swap(:, :, :)
    REAL(r_kind) :: h1, h2, h3, h4, a1, a2, a3, a4, as(4, this%sg%vLevel)
    REAL(r_kind) ::  dt

    ! IF(this%mpddSub%isBaseProc()) PRINT *, 'this%scalePara', this%scaleParaX, this%scaleParaY, this%scaleParaZ

    ! DEALLOCATE(this%sigma)
    !           ALLOCATE(this%sigma(this%sg%vLevel, this%sg%num_cell, this%sg%tSlots))
    !           this%sigma = 1.0

    ! PRINT*, 'sigma is ', this%sigma(1, 1, 1)
    ALLOCATE (swap(this%sg%vLevel, this%sg%num_cell, this%sg%tSlots))
    ASSOCIATE (sg => this%sg)

      ! 1/sigma*X
      swap = field%DATA / this%sigma
      CALL sg%ExchangeMatOnHaloForFieldGrid(sg%tSlots, sg%vLevel, swap)

      BLOCK
        REAL(r_kind) :: gaussian(3,9,3) ! vLevel x num_cell x tSlots  
        REAL(r_kind) :: sigma_x = 300e3, sigma_y = 300e3, sigma_h = 2e3, sigma_t = 600, dx, dy
        sigma_x = sigma_x * sigma_x
        sigma_y = sigma_y * sigma_y
        sigma_h = sigma_h * sigma_h
        sigma_t = sigma_t * sigma_t

      field%DATA = 0.0D0
      DO k = 1, sg%vLevel
        DO j = 1, sg%num_icell
          DO t = 1, sg%tSlots
            
            gaussian = 1.0D0

            IF(k > 1 .AND. k < sg%vLevel) THEN
              gaussian(1, :, :) = gaussian(1, :, :) * exp((sg%sigma(k) - sg%sigma(k-1))**2/ sigma_h)
              gaussian(2, :, :) = gaussian(2, :, :) * exp((sg%sigma(k) - sg%sigma(k))**2/ sigma_h)
              gaussian(3, :, :) = gaussian(3, :, :) * exp((sg%sigma(k) - sg%sigma(k+1))**2/ sigma_h)
            ELSE IF(k == 1)THEN
              gaussian(1, :, :) = gaussian(1, :, :) * exp((sg%sigma(k) - sg%sigma(k+1))**2/ sigma_h)
              gaussian(2, :, :) = gaussian(2, :, :) * exp((sg%sigma(k) - sg%sigma(k))**2/ sigma_h)
              gaussian(3, :, :) = gaussian(3, :, :) * exp((sg%sigma(k) - sg%sigma(k+1))**2/ sigma_h)
            ELSE IF(k == sg%vLevel)THEN
              gaussian(1, :, :) = gaussian(1, :, :) * exp((sg%sigma(k) - sg%sigma(k-1))**2/ sigma_h)
              gaussian(2, :, :) = gaussian(2, :, :) * exp((sg%sigma(k) - sg%sigma(k))**2/ sigma_h)
              gaussian(3, :, :) = gaussian(3, :, :) * exp((sg%sigma(k) - sg%sigma(k-1))**2/ sigma_h)
            END IF
            
            IF(t > 1 .AND. t < sg%tSlots) THEN
              gaussian(:, :, 1) = gaussian(:, :, 1) * exp((sg%tt(t) - sg%tt(t-1))**2/ sigma_t)
              gaussian(:, :, 2) = gaussian(:, :, 2) * exp((sg%tt(t) - sg%tt(t))**2/ sigma_t)
              gaussian(:, :, 3) = gaussian(:, :, 3) * exp((sg%tt(t) - sg%tt(t+1))**2/ sigma_t)
            ELSE IF(t == 1)THEN
              gaussian(:, :, 1) = gaussian(:, :, 1) * exp((sg%tt(t) - sg%tt(t+1))**2/ sigma_t)
              gaussian(:, :, 2) = gaussian(:, :, 2) * exp((sg%tt(t) - sg%tt(t))**2/ sigma_t)
              gaussian(:, :, 3) = gaussian(:, :, 3) * exp((sg%tt(t) - sg%tt(t+1))**2/ sigma_t)
            ELSE IF(t == sg%tSlots)THEN
              gaussian(:, :, 1) = gaussian(:, :, 1) * exp((sg%tt(t) - sg%tt(t-1))**2/ sigma_t)
              gaussian(:, :, 2) = gaussian(:, :, 2) * exp((sg%tt(t) - sg%tt(t))**2/ sigma_t)
              gaussian(:, :, 3) = gaussian(:, :, 3) * exp((sg%tt(t) - sg%tt(t-1))**2/ sigma_t)
            END IF

            IF(sg%cell_dist(2, j) == 0) THEN 
              dy = sg%cell_dist(8, j)
            ELSE
              dy = sg%cell_dist(2, j)
            ENDIF

            IF(sg%cell_dist(4, j) == 0) THEN 
              dx = sg%cell_dist(6, j)
            ELSE
              dx = sg%cell_dist(4, j)
            ENDIF

            gaussian(:, 1, :) = gaussian(:, 1, :) * exp( dx ** 2 / sigma_x + dy **2/sigma_y)
            gaussian(:, 2, :) = gaussian(:, 2, :) * exp( dy ** 2 / sigma_y)
            gaussian(:, 3, :) = gaussian(:, 3, :) * exp( dx ** 2 / sigma_x + dy **2/sigma_y)
            gaussian(:, 4, :) = gaussian(:, 4, :) * exp( dx ** 2 / sigma_x)
            gaussian(:, 6, :) = gaussian(:, 6, :) * exp( dx ** 2 / sigma_x)
            gaussian(:, 7, :) = gaussian(:, 7, :) * exp( dx ** 2 / sigma_x + dy **2/sigma_y)
            gaussian(:, 8, :) = gaussian(:, 8, :) * exp( dy ** 2 / sigma_y)
            gaussian(:, 9, :) = gaussian(:, 9, :) * exp( dx ** 2 / sigma_x + dy **2/sigma_y)

            gaussian = gaussian/SUM(gaussian)
            gaussian (2, 5, 2) = gaussian(2, 5, 2) - 1.0

            DO kk = k-1, k+1
              DO jj = 1, 9
                DO tt = t-1, t+1
                  IF (kk > 0 .AND. kk <= sg%vLevel .AND. sg%cell_stcl(jj, j) /= 0 .AND. tt > 0 .AND. tt <= sg%tSlots) THEN
                    field%DATA(k, j, t) = field%DATA(k, j, t) + swap(kk, sg%cell_stcl(jj, j), tt) * gaussian(kk-k+2, jj, tt-t+2)
                  END IF
                END DO
              END DO
            END DO

          END DO
        END DO
      END DO
      END BLOCK

      CALL sg%ExchangeMatOnHaloForFieldGrid(sg%tSlots, sg%vLevel, field%DATA)! Exchange hale

    END ASSOCIATE
    DEALLOCATE (swap)

  END SUBROUTINE sqrt_inverse_multiply

END MODULE BFieldLaplace_m
