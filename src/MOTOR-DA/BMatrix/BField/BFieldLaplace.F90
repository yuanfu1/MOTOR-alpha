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
    REAL(r_kind) :: RelativeWeightJb2Jo, temp

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

    ! Yuanfu Xie modified the yaml file for scalePara to reduce 0.1 instead of coding them here: : 2023/10/20
    ! Yuanfu Xie changed this to factor in the grid XY ratio (see the docx:
    ! How to Scale the Laplacian Operator for Background Covariance under $HOME/developments/da/motor_training/test/xie)
    temp = this%scaleParaX + this%scaleParaY
    this%scaleParaY = this%scaleParaY * (sg%dlon / sg%dlat)**2
    this%scaleParaX = this%scaleParaX / (this%scaleParaX + this%scaleParaY) * temp
    this%scaleParaY = temp - this%scaleParaX

    ! PRINT*, this%scaleParaY, this%scaleParaX
    ! STOP

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
    IF (this%scaleParaX + this%scaleParaY + this%scaleParaZ + this%scaleParaT == 0) RelativeWeightJb2Jo = 0.0

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

    ! IF (TRIM(varName) == 'uwnd') RelativeWeightJb2Jo = RelativeWeightJb2Jo / 5
    ! IF (TRIM(varName) == 'vwnd') RelativeWeightJb2Jo = RelativeWeightJb2Jo / 5


    this%scaleParaX = this%scaleParaX * RelativeWeightJb2Jo
    this%scaleParaY = this%scaleParaY * RelativeWeightJb2Jo
    this%scaleParaZ = this%scaleParaZ * RelativeWeightJb2Jo
    this%scaleParaT = this%scaleParaT * RelativeWeightJb2Jo
    PRINT *, 'BFieldLaplace_t: ', TRIM(varName), RelativeWeightJb2Jo, ' this%scaleParaX, this%scaleParaY, this%scaleParaZ, this%scaleParaT', this%scaleParaX, this%scaleParaY, this%scaleParaZ, this%scaleParaT

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
      IF (sg%gLevel < mg%mg_finest) THEN
        DO i = mg%mg_finest, sg%gLevel + 1, -1
          this%scaleParaX = this%scaleParaX / 4.0
          this%scaleParaY = this%scaleParaY / 4.0
          ! Yuanfu Xie added this check to ensure the vertical scaling is set properly on 2025-02-13:
          IF (mg%sg(i)%vLevel .GT. 1) THEN
            IF (sg%gLevel > 1) this%scaleParaZ = this%scaleParaZ / &
              (((mg%sg(i)%vLevel - 1) / (mg%sg(i - 1)%vLevel - 1))**2.0D0)
          ELSE
            this%scaleParaZ = this%scaleParaZ ! Temporarily set to default by Yuanfu Xie on 2025-02-13
          END IF
          this%scaleParaT = this%scaleParaT / (((mg%sg(i)%tSlots - 1) / (mg%sg(i - 1)%tSlots - 1))**2.0D0)
        END DO
      END IF
    END SELECT

    ! Yuanfu Xie modified the print statement for more info: : 2023/10/20
    WRITE (*, 1) sg%dlon_dis, sg%dlat_dis, sg%gLevel, sg%mpddInfo_sg%myrank
1   FORMAT('Laplace B constructor- gridSizes in X & Y:', 2E12.4, ' Glvl: ', I2, ' pc:', I3)

    ! ---------------------------------------- For testing ----------------------------------------
    IF (TRIM(varname) == 'wwnd') this%scaleParaX = 0
    IF (TRIM(varname) == 'wwnd') this%scaleParaY = 0
    IF (TRIM(varname) == 'wwnd') this%scaleParaZ = 0
    IF (TRIM(varname) == 'wwnd') this%scaleParaT = 0
    ! ---------------------------------------- For testing ----------------------------------------

    IF (PRESENT(Y)) THEN
      this%Y => Y
    END IF

  END SUBROUTINE Initialize

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
    INTEGER :: i, ti, ni, k, it
    REAL(r_kind), ALLOCATABLE :: swap(:, :, :)
    REAL(r_kind) :: h1, h2, h3, h4, a1, a2, a3, a4, as(4, this%sg%vLevel)

    ! IF(this%mpddSub%isBaseProc()) PRINT *, 'this%scalePara', this%scaleParaX, this%scaleParaY, this%scaleParaZ

    ALLOCATE (swap(this%sg%vLevel, this%sg%num_cell, this%sg%tSlots))
    ASSOCIATE (sg => this%sg)

      ! 1/sigma*X
      swap = field%DATA / this%sigma
      CALL sg%ExchangeMatOnHaloForFieldGrid(sg%tSlots, sg%vLevel, swap)

      field%DATA = 0.0D0
      ! Laplacian x direction
      IF (field%maskHorizontal .EQ. 1) THEN
        ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
        DO i = 1, this%sg%num_icell
          !        IF (sg%cell_type(i) .EQ. 0) THEN
          IF ((sg%cell_stcl(4, i) .NE. 0) .AND. (sg%cell_stcl(6, i) .NE. 0)) THEN

            ! field%data(:, i, :) = field%data(:, i, :) + ((swap(:, sg%cell_stcl(4, i), :) + swap(:, sg%cell_stcl(6, i), :) - &
            !                                               2*swap(:, i, :))/(sg%cell_dist(4, i)**2))*this%ScaleParaX
            ! Yuanfu Xie added the horizonSimilarity parameters: 2023/10/20
            ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
            DO it = 1, sg%tSlots
              IF (field%maskTemporal(it) .EQ. 1) THEN
                DO k = 1, sg%vLevel
                  IF (field%maskVertical(k) .EQ. 1) THEN
                    field%DATA(k, sg%cell_stcl(4, i), it) = field%DATA(k, sg%cell_stcl(4, i), it) + &
                                                            swap(k, i, it) * this%scaleParaX * &
                                                            (sg%horizonSimilarity(i) * (1.0D0 - sg%vertcalSimilarity(k)) + &
                                                             sg%vertcalSimilarity(k)) * this%c1
                    field%DATA(k, sg%cell_stcl(6, i), it) = field%DATA(k, sg%cell_stcl(6, i), it) + &
                                                            swap(k, i, it) * this%scaleParaX * &
                                                            (sg%horizonSimilarity(i) * (1.0D0 - sg%vertcalSimilarity(k)) + &
                                                             sg%vertcalSimilarity(k)) * this%c3
                    field%DATA(k, i, it) = field%DATA(k, i, it) + swap(k, i, it) * this%scaleParaX * &
                                           (sg%horizonSimilarity(i) * (1.0D0 - sg%vertcalSimilarity(k)) + &
                                            sg%vertcalSimilarity(k)) * this%c2
                  END IF
                END DO
              END IF
            END DO
          END IF

          IF (this%useNeumCondOnBndy) THEN
            IF ((sg%cell_stcl(4, i) .EQ. 0)) THEN
              ! Yuanfu Xie added the horizonSimilarity parameters: 2023/10/20
              ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
              DO it = 1, sg%tSlots
                IF (field%maskTemporal(it) .EQ. 1) THEN
                  DO k = 1, sg%vLevel
                    IF (field%maskVertical(k) .EQ. 1) THEN
                      ! field%data(:, i, :) = field%data(:, i, :) + ((swap(:, sg%cell_stcl(6, i), :) - swap(:, i, :))/(sg%cell_dist(6, i)))*(this%ScaleParaX)/(sg%cell_dist(6, i))
                      field%DATA(k, sg%cell_stcl(6, i), it) = field%DATA(k, sg%cell_stcl(6, i), it) &
                                                              + swap(k, i, it) * this%scaleParaX * &
                                                              (sg%horizonSimilarity(i) * (1.0D0 - sg%vertcalSimilarity(k)) + &
                                                               sg%vertcalSimilarity(k))
                      field%DATA(k, i, it) = field%DATA(k, i, it) &
                                             - swap(k, i, it) * this%scaleParaX * &
                                             (sg%horizonSimilarity(i) * (1.0D0 - sg%vertcalSimilarity(k)) + &
                                              sg%vertcalSimilarity(k))
                    END IF
                  END DO
                END IF
              END DO
            END IF

            IF ((sg%cell_stcl(6, i) .EQ. 0)) THEN
              ! Yuanfu Xie added the horizonSimilarity parameters: 2023/10/20
              ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
              DO it = 1, sg%tSlots
                IF (field%maskTemporal(it) .EQ. 1) THEN
                  DO k = 1, sg%vLevel
                    IF (field%maskVertical(k) .EQ. 1) THEN
                      ! field%data(:, i, :) = field%data(:, i, :) + ((swap(:, sg%cell_stcl(4, i), :) - swap(:, i, :))/(sg%cell_dist(4, i)))*(this%ScaleParaX)/(sg%cell_dist(4, i))
                      field%DATA(k, sg%cell_stcl(4, i), it) = field%DATA(k, sg%cell_stcl(4, i), it) &
                                                              + swap(k, i, it) * this%scaleParaX * &
                                                              (sg%horizonSimilarity(i) * (1.0D0 - sg%vertcalSimilarity(k)) + &
                                                               sg%vertcalSimilarity(k))
                      field%DATA(k, i, it) = field%DATA(k, i, it) &
                                             - swap(k, i, it) * this%scaleParaX * &
                                             (sg%horizonSimilarity(i) * (1.0D0 - sg%vertcalSimilarity(k)) + &
                                              sg%vertcalSimilarity(k))
                    END IF
                  END DO
                END IF
              END DO
            END IF
          END IF

          ! IF ((sg%cell_stcl(4, i) .EQ. 0) .OR. (sg%cell_stcl(6, i) .EQ. 0)) THEN
          !   field%data(:, i, :) = swap(:, i, :)
          ! END IF

          IF ((sg%cell_stcl(2, i) .NE. 0) .AND. (sg%cell_stcl(8, i) .NE. 0)) THEN
            ! Yuanfu Xie added the horizonSimilarity parameters: 2023/10/20
            ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
            ! field%data(:, i, :) = field%data(:, i, :) + ((swap(:, sg%cell_stcl(2, i), :) + swap(:, sg%cell_stcl(8, i), :) - &
            !                                               2*swap(:, i, :))/(sg%cell_dist(2, i)**2))*this%ScaleParaY
            DO it = 1, sg%tSlots
              IF (field%maskTemporal(it) .EQ. 1) THEN
                DO k = 1, sg%vLevel
                  IF (field%maskVertical(k) .EQ. 1) THEN
                    field%DATA(k, sg%cell_stcl(2, i), it) = field%DATA(k, sg%cell_stcl(2, i), it) + &
                                                            swap(k, i, it) * this%scaleParaY * &
                                                            (sg%horizonSimilarity(i) * (1.0D0 - sg%vertcalSimilarity(k)) + &
                                                             sg%vertcalSimilarity(k)) * this%c1
                    field%DATA(k, sg%cell_stcl(8, i), it) = field%DATA(k, sg%cell_stcl(8, i), it) + &
                                                            swap(k, i, it) * this%scaleParaY * &
                                                            (sg%horizonSimilarity(i) * (1.0D0 - sg%vertcalSimilarity(k)) + &
                                                             sg%vertcalSimilarity(k)) * this%c3
                    field%DATA(k, i, it) = field%DATA(k, i, it) + swap(k, i, it) * this%scaleParaY * &
                                           (sg%horizonSimilarity(i) * (1.0D0 - sg%vertcalSimilarity(k)) + &
                                            sg%vertcalSimilarity(k)) * this%c2
                  END IF
                END DO
              END IF
            END DO
          END IF

          ! IF ((sg%cell_stcl(2, i) .EQ. 0) .OR. (sg%cell_stcl(8, i) .EQ. 0)) THEN
          !   field%data(:, i, :) = swap(:, i, :)
          ! END IF

          IF (this%useNeumCondOnBndy) THEN
            IF ((sg%cell_stcl(2, i) .EQ. 0)) THEN
              ! Yuanfu Xie added the horizonSimilarity parameters: 2023/10/20
              ! field%data(:, i, :) = field%data(:, i, :) + ((swap(:, sg%cell_stcl(8, i), :) - swap(:, i, :))/(sg%cell_dist(8, i)))*(this%ScaleParaY)/(sg%cell_dist(8, i))
              ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
              DO it = 1, sg%tSlots
                IF (field%maskTemporal(it) .EQ. 1) THEN
                  DO k = 1, sg%vLevel
                    IF (field%maskVertical(k) .EQ. 1) THEN
                      field%DATA(k, sg%cell_stcl(8, i), it) = field%DATA(k, sg%cell_stcl(8, i), it) &
                                                              + swap(k, i, it) * this%scaleParaY * &
                                                              (sg%horizonSimilarity(i) * (1.0D0 - sg%vertcalSimilarity(k)) + &
                                                               sg%vertcalSimilarity(k))
                      field%DATA(k, i, it) = field%DATA(k, i, it) &
                                             - swap(k, i, it) * this%scaleParaY * &
                                             (sg%horizonSimilarity(i) * (1.0D0 - sg%vertcalSimilarity(k)) + &
                                              sg%vertcalSimilarity(k))
                    END IF
                  END DO
                END IF
              END DO
            END IF

            IF ((sg%cell_stcl(8, i) .EQ. 0)) THEN
              ! Yuanfu Xie added the horizonSimilarity parameters: 2023/10/20
              ! field%data(:, i, :) = field%data(:, i, :) + ((swap(:, sg%cell_stcl(2, i), :) - swap(:, i, :))/(sg%cell_dist(2, i)))*(this%ScaleParaY)/(sg%cell_dist(2, i))
              ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
              DO it = 1, sg%tSlots
                IF (field%maskTemporal(it) .EQ. 1) THEN
                  DO k = 1, sg%vLevel
                    IF (field%maskVertical(k) .EQ. 1) THEN
                      field%DATA(k, sg%cell_stcl(2, i), it) = field%DATA(k, sg%cell_stcl(2, i), it) &
                                                              + swap(k, i, it) * this%scaleParaY * &
                                                              (sg%horizonSimilarity(i) * (1.0D0 - sg%vertcalSimilarity(k)) + &
                                                               sg%vertcalSimilarity(k))
                      field%DATA(k, i, it) = field%DATA(k, i, it) &
                                             - swap(k, i, it) * this%scaleParaY * &
                                             (sg%horizonSimilarity(i) * (1.0D0 - sg%vertcalSimilarity(k)) + &
                                              sg%vertcalSimilarity(k))
                    END IF
                  END DO
                END IF
              END DO
            END IF
          END IF
        END DO
      END IF

      CALL sg%ExchangeMatOnHaloReverseSumForFieldGrid(sg%tSlots, sg%vLevel, field%DATA)! Exchange hale
      CALL sg%ExchangeMatOnHaloForFieldGrid(sg%tSlots, sg%vLevel, field%DATA)! Exchange hale

      ! Laplacian z direction
      ! This second-order scheme, refer to the docs by Xie and Qin, which will be added to the MOTOR-DA sooner.
      IF (this%sg%vLevel .GT. 3 .AND. (TRIM(field%Get_Name()) /= 'pres') .AND. (TRIM(field%Get_Name()) /= 'psl')) THEN
        ! Middle points, use i-1, i, i+1 and i+2

        IF (this%ScaleParaZ >= 500000) THEN
          DO i = 2, this%sg%vLevel - 2
            h1 = this%sg%sigma(i - 1) - this%sg%sigma(i)
            h2 = 0
            h3 = this%sg%sigma(i + 1) - this%sg%sigma(i)
            h4 = this%sg%sigma(i + 2) - this%sg%sigma(i)

            as(1, i) = (-2 * (h2 + h3 + h4) / ((h1 - h2) * (h1 - h3) * (h1 - h4)))
            as(2, i) = (2 * (h1 + h3 + h4) / ((h1 - h2) * (h2 - h3) * (h2 - h4)))
            as(3, i) = (2 * (h1 + h2 + h4) / ((h1 - h3) * (h3 - h2) * (h3 - h4)))
            as(4, i) = (2 * (h1 + h2 + h3) / ((h1 - h4) * (h4 - h2) * (h4 - h3)))

            ! PRINT *, as(1, i), as(2, i), as(3, i), as(4, i)
          END DO
          ! Last point, use i-2, i-1, i and i+1
          i = this%sg%vLevel - 1
          h1 = this%sg%sigma(i - 2) - this%sg%sigma(i)
          h2 = this%sg%sigma(i - 1) - this%sg%sigma(i)
          h3 = 0
          h4 = this%sg%sigma(i + 1) - this%sg%sigma(i)

          a1 = (-2 * (h2 + h3 + h4) / ((h1 - h2) * (h1 - h3) * (h1 - h4)))
          a2 = (2 * (h1 + h3 + h4) / ((h1 - h2) * (h2 - h3) * (h2 - h4)))
          a3 = (2 * (h1 + h2 + h4) / ((h1 - h3) * (h3 - h2) * (h3 - h4)))
          a4 = (2 * (h1 + h2 + h3) / ((h1 - h4) * (h4 - h2) * (h4 - h3)))

        ELSE
          as(1, :) = this%c1
          as(2, :) = this%c2
          as(3, :) = this%c3
          as(4, :) = 0.0D0
          a1 = 0.0D0
          a2 = this%c1
          a3 = this%c2
          a4 = this%c3
        END IF

        DO i = 2, this%sg%vLevel - 2
          IF (field%maskVertical(i) .EQ. 1) THEN
            DO it = 1, sg%tSlots
              ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
              IF (field%maskTemporal(it) .EQ. 1) THEN
                IF (field%maskHorizontal .EQ. 1) THEN
                  IF (i > 2) &
                    field%DATA(i - 1, :, it) = field%DATA(i - 1, :, it) + as(1, i) * swap(i, :, it) * this%ScaleParaZ

                  field%DATA(i, :, it) = field%DATA(i, :, it) + as(2, i) * swap(i, :, it) * this%ScaleParaZ
                  field%DATA(i + 1, :, it) = field%DATA(i + 1, :, it) + as(3, i) * swap(i, :, it) * this%ScaleParaZ
                  field%DATA(i + 2, :, it) = field%DATA(i + 2, :, it) + as(4, i) * swap(i, :, it) * this%ScaleParaZ
                END IF
              END IF
            END DO
          END IF
        END DO

        ! Last point, use i-2, i-1, i and i+1
        i = this%sg%vLevel - 1

        ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
        IF (field%maskVertical(i) .EQ. 1) THEN
          DO it = 1, sg%tSlots
            ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
            IF (field%maskTemporal(it) .EQ. 1) THEN
              IF (field%maskHorizontal .EQ. 1) THEN
                field%DATA(i - 2, :, it) = field%DATA(i - 2, :, it) + a1 * swap(i, :, it) * this%ScaleParaZ
                field%DATA(i - 1, :, it) = field%DATA(i - 1, :, it) + a2 * swap(i, :, it) * this%ScaleParaZ
                field%DATA(i, :, it) = field%DATA(i, :, it) + a3 * swap(i, :, it) * this%ScaleParaZ
                field%DATA(i + 1, :, it) = field%DATA(i + 1, :, it) + a4 * swap(i, :, it) * this%ScaleParaZ
              END IF
            END IF
          END DO
        END IF

        ! This set the open condition at the upper and bottom boundary, which may lead the unconvergent of the minization.
        ! i = 1
        ! IF (this%ScaleParaZ >= 500000) THEN
        ! field%data(1, :, :) = field%data(1, :, :) - swap(i, :, :)/(this%sg%sigma(2)-this%sg%sigma(1))*(this%ScaleParaZ)/ABS(this%sg%sigma(1)-this%sg%sigma(2))
        ! field%data(2, :, :) = field%data(2, :, :) + swap(i, :, :)/(this%sg%sigma(2)-this%sg%sigma(1))*(this%ScaleParaZ)/ABS(this%sg%sigma(1)-this%sg%sigma(2))
        ! ELSE
        !   field%data(1, :, :) = field%data(1, :, :) - swap(i, :, :)*(this%ScaleParaZ)
        !   field%data(2, :, :) = field%data(2, :, :) + swap(i, :, :)*(this%ScaleParaZ)
        ! END IF

        IF (this%useNeumCondOnBndy) THEN
          i = this%sg%vLevel
          ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
          IF (field%maskVertical(i) .EQ. 1) THEN
            DO it = 1, sg%tSlots
              ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
              IF (field%maskTemporal(it) .EQ. 1) THEN
                IF (field%maskHorizontal .EQ. 1) THEN
                  IF (this%ScaleParaZ >= 500000) THEN
                    field%DATA(this%sg%vLevel, :, it) = field%DATA(this%sg%vLevel, :, it) - swap(i, :, it) &
                                                        / (this%sg%sigma(this%sg%vLevel) - this%sg%sigma(this%sg%vLevel - 1)) * (this%ScaleParaZ) / ABS(this%sg%sigma(this%sg%vLevel) - this%sg%sigma(this%sg%vLevel - 1))

                    field%DATA(this%sg%vLevel - 1, :, it) = field%DATA(this%sg%vLevel - 1, :, it) + swap(i, :, it) &
                                                            / (this%sg%sigma(this%sg%vLevel) - this%sg%sigma(this%sg%vLevel - 1)) * (this%ScaleParaZ) / ABS(this%sg%sigma(this%sg%vLevel) - this%sg%sigma(this%sg%vLevel - 1))

                  ELSE
                    field%DATA(this%sg%vLevel, :, it) = field%DATA(this%sg%vLevel, :, it) - swap(i, :, it) &
                                                        * this%ScaleParaZ

                    field%DATA(this%sg%vLevel - 1, :, it) = field%DATA(this%sg%vLevel - 1, :, it) + swap(i, :, it) &
                                                            * this%ScaleParaZ
                  END IF
                END IF
              END IF
            END DO
          END IF
        END IF
      END IF

      ! Laplacian T direction
      IF (this%sg%tSlots .GT. 2) THEN
        BLOCK
          REAL(r_kind) :: dt

          dt = this%sg%tt(2) - this%sg%tt(1)
          DO it = 2, this%sg%tSlots - 1
            ! field%data(:, :, i) = field%data(:, :, i) + ((swap(:, :, i + 1) + swap(:, :, i - 1) - &
            !                                               2*swap(:, :, i)))/()/dt/dt*this%scaleParaT
            ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
            IF (field%maskTemporal(it) .EQ. 1) THEN
              DO k = 1, sg%vLevel
                IF (field%maskVertical(k) .EQ. 1) THEN
                  IF (field%maskHorizontal .EQ. 1) THEN
                    field%DATA(k, :, it) = field%DATA(k, :, it) + swap(k, :, it) * this%scaleParaT * this%c2
                    field%DATA(k, :, it + 1) = field%DATA(k, :, it + 1) + swap(k, :, it) * this%scaleParaT * this%c1
                    field%DATA(k, :, it - 1) = field%DATA(k, :, it - 1) + swap(k, :, it) * this%scaleParaT * this%c3
                  END IF
                END IF
              END DO
            END IF
          END DO

          IF (this%useNeumCondOnBndy) THEN
            ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
            IF (field%maskTemporal(1) .EQ. 1) THEN
              DO k = 1, sg%vLevel
                IF (field%maskVertical(k) .EQ. 1) THEN
                  IF (field%maskHorizontal .EQ. 1) THEN
                    field%DATA(k, :, 1) = field%DATA(k, :, 1) - swap(k, :, 1) * this%scaleParaT
                    field%DATA(k, :, 2) = field%DATA(k, :, 2) + swap(k, :, 1) * this%scaleParaT
                  END IF
                END IF
              END DO
            END IF

            ! Yuanfu Xie added the mask to the Laplacian operator: 2024/04/18
            IF (field%maskTemporal(this%sg%tSlots) .EQ. 1) THEN
              DO k = 1, sg%vLevel
                IF (field%maskVertical(k) .EQ. 1) THEN
                  IF (field%maskHorizontal .EQ. 1) THEN
                    field%DATA(k, :, this%sg%tSlots) = field%DATA(k, :, this%sg%tSlots) - &
                                                       swap(k, :, this%sg%tSlots) * this%scaleParaT

                    field%DATA(k, :, this%sg%tSlots - 1) = field%DATA(k, :, this%sg%tSlots - 1) + &
                                                           swap(k, :, this%sg%tSlots) * this%scaleParaT
                  END IF
                END IF
              END DO
            END IF
          END IF
        END BLOCK
      END IF
      ! CALL sg%ExchangeMatOnHaloForFieldGrid(sg%tSlots, sg%vLevel, field%data)! Exchange hale

      IF (field%Get_Name() == 'wwnd') THEN
        field%DATA(1, :, :) = 0.0D0; 
      END IF

      IF (field%Get_Name() == 'wwnd') THEN
        field%DATA(this%sg%vLevel, :, :) = 0.0D0; 
      END IF
      ! DO i = LBOUND(this%Y%ObsFields, 1), UBOUND(this%Y%ObsFields, 1)
      !   IF ((TRIM(this%Y%ObsFields(i)%Get_Name())) .EQ. (TRIM(field%Get_Name()))) THEN
      !     DO k = LBOUND(this%Y%ObsFields(i)%idx, 1), UBOUND(this%Y%ObsFields(i)%idx, 1)
      !       !  PRINT *,'X%fields(j)%Get_Name(): ',X%fields(j)%Get_Name(), sg%gLevel
      !       CALL field%Set_Value(this%Y%ObsFields(i)%idx(k), 0.0D0)
      !     END DO
      !   END IF
      ! END DO

      ! Reset all non-control gradient to zero by Yuanfu Xie 2024-07-09:
      IF (field%maskHorizontal .NE. 1) THEN
        field%DATA(:, :, :) = 0.0D0
      ELSE
        DO k = 1, this%sg%vLevel
          IF (field%maskVertical(k) .NE. 1) THEN
            field%DATA(k, :, :) = 0.0D0
          ELSE
            DO it = 1, this%sg%tSlots
              IF (field%maskTemporal(it) .NE. 1) &
                field%DATA(k, :, it) = 0.0D0
            END DO
          END IF
        END DO
      END IF

    END ASSOCIATE
    DEALLOCATE (swap)

  END SUBROUTINE sqrt_inverse_multiply_adjoint

  SUBROUTINE sqrt_inverse_multiply(this, field)
    IMPLICIT NONE
    CLASS(BFieldLaplace_t) :: this
    TYPE(Field_t), INTENT(INOUT) :: field
    INTEGER :: i, it, ni, k
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

      field%DATA = 0.0D0
      ! Laplacian x direction
      DO i = 1, this%sg%num_icell
!        IF (sg%cell_type(i) .EQ. 0) THEN
        IF ((sg%cell_stcl(4, i) .NE. 0) .AND. (sg%cell_stcl(6, i) .NE. 0)) THEN
          field%DATA(:, i, :) = field%DATA(:, i, :) + ((swap(:, sg%cell_stcl(4, i), :) * this%c1 + swap(:, sg%cell_stcl(6, i), :) * this%c3 + &
                                                        this%c2 * swap(:, i, :))) * this%scaleParaX
        END IF

        IF (this%useNeumCondOnBndy) THEN
          IF ((sg%cell_stcl(4, i) .EQ. 0)) THEN
            field%DATA(:, i, :) = field%DATA(:, i, :) + ((swap(:, sg%cell_stcl(6, i), :) - swap(:, i, :))) * this%scaleParaX * 1
          END IF

          IF ((sg%cell_stcl(6, i) .EQ. 0)) THEN
            field%DATA(:, i, :) = field%DATA(:, i, :) + ((swap(:, sg%cell_stcl(4, i), :) - swap(:, i, :))) * this%scaleParaX * 1
          END IF
        END IF

        ! Laplacian y direction
!        IF (sg%cell_type(i) .EQ. 0) THEN
        IF ((sg%cell_stcl(2, i) .NE. 0) .AND. (sg%cell_stcl(8, i) .NE. 0)) THEN
          field%DATA(:, i, :) = field%DATA(:, i, :) + ((swap(:, sg%cell_stcl(2, i), :) * this%c1 + swap(:, sg%cell_stcl(8, i), :) * this%c3 + &
                                                        this%c2 * swap(:, i, :))) * this%scaleParaY
        END IF

        IF (this%useNeumCondOnBndy) THEN
          IF ((sg%cell_stcl(2, i) .EQ. 0)) THEN
            field%DATA(:, i, :) = field%DATA(:, i, :) + ((swap(:, sg%cell_stcl(8, i), :) - swap(:, i, :))) * this%scaleParaY * 1
          END IF

          IF ((sg%cell_stcl(8, i) .EQ. 0)) THEN
            field%DATA(:, i, :) = field%DATA(:, i, :) + ((swap(:, sg%cell_stcl(2, i), :) - swap(:, i, :))) * this%scaleParaY * 1
          END IF
        END IF
      END DO

      ! Laplacian z direction
      ! This second-order scheme, refer to the docs by Xie and Qin, which will be added to the MOTOR-DA sooner.
      IF (this%sg%vLevel .GT. 3 .AND. (TRIM(field%Get_Name()) /= 'pres') .AND. (TRIM(field%Get_Name()) /= 'psl')) THEN
        ! Middle points, use i-1, i, i+1 and i+2

        IF (this%ScaleParaZ >= 500000) THEN
          DO i = 2, this%sg%vLevel - 2
            h1 = this%sg%sigma(i - 1) - this%sg%sigma(i)
            h2 = 0
            h3 = this%sg%sigma(i + 1) - this%sg%sigma(i)
            h4 = this%sg%sigma(i + 2) - this%sg%sigma(i)

            as(1, i) = (-2 * (h2 + h3 + h4) / ((h1 - h2) * (h1 - h3) * (h1 - h4)))
            as(2, i) = (2 * (h1 + h3 + h4) / ((h1 - h2) * (h2 - h3) * (h2 - h4)))
            as(3, i) = (2 * (h1 + h2 + h4) / ((h1 - h3) * (h3 - h2) * (h3 - h4)))
            as(4, i) = (2 * (h1 + h2 + h3) / ((h1 - h4) * (h4 - h2) * (h4 - h3)))

            ! PRINT *, as(1, i), as(2, i), as(3, i), as(4, i)
          END DO
          ! Last point, use i-2, i-1, i and i+1
          i = this%sg%vLevel - 1
          h1 = this%sg%sigma(i - 2) - this%sg%sigma(i)
          h2 = this%sg%sigma(i - 1) - this%sg%sigma(i)
          h3 = 0
          h4 = this%sg%sigma(i + 1) - this%sg%sigma(i)

          a1 = (-2 * (h2 + h3 + h4) / ((h1 - h2) * (h1 - h3) * (h1 - h4)))
          a2 = (2 * (h1 + h3 + h4) / ((h1 - h2) * (h2 - h3) * (h2 - h4)))
          a3 = (2 * (h1 + h2 + h4) / ((h1 - h3) * (h3 - h2) * (h3 - h4)))
          a4 = (2 * (h1 + h2 + h3) / ((h1 - h4) * (h4 - h2) * (h4 - h3)))

        ELSE
          as(1, :) = this%c1
          as(2, :) = this%c2
          as(3, :) = this%c3
          as(4, :) = 0.0D0
          a1 = 0.0D0
          a2 = this%c1
          a3 = this%c2
          a4 = this%c3
        END IF

        DO i = 2, this%sg%vLevel - 2
          field%DATA(i, :, :) = field%DATA(i, :, :) + (as(1, i) * swap(i - 1, :, :) + as(2, i) * swap(i, :, :) + &
                                                       as(3, i) * swap(i + 1, :, :) + as(4, i) * swap(i + 2, :, :)) * this%ScaleParaZ
        END DO

        i = this%sg%vLevel - 1
        ! PRINT *, 'Top', a2, a3, a4, h2
        field%DATA(i, :, :) = field%DATA(i, :, :) + (a1 * swap(i - 2, :, :) + a2 * swap(i - 1, :, :) + &
                                                     a3 * swap(i, :, :) + a4 * swap(i + 1, :, :)) * this%ScaleParaZ

        ! This set the open condition at the upper and bottom boundary, which may lead the unconvergent of the minization.
        ! i = 1
        ! IF (this%ScaleParaZ >= 500000) THEN
        ! field%data(i, :, :) = field%data(i, :, :) + (swap(2, :, :)-swap(1, :, :))/(this%sg%sigma(2)-this%sg%sigma(1))*(this%ScaleParaZ)/ABS(this%sg%sigma(1)-this%sg%sigma(2))
        ! ELSE
        !   field%data(i, :, :) = field%data(i, :, :) + (swap(2, :, :) - swap(1, :, :))*this%ScaleParaZ
        ! END IF
        ! !   !                                               ! print *, shape(field%data)

        IF (this%useNeumCondOnBndy) THEN
          i = this%sg%vLevel
          IF (this%ScaleParaZ >= 500000) THEN
            field%DATA(i, :, :) = field%DATA(i, :, :) + (swap(this%sg%vLevel - 1, :, :) - swap(this%sg%vLevel, :, :)) &
                                  / (this%sg%sigma(this%sg%vLevel) - this%sg%sigma(this%sg%vLevel - 1)) * (this%ScaleParaZ) / ABS(this%sg%sigma(this%sg%vLevel) - this%sg%sigma(this%sg%vLevel - 1))
          ELSE
            field%DATA(i, :, :) = field%DATA(i, :, :) + (swap(this%sg%vLevel - 1, :, :) - swap(this%sg%vLevel, :, :)) * this%ScaleParaZ
          END IF
        END IF

      END IF

      ! Laplacian T direction
      IF (this%sg%tSlots .GT. 2) THEN
        dt = this%sg%tt(2) - this%sg%tt(1)
        DO i = 2, this%sg%tSlots - 1
          field%DATA(:, :, i) = field%DATA(:, :, i) + this%c1 * ((swap(:, :, i + 1) + this%c3 * swap(:, :, i - 1) + &
                                                                  this%c2 * swap(:, :, i))) * this%scaleParaT
        END DO

        field%DATA(:, :, 1) = field%DATA(:, :, 1) + (swap(:, :, 2) - swap(:, :, 1)) * this%scaleParaT
        field%DATA(:, :, this%sg%tSlots) = field%DATA(:, :, this%sg%tSlots) + &
                                           (swap(:, :, this%sg%tSlots - 1) - swap(:, :, this%sg%tSlots)) * this%scaleParaT
      END IF

      ! IF (field%Get_Name() == 'wwnd') THEN
      !   field%DATA(1, :, :) = 0.0D0;
      ! END IF

      ! FORALL(i = 1:this%sg%num_icell, this%sg%cell_type(i)/=0) field%data(:, i, :) = 0.0D0

      ! DO i = LBOUND(this%Y%ObsFields, 1), UBOUND(this%Y%ObsFields, 1)
      !   IF ((TRIM(this%Y%ObsFields(i)%Get_Name())) .EQ. (TRIM(field%Get_Name()))) THEN
      !     DO k = LBOUND(this%Y%ObsFields(i)%idx, 1), UBOUND(this%Y%ObsFields(i)%idx, 1)
      !       !  PRINT *,'X%fields(j)%Get_Name(): ',X%fields(j)%Get_Name(), sg%gLevel
      !       CALL field%Set_Value(this%Y%ObsFields(i)%idx(k), 0.0D0)
      !     END DO
      !   END IF
      ! END DO

      ! Reset all non-control gradient to zero by Yuanfu Xie 2024-07-09:
      IF (field%maskHorizontal .NE. 1) THEN
        field%DATA(:, :, :) = 0.0D0
      ELSE
        DO k = 1, this%sg%vLevel
          IF (field%maskVertical(k) .NE. 1) THEN
            field%DATA(k, :, :) = 0.0D0
          ELSE
            DO it = 1, this%sg%tSlots
              IF (field%maskTemporal(it) .NE. 1) &
                field%DATA(k, :, it) = 0.0D0
            END DO
          END IF
        END DO
      END IF

      CALL sg%ExchangeMatOnHaloForFieldGrid(sg%tSlots, sg%vLevel, field%DATA)! Exchange hale

    END ASSOCIATE
    DEALLOCATE (swap)

  END SUBROUTINE sqrt_inverse_multiply

END MODULE BFieldLaplace_m
