!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.BMatrix.BField
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong CHEN, Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/4/23, @GBA-MWF, Shenzhen
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/10/19, @GBA-MWF, Shenzhen
!     for adding an option currentField preparing to add a conversion of pressure obs to log(pressure).
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE BFieldLaplace_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE Field_m, ONLY: Field_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE AuxTypeSG_m, ONLY: AuxTypeSG_t
  USE BFieldBase_m, ONLY: BFieldBase_t
  USE GenSgDiffCoef_m, ONLY: GenSgDiffCoef_t
  USE ObsSet_m, ONLY: ObsSet_t
  IMPLICIT NONE

  TYPE, EXTENDS(BFieldBase_t):: BFieldLaplace_t

    REAL(r_kind) :: scaleParaH

  CONTAINS
    PROCEDURE, PRIVATE :: loadBMatFiles

    PROCEDURE, PUBLIC :: sqrt_inverse_multiply
    PROCEDURE, PUBLIC :: sqrt_inverse_multiply_adjoint

    PROCEDURE :: sqrtInvMul => sqrt_inverse_multiply
    PROCEDURE :: sqrtInvMulAdj => sqrt_inverse_multiply_adjoint

    PROCEDURE, PUBLIC :: inverse_multiply

    PROCEDURE, PUBLIC :: Initialize

    PROCEDURE :: Laplacia_Terrain
    PROCEDURE :: FirstOrderDerivative
    PROCEDURE :: SecondOrderDerivative
    PROCEDURE :: Laplacia
    PROCEDURE :: Divergen
    PROCEDURE :: LAPLACIA_TERRAIN_B
    PROCEDURE :: FIRSTORDERDERIVATIVE_B
    PROCEDURE :: SECONDORDERDERIVATIVE_B
    PROCEDURE :: LAPLACIA_B
    PROCEDURE :: DIVERGEN_B

    FINAL :: destructor
  END TYPE BFieldLaplace_t

CONTAINS

  SUBROUTINE Initialize(this, configFile, sg, varName, Y)
    !USE NMLRead_m
    USE YAMLRead_m
    IMPLICIT NONE

    CLASS(BFieldLaplace_t) :: this
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    TYPE(GenSgDiffCoef_t) :: GenSgDiffCoef

    INTEGER(i_kind) :: ifile, numCell, mgStart
    TYPE(ObsSet_t), TARGET, OPTIONAL :: Y

    CHARACTER(Len=1024), INTENT(IN) :: configFile
    CHARACTER(Len=10), INTENT(IN) :: varName

    CALL this%AuxTypeSG_t%aux_initialize(sg)
    this%name = varName

    this%scaleParaX = 1.0D0
    this%scaleParaY = 1.0D0
    this%scaleParaZ = 1.0D0
    this%scaleParaT = 1.0D0
    this%scaleParaH = 1.0D0

    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaX', this%scaleParaX)
    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaY', this%scaleParaY)
    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaZ', this%scaleParaZ)
    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaT', this%scaleParaT)
    ifile = yaml_get_var(TRIM(configFile), 'BMat', 'ScaleParaH', this%scaleParaH)

    ! IF (this%scaleParaX >= 50) this%scaleParaX = this%scaleParaX/sg%dlon_dis/sg%dlon_dis
    ! IF (this%scaleParaY >= 50) this%scaleParaY = this%scaleParaY/sg%dlat_dis/sg%dlat_dis
    ! IF (this%ScaleParaT >= 50) THEN
    ! IF (sg%tSlots > 1) this%ScaleParaT = this%ScaleParaT/(this%sg%tt(2) - this%sg%tt(1))/(this%sg%tt(2) - this%sg%tt(1))
    ! END IF
    CALL this%loadBMatFiles(varName)

    ! this%scaleParaX = this%scaleParaX*0.1
    ! this%scaleParaY = this%scaleParaY*0.1
    ! this%scaleParaZ = this%scaleParaZ*0.1
    ! this%scaleParaT = this%scaleParaT*0.1
    ! this%scaleParaH = this%scaleParaH*0.1

    ! IF (sg%gLevel > 3) THEN
    !    this%scaleParaX = this%scaleParaX*(1.7**(sg%gLevel - 3))
    !    this%scaleParaY = this%scaleParaY*(1.7**(sg%gLevel - 3))
    !    this%scaleParaZ = this%scaleParaZ*(1.7**(sg%gLevel - 3))
    !    this%scaleParaT = this%scaleParaT*(1.7**(sg%gLevel - 3))
    !    this%scaleParaH = this%scaleParaH*(1.7**(sg%gLevel - 3))
    ! END IF

    CALL GenSgDiffCoef%GenFirstOrderFull(this%sg)
    CALL GenSgDiffCoef%GenSecondOrder(this%sg)
    CALL GenSgDiffCoef%GenHzParams(this%sg)
    CALL GenSgDiffCoef%GenLapaceParams(this%sg, 'DA')
! this%ScaleParaT = 0
    ! IF (TRIM(varname) == 'qvapor') this%scaleParaZ = this%scaleParaZ*2
    PRINT *, 'constructor', MINVAL(this%sigma)

  END SUBROUTINE

  SUBROUTINE loadBMatFiles(this, varName)
    CLASS(BFieldLaplace_t) :: this
    CHARACTER(LEN=10), INTENT(IN) :: varName
    INTEGER(i_kind) :: i, k

    ! Here is the code of monk for simulating the fileinput of B mat files
    IF (ALLOCATED(this%sigma)) DEALLOCATE (this%sigma)
    ALLOCATE (this%sigma(this%sg%vLevel, this%sg%num_cell, this%sg%tSlots))
    ! Set all the sigma to 1.0D0
    this%sigma = 1.0D0
    SELECT CASE (TRIM(varName))
    CASE ('pres')
      this%sigma = 100.0D0
    CASE ('qvapor')
      this%sigma = 0.001D0
      ! DO i = 1, this%sg%tSlots
      !   this%sigma(:, :, i) = this%sigma(: ,:, i)*this%sg%s1
      ! END DO

      DO i = 1, this%sg%tSlots
        DO k = 1, this%sg%num_cell
          this%sigma(:, k, i) = this%sg%s_qvapor
        END DO
      END DO

    END SELECT
    PRINT *, 'loadBMatFiles', MINVAL(this%sigma)

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

  SUBROUTINE sqrt_inverse_multiply(this, field)
    IMPLICIT NONE
    CLASS(BFieldLaplace_t) :: this
    TYPE(Field_t), INTENT(INOUT) :: field
    INTEGER :: i, ti, ni
    REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE :: swap, field_h, field_v, field_time

    ! IF(this%mpddSub%isBaseProc()) PRINT *, 'this%scalePara', this%scaleParaX, this%scaleParaY, this%scaleParaZ

    ! DEALLOCATE(this%sigma)
    !           ALLOCATE(this%sigma(this%sg%vLevel, this%sg%num_cell, this%sg%tSlots))
    !           this%sigma = 1.0

    ! PRINT*, 'sigma is ', this%sigma(1, 1, 1)

    ASSOCIATE (sg => this%sg, &
               tSlots => this%sg%tSlots, &
               num_icell => this%sg%num_icell, &
               num_cell => this%sg%num_cell, &
               vLevel => this%sg%vLevel)
      !            cell_type => sg%cell_type, &
      !            num_edge => sg%num_edge, &
      !            numQuadPerEdge => sg%numQuadPerEdge, &
      !            edge_stcl => sg%edge_stcl, &
      !            coef_norm => sg%coef_norm, &
      !            coef_norm => sg%coef_gl, &
      !            edge_leng => sg%edge_leng, &
      !            cell_area => sg%cell_area)

      ! 1/sigma*X
      ALLOCATE (swap(vLevel, num_cell, tSlots), field_h(vLevel, num_cell, tSlots), &
                field_v(vLevel, num_cell, tSlots), field_time(vLevel, num_cell, tSlots))
      swap = field%DATA / this%sigma
      CALL sg%ExchangeMatOnHaloForFieldGrid(sg%tSlots, sg%vLevel, swap)

      field%DATA = 0.0D0
      field_h = 0.0D0
      field_v = 0.0D0
      field_time = 0.0D0

      ! horizontal Laplace
      DO i = 1, tSlots
        CALL this%Laplacia_Terrain(swap(:, :, i), field_h(:, :, i))
      END DO
      CALL sg%ExchangeMatOnHaloForFieldGrid(sg%tSlots, sg%vLevel, field_h)

      ! Vertical Laplacian
      ! This second-order scheme, refer to the docs by Xie and Qin, which will be added to the MOTOR-DA sooner.
      IF (vLevel .GT. 4 .AND. (TRIM(field%Get_Name()) /= 'pres')) THEN
        ! Middle points, use i-1, i, i+1 and i+2
        DO i = 1, tSlots
          CALL this%SecondOrderDerivative(swap(:, :, i), sg%coef_scddif, field_v(:, :, i))
        END DO

      END IF

      ! Laplacian T direction
      IF (tSlots .GT. 2) THEN
        DO i = 2, tSlots - 1
          field_time(:, :, i) = ((swap(:, :, i + 1) + swap(:, :, i - 1) - &
                                  2 * swap(:, :, i)))
        END DO
      END IF

      IF (tSlots .GT. 3) THEN
        field_time(:, :, 1) = 2.0D0 * swap(:, :, 1) - 5.0D0 * swap(:, :, 2) + 4.0D0 * swap(:, :, 3) &
                              - swap(:, :, 4)
        field_time(:, :, tSlots) = 2.0D0 * swap(:, :, tSlots) - 5.0D0 * swap(:, :, tSlots - 1) &
                                   + 4.0D0 * swap(:, :, tSlots - 2) - swap(:, :, tSlots - 3)
        ! ELSE
        !    field_time(:, :, 1) = field_time(:, :, 2)
        !    field_time(:, :, tSlots) = field_time(:, :, 2)
      END IF

      field%DATA = field_h * this%scaleParaH + field_v * this%scaleParaZ + field_time * this%scaleParaT

      IF (field%Get_Name() == 'wwnd') THEN
        field%DATA(1, :, :) = 0.0D0; 
      END IF

      CALL sg%ExchangeMatOnHaloForFieldGrid(sg%tSlots, sg%vLevel, field%DATA)! Exchange hale

    END ASSOCIATE
    DEALLOCATE (swap, field_h, field_v, field_time)

  END SUBROUTINE sqrt_inverse_multiply

  SUBROUTINE Laplacia_Terrain(this, opr, value_out)

    IMPLICIT NONE

    CLASS(BFieldLaplace_t) :: this
    REAL(r_kind), INTENT(IN) :: opr(:, :)
    ! These operands are in cell space:
    REAL(r_kind), INTENT(OUT) :: value_out(:, :)
    REAL(r_kind), ALLOCATABLE :: Lopr(:, :), &
                                 F_oprsigma_z(:, :), &
                                 paropr_sigma_1rst(:, :), &
                                 paropr_sigma_2end(:, :)
    INTEGER(i_kind) :: i

    ASSOCIATE (num_icell => this%sg%num_icell, &
               num_cell => this%sg%num_cell, &
               cell_type => this%sg%cell_type, &
               vLevel => this%sg%vLevel, &
               coef_fstdif => this%sg%coef_fstdif, &
               coef_scddif => this%sg%coef_scddif, &
               z => this%sg%zHght, &
               hz => this%sg%hz, &
               parhz_sigma => this%sg%parHz_parsigma, &
               Lz => this%sg%Lz, &
               F_z_z => this%sg%F_z_z, &
               F_Hz_z => this%sg%F_Hz_z, &
               F_invHz_z => this%sg%F_invHz_z)

      ALLOCATE (Lopr(vLevel, num_cell), &
                F_oprsigma_z(vLevel, num_cell), &
                paropr_sigma_1rst(vLevel, num_cell), &
                paropr_sigma_2end(vLevel, num_cell))

      Lopr = 0.0D0; 
      F_oprsigma_z = 0.0D0
      paropr_sigma_1rst = 0.0D0
      paropr_sigma_2end = 0.0D0

      CALL this%FirstOrderDerivative(opr, coef_fstdif, paropr_sigma_1rst)
      CALL this%SecondOrderDerivative(opr, coef_scddif, paropr_sigma_2end)

      CALL this%Divergen(paropr_sigma_1rst, z, F_oprsigma_z)

      DO i = 1, num_icell
        IF (cell_type(i) .EQ. 0) THEN
          value_out(:, i) = Lopr(:, i) &
                            - Hz(:, i) * (z(:, i) * Hz(:, i) &
                                          * paropr_sigma_2end(:, i) &
                                          - paropr_sigma_1rst(:, i) &
                                          + z(:, i) &
                                          * paropr_sigma_1rst(:, i) &
                                          * parHz_sigma(:, i)) * Lz(:, i) &
                            + Hz(:, i) * (Hz(:, i) * paropr_sigma_2end(:, i) &
                                          + paropr_sigma_1rst(:, i) &
                                          * parHz_sigma(:, i)) * F_z_z(:, i) &
                            - 2 * Hz(:, i) * F_oprsigma_z(:, i) &
                            + Hz(:, i)**2 * paropr_sigma_1rst(:, i) * F_invHz_z(:, i) &
                            - paropr_sigma_1rst(:, i) * F_Hz_z(:, i)
        END IF
      END DO
      ! CALL sg%ExchangeMatOnHalo2D(vLevel, value_out)
    END ASSOCIATE
    DEALLOCATE (Lopr, &
                F_oprsigma_z, &
                paropr_sigma_1rst, &
                paropr_sigma_2end)
  END SUBROUTINE Laplacia_Terrain

  SUBROUTINE FirstOrderDerivative(this, A, coef, parA_sigma)

    IMPLICIT NONE

    CLASS(BFieldLaplace_t) :: this
    REAL(8), INTENT(IN) :: coef(:, :, :), A(:, :)
    REAL(8), INTENT(OUT) :: parA_sigma(:, :)
    INTEGER(4) :: i, k

    ASSOCIATE (num_icell => this%sg%num_icell, &
               num_cell => this%sg%num_cell, &
               cell_type => this%sg%cell_type, &
               vLevel => this%sg%vLevel)

      parA_sigma = 0.0D0
      DO i = 1, num_icell
        IF (cell_type(i) .EQ. 0) THEN
          DO k = 1, vLevel
            IF (k .EQ. 1) THEN
              parA_sigma(k, i) = SUM(coef(k, i, :) * A(k:k + 2, i))
            ELSEIF (k .EQ. vLevel) THEN
              parA_sigma(k, i) = SUM(coef(k, i, :) * A(k - 2:k, i))
            ELSE
              parA_sigma(k, i) = SUM(coef(k, i, :) * A(k - 1:k + 1, i))
            END IF
          END DO
        END IF
      END DO
      CALL this%sg%ExchangeMatOnHalo2D(vLevel, parA_sigma)

    END ASSOCIATE

  END SUBROUTINE FirstOrderDerivative

  SUBROUTINE SecondOrderDerivative(this, A, coef, parA_sigma_2nd)
    IMPLICIT NONE

    CLASS(BFieldLaplace_t) :: this

    REAL(8), INTENT(IN) :: coef(:, :, :), A(:, :)
    REAL(8), INTENT(INOUT) :: parA_sigma_2nd(:, :)
    INTEGER(4) :: i, k

    ASSOCIATE (num_icell => this%sg%num_icell, &
               num_cell => this%sg%num_cell, &
               cell_type => this%sg%cell_type, &
               vLevel => this%sg%vLevel)
      parA_sigma_2nd = 0.0D0
      PRINT *, 'second order is in'
      PRINT *, 'size of coef is ', SIZE(coef, dim=1), SIZE(coef, dim=2), SIZE(coef, dim=3)
      PRINT *, 'size of A is ', SIZE(A, dim=1), SIZE(A, dim=2)
      DO i = 1, num_icell
        IF (cell_type(i) .EQ. 0) THEN
          !  PRINT *, 'i is', i
          DO k = 1, vLevel
            IF (k .EQ. 1) THEN
              parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k:k + 4, i))
            ELSEIF (k .EQ. 2) THEN
              parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k - 1:k + 3, i))
            ELSEIF (k .EQ. vLevel) THEN
              parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k - 4:k, i))
            ELSEIF (k .EQ. vLevel - 1) THEN
              parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k - 3:k + 1, i))
            ELSE
              parA_sigma_2nd(k, i) = SUM(coef(k, i, :) * A(k - 2:k + 2, i))
            END IF
          END DO
        END IF
      END DO
      CALL this%sg%ExchangeMatOnHalo2D(vLevel, parA_sigma_2nd)
    END ASSOCIATE
    PRINT *, 'second order closed'
  END SUBROUTINE

  SUBROUTINE Laplacia(this, oprand, VALUE)

    IMPLICIT NONE
    CLASS(BFieldLaplace_t) :: this

    REAL(8), INTENT(IN) :: oprand(:, :)
    REAL(8), INTENT(OUT) :: VALUE(:, :)

    ! Local variables:
    INTEGER(4) :: i, ic, ie, ig
    REAL(8), ALLOCATABLE :: edge_value(:), quad_value(:)

    ASSOCIATE (sg => this%sg, &
               num_icell => this%sg%num_icell, &
               num_cell => this%sg%num_cell, &
               cell_type => this%sg%cell_type, &
               vLevel => this%sg%vLevel, &
               num_edge => this%sg%num_edge, &
               numquadperedge => this%sg%numquadperedge, &
               edge_stcl => this%sg%edge_stcl, &
               coef_norm => this%sg%coef_norm, &
               coef_gl => this%sg%coef_gl, &
               edge_leng => this%sg%edge_leng, &
               cell_area => this%sg%cell_area)

      ALLOCATE (edge_value(vLevel), quad_value(vLevel))

      ! Loop through all cells:
      DO ic = 1, num_icell  ! Apply to inner points only, icell has halo points etc.
        IF (cell_type(ic) .NE. 0) CYCLE ! Apply to interior points only
        VALUE(:, ic) = 0.0D0
        DO ie = 1, num_edge
          edge_value = 0.0D0
          DO ig = 1, numQuadPerEdge
            DO i = 1, UBOUND(edge_stcl, 1)
              quad_value = quad_value + &
                           oprand(:, edge_stcl(i, ie, ic)) * coef_norm(i, ig, ie, ic)
            END DO
            edge_value = edge_value + coef_gl(ig) * quad_value
          END DO
          VALUE(:, ic) = VALUE(:, ic) + edge_value * edge_leng(ie, ic)
        END DO
        VALUE(:, ic) = VALUE(:, ic) / cell_area(ic)
      END DO

      DEALLOCATE (edge_value, quad_value)

    END ASSOCIATE

  END SUBROUTINE Laplacia

  SUBROUTINE Divergen(this, opr_left, opr_right, VALUE)

    IMPLICIT NONE

    CLASS(BFieldLaplace_t) :: this
    REAL(8), INTENT(IN) :: opr_left(:, :), opr_right(:, :)
    REAL(8), INTENT(OUT) :: VALUE(:, :)

    ! Local variables:
    INTEGER(4) :: i, ic, ie, ig
    REAL(8), ALLOCATABLE :: edge_value(:), left_value(:), rightvalue(:)

    ASSOCIATE (sg => this%sg, &
               num_icell => this%sg%num_icell, &
               num_cell => this%sg%num_cell, &
               cell_type => this%sg%cell_type, &
               vLevel => this%sg%vLevel, &
               num_edge => this%sg%num_edge, &
               numquadperedge => this%sg%numquadperedge, &
               edge_stcl => this%sg%edge_stcl, &
               coef_norm => this%sg%coef_norm, &
               coef_gl => this%sg%coef_gl, &
               edge_leng => this%sg%edge_leng, &
               cell_area => this%sg%cell_area, &
               coef_func => this%sg%coef_func)

      ALLOCATE (edge_value(vLevel), left_value(vLevel), rightvalue(vLevel))

      ! Loop through all cells:
      DO ic = 1, num_icell  ! Apply to inner points only, icell has halo points etc.
        IF (cell_type(ic) .NE. 0) CYCLE ! Apply to interior points only
        VALUE(:, ic) = 0.0D0
        DO ie = 1, num_edge
          edge_value = 0.0D0
          left_value = 0.0D0
          rightvalue = 0.0D0
          DO ig = 1, numQuadPerEdge
            DO i = 1, UBOUND(edge_stcl, 1)
              left_value = left_value + &
                           opr_left(:, edge_stcl(i, ie, ic)) * coef_func(i, ig, ie, ic)
              rightvalue = rightvalue + &
                           opr_right(:, edge_stcl(i, ie, ic)) * coef_norm(i, ig, ie, ic)
            END DO
            edge_value = edge_value + coef_gl(ig) * left_value * rightvalue
          END DO
          VALUE(:, ic) = VALUE(:, ic) + edge_value * edge_leng(ie, ic)
        END DO
        VALUE(:, ic) = VALUE(:, ic) / cell_area(ic)
      END DO

    END ASSOCIATE

    DEALLOCATE (edge_value, left_value, rightvalue)

  END SUBROUTINE Divergen

  SUBROUTINE sqrt_inverse_multiply_adjoint(this, field, currentfield)
    IMPLICIT NONE

    CLASS(BFieldLaplace_t) :: this
    TYPE(Field_t), INTENT(INOUT) :: field
    TYPE(Field_t), OPTIONAL, INTENT(IN) :: currentfield
    ! REAL*8, intent(inout) :: field1b(:,:,:)
    ! REAL*8, INTENT(INOUT) :: field2b(:,:,:)
    INTEGER :: i, it
    REAL*8, DIMENSION(:, :, :), ALLOCATABLE :: swapb, field_hb, &
                                               field_vb, field_timeb, &
                                               field1b
    INTEGER*4 :: branch

    ASSOCIATE (sg => this%sg, &
               tSlots => this%sg%tSlots, &
               num_icell => this%sg%num_icell, &
               num_cell => this%sg%num_cell, &
               coef_scddif => this%sg%coef_scddif, &
               vLevel => this%sg%vLevel, &
               field2b => field%DATA)

      ALLOCATE (swapb(vlevel, num_cell, tslots), field_hb(vlevel, num_cell, tslots), &
                field_vb(vlevel, num_cell, tslots), field_timeb(vlevel, num_cell, tslots), &
                field1b(vlevel, num_cell, tslots))

! Vertical Laplacian
      IF (vLevel .GT. 4 .AND. (TRIM(field%Get_Name()) /= 'pres')) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
! Laplacian T direction
      IF (tSlots .GT. 2) THEN
        ! i = tslots
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (tslots .GT. 3) THEN
        CALL PUSHCONTROL1B(0)
      ELSE
        CALL PUSHCONTROL1B(1)
      END IF
      IF (field%Get_Name() == 'wwnd') field2b(1, :, :) = 0.0_8
      field_hb = 0.0_8
      field_vb = 0.0_8
      field_timeb = 0.0_8
      field_hb = this%scaleParaH * field2b
      field_vb = this%scaleparaz * field2b
      field_timeb = this%scaleparat * field2b
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        swapb = 0.0_8
        swapb(:, :, tslots) = swapb(:, :, tslots) + 2.0D0 * field_timeb(:, :, &
    &     tslots)
        swapb(:, :, tslots - 1) = swapb(:, :, tslots - 1) - 5.0D0 * field_timeb(:&
    &     , :, tslots)
        swapb(:, :, tslots - 2) = swapb(:, :, tslots - 2) + 4.0D0 * field_timeb(:&
    &     , :, tslots)
        swapb(:, :, tslots - 3) = swapb(:, :, tslots - 3) - field_timeb(:, :, &
    &     tslots)
        field_timeb(:, :, tslots) = 0.0_8
        swapb(:, :, 1) = swapb(:, :, 1) + 2.0D0 * field_timeb(:, :, 1)
        swapb(:, :, 2) = swapb(:, :, 2) - 5.0D0 * field_timeb(:, :, 1)
        swapb(:, :, 3) = swapb(:, :, 3) + 4.0D0 * field_timeb(:, :, 1)
        swapb(:, :, 4) = swapb(:, :, 4) - field_timeb(:, :, 1)
        field_timeb(:, :, 1) = 0.0_8
      ELSE
        swapb = 0.0_8
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO i = tslots - 1, 2, -1
          swapb(:, :, i + 1) = swapb(:, :, i + 1) + field_timeb(:, :, i)
          swapb(:, :, i - 1) = swapb(:, :, i - 1) + field_timeb(:, :, i)
          swapb(:, :, i) = swapb(:, :, i) - 2 * field_timeb(:, :, i)
          field_timeb(:, :, i) = 0.0_8
        END DO
      END IF
      CALL POPCONTROL1B(branch)
      IF (branch .EQ. 0) THEN
        DO it = tslots, 1, -1
          CALL this%SECONDORDERDERIVATIVE_B(swapb(:, :, it), &
    &                            coef_scddif, field_vb(:, :, it))
          field_vb(:, :, it) = 0.0_8
        END DO
      END IF
      DO it = tslots, 1, -1
        CALL this%LAPLACIA_TERRAIN_B(swapb(:, :, it), field_hb(:, :, it))
      END DO
      field1b = 0.0_8
      field1b = swapb / this%sigma
      field2b = field1b
    END ASSOCIATE
    DEALLOCATE (swapb, field_hb, &
                field_vb, field_timeb, field1b)
  END SUBROUTINE sqrt_inverse_multiply_adjoint

  SUBROUTINE LAPLACIA_TERRAIN_B(this, oprb, value_outb)
    IMPLICIT NONE

    CLASS(BFieldLaplace_t) :: this

    REAL*8, INTENT(INOUT) :: oprb(:, :)
! These operands are in cell space:
    REAL*8 :: value_outb(:, :)
    REAL*8, ALLOCATABLE, DIMENSION(:, :) :: loprb, f_oprsigma_zb, &
  & paropr_sigma_1rstb, paropr_sigma_2endb

    INTEGER*4 :: i
    REAL*8, DIMENSION(:), ALLOCATABLE :: tempb
    REAL*8, DIMENSION(:), ALLOCATABLE :: tempb0
    INTEGER*4 :: branch

    ASSOCIATE (tSlots => this%sg%tSlots, &
               num_icell => this%sg%num_icell, &
               num_cell => this%sg%num_cell, &
               vLevel => this%sg%vLevel, &
               cell_type => this%sg%cell_type, &
               coef_fstdif => this%sg%coef_fstdif, &
               coef_scddif => this%sg%coef_scddif, &
               z => this%sg%zHght, &
               hz => this%sg%hz, &
               parhz_sigma => this%sg%parHz_parsigma, &
               Lz => this%sg%Lz, &
               F_z_z => this%sg%F_z_z, &
               F_Hz_z => this%sg%F_Hz_z, &
               F_invHz_z => this%sg%F_invHz_z)

      ALLOCATE (loprb(vLevel, num_cell), f_oprsigma_zb(vLevel, num_cell), &
                paropr_sigma_1rstb(vLevel, num_cell), paropr_sigma_2endb(vLevel, num_cell), &
                tempb(vLevel), tempb0(vLevel))

      DO i = 1, num_icell
        IF (cell_type(i) .EQ. 0) THEN
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
      paropr_sigma_2endb = 0.0_8
      f_oprsigma_zb = 0.0_8
      loprb = 0.0_8
      paropr_sigma_1rstb = 0.0_8
      DO i = num_icell, 1, -1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          loprb(:, i) = loprb(:, i) + value_outb(:, i)
          tempb = -(hz(:, i) * lz(:, i) * value_outb(:, i))
          tempb0 = hz(:, i) * f_z_z(:, i) * value_outb(:, i)
          paropr_sigma_1rstb(:, i) = paropr_sigma_1rstb(:, i) + (f_invhz_z(:&
    &       , i) * hz(:, i)**2 - f_hz_z(:, i)) * value_outb(:, i) + parhz_sigma(:&
    &       , i) * tempb0 + (z(:, i) * parhz_sigma(:, i) - 1.0) * tempb
          f_oprsigma_zb(:, i) = f_oprsigma_zb(:, i) - hz(:, i) * 2 * value_outb(&
    &       :, i)
          value_outb(:, i) = 0.0_8
          paropr_sigma_2endb(:, i) = paropr_sigma_2endb(:, i) + hz(:, i) *&
    &       tempb0 + z(:, i) * hz(:, i) * tempb
        END IF
      END DO
      CALL this%DIVERGEN_B(paropr_sigma_1rstb, z, f_oprsigma_zb)
      CALL this%LAPLACIA_B(oprb, loprb)
      CALL this%SECONDORDERDERIVATIVE_B(oprb, coef_scddif, paropr_sigma_2endb)
      CALL this%FIRSTORDERDERIVATIVE_B(oprb, coef_fstdif, paropr_sigma_1rstb)
    END ASSOCIATE

    DEALLOCATE (loprb, f_oprsigma_zb, &
                paropr_sigma_1rstb, paropr_sigma_2endb)
  END SUBROUTINE LAPLACIA_TERRAIN_B

  SUBROUTINE FIRSTORDERDERIVATIVE_B(this, ab, coef, para_sigmab)
    IMPLICIT NONE

    CLASS(BFieldLaplace_t) :: this
    REAL*8, INTENT(IN) :: coef(:, :, :)
    REAL*8 :: ab(:, :)
    REAL*8 :: para_sigmab(:, :)
    INTEGER*4 :: i, k
    INTRINSIC SUM
    INTEGER*4 :: branch

    ASSOCIATE (num_icell => this%sg%num_icell, &
               num_cell => this%sg%num_cell, &
               cell_type => this%sg%cell_type, &
               vLevel => this%sg%vLevel)

      DO i = 1, num_icell
        IF (cell_type(i) .EQ. 0) THEN
          DO k = 1, vlevel
            IF (k .EQ. 1) THEN
              CALL PUSHCONTROL2B(2)
            ELSE IF (k .EQ. vlevel) THEN
              CALL PUSHCONTROL2B(1)
            ELSE
              CALL PUSHCONTROL2B(0)
            END IF
          END DO
          CALL PUSHCONTROL1B(1)
        ELSE
          CALL PUSHCONTROL1B(0)
        END IF
      END DO
      DO i = num_icell, 1, -1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          DO k = vlevel, 1, -1
            CALL POPCONTROL2B(branch)
            IF (branch .EQ. 0) THEN
              ab(k - 1:k + 1, i) = ab(k - 1:k + 1, i) + coef(k, i, :) * para_sigmab(k&
    &           , i)
              para_sigmab(k, i) = 0.0_8
            ELSE IF (branch .EQ. 1) THEN
              ab(k - 2:k, i) = ab(k - 2:k, i) + coef(k, i, :) * para_sigmab(k, i)
              para_sigmab(k, i) = 0.0_8
            ELSE
              ab(k:k + 2, i) = ab(k:k + 2, i) + coef(k, i, :) * para_sigmab(k, i)
              para_sigmab(k, i) = 0.0_8
            END IF
          END DO
        END IF
      END DO

    END ASSOCIATE
  END SUBROUTINE FIRSTORDERDERIVATIVE_B

  SUBROUTINE SECONDORDERDERIVATIVE_B(this, ab, coef, para_sigma_2ndb)
    IMPLICIT NONE

    CLASS(BFieldLaplace_t) :: this
    REAL*8, INTENT(IN) :: coef(:, :, :)
    REAL*8 :: ab(:, :)
    REAL*8 :: para_sigma_2ndb(:, :)
    INTEGER*4 :: i, k
    INTRINSIC SUM
    INTEGER*4 :: branch

    ASSOCIATE (num_icell => this%sg%num_icell, &
               num_cell => this%sg%num_cell, &
               cell_type => this%sg%cell_type, &
               vLevel => this%sg%vLevel)
    DO i = 1, num_icell
      IF (cell_type(i) .EQ. 0) THEN
        DO k = 1, vlevel
          IF (k .EQ. 1) THEN
            CALL PUSHCONTROL3B(4)
          ELSE IF (k .EQ. 2) THEN
            CALL PUSHCONTROL3B(3)
          ELSE IF (k .EQ. vlevel) THEN
            CALL PUSHCONTROL3B(2)
          ELSE IF (k .EQ. vlevel - 1) THEN
            CALL PUSHCONTROL3B(1)
          ELSE
            CALL PUSHCONTROL3B(0)
          END IF
        END DO
        CALL PUSHCONTROL1B(1)
      ELSE
        CALL PUSHCONTROL1B(0)
      END IF
    END DO
    DO i = num_icell, 1, -1
      CALL POPCONTROL1B(branch)
      IF (branch .NE. 0) THEN
        DO k = vlevel, 1, -1
          CALL POPCONTROL3B(branch)
          IF (branch .LT. 2) THEN
            IF (branch .EQ. 0) THEN
              ab(k - 2:k + 2, i) = ab(k - 2:k + 2, i) + coef(k, i, :) *&
  &             para_sigma_2ndb(k, i)
              para_sigma_2ndb(k, i) = 0.0_8
            ELSE
              ab(k - 3:k + 1, i) = ab(k - 3:k + 1, i) + coef(k, i, :) *&
  &             para_sigma_2ndb(k, i)
              para_sigma_2ndb(k, i) = 0.0_8
            END IF
          ELSE IF (branch .EQ. 2) THEN
            ab(k - 4:k, i) = ab(k - 4:k, i) + coef(k, i, :) * para_sigma_2ndb(k&
  &           , i)
            para_sigma_2ndb(k, i) = 0.0_8
          ELSE IF (branch .EQ. 3) THEN
            ab(k - 1:k + 3, i) = ab(k - 1:k + 3, i) + coef(k, i, :) *&
  &           para_sigma_2ndb(k, i)
            para_sigma_2ndb(k, i) = 0.0_8
          ELSE
            ab(k:k + 4, i) = ab(k:k + 4, i) + coef(k, i, :) * para_sigma_2ndb(k&
  &           , i)
            para_sigma_2ndb(k, i) = 0.0_8
          END IF
        END DO
      END IF
    END DO

    END ASSOCIATE
  END SUBROUTINE SECONDORDERDERIVATIVE_B

  SUBROUTINE LAPLACIA_B(this, oprandb, valueb)
    IMPLICIT NONE

    CLASS(BFieldLaplace_t) :: this
    REAL*8 :: oprandb(:, :)
    REAL*8 :: valueb(:, :)
! Local variables:
    INTEGER*4 :: i, ic, ie, ig
    REAL*8, ALLOCATABLE :: edge_valueb(:), quad_valueb(:)
    INTRINSIC UBOUND
    INTEGER*4 :: ad_to
    INTEGER*4 :: branch

    ASSOCIATE (num_icell => this%sg%num_icell, &
               num_cell => this%sg%num_cell, &
               cell_type => this%sg%cell_type, &
               vLevel => this%sg%vLevel, &
               num_edge => this%sg%num_edge, &
               numquadperedge => this%sg%numquadperedge, &
               edge_stcl => this%sg%edge_stcl, &
               coef_norm => this%sg%coef_norm, &
               coef_gl => this%sg%coef_gl, &
               edge_leng => this%sg%edge_leng, &
               cell_area => this%sg%cell_area)

      ALLOCATE (edge_valueb(vLevel), quad_valueb(vLevel))

! Loop through all cells:
! Apply to inner points only, icell has halo points etc.
      DO ic = 1, num_icell
! Apply to interior points only
        IF (cell_type(ic) .NE. 0) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          DO ie = 1, num_edge
            DO ig = 1, numquadperedge
              DO i = 1, UBOUND(edge_stcl, 1)

              END DO
              CALL PUSHINTEGER4(i - 1)
            END DO
          END DO
          CALL PUSHCONTROL1B(1)
        END IF
      END DO
      quad_valueb = 0.0_8
      DO ic = num_icell, 1, -1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          valueb(:, ic) = valueb(:, ic) / cell_area(ic)
          DO ie = num_edge, 1, -1
            edge_valueb = 0.0_8
            edge_valueb = edge_leng(ie, ic) * valueb(:, ic)
            DO ig = numquadperedge, 1, -1
              quad_valueb = quad_valueb + coef_gl(ig) * edge_valueb
              CALL POPINTEGER4(ad_to)
              DO i = ad_to, 1, -1
                oprandb(:, edge_stcl(i, ie, ic)) = oprandb(:, edge_stcl(i, &
    &             ie, ic)) + coef_norm(i, ig, ie, ic) * quad_valueb
              END DO
            END DO
          END DO
          valueb(:, ic) = 0.0_8
        END IF
      END DO

      DEALLOCATE (edge_valueb, quad_valueb)

    END ASSOCIATE
  END SUBROUTINE LAPLACIA_B

  SUBROUTINE DIVERGEN_B(this, opr_leftb, opr_right, valueb)
    IMPLICIT NONE

    CLASS(BFieldLaplace_t) :: this
    REAL*8, INTENT(IN) :: opr_right(:, :)
    REAL*8, INTENT(INOUT) :: opr_leftb(:, :)
    REAL*8, INTENT(INOUT) :: valueb(:, :)
! Local variables:
    INTEGER*4 :: i, ic, ie, ig
    REAL*8, ALLOCATABLE :: rightvalue(:), edge_valueb(:), left_valueb(:)
    INTRINSIC UBOUND
    INTEGER*4 :: ad_to
    INTEGER*4 :: branch

    ASSOCIATE (num_icell => this%sg%num_icell, &
               num_cell => this%sg%num_cell, &
               cell_type => this%sg%cell_type, &
               vLevel => this%sg%vLevel, &
               num_edge => this%sg%num_edge, &
               numquadperedge => this%sg%numquadperedge, &
               edge_stcl => this%sg%edge_stcl, &
               coef_norm => this%sg%coef_norm, &
               coef_gl => this%sg%coef_gl, &
               edge_leng => this%sg%edge_leng, &
               cell_area => this%sg%cell_area, &
               coef_func => this%sg%coef_func)

      ALLOCATE (rightvalue(vlevel), edge_valueb(vlevel), left_valueb(vlevel))

! Loop through all cells:
! Apply to inner points only, icell has halo points etc.
      DO ic = 1, num_icell
! Apply to interior points only
        IF (cell_type(ic) .NE. 0) THEN
          CALL PUSHCONTROL1B(0)
        ELSE
          DO ie = 1, num_edge
            CALL PUSHREAL8ARRAY(rightvalue, vlevel)
            rightvalue = 0.0D0
            DO ig = 1, numquadperedge
              DO i = 1, UBOUND(edge_stcl, 1)
                CALL PUSHREAL8ARRAY(rightvalue, vlevel)
                rightvalue = rightvalue + opr_right(:, edge_stcl(i, ie, ic))&
                            &             * coef_norm(i, ig, ie, ic)
              END DO
              CALL PUSHINTEGER4(i - 1)
            END DO
          END DO
          CALL PUSHCONTROL1B(1)
        END IF
      END DO
      DO ic = num_icell, 1, -1
        CALL POPCONTROL1B(branch)
        IF (branch .NE. 0) THEN
          valueb(:, ic) = valueb(:, ic) / cell_area(ic)
          DO ie = num_edge, 1, -1
            edge_valueb = 0.0_8
            edge_valueb = edge_leng(ie, ic) * valueb(:, ic)
            left_valueb = 0.0_8
            DO ig = numquadperedge, 1, -1
              left_valueb = left_valueb + rightvalue * coef_gl(ig) * edge_valueb
              CALL POPINTEGER4(ad_to)
              DO i = ad_to, 1, -1
                CALL POPREAL8ARRAY(rightvalue, vlevel)
                opr_leftb(:, edge_stcl(i, ie, ic)) = opr_leftb(:, edge_stcl(&
    &             i, ie, ic)) + coef_func(i, ig, ie, ic) * left_valueb
              END DO
            END DO
            CALL POPREAL8ARRAY(rightvalue, vlevel)
          END DO
          valueb(:, ic) = 0.0_8
        END IF
      END DO
      DEALLOCATE (rightvalue, edge_valueb, left_valueb)
    END ASSOCIATE
  END SUBROUTINE DIVERGEN_B

END MODULE BFieldLaplace_m

