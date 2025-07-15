!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yongjian Huang
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yongjian Huang, 2024-09-13, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module forwards temprature, pressure, specific humility of a model state
!! toward refractivity in observation space
MODULE OprGNSSRefrac_m
  USE State_m, ONLY: State_t
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE ObsSet_m, ONLY: ObsSet_t
  USE M2OBase_m, ONLY: M2OBase_t
  USE parameters_m, ONLY: degree2radian
  USE YAMLRead_m, ONLY: yaml_get_var

  TYPE, EXTENDS(M2OBase_t) :: OprGNSSRefrac_t
    TYPE(ObsSet_t), POINTER :: Y
    TYPE(State_t), POINTER :: Xb

    REAL(r_kind) :: k1, k2, k3
    REAL(r_kind) :: R_dry, R_vap, epsilon_water
    LOGICAL :: use_qvapor_ctl, use_pres_ctl
    !> gnss_adj_mode: 1: update temp, 2: update temp and vapor, others: update all
    INTEGER(i_kind) :: gnss_adj_mode
  CONTAINS

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
    FINAL :: destructor
  END TYPE OprGNSSRefrac_t

  INTERFACE OprGNSSRefrac_t
    PROCEDURE :: constructor
  END INTERFACE OprGNSSRefrac_t

CONTAINS

  FUNCTION constructor(configFile, X, Y) RESULT(this)
    IMPLICIT NONE
    TYPE(OprGNSSRefrac_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X
    TYPE(ObsSet_t), OPTIONAL, TARGET, INTENT(IN) :: Y

    INTEGER(i_kind) :: istatus

    this%Xb => X

    BLOCK
    INTEGER(i_kind) :: i

    ! DO i = LBOUND(X%Fields, 1), UBOUND(X%Fields, 1)
    !   PRINT*, 'HYJ+++', X%Fields(i)%Get_Name(), MAXVAL(X%Fields(i)%DATA), MINVAL(X%Fields(i)%DATA)
    ! END DO

    END BLOCK
    
    IF (PRESENT(Y)) THEN
      this%Y => Y

      !! constants for ideal gas refractivity equation
      !! N = k1 * pdry/T + k2* pwvp/T^2 + k3 * pwvp/T
      this%k1 = 0.7760D0
      this%k2 = 3730.0D0
      this%k3 = 0.7760D0

      this%R_dry = 287.0597D0
      this%R_vap = 461.5250D0

      WRITE (*, 1) this%k1, this%k2, this%k3, this%R_dry, this%R_vap
  1   FORMAT("GNSS operator parameters [k1, k2, k3, R_dry, R_vap]: ", 5f5.5)

      istatus = yaml_get_var(TRIM(configFile), 'IO', 'GNSSRO_adj_mode', this%gnss_adj_mode)
      IF (istatus .EQ. 0) THEN
        WRITE (*, 4) this%gnss_adj_mode
    4   FORMAT('gnss ro mode is set to ', I3)
      ELSE
        this%gnss_adj_mode = 2
        WRITE (*, 5) this%gnss_adj_mode
    5   FORMAT('gnss ro mode is set to DEFAULT ', I3)
      END IF 

      this%epsilon_water = this%R_dry / this%R_vap
    ELSE
      WRITE (*, 99)
99     FORMAT('OprGNSSRefrac - transFwdTanLinear_opr requires OprRefrac is constructed by option Y but missing! STOP')
      STOP
    END IF
  END FUNCTION constructor

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    TYPE(OprGNSSRefrac_t), INTENT(INOUT) :: this
    INTEGER(i_kind) :: i
  END SUBROUTINE destructor

  SUBROUTINE transFwdNonLinear(this, X, Y)
    IMPLICIT NONE
    CLASS(OprGNSSRefrac_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(ObsSet_t), INTENT(INOUT) :: Y

    INTEGER(i_kind) :: refrac_ind, pres_ind, temp_ind, shum_ind, k, pres_ind_bkg, shum_ind_bkg
    REAL(r_kind) :: pwvp_val, pdry_val, pres_val, temp_val, shum_val, refrac_val

    refrac_ind = Y%getObsIdx('refractivity')

    pres_ind = X%getVarIdx('pres')
    temp_ind = X%getVarIdx('temp')
    shum_ind = X%getVarIdx('qvapor')


    pres_ind_bkg = this%Xb%getVarIdx('pres')
    shum_ind_bkg = this%Xb%getVarIdx('qvapor')



    ! If needed variables are not present:
    IF ((refrac_ind .EQ. 0) .OR. &
        (pres_ind .EQ. 0) .OR. &
        (temp_ind .EQ. 0) .OR. &
        (shum_ind .EQ. 0)) THEN

      PRINT *, 'OprGNSSRefrac%transFwdNonLinear variables not present, SKIP!'

      RETURN
    END IF

    ASSOCIATE (pres => X%Fields(pres_ind), &
               temp => X%Fields(temp_ind), &
               shum => X%Fields(shum_ind))

      ASSOCIATE (refrac => Y%ObsFields(refrac_ind))

        DO k = LBOUND(refrac%idx, 1), UBOUND(refrac%idx, 1)


          IF (this%gnss_adj_mode == 1) THEN      
            shum_val = this%Xb%Fields(shum_ind_bkg)%Get_Value(refrac%idx(k))
            pres_val = this%Xb%Fields(pres_ind_bkg)%Get_Value(refrac%idx(k)) 
            temp_val = temp%Get_Value(refrac%idx(k))
            ! shum_val = shum%Get_Value(refrac%idx(k))
            ! pres_val = pres%Get_Value(refrac%idx(k))
          ELSE IF (this%gnss_adj_mode == 2) THEN
            shum_val = shum%Get_Value(refrac%idx(k))
            temp_val = temp%Get_Value(refrac%idx(k))
            pres_val = this%Xb%Fields(pres_ind_bkg)%Get_Value(refrac%idx(k)) 
          ELSE
            temp_val = temp%Get_Value(refrac%idx(k))
            shum_val = shum%Get_Value(refrac%idx(k))
            pres_val = pres%Get_Value(refrac%idx(k))
          END IF
          
          ! temp_val = temp%Get_Value(refrac%idx(k))
          ! shum_val = shum%Get_Value(refrac%idx(k))
          ! pres_val = pres%Get_Value(refrac%idx(k))

          ! PRINT*, 'HYJ+++: OPR forward: ', pres_val, temp_val, shum_val

          ! dry pres and partial pres of vapor
          pwvp_val = pres_val * shum_val / (this%epsilon_water + (1.0D0 - this%epsilon_water) * shum_val)
          pdry_val = pres_val - pwvp_val

          ! cal refractivity
          refrac_val = this%k1 * (pdry_val / temp_val) + this%k2 * (pwvp_val / (temp_val**2.0D0)) + this%k3 * (pwvp_val / (temp_val))

          refrac%values(k) = refrac_val

        END DO

      END ASSOCIATE
    END ASSOCIATE
  END SUBROUTINE transFwdNonLinear

  FUNCTION transFwdNonLinear_opr(this, X) RESULT(Y)
    IMPLICIT NONE
    CLASS(OprGNSSRefrac_t) :: this
    TYPE(State_t), INTENT(in) :: X
    TYPE(ObsSet_t) :: Y
    INTEGER(i_kind) :: i, j, k

    Y = this%Y%zeroCopy()
    CALL this%transFwdNonLinear(X, Y)
  END FUNCTION

  SUBROUTINE transAdjMultiply(this, D, dX, X)
    IMPLICIT NONE
    !! This multiplies D by the adjoint operator to dx:
    !! Assume the Jacobian is L, dX = L^T D:
    CLASS(OprGNSSRefrac_t) :: this
    TYPE(ObsSet_t), INTENT(IN) :: D
    TYPE(State_t), INTENT(INOUT) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    INTEGER(i_kind) :: refrac_ad_ind, dpres_ind, dtemp_ind, dshum_ind
    INTEGER(i_kind) :: pres_ind, temp_ind, shum_ind, k, pres_ind_bkg, shum_ind_bkg

    REAL(r_kind) :: pres_val, shum_val, temp_val, refrac_val, refrac_ad_val, pdry_ad_val, pwvp_ad_val
    REAL(r_kind) :: pdry_val, pwvp_val, temp_ad_val, pres_ad_val, shum_ad_val

    ! This is a nonlinear operator and the current state is required:
    IF (.NOT. PRESENT(X)) THEN
      WRITE (*, 2)
2     FORMAT('OprGNSSRefrac - transAdjMultiply requires a state but it is missing, check and rerun!')
      STOP
    END IF

    refrac_ad_ind = D%getObsIdx('refractivity')

    dpres_ind = dX%getVarIdx('pres')
    dtemp_ind = dX%getVarIdx('temp')
    dshum_ind = dX%getVarIdx('qvapor')



    pres_ind = X%getVarIdx('pres')
    temp_ind = X%getVarIdx('temp')
    shum_ind = X%getVarIdx('qvapor')


    pres_ind_bkg = this%Xb%getVarIdx('pres')
    shum_ind_bkg = this%Xb%getVarIdx('qvapor')

    ! First clean the dX halo 
    CALL dX%cleanHalo()
    ! Then rev sum halo, so that the halos are filled by innercell value
    CALL dX%exHaloRevSum()

    IF ((refrac_ad_ind .EQ. 0) .OR. &
        (pres_ind .EQ. 0) .OR. &
        (temp_ind .EQ. 0) .OR. &
        (shum_ind .EQ. 0)) THEN

      WRITE (*, 1)
1     FORMAT('Error use the OprGNSSRefrac_t in transAdjMultiply', /, &
             ' check your yaml for pres temp qvapor for model state and refractivity!')
      STOP
    END IF

    ASSOCIATE (pres => X%Fields(pres_ind), &
               dpres => dX%Fields(dpres_ind), &
               temp => X%Fields(temp_ind), &
               dtemp => dX%Fields(dtemp_ind), &
               shum => X%Fields(shum_ind), &
               dshum => dX%Fields(dshum_ind) &
               )

      ASSOCIATE (refrac_ad => D%ObsFields(refrac_ad_ind))

        DO k = LBOUND(refrac_ad%idx, 1), UBOUND(refrac_ad%idx, 1)

          pres_ad_val = 0.0D0
          temp_ad_val = 0.0D0
          shum_ad_val = 0.0D0
          pwvp_ad_val = 0.0D0
          pdry_ad_val = 0.0D0
          refrac_ad_val = 0.0D0

          IF (this%gnss_adj_mode == 1) THEN        
            shum_val = this%Xb%Fields(shum_ind_bkg)%Get_Value(refrac_ad%idx(k))
            pres_val = this%Xb%Fields(pres_ind_bkg)%Get_Value(refrac_ad%idx(k)) 
            temp_val = temp%Get_Value(refrac_ad%idx(k))
          ELSE IF (this%gnss_adj_mode == 2) THEN   
            shum_val = shum%Get_Value(refrac_ad%idx(k))
            temp_val = temp%Get_Value(refrac_ad%idx(k))
            pres_val = this%Xb%Fields(pres_ind_bkg)%Get_Value(refrac_ad%idx(k)) 
          ELSE
            shum_val = shum%Get_Value(refrac_ad%idx(k))
            pres_val = pres%Get_Value(refrac_ad%idx(k)) 
            temp_val = temp%Get_Value(refrac_ad%idx(k))
          END IF

          ! shum_val = shum%Get_Value(refrac_ad%idx(k))
          ! pres_val = pres%Get_Value(refrac_ad%idx(k)) 
          ! temp_val = temp%Get_Value(refrac_ad%idx(k))

          ! PRINT*, 'HYJ+++, TEST CTL', pres_val, temp_val, shum_val

          pwvp_val = pres_val * shum_val / (this%epsilon_water + (1.0D0 - this%epsilon_water) * shum_val)
          pdry_val = pres_val - pwvp_val

          ! cal refractivity
          refrac_val = this%k1 * (pdry_val / temp_val) + this%k2 * (pwvp_val / (temp_val**2.0D0)) + this%k3 * (pwvp_val / temp_val)

          !! get obs adjoint
          refrac_ad_val = refrac_ad%values(k)

          pdry_ad_val = pdry_ad_val + refrac_ad_val * this%k1 / temp_val
          pwvp_ad_val = pwvp_ad_val + refrac_ad_val * (this%k2 / (temp_val**2.0D0) + this%k3 / temp_val)

          temp_ad_val = temp_ad_val - refrac_ad_val * (this%k1 * pdry_val / (temp_val**2.0D0) + 2 * this%k2 * pwvp_val / (temp_val**3) + this%k3 * pwvp_val / (temp_val**2.0D0))

          refrac_ad_val = 0

          pres_ad_val = pres_ad_val + pdry_ad_val
          pwvp_ad_val = pwvp_ad_val - pdry_ad_val

          pdry_ad_val = 0

          pres_ad_val = pres_ad_val + pwvp_ad_val * (shum_val / (this%epsilon_water + (1.0D0 - this%epsilon_water) * shum_val))

          ! PRINT*, shum_ad_val, shum_val, pwvp_ad_val, pres_val, pres_val / (this%epsilon_water + (1.0D0 - this%epsilon_water) * shum_val), pres_val * shum_val * (1.0D0 - this%epsilon_water) / ((this%epsilon_water + (1.0D0 - this%epsilon_water) * shum_val) * &
          ! (this%epsilon_water + (1.0D0 - this%epsilon_water) * shum_val))

          shum_ad_val = shum_ad_val + pwvp_ad_val * (pres_val  / (this%epsilon_water + (1.0D0 - this%epsilon_water) * shum_val) - &
                                                     pres_val  * shum_val * (1.0D0 - this%epsilon_water) / ((this%epsilon_water + (1.0D0 - this%epsilon_water) * shum_val) * &
                                                                                                           (this%epsilon_water + (1.0D0 - this%epsilon_water) * shum_val)))

          ! shum_ad_val = shum_ad_val 
          ! print *, 'HYJ+++: AD : pres_ad_val, temp_ad_val,shum_ad_val', pres_ad_val,temp_ad_val,shum_ad_val

          IF (this%gnss_adj_mode == 1) THEN                                                                                  
            ! CALL dpres%Add_Value(refrac_ad%idx(k), 0.0D0)
            CALL dtemp%Add_Value(refrac_ad%idx(k), temp_ad_val)
            ! CALL dshum%Add_Value(refrac_ad%idx(k), 0.0D0)
          ELSE IF (this%gnss_adj_mode == 2) THEN 
            CALL dtemp%Add_Value(refrac_ad%idx(k), temp_ad_val)
            CALL dshum%Add_Value(refrac_ad%idx(k), shum_ad_val)
          ELSE 
            CALL dpres%Add_Value(refrac_ad%idx(k), pres_ad_val)
            CALL dtemp%Add_Value(refrac_ad%idx(k), temp_ad_val)
            CALL dshum%Add_Value(refrac_ad%idx(k), shum_ad_val)
          END IF 
        END DO
      END ASSOCIATE

    END ASSOCIATE

    CALL dX%exHalo()

  END SUBROUTINE transAdjMultiply

  FUNCTION transAdjMultiply_opr(this, D, X) RESULT(dX)
    IMPLICIT NONE
    CLASS(OprGNSSRefrac_t) :: this
    TYPE(ObsSet_t), INTENT(IN) :: D
    TYPE(State_t) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    dX = X%zeroCopy()
    CALL this%transAdjMultiply(D, dX, X)
  END FUNCTION

  SUBROUTINE transFwdTanLinear(this, dX, dY, X)
    IMPLICIT NONE
    CLASS(OprGNSSRefrac_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(ObsSet_t), INTENT(INOUT) :: dY
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X

    ! local variables
    INTEGER(i_kind) :: dpres_ind, dtemp_ind, dshum_ind, pres_ind, temp_ind, shum_ind, drefrac_ind, k
    REAL(r_kind) :: pwvp_val, dpwvp_val, pdry_val, dpdry_val, refrac_val, drefrac_val, pwvp_tl_val2, refrac_tl_val2
    REAL(r_kind) :: pdry_tl_val, pwvp_tl_val, dpres_val, dtemp_val, dshum_val, pres_val, temp_val, shum_val, refrac_tl_val

    TYPE(ObsSet_t) :: dZ
    ! For a nonlinear operator, X is required:
    IF (.NOT. PRESENT(X)) THEN
      WRITE (*, 2)
2     FORMAT('OprGNSSRefrac - transFwdTanLinear requires a state but it is missing, check and rerun!')
      STOP
    END IF

    dpres_ind = dX%getVarIdx('pres')
    dtemp_ind = dX%getVarIdx('temp')
    dshum_ind = dX%getVarIdx('qvapor')


    pres_ind = X%getVarIdx('pres')
    temp_ind = X%getVarIdx('temp')
    shum_ind = X%getVarIdx('qvapor')

    drefrac_ind = dY%getObsIdx('refractivity')


    ! CALL this%transFwdNonLinear(X, this%Y)

    ASSOCIATE (dpres => dX%Fields(dpres_ind), dtemp => dX%Fields(dtemp_ind), &
               dshum => dX%Fields(dshum_ind), pres => X%Fields(pres_ind), &
               temp => X%Fields(temp_ind), shum => X%Fields(shum_ind), &
               drefrac => dY%ObsFields(drefrac_ind))

      DO k = LBOUND(drefrac%idx, 1), UBOUND(drefrac%idx, 1)

        !! todoL check the idx k
        dpres_val = dpres%Get_Value(drefrac%idx(k))
        dtemp_val = dtemp%Get_Value(drefrac%idx(k))
        dshum_val = dshum%Get_Value(drefrac%idx(k))


        pres_val = pres%Get_Value(drefrac%idx(k))
        shum_val = shum%Get_Value(drefrac%idx(k))
        temp_val = temp%Get_Value(drefrac%idx(k))

        ! print*, 'HYJ+++: debug tl'
        ! PRINT*, 'pres_val: ', pres_val
        ! PRINT*, 'temp_val: ', temp_val
        ! PRINT*, 'shum_val: ', shum_val
        ! PRINT*, 'dpres_val: ', dpres_val
        ! PRINT*, 'dtemp_val: ', dtemp_val
        ! PRINT*, 'dshum_val: ', dshum_val
        ! PRINT*, 'drefrac_val: ', drefrac_val
        ! print*, 'HYJ+++: *******'

        IF (this%gnss_adj_mode == 1) THEN       
          
          pwvp_val = pres_val * shum_val / (this%epsilon_water + (1.0D0 - this%epsilon_water) * shum_val)
          pdry_val = pres_val - pwvp_val

          refrac_tl_val = -(this%k1 * pdry_val / (temp_val ** 2.0D0) &
                    + 2.0D0 * this%k2 * pwvp_val / (temp_val**3.0D0) & 
                    + this%k3 * pwvp_val / (temp_val**2.0D0)) * dtemp_val

        ELSE IF (this%gnss_adj_mode == 2) THEN       
          ! pwvp_val = pres_val * shum_val / (this%epsilon_water + (1.0D0 - this%epsilon_water) * shum_val)

  
          pwvp_tl_val = pres_val * this%epsilon_water * dshum_val / ((this%epsilon_water + (1.0D0 - this%epsilon_water) * shum_val) ** 2)

          pdry_tl_val = -pwvp_tl_val

          refrac_tl_val = this%k1 * pdry_tl_val / temp_val + &
                          this%k2 * pwvp_tl_val / temp_val**2.0D0 + &
                          this%k3 * pwvp_tl_val / temp_val - &
                          (this%k1 * pdry_val / temp_val**2.0D0 + &
                          2.0D0 * this%k2 * pwvp_val / temp_val**3.0D0 + &
                          this%k3 * pwvp_val / temp_val**2.0D0) * dtemp_val
        ELSE

          pwvp_val = pres_val * shum_val / (this%epsilon_water + (1.0D0 - this%epsilon_water) * shum_val)

          pwvp_tl_val = pwvp_val * (dpres_val / pres_val + dshum_val / shum_val - dshum_val * pwvp_val * (1.0D0 - this%epsilon_water) / (pres_val * shum_val))

          pdry_val = pres_val - pwvp_val

          pdry_tl_val = dpres_val - pwvp_tl_val

          refrac_tl_val = this%k1 * pdry_tl_val / temp_val + &
                          this%k2 * pwvp_tl_val / temp_val**2.0D0 + &
                          this%k3 * pwvp_tl_val / temp_val - &
                          (this%k1 * pdry_val / temp_val**2.0D0 + &
                          2.0D0 * this%k2 * pwvp_val / temp_val**3.0D0 + &
                          this%k3 * pwvp_val / temp_val**2.0D0) * dtemp_val
        END IF 
        ! print*, 'HYJ+++: debug tl'
        ! PRINT*, 'pres_val: ', pres_val
        ! PRINT*, 'temp_val: ', temp_val
        ! PRINT*, 'shum_val: ', shum_val
        ! PRINT*, 'dpres_val: ', dpres_val
        ! PRINT*, 'dtemp_val: ', dtemp_val
        ! PRINT*, 'dshum_val: ', dshum_val
        ! PRINT*, 'pdry_val: ', pdry_val
        ! PRINT*, 'pwvp_val: ', pwvp_val
        ! PRINT*, 'pdry_tl_val: ', pdry_tl_val
        ! PRINT*, 'pwvp_tl_val: ', pwvp_tl_val
        ! PRINT*, 'refrac_tl_val: ', refrac_tl_val
        ! dZ = (this%fwdNL_opr(X + dX/2.0D0) - this%fwdNL_opr(X- dX/2.0D0))
        ! PRINT*, 'fwd :', dZ%ObsFields(drefrac_ind)%values(k)
        ! print*, 'HYJ+++: *******'
        drefrac%values(k) = refrac_tl_val
      END DO

    END ASSOCIATE

    ! dY = dY + this%Y

  END SUBROUTINE transFwdTanLinear

  FUNCTION transFwdTanLinear_opr(this, dX, X) RESULT(dY)
    CLASS(OprGNSSRefrac_t) :: this
    TYPE(State_t), INTENT(in) :: dX
    TYPE(State_t), OPTIONAL, INTENT(IN) :: X
    TYPE(ObsSet_t) :: dY

    dY = this%Y%zeroCopy()
    CALL this%transFwdTanLinear(dX, dY, X)
  END FUNCTION

END MODULE OprGNSSRefrac_m
