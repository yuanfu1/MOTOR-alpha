!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie, Jiongming Pang
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!   Modified by Zilong Qin (zilong.qin@gmail.com), 2022/3/13, @GBA-MWF, Shenzhen
!     Update the R matrix into the J function
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/10/07, @GBA-MWF, Shenzhen
!     Add a function to print out || H(x)-Y ||^2. for debugging purpose
!   Modified by Yuanfu Xie (yuanfu_xie@yahoo.com), 2022/10/16, @GBA-MWF, Shenzhen
!     Add a geostrophic weak constraint
!   Modified by Zilong Qin, add Jc of InCompressiable weak constrain.
!   Modified by Jiongming Pang (pang.j.m@hotmail.com), 2023/05/18, @GBA-MWF, Shenzhen
!     Add EnLoc as weak constraint, like hybrid
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE JFunc_m
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE C2O_m, ONLY: C2O_t
  USE BMatrix_m, ONLY: BMatrix_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE RMatrix_m, ONLY: RMatrix_t
  USE UV2Divg_m, ONLY: UV2Divg_t
  USE YAMLRead_m, ONLY: yaml_get_var
  IMPLICIT NONE

! #define TRACE_PRESSURE
! C1 : InCompressiable weak constrain

  TYPE JFunc_t
    TYPE(State_t), POINTER :: X
    TYPE(ObsSet_t), POINTER :: Y
    TYPE(ObsSet_t) :: INV
    TYPE(C2O_t), POINTER :: H
    TYPE(BMatrix_t), POINTER :: B
    TYPE(BMatrix_t), POINTER :: B_e
    TYPE(RMatrix_t), POINTER :: R
    TYPE(SingleGrid_t), POINTER :: sg
    TYPE(State_t) :: Xb, Xb_inc
    TYPE(State_t), DIMENSION(:), ALLOCATABLE :: Xb_enperts
    LOGICAL :: Use_JcTerm = .FALSE.
    REAL(r_kind) :: Weight_Jc = 1.0
    TYPE(UV2Divg_t) :: C1
    REAL(r_kind) :: C1Coef
    LOGICAL :: hasJc = .FALSE.
    CHARACTER(len=30) :: framework = 'FullState'
    LOGICAL :: LineSearch = .TRUE.
    CHARACTER(len=20) :: BoundaryCheck = 'inner'
    LOGICAL :: Use_JcTerm_InCompres = .FALSE.
    CHARACTER(LEN=20) :: Type_Jc_InCompres
    LOGICAL :: Hybrid
    LOGICAL :: HybInc
    ! Yuanfu Xie temporarily added a variable to indicate if this%Xb has been initialized on 2025-02-13:
    LOGICAL :: xbPresent = .FALSE.
    INTEGER(i_kind) :: ensNum
    REAL(r_kind)    :: beta_s, beta_e

  CONTAINS
    PROCEDURE, PUBLIC :: initialize
    FINAL :: destructor
    PROCEDURE, PUBLIC :: get_J0_value
    PROCEDURE, PUBLIC :: get_J_value
    PROCEDURE, PUBLIC :: get_grad_J_vec
    PROCEDURE, PUBLIC :: get_J_and_grad_J_vec
    !PROCEDURE, PUBLIC :: get_A_mat_mul             !< Hessian matrix multiply
    PROCEDURE, PUBLIC :: obsResidue
    PROCEDURE, PUBLIC :: normalized
    PROCEDURE, PUBLIC :: get_fullstate_CtlVar
  END TYPE JFunc_t

CONTAINS

  SUBROUTINE initialize(this, configFile, X, Y, H, B, R, sg, Xb, B_e)
    IMPLICIT NONE
    CLASS(JFunc_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET :: X
    TYPE(ObsSet_t), TARGET :: Y
    TYPE(C2O_t), TARGET :: H
    TYPE(BMatrix_t), TARGET :: B
    TYPE(RMatrix_t), TARGET :: R
    TYPE(SingleGrid_t), TARGET :: sg
    TYPE(State_t), OPTIONAL :: Xb
    TYPE(BMatrix_t), TARGET, OPTIONAL :: B_e

    INTEGER :: ifile
    INTEGER(i_kind) :: istatus

    istatus = yaml_get_var(TRIM(configFile), 'BMat', 'Hybrid', this%Hybrid)
    IF (istatus .NE. 0) THEN
      this%Hybrid = .FALSE.
    END IF

    istatus = yaml_get_var(TRIM(configFile), 'BMat', 'ensNum', this%ensNum)
    IF (istatus .NE. 0) THEN
      this%ensNum = 0
    END IF

    istatus = yaml_get_var(TRIM(configFile), 'BMat', 'beta_s', this%beta_s)
    IF (istatus .NE. 0) THEN
      this%beta_s = 1.0D0
    END IF

    istatus = yaml_get_var(TRIM(configFile), 'BMat', 'beta_e', this%beta_e)
    IF (istatus .NE. 0) THEN
      this%beta_e = 0.0D0
    END IF

    this%X => X
    this%Y => Y
    this%H => H
    this%B => B
    this%R => R
    this%sg => sg

    ! Check from the configfile whether to use the imcompressible weak constrain
    istatus = yaml_get_var(configFile, 'Minimization', 'Use_JcTerm_InCompres', this%Use_JcTerm_InCompres)
    PRINT *, 'Use_JcTerm_InCompres read status: ', istatus, this%Use_JcTerm_InCompres
    IF (istatus .NE. 0) this%Use_JcTerm_InCompres = .FALSE.  ! If users do not specify fill option, it does not fill

    IF (this%Use_JcTerm_InCompres) THEN
      this%C1 = UV2Divg_t(configFile, X)
      istatus = yaml_get_var(configFile, 'Minimization', 'Weight_Jc_InCompres', this%C1Coef)
      IF (istatus .NE. 0) this%C1Coef = 10.0  ! If users do not specify fill option, it does not fill
      this%C1Coef = this%C1Coef * 1.0**(sg%gLevel - 4) * 1.0D8
      istatus = yaml_get_var(configFile, 'Minimization', 'Type_Jc_InCompres', this%Type_Jc_InCompres)
    END IF

    ! Yuanfu Xie added a check to ensure Xb present on 2025-02-13:
    IF (PRESENT(Xb)) THEN
      this%Xb = Xb
      this%xbPresent = .TRUE.
    END IF

    ! Jiongming Pang added: 1) EnLoc as weak constraint like Hybrid; 2) HybInc (Hybrid in incremental analysis) referred to GSI
    IF (PRESENT(B_e)) THEN
      this%B_e => B_e
    END IF

    ifile = yaml_get_var(TRIM(configFile), 'Minimization', 'Use_JcTerm', this%Use_JcTerm)
    PRINT *, 'Reading Use_JcTerm: ', ifile, this%Use_JcTerm, X%sg%mpddInfo_sg%myrank
    ifile = yaml_get_var(TRIM(configFile), 'Minimization', 'Weight_Jc', this%Weight_Jc)
    ifile = yaml_get_var(TRIM(configFile), 'RunMode', 'Framework', this%framework)
    ifile = yaml_get_var(TRIM(configFile), 'RunMode', 'LineSearch', this%LineSearch)
    ifile = yaml_get_var(TRIM(configFile), 'RunMode', 'BoundaryCheck', this%BoundaryCheck)

    PRINT *, 'Calculate innovation vectors for the incremental scheme'

    ! Yuanfu Xie added this check to handle Xb is missing on 2025-02-13:
    IF (PRESENT(Xb)) THEN
      this%INV = this%H%fwdNL_opr(this%Xb) - this%Y
      IF (ifile .NE. 0) THEN
        this%Use_JcTerm = .FALSE.
        this%Weight_Jc = 1.0
      END IF
    ELSE
      this%INV = this%Y
      this%INV = this%INV%zeroCopy()
      this%INV = this%INV - this%Y  ! Yuanfu Xie added to retreat background state is zero if it is supplied on 2025-02-13
    END IF

  END SUBROUTINE

  SUBROUTINE get_fullstate_CtlVar(this, X_Inc2FS)
    IMPLICIT NONE
    CLASS(JFunc_t) :: this
    TYPE(State_t), INTENT(INOUT) :: X_Inc2FS

    IF (TRIM(this%framework) .EQ. 'Incremental') THEN
      X_Inc2FS = X_Inc2FS + this%Xb
    END IF

  END SUBROUTINE

  SUBROUTINE get_J0_value(this, d0, X, J0)
    CLASS(JFunc_t) :: this
    REAL(r_kind), INTENT(INOUT) :: J0
    TYPE(State_t), INTENT(IN) :: d0
    TYPE(State_t)  :: XX, X, XX_e, CC
    TYPE(ObsSet_t) :: D
    REAL(r_kind)   :: Jo

    TYPE(State_t), ALLOCATABLE :: XX_einc(:), sigmaD_einc(:)
    TYPE(ObsSet_t), ALLOCATABLE :: D_e(:)
    TYPE(State_t)   :: Dsg_tmp, IXXeinc_tmp, d0e
    INTEGER(i_kind) :: i, t

    ! Update function value
    XX = this%B%sqrt_inverse_multiply_tl(d0, X)

    ! Jiongming Pang added: 1) EnLoc as weak constraint like Hybrid; 2) HybInc (Hybrid in incremental analysis) referred to GSI
    IF (this%Hybrid) THEN
      XX_e = this%B_e.SQRTINVMUL.d0
    END IF

    ! Jiongming Pang added HybInc (Hybrid in incremental analysis) referred to GSI

    SELECT CASE (TRIM(this%framework))
    CASE ('FullState')
      D = this%R.SQRTINVMUL. (this%H%fwdTL_opr(d0, X))
    CASE ('Incremental')
      D = this%R.SQRTINVMUL. (this%H%fwdTL_opr(d0, this%Xb))
    END SELECT

    J0 = (XX.DOT.XX)
    Jo = (D.DOT.D)

    ! Jiongming Pang added: 1) EnLoc as weak constraint like Hybrid; 2) HybInc (Hybrid in incremental analysis) referred to GSI
    IF (this%Hybrid) THEN
      J0 = J0 * this%beta_s + (XX_e.DOT.XX_e) * this%beta_e
    END IF

    J0 = J0 + Jo

    IF (this%Use_JcTerm) THEN
      J0 = J0 + (d0.DOT.d0) * this%Weight_Jc
    END IF

    ! Has Jc, weak constrain - Incompressible
    IF (this%Use_JcTerm_InCompres) THEN

      IF (this%Type_Jc_InCompres == "UV2W-UVW2Divg") THEN
        CC = this%C1%fwdNL_opr(this%H%UV2W%fwdNL_opr(d0))
      ELSE IF (this%Type_Jc_InCompres == "UV2Divg") THEN
        CC = this%C1%fwdNL_opr(d0)
      ELSE
        PRINT *, 'Error use the transFwdNonLinear in UV2Divg_t, STOP! '
        STOP
      END IF

      J0 = J0 + (CC.DOT.CC) * this%C1Coef
    END IF

    ! Yuanfu Xie added Geostrophic balance:
    ! Note: 1. This is used to calculate d^T A d for FRCG lambda
    !       2. Currently implementation using u and v as control only
    IF (this%H%GeosBal%weight .GT. 0.0D0 .AND. this%H%GeosBal%iu .GT. 0) &
      J0 = J0 + this%H%GeosBal%weight * (( &
                                         d0%fields(this%H%GeosBal%iu) .DOT. &
                                         d0%fields(this%H%GeosBal%iu)) + ( &
                                         d0%fields(this%H%GeosBal%iu) .DOT. &
                                         d0%fields(this%H%GeosBal%iu)))

  END SUBROUTINE get_J0_value

  SUBROUTINE get_J_value(this, X, J)
    CLASS(JFunc_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    REAL(r_kind), INTENT(INOUT) :: J
    TYPE(State_t) :: XX, YY, CC, XX1, XX_e, Xe
    TYPE(ObsSet_t) :: DD
    INTEGER(r_kind) :: i, t

    ! Update function value
    SELECT CASE (TRIM(this%framework))
    CASE ('FullState')
      XX = this%B.SQRTINVMUL. (X - this%Xb)
      DD = this%R.SQRTINVMUL. ((this%H%fwdNL_opr(X)) - this%Y)
    CASE ('Incremental')
      XX = this%B.SQRTINVMUL.X
      DD = this%R.SQRTINVMUL. (this%INV + this%H%fwdTL_opr(X, this%Xb))
    END SELECT

    ! XX = this%B.SQRTINVMUL. (X - this%Xb)
    ! DD = this%R.SQRTINVMUL. ((this%H%fwdNL_opr(X)) - this%Y)

    ! Jiongming Pang added HybInc (RF localization refer to GSI) and Hybrid (EnLoc)
    IF (this%Hybrid) THEN
      Xe = X
      XX_e = this%B_e.SQRTINVMUL. (Xe - this%Xb)
      J = (XX.DOT.XX) * this%beta_s + (XX_e.DOT.XX_e) * this%beta_e + (DD.DOT.DD)
    ELSE
      J = (XX.DOT.XX) + (DD.DOT.DD)
    END IF

    IF (this%Use_JcTerm) THEN
      SELECT CASE (TRIM(this%framework))
      CASE ('FullState')
        XX1 = X - this%Xb
      CASE ('Incremental')
        XX1 = X
      END SELECT
      J = J + (XX1.DOT.XX1) * this%Weight_Jc
    END IF

    ! Has Jc, weak constrain
    IF (this%Use_JcTerm_InCompres) THEN

      IF (this%Type_Jc_InCompres == "UV2W-UVW2Divg") THEN
        CC = this%C1%fwdNL_opr(this%H%UV2W%fwdNL_opr(X - this%Xb))
      ELSE IF (this%Type_Jc_InCompres == "UV2Divg") THEN
        CC = this%C1%fwdNL_opr(X - this%Xb)
      ELSE
        PRINT *, 'Error use the transFwdNonLinear in UV2Divg_t, STOP! '
        STOP
      END IF

      J = J + (CC.DOT.CC) * this%C1Coef
    END IF

    ! Yuanfu Xie added Geostrophic balance:
    ! Note: 1. This is used to calculate d^T A d for FRCG lambda
    !       2. Currently implementation using u and v as control only
    IF (this%H%GeosBal%weight .GT. 0.0D0 .AND. this%H%GeosBal%iu .GT. 0) THEN
      XX = X
      CALL this%H%Ctl2State%transFwdNonLinear(XX)
      CALL this%H%GeosBal%dJcDUV(XX, YY, 0)

      J = J + 0.5D0 * this%H%GeosBal%weight * &
          ((YY%fields(this%H%GeosBal%iv) .DOT. &
            YY%fields(this%H%GeosBal%iv)) + &
           (YY%fields(this%H%GeosBal%iu) .DOT. &
            YY%fields(this%H%GeosBal%iu)))
    END IF

  END SUBROUTINE get_J_value

  SUBROUTINE get_grad_J_vec(this, X, grad_J)
    CLASS(JFunc_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t) :: XX, GG, CC, XX1, XX_e, X_e, GG_s, GG_e
    TYPE(ObsSet_t) :: DD
    TYPE(State_t), INTENT(INOUT) :: grad_J

    SELECT CASE (TRIM(this%framework))
    CASE ('FullState')
      XX = this%B.SQRTINVMUL. (X - this%Xb)
      DD = this%R.SQRTINVMUL. ((this%H%fwdNL_opr(X)) - this%Y)
    CASE ('Incremental')
      XX = this%B.SQRTINVMUL.X
      DD = this%R.SQRTINVMUL. (this%INV + this%H%fwdTL_opr(X, this%Xb))
    END SELECT

    ! Jiongming Pang added EnLoc as weak constraint like Hybrid
    IF (this%Hybrid) THEN
      X_e = X
    END IF

    GG_s = this%B%sqrt_inverse_multiply_adjoint(XX, X)

    ! Jiongming Pang added EnLoc as weak constraint like Hybrid
    IF (this%Hybrid) THEN
      XX_e = this%B_e.SQRTINVMUL. (X_e - this%Xb)
      GG_e = this%B_e%sqrt_inverse_multiply_adjoint(XX_e, X_e)
      GG = GG_s * this%beta_s + GG_e * this%beta_e
    ELSE
      GG = GG_s
    END IF

    SELECT CASE (TRIM(this%framework))
    CASE ('FullState')
      grad_J = GG * 2.0D0 &
               + (this%H%adjMul_opr(this%R.SQRTINVMULADJ.DD, X)) * 2.0D0
    CASE ('Incremental')
      grad_J = GG * 2.0D0 &
               + (this%H%adjMul_opr(this%R.SQRTINVMULADJ.DD, this%Xb)) * 2.0D0
    END SELECT

    IF (this%Use_JcTerm) THEN
      SELECT CASE (TRIM(this%framework))
      CASE ('FullState')
        XX1 = X - this%Xb
      CASE ('Incremental')
        XX1 = X
      END SELECT
      grad_J = grad_J + XX1 * 2.0D0 * this%Weight_Jc
    END IF

    ! Has Jc, weak constrain
    IF (this%Use_JcTerm_InCompres) THEN

      IF (this%Type_Jc_InCompres == "UV2W-UVW2Divg") THEN
        CC = this%C1%fwdNL_opr(this%H%UV2W%fwdNL_opr(X - this%Xb))
        grad_J = grad_J + this%H%UV2W%adjMul_opr(this%C1%adjMul_opr(CC)) * this%C1Coef * 2.0D0
      ELSE IF (this%Type_Jc_InCompres == "UV2Divg") THEN
        CC = this%C1%fwdNL_opr(X - this%Xb)
        grad_J = grad_J + this%C1%adjMul_opr(CC) * this%C1Coef * 2.0D0
      ELSE
        PRINT *, 'Error use the transFwdNonLinear in UV2Divg_t, STOP! '
        STOP
      END IF

    END IF

    ! Yuanfu Xie added Geostrophic balance:
    ! Note: 1. This is used to calculate d^T A d for FRCG lambda
    !       2. Currently implementation using u and v as control only
    IF (this%H%GeosBal%weight .GT. 0.0D0 .AND. this%H%GeosBal%iu .GT. 0) THEN
      XX = X
      CALL this%H%Ctl2State%transFwdNonLinear(XX)
      CALL this%H%GeosBal%dJcDUV(XX, GG, 0)
      grad_J = grad_J + GG
    END IF
  END SUBROUTINE get_grad_J_vec

  SUBROUTINE get_J_and_grad_J_vec(this, X, J, grad_J)
    CLASS(JFunc_t) :: this
    TYPE(State_t), INTENT(IN)   :: X
    REAL(r_kind), INTENT(INOUT)    :: J
    TYPE(State_t), INTENT(INOUT)   :: grad_J

    TYPE(State_t)     :: XX, XX1, grad_Jc, grad_Jb, grad_Jo, XTemp, XX_e, X_e
    ! TYPE(State_t) :: XmXb
    TYPE(ObsSet_t)    :: DD, D, YY
    REAL(r_kind)      :: Jb, Jo, Jc
    TYPE(State_t)     :: GG, HH, CC, GG_s, GG_e
    INTEGER(i_kind)   :: k, i, t

    ! REAL(r_kind) :: sumFieldVal, sumFieldVal_loc

    SELECT CASE (TRIM(this%framework))
    CASE ('FullState')
      ! PRINT *, 'get_J_and_grad_J_vec on the FullState framework'
      XX = this%B.SQRTINVMUL. (X - this%Xb)
      DD = this%R.SQRTINVMUL. ((this%H%fwdNL_opr(X)) - this%Y)
    CASE ('Incremental')
      ! PRINT *, 'get_J_and_grad_J_vec on the Incremental framework'
      XX = this%B.SQRTINVMUL.X
      DD = this%R.SQRTINVMUL. (this%INV + this%H%fwdTL_opr(X, this%Xb))
    END SELECT

#ifdef TRACE_PRESSURE
    PRINT *, 'XX is computed...', this%Hybrid
#endif

    ! Jiongming Pang added: 1) EnLoc as weak constraint like Hybrid; 2) HybInc (Hybrid in incremental analysis) referred to GSI
    IF (this%Hybrid) THEN
      X_e = X
    END IF

    IF (this%Hybrid) THEN
      XX = this%B.SQRTINVMUL. (X - this%Xb)
      XX_e = this%B_e.SQRTINVMUL. (X_e - this%Xb)
      Jb = (XX.DOT.XX) * this%beta_s + (XX_e.DOT.XX_e) * this%beta_e
    ELSE
      ! Yuanfu Xie added a check to ensure Xb is present: 2025-02-13
      IF (this%xbPresent) THEN
        XX = this%B.SQRTINVMUL. (X - this%Xb)
      ELSE
        XX = this%B.SQRTINVMUL.X
      END IF
      ! XX = this%B.SQRTINVMUL. (X - this%Xb)
      Jb = (XX.DOT.XX)
    END IF

    ! DD = this%R.SQRTINVMUL. ((this%H%fwdNL_opr(X)) - this%Y)
    Jo = (DD.DOT.DD)

    ! return
    J = Jb + Jo

#ifdef TRACE_PRESSURE
    PRINT *, 'J is computed---', Jb, Jo
#endif

    ! Jiongming Pang added: 1) EnLoc as weak constraint like Hybrid; 2) HybInc (Hybrid in incremental analysis) referred to GSI
    IF (this%Hybrid) THEN
      GG_s = this%B%sqrt_inverse_multiply_adjoint(XX, X)
      GG_e = this%B_e%sqrt_inverse_multiply_adjoint(XX_e, X_e)
      GG = GG_s * this%beta_s + GG_e * this%beta_e
    ELSE
      GG = this%B%sqrt_inverse_multiply_adjoint(XX, X)
    END IF

#ifdef TRACE_PRESSURE
    PRINT *, 'GG is computed', (GG.DOT.GG)
#endif

    YY = this%R.SQRTINVMULADJ.DD

#ifdef TRACE_PRESSURE
    PRINT *, 'YY is computed'
#endif

    ! BLOCK
    !   TYPE(State_t) :: GGo, Xfs
    !   TYPE(ObsSet_t) :: BB
    !   INTEGER :: iji=7, tIdx, vIdx, hIdx, iobs, ih, it

    !   SELECT CASE (TRIM(this%framework))
    !   CASE ('FullState')
    !     GGo = this%H%adjMul_opr(YY, X)
    !     Xfs = X
    !   CASE ('Incremental')
    !     GGo = this%H%adjMul_opr(YY, this%Xb)
    !     Xfs = X
    !     CALL this%get_fullstate_CtlVar(Xfs)
    !   END SELECT

    !   BB = this%H%fwdNL_opr(X)

    !   ! iobs = 4220
    !   ! ih = this%Y%ObsFields(2)%idx(iobs)%hIdx
    !   ! it = this%Y%ObsFields(2)%idx(iobs)%tIdx
    !   ! PRINT *,  'O: ', iobs, this%Y%ObsFields(2)%values(iobs), 'B: ', BB%ObsFields(2)%values(iobs)
    !   ! DO vIdx = 20,40
    !   !   PRINT *,  'BKG: ',Xfs%Fields(iji)%data(vIdx, ih, it), '        GRAD: ', GGo%Fields(iji)%data(vIdx, ih, it)
    !   ! END DO

    !   DO vIdx=1, X%sg%vLevel
    !     DO hIdx=1,X%sg%num_cell
    !       ! IF (Xfs%Fields(iji)%data(vIdx, hIdx, X%sg%tSlots) < 1.0D-16) THEN
    !       IF (Xfs%Fields(iji)%data(vIdx, hIdx, X%sg%tSlots) < 1.0D-16 .AND. ABS(GGo%Fields(iji)%data(vIdx, hIdx, X%sg%tSlots)) > 0.0D0) THEN
    !         ! PRINT *,  'BKG: ',Xfs%Fields(iji)%data(vIdx, hIdx, X%sg%tSlots), 'GAD:', GGo%Fields(iji)%data(vIdx, hIdx, X%sg%tSlots)
    !         DO iobs = 1, SIZE(DD%obsFields(2)%idx)
    !           IF (abs(DD%obsFields(2)%idx(iobs)%hIdx-hIdx)<1 .AND. abs(DD%obsFields(2)%idx(iobs)%tIdx-X%sg%tSlots)<1) THEN
    !             ! PRINT *, 'iobs = ', DD%obsFields(2)%idx(iobs)%hIdx,hIdx,DD%obsFields(2)%idx(iobs)%vIdx,vIdx,DD%obsFields(2)%idx(iobs)%tIdx,X%sg%tSlots
    !             PRINT *,  'BKG: ',X%Fields(iji)%data(vIdx, hIdx, X%sg%tSlots), 'GAD:', GGo%Fields(iji)%data(vIdx, hIdx, X%sg%tSlots)
    !             PRINT *,  'O: ', iobs, this%Y%ObsFields(2)%values(iobs), 'B: ', BB%ObsFields(2)%values(iobs)
    !           END IF
    !         END DO
    !       END IF
    !     END DO
    !   END DO

    ! END BLOCK

    SELECT CASE (TRIM(this%framework))
    CASE ('FullState')
      HH = this%H%adjMul_opr(YY, X)
    CASE ('Incremental')
      HH = this%H%adjMul_opr(YY, this%Xb)
    END SELECT

#ifdef TRACE_PRESSURE
    PRINT *, 'HH is computed', J, this%H%GeosBal%weight
#endif

    grad_Jb = GG * 2.0D0
    grad_Jo = HH * 2.0D0
    grad_J = grad_Jb + grad_Jo

    IF (this%Use_JcTerm) THEN
      SELECT CASE (TRIM(this%framework))
      CASE ('FullState')
        XX1 = X - this%Xb
      CASE ('Incremental')
        XX1 = X
      END SELECT
      J = J + (XX1.DOT.XX1) * this%Weight_Jc
      grad_J = grad_J + XX1 * 2.0D0 * this%Weight_Jc
    END IF

    ! has Jc, weak constrain
    IF (this%Use_JcTerm_InCompres) THEN
      IF (this%Type_Jc_InCompres == "UV2W-UVW2Divg") THEN
        CC = this%C1%fwdNL_opr(this%H%UV2W%fwdNL_opr(X - this%Xb))
        Jc = (CC.DOT.CC) * this%C1Coef
        IF (this%sg%isBaseProc()) PRINT *, 'Jc-inCompress: ', Jc, Jc / this%C1Coef
        J = J + Jc
        grad_J = grad_J + this%H%UV2W%adjMul_opr(this%C1%adjMul_opr(CC)) * this%C1Coef * 2.0D0
      ELSE IF (this%Type_Jc_InCompres == "UV2Divg") THEN
        CC = this%C1%fwdNL_opr(X - this%Xb)
        Jc = (CC.DOT.CC) * this%C1Coef
        IF (this%sg%isBaseProc()) PRINT *, 'Jc-inCompress: ', Jc, Jc / this%C1Coef
        J = J + Jc
        grad_J = grad_J + this%C1%adjMul_opr(CC) * this%C1Coef * 2.0D0
      ELSE
        PRINT *, 'Error use the transFwdNonLinear in UV2Divg_t, STOP! '
        STOP
      END IF

    END IF

    IF (this%sg%isBaseProc()) THEN
      IF (this%Use_JcTerm_InCompres) THEN
        WRITE (*, 12) Jb, Jo, Jc
12      FORMAT('J-G: Jb', D16.6, ' Jo ', D16.6, ' Jc ', D16.6)
      ELSE
        WRITE (*, 10) Jb, Jo
10      FORMAT('J-G: Jb', D16.6, ' Jo ', D16.6)
      END IF
    END IF

    ! Yuanfu Xie added Geostrophic balance:
    ! Note: 1. This is used to calculate d^T A d for FRCG lambda
    !       2. Currently implementation using u and v as control only
    IF (this%H%GeosBal%weight .GT. 0.0D0 .AND. this%H%GeosBal%iu .GT. 0) THEN
      XX = X
      CALL this%H%Ctl2State%transFwdNonLinear(XX)
      CALL this%H%GeosBal%dJcDUV(XX, GG, 0)

      J = J + 0.5D0 * this%H%GeosBal%weight * &
          ((GG%fields(this%H%GeosBal%iv) .DOT. &
            GG%fields(this%H%GeosBal%iv)) + &
           (GG%fields(this%H%GeosBal%iu) .DOT. &
            GG%fields(this%H%GeosBal%iu)))

      WRITE (*, 2) MINVAL(GG%fields(this%H%GeosBal%iu)%DATA), &
        MAXVAL(GG%fields(this%H%GeosBal%iu)%DATA), &
        MINVAL(GG%fields(this%H%GeosBal%iv)%DATA), &
        MAXVAL(GG%fields(this%H%GeosBal%iv)%DATA), J, this%H%GeosBal%weight
2     FORMAT('get_J_and_grad_J_vec: uv mx/mn: ', 4D12.4, ' J: ', D12.4, ' w', D12.4)

      grad_J = grad_J + GG
    END IF

    ! ! Diagnose only
    ! BLOCK
    !   USE parameters_m
    !   INTEGER :: iloc, iv, ivar, ih, it, ilev
    !   REAL :: lat, lon, dlat, dlon
    !   CHARACTER(len=100) :: filename, stepstr

    !   dlat = X%sg%cell_cntr(1, X%sg%cell_stcl(8, 1)) - X%sg%cell_cntr(1, 1)
    !   dlon = X%sg%cell_cntr(2, X%sg%cell_stcl(6, 1)) - X%sg%cell_cntr(2, 1)
    !   dlat = dlat / degree2radian * 0.5
    !   dlon = dlon / degree2radian * 0.5

    !   !   DO iloc = 1, X%sg%num_icell
    !   !     IF (ABS(X%sg%cell_cntr(1,iloc) - 30.2) < 0.1 .AND. ABS(X%sg%cell_cntr(2,iloc) - 101.5) < 0.1) THEN
    !   !       PRINT *,
    !   IF (PRESENT(debug)) THEN
    !     DO iv = 1, SIZE(DD%ObsFields)
    !       IF (DD%ObsFields(iv)%name .EQ. 'tbb') THEN
    !         DO iloc = 1, SIZE(DD%ObsFields(iv)%values)
    !           ih = DD%ObsFields(iv)%idx(iloc)%hIdx
    !           it = DD%ObsFields(iv)%idx(iloc)%tIdx
    !           lat = DD%obsFields(iv)%ObsAttrSat%latitude(iloc) / degree2radian
    !           lon = DD%obsFields(iv)%ObsAttrSat%longitude(iloc) / degree2radian

    !           ! IF (ABS(lat - 16.2) < dlat .AND. ABS(lon - 116.5) < dlon) THEN  ! a total clear location (ocean)
    !           IF (ABS(lat - 16.84) < dlat .AND. ABS(lon - 116.5) < dlon) THEN  ! a clear location (ocean)
    !               ! PRINT *, 'check if clear: ', SUM(X%fields(X%getVarIdx('qcloud'))%data(:,ih,it)),SUM(X%fields(X%getVarIdx('qice'))%data(:,ih,it))
    !               PRINT *, 'check J at a clear1 location: ichan = ',iv, 'rank = ', X%sg%mpddInfo_sg%myrank, ' step = ', step, &
    !               ' Jo = ', DD%obsFields(iv)%values(iloc) * DD%obsFields(iv)%values(iloc), &
    !               ' grad_Jo = ', SUM(HH%fields(HH%getVarIdx('qvapor_ctl'))%data(:,ih,it)) * 2.0, SUM(HH%fields(HH%getVarIdx('qcloud'))%data(:,ih,it)) * 2.0, &
    !                             SUM(HH%fields(HH%getVarIdx('qice'))%data(:,ih,it)) * 2.0
    !               ! ' grad_Jo = ', MAXVAL(HH%fields(HH%getVarIdx('qvapor_ctl'))%data(:,ih,it)) * 2.0
    !               ! DO ivar = 1, SIZE(XX%fields)
    !               !   print *, 'xx field name: ', ivar, XX%fields(ivar)%name
    !               ! END DO
    !               PRINT *, 'check J at a clear1 location: ichan = ',iv, 'rank = ', X%sg%mpddInfo_sg%myrank, ' step = ', step, &
    !               ' Jb = ', SUM(XX%fields(XX%getVarIdx('qvapor'))%data(:,ih,it) *  &
    !                             XX%fields(XX%getVarIdx('qvapor'))%data(:,ih,it)), &
    !               ' grad_Jb = ', MAXVAL(GG%fields(GG%getVarIdx('qvapor_ctl'))%data(:,ih,it)) * 2.0
    !               IF (step < 10) THEN
    !                 WRITE(stepstr,'(i1)') step
    !               ELSE
    !                 IF (step < 100) THEN
    !                   WRITE(stepstr,'(i2)') step
    !                 ELSE
    !                   WRITE(stepstr,'(i3)') step
    !                 END IF
    !               END IF
    !               filename = '/Users/yaliwu/Desktop/MOTOR/MOTOR/output/test_dpp_cloudy/clear1_'//TRIM(stepstr)//'_gradJo.txt'
    !               OPEN(unit=10, file=filename, status='replace')
    !               DO ilev = 1, XX%sg%vLevel
    !                 WRITE(10, *) HH%fields(HH%getVarIdx('qvapor_ctl'))%data(ilev,ih,it) * 2.0, HH%fields(HH%getVarIdx('qcloud'))%data(ilev,ih,it) * 2.0, &
    !                             HH%fields(HH%getVarIdx('qice'))%data(ilev,ih,it) * 2.0
    !               END DO
    !               CLOSE(10)
    !           END IF

    !           ! IF (ABS(lat - 29.0) < dlat .AND. ABS(lon - 97.8) < dlon) THEN  ! a total clear location (mountain)
    !           IF (ABS(lat - 29.6) < dlat .AND. ABS(lon - 97.5) < dlon) THEN  ! a clear location (mountain)
    !             ! PRINT *, 'check if clear: ', SUM(X%fields(X%getVarIdx('qcloud'))%data(:,ih,it)),SUM(X%fields(X%getVarIdx('qice'))%data(:,ih,it))
    !             PRINT *, 'check J at a clear2 location: ichan = ',iv, 'rank = ', X%sg%mpddInfo_sg%myrank, ' step = ', step, &
    !             ' Jo = ', DD%obsFields(iv)%values(iloc) * DD%obsFields(iv)%values(iloc), &
    !             ' grad_Jo = ', SUM(HH%fields(HH%getVarIdx('qvapor_ctl'))%data(:,ih,it)) * 2.0, SUM(HH%fields(HH%getVarIdx('qcloud'))%data(:,ih,it)) * 2.0, &
    !                           SUM(HH%fields(HH%getVarIdx('qice'))%data(:,ih,it)) * 2.0
    !             ! ' grad_Jo = ', MAXVAL(HH%fields(HH%getVarIdx('qvapor_ctl'))%data(:,ih,it)) * 2.0
    !             ! DO ivar = 1, SIZE(XX%fields)
    !             !   print *, 'xx field name: ', ivar, XX%fields(ivar)%name
    !             ! END DO
    !             PRINT *, 'check J at a clear2 location: ichan = ',iv, 'rank = ', X%sg%mpddInfo_sg%myrank, ' step = ', step, &
    !             ' Jb = ', SUM(XX%fields(XX%getVarIdx('qvapor'))%data(:,ih,it) *  &
    !                           XX%fields(XX%getVarIdx('qvapor'))%data(:,ih,it)), &
    !             ' grad_Jb = ', MAXVAL(GG%fields(GG%getVarIdx('qvapor_ctl'))%data(:,ih,it)) * 2.0

    !             IF (step < 10) THEN
    !               WRITE(stepstr,'(i1)') step
    !             ELSE
    !               IF (step < 100) THEN
    !                 WRITE(stepstr,'(i2)') step
    !               ELSE
    !                 WRITE(stepstr,'(i3)') step
    !               END IF
    !             END IF
    !             filename = '/Users/yaliwu/Desktop/MOTOR/MOTOR/output/test_dpp_cloudy/clear2_'//TRIM(stepstr)//'_gradJo.txt'
    !             OPEN(unit=10, file=filename, status='replace')
    !             DO ilev = 1, XX%sg%vLevel
    !               WRITE(10, *) HH%fields(HH%getVarIdx('qvapor_ctl'))%data(ilev,ih,it) * 2.0, HH%fields(HH%getVarIdx('qcloud'))%data(ilev,ih,it) * 2.0, &
    !                            HH%fields(HH%getVarIdx('qice'))%data(ilev,ih,it) * 2.0
    !             END DO
    !             CLOSE(10)
    !           END IF

    !           IF (ABS(lat - 24.0) < dlat .AND. ABS(lon - 109.9) < dlon) THEN  ! a cloudy location
    !             PRINT *, 'check J at a cloudy location: ichan = ',iv, 'rank = ', X%sg%mpddInfo_sg%myrank, ' step = ', step, &
    !             ' Jo = ', DD%obsFields(iv)%values(iloc) * DD%obsFields(iv)%values(iloc), &
    !             ' grad_Jo = ', SUM(HH%fields(HH%getVarIdx('qvapor_ctl'))%data(:,ih,it)) * 2.0, SUM(HH%fields(HH%getVarIdx('qcloud'))%data(:,ih,it)) * 2.0, &
    !                           SUM(HH%fields(HH%getVarIdx('qice'))%data(:,ih,it)) * 2.0
    !             ! ' grad_Jo = ', MAXVAL(HH%fields(HH%getVarIdx('qcloud'))%data(:,ih,it)) * 2.0, &
    !             !               MAXVAL(HH%fields(HH%getVarIdx('qice'))%data(:,ih,it)) * 2.0
    !             PRINT *, 'check J at a cloudy location: ichan = ',iv, 'rank = ', X%sg%mpddInfo_sg%myrank, ' step = ', step, &
    !             ' Jb = ', SUM(GG%fields(GG%getVarIdx('qcloud'))%data(:,ih,it) * &
    !                           GG%fields(GG%getVarIdx('qcloud'))%data(:,ih,it)), SUM(GG%fields(GG%getVarIdx('qcloud'))%data(:,ih,it) * &
    !                           GG%fields(GG%getVarIdx('qcloud'))%data(:,ih,it)), &
    !             ' grad_Jb = ', MAXVAL(GG%fields(GG%getVarIdx('qcloud'))%data(:,ih,it)) * 2.0, &
    !                           MAXVAL(GG%fields(GG%getVarIdx('qice'))%data(:,ih,it)) * 2.0
    !             IF (step < 10) THEN
    !               WRITE(stepstr,'(i1)') step
    !             ELSE
    !               IF (step < 100) THEN
    !                 WRITE(stepstr,'(i2)') step
    !               ELSE
    !                 WRITE(stepstr,'(i3)') step
    !               END IF
    !             END IF
    !             filename = '/Users/yaliwu/Desktop/MOTOR/MOTOR/output/test_dpp_cloudy/cloudy_'//TRIM(stepstr)//'_gradJo.txt'
    !             OPEN(unit=10, file=filename, status='replace')
    !             DO ilev = 1, XX%sg%vLevel
    !               WRITE(10, *) HH%fields(HH%getVarIdx('qvapor_ctl'))%data(ilev,ih,it) * 2.0, HH%fields(HH%getVarIdx('qcloud'))%data(ilev,ih,it) * 2.0, &
    !                             HH%fields(HH%getVarIdx('qice'))%data(ilev,ih,it) * 2.0
    !             END DO
    !             CLOSE(10)
    !           END IF

    !         END DO
    !       END IF
    !     END DO
    !   END IF
    ! END BLOCK

  END SUBROUTINE get_J_and_grad_J_vec

  ! A = B^(-1)+H^(T)*R^(-1)*H
  ! SUBROUTINE get_A_mat_mul(this, Xi, Xo)
  !   CLASS(JFunc_t) :: this
  !   TYPE(State_t), INTENT(IN) :: Xi
  !   TYPE(State_t), INTENT(INOUT) :: Xo
  !   TYPE(State_t) :: XX

  !   ! Xo = this%B.SQRTINVMULADJ. (this%B.SQRTINVMUL.Xi)
  !   Xo = this%B%sqrt_inverse_multiply_adjoint(Xi,)
  !   Xo = Xo + (this%H%adjMul_opr(this%R.SQRTINVMULADJ. (this%R.SQRTINVMUL. (this%H%fwdNL_opr(Xi))), Xi))

  ! END SUBROUTINE get_A_mat_mul

  SUBROUTINE obsResidue(this, X, itr, MaxOptStep, costJo)
    CLASS(JFunc_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    INTEGER(i_kind), INTENT(IN) :: itr, MaxOptStep
    REAL(r_kind), INTENT(INOUT) :: costJo(MaxOptStep, 49)

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, l
    REAL(r_kind) :: totalJo
    TYPE(ObsSet_t) :: DD, RD

    DD = (this%H%fwdNL_opr(X)) - this%Y
    RD = this%R.SQRTINVMUL. ((this%H%fwdNL_opr(X)) - this%Y)
    totalJo = (RD.DOT.RD)

#ifdef TRACE_PRESSURE
    WRITE (*, 5) X%sg%glevel, itr, totalJo
5   FORMAT('TotalJo at G', I1, ' itr: ', I4, ' total: ', E14.6)
#endif

    ! For all obs:
    k = LBOUND(DD%ObsFields, 1) ! Start with the first obsFields
    j = LBOUND(DD%ObsFields, 1) - 1 ! Counter
#ifdef TRACE_PRESSURE
    WRITE (*, 2) (TRIM(DD%ObsFields(i)%Get_ObsType()), i=LBOUND(DD%ObsFields, 1), UBOUND(DD%ObsFields, 1))
2   FORMAT('OBSFields obsType: ', 20(1X, A))
#endif
    DO i = LBOUND(DD%ObsFields, 1), UBOUND(DD%ObsFields, 1)
      ! Normalizing factors:
      costJo(itr, i) = (DD%ObsFields(i) .DOT.DD%ObsFields(i))
      IF (costJo(itr, i) .LE. 0.0D0) costJo(itr, i) = 1.0D0
      IF (TRIM(DD%ObsFields(i)%Get_ObsType()) .EQ. &
          TRIM(DD%ObsFields(k)%Get_ObsType())) THEN
        j = j + 1
      ELSE
        ! Output || H(x)-Y || by components:
        WRITE (*, 1) X%sg%glevel, itr, TRIM(DD%ObsFields(k)%Get_ObsType()), &
          (TRIM(DD%ObsFields(l)%Get_Name()), &
           (DD%ObsFields(l) .DOT.DD%ObsFields(l)), &
           l=k, j)
1       FORMAT('||H(X)-Y||2 G', I1, ' itr:', I4, ' obsType ', A, 20(1X, A, E14.6))

        ! Output (H(x)-Y)^T R^{-1} (H(x)-Y) by components:
        WRITE (*, 3) X%sg%glevel, itr, TRIM(RD%ObsFields(k)%Get_ObsType()), &
          (TRIM(RD%ObsFields(l)%Get_Name()), &
           (RD%ObsFields(l) .DOT.RD%ObsFields(l)), &
           l=k, j)
3       FORMAT('Jo-item G', I1, ' itr:', I4, ' obsType ', A, ' components', 20(1X, A, E14.6))

        k = i
        j = i

      END IF
    END DO
    ! Output || H(x)-Y || by components: Last part
    WRITE (*, 1) X%sg%glevel, itr, TRIM(DD%ObsFields(k)%Get_ObsType()), &
      (TRIM(DD%ObsFields(l)%Get_Name()), &
       (DD%ObsFields(l) .DOT.DD%ObsFields(l)), &
       l=k, j)

    ! Output (H(x)-Y)^T R^{-1} (H(x)-Y) by components: last part
    WRITE (*, 3) X%sg%glevel, itr, TRIM(RD%ObsFields(k)%Get_ObsType()), &
      (TRIM(RD%ObsFields(l)%Get_Name()), &
       (RD%ObsFields(l) .DOT.RD%ObsFields(l)), &
       l=k, j)
  END SUBROUTINE obsResidue

  SUBROUTINE normalized(this, X, itr, MaxOptStep, costJo)
    CLASS(JFunc_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    INTEGER(i_kind), INTENT(IN) :: itr, MaxOptStep
    REAL(r_kind), INTENT(INOUT) :: costJo(MaxOptStep, 50)

    ! Local variables:
    INTEGER(i_kind) :: i, j, k, m, l
    REAL(r_kind) :: maxmin(50)

    ! The ranges of the costs:
    maxmin(1) = (MAXVAL(costJo(1:itr, 1)) - MINVAL(costJo(1:itr, 1)))
    DO i = LBOUND(this%Y%ObsFields, 1), UBOUND(this%Y%ObsFields, 1)
      maxmin(i + 1) = MAXVAL(costJo(1:itr, i + 1) - MINVAL(costJo(1:itr, i + 1)))
      IF (maxmin(i + 1) .LE. 0.0D0) maxmin(i + 1) = 1.0D0
      WRITE (*, 3) this%Y%ObsFields(i)%Get_ObsType(), this%Y%ObsFields(i)%Get_Name(), &
        MINVAL(costJo(1:itr, i + 1)), MAXVAL(costJo(1:itr, i + 1)), maxmin(i + 1), X%sg%glevel
3     FORMAT('ObsTYPE:', 2(1X, A), ' min.max:', 3E14.6, ' G', I1)
    END DO

    ! For iterations:
    DO m = 1, itr
      WRITE (*, 2) X%sg%glevel, m, (costJo(m, 1) - MINVAL(costJo(1:itr, 1))) / maxmin(1)
2     FORMAT('Normalized J G', I1, ' itr:', I4, ' totalJo ', E14.6)
      WRITE (*, 4) X%sg%glevel, m, costJo(m, 1)
4     FORMAT('Un-normal  J G', I1, ' itr:', I4, ' totalJo ', E14.6)
      ! For all obs:
      k = LBOUND(this%Y%ObsFields, 1) ! Start with the first obsFields
      j = LBOUND(this%Y%ObsFields, 1) - 1 ! Counter
      DO i = LBOUND(this%Y%ObsFields, 1), UBOUND(this%Y%ObsFields, 1)
        IF (TRIM(this%Y%ObsFields(i)%Get_ObsType()) .EQ. &
            TRIM(this%Y%ObsFields(k)%Get_ObsType())) THEN
          j = j + 1
        ELSE
          WRITE (*, 1) X%sg%glevel, m, TRIM(this%Y%ObsFields(k)%Get_ObsType()), &
            (TRIM(this%Y%ObsFields(l)%Get_Name()), &
             (costJo(m, l + 1) - MINVAL(costJo(1:itr, l + 1))) / maxmin(l + 1), &
             l=k, j)
1         FORMAT('Normalized Jos G', I1, ' itr:', I4, ' obsType ', A, 20(1X, A, E14.6))

          k = i
          j = i

        END IF
      END DO
      WRITE (*, 1) X%sg%glevel, m, TRIM(this%Y%ObsFields(k)%Get_ObsType()), &
        (TRIM(this%Y%ObsFields(l)%Get_Name()), &
         (costJo(m, l + 1) - MINVAL(costJo(1:itr, l + 1))) / maxmin(l + 1), &
         l=k, j)
    END DO
  END SUBROUTINE normalized

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(JFunc_t), INTENT(INOUT) :: this

    PRINT *, 'End of destructor of JFunc.'
  END SUBROUTINE destructor

END MODULE JFunc_m
