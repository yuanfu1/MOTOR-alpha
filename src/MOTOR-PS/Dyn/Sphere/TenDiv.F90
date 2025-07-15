MODULE TenDiv_m
  USE gzm_m, ONLY: gzm_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE InterpValue_m, ONLY: InterpValue_t
  USE CalVerDer_m, ONLY: CalVerDer_t
  USE PreCal_m, ONLY: PreCal_t
  USE State_m, ONLY: State_t
  USE parameters_m, ONLY: EarthRadius

  IMPLICIT NONE

  TYPE :: TenDiv_t
    TYPE(SingleGrid_t), POINTER :: sg
    TYPE(InterpValue_t) :: InterpValue
    TYPE(gzm_t) :: gzm
    TYPE(CalVerDer_t) :: CalVerDer
    ! TYPE(PreCal_t), pointer :: PreCal
  CONTAINS
    PROCEDURE :: Tendency
    PROCEDURE :: destroy
    FINAL :: destructor
  END TYPE TenDiv_t

  INTERFACE TenDiv_t
    PROCEDURE :: constructor
  END INTERFACE TenDiv_t

CONTAINS

  FUNCTION constructor(sg) RESULT(this)
    TYPE(TenDiv_t) :: this
    TYPE(SingleGrid_t), TARGET :: sg
    TYPE(PreCal_t), TARGET :: PreCal

    this%sg => sg
    ! this%PreCal => PreCal
    ! PRINT *, 'max of DM in TenDiv initial function is: ', MAXVAL(this%precal%DM), MINVAL(this%precal%DM), SIZE(this%precal%DM)

    this%gzm = gzm_t(this%sg)
    this%InterpValue = InterpValue_t(this%sg)

  END FUNCTION constructor

  SUBROUTINE Tendency(this, div, vor, psi, chi, w, rho, P, PreCal, ten_div, X1)
    IMPLICIT NONE

    CLASS(TenDiv_t), INTENT(INOUT) :: this
    TYPE(State_t), INTENT(INOUT), OPTIONAL :: X1
    TYPE(PreCal_t), INTENT(IN) :: PreCal
    REAL(r_kind), DIMENSION(:, :), INTENT(IN) :: div, vor, psi, chi, w, rho, P
    REAL(r_kind), INTENT(INOUT) :: ten_div(:, :)
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: w2, omg, &
                                                  F_vor_psi, J_vor_chi, &
                                                  J_DM_psisigma, F_DM_chisigma, &
                                                  J_omg_psisigma, F_omg_chisigma, &
                                                  F_invRho_P, F_hzrhoP_z, &
                                                  F_f_psi, J_f_chi, F_psi_psi, &
                                                  F_chi_chi, J_psi_chi, KE, LKE, J_DM_psi, &
                                                  F_DM_chi, J_omg_psi, F_omg_chi

    REAL(r_kind) :: start_time_t, end_time_t, start_interp, end_interp, &
                    start_gzm, end_gzm, &
                    start_time_al, end_time_al

    REAL(r_kind) :: r_sg
    INTEGER(i_kind) :: k

    ASSOCIATE (num_cell => this%sg%num_cell, &
               num_iCell => this%sg%num_icell, &
               vLevel => this%sg%vLevel, &
               sigma => this%sg%sigma)

      CALL CPU_TIME(start_time_t)

      ALLOCATE (w2(vLevel, num_cell), omg(vLevel, num_cell), &
                F_vor_psi(vLevel, num_cell), J_vor_chi(vLevel, num_cell), &
                KE(vLevel, num_cell), LKE(vLevel, num_cell), &
                F_psi_psi(vLevel, num_cell), F_chi_chi(vLevel, num_cell), &
                J_DM_psisigma(vLevel, num_cell), F_DM_chisigma(vLevel, num_cell), &
                J_omg_psisigma(vLevel, num_cell), F_omg_chisigma(vLevel, num_cell), &
                F_invRho_P(vLevel, num_cell), F_hzrhoP_z(vLevel, num_cell), &
                F_f_psi(vLevel, num_cell), J_f_chi(vLevel, num_cell), &
                J_DM_psi(vLevel, num_cell), F_DM_chi(vLevel, num_cell), &
                J_omg_psi(vLevel, num_cell), F_omg_chi(vLevel, num_cell), &
                J_psi_chi(vLevel, num_cell))

      w2 = 0.0D0
      omg = 0.0D0
      F_vor_psi = 0.0D0
      J_vor_chi = 0.0D0
      J_DM_psisigma = 0.0D0
      F_DM_chisigma = 0.0D0
      J_omg_psisigma = 0.0D0
      F_omg_chisigma = 0.0D0
      F_invRho_P = 0.0D0
      F_hzrhoP_z = 0.0D0
      F_f_psi = 0.0D0
      J_f_chi = 0.0D0
      F_psi_psi = 0.0D0
      F_chi_chi = 0.0D0
      J_psi_chi = 0.0D0
      KE = 0.0D0
      LKE = 0.0D0
      J_DM_psi = 0.0D0
      F_DM_chi = 0.0D0
      J_omg_psi = 0.0D0
      F_omg_chi = 0.0D0
      CALL CPU_TIME(end_time_al)
      PRINT *, 'allocate time in tendiv is: ', end_time_al - start_time_t, 'S'
      !  CALL interpvalue(sigma1, sigma2, w(:, 1:num_iCell), w2(:, 1:num_iCell))
      CALL CPU_TIME(start_interp)
      CALL this%InterpValue%InterpValue(w, w2)
      ! w2(:, 1:num_iCell) = (w(1:vLevel, num_iCell) + w(2:vLevel + 1, num_iCell))/2.0D0
      omg(:, 1:num_iCell) = this%sg%Hz(:, 1:num_iCell) * w2(:, 1:num_iCell)
      CALL CPU_TIME(end_interp)
      PRINT *, 'interp time in tendiv is: ', end_interp - start_interp, 'S'

      ! PRINT *, 'max of DM in TenDiv is: ', MAXVAL(PreCal%DM), MINVAL(PreCal%DM), SIZE(PreCal%DM)
      ! PRINT *, 'max of omg in TenDiv is: ', MAXVAL(omg), MINVAL(omg), SIZE(omg)
      CALL CPU_TIME(start_gzm)
      CALL CPU_TIME(start_time_al)
      CALL this%gzm%Divergen(vor, psi, F_vor_psi)
      CALL CPU_TIME(end_gzm)
      PRINT *, 'divergence time in tendiv is: ', end_gzm - start_gzm, 'S'
      ! IF (MAXVAL(F_vor_psi) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in F_vor_psi.', MAXVAL(psi), MAXVAL(vor)
      !    stop
      ! END IF
      CALL CPU_TIME(start_gzm)
      CALL this%gzm%Jacobian(vor, chi, J_vor_chi)
      CALL CPU_TIME(end_gzm)
      PRINT *, 'jacobian time in tendiv is: ', end_gzm - start_gzm, 'S'
      ! IF (MAXVAL(J_vor_chi) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in J_vor_chi.', MAXVAL(chi), MAXVAL(vor)
      !    stop
      ! END IF
      CALL this%gzm%Divergen(psi, psi, F_psi_psi)
      ! IF (MAXVAL(F_psi_psi) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in F_psi_psi.', MAXVAL(psi)
      !    stop
      ! END IF
      CALL this%gzm%Divergen(chi, chi, F_chi_chi)
      ! IF (MAXVAL(F_chi_chi) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in F_chi_chi.', MAXVAL(chi)
      !    stop
      ! END IF
      CALL this%gzm%Jacobian(psi, chi, J_psi_chi)
      ! IF (MAXVAL(J_psi_chi) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in J_psi_chi.', MAXVAL(psi), MAXVAL(chi)
      !    stop
      ! END IF
      CALL this%gzm%Jacobian(PreCal%DM, PreCal%par_psi_sigma, J_DM_psisigma)
      ! IF (MAXVAL(J_DM_psisigma) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in J_DM_psisigma.', MAXVAL(PreCal%DM), MAXVAL(PreCal%par_psi_sigma)
      !    stop
      ! END IF
      CALL this%gzm%Divergen(PreCal%DM, PreCal%par_chi_sigma, F_DM_chisigma)
      ! IF (MAXVAL(F_DM_chisigma) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in F_DM_chisigma.', MAXVAL(PreCal%DM), MAXVAL(PreCal%par_chi_sigma)
      !    stop
      ! END IF
      CALL this%gzm%Jacobian(omg, PreCal%par_psi_sigma, J_omg_psisigma)
      ! IF (MAXVAL(J_omg_psisigma) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in J_omg_psisigma.', MAXVAL(omg), MAXVAL(PreCal%par_psi_sigma)
      !    stop
      ! END IF
      CALL this%gzm%Divergen(omg, PreCal%par_chi_sigma, F_omg_chisigma)
      CALL this%gzm%Divergen(PreCal%invRho, P, F_invRho_P)
      ! IF (MAXVAL(F_invRho_P) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in F_invRho_P.', MAXVAL(PreCal%invRho), MAXVAL(P)
      !    stop
      ! END IF
      CALL this%gzm%Divergen(PreCal%HzInvRhoParPsigma, this%sg%zHght, F_hzrhoP_z)
      ! IF (MAXVAL(F_hzrhoP_z) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in F_hzrhoP_z.', MAXVAL(PreCal%HzInvRhoParPsigma), MAXVAL(this%sg%zHght)
      !    stop
      ! END IF
      CALL this%gzm%Divergen(this%sg%f, psi, F_f_psi)
      ! IF (MAXVAL(F_f_psi) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in F_f_psi.', MAXVAL(this%sg%f), MAXVAL(psi)
      !    stop
      ! END IF
      CALL this%gzm%Jacobian(this%sg%f, chi, J_f_chi)
      ! IF (MAXVAL(J_f_chi) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in J_f_chi.', MAXVAL(this%sg%f), MAXVAL(chi)
      !    stop
      ! END IF
      CALL this%gzm%Divergen(PreCal%DM, chi, F_DM_chi)
      ! IF (MAXVAL(F_DM_chi) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in F_DM_chi.', MAXVAL(PreCal%DM), MAXVAL(chi), MAXVAL(F_DM_chi)
      !    stop
      ! END IF
      CALL this%gzm%Jacobian(PreCal%DM, psi, J_DM_psi)
      ! IF (MAXVAL(J_DM_psi) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in J_DM_psi.', MAXVAL(PreCal%DM), MAXVAL(psi)
      !    stop
      ! END IF
      CALL this%gzm%Jacobian(omg, psi, J_omg_psi)

      CALL this%gzm%Divergen(omg, chi, F_omg_chi)
      CALL CPU_TIME(end_gzm)
      PRINT *, 'total gzm time in tendiv is: ', end_gzm - start_time_al, 'S'

      KE(:, 1:num_iCell) = (F_psi_psi(:, 1:num_iCell) - psi(:, 1:num_iCell) * vor(:, 1:num_iCell) &
                            + F_chi_chi(:, 1:num_iCell) - chi(:, 1:num_iCell) * div(:, 1:num_iCell)) / 2.0D0 &
                           + J_psi_chi(:, 1:num_iCell)
      CALL this%gzm%Laplacia(KE, LKE)

      ! Calculate Tendency of divergence
      DO k = 1, vLevel
        r_sg = EarthRadius + sigma(k)
        ten_div(k, 1:num_iCell) = F_vor_psi(k, 1:num_iCell) + J_vor_chi(k, 1:num_iCell) &
                                  - LKE(k, 1:num_iCell) &
                                  - J_DM_psisigma(k, 1:num_iCell) + F_DM_chisigma(k, 1:num_iCell) &
                                  + J_omg_psisigma(k, 1:num_iCell) - F_omg_chisigma(k, 1:num_iCell) &
                                  + 1.0D0 / r_sg * (J_DM_psi(k, 1:num_iCell) - F_DM_chi(k, 1:num_iCell) &
                                                    + PreCal%DM(k, 1:num_iCell) * div(k, 1:num_iCell) &
                                                    - J_omg_psi(k, 1:num_iCell) + F_omg_chi(k, 1:num_iCell) &
                                                    - omg(k, 1:num_iCell) * div(k, 1:num_iCell)) &
                                  - F_invRho_P(k, 1:num_iCell) + F_hzrhoP_z(k, 1:num_iCell) &
                                  + F_f_psi(k, 1:num_iCell) + J_f_chi(k, 1:num_iCell)

      END DO
    END ASSOCIATE

    ! IF (PRESENT(X1)) THEN
    !    X1%fields(X1%getVarIdx('J_vor_chi'))%data(:, :, 1) = J_vor_chi
    !    X1%fields(X1%getVarIdx('F_vor_psi'))%data(:, :, 1) = F_vor_psi
    !    X1%fields(X1%getVarIdx('LKE'))%data(:, :, 1) = LKE
    !    X1%fields(X1%getVarIdx('J_DM_psisigma'))%data(:, :, 1) = J_DM_psisigma
    !    X1%fields(X1%getVarIdx('F_DM_chisigma'))%data(:, :, 1) = F_DM_chisigma
    !    X1%fields(X1%getVarIdx('J_DM_psi'))%data(:, :, 1) = J_DM_psi
    !    X1%fields(X1%getVarIdx('F_DM_chi'))%data(:, :, 1) = F_DM_chi
    !    X1%fields(X1%getVarIdx('DM'))%data(:, :, 1) = PreCal%DM
    !    X1%fields(X1%getVarIdx('psi_sigma'))%data(:, :, 1) = PreCal%par_psi_sigma
    !    X1%fields(X1%getVarIdx('F_invRho_P'))%data(:, :, 1) = F_invRho_P
    !    X1%fields(X1%getVarIdx('F_HzinvRhoPsigma_z'))%data(:, :, 1) = F_hzrhoP_z
    !    X1%fields(X1%getVarIdx('F_f_psi'))%data(:, :, 1) = F_f_psi
    !    X1%fields(X1%getVarIdx('J_f_chi'))%data(:, :, 1) = J_f_chi
    ! END IF

    DEALLOCATE (w2, omg, &
                F_vor_psi, J_vor_chi, &
                KE, LKE, &
                F_psi_psi, F_chi_chi, &
                J_DM_psisigma, F_DM_chisigma, &
                J_omg_psisigma, F_omg_chisigma, &
                F_invRho_P, F_hzrhoP_z, &
                F_f_psi, J_f_chi, &
                J_DM_psi, F_DM_chi, &
                J_omg_psi, F_omg_chi, &
                J_psi_chi)

  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(TenDiv_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    ! IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destructor

  SUBROUTINE destroy(this)
    IMPLICIT NONE
    CLASS(TenDiv_t), INTENT(INOUT) :: this

    IF (ASSOCIATED(this%sg)) NULLIFY (this%sg)
    ! IF (ASSOCIATED(this%PreCal)) NULLIFY (this%PreCal)

  END SUBROUTINE destroy

END MODULE
