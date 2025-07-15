MODULE PreCal_m
  USE State_m, ONLY: State_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind
  USE gzm_m, ONLY: gzm_t
  USE CalVerDer_m, ONLY: CalVerDer_t

  IMPLICIT NONE

  TYPE :: PreCal_t
    TYPE(gzm_t) :: gzm
    TYPE(CalVerDer_t) :: derivatives
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: DM, par_psi_sigma, &
                                                  par_chi_sigma, invRho, &
                                                  HzInvRhoParPsigma

  CONTAINS
    PROCEDURE :: destroy
    PROCEDURE :: PreCals
    FINAL :: destructor
  END TYPE

  INTERFACE PreCal_t
    PROCEDURE :: constructor
  END INTERFACE PreCal_t

CONTAINS

  SUBROUTINE PreCals(this, X)
    IMPLICIT NONE
    CLASS(PreCal_t), INTENT(INOUT) :: this
    TYPE(State_t), INTENT(IN) :: X
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE ::par_P_sigma, &
                                                  J_z_psi, F_z_chi

    ASSOCIATE (psi => X%Fields(X%getVarIdx('psi'))%DATA(:, :, 1), &
               chi => X%Fields(X%getVarIdx('chi'))%DATA(:, :, 1), &
               rho => X%Fields(X%getVarIdx('rho'))%DATA(:, :, 1), &
               div => X%Fields(X%getVarIdx('div'))%DATA(:, :, 1), &
               pres => X%Fields(X%getVarIdx('pres'))%DATA(:, :, 1), &
               z => X%sg%zHght, &
               Hz => X%sg%Hz, &
               vLevel => X%sg%vLevel, &
               num_cell => X%sg%num_cell, &
               num_icell => X%sg%num_icell &
               )

      ALLOCATE (F_z_chi(vLevel, num_cell), J_z_psi(vLevel, num_cell), &
                par_P_sigma(vLevel, num_cell))
      F_z_chi = 0.0D0
      J_z_psi = 0.0D0
      par_P_sigma = 0.0D0

      CALL this%derivatives%FirstOrder(psi, this%par_psi_sigma)
      CALL this%derivatives%FirstOrder(chi, this%par_chi_sigma)
      CALL this%derivatives%FirstOrder(pres, par_P_sigma)
      CALL this%gzm%Jacobian(z, psi, J_z_psi)
      ! IF (MAXVAL(J_z_psi) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in J_z_psi in precal.', MAXVAL(z), MAXVAL(psi)
      !    stop
      ! END IF
      ! PRINT *, 'F_z_chi debug starts: '
      CALL this%gzm%Divergen(z, chi, F_z_chi)
      ! IF (MAXVAL(F_z_chi) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in F_z_chi in precal.', MAXVAL(z), MAXVAL(chi), MAXVAL(F_z_chi)
      !    stop
      ! END IF
      ! IF (MAXVAL(div) .GE. 1.0D9) THEN
      !    PRINT *, 'there are nan in div in precal.', MAXVAL(div), MAXVAL(z)
      !    stop
      ! END IF

      this%invRho(:, 1:num_cell) = 1.0D0 / rho(:, 1:num_cell)
      this%HzInvRhoParPsigma(:, 1:num_cell) = Hz(:, 1:num_cell) * this%invRho(:, 1:num_cell) &
                                              * par_P_sigma(:, 1:num_cell)
      this%DM(:, 1:num_cell) = Hz(:, 1:num_cell) * (F_z_chi(:, 1:num_cell) - J_z_psi(:, 1:num_cell) &
                                                    - z(:, 1:num_cell) * div(:, 1:num_cell))

    END ASSOCIATE
    DEALLOCATE (par_P_sigma, &
                J_z_psi, F_z_chi)
  END SUBROUTINE

  FUNCTION constructor(sg) RESULT(this)
    IMPLICIT NONE
    TYPE(PreCal_t) :: this
    TYPE(SingleGrid_t) :: sg

    ASSOCIATE (vLevel => sg%vLevel, &
               num_cell => sg%num_cell, &
               num_icell => sg%num_icell &
               )
      ALLOCATE (this%DM(vLevel, num_cell), this%par_psi_sigma(vLevel, num_cell), &
                this%par_chi_sigma(vLevel, num_cell), this%invRho(vLevel, num_cell), &
                this%HzInvRhoParPsigma(vLevel, num_cell))

      this%DM = 0.0D0
      this%par_psi_sigma = 0.0D0
      this%par_chi_sigma = 0.0D0
      this%invRho = 0.0D0
      this%HzInvRhoParPsigma = 0.0D0

      this%derivatives = CalVerDer_t(sg)
      this%gzm = gzm_t(sg)
    END ASSOCIATE

  END FUNCTION constructor

  SUBROUTINE destroy(this)
    IMPLICIT NONE
    CLASS(PreCal_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%DM)) THEN
      DEALLOCATE (this%DM, this%par_psi_sigma, &
                  this%par_chi_sigma, this%invRho, &
                  this%HzInvRhoParPsigma)
    END IF

  END SUBROUTINE destroy

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(PreCal_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%DM)) THEN
      DEALLOCATE (this%DM, this%par_psi_sigma, &
                  this%par_chi_sigma, this%invRho, &
                  this%HzInvRhoParPsigma)
    END IF

  END SUBROUTINE destructor

END MODULE PreCal_m
