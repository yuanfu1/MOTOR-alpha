PROGRAM Test_CalVerDer_tlm
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE CalVerDer_tlm_m, ONLY: CalVerDer_tlm_t
  IMPLICIT NONE

  TYPE(SingleGrid_t), TARGET :: sg
  TYPE(CalVerDer_tlm_t) :: calverder_tlm

  REAL(r_kind), ALLOCATABLE :: A(:, :), A_tlm(:, :), A_pert(:, :), parA_sigma(:, :), parA_sigma_tlm(:, :), parA_sigma_pert(:, :)
  REAL(r_kind) :: epsilon

  CHARACTER(LEN=1024) :: configFile

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testCalVerDer_TLAD.yaml"

  ! Initialize SingleGrid
  CALL sg%initialize(configFile, mpddGlob)

  ! Initialize CalVerDer_tlm
  CALL calverder_tlm%initialize(sg)

  ! Allocate and initialize arrays
  ALLOCATE (A(sg%vLevel, sg%num_cell), A_tlm(sg%vLevel, sg%num_cell), A_pert(sg%vLevel, sg%num_cell))
  ALLOCATE (parA_sigma(sg%vLevel, sg%num_cell), parA_sigma_tlm(sg%vLevel, sg%num_cell), parA_sigma_pert(sg%vLevel, sg%num_cell))
  A = 0.05_R_KIND
  A_tlm = 0.01_R_KIND
  epsilon = 1.0E-6_R_KIND

  ! Test FirstOrder_tlm
  CALL Test_Subroutine("FirstOrder_tlm", calverder_tlm%FirstOrder_tlm, A, parA_sigma, A_tlm, parA_sigma_tlm, A_pert, parA_sigma_pert, epsilon)

  ! Test SecondOrder_tlm
  CALL Test_Subroutine("SecondOrder_tlm", calverder_tlm%SecondOrder_tlm, A, parA_sigma, A_tlm, parA_sigma_tlm, A_pert, parA_sigma_pert, epsilon)

  ! Test FirstOrder_half_tlm
  CALL Test_Subroutine("FirstOrder_half_tlm", calverder_tlm%FirstOrder_half_tlm, A, parA_sigma, A_tlm, parA_sigma_tlm, A_pert, parA_sigma_pert, epsilon)

CONTAINS

  SUBROUTINE Test_Subroutine(name, SUBROUTINE, A, parA_sigma, A_tlm, parA_sigma_tlm, A_pert, parA_sigma_pert, epsilon)
    CHARACTER(LEN=*), INTENT(IN) :: name
    PROCEDURE(subroutine_interface), POINTER :: SUBROUTINE
    REAL(r_kind), INTENT(IN) :: A(:, :), A_tlm(:, :), epsilon
    REAL(r_kind), INTENT(OUT) :: parA_sigma(:, :), parA_sigma_tlm(:, :), parA_sigma_pert(:, :)

    ! Call the subroutine with the original input
    CALL SUBROUTINE(A, parA_sigma)

    ! Call the subroutine with the tangent linear model input
    A_pert = A + epsilon * A_tlm
    CALL SUBROUTINE(A_pert, parA_sigma_pert)

    ! Calculate the difference
    parA_sigma_tlm = (parA_sigma_pert - parA_sigma) / epsilon

    PRINT *, name, " - Relative Difference (should be close to 1): ", SUM(ABS(parA_sigma_tlm - A_tlm)) / SUM(ABS(A_tlm))
  END SUBROUTINE Test_Subroutine

  ABSTRACT INTERFACE
    SUBROUTINE subroutine_interface(A, parA_sigma)
      USE kinds_m, ONLY: r_kind
      REAL(r_kind), INTENT(IN) :: A(:, :)
      REAL(r_kind), INTENT(OUT) :: parA_sigma(:, :)
    END SUBROUTINE subroutine_interface
  END INTERFACE

END PROGRAM Test_CalVerDer_tlm
