MODULE GenSgDiffCoef_TL_AD_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE gzm_tlm_m, ONLY: gzm_tlm_t, Divergen_TLM, Jacobian_TLM, Laplacia_TLM
  USE gzm_adj_m, ONLY: gzm_adj_t, Divergen_AD, Jacobian_AD, Laplacia_AD
  IMPLICIT NONE

  TYPE :: GenSgDiffCoef_TL_AD_t
    TYPE(gzm_tlm_t), POINTER :: gzm_tlm
    TYPE(gzm_adj_t), POINTER :: gzm_adj
  CONTAINS
    PROCEDURE :: GenLapaceParams_TL_AD
  END TYPE GenSgDiffCoef_TL_AD_t

CONTAINS

  SUBROUTINE GenLapaceParams_TL_AD(this, sg)
    IMPLICIT NONE
    CLASS(GenSgDiffCoef_TL_AD_t), INTENT(INOUT) :: this
    TYPE(SingleGrid_t), INTENT(INOUT) :: sg

    ALLOCATE (sg%F_z_z(sg%vLevel, sg%num_cell))
    ALLOCATE (sg%Lz(sg%vLevel, sg%num_cell))
    ALLOCATE (sg%F_Hz_z(sg%vLevel, sg%num_cell))
    ALLOCATE (sg%F_invHz_z(sg%vLevel, sg%num_cell))

    sg%F_Hz_z = 0.0D0
    sg%F_z_z = 0.0D0
    sg%Lz = 0.0D0
    sg%F_invHz_z = 0.0D0

    ! Use the TLM and AD versions of the gzm subroutines
    CALL Divergen_TLM(this%gzm_tlm, sg%Hz, sg%zHght, sg%F_Hz_z)
    CALL Divergen_TLM(this%gzm_tlm, sg%zHght, sg%zHght, sg%F_z_z)
    CALL Divergen_TLM(this%gzm_tlm, sg%parz_parsigma, sg%zHght, sg%F_invHz_z)
    CALL Laplacia_TLM(this%gzm_tlm, sg%zHght, sg%Lz)

    CALL sg%ExchangeMatOnHalo2D(sg%vLevel, sg%F_Hz_z)
    CALL sg%ExchangeMatOnHalo2D(sg%vLevel, sg%F_z_z)
    CALL sg%ExchangeMatOnHalo2D(sg%vLevel, sg%F_invHz_z)
    CALL sg%ExchangeMatOnHalo2D(sg%vLevel, sg%Lz)
  END SUBROUTINE GenLapaceParams_TL_AD

END MODULE GenSgDiffCoef_TL_AD_m
