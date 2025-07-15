!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather
! Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong CHEN
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong CHEN (jchen@link.cuhk.edu.hk), 2021/1/26, @GBA-MWF,
!   Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!

MODULE RK4TL_m
  USE State_m, ONLY: State_t
  USE parameters_m, ONLY: g, cp, Lv, k_d, r_k_d, epsw_r, gamma_d, r_d
  USE kinds_m, ONLY: i_kind, r_kind
  USE PhyCon_m, ONLY: PhyCon_t
  USE RK4D_m, ONLY: RK4_D, CUMFWD_D, DIVERGENCE_D, TENDENCY_D, INTERP1D_D
  ! USE CumFwd_m
  ! USE RK4B_m
  ! USE DataPre_m, ONLY: DataPre

  TYPE CumTL_t
    CHARACTER(len=1024) :: configFile
    REAL(r_kind), ALLOCATABLE :: uwnd_d(:, :, :), vwnd_d(:, :, :), wwnd_d(:, :, :), &
                                 pres_d(:, :, :), theta_d(:, :, :), qvapor_d(:, :, :), &
                                 precipd(:, :, :)

  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: TLDataPre
    PROCEDURE, PUBLIC, PASS(this) :: RK4_TL

    FINAL :: destroy
  END TYPE CumTL_t
  INTERFACE CumTL_t
    PROCEDURE :: constructor
  END INTERFACE CumTL_t

CONTAINS
  SUBROUTINE TLDataPre(this, PhyCon)

    CLASS(CumTL_t), INTENT(INOUT) :: this
    TYPE(PhyCon_t), INTENT(IN) :: PhyCon

    INTEGER(i_kind) :: size_dim1, size_dim2, size_dim3

    ASSOCIATE (uwnd_out => PhyCon%uwnd_out)

      size_dim1 = SIZE(uwnd_out, dim=1)
      size_dim2 = SIZE(uwnd_out, dim=2)
      size_dim3 = SIZE(uwnd_out, dim=3)

      ALLOCATE (this%uwnd_d(size_dim1, size_dim2, size_dim3), &
                this%vwnd_d(size_dim1, size_dim2, size_dim3), &
                this%wwnd_d(size_dim1, size_dim2, size_dim3), &
                this%theta_d(size_dim1, size_dim2, size_dim3), &
                this%qvapor_d(size_dim1, size_dim2, size_dim3), &
                this%pres_d(size_dim1, size_dim2, size_dim3), &
                this%precipd(1, size_dim2, size_dim3))
      this%uwnd_d = 0.0D0
      this%vwnd_d = 0.0D0
      this%wwnd_d = 0.0D0
      this%theta_d = 0.0D0
      this%qvapor_d = 0.0D0
      this%pres_d = 0.0D0
      this%precipd = 0.0D0

    END ASSOCIATE

  END SUBROUTINE

  SUBROUTINE RK4_TL(this, PhyCon, X)

    CLASS(CumTL_t), INTENT(INOUT) :: this
    TYPE(PhyCon_t), INTENT(INOUT) :: PhyCon
    TYPE(state_t), INTENT(IN) :: X

    INTEGER(i_kind) :: size_q_dim1, size_q_dim2, size_q_dim3

    ASSOCIATE (cell_type => X%sg%cell_type, &
               cell_stcl => X%sg%cell_stcl, &
               cell_dist => X%sg%cell_dist, &
               vLevel => X%sg%vLevel, &
               num_cell => X%sg%num_cell, &
               num_icell => X%sg%num_icell, &
               gph => X%sg%zHght, &
               sigma => X%sg%sigma, &
               topo => X%sg%topo, &
               uwnd_out => PhyCon%uwnd_out, &
               vwnd_out => PhyCon%vwnd_out, &
               wwnd_out => PhyCon%wwnd_out, &
               theta_out => PhyCon%theta_out, &
               qvapor_out => PhyCon%qvapor_out, &
               pres_out => PhyCon%pres_out, &
               uwnd_outd => this%uwnd_d, &
               vwnd_outd => this%vwnd_d, &
               wwnd_outd => this%wwnd_d, &
               theta_outd => this%theta_d, &
               qvapor_outd => this%qvapor_d, &
               pres_outd => this%pres_d, &
               precip => PhyCon%precip, &
               precipd => this%precipd, &
               ztop => PhyCon%ztop, &
               cu_dt => PhyCon%cu_dt, &
               size_3m => PhyCon%size_3m)

      size_q_dim1 = SIZE(qvapor_out, dim=1)
      size_q_dim2 = SIZE(qvapor_out, dim=2)
      size_q_dim3 = SIZE(qvapor_out, dim=3)

      CALL RK4_D(uwnd_out, uwnd_outd, vwnd_out, vwnd_outd, wwnd_out, &
      & wwnd_outd, theta_out, theta_outd, qvapor_out, qvapor_outd, pres_out, &
      & pres_outd, precip, precipd, size_q_dim1, size_q_dim2, size_q_dim3, &
      & cu_dt, ztop, cell_dist, cell_stcl, gph, cell_type, sigma, topo, &
      & num_cell, num_icell, size_3m)

    END ASSOCIATE

  END SUBROUTINE

  IMPURE ELEMENTAL SUBROUTINE destroy(this)
    IMPLICIT NONE
    TYPE(CumTL_t), INTENT(INOUT) :: this

    IF (ALLOCATED(this%uwnd_d)) DEALLOCATE (this%uwnd_d, &
                                            this%vwnd_d, this%wwnd_d, this%theta_d, &
                                            this%qvapor_d, this%pres_d, this%precipd)
  END SUBROUTINE

  FUNCTION constructor(configFile) RESULT(this)
    IMPLICIT NONE
    TYPE(CumTL_t) :: this
    CHARACTER(len=1024) :: configFile

    this%configFile = configFile

  END FUNCTION constructor

END MODULE RK4TL_m
