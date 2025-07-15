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
SUBROUTINE RK4_AD(this, X)

  USE RK4B_m

  IMPLICIT NONE

  CLASS(PhyCon_t), INTENT(INOUT) :: this
  TYPE(State_t), INTENT(IN) :: X

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
             tt_in => X%sg%tt, &
             uwnd_out => this%uwnd_out, &
             vwnd_out => this%vwnd_out, &
             wwnd_out => this%wwnd_out, &
             theta_out => this%theta_out, &
             qvapor_out => this%qvapor_out, &
             pres_out => this%pres_out, &
             uwnd_outb => this%uwnd_outb, &
             vwnd_outb => this%vwnd_outb, &
             wwnd_outb => this%wwnd_outb, &
             theta_outb => this%theta_outb, &
             qvapor_outb => this%qvapor_outb, &
             pres_outb => this%pres_outb, &
             !  precip => this%precip, &
             precipb => this%precipb, &
             dthdt => this%dthdt, &
             dqdt => this%dqdt, &
             ztop => this%ztop, &
             cu_dt => this%cu_dt, &
             tt_out => this%tt_out, &
             kt => this%kt, &
             size_3m => this%size_3m)

    size_q_dim1 = SIZE(qvapor_out, dim=1)
    size_q_dim2 = SIZE(qvapor_out, dim=2)
    size_q_dim3 = SIZE(qvapor_out, dim=3)

    CALL RK4_B(uwnd_out, uwnd_outb, vwnd_out, vwnd_outb, wwnd_out, &
    & wwnd_outb, theta_out, theta_outb, qvapor_out, qvapor_outb, pres_out, &
    & pres_outb, precipb, dthdt, dqdt, size_q_dim1, size_q_dim2, size_q_dim3, &
    & cu_dt, ztop, cell_dist, cell_stcl, gph, cell_type, sigma, topo, &
    & num_cell, num_icell, size_3m)

  END ASSOCIATE

END SUBROUTINE RK4_AD
