!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.Constrants.Cumulus.Prepare data for Adjoint model.
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research
! Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong Chen
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong Chen (jchen@link.cuhk.edu.hk), 2022/2/18, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
  SUBROUTINE AdjDataPre(this, X, dX)

    CLASS(PhyCon_t), INTENT(INOUT) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t), INTENT(IN) :: dX
    INTEGER :: it_in, size_tt_out, it_out

    ASSOCIATE (sg => X%sg, &
               duwnd => dX%Fields(dX%getVarIdx('uwnd'))%DATA, &
               dvwnd => dX%Fields(dX%getVarIdx('vwnd'))%DATA, &
               dwwnd => dX%Fields(dX%getVarIdx('wwnd'))%DATA, &
               dtheta => dX%Fields(dX%getVarIdx('theta'))%DATA, &
               dqvapor => dX%Fields(dX%getVarIdx('qvapor'))%DATA, &
               dpres => dX%Fields(dX%getVarIdx('pres'))%DATA, &
               tt_in => X%sg%tt)

      size_tt_out = (tt_in(SIZE(tt_in)) - tt_in(1)) / this%cu_dt + 1.0D0

      ALLOCATE (this%uwnd_outb(sg%vLevel, sg%num_cell, size_tt_out), &
                this%vwnd_outb(sg%vLevel, sg%num_cell, size_tt_out), &
                this%wwnd_outb(sg%vLevel, sg%num_cell, size_tt_out), &
                this%theta_outb(sg%vLevel, sg%num_cell, size_tt_out), &
                this%qvapor_outb(sg%vLevel, sg%num_cell, size_tt_out), &
                this%pres_outb(sg%vLevel, sg%num_cell, size_tt_out), &
                this%precipb(1, sg%num_cell, size_tt_out) &
                )

      this%uwnd_outb = 0.0D0
      this%vwnd_outb = 0.0D0
      this%wwnd_outb = 0.0D0
      this%theta_outb = 0.0D0
      this%qvapor_outb = 0.0D0
      this%pres_outb = 0.0D0
      this%precipb = 0.0D0

      it_in = 1
      DO it_out = 1, size_tt_out
        IF (MOD(it_out, this%kt) .EQ. 1) THEN
          this%uwnd_outb(:, :, it_out) = duwnd(:, :, it_in)
          this%vwnd_outb(:, :, it_out) = dvwnd(:, :, it_in)
          this%wwnd_outb(:, :, it_out) = dwwnd(:, :, it_in)
          this%theta_outb(:, :, it_out) = dtheta(:, :, it_in)
          this%qvapor_outb(:, :, it_out) = dqvapor(:, :, it_in)
          this%pres_outb(:, :, it_out) = dpres(:, :, it_in)
          it_in = it_in + 1
        END IF
      END DO

    END ASSOCIATE
  END SUBROUTINE
