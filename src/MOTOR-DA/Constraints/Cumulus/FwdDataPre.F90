!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.Constrants.Cumulus.Prepare data for foreward model.
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research
! Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Jilong Chen
! VERSION           : V 0.0
! HISTORY           :
!   Created by Jilong Chen (jchen@link.cuhk.edu.hk), 2022/2/18, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

  SUBROUTINE FwdDataPre(this, X)

    CLASS(PhyCon_t), INTENT(INOUT) :: this
    TYPE(State_t), INTENT(IN) :: X
    INTEGER :: i, size_tt_in, it_in, size_tt_out, it_out, vLevel
    INTEGER(i_kind), ALLOCATABLE :: cell_eage(:)

    ASSOCIATE (sg => X%sg, &
               uwnd => X%Fields(X%getVarIdx('uwnd'))%DATA, &
               vwnd => X%Fields(X%getVarIdx('vwnd'))%DATA, &
               wwnd => X%Fields(X%getVarIdx('wwnd'))%DATA, &
               theta => X%Fields(X%getVarIdx('theta'))%DATA, &
               qvapor => X%Fields(X%getVarIdx('qvapor'))%DATA, &
               pres => X%Fields(X%getVarIdx('pres'))%DATA, &
               tt_in => X%sg%tt)

      size_tt_in = SIZE(tt_in)
      size_tt_out = (tt_in(SIZE(tt_in)) - tt_in(1)) / this%cu_dt + 1.0D0
      this%kt = (tt_in(2) - tt_in(1)) / this%cu_dt
      vLevel = sg%vLevel
      IF (.NOT. ALLOCATED(this%uwnd_out)) THEN
        ALLOCATE (this%tt_out(size_tt_out), &
                  this%uwnd_out(sg%vLevel, sg%num_cell, size_tt_out), &
                  this%vwnd_out(sg%vLevel, sg%num_cell, size_tt_out), &
                  this%wwnd_out(sg%vLevel, sg%num_cell, size_tt_out), &
                  this%theta_out(sg%vLevel, sg%num_cell, size_tt_out), &
                  this%qvapor_out(sg%vLevel, sg%num_cell, size_tt_out), &
                  this%pres_out(sg%vLevel, sg%num_cell, size_tt_out), &
                  this%dthdt(sg%vLevel, sg%num_cell, 4, size_tt_out), &
                  this%dqdt(sg%vLevel, sg%num_cell, 4, size_tt_out), &
                  this%precip(1, sg%num_cell, size_tt_out) &
                  )
      END IF

      this%uwnd_out = 0.0D0
      this%vwnd_out = 0.0D0
      this%wwnd_out = 0.0D0
      this%theta_out = 0.0D0
      this%qvapor_out = 0.0D0
      this%pres_out = 0.0D0
      this%precip = 0.0D0
      this%dthdt = 0.0D0
      this%dqdt = 0.0D0

      DO it_out = 1, size_tt_out
        this%tt_out(it_out) = tt_in(1) + this%cu_dt * (REAL(it_out, 8) - 1.0D0)
      END DO

      cell_eage = PACK([(i, i=1, sg%num_cell)], sg%cell_type == 1 .OR. sg%cell_type == 2)

      this%uwnd_out(:, :, 1) = uwnd(:, :, 1)
      this%vwnd_out(:, :, 1) = vwnd(:, :, 1)
      this%wwnd_out(:, :, 1) = wwnd(:, :, 1)
      this%theta_out(:, :, 1) = theta(:, :, 1)
      this%pres_out(:, :, 1) = pres(:, :, 1)
      this%qvapor_out(:, :, 1) = qvapor(:, :, 1)

      DO it_out = 2, size_tt_out
        DO it_in = 1, size_tt_in
          IF (dabs(tt_in(it_in) - this%tt_out(it_out)) < 1.0D-14) THEN
            this%uwnd_out(:, :, it_out) = uwnd(:, :, it_in)
            this%vwnd_out(:, :, it_out) = vwnd(:, :, it_in)
            this%wwnd_out(:, :, it_out) = wwnd(:, :, it_in)
            this%pres_out(:, :, it_out) = pres(:, :, it_in)
            this%theta_out(2:vLevel, cell_eage, it_out) = theta(2:vLevel, cell_eage, it_in)
            this%qvapor_out(2:vLevel, cell_eage, it_out) = qvapor(2:vLevel, cell_eage, it_in)
            this%theta_out(1, :, it_out) = theta(1, :, it_in)
            this%qvapor_out(1, :, it_out) = qvapor(1, :, it_in)
            EXIT
          ELSE IF (tt_in(it_in) > this%tt_out(it_out) .OR. it_in == size_tt_in) THEN
            this%uwnd_out(:, :, it_out) = uwnd(:, :, it_in - 1) + &
                                          (uwnd(:, :, it_in) - uwnd(:, :, it_in - 1)) / &
                                          (tt_in(it_in) - tt_in(it_in - 1)) * &
                                          (this%tt_out(it_out) - tt_in(it_in - 1))
            this%vwnd_out(:, :, it_out) = vwnd(:, :, it_in - 1) + &
                                          (vwnd(:, :, it_in) - vwnd(:, :, it_in - 1)) / &
                                          (tt_in(it_in) - tt_in(it_in - 1)) * &
                                          (this%tt_out(it_out) - tt_in(it_in - 1))
            this%wwnd_out(:, :, it_out) = wwnd(:, :, it_in - 1) + &
                                          (wwnd(:, :, it_in) - wwnd(:, :, it_in - 1)) / &
                                          (tt_in(it_in) - tt_in(it_in - 1)) * &
                                          (this%tt_out(it_out) - tt_in(it_in - 1))
            this%pres_out(:, :, it_out) = pres(:, :, it_in - 1) + &
                                          (pres(:, :, it_in) - pres(:, :, it_in - 1)) / &
                                          (tt_in(it_in) - tt_in(it_in - 1)) * &
                                          (this%tt_out(it_out) - tt_in(it_in - 1))
            this%theta_out(2:vLevel, cell_eage, it_out) = theta(2:vLevel, cell_eage, it_in - 1) + &
                                                          (theta(2:vLevel, cell_eage, it_in) - theta(2:vLevel, cell_eage, it_in - 1)) / &
                                                          (tt_in(it_in) - tt_in(it_in - 1)) * &
                                                          (this%tt_out(it_out) - tt_in(it_in - 1))
            this%qvapor_out(2:vLevel, cell_eage, it_out) = qvapor(2:vLevel, cell_eage, it_in - 1) + &
                                                           (qvapor(2:vLevel, cell_eage, it_in) - qvapor(2:vLevel, cell_eage, it_in - 1)) / &
                                                           (tt_in(it_in) - tt_in(it_in - 1)) * &
                                                           (this%tt_out(it_out) - tt_in(it_in - 1))
            this%theta_out(1, :, it_out) = theta(1, :, it_in - 1) + &
                                           (theta(1, :, it_in) - theta(1, :, it_in - 1)) / &
                                           (tt_in(it_in) - tt_in(it_in - 1)) * &
                                           (this%tt_out(it_out) - tt_in(it_in - 1))
            this%qvapor_out(1, :, it_out) = qvapor(1, :, it_in - 1) + &
                                            (qvapor(1, :, it_in) - qvapor(1, :, it_in - 1)) / &
                                            (tt_in(it_in) - tt_in(it_in - 1)) * &
                                            (this%tt_out(it_out) - tt_in(it_in - 1))
            EXIT
          END IF
        END DO
      END DO

    END ASSOCIATE
  END SUBROUTINE
