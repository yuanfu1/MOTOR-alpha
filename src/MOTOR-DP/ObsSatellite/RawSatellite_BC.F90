!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC.RawSateObs_BC
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu (wuyali@gbamwf.com), 2022/01/17, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements bias correction (BC) for satellite radiances.

MODULE RawSatellite_BC_m
  USE kinds_m
  USE parameters_m
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE BcPredictors_m
  USE State_m, ONLY: State_t
  USE mpObs_m, ONLY: mpObs_t
  USE YAMLRead_m
  USE FLog_m, ONLY: logger

CONTAINS

  SUBROUTINE Get_pdf_mode(sg, nbins, tb_inv, ObsBias)
    IMPLICIT NONE
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    INTEGER(i_kind), INTENT(in) :: nbins
    REAL(r_kind), INTENT(in) :: tb_inv(:)
    REAL(r_kind), INTENT(INOUT) :: ObsBias
    INTEGER(i_kind) :: nobs, nchans
    INTEGER(i_kind) :: iobs, ichan, ibin, icol, iline
    REAL(r_kind), ALLOCATABLE :: hist(:)
    ! INTEGER(i_kind), PARAMETER    :: nbins = 200                       ! Number of Hist bins.
    REAL(r_kind), PARAMETER       :: maxhist = 10.0                     ! Maximum bin value.
    REAL(r_kind)                  :: zbin                              ! Hist bin width.
    REAL(r_kind) :: mode, modetmp(1), avg
    ! For mpi
    INTEGER(i_kind), ALLOCATABLE :: ncount_group(:), disp_group(:), Cbin(:)
    REAL(r_kind), ALLOCATABLE :: idxArray(:)
    INTEGER(i_kind) :: i, icount

    INCLUDE "mpif.h"

    ALLOCATE (hist(nbins))
    ALLOCATE (Cbin(nbins))
    hist(:) = ZERO
    mode = ZERO
    zbin = 2 * maxhist / REAL(nbins)
    Cbin = 0

    ! Accumulate statistics for histogram
    ! -----------------------------------
    nobs = UBOUND(tb_inv, 1)
    ! nchans = UBOUND(tb_inv, 2)
    ! ALLOCATE(ObsBias(nchans))
    ! PRINT *, 'BC: nobs = ', nobs, MAXVAL(tb_inv), MINVAL(tb_inv)
    !PRINT *, 'BC: nchans = ', nchans

    ALLOCATE (ncount_group(sg%mpddInfo_sg%nproc))
    ALLOCATE (disp_group(sg%mpddInfo_sg%nproc))

    ichan = 1
    ! DO ichan = 5, nchans

    ncount_group = 0
    disp_group = 0

    PRINT *, 'nproc BC: ', sg%mpddInfo_sg%nProc
    ! 将不同进程上的观测个数汇总到主线程 0
    CALL MPI_GATHER(nobs, 1, MPI_INTEGER4, &
                    ncount_group, 1, MPI_INTEGER4, 0, &
                    sg%mpddInfo_sg%comm, sg%mpddInfo_sg%ierr)

    print *, 'isBaseProc BC = ', sg%isBaseProc()
    IF (sg%isBaseProc()) THEN
      PRINT *, "ncount_group BC is ", ncount_group
      ALLOCATE (idxArray(SUM(ncount_group)))
      idxArray = ZERO
      disp_group = 0
      FORALL (i=2:sg%mpddInfo_sg%nproc) disp_group(i) = SUM(ncount_group(1:i - 1)) !Gatherv 时，每个进程所需的位移
      PRINT *, 'check out BC:', sg%mpddInfo_sg%myrank, disp_group, ncount_group
    END IF

    CALL sg%mpddInfo_sg%bcast(ncount_group)
    CALL sg%mpddInfo_sg%bcast(disp_group)

    ! 将所有进程的innovation值汇总，存在idxArray中
    CALL MPI_GATHERV(tb_inv(:), SIZE(tb_inv, 1), MPI_DOUBLE_PRECISION, &
                     idxArray, ncount_group, disp_group, MPI_DOUBLE_PRECISION, 0, &
                     sg%mpddInfo_sg%comm, sg%mpddInfo_sg%ierr)

    ! IF (sg%isBaseProc()) print *, 'check idxArray = ', MAXVAL(idxArray), MINVAL(idxArray)
    IF (sg%isBaseProc()) PRINT *, 'check idxArray = ', SHAPE(idxArray)

    IF (sg%isBaseProc()) THEN
      hist = ZERO
      icount = 0
      avg = ZERO
      DO iobs = 1, SIZE(idxArray, 1)

        ! print *, idxArray(iobs), missing, idxArray(iobs) - missing
        ! IF (ABS(idxArray(iobs) - missing ) < 10.0 ) CYCLE
        ! IF (ABS(idxArray(iobs)) .GT. 50.0 ) CYCLE
        IF (ABS(idxArray(iobs)) .GT. maxhist + 1) CYCLE
        ! PRINT *, 'iobs = ', iobs, ' ', ABS(idxArray(iobs))
        icount = icount + 1
        ibin = NINT((idxArray(iobs) + maxhist) / zbin)
        IF ((ibin > 0) .AND. (ibin <= nbins)) THEN
          hist(ibin) = hist(ibin) + 1
          Cbin(ibin) = Cbin(ibin) + 1
        END IF
        avg = avg + idxArray(iobs)

      END DO
      ! PRINT *, 'Total num of obs used for BC (not counting missing) ', icount
      avg = avg / REAL(icount - 1)

      ! Determine mode of Histogram
      !----------------------------
      ! PRINT *, 'BC: hist ', MAXVAL(hist), MINVAL(hist), SUM(hist), MAXLOC(hist)
      ! PRINT *, 'BC: hist ', Cbin
      IF (SUM(hist(:)) > 0) THEN
        ! PRINT *, 'check BC values: ', MAXLOC(hist(:)), zbin, maxhist
        modetmp(1:1) = MAXLOC(hist(:))*zbin - maxhist
        ObsBias = modetmp(1)
      END IF
      PRINT *, 'mode/ObsBias = ', ObsBias, avg
    END IF
    ! ObsBias(ichan) = avg
    call sg%mpddInfo_sg%bcast(ObsBias)

    IF (ALLOCATED(idxArray)) DEALLOCATE (idxArray)

    ! END DO

    ! tb_inv has a intent(in) attribute
    ! No need to scatterv

    ! BLOCK
    !   REAL(r_kind) :: llsrcEach(size(tb_inv,1))

    !   ! 将经纬度分发给各个线程
    !   CALL MPI_SCATTERV(idxArray, ncount_group, disp_group, MPI_DOUBLE_PRECISION, &
    !                     llsrcEach, size(tb_inv, 1), MPI_DOUBLE_PRECISION, 0, &
    !                     sg%mpddInfo_sg%comm, sg%mpddInfo_sg%ierr)

    ! END BLOCK
    DEALLOCATE (hist, Cbin)
    DEALLOCATE (ncount_group, disp_group)
    ! PRINT *, 'Get_pdf_mode is successfully called'

  END SUBROUTINE Get_pdf_mode

  SUBROUTINE GetBias_ScanAirs(num_rad, PassDomainCheck, npred, zpred, iscanpos, zlat, BCOEF, bias, SCOEF)
    !---------------------------------------         
    !     INPUT.
    !     num_rad  - number of obs
    !     scanpos  - scan postion of obs
    !     zlat     - latitude of obs
    !     zred     - predictor of at each scan position
    !     OUTPUT.
    !     -------
    !     bias    - bias corrections
    !
    !     AUTOHR(S)         : Hua Zhang, 2022/11/18,  Shenzhen
    !---------------------------------------
    IMPLICIT NONE
    
    INTEGER,          INTENT(IN)  :: num_rad, npred
    LOGICAL, INTENT(IN), OPTIONAL       :: PassDomainCheck(:)
    REAL(r_kind),     INTENT(IN)  :: zpred(:,:),zlat(:)
    INTEGER(i_kind),  INTENT(IN)  :: iscanpos(:)
    REAL(r_kind),     INTENT(IN)  :: BCOEF(:)
    REAL(r_kind),     INTENT(IN), OPTIONAL :: SCOEF(:,:)
    REAL(r_kind),     INTENT(INOUT) :: bias(:)
    INTEGER :: io, ipred
    REAL(r_kind), ALLOCATABLE :: pred(:), BIAS8(:)
    REAL(r_kind)    :: SCORR, BCOR = ZERO

    PRINT *, "BC: num_rad = ", num_rad
    PRINT *, "BC coeffs are ", BCOEF(:)
    PRINT *, 'BC predictors are ', MAXVAL(zpred), MINVAL(zpred)
    
   ! ALLOCATE(SCORR(1:nchanl))
    ALLOCATE(pred(1:npred))
    DO io=1,num_rad
  
      IF (PassDomainCheck(io)) THEN
        ! 3.0 Get scan position and scan lat band
        !---------------------------------

          !ISCAN = iscanpos(io)
          !  IF ((ISCAN > maxscan(inst,jsat)) .OR. (ISCAN < 1)) THEN
              ! WRITE(*,*) 'SCAN POSITION: ISCAN ', iscanpos(io),zlat(io)/degree2radian
          !    RETURN
          !  ENDIF

          IF ( PRESENT(SCOEF)) THEN
            CALL GETSCORR(SCORR,zlat(io)/degree2radian,SCOEF(:,:),iscanpos(io))
          ELSE
            SCORR = 0.0D0
          END IF
        !  4.0 Calculate bias correction
        !-------------------------------     
          pred(1:npred)=zpred(io,1:npred)

          IF ( PRESENT(SCOEF)) THEN
            BCOR =  SCORR + SUM(BCOEF(1:NPRED)*PRED(1:NPRED)) + BCOEF(NPRED+1)
          ELSE
            BCOR =  SUM(BCOEF(1:NPRED)*PRED(1:NPRED))
          END IF

          ALLOCATE(BIAS8(NPRED))
          BIAS8(1:NPRED) = BCOEF(1:NPRED) * PRED(1:NPRED)
          ! PRINT *, 'why nan? bcoef: ', bcoef
          ! PRINT *, 'why nan? pred: ', pred

          bias(io)=BCOR
          ! PRINT *, 'Check output BIAS values of each predictor: ', BIAS8
          DEALLOCATE(BIAS8)
        END IF
        
      ENDDO
      ! PRINT *, 'check bias inside: ', maxval(bias), minval(bias)

    DEALLOCATE (pred)
  END SUBROUTINE GetBias_ScanAirs

  SUBROUTINE Static_BC(sg, pred, beta, ObsBias)
    IMPLICIT NONE
    TYPE(SingleGrid_t), INTENT(IN) :: sg
    REAL(r_kind), INTENT(IN) :: pred(:,:), beta(:,:)
    REAL(r_kind), INTENT(INOUT) :: ObsBias
    INTEGER(i_kind) :: nobs, nchans
    INTEGER(i_kind) :: iobs, ichan
    INTEGER(i_kind) :: i, icount

    ObsBias = 0.0D0
    IF (sg%isBaseProc()) THEN

      PRINT *, 'Static_BC for each channel'
      ! DO ichan = 1, nchans
      !   DO ipred = 1, npred
      !     ObsBias(iobs, ichan) = ObsBias(iobs, ichan) + pred(iobs, ipred)*beta(iobs, ipred)
      !   END DO
      ! END DO

    END IF
    ! call sg%mpddInfo_sg%bcast(ObsBias(ichan))

    PRINT *, 'Static_BC is successfully called'

  END SUBROUTINE Static_BC

  ! SUBROUTINE Var_BC(this)
  !   IMPLICIT NONE
  !   CLASS(RawSatellite_BC_t) :: this
  !   TYPE(State_t), INTENT(IN) :: X
  !   TYPE(ObsSet_t), INTENT(IN) :: Y

  !   PRINT *, 'Nothing yet'
  !   CALL BC_predictors(this, X, Y)

  !   DO ipred = 1, npred
  !     this%ObsBias = predictors(ipred)*beta(ipred)
  !   END DO

  !   ! Re-calculate Jo with bias considered
  !   Jo_w_bias = y - HX - this%ObsBias

  !   ! Calculate Jb with coefs as extended control variables
  !   Jb_w_bias = this%B_beta.SQRTINVMULT. (beta - beta_b)

  !   ! Re-calculate J
  !   Jfunc = Jfunc + Jb_w_bias.DOT.Jb_w_bias + Jo_w_bias.DOT.Jo_w_bias

  ! END SUBROUTINE Var_BC

  ! SUBROUTINE var_BC_grad(this)
  !   IMPLICIT NONE
  !   CLASS(RawSateObs_BC_t) :: this

  !   ! Re-calculate grad_Jo with bias considered
  !   grad_Jo_w_bias = y - HX - beta

  !   ! Calculate grad_Jb with coefs as extended control variables
  !   grad_Jb_w_bias = this%B_beta.SQRTINVMULT.beta

  !   ! Re-calculate grad_J
  !   grad_J = grad_J + grad_Jb_w_bias*2.0D0 + grad_Jo_w_bias*2.0D0

  ! END SUBROUTINE var_BC_grad

  ! SUBROUTINE var_BC_tl()
  !   IMPLICIT NONE
  !   CLASS(RawSateObs_BC_t) :: this

  !   PRINT *, 'Nothing yet'

  ! END SUBROUTINE var_BC_tl

  ! SUBROUTINE var_BC_ad()
  !   IMPLICIT NONE
  !   CLASS(RawSateObs_BC_t) :: this

  !   PRINT *, 'Nothing yet'

  ! END SUBROUTINE var_BC_ad

END MODULE RawSatellite_BC_m
