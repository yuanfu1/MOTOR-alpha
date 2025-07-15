!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-QC. Conventional Obs QC
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for
!                     Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Hua Zhang
! VERSION           : V 0.0
! HISTORY           : Origionally from RCNMP
!
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This module implements a conventional Obs QC.
MODULE GtsQC_m
  USE Meteoro_Constants, ONLY: threshold
  USE Meteoro_type_define
  USE extreme_value_m
  USE interp_m
  USE qc_flag_set_m

  IMPLICIT NONE

!TYPE GtsQC_t

!CONTAINS
  !PROCEDURE, PUBLIC  :: p_extreme_qc
  !PROCEDURE, PUBLIC  :: h_extreme_qc
  !PROCEDURE, PUBLIC  :: t_extreme_qc
  !PROCEDURE, PUBLIC  :: w_extreme_qc
  !PROCEDURE, PUBLIC  :: w_extreme_qc_with_h
  !PROCEDURE, PUBLIC  :: ps_extreme_qc
  !PROCEDURE, PUBLIC  :: logp_interp
  !PROCEDURE, PUBLIC  :: logp_interh
  !PROCEDURE, PUBLIC  :: comp_extr
  !PROCEDURE, PUBLIC  :: comp_extrw
  !PROCEDURE, PUBLIC  :: t_td_consist_check
  !PROCEDURE, PUBLIC  :: wind_consist_check
  !PROCEDURE, PUBLIC  :: ddff_uv
  !PROCEDURE, PUBLIC  :: FSIGTOBT
  !PROCEDURE, PUBLIC  :: FSIGTOBH
  !PROCEDURE, PUBLIC  :: PTOSIG
  !PROCEDURE, PUBLIC  :: complex_dma

!END TYPE GtsQC_t

CONTAINS

  SUBROUTINE p_extreme_qc(p, flag)
!  purpose:
!  this is the  programm of extreme extreme (or gross) check for
!  for pressure data of radiosonding  observasion
!
!-------------------------------------------------------------------
! the input variables
    IMPLICIT NONE

! Subroutine arguments:
    REAL                                 :: P, robsmissing
    INTEGER :: flag
!----------------------------------
    IF (p .LT. 0 .OR. p .GT. 1100) THEN
      CALL flag_replace(flag, extreme_flag, 2)
      CALL flag_replace(flag, DMA_flag, 2)

!      if(logic_temp_extreme_print) print '(a,10f10.2)','p error ',temp(num_stn)%header_info%stn_id,p

    END IF
    RETURN
  END SUBROUTINE p_extreme_qc

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE h_extreme_qc(p, lat, h, flag)
!_____________________________________________________________________
! purpose:
!     this is the  programm of extreme extreme (or gross) check
!     for height data of radiosonding  observasion
!_______________________________________________________________________
! method:
!   1. compare the observation data with it's extreme extreme
!      if exceed the extreme 2  ,it is considered wrong
!      if exceed the extreme 1 ,it is considered suspect
!______________________________________________________________________
!______________________________________________________________________
! the input variables
    IMPLICIT NONE

! Subroutine arguments:
! Locale arguments
    INTEGER                         :: iol, lll, flag, iflag
    INTEGER, PARAMETER              :: check_id = 10**2
    REAL                            :: hmax1, hmax2, hmin1, hmin2, p, h, lat

    IF (ABS(lat) .LT. 45) THEN
      CALL logp_interp(h_min_mid_sus, p, hmin1)
      CALL logp_interp(h_min_mid_err, p, hmin2)
      CALL logp_interp(h_max_mid_sus, p, hmax1)
      CALL logp_interp(h_max_mid_err, p, hmax2)
    ELSE
      CALL logp_interp(h_min_hig_sus, p, hmin1)
      CALL logp_interp(h_min_hig_err, p, hmin2)
      CALL logp_interp(h_max_hig_sus, p, hmax1)
      CALL logp_interp(h_max_hig_err, p, hmax2)
    END IF

    CALL comp_extr(hmin1, hmax1, hmin2, hmax2, h, iflag)
    IF (flag == 2) THEN
      CALL flag_replace(flag, DMA_flag, iflag)
      CALL flag_replace(flag, extreme_flag, iflag)
!        if(logic_temp_extreme_print)&
!        print'(a,10f10.2)','hpi error  ',temp(num_stn)%header_info%stn_id,temp(num_stn)%level_data(iol)%p%value,h

    ELSE
      CALL flag_replace(flag, extreme_flag, iflag)
    END IF
    RETURN
  END SUBROUTINE h_extreme_qc

!____________________________________________________________
! Temperature extreme extreme check
!_____________________________________________________________
  SUBROUTINE t_extreme_qc(p, lat, t, flag)
!_____________________________________________________________________
! purpose:
!     this is the  programm of extreme extreme (or gross) check
!     for height data of radiosonding  observasion
!_______________________________________________________________________
! method:
!   1. compare the observation data with it's extreme extreme
!      if exceed the extreme 2  ,it is considered wrong
!      if exceed the extreme 1 ,it is considered suspect
!______________________________________________________________________
!______________________________________________________________________
! the input variables
    IMPLICIT NONE

! Subroutine arguments:
! Locale arguments
    INTEGER                         :: iol, lll, iflag, flag
    INTEGER, PARAMETER              :: check_id = 10**2
    REAL                            :: tmax1, tmax2, tmin1, tmin2, p, t, lat

    IF (ABS(lat) .LT. 45) THEN
      CALL logp_interp(t_min_mid_sus, p, tmin1)
      CALL logp_interp(t_min_mid_err, p, tmin2)
      CALL logp_interp(t_max_mid_sus, p, tmax1)
      CALL logp_interp(t_max_mid_err, p, tmax2)
    ELSE
      CALL logp_interp(t_min_hig_sus, p, tmin1)
      CALL logp_interp(t_min_hig_err, p, tmin2)
      CALL logp_interp(t_max_hig_sus, p, tmax1)
      CALL logp_interp(t_max_hig_err, p, tmax2)
    END IF

    CALL comp_extr(tmin1, tmax1, tmin2, tmax2, t, iflag)
    IF (iflag == 2) THEN
      CALL flag_replace(flag, DMA_flag, iflag)
      CALL flag_replace(flag, extreme_flag, iflag)
!        if(logic_temp_extreme_print)&
!        print'(a,10f10.2)','hpi error  ',p,t

    ELSE
      CALL flag_replace(flag, extreme_flag, iflag)
    END IF
    RETURN
  END SUBROUTINE t_extreme_qc

!____________________________________________________________________
! Wind extreme extreme check
!____________________________________________________________________

  SUBROUTINE w_extreme_qc(p, lat, w, flag)

    IMPLICIT NONE

! Locale arguments
    INTEGER                         :: iol, lll, flag, iflag
    INTEGER, PARAMETER              :: check_id = 10**2
    REAL                            :: wmax, p, w, lat
!-----------------------------------------------------
    IF (w .LT. 0) THEN
      iflag = 2
      CALL flag_replace(flag, DMA_flag, iflag)
      CALL flag_replace(flag, extreme_flag, iflag)
    ELSE
      CALL logp_interp(v_max_sus, p, wmax)
      CALL comp_extrw(wmax, w, iflag)
      IF (iflag == 2) THEN
        CALL flag_replace(flag, DMA_flag, iflag)
        CALL flag_replace(flag, extreme_flag, iflag)
      ELSE
        CALL flag_replace(flag, extreme_flag, iflag)
      END IF
    END IF

    RETURN
!----------------------------------------------------------------------
  END SUBROUTINE w_extreme_qc

!____________________________________________________________________

  SUBROUTINE w_extreme_qc_with_h(h, lat, w, flag)

    IMPLICIT NONE

! Locale arguments
    INTEGER                         :: iol, lll, flag, iflag
    INTEGER, PARAMETER              :: check_id = 10**2
    REAL                            :: wmax, p, w, lat, h
!-----------------------------------------------------
    IF (w .LT. 0) THEN
      iflag = 2
      CALL flag_replace(flag, DMA_flag, iflag)
      CALL flag_replace(flag, extreme_flag, iflag)
    ELSE
      CALL logp_interh(v_max_sus, h, wmax)
      CALL comp_extrw(wmax, w, iflag)
      CALL flag_replace(flag, extreme_flag, iflag)
      IF (iflag == 2) THEN
        CALL flag_replace(flag, DMA_flag, iflag)
      END IF
    END IF
    RETURN
!----------------------------------------------------------------------
  END SUBROUTINE w_extreme_qc_with_h

  SUBROUTINE ps_extreme_qc(ps, flag)

    IMPLICIT NONE
    ! CLASS(GtsQC_t) :: this
    ! TYPE(qc_flag_set_t) :: qc_flag
! Subroutine arguments:
    REAL                            :: ps
    INTEGER :: flag

    IF (ps > 1100 .OR. ps < 850) THEN
      CALL flag_replace(flag, DMA_flag, 2)
      CALL flag_replace(flag, extreme_flag, 2)
    END IF
    RETURN
  END SUBROUTINE ps_extreme_qc

!-----------------------------------------------------------------------------------------------

  SUBROUTINE logp_interp(ax1, p, amax)
!_______________________________________________________________________
! purpose:
!      calculate the p level extreme extreme from ax1
!_______________________________________________________________________
! the input variables
!    ax1 the array of the maximum extreme for standard levels
!    p: observation pressure
!______________________________________________________________________
! the output variables are the same as the input variables
!    amax: the maximum extreme for p level
!______________________________________________________________________
    !  USE standard_pressure

    IMPLICIT NONE
    ! CLASS(GtsQC_t) :: this
! Local arguments:
    REAL            ::  ax1(maxstl)
    REAL            :: plin(maxstl), p, alogp, amax
    INTEGER         :: i, istl, l1, l2

    DO i = 1, maxstl
      plin(i) = alog(stdp(i))
    END DO

    IF (p .GT. stdp(1)) THEN
      l1 = 1
      l2 = 2
      alogp = alog(p)
      amax = ax1(l1) + (ax1(l1) - ax1(l2)) * (plin(l1) - alogp) / (plin(l2) - plin(l1))
      RETURN
    END IF

    DO istl = 1, maxstl
      IF (ABS(p - stdp(istl)) .LT. 0.1) THEN
        amax = ax1(istl)
        RETURN
      END IF
      IF (p .GT. stdp(istl)) THEN
        l1 = istl - 1
        l2 = istl
        alogp = alog(p)
        amax = ax1(l1) + (ax1(l2) - ax1(l1)) * (plin(l1) - alogp) / (plin(l1) - plin(l2))
        RETURN
      END IF
    END DO
    RETURN

  END SUBROUTINE logp_interp

!-----------------------------------------------------------------------------------------------

  SUBROUTINE logp_interh(ax1, h, amax)
!_______________________________________________________________________
! purpose:
!      calculate the h level extreme extreme from ax1
!_______________________________________________________________________
! the input variables
!    ax1 the array of the maximum extreme for standard levels
!    p: observation pressure
!______________________________________________________________________
! the output variables are the same as the input variables
!    amax: the maximum extreme for p level
!______________________________________________________________________
    ! USE standard_pressure

    IMPLICIT NONE
    ! CLASS(GtsQC_t) :: this
! Local arguments:
    REAL            ::  ax1(maxstl)
    REAL            :: h, amax
    INTEGER         :: i, istl, l1, l2

    IF (h .LT. stdh(1)) THEN
      l1 = 1
      l2 = 2
      amax = ax1(l1) + (ax1(l1) - ax1(l2)) * (stdh(l1) - h) / (stdh(l2) - stdh(l1))
      RETURN
    END IF

    DO istl = 1, maxstl
      IF (ABS(h - stdh(istl)) .LT. 0.1) THEN
        amax = ax1(istl)
        RETURN
      END IF
      IF (h .LT. stdh(istl)) THEN
        l1 = istl + 1
        l2 = istl
        amax = ax1(l1) + (ax1(l2) - ax1(l1)) * (stdh(l1) - h) / (stdh(l1) - stdh(l2))
        RETURN
      END IF
    END DO
  END SUBROUTINE logp_interh

  SUBROUTINE comp_extr(hmin1, hmax1, hmin2, hmax2, h, iflag)
!_______________________________________________________________________
! purpose:
!    compare the observation data with it's extreme extreme
!    set the ifalg for qc mark
!_______________________________________________________________________
! method:

!   1. compare the observation data with it's extreme extreme
!      if exceed the extreme 2  ,it is considered wrong
!      if exceed the extreme 1 ,it is considered suspect
!______________________________________________________________________
! the input variables
!    hmin1: the array of the minimum extreme 1 for the observation
!    hmax1: the array of the maximum extreme 1 for the observation
!    hmin2: the array of the minimum extreme 2 for the observation
!    hmax2: the array of the maximum extreme 2 for the observation
!    h   : the observation value
!______________________________________________________________________
! the output variables:
!     iflag:the qc flag for the ob.
!     the other output variables are the same meaning as the input variables
!______________________________________________________________________
    !  USE standard_pressure
    IMPLICIT NONE
    !CLASS(GtsQC_t) :: this
    REAL         :: hmax1, hmax2, hmin1, hmin2, h
    INTEGER      :: iflag

    iflag = 0
    IF (h .LT. hmin1 .OR. h .GT. hmax1) THEN
      iflag = 1
      IF (h .LT. hmin2 .OR. h .GT. hmax2) THEN
        iflag = 2
      END IF
    END IF
    RETURN

  END SUBROUTINE comp_extr

  SUBROUTINE comp_extrw(wmax, w, iflag)
!_______________________________________________________________________
! purpose:
!    compare the wind speed observation data with it's extreme extreme
!    set the ifalg for qc mark
!_______________________________________________________________________
! method:
!   1. compare the observation data with it's extreme extreme
!      if exceed the extreme 1.2 time ,it is considered wrong
!      if exceed the extreme 1.0 time ,it is considered suspect
!______________________________________________________________________
! the input variables
!    wmax: the array of the maximum extreme for the wind speed observation
!    w   : the wind speed observation value
!______________________________________________________________________
! the output variables
!     iflag:the qc flag for the ob.
!     the other output variables are the same as the input variables
!______________________________________________________________________
    IMPLICIT NONE
    !CLASS(GtsQC_t) :: this
    REAL         :: wmax1, wmax, w
    INTEGER      :: iflag
    iflag = 0
    wmax1 = wmax * 1.2
    IF (w .GT. wmax) THEN
      iflag = 1
      IF (w .GT. wmax1) THEN
        iflag = 2
      END IF
    END IF
    RETURN

  END SUBROUTINE comp_extrw

  SUBROUTINE t_td_consist_check(t, td, iq)
!********************************************************************
!     this routine checks temperature and dew-point for mutual
!     consistency, and flag is set.
! 0 is right, 1 is suspect, 2 is wrong
!********************************************************************

    IMPLICIT NONE
    !CLASS(GtsQC_t) :: this
    INTEGER :: iq
    REAL    :: it, t, td

    iq = 9

    IF (ABS(td - robsmissing) .LT. threshold .OR. ABS(t - robsmissing) .LT. threshold) RETURN

    iq = 0
    it = t - td
    IF (it .LT. 0 .OR. it .GE. 55) THEN
      iq = 2
      RETURN
    END IF

    IF (it .GE. 50) THEN
      iq = 1
      RETURN
    END IF

    RETURN
  END SUBROUTINE t_td_consist_check

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  SUBROUTINE wind_consist_check(u, v, iq)

    IMPLICIT NONE
    REAL              :: u, v
    INTEGER           :: iq

!===================================================
    iq = 0
    IF (u == 620) THEN
      iq = 3
      RETURN
    END IF
    IF (ABS(u - robsmissing) < threshold .AND. ABS(v - robsmissing) > threshold) iq = 2
    IF (ABS(u - robsmissing) > threshold .AND. ABS(v - robsmissing) < threshold) iq = 2
    IF (u < 0 .OR. u > 360) THEN
      iq = 2
    ELSE IF (u /= 0 .AND. v == 0) THEN
      iq = 2
    ELSE IF (u == 0 .AND. v /= 0) THEN
      iq = 2
    END IF
    RETURN
  END SUBROUTINE wind_consist_check

  SUBROUTINE ddff_uv(dd, ff, uu, vv, iflag, robsmissing)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Description:
!  purpose: computer the u v from dd ff
!
! method :
!        ANG=450-(DD+180)
!          u=FF*cos(ANG)
!          v=FF*sin(ANG)
!
! Current code owner: NMC
! History:
! Version   Date         Comment
! 1.0       201203    Writen by tianwh
! modified :
! input : dd : wind direction
!         ff : wind speed
! output: uu : wind U
!         vv : wind v
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! use Meteoro_constants
    USE kinds_m, ONLY: r_kind

    IMPLICIT NONE
    !CLASS(GtsQC_t) :: this
    INTEGER :: iflag
    REAL                        :: dd, ff, uu, vv, ang
    REAL(r_kind) :: robsmissing
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    uu = robsmissing
    vv = robsmissing
    IF (ABS(dd - robsmissing) < 0.01) RETURN
    IF (ABS(ff - robsmissing) < 0.01) RETURN
    IF (dd > 360 .AND. dd <= 0) RETURN
    IF (MOD(iflag, 10) == 2) RETURN
    IF (MOD(iflag, 10) == 3) RETURN
    ANG = 270 - dd
    uu = ff * COS(pi180 * ANG)  ! dd --> u
    vv = ff * SIN(pi180 * ANG)  ! ff --> v
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    RETURN
  END SUBROUTINE ddff_uv

  SUBROUTINE FSIGTOBT(SIGP, FSIGH, FSIGZ, POB, POBF, NOL, LP)
    USE Meteoro_constants
    !use obs_size_constants
    IMPLICIT NONE
    !CLASS(GtsQC_t) :: this

    INTEGER        :: NOL, LP, LL, L1, L2, L, LLS
    REAL           :: FSIGZ(LP)
    REAL           :: FSIGH(LP), SIGP(LP)
    REAL           :: POB(NOL), PLIN(LP), POBF(NOL)

    REAL           :: alph, t0, pplin, tplat, tdt, hmod, tstar
    REAL           :: plinc, hr, tmod
    REAL          :: FMIS, tbar, po
    REAL           :: ROG
!========================================================
    ROG = rd / Gravity_s
    FMIS = robsmissing
    POBF(:) = FMIS
    DO L = 1, LP
      PLIN(L) = LOG(SIGP(L))
    END DO

    DO L = 1, NOL
      PO = POB(L)
      PLINC = LOG(PO)
!      print*,'PO',PO,SIGP(LP)

      IF (PO .GT. SIGP(1)) THEN
! -----------------------------------------------------------------
!  OB. LEVEL PRESSURE BELOW THE FIRST MANDATOR LEVEL :P><1000.PHA
!  FIRST, HORIZONTAL INTERPOLATE IN FIRST LEVEL AND SECOND LEVEL BE MADE
!  THEN , VERTICAL EXTRPOLATE WITH HR1 AND HR2
        L1 = 1
        L2 = 2
!        :
!     _______   2
!        :
!        :
!     --------  1
!        :
!        .      CURRENT  OBSERVATION LEVEL
        TMOD = FSIGH(1)
        HMOD = FSIGZ(1)

        TDT = ROG * TMOD * (PO / SIGP(1) - 1)

!     print*,ROG,TMOD,ps,SIGP(1)

!     print*,PO ,TMOD,HMOD,TDT

        IF (TDT .GT. 800) RETURN
        PPLIN = PLINC - LOG(SIGP(1))
        ALPH = 0.0065 * ROG
        HR = TMOD * (1 + 0.5 * ALPH * PPLIN + 1 / 2 * (ALPH * PPLIN)**2 &
                     + 1 / 6 * (ALPH * PPLIN)**3)
        POBF(L) = HR
!       print*,HR

      ELSEIF (PO .LT. SIGP(LP)) THEN
        POBF(L) = FMIS
      ELSE
        DO LL = 1, LP
          LLS = LL
          IF (ABS(PO - SIGP(LL)) .LT. 0.01) THEN
            HR = FSIGH(LLS)
            POBF(L) = HR
            EXIT
          ELSEIF (PO .GT. SIGP(LL)) THEN
!      print*,PO,SIGP(LL)

            L1 = LL
            L2 = LL + 1

! -----------------------------------------------------------------
!  OB. LEVEL PRESSURE P<1000.PHA AND P>P(LP)
!        :
!     _______   2
!        :
!        .      CURRENT  OBSERVATION LEVEL
!        :
!        :
!     --------  1
            HR = FSIGH(L1) - (FSIGH(L1) - FSIGH(L2)) * (PLIN(L1) - PLINC) &
                 / (PLIN(L1) - PLIN(L2))
            POBF(L) = HR
            EXIT
          END IF
        END DO  ! the loop of the model layers
! -----------------------------------------------------------------
      END IF
    END DO  ! the loop of the obs layers
!      print*,POBF(1:LP)
    RETURN
  END SUBROUTINE FSIGTOBT

!--------------------------------------------------------------
  SUBROUTINE FSIGTOBH(SIGP, FSIGH, FSIGT, &
                      POB, ZOB, POBF, NOL, LP)
    USE Meteoro_constants
    ! use obs_size_constants
    IMPLICIT NONE
    !CLASS(GtsQC_t) :: this
    REAL          :: FMIS
    INTEGER       :: NOL, LP
    REAL          :: FSIGH(LP), FSIGT(LP), SIGP(LP), POB(NOL)
    REAL          :: PLIN(LP), ZOB(NOL), POBF(NOL)
    INTEGER       :: L, kp, k
    REAL          :: dh, hmh3, z3, tsig, zsig, alnp, nsig, rog, alnpo
    REAL          :: tbar, tplat, polin, pplin, hmod, alph, t0, pa, po, tstar, relp
    REAL          :: hr
    INTEGER       :: ls, lss
!-------------------------------------------------------------
    ROG = rd / Gravity_s
    FMIS = robsmissing
    POBF = FMIS

    DO L = 1, LP
      PLIN(L) = LOG(SIGP(L))
    END DO

    DO L = 1, NOL
      PO = POB(L)
      IF (PO .LT. SIGP(LP)) RETURN

      IF (PO .GT. SIGP(1)) THEN
        PA = LOG(PO)
        RELP = 1 / (1 + 5 * (PA - PLIN(1))**2)
        IF (RELP .LT. 0.995) CYCLE
        IF (ZOB(L) .GE. 90000.) CYCLE
!     TMOD=FSIGT(LP)
!     HR=FSIGH(LP)+29.713*TMOD*(PLIN(LP)-PA)
!     TMOD=FSIGT(LP)+0.00325*(FSIGH(LP)-HR)
!     HR=FSIGH(LP)+29.713*TMOD*(PLIN(LP)-PA)
!     POBF(L)=HR
!     1998.2.10
        TSTAR = FSIGT(1)
        HMOD = FSIGH(1)
        PPLIN = PA - PLIN(1)
        ALPH = 0.0065 * ROG
        HR = HMOD - ROG * TSTAR * PPLIN * (1 + 0.5 * ALPH * PPLIN + 1 / 6 * ALPH * PPLIN**2)
        POBF(L) = HR
      ELSE
        DO LS = 1, LP
          LSS = LS
          IF ((ABS(PO - SIGP(LS))) .LT. 0.01) THEN
            POBF(L) = FSIGH(LSS)
          ELSE
            POLIN = ALOG(PO)
            CALL PTOSIG(PLIN, FSIGH, FSIGT, LP, POLIN, HR)
            POBF(L) = HR
          END IF
        END DO
      END IF
    END DO
    RETURN

  END SUBROUTINE FSIGTOBH

  SUBROUTINE PTOSIG(ALNP, ZSIG, TSIG, NSIG, ALNPO, ZOB)
    ! use Meteoro_constants
    IMPLICIT NONE
    !CLASS(GtsQC_t) :: this
    INTEGER  :: NSIG, kp, k
    REAL     :: dh, A, B, Z3, HMH3, ZOB, ROG, ALNPO

    REAL     :: ALNP(NSIG), ZSIG(NSIG), TSIG(NSIG)
!H
! ABSTRACT: FINDS VALUES OF TEMPERATURE, WIND AND RELATIVE HUMIDITY
!           FOR EACH OF THE MODEL'S SIGMA COORDINATE LEVELS.
!C           IN ORDER TO CONSTRUCT THE TEMPERATURE, IT FIRST DERIVES
!C           GEOPOTENTIAL HEIGHTS ON THE SIGMA LEVEL INTERFACES
!C           AND THEN DERIVES LEVEL MEAN TEMPERATURE FROM THE
!C           HYDROSTATIC EQUATION
!CH
!      ROG= rd/Gravity_s
!CH
!C     HERE FIND HEIGHTS OF SIGMA SURFACES BY SAME HYDROSTATIC FORMULA
!C     AND FIND WINDS IN THE D SIGMA LAYERS BY INTERPOLATION FROM
!C     MANDATORY LEVELS - INTERPOLATION SAYS WIND COMPONENTS ARE
!C     LINEAR WITH LN(PRESSURE) TO CENTERS OF SIGMA LAYERS.
!C     LAYER DEFINED AS AT LN(PBAR*SIG)
!CH
!C     FINDING PRESSURE LEVELS WHICH BRACKET KS SIG LEVEL
    DO K = 2, NSIG
      KP = K
      IF (ALNP(K) .LT. ALNPO) EXIT
    END DO
    IF (KP .EQ. 2) THEN
      DH = ALNP(KP) - ALNPO
      ZOB = ZSIG(KP) + ROG * TSIG(KP) * DH
      RETURN
    END IF
!C     HEIGHTS - HYDROSTATIC INTERPOLATION
    DH = ALNP(KP - 1) - ALNP(KP)
    A = -(ZSIG(KP - 1) - ZSIG(KP)) / DH
    B = ROG * (TSIG(KP - 1) - TSIG(KP)) / DH
    Z3 = 0.5 * (ZSIG(KP - 1) + ZSIG(KP)) + 0.125 * B * (DH**2)
    HMH3 = ALNPO - 0.5 * (ALNP(KP - 1) + ALNP(KP))
    ZOB = Z3 - (A + 0.5 * B * HMH3) * HMH3
    RETURN
  END SUBROUTINE PTOSIG

  SUBROUTINE complex_dma(iflag)

    IMPLICIT NONE
    INTEGER :: iflag    !  flag before complex_dma and after complex_dma
    INTEGER :: n        !  number of time s to decompse
    INTEGER :: idma     !  final DMA
    INTEGER :: i, ic, m
    idma = 0
    ic = 0
    !        print *,'iflag in ',iflag
    DO i = 2, num_flag
      CALL flag_decompse(iflag, i, m)
      IF (m > 2) CYCLE
      idma = idma + m
      ic = ic + 1
    END DO

    IF (ic < 1) RETURN
    IF (idma > 2) idma = 2
    CALL flag_decompse(iflag, DMA_flag, m)
    IF (m < 3 .AND. idma < m) RETURN
    IF (m < 2) idma = MAX(idma, m)
    CALL flag_replace(iflag, DMA_flag, idma)
    !        print *,'iflag out',iflag
    RETURN
  END SUBROUTINE complex_dma

END MODULE GtsQC_m
