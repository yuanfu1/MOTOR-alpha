!!--------------------------------------------------------------------------------------------------
!! PROJECT           : Multi-Grid generation: test edge function, normal and tangential derivative
!!                    calculation.
!! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring 
!!                    Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
!! AUTOHR(S)         : Yuanfu Xie
!! VERSION           : Beta 0.0
!! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2024/12/24, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------
PROGRAM testEdgeSetting
  USE kinds_m, ONLY: i_kind, r_kind

  IMPLICIT NONE

  INTEGER(i_kind) :: ij(2),nm(2),numStencil,info,idxStencil(7),istatus
  REAL(r_kind) :: coeff_func(7),coeff_norm(7),coeff_tang(7),ll(2),dll(2)

  numStencil = 7
  nm = 41
  ij = nm
  ! ij(1) = nm(1)
  ! ij = 2
  ij(2) = 1
  ll = 0.0D0
  dll = 1.0D0
  DO info =1,4
    CALL edgeSetting_s(ij,nm,ll,dll,numStencil,info,idxStencil,coeff_func,coeff_norm,coeff_tang,istatus)
    WRITE(*,10) info,ij
  10 FORMAT('Cell info: ',I2,' ij: ',2I3)
    WRITE(*,1) idxStencil
1   FORMAT('Stencil: ',3I5)
    WRITE(*,2) coeff_func
2   FORMAT('Func cf: ',3E12.4)
    WRITE(*,3) coeff_norm
3   FORMAT('Norm cf: ',3E12.4)
    WRITE(*,4) coeff_tang
4   FORMAT('Tang cf: ',3E12.4)
  END DO
END PROGRAM testEdgeSetting
