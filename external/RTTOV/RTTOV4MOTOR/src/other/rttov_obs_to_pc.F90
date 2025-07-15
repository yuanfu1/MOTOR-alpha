! Description:
!> @file
!!   Example which converts obs to PC-space.
!
!> @brief
!!   Example program demonstrating how observations may
!!   be converted to PC-space for assimilation/retrieval
!!   applications using PC-RTTOV.
!!
!! @details
!!   For usage details see user guide or run:
!!   $ rttov_obs_to_pc.exe
!!
!!   See the source code for more details: this is intended
!!   as example code which should be modified for your own
!!   purposes.
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 7 December 2016, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, DWD and MeteoFrance.
!
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
PROGRAM rttov_obs_to_pc

#include "throw.h"

  USE rttov_types, ONLY : rttov_options, rttov_coefs
  USE parkind1, ONLY : jplm, jpim, jprb
  USE rttov_getoptions
  USE rttov_unix_env, ONLY : rttov_exit
  USE rttov_math_mod, ONLY : planck

  IMPLICIT NONE

#include "rttov_errorreport.interface"
#include "rttov_read_coefs.interface"
#include "rttov_dealloc_coefs.interface"

  INTEGER(jpim), PARAMETER :: lu = 70

  TYPE(rttov_options)      :: opts
  TYPE(rttov_coefs)        :: coefs

  CHARACTER(256)           :: rtcoef_file
  CHARACTER(256)           :: pccoef_file
  CHARACTER(256)           :: obs_file

  INTEGER(jpim)            :: err
  INTEGER(jpim)            :: i, j
  INTEGER(jpim)            :: nchannels
  INTEGER(jpim)            :: npcscores, ipcbnd, ipcreg
  LOGICAL(jplm)            :: obs_bt
  LOGICAL(jplm)            :: exists

  REAL(jprb), ALLOCATABLE  :: obs(:), pcscores(:), bt_eff(:)

! --End of header--------------------------------------------------------------

TRY

  ! ===========================================================================
  ! This is demonstration code: you would typically want to adapt this for
  ! your own uses.

  ! This code reads in a single observation from an ASCII file: this consists
  ! of radiances or BTs for *all* instrument channels separated by white-space.

  ! You must specify the same PC band (ipcbnd) and regression predictor set
  ! (ipcreg) used in your PC-RTTOV calculations.

  ! The number of PC scores required is specified on the commandline.

  ! The computed PC scores are written to the commandline.
  ! ===========================================================================

  ! ----------------------------------------
  ! Handle program arguments
  ! ----------------------------------------

  CALL initoptions("This program demonstrates conversion of observations to PC space")

  CALL getoption("--rtcoef_file", rtcoef_file, mnd=.TRUE._jplm, use="Input rtcoef file name")
  CALL getoption("--pccoef_file", pccoef_file, mnd=.TRUE._jplm, use="Input pccoef file name")
  CALL getoption("--obs_file", obs_file, mnd=.TRUE._jplm, &
      use="Input ASCII file containing observations (white-space separated)")

  CALL getoption("--ipcbnd", ipcbnd, mnd=.TRUE._jplm, use="PC-RTTOV band (usually 1)")
  CALL getoption("--ipcreg", ipcreg, mnd=.TRUE._jplm, use="PC-RTTOV predictor set index")
  CALL getoption("--npcscores", npcscores, mnd=.TRUE._jplm, use="Number of PC scores to output")

  obs_bt = .FALSE.
  CALL getoption( "--obs_bt", obs_bt, use="If present obs units are Kelvin, otherwise mW/m-2/sr-1/cm-1")

  CALL checkoptions()

  INQUIRE(FILE=rtcoef_file, EXIST=exists)
  IF (.NOT. exists) THEN
    PRINT *, 'Cannot find rtcoef file: '//TRIM(rtcoef_file)
    STOP
  ENDIF
  INQUIRE(FILE=pccoef_file, EXIST=exists)
  IF (.NOT. exists) THEN
    PRINT *, 'Cannot find pccoef file: '//TRIM(pccoef_file)
    STOP
  ENDIF
  INQUIRE(FILE=obs_file, EXIST=exists)
  IF (.NOT. exists) THEN
    PRINT *, 'Cannot find observations file: '//TRIM(obs_file)
    STOP
  ENDIF

  ! ----------------------------------------
  ! Set options
  ! ----------------------------------------

  opts%rt_ir%pc%addpc = .TRUE.
  opts%rt_ir%pc%ipcbnd = ipcbnd
  opts%rt_ir%pc%ipcreg = ipcreg

  ! ----------------------------------------
  ! Read coefficients and observations
  ! ----------------------------------------

  CALL rttov_read_coefs(err, coefs, opts, file_coef=rtcoef_file, file_pccoef=pccoef_file)
  THROWM(err .NE. 0, 'Error reading coefficients')

  nchannels = coefs%coef%fmv_chn
  IF (nchannels /= MAXVAL(coefs%coef%ff_ori_chn)) THEN
    err = errorstatus_fatal
    THROWM(err .NE. 0, 'The full coefficient file containing all instrument channels is required')
  ENDIF

  OPEN(lu, FILE=obs_file, IOSTAT=err)
  THROWM(err .NE. 0, 'Error opening observations file')

  ALLOCATE(obs(nchannels))
  READ(lu, *, IOSTAT=err) obs
  THROWM(err .NE. 0, 'Error reading observations')

  CLOSE(lu)

  ! ----------------------------------------
  ! Convert BTs to radiances (if necessary)
  ! ----------------------------------------

  ! Coefficient structure contains correction scaling and offset for BTs
  ! and constants for use in Planck function for each channel

  IF (obs_bt) THEN
    ALLOCATE(bt_eff(nchannels))
    DO j = 1, nchannels
      bt_eff(j) = coefs%coef%ff_bco(j) + coefs%coef%ff_bcs(j) * obs(j)
      CALL planck(coefs%coef%planck1(j), coefs%coef%planck2(j), bt_eff(j), obs(j))
    ENDDO
    DEALLOCATE(bt_eff)
  ENDIF

  ! ----------------------------------------
  ! Convert observations to PCs
  ! ----------------------------------------

  ALLOCATE(pcscores(npcscores))
  pcscores = 0._jprb

  DO i = 1, npcscores
    DO j = 1, nchannels
      pcscores(i)= pcscores(i) + &
        coefs%coef_pccomp%eigen(ipcbnd)%eigenvectors(j,i) * &
        obs(j) / coefs%coef_pccomp%noise_in(j)
    ENDDO
  ENDDO

  ! ----------------------------------------
  ! Write out the PC scores
  ! ----------------------------------------

  WRITE(*,'(a)') 'PC scores:'
  WRITE(*,'(6e13.5)') pcscores

  ! ----------------------------------------
  ! Tidy up
  ! ----------------------------------------

  IF (ALLOCATED(obs))      DEALLOCATE(obs)
  IF (ALLOCATED(pcscores)) DEALLOCATE(pcscores)

  CALL rttov_dealloc_coefs(err, coefs)
  THROWM(err .NE. 0, 'Error deallocating coefficients')

PCATCH
END PROGRAM rttov_obs_to_pc
