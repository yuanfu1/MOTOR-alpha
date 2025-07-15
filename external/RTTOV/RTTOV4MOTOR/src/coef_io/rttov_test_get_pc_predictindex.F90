! Description:
!> @file
!!   Test program for running the rttov_get_pc_predictindex
!!   subroutine on the commandline.
!
!> @brief
!!   Test program for running the rttov_get_pc_predictindex
!!   subroutine on the commandline.
!!
!! @details
!!   This prints out the predictor channels for a specified PC band
!!   (ipcbnd) and predictor regression set (ipcreg) for the given
!!   PC-RTTOV coefficient file.
!!
!!   For usage details run:
!!   $ rttov_test_get_pc_predictindex.exe \-\-help
!!
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
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!
PROGRAM rttov_test_get_pc_predictindex

#include "throw.h"

  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_options
  USE rttov_getoptions, ONLY : initoptions, getoption, checkoptions
  USE rttov_unix_env, ONLY : rttov_exit

  IMPLICIT NONE

#include "rttov_errorreport.interface"
#include "rttov_get_pc_predictindex.interface"

  TYPE(rttov_options)  :: opts
  INTEGER(KIND=jpim)   :: err
  INTEGER(KIND=jpim)   :: ipcbnd, ipcreg
  CHARACTER(LEN=256)   :: f_pccoef_in = ""
  INTEGER(KIND=jpim), POINTER ::predictindex(:)

  !- End of header --------------------------------------------------------
  TRY

  CALL initoptions( "This program writes out the selected PC-RTTOV predictor channel set")

  CALL getoption( "--pccoef-in",  f_pccoef_in, mnd=.TRUE._jplm, use="PC-RTTOV coefficient file")
  CALL getoption( "--ipcbnd",  ipcbnd, mnd=.TRUE._jplm, use="PC-RTTOV band (usually 1)")
  CALL getoption( "--ipcreg",  ipcreg, mnd=.TRUE._jplm, use="PC-RTTOV predictor set index")

  CALL checkoptions()

  opts%rt_ir%pc%addpc  = f_pccoef_in .NE. ""
  opts%rt_ir%pc%ipcbnd = ipcbnd
  opts%rt_ir%pc%ipcreg = ipcreg

  CALL rttov_get_pc_predictindex( &
            & err,           &
            & opts,          &
            & predictindex,  &
            & file_pccoef   = f_pccoef_in)
  THROW(err.NE.0)

  WRITE(*,*) "Number of channels ", SIZE(predictindex)
  WRITE(*,*) predictindex
  DEALLOCATE(predictindex)

  PCATCH
END PROGRAM rttov_test_get_pc_predictindex
