! Description:
!> @file
!!   Extract data for given channel list from a cloud coefficients
!!   structure.
!
!> @brief
!!   Extract data for given channel list from a cloud coefficients
!!   structure.
!!
!! @details
!!   This is used by HDF5 I/O code to read in a subset of channels from a
!!   coefficient file. The first coef arguments contain the coefficients
!!   from the file. The second arguments are uninitialised structures
!!   which contain the extracted coefficients on exit.
!!
!! @param[out]     err           status on exit
!! @param[in]      coef_scatt1   input rttov_coef_scatt structure read from file
!! @param[in,out]  coef_scatt2   output rttov_coef_scatt structure, uninitialised on entry
!! @param[in]      channels      list of channels to extract
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
!    Copyright 2018, EUMETSAT, All Rights Reserved.
!
SUBROUTINE rttov_channel_extract_sccldcoef(err, coef_scatt1, coef_scatt2, channels)
!INTF_OFF
#include "throw.h"
  USE rttov_cldaer_io_mod, ONLY : channel_extract_optp
!INTF_ON
  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY : rttov_coef_scatt

  IMPLICIT NONE

  INTEGER(jpim),          INTENT(OUT)   :: err
  TYPE(rttov_coef_scatt), INTENT(IN)    :: coef_scatt1
  TYPE(rttov_coef_scatt), INTENT(INOUT) :: coef_scatt2
  INTEGER(jpim),          INTENT(IN)    :: channels(:)
!INTF_END

#include "rttov_errorreport.interface"

  TRY

  CALL channel_extract_optp(err, coef_scatt1%optp_wcl_opac, coef_scatt2%optp_wcl_opac, channels)
  THROW(err.NE.0)

  CALL channel_extract_optp(err, coef_scatt1%optp_wcl_deff, coef_scatt2%optp_wcl_deff, channels)
  THROW(err.NE.0)

  CALL channel_extract_optp(err, coef_scatt1%optp_icl_baum, coef_scatt2%optp_icl_baum, channels)
  THROW(err.NE.0)

  CATCH
END SUBROUTINE rttov_channel_extract_sccldcoef
