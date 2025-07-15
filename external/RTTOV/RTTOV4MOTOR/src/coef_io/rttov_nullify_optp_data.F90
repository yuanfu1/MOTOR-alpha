! Description:
!> @file
!!   Nullify/zero a VIS/IR aerosol/cloud optical property data structure.
!
!> @brief
!!   Nullify/zero a VIS/IR aerosol/cloud optical property data structure.
!!
!!
!! @param[in,out]  data  the cloud/aerosol optical property data structure to nullify/zero
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
SUBROUTINE rttov_nullify_optp_data(data)
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  USE rttov_types, ONLY : rttov_optp_data
  IMPLICIT NONE
  TYPE(rttov_optp_data), INTENT(INOUT) :: data
!INTF_END

  data%name = ''
  data%confac = 0._jprb
  data%nrelhum = 0
  data%ndeff = 0

  NULLIFY(data%relhum, &
          data%deff, &
          data%abs, &
          data%sca, &
          data%bpr, &
          data%nmom, &
          data%legcoef, &
          data%pha)

END SUBROUTINE rttov_nullify_optp_data
