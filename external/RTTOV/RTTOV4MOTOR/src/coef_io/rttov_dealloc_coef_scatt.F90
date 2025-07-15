! Description:
!> @file
!!   Deallocate a VIS/IR cloud/aerosol optical properties structure.
!!   This should usually be called via rttov_dealloc_coefs.
!
!> @brief
!!   Deallocate a VIS/IR cloud/aerosol optical properties structure.
!!   This should usually be called via rttov_dealloc_coefs.
!!
!! @param[out]     err          status on exit
!! @param[in,out]  coef_scatt   the cloud/aerosol optical properties structure to deallocate
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
SUBROUTINE rttov_dealloc_coef_scatt(err, coef_scatt)
!INTF_OFF
#include "throw.h"
!INTF_ON
  USE rttov_types, ONLY : rttov_coef_scatt
  USE parkind1, ONLY : jpim

  IMPLICIT NONE

  INTEGER(jpim),          INTENT(OUT)   :: err
  TYPE(rttov_coef_scatt), INTENT(INOUT) :: coef_scatt
!INTF_END
#include "rttov_errorreport.interface"
#include "rttov_nullify_coef_scatt.interface"
#include "rttov_alloc_phfn_int.interface"

!- End of header --------------------------------------------------------
  TRY

  CALL deallocate_optp(err, coef_scatt%optp_aer)
  THROW(err.NE.0)
  CALL deallocate_optp(err, coef_scatt%optp_wcl_opac)
  THROW(err.NE.0)
  CALL deallocate_optp(err, coef_scatt%optp_wcl_deff)
  THROW(err.NE.0)
  CALL deallocate_optp(err, coef_scatt%optp_icl_baum)
  THROW(err.NE.0)
  CALL deallocate_optp_baran(err, coef_scatt%optp_icl_baran2014)
  THROW(err.NE.0)
  CALL deallocate_optp_baran(err, coef_scatt%optp_icl_baran2018)
  THROW(err.NE.0)

  CALL rttov_nullify_coef_scatt(coef_scatt)

  CATCH
CONTAINS

  SUBROUTINE deallocate_optp_data(err, data)
    USE rttov_types, ONLY : rttov_optp_data

    INTEGER(jpim),         INTENT(OUT)   :: err
    TYPE(rttov_optp_data), INTENT(INOUT) :: data

    TRY

    IF (ASSOCIATED(data%relhum))  DEALLOCATE(data%relhum, STAT=err)
    THROW(err.NE.0)
    IF (ASSOCIATED(data%deff))    DEALLOCATE(data%deff, STAT=err)
    THROW(err.NE.0)
    IF (ASSOCIATED(data%abs))     DEALLOCATE(data%abs, STAT=err)
    THROW(err.NE.0)
    IF (ASSOCIATED(data%sca))     DEALLOCATE(data%sca, STAT=err)
    THROW(err.NE.0)
    IF (ASSOCIATED(data%bpr))     DEALLOCATE(data%bpr, STAT=err)
    THROW(err.NE.0)
    IF (ASSOCIATED(data%nmom))    DEALLOCATE(data%nmom, STAT=err)
    THROW(err.NE.0)
    IF (ASSOCIATED(data%legcoef)) DEALLOCATE(data%legcoef, STAT=err)
    THROW(err.NE.0)
    IF (ASSOCIATED(data%pha))     DEALLOCATE(data%pha, STAT=err)
    THROW(err.NE.0)

    CATCH
  END SUBROUTINE deallocate_optp_data

  SUBROUTINE deallocate_optp(err, optp)
    USE rttov_types, ONLY : rttov_optp

    INTEGER(jpim),    INTENT(INOUT) :: err
    TYPE(rttov_optp), INTENT(INOUT) :: optp

    INTEGER(jpim) :: i

    TRY

    IF (ASSOCIATED(optp%chan_pha))       DEALLOCATE(optp%chan_pha, STAT=err)
    THROW(err.NE.0)
    IF (ASSOCIATED(optp%chan_pha_index)) DEALLOCATE(optp%chan_pha_index, STAT=err)
    THROW(err.NE.0)
    IF (ASSOCIATED(optp%data)) THEN
      DO i = 1, optp%ntypes
        CALL deallocate_optp_data(err, optp%data(i))
        THROW(err.NE.0)
      ENDDO

      DEALLOCATE(optp%data, STAT=err)
      THROW(err.NE.0)
    ENDIF
    IF (ASSOCIATED(optp%phangle)) THEN
      CALL rttov_alloc_phfn_int(err, optp%phangle, optp%phfn_int, 0_jpim)
      THROW(err.NE.0)

      DEALLOCATE(optp%phangle, STAT=err)
      THROW(err.NE.0)
    ENDIF

    CATCH
  END SUBROUTINE deallocate_optp

  SUBROUTINE deallocate_optp_baran(err, optp)
    USE rttov_types, ONLY : rttov_optp_baran
    USE rttov_const, ONLY : phangle_hires

    INTEGER(jpim),          INTENT(INOUT) :: err
    TYPE(rttov_optp_baran), INTENT(INOUT) :: optp

    TRY

    IF (ASSOCIATED(optp%iwn))    DEALLOCATE(optp%iwn, STAT=err)
    THROW(err.NE.0)
    IF (ASSOCIATED(optp%jwn))    DEALLOCATE(optp%jwn, STAT=err)
    THROW(err.NE.0)
    IF (ASSOCIATED(optp%dx_dwn)) DEALLOCATE(optp%dx_dwn, STAT=err)
    THROW(err.NE.0)
    IF (ASSOCIATED(optp%q))      DEALLOCATE(optp%q, STAT=err)
    THROW(err.NE.0)
    IF (ASSOCIATED(optp%w))      DEALLOCATE(optp%w, STAT=err)
    THROW(err.NE.0)
    CALL rttov_alloc_phfn_int(err, phangle_hires, optp%phfn_int, 0_jpim)
    THROW(err.NE.0)

    CATCH
  END SUBROUTINE deallocate_optp_baran

END SUBROUTINE rttov_dealloc_coef_scatt
