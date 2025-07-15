! Description:
!> @file
!!   Allocate/deallocate pccomp structure for PC-RTTOV and HTFRTC.
!
!> @brief
!!   Allocate/deallocate pccomp structure for PC-RTTOV and HTFRTC.
!!
!! @details
!!   The pccomp structure contains the output PC scores and
!!   reconstructed radiances, BTs for the PC-RTTOV direct model,
!!   the corresponding perturbations for the TL model, and the
!!   input gradients and perturbations for the AD and K models.
!!
!!   For HTFRTC there are additional outputs for clear, cloudy
!!   and overcast radiances and BTs. When allocating this structure
!!   for HTFRTC simulations, the relevant options must be set in the
!!   options structure opts.
!!
!! @param[out]    err            status on exit
!! @param[in,out] pccomp         PC scores and radiances for PC-RTTOV
!! @param[in]     npcscores      number of PC scores to simulate for all profiles
!! @param[in]     asw            1_jpim => allocate; 0_jpim => deallocate
!! @param[in]     init           set .TRUE. to initialise newly allocated structures, optional
!! @param[in]     nchannels_rec  total number of reconstructed radiances for all profiles, optional
!! @param[in]     opts           RTTOV options, optional (required for HTFRTC)
!! @param[in]     nlevels        number of pressure levels in the input profiles
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
SUBROUTINE rttov_alloc_pccomp( &
              err,           &
              pccomp,        &
              npcscores,     &
              asw,           &
              init,          &
              nchannels_rec, &
              opts,          &
              nlevels)
!INTF_OFF
#include "throw.h"
!INTF_ON

  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_pccomp, rttov_options

  IMPLICIT NONE

  INTEGER(KIND=jpim),  INTENT(OUT)            :: err
  TYPE(rttov_pccomp),  INTENT(INOUT)          :: pccomp
  INTEGER(KIND=jpim),  INTENT(IN)             :: npcscores
  INTEGER(KIND=jpim),  INTENT(IN)             :: asw
  LOGICAL(KIND=jplm),  INTENT(IN),   OPTIONAL :: init
  INTEGER(KIND=jpim),  INTENT(IN),   OPTIONAL :: nchannels_rec
  TYPE(rttov_options), INTENT(IN),   OPTIONAL :: opts
  INTEGER(KIND=jpim),  INTENT(IN),   OPTIONAL :: nlevels
!INTF_END

#include "rttov_errorreport.interface"
#include "rttov_init_pccomp.interface"

  LOGICAL(KIND=jplm) :: init1, htfrtc
!- End of header --------------------------------------------------------

  TRY
  IF (asw == 1) THEN
    init1 = .FALSE.
    IF (PRESENT(init)) init1 = init
    htfrtc = .FALSE.
    IF (PRESENT(opts)) htfrtc = opts%htfrtc_opts%htfrtc
    IF (htfrtc) THEN
      IF (opts%htfrtc_opts%overcast .AND. .NOT. PRESENT(nlevels)) THEN
        err = errorstatus_fatal
        THROWM(err .NE. 0, "nlevels argument is mandatory for HTFRTC overcast simulations")
      ENDIF
    ENDIF

    CALL nullify_struct(pccomp)

    IF (htfrtc) THEN
      ! Arrays for HTFRTC
      ALLOCATE (pccomp%total_pcscores(npcscores), pccomp%clear_pcscores(npcscores), STAT = err)
      THROWM(err .NE. 0, "allocation of pcscores")

      IF (PRESENT(nchannels_rec)) THEN
        IF (nchannels_rec > 0) THEN
          ALLOCATE (pccomp%bt_clear_pccomp(nchannels_rec), pccomp%clear_pccomp(nchannels_rec), STAT = err)
          THROWM(err .NE. 0, "allocation of bt_clear_pccomp, clear_pccomp")
          ALLOCATE (pccomp%bt_pccomp(nchannels_rec), pccomp%total_pccomp(nchannels_rec), STAT = err)
          THROWM(err .NE. 0, "allocation of bt_pccomp, total_pccomp")
        ENDIF
      ENDIF

      IF (opts%htfrtc_opts%overcast) THEN
        ALLOCATE (pccomp%overcast_pcscores(nlevels-1,npcscores), STAT = err)
        THROWM(err .NE. 0, "allocation of overcast_pcscores")
        IF (PRESENT(nchannels_rec)) THEN
          IF (nchannels_rec > 0) THEN
            ALLOCATE (pccomp%overcast_pccomp(nlevels-1,nchannels_rec), STAT = err)
            THROWM(err .NE. 0, "allocation of overcast_pccomp")
          ENDIF
        ENDIF
      ENDIF

      IF (opts%htfrtc_opts%simple_cloud) THEN
        ALLOCATE (pccomp%cloudy_pcscores(npcscores), STAT = err)
        THROWM(err .NE. 0, "allocation of cloudy_pcscores")
        IF (PRESENT(nchannels_rec)) THEN
          IF (nchannels_rec > 0) THEN
            ALLOCATE (pccomp%cloudy_pccomp(nchannels_rec),  STAT = err)
            THROWM(err .NE. 0, "allocation of cloudy_pccomp")
          ENDIF
        ENDIF
      ENDIF
    ELSE
      ! Arrays for PC-RTTOV
      ALLOCATE (pccomp%total_pcscores(npcscores), STAT = err)
      THROWM(err .NE. 0, "allocation of total_pcscores")
      IF (PRESENT(nchannels_rec)) THEN
        IF (nchannels_rec > 0) THEN
          ALLOCATE (pccomp%bt_pccomp(nchannels_rec), pccomp%total_pccomp(nchannels_rec), STAT = err)
          THROWM(err .NE. 0, "allocation of bt_pccomp, total_pccomp")
        ENDIF
      ENDIF
    ENDIF
    IF (init1) CALL rttov_init_pccomp(pccomp)
  ENDIF

  IF (asw == 0) THEN
    IF (ASSOCIATED(pccomp%clear_pcscores)) THEN
      DEALLOCATE (pccomp%clear_pcscores, STAT = err)
      THROWM(err .NE. 0, "deallocation of clear_pcscores")
    ENDIF
    IF (ASSOCIATED(pccomp%total_pcscores)) THEN
      DEALLOCATE (pccomp%total_pcscores, STAT = err)
      THROWM(err .NE. 0, "deallocation of total_pcscores")
    ENDIF
    IF (ASSOCIATED(pccomp%overcast_pcscores)) THEN
      DEALLOCATE (pccomp%overcast_pcscores, STAT = err)
      THROWM(err .NE. 0, "deallocation of overcast_pcscores")
    ENDIF
    IF (ASSOCIATED(pccomp%cloudy_pcscores)) THEN
      DEALLOCATE (pccomp%cloudy_pcscores, STAT = err)
      THROWM(err .NE. 0, "deallocation of cloudy_pcscores")
    ENDIF
    IF (ASSOCIATED(pccomp%clear_pccomp)) THEN
      DEALLOCATE (pccomp%clear_pccomp, STAT = err)
      THROWM(err .NE. 0, "deallocation of clear_pccomp")
    ENDIF
    IF (ASSOCIATED(pccomp%total_pccomp)) THEN
      DEALLOCATE (pccomp%total_pccomp, STAT = err)
      THROWM(err .NE. 0, "deallocation of total_pccomp")
    ENDIF
    IF (ASSOCIATED(pccomp%overcast_pccomp)) THEN
      DEALLOCATE (pccomp%overcast_pccomp, STAT = err)
      THROWM(err .NE. 0, "deallocation of overcast_pccomp")
    ENDIF
    IF (ASSOCIATED(pccomp%cloudy_pccomp)) THEN
      DEALLOCATE (pccomp%cloudy_pccomp, STAT = err)
      THROWM(err .NE. 0, "deallocation of cloudy_pccomp")
    ENDIF
    IF (ASSOCIATED(pccomp%bt_pccomp)) THEN
      DEALLOCATE (pccomp%bt_pccomp, STAT = err)
      THROWM(err .NE. 0, "deallocation of bt_pccomp")
    ENDIF
    IF (ASSOCIATED(pccomp%bt_clear_pccomp)) THEN
      DEALLOCATE (pccomp%bt_clear_pccomp, STAT = err)
      THROWM(err .NE. 0, "deallocation of bt_clear_pccomp")
    ENDIF
    CALL nullify_struct(pccomp)
  ENDIF
  CATCH
CONTAINS
  SUBROUTINE nullify_struct(pccomp)
    TYPE(rttov_pccomp), INTENT(INOUT) :: pccomp
    NULLIFY(pccomp%clear_pcscores,    &
            pccomp%total_pcscores,    &
            pccomp%overcast_pcscores, &
            pccomp%cloudy_pcscores,   &
            pccomp%clear_pccomp,      &
            pccomp%total_pccomp,      &
            pccomp%overcast_pccomp,   &
            pccomp%cloudy_pccomp,     &
            pccomp%bt_pccomp,         &
            pccomp%bt_clear_pccomp)
  END SUBROUTINE nullify_struct
END SUBROUTINE 
