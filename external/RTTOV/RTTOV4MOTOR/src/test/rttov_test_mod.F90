! Description:
!> @file
!!   Helper subroutines for the test suite for writing output data to files
!!   and converting profile units.
!
!> @brief
!!   Helper subroutines for the test suite for writing output data to files
!!   and converting profile units.
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
!    Copyright 2017, EUMETSAT, All Rights Reserved.
!
MODULE rttov_test_mod

#include "throw.h"

  USE parkind1, ONLY : jprb, jpim

  USE rttov_const, ONLY : &
      gas_id_watervapour, &
      gas_id_ozone,       &
      gas_id_co2,         &
      gas_id_n2o,         &
      gas_id_co,          &
      gas_id_ch4,         &
      gas_id_so2,         &
      mair,               &
      gas_mass

  USE rttov_types, ONLY :  &
      rttov_profile,       &
      rttov_profile_cloud, &
      rttov_opt_param,     &
      rttov_radiance,      &
      rttov_radiance2,     &
      rttov_reflectivity,  &
      rttov_transmission,  &
      rttov_pccomp,        &
      rttov_scatt_emis_retrieval_type

  USE rttov_lun, ONLY : rttov_get_lun, rttov_put_lun

  IMPLICIT NONE

#include "rttov_errorreport.interface"

  PRIVATE
  PUBLIC :: setGformat, &
            printprofiles, printcldprofiles, printoptparam, &
            printradiance, printtransmission, &
            printreflectivity, printemisterms, &
            printsurfvar, printpcscores, &
            ppmvdry2ppmvwet, ppmvdry2kgkgwet, ppmvwet2ppmvdry, &
            ppmvwet2kgkgwet, kgkgwet2ppmvdry, kgkgwet2ppmvwet

  CHARACTER(LEN=8)            :: Gformat = 'G15.6E3'
  CHARACTER(LEN=*), PARAMETER :: Iformat = 'I10,5X'

CONTAINS

  SUBROUTINE setGformat(realprec)
    INTEGER(jpim), INTENT(IN) :: realprec
    IF (realprec < 1) THEN
      Gformat = 'G15.6E3'
    ELSEIF (realprec < 10) THEN
      WRITE(Gformat,'(A,I2,A,I1,A)') 'G', realprec+9, '.', realprec, 'E3'
    ELSE
      WRITE(Gformat,'(A,I2,A,I2)') 'E', realprec+9, '.', realprec
    ENDIF
  END SUBROUTINE setGformat

  SUBROUTINE printprofiles(err, profiles, name, f)
    USE rttov_chain, ONLY : chain, chain_rttov_profile, print_chain, free_chain

    INTEGER(jpim),       INTENT(OUT)        :: err
    TYPE(rttov_profile), INTENT(IN), TARGET :: profiles(:) ! target because of chain_rttov_profile
    CHARACTER(LEN=*),    INTENT(IN)         :: name, f

    TYPE(chain) :: chain_prof
    INTEGER(jpim) :: file_id

    TRY

    CALL rttov_get_lun(file_id)
    OPEN(file_id, file = f, form = 'formatted', iostat = err)
    THROWM(err.NE.0,"Cannot open "//TRIM(f))

    CALL chain_rttov_profile(err, chain_prof, name, a1 = profiles)
    THROWM(err.NE.0,"Cannot chain profiles")

    CALL print_chain(INT(file_id,jpim), chain_prof, TRIM(Gformat))

    CALL free_chain(chain_prof)

    CLOSE(file_id, iostat = err)
    THROWM(err.NE.0,"Cannot close "//TRIM(f))
    CALL rttov_put_lun(file_id)

    CATCH
  END SUBROUTINE

  SUBROUTINE printcldprofiles(err, cld_profiles, name, f)
    USE rttov_chain, ONLY : chain, chain_rttov_profile_cloud, print_chain, free_chain

    INTEGER(jpim),             INTENT(OUT)        :: err
    TYPE(rttov_profile_cloud), INTENT(IN), TARGET :: cld_profiles(:) ! target because of chain_rttov_cld_profile
    CHARACTER(LEN=*),          INTENT(IN)         :: name, f

    TYPE(chain) :: chain_prof
    INTEGER(jpim) :: file_id

    TRY

    CALL rttov_get_lun(file_id)
    OPEN(file_id, file = f, form = 'formatted', iostat = err)
    THROWM(err.NE.0,"Cannot open "//TRIM(f))

    CALL chain_rttov_profile_cloud(err, chain_prof, name, a1 = cld_profiles)
    THROWM(err.NE.0,"Cannot chain cld_profiles")

    CALL print_chain(INT(file_id,jpim), chain_prof, TRIM(Gformat))

    CALL free_chain(chain_prof)

    CLOSE(file_id, iostat = err)
    THROWM(err.NE.0,"Cannot close "//TRIM(f))
    CALL rttov_put_lun(file_id)

    CATCH
  END SUBROUTINE

  SUBROUTINE printoptparam(err, opt_param, name, f)
    USE rttov_chain, ONLY : print_array

    INTEGER(jpim),         INTENT(OUT) :: err
    TYPE(rttov_opt_param), INTENT(IN)  :: opt_param
    CHARACTER(LEN=*),      INTENT(IN)  :: name, f

    INTEGER(jpim) :: file_id

    TRY

    CALL rttov_get_lun(file_id)
    OPEN(file_id, file = f, form = 'formatted', iostat = err)
    THROWM(err.NE.0,"Cannot open "//TRIM(f))

    WRITE(file_id, '(A," = (")') name//'%ABS'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = opt_param%abs)
    WRITE(file_id, *) ')'

    WRITE(file_id, '(A," = (")') name//'%SCA'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = opt_param%sca)
    WRITE(file_id, *) ')'

    IF (ANY(opt_param%bpr /= 0._jprb)) THEN
      WRITE(file_id, '(A," = (")') name//'%BPR'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = opt_param%bpr)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ANY(opt_param%legcoef /= 0._jprb)) THEN
      WRITE(file_id, '(A," = (")') name//'%LEGCOEF'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A3 = opt_param%legcoef)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ANY(opt_param%pha /= 0._jprb)) THEN
      WRITE(file_id, '(A," = (")') name//'%PHA'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A3 = opt_param%pha)
      WRITE(file_id, *) ')'
    ENDIF

    CLOSE(file_id, iostat = err)
    THROWM(err.NE.0,"Cannot close "//TRIM(f))
    CALL rttov_put_lun(file_id)

    CATCH
  END SUBROUTINE

  SUBROUTINE printradiance(err, radiance, name, f, radiance2)
    USE rttov_chain, ONLY : print_array

    INTEGER(jpim),         INTENT(OUT)          :: err
    TYPE(rttov_radiance),  INTENT(IN)           :: radiance
    CHARACTER(LEN=*),      INTENT(IN)           :: name, f
    TYPE(rttov_radiance2), INTENT(IN), OPTIONAL :: radiance2

    INTEGER(jpim) :: file_id

    TRY

    CALL rttov_get_lun(file_id)
    OPEN(file_id, file = f, form = 'formatted', iostat = err)
    THROWM(err.NE.0,"Cannot open "//TRIM(f))

    WRITE(file_id, '(A," = (")') name//'%TOTAL'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = radiance%total)
    WRITE(file_id, *) ')'

    IF (ANY(radiance%bt /= 0._jprb)) THEN
      WRITE(file_id, '(A," = (")') name//'%BT'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = radiance%bt)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ANY(radiance%refl /= 0._jprb)) THEN
      WRITE(file_id, '(A," = (")') name//'%REFL'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = radiance%refl)
      WRITE(file_id, *) ')'
    ENDIF

    WRITE(file_id, '(A," = (")') name//'%CLEAR'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = radiance%clear)
    WRITE(file_id, *) ')'

    IF (ANY(radiance%bt_clear /= 0._jprb)) THEN
      WRITE(file_id, '(A," = (")') name//'%BT_CLEAR'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = radiance%bt_clear)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ANY(radiance%refl_clear /= 0._jprb)) THEN
      WRITE(file_id, '(A," = (")') name//'%REFL_CLEAR'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = radiance%refl_clear)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ANY(radiance%overcast /= 0._jprb)) THEN
      WRITE(file_id, '(A," = (")') name//'%OVERCAST'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = radiance%overcast)
      WRITE(file_id, *) ')'
    ENDIF

    WRITE(file_id, '(A," = (")') name//'%CLOUDY100%'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = radiance%cloudy)
    WRITE(file_id, *) ')'

    WRITE(file_id, '(A," = (")') name//'%QUALITY'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Iformat) // '))', AI1 = radiance%quality)
    WRITE(file_id, *) ')'

    WRITE(file_id, '(A," = (")') name//'%GEOMETRIC_HEIGHT'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = radiance%geometric_height)
    WRITE(file_id, *) ')'

    IF (PRESENT(radiance2)) THEN

      WRITE(file_id, '(A," = (")') name//'%UPCLEAR'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = radiance2%upclear)
      WRITE(file_id, *) ')'

      WRITE(file_id, '(A," = (")') name//'%DNCLEAR'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = radiance2%dnclear)
      WRITE(file_id, *) ')'

      WRITE(file_id, '(A," = (")') name//'%REFLDNCLEAR'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = radiance2%refldnclear)
      WRITE(file_id, *) ')'

      WRITE(file_id, '(A," = (")') name//'%UP'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = radiance2%up)
      WRITE(file_id, *) ')'

      WRITE(file_id, '(A," = (")') name//'%DOWN'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = radiance2%down)
      WRITE(file_id, *) ')'

      WRITE(file_id, '(A," = (")') name//'%SURF'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = radiance2%surf)
      WRITE(file_id, *) ')'

    ENDIF

    CLOSE(file_id, iostat = err)
    THROWM(err.NE.0,"Cannot close "//TRIM(f))
    CALL rttov_put_lun(file_id)

    CATCH
  END SUBROUTINE

  SUBROUTINE printtransmission(err, transmission, name, f)
    USE rttov_chain, ONLY : print_array

    INTEGER(jpim),            INTENT(OUT) :: err
    TYPE(rttov_transmission), INTENT(IN) :: transmission
    CHARACTER(LEN=*),         INTENT(IN) :: name, f

    INTEGER(jpim) :: file_id

    TRY

    CALL rttov_get_lun(file_id)
    OPEN(file_id, file = f, form = 'formatted', iostat = err)
    THROWM(err.NE.0,"Cannot open "//TRIM(f))

    WRITE(file_id, '(A," = (")') name//'%TAU_TOTAL'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = transmission%tau_total)
    WRITE(file_id, *) ')'

    IF (ANY(transmission%tausun_total_path1 /= 0._jprb)) THEN
      WRITE(file_id, '(A," = (")') name//'%TAUSUN_TOTAL_PATH1'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = transmission%tausun_total_path1)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ANY(transmission%tausun_total_path2 /= 0._jprb)) THEN
      WRITE(file_id, '(A," = (")') name//'%TAUSUN_TOTAL_PATH2'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = transmission%tausun_total_path2)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ANY(transmission%tau_total_cld /= 0._jprb)) THEN
      WRITE(file_id, '(A," = (")') name//'%TAU_TOTAL_CLD'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = transmission%tau_total_cld)
      WRITE(file_id, *) ')'
    ENDIF

    WRITE(file_id, '(A," = (")') name//'%TAU_LEVELS'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = transmission%tau_levels)
    WRITE(file_id, *) ')'

    IF (ANY(transmission%tausun_levels_path1 /= 0._jprb)) THEN
      WRITE(file_id, '(A," = (")') name//'%TAUSUN_LEVELS_PATH1'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = transmission%tausun_levels_path1)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ANY(transmission%tausun_levels_path2 /= 0._jprb)) THEN
      WRITE(file_id, '(A," = (")') name//'%TAUSUN_LEVELS_PATH2'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = transmission%tausun_levels_path2)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ANY(transmission%tau_levels_cld /= 0._jprb)) THEN
      WRITE(file_id, '(A," = (")') name//'%TAU_LEVELS_CLD'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = transmission%tau_levels_cld)
      WRITE(file_id, *) ')'
    ENDIF

    CLOSE(file_id, iostat = err)
    THROWM(err.NE.0,"Cannot close "//TRIM(f))
    CALL rttov_put_lun(file_id)

    CATCH
  END SUBROUTINE

  SUBROUTINE printreflectivity(err, reflectivity, name, f)
    USE rttov_chain, ONLY : print_array

    INTEGER(jpim),            INTENT(OUT) :: err
    TYPE(rttov_reflectivity), INTENT(IN)  :: reflectivity
    CHARACTER(LEN=*),         INTENT(IN)  :: name, f

    INTEGER(jpim) :: file_id

    TRY

    CALL rttov_get_lun(file_id)
    OPEN(file_id, file = f, form = 'formatted', iostat = err)
    THROWM(err.NE.0,"Cannot open "//TRIM(f))

    WRITE(file_id, '(A," = (")') name//'%ZEF'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = reflectivity%zef)
    WRITE(file_id, *) ')'

    WRITE(file_id, '(A," = (")') name//'%AZEF'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = reflectivity%azef)
    WRITE(file_id, *) ')'

    CLOSE(file_id, iostat = err)
    THROWM(err.NE.0,"Cannot close "//TRIM(f))
    CALL rttov_put_lun(file_id)

    CATCH
  END SUBROUTINE

  SUBROUTINE printemisterms(err, emis_terms, name, f)
    USE rttov_chain, ONLY : print_array

    INTEGER(jpim),                         INTENT(OUT) :: err
    TYPE(rttov_scatt_emis_retrieval_type), INTENT(IN)  :: emis_terms
    CHARACTER(LEN=*),                      INTENT(IN)  :: name, f

    INTEGER(jpim) :: file_id

    TRY

    CALL rttov_get_lun(file_id)
    OPEN(file_id, file = f, form = 'formatted', iostat = err)
    THROWM(err.NE.0,"Cannot open "//TRIM(f))

    WRITE(file_id, '(A," = (")') name//'%CFRAC'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = emis_terms%cfrac)
    WRITE(file_id, *) ')'

    WRITE(file_id, '(A," = (")') name//'%BSFC'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = emis_terms%bsfc)
    WRITE(file_id, *) ')'

    WRITE(file_id, '(A," = (")') name//'%TAU_CLD'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = emis_terms%tau_cld)
    WRITE(file_id, *) ')'

    WRITE(file_id, '(A," = (")') name//'%UP_CLD'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = emis_terms%up_cld)
    WRITE(file_id, *) ')'

    WRITE(file_id, '(A," = (")') name//'%DOWN_CLD'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = emis_terms%down_cld)
    WRITE(file_id, *) ')'

    WRITE(file_id, '(A," = (")') name//'%TAU_CLR'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = emis_terms%tau_clr)
    WRITE(file_id, *) ')'

    WRITE(file_id, '(A," = (")') name//'%UP_CLR'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = emis_terms%up_clr)
    WRITE(file_id, *) ')'

    WRITE(file_id, '(A," = (")') name//'%DOWN_CLR'
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = emis_terms%down_clr)
    WRITE(file_id, *) ')'


    CLOSE(file_id, iostat = err)
    THROWM(err.NE.0,"Cannot close "//TRIM(f))
    CALL rttov_put_lun(file_id)

    CATCH
  END SUBROUTINE

  SUBROUTINE printsurfvar(err, var, name, f)
    USE rttov_chain, ONLY : print_array

    INTEGER(jpim), INTENT(OUT) :: err
    REAL(jprb),    INTENT(IN)  :: var(:)
    CHARACTER(LEN=*),   INTENT(IN)  :: name, f

    INTEGER(jpim) :: file_id

    TRY

    CALL rttov_get_lun(file_id)
    OPEN(file_id, file = f, form = 'formatted', iostat = err)
    THROWM(err.NE.0,"Cannot open "//TRIM(f))

    WRITE(file_id, '(A," = (")') name
    CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = var)
    WRITE(file_id, *) ')'

    CLOSE(file_id, iostat = err)
    THROWM(err.NE.0,"Cannot close "//TRIM(f))
    CALL rttov_put_lun(file_id)

    CATCH
  END SUBROUTINE

  SUBROUTINE printpcscores(err, pccomp, name, f)
    USE rttov_chain, ONLY : print_array

    INTEGER(jpim),      INTENT(OUT) :: err
    TYPE(rttov_pccomp), INTENT(IN)  :: pccomp
    CHARACTER(LEN=*),   INTENT(IN)  :: name, f

    INTEGER(jpim) :: file_id

    TRY

    CALL rttov_get_lun(file_id)
    OPEN(file_id, file = f, form = 'formatted', iostat = err)
    THROWM(err.NE.0,"Cannot open "//TRIM(f))

    IF (ASSOCIATED(pccomp%total_pcscores)) THEN
      WRITE(file_id, '(A," = (")') name//'%TOTAL_PCSCORES'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = pccomp%total_pcscores)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ASSOCIATED(pccomp%total_pccomp)) THEN
      WRITE(file_id, '(A," = (")') name//'%TOTAL_PCCOMP'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = pccomp%total_pccomp)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ASSOCIATED(pccomp%bt_pccomp)) THEN
      WRITE(file_id, '(A," = (")') name//'%BT_PCCOMP'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = pccomp%bt_pccomp)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ASSOCIATED(pccomp%clear_pcscores)) THEN
      WRITE(file_id, '(A," = (")') name//'%CLEAR_PCSCORES'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = pccomp%clear_pcscores)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ASSOCIATED(pccomp%clear_pccomp)) THEN
      WRITE(file_id, '(A," = (")') name//'%CLEAR_PCCOMP'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = pccomp%clear_pccomp)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ASSOCIATED(pccomp%bt_clear_pccomp)) THEN
      WRITE(file_id, '(A," = (")') name//'%BT_CLEAR_PCCOMP'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = pccomp%bt_clear_pccomp)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ASSOCIATED(pccomp%cloudy_pcscores)) THEN
      WRITE(file_id, '(A," = (")') name//'%CLOUDY_PCSCORES'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = pccomp%cloudy_pcscores)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ASSOCIATED(pccomp%cloudy_pccomp)) THEN
      WRITE(file_id, '(A," = (")') name//'%CLOUDY_PCCOMP'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A1 = pccomp%cloudy_pccomp)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ASSOCIATED(pccomp%overcast_pcscores)) THEN
      WRITE(file_id, '(A," = (")') name//'%OVERCAST_PCSCORES'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = pccomp%overcast_pcscores)
      WRITE(file_id, *) ')'
    ENDIF

    IF (ASSOCIATED(pccomp%overcast_pccomp)) THEN
      WRITE(file_id, '(A," = (")') name//'%OVERCAST_PCCOMP'
      CALL print_array(INT(file_id,jpim), '(10000(' // TRIM(Gformat) // '))', A2 = pccomp%overcast_pccomp)
      WRITE(file_id, *) ')'
    ENDIF

    CLOSE(file_id, iostat = err)
    THROWM(err.NE.0,"Cannot close "//TRIM(f))
    CALL rttov_put_lun(file_id)

    CATCH
  END SUBROUTINE

  SUBROUTINE ppmvdry2ppmvwet(prof)
    TYPE(rttov_profile), INTENT(INOUT) :: prof

    ! Convert q first
    prof%q =  prof%q / (1._jprb + 1.E-06_jprb * prof%q)
    prof%s2m%q =  prof%s2m%q / (1._jprb + 1.E-06_jprb * prof%s2m%q)

    ! q is now ppmv wet
    IF (ASSOCIATED(prof%o3)) THEN
      prof%o3 = prof%o3 * (1._jprb - 1.E-06_jprb * prof%q)
      prof%s2m%o = prof%s2m%o * (1._jprb - 1.E-06_jprb * prof%s2m%q)
    ENDIF
    IF (ASSOCIATED(prof%co2)) THEN
      prof%co2 = prof%co2 * (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%co)) THEN
      prof%co = prof%co * (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%ch4)) THEN
      prof%ch4 = prof%ch4 * (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%so2)) THEN
      prof%so2 = prof%so2 * (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%n2o)) THEN
      prof%n2o = prof%n2o * (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
  END SUBROUTINE

  SUBROUTINE ppmvdry2kgkgwet(prof)
    TYPE(rttov_profile), INTENT(INOUT) :: prof

    IF (ASSOCIATED(prof%o3)) THEN
      prof%o3 = gas_mass(gas_id_ozone) * prof%o3 / (Mair * 1.E06_jprb + gas_mass(gas_id_watervapour) * prof%q)
      prof%s2m%o = gas_mass(gas_id_ozone) * prof%s2m%o / (Mair * 1.E06_jprb + gas_mass(gas_id_watervapour) * prof%s2m%q)
    ENDIF
    IF (ASSOCIATED(prof%co2)) THEN
      prof%co2 = gas_mass(gas_id_co2) * prof%co2 / (Mair * 1.E06_jprb + gas_mass(gas_id_watervapour) * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%co)) THEN
      prof%co = gas_mass(gas_id_co) * prof%co / (Mair * 1.E06_jprb + gas_mass(gas_id_watervapour) * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%ch4)) THEN
      prof%ch4 = gas_mass(gas_id_ch4) * prof%ch4 / (Mair * 1.E06_jprb + gas_mass(gas_id_watervapour) * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%so2)) THEN
      prof%so2 = gas_mass(gas_id_so2) * prof%so2 / (Mair * 1.E06_jprb + gas_mass(gas_id_watervapour) * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%n2o)) THEN
      prof%n2o = gas_mass(gas_id_n2o) * prof%n2o / (Mair * 1.E06_jprb + gas_mass(gas_id_watervapour) * prof%q)
    ENDIF

    ! Convert q last
    prof%q = gas_mass(gas_id_watervapour) * prof%q / (Mair * 1.E06_jprb  + gas_mass(gas_id_watervapour) * prof%q)
    prof%s2m%q = gas_mass(gas_id_watervapour) * prof%s2m%q / (Mair * 1.E06_jprb  + gas_mass(gas_id_watervapour) * prof%s2m%q)

  END SUBROUTINE

  SUBROUTINE ppmvwet2ppmvdry(prof)
    TYPE(rttov_profile), INTENT(INOUT) :: prof

    ! q is ppmv wet
    IF (ASSOCIATED(prof%o3)) THEN
      prof%o3 = prof%o3 / (1._jprb - 1.E-06_jprb * prof%q)
      prof%s2m%o = prof%s2m%o / (1._jprb - 1.E-06_jprb * prof%s2m%q)
    ENDIF
    IF (ASSOCIATED(prof%co2)) THEN
      prof%co2 = prof%co2 / (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%co)) THEN
      prof%co = prof%co / (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%ch4)) THEN
      prof%ch4 = prof%ch4 / (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%so2)) THEN
      prof%so2 = prof%so2 / (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF
    IF (ASSOCIATED(prof%n2o)) THEN
      prof%n2o = prof%n2o / (1._jprb - 1.E-06_jprb * prof%q)
    ENDIF

    ! Convert q last
    prof%q =  prof%q / (1._jprb - 1.E-06_jprb * prof%q)
    prof%s2m%q =  prof%s2m%q / (1._jprb - 1.E-06_jprb * prof%s2m%q)

  END SUBROUTINE

  SUBROUTINE ppmvwet2kgkgwet(prof)
    TYPE(rttov_profile), INTENT(INOUT) :: prof

    CALL ppmvwet2ppmvdry(prof)
    CALL ppmvdry2kgkgwet(prof)

  END SUBROUTINE

  SUBROUTINE kgkgwet2ppmvdry(prof)
    TYPE(rttov_profile), INTENT(INOUT) :: prof

    IF (ASSOCIATED(prof%o3)) THEN
      prof%o3 = (1.E06_jprb * Mair / gas_mass(gas_id_ozone)) * prof%o3 / (1._jprb - prof%q)
      prof%s2m%o = (1.E06_jprb * Mair / gas_mass(gas_id_ozone)) * prof%s2m%o / (1._jprb - prof%s2m%q)
    ENDIF
    IF (ASSOCIATED(prof%co2)) THEN
      prof%co2 = (1.E06_jprb * Mair / gas_mass(gas_id_co2)) * prof%co2 / (1._jprb - prof%q)
    ENDIF
    IF (ASSOCIATED(prof%co)) THEN
      prof%co = (1.E06_jprb * Mair / gas_mass(gas_id_co)) * prof%co / (1._jprb - prof%q)
    ENDIF
    IF (ASSOCIATED(prof%ch4)) THEN
      prof%ch4 = (1.E06_jprb * Mair / gas_mass(gas_id_ch4)) * prof%ch4 / (1._jprb - prof%q)
    ENDIF
    IF (ASSOCIATED(prof%so2)) THEN
      prof%so2 = (1.E06_jprb * Mair / gas_mass(gas_id_so2)) * prof%so2 / (1._jprb - prof%q)
    ENDIF
    IF (ASSOCIATED(prof%n2o)) THEN
      prof%n2o = (1.E06_jprb * Mair / gas_mass(gas_id_n2o)) * prof%n2o / (1._jprb - prof%q)
    ENDIF

    ! Convert q last
    prof%q = (1.E06_jprb * Mair / gas_mass(gas_id_watervapour)) * prof%q / (1._jprb - prof%q)
    prof%s2m%q = (1.E06_jprb * Mair / gas_mass(gas_id_watervapour)) * prof%s2m%q / (1._jprb - prof%s2m%q)

  END SUBROUTINE

  SUBROUTINE kgkgwet2ppmvwet(prof)
    TYPE(rttov_profile), INTENT(INOUT) :: prof

    CALL kgkgwet2ppmvdry(prof)
    CALL ppmvdry2ppmvwet(prof)

  END SUBROUTINE

END MODULE rttov_test_mod
