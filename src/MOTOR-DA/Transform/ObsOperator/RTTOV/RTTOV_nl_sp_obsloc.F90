!!--------------------------------------------------------------------------------------------------
! PROJECT           : rttov_nl_sp
! AFFILIATION       : GBA-MWF(SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu, 2021/8/27
! Information       : Add interfaces to call rttov NL..
!                   : Only the clear-sky capability is enabled so far.
!                   : Requirement: rttov V13
! Function          : Calculate J_radiance and grad_J_radiance wrt CV variables
!
! The usual steps to take when running RTTOV NL are as follows:
!   1. Specify required RTTOV options
!   2. Read coefficients
!   3. Allocate RTTOV input and output structures
!   4. Set up the chanprof array with the channels/profiles to simulate
!   5. Read input profile(s)
!   6. Set up surface emissivity and/or reflectance
!   7. Call rttov_direct/rttov_k and store results
!   8. Deallocate all structures and arrays
! For rttov nl model,
! y=H(x), where H is the RTTOV forward (or direct) model.
! Input: model profiles + surface variables + satellite info (angles, channels, etc)
! Output:
! The TOA radiance outputs are stored in the rttov_radiance structure.
! For channels with significant thermal component (lamda>3um), RTTOV computes BTs;
! For solar-affected channels (lamdda<5um), RTTOV computes reflectances;
! For all channels, RTTOV computes radiances. Each array has size nchanprof:
! >radiance%clear(:), unit: mW/m^2/sr/cm^-1
! >radiance%total(:), unit: mW/m^2/sr/cm^-1
! >radiance%bt_clear(:), unit: K
! >radiance%bt(:), unit: K
! >radiance%refl_clear(:), unit: none
! >radaicne%refl(:), unit: none
! RTTOV K (Jacobian) model is useful for computing the sensitivity of the TOA radainces to
! each input variable.
! Note:
! (1) all visible/IR scattering-related inputs are provided on layers rather than levels
!     where layer i is bounded by input pressure levels i and i+1.

!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE rttov_nl_sp_obsloc_m
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: r_kind, i_kind
  ! Self defined rttov utilities and types
  USE rttov_typedefine_m
  USE rttov_utils_m, ONLY: rttov_utils_t
  USE rttov_TYPES, ONLY: rttov_options, rttov_profile, rttov_coefs, &
                         rttov_radiance, rttov_transmission, rttov_emissivity, &
                         rttov_chanprof
  USE parkind1, ONLY: jpim
  USE rttov_unix_env, ONLY: rttov_exit

  TYPE rttov_nl_sp_obsloc_t
    TYPE(State_t), POINTER :: X
    TYPE(rttov_utils_t) :: rttov_utils  ! Pointer or not?
    CHARACTER(LEN=1024) :: configFile

  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: initialize
    FINAL :: destructor
    PROCEDURE, PUBLIC :: simobsloc
  END TYPE rttov_nl_sp_obsloc_t

CONTAINS

  SUBROUTINE initialize(this, configFile, X, inst_name, platform_name)
    IMPLICIT NONE
    CLASS(rttov_nl_sp_obsloc_t) :: this
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(State_t), TARGET, INTENT(IN) :: X
    CHARACTER(len=*), INTENT(IN)   :: inst_name, platform_name

    CALL this%rttov_utils%initialize(configFile, X)
    this%X => X

    !   1. Specify required RTTOV options
    CALL this%rttov_utils%rttov_opts_all()
    CALL this%rttov_utils%rttov_opts_ir()
    CALL this%rttov_utils%rttov_opts_mw()
    CALL this%rttov_utils%rttov_opts_vis_nir()
    !CALL this%rttov_utils%rttov_GetMVars(X)
    CALL this%rttov_utils%rttov_GetOption(TRIM(inst_name), TRIM(platform_name))
    !   2. Read coefficients
    CALL this%rttov_utils%read_coefs()

  END SUBROUTINE initialize

  SUBROUTINE simobsloc(this, vLevel, ichan, pres, t, q, u, v, tskin, landmask, elevation, soil_type, snowc, &
                                lat, lon, SatAngles, BT, BTclr, diag, qc, qr, qi, qs, qg)
    USE rttov_typedefine_m, ONLY : rttov_diag_output
    CLASS(rttov_nl_sp_obsloc_t) :: this
    INTEGER(i_kind), INTENT(IN)  :: vLevel
    REAL(r_kind), INTENT(IN) :: pres(:), t(:), q(:), u(:), v(:)
    REAL(r_kind), INTENT(IN) :: tskin, landmask, elevation, soil_type, snowc, lat, lon
    INTEGER(i_kind), INTENT(IN) :: ichan
    REAL(r_kind), INTENT(IN) :: SatAngles(4)
    REAL(r_kind) :: InLatLon(2)
    REAL(r_kind), INTENT(INOUT) :: BT, BTclr
    TYPE(rttov_diag_output), INTENT(INOUT), OPTIONAL :: diag
    REAL(r_kind), INTENT(IN), OPTIONAL :: qc(:), qr(:), qi(:), qs(:), qg(:)

    ! Local variables:
    TYPE(ObsSet_t)  :: Y
    INTEGER(i_kind) :: i_prof, it, nchans, nlevels, nch
    INTEGER(i_kind) :: nprofs, nsimobs
    INTEGER(i_kind) :: asw
    INTEGER(kind=jpim)  :: errstatus
    INTEGER(i_kind) :: jchan, l, i, j, k
    REAL(r_kind) :: t1, t2, t3, t4

    INCLUDE 'rttov_user_options_checkinput.interface'
    INCLUDE 'rttov_print_profile.interface'
    INCLUDE 'rttov_k.interface'
    INCLUDE 'rttov_direct.interface'
    INCLUDE 'rttov_init_emis_refl.interface'
    INCLUDE 'rttov_dealloc_coefs.interface'
    INCLUDE 'rttov_alloc_transmission.interface'
    INCLUDE 'rttov_alloc_rad.interface'

    IF (ANY(t<0.1)) THEN  ! exclude invalid profiles
      BT = missing
      BTclr = missing
    ELSE
      InLatLon(1) = lat
      InLatLon(2) = lon

      nlevels = vLevel
      IF (this%rttov_utils%opts%rt_ir%addclouds .AND. PRESENT(qc)) THEN
        CALL this%rttov_utils%rttov_GetMVars_obsloc(pres, t, q, u, v, tskin, landmask, elevation, soil_type, snowc, &
                                                    qc, qr, qi, qs, qg)
      ELSE
        CALL this%rttov_utils%rttov_GetMVars_obsloc(pres, t, q, u, v, tskin, landmask, elevation, soil_type, snowc)
      END IF
      nchans = 1
      nprofs = 1
      nsimobs = nprofs * nchans

      !   3. Allocate RTTOV input and output structures
      !   4. Set up the chanprof array with the channels/profiles to simulate

      asw = 1 !allocate

      ! CALL this%rttov_utils%alloc_profs(asw, nprofs, nchans, nlevels)
      IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%use_kmodel .or. &
          this%rttov_utils%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%write_jac) THEN
        ! CALL this%rttov_utils%alloc_profs_k(asw, nprofs, nchans, nlevels)
        CALL this%rttov_utils%alloc_k(asw, nprofs, nchans, nlevels, ichan)
      ELSE
        CALL this%rttov_utils%alloc_direct(asw, nprofs, nchans, nlevels, ichan)
      END IF

      !   5. Read input profile(s)
      CALL this%rttov_utils%rttov_MVar2RttovVar(nlevels, SatAngles, InLatLon)

      IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%debug_dev) THEN
          
          i_prof = 1
          
          CALL rttov_print_opts(this%rttov_utils%opts)
          PRINT *, 'rttov_print_profile'
          CALL rttov_print_profile(this%rttov_utils%rttov_in_calc%profs(i_prof))

          ! PRINT *, 'rttov_user_profile_checkinput'
          ! call rttov_user_profile_checkinput(errstatus, &
          !                                 this%rttov_utils%rttov_setup_input%rttov_in_userdefs%opts, &
          !                                 this%rttov_utils%coefs, &
          !                                 this%rttov_utils%rttov_in_calc%profs(i_prof))

          PRINT *, 'rttov_user_options_checkinput'
          CALL rttov_user_options_checkinput(errstatus, this%rttov_utils%opts, &
                                            this%rttov_utils%coefs)
      END IF

      !   6. Set up surface emissivity and/or reflectance
      ! Method 1: specify them in RTTOV_utils.F90
      ! CALL this%rttov_utils%init_emiss(nsimobs)
      ! Method 2: initialise all inputs to zero, and calculate them within RTTOV
      ! Calculate emissivity within RTTOV where the input emissivity value is
      ! zero or less
      ! this%rttov_utils%rttov_in_calc%emis(:)%emis_in = 0.0
      ! Method 3: Method 1 + Method 2, that is, if surftype_sea, calculate within RTTOV, otherwise, specify in RTTOV_utils.
      this%rttov_utils%rttov_in_calc%emis(:)%emis_in = 0.0
      CALL this%rttov_utils%init_emiss(nsimobs)

      ! CALL rttov_init_emis_refl(this%rttov_utils%rttov_in_calc%emis, &
      !                           this%rttov_utils%rttov_in_calc%reflectance)

      this%rttov_utils%rttov_in_calc%calcemis(:) = &
        (this%rttov_utils%rttov_in_calc%emis(:)%emis_in <= 0.0)
      this%rttov_utils%rttov_in_calc%calcrefl(:) = &
        (this%rttov_utils%rttov_in_calc%reflectance(:)%refl_in <= 0.0)

      !   7. Call rttov_direct and store results
      asw = 1 ! allocate

      IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%use_kmodel) THEN

        ! The input/output K variables must be initialised to zero
        ! before every call to rttov_k:

        ! Initialise RTTOV Jacobian structures to zero
        CALL rttov_k(errstatus, &! out   error flag
                      this%rttov_utils%rttov_in_calc%chanprof, &! in LOCAL channel and profile index structure
                      this%rttov_utils%opts, &! in    options structure
                      this%rttov_utils%rttov_in_calc%profs, &! in    profile array
                      this%rttov_utils%rttov_in_calc%profs_k, &! in    profile array
                      this%rttov_utils%coefs, &! in    coefficients structure
                      this%rttov_utils%rttov_in_calc%transm, &! inout computed transmittances
                      this%rttov_utils%rttov_in_calc%transm_k, &! inout computed transmittances
                      this%rttov_utils%rttov_in_calc%radiance, &! inout computed radiances
                      this%rttov_utils%rttov_in_calc%radiance_k, &! inout computed radiances
                      calcemis=this%rttov_utils%rttov_in_calc%calcemis, &! in    flag for internal emissivity calcs
                      emissivity=this%rttov_utils%rttov_in_calc%emis, &!, &! inout input/output emissivities per channel
                      emissivity_k=this%rttov_utils%rttov_in_calc%emis_k)                    ! inout input/output emissivities per channel
        !cld_opt_param = .False.)

      ELSE

      ! PRINT *, 'chanprof = ', this%rttov_utils%rttov_in_calc%chanprof(1)%chan, this%rttov_utils%coefs%coef%fmv_chn
      CALL rttov_direct(errstatus, &! out   error flag
                        this%rttov_utils%rttov_in_calc%chanprof, &! in LOCAL channel and profile index structure
                        this%rttov_utils%opts, &! in    options structure
                        this%rttov_utils%rttov_in_calc%profs, &! in    profile array
                        this%rttov_utils%coefs, &! in    coefficients structure
                        this%rttov_utils%rttov_in_calc%transm, &! inout computed transmittances
                        this%rttov_utils%rttov_in_calc%radiance, &! inout computed radiances
                        calcemis=this%rttov_utils%rttov_in_calc%calcemis, &! in    flag for internal emissivity calcs
                        emissivity=this%rttov_utils%rttov_in_calc%emis)!, &
      END IF

      ! TODO: replace the following by storing the outputs
      BT = this%rttov_utils%rttov_in_calc%radiance%bt(1)
      BTclr = this%rttov_utils%rttov_in_calc%radiance%bt_clear(1)
      ! PRINT *, 'BT = ', BT
      IF (PRESENT(diag)) &
      diag%bt = this%rttov_utils%rttov_in_calc%radiance%bt(1)
      IF (PRESENT(diag)) &
      diag%btclr = this%rttov_utils%rttov_in_calc%radiance%bt_clear(1)
      
      IF ( .NOT. this%rttov_utils%opts%rt_ir%addclouds .AND. &
      PRESENT(diag) .AND. PRESENT(diag) .AND. bt /= btclr) THEN
        PRINT*, 'FATAL ERROR: bt should be identical to btclr in clear skies'
      END IF 

      !   8. Deallocate all structures and arrays
      asw = 0 !deallocate

      ! CALL this%rttov_utils%alloc_profs(asw, nprofs, nchans, nlevels)
      IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%use_kmodel .or. &
          this%rttov_utils%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%write_jac) THEN

        ! CALL this%rttov_utils%zero_k()
        CALL this%rttov_utils%alloc_k(asw, nprofs, nchans, nlevels, ichan)

        ! CALL this%rttov_utils%alloc_profs_k(asw, nprofs, nchans, nlevels)
      ELSE
        CALL this%rttov_utils%alloc_direct(asw, nprofs, nchans, nlevels, ichan)
      END IF
      CALL this%rttov_utils%rttov_dealloc_array()
      ! PRINT *, 'rttov_nl_sp_simobs OVER'
    END IF
  END SUBROUTINE simobsloc

  IMPURE ELEMENTAL SUBROUTINE destructor(this)
    IMPLICIT NONE
    TYPE(rttov_nl_sp_obsloc_t), INTENT(INOUT) :: this

    CALL this%rttov_utils%rttov_dealloc_types()

  END SUBROUTINE destructor

END MODULE rttov_nl_sp_obsloc_m
