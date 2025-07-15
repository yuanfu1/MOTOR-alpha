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
MODULE rttov_nl_sp_m
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

  TYPE rttov_nl_sp_t
    TYPE(State_t), POINTER :: X
    TYPE(rttov_utils_t) :: rttov_utils
    CHARACTER(LEN=1024) :: configFile

  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: initialize
    FINAL :: destructor_nl_sp
    PROCEDURE, PUBLIC :: rttov_nl_sp_simobs
  END TYPE rttov_nl_sp_t

CONTAINS

  SUBROUTINE initialize(this, configFile, X, inst_name, platform_name)
    IMPLICIT NONE
    CLASS(rttov_nl_sp_t) :: this
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
    ! Note that coef%id_sensor should be determined according to inst_name and platform_name
    ! No need to explicitly set id_sensor = sensor_id_ir or sensor_id_mw or sensor_id_hi
    CALL this%rttov_utils%rttov_GetOption(TRIM(inst_name), TRIM(platform_name))
    !   2. Read coefficients
    CALL this%rttov_utils%read_coefs()
    PRINT *, 'check init of inst_name and platform_name: ', TRIM(inst_name), '   ', TRIM(platform_name)

  END SUBROUTINE initialize

  SUBROUTINE rttov_nl_sp_simobs(this, X, ichan, locX, SatAngles, Hx_Sp, diag)
    USE rttov_typedefine_m, ONLY: rttov_diag_output
    CLASS(rttov_nl_sp_t) :: this
    TYPE(State_t), INTENT(IN) ::  X
    INTEGER(i_kind), INTENT(IN) :: ichan
    INTEGER(i_kind), INTENT(IN) :: locX(2)
    REAL(r_kind), INTENT(IN) :: SatAngles(4)
    REAL(r_kind) :: InLatLon(2)
    REAL(r_kind), INTENT(INOUT) :: Hx_Sp
    TYPE(rttov_diag_output), INTENT(INOUT), OPTIONAL :: diag

    ! Local variables:
    TYPE(ObsSet_t)  :: Y
    INTEGER(i_kind) :: i_prof, it, nchans, nlevels, vLevel, nch
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
    INCLUDE 'rttov_print_opts.interface'
    INCLUDE 'rttov_user_profile_checkinput.interface'

    vLevel = X%sg%vLevel
    InLatLon(1) = X%sg%cell_cntr(1, locX(1))
    InLatLon(2) = X%sg%cell_cntr(2, locX(1))

    nlevels = vLevel

    CALL this%rttov_utils%rttov_GetMVars(X, locX)
    nchans = 1
    nprofs = 1
    nsimobs = nprofs * nchans

    !   3. Allocate RTTOV input and output structures
    !   4. Set up the chanprof array with the channels/profiles to simulate

    asw = 1 !allocate

    ! CALL this%rttov_utils%alloc_profs(asw, nprofs, nchans, nlevels)
    IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%use_kmodel .OR. &
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
        
        ! PRINT *, 'rttov_print_opts'
        ! CALL rttov_print_opts(this%rttov_utils%opts)
        ! PRINT *, 'rttov_print_profile'
        CALL rttov_print_profile(this%rttov_utils%rttov_in_calc%profs(i_prof))

        ! PRINT *, 'rttov_user_profile_checkinput'
        ! call rttov_user_profile_checkinput(errstatus, &
        !                                 this%rttov_utils%opts, &
        !                                 this%rttov_utils%coefs, &
        !                                 this%rttov_utils%rttov_in_calc%profs(i_prof))

        ! PRINT *, 'rttov_user_options_checkinput'
        ! CALL rttov_user_options_checkinput(errstatus, this%rttov_utils%opts, &
        !                                   this%rttov_utils%coefs)
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
                   radiance2=this%rttov_utils%rttov_in_calc%radiance2, &
                   calcemis=this%rttov_utils%rttov_in_calc%calcemis, &! in    flag for internal emissivity calcs
                   emissivity=this%rttov_utils%rttov_in_calc%emis, &!, &! inout input/output emissivities per channel
                   emissivity_k=this%rttov_utils%rttov_in_calc%emis_k)                    ! inout input/output emissivities per channel
      !cld_opt_param = .False.)

    ELSE

      CALL rttov_direct(errstatus, &! out   error flag
                        this%rttov_utils%rttov_in_calc%chanprof, &! in LOCAL channel and profile index structure
                        this%rttov_utils%opts, &! in    options structure
                        this%rttov_utils%rttov_in_calc%profs, &! in    profile array
                        this%rttov_utils%coefs, &! in    coefficients structure
                        this%rttov_utils%rttov_in_calc%transm, &! inout computed transmittances
                        this%rttov_utils%rttov_in_calc%radiance, &! inout computed radiances
                        radiance2=this%rttov_utils%rttov_in_calc%radiance2, & ! out for up/down welling radiances
                        calcemis=this%rttov_utils%rttov_in_calc%calcemis, &! in    flag for internal emissivity calcs
                        emissivity=this%rttov_utils%rttov_in_calc%emis)!, &
    END IF

    ! TODO: replace the following by storing the outputs
    ! BLOCK
    !   INTEGER :: k 
    !   PRINT *, 'profiles_k%cloud(1:6), profiles%cloud(1), profiles%cloud(6) '
    !   DO k = 1, nlevels-1
    !     PRINT *, K, this%rttov_utils%rttov_in_calc%profs_k(1)%cloud(1:6,k), this%rttov_utils%rttov_in_calc%profs(1)%cloud(1,k), this%rttov_utils%rttov_in_calc%profs(1)%cloud(6,k)
    !     ! PRINT *, 'rttov_k ice : ', this%rttov_utils%rttov_in_calc%profs_k(1)%cloud(6,1:nlevels-1:1)
    !     ! PRINT *, 'rttov ice:', this%rttov_utils%rttov_in_calc%profs(1)%cloud(6,1:nlevels-1:1)
    !     ! PRINT *, 'rttov_k Q : ', this%rttov_utils%rttov_in_calc%profs_k(1)%q(1:nlevels:1)
    !     ! PRINT *, 'rttov Q:', this%rttov_utils%rttov_in_calc%profs(1)%q(1:nlevels:1)
    !   END DO 
    ! END BLOCK

    Hx_Sp = this%rttov_utils%rttov_in_calc%radiance%bt(1)
    ! PRINT *, 'Hx_Sp = ', Hx_Sp
    IF (PRESENT(diag)) &
      diag%bt = this%rttov_utils%rttov_in_calc%radiance%bt(1)
    IF (PRESENT(diag)) &
      diag%btclr = this%rttov_utils%rttov_in_calc%radiance%bt_clear(1)

    IF (.NOT. this%rttov_utils%opts%rt_ir%addclouds .AND. &
        PRESENT(diag) .AND. PRESENT(diag) .AND. bt /= btclr) THEN
      PRINT *, 'FATAL ERROR: bt should be identical to btclr in clear skies'
    END IF

    !!!---------- Output RTTOV FWD results -----------!!!
    ! Write for diagnosis
    IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%use_kmodel &
        .AND. present(diag) ) &
      CALL this%rttov_utils%Convert_to_diag_out(diag)

    IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%debug_dev) THEN
      IF (PRESENT(diag)) PRINT *, 'Brightness Temperature (K) = ', diag%bt

      PRINT *, 'Radiance overcast (mW/m2/sr/cm-1) = '
      PRINT *, this%rttov_utils%rttov_in_calc%radiance%cloudy

      ! PRINT *, 'Transmittance 0-1 = '
      ! PRINT *, 'tau_total: ', this%rttov_utils%rttov_in_calc%transm%tau_total
      ! PRINT *, 'tau_levels: ',this%rttov_utils%rttov_in_calc%transm%tau_levels
      ! PRINT *, 'tausun_total_path2: ', this%rttov_utils%rttov_in_calc%transm%tausun_total_path2
      ! PRINT *, 'tausun_levels_path2: ', this%rttov_utils%rttov_in_calc%transm%tausun_levels_path2
      ! PRINT *, 'tausun_total_path1: ', this%rttov_utils%rttov_in_calc%transm%tausun_total_path1
      ! PRINT *, 'tausun_levels_path1: ', this%rttov_utils%rttov_in_calc%transm%tausun_levels_path1
      ! PRINT *, 'tau_total_cld: ', this%rttov_utils%rttov_in_calc%transm%tau_total_cld
      ! PRINT *, 'tau_levels_cld: ', this%rttov_utils%rttov_in_calc%transm%tau_levels_cld

      ! PRINT *, 'Radiance2 (mW/m2/sr/cm-1) = '
      ! PRINT *, 'upclear: ', this%rttov_utils%rttov_in_calc%radiance2%upclear
      ! PRINT *, 'dnclear: ', this%rttov_utils%rttov_in_calc%radiance2%dnclear
      ! PRINT *, 'refldnclear: ', this%rttov_utils%rttov_in_calc%radiance2%refldnclear
      ! PRINT *, 'up: ', this%rttov_utils%rttov_in_calc%radiance2%up
      ! PRINT *, 'down: ', this%rttov_utils%rttov_in_calc%radiance2%down
      ! PRINT *, 'surf: ', this%rttov_utils%rttov_in_calc%radiance2%surf

      IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%use_kmodel) THEN
        PRINT *, '-------------Print profile_k-----------'

        !CALL rttov_print_profile(this%rttov_utils%rttov_in_calc%profs_k(i_prof))
        PRINT *, '================ Print from NL ==================='
        PRINT *, "Level", "Pressure", "JAC(Temp) (K/K)", "JAC(WV) ", this%rttov_utils%rttov_in_calc%profs(1)%gas_units
        jch = 1
        DO l = 1, this%rttov_utils%rttov_in_calc%profs_k(jch)%nlevels, 10
          PRINT *, 'Level = ', l
          PRINT *, 'P: ', this%rttov_utils%rttov_in_calc%profs(1)%p(l)
          PRINT *, 't, k_t: ', this%rttov_utils%rttov_in_calc%profs(jch)%t(l), this%rttov_utils%rttov_in_calc%profs_k(jch)%t(l)
          PRINT *, 'q, q_t: ', this%rttov_utils%rttov_in_calc%profs(jch)%q(l), this%rttov_utils%rttov_in_calc%profs_k(jch)%q(l)
          ! PRINT *, 'transm_k: ', this%rttov_utils%rttov_in_calc%transm_k%tau_total, this%rttov_utils%rttov_in_calc%transm_k%tau_levels
          ! PRINT *, 'radiance_k: ', this%rttov_utils%rttov_in_calc%radiance_k%bt
          ! PRINT *, 'emis_k: ', this%rttov_utils%rttov_in_calc%emis_k(1)%emis_out, this%rttov_utils%rttov_in_calc%emis_k(1)%emis_in
        END DO
      END IF
    END IF

    !   8. Deallocate all structures and arrays
    asw = 0 !deallocate

    ! CALL this%rttov_utils%alloc_profs(asw, nprofs, nchans, nlevels)
    IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%use_kmodel .OR. &
        this%rttov_utils%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%write_jac) THEN

      ! CALL this%rttov_utils%zero_k()
      CALL this%rttov_utils%alloc_k(asw, nprofs, nchans, nlevels, ichan)

      ! CALL this%rttov_utils%alloc_profs_k(asw, nprofs, nchans, nlevels)
      ! we will re-calculate profs_k in RTTOV_tlad.F90
    ELSE
      CALL this%rttov_utils%alloc_direct(asw, nprofs, nchans, nlevels, ichan)
    END IF
    CALL this%rttov_utils%rttov_dealloc_array()
    ! PRINT *, 'rttov_nl_sp_simobs OVER'
  END SUBROUTINE rttov_nl_sp_simobs

  IMPURE ELEMENTAL SUBROUTINE destructor_nl_sp(this)
    IMPLICIT NONE
    TYPE(rttov_nl_sp_t), INTENT(INOUT) :: this

    CALL this%rttov_utils%rttov_dealloc_types()
  END SUBROUTINE destructor_nl_sp

END MODULE rttov_nl_sp_m
