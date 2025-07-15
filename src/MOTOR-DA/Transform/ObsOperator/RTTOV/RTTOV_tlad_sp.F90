!!--------------------------------------------------------------------------------------------------
! PROJECT           : rttov_tl_spad
! AFFILIATION       : GBA-MWF(SIMI)
! AUTOHR(S)         : Yali Wu
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yali Wu, 2021/8/27
! Information       : Add interfaces to call rttov TL/AD..
!                   : Only the clear-sky capability is enabled so far.
!                   : Requirement: rttov V13
! Function          : Calculate J_radiance and grad_J_radiance wrt CV variables
!
! Use K model for the calculation of TL/AD, which is faster than directly call TL/AD models.
! For the rttov TL model,
! δy = H1(x0)δx
! δy: the change in radiance
! H1(x0): a linear relationship about a given atmospheric state x0
! δx: a change of the state vector
! The Jacobian is a Nxm matrix
! The elements of H1 contain the partial derivatives \frac{\partial y_i}{\partial x_j}
! where the subscript i refers to channel number and j to position in state vector.
! rttov_k_sp computes the H(x0) matrix for each input profile.
! We always need the tangent linear value δy, but not always need the full Jacobian matrix H.
! TL input:
!         Profile(x), channel specifications, observation geometry,
!         increments in profile variables (δx)
! TL output: δy, i.e., increments in radiances (δBT)
! The TL is a mx1 vector
!
! For the rttov AD model,
! Conversely the adjoint routines compute the change in any scalar quantity to the state vector δx
! for an assumed atmospheric state, x0, given a change in the radiances, δy.
! δx = H2(x0)δy
! δy: the change in radiance
! H2(x0): adjoint, \frac{\partial x_i}{\partial y_j}
! δx: a change of the state vector
! AD input:
!         Profile(x), channel specifications, observation geometry,
!         normalised departures: R^-1.d=R^-1.(y-H(x))=R^-1.(y-H(xb)-H1δx)
! AD output: gradients of cost wrt profile variables
! AD transforms gradients from observation space to model space
! The AD is a Nx1 vector
!> @brief
!!
MODULE rttov_tlad_sp_m
  USE State_m, ONLY: State_t
  USE ObsSet_m, ONLY: ObsSet_t
  USE ObsField_m, ONLY: ObsField_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: r_kind, i_kind
  ! Self defined rttov utilities and types
  USE rttov_typedefine_m
  USE rttov_utils_m, ONLY: rttov_utils_t
  USE rttov_TYPES, ONLY: rttov_options, rttov_profile, rttov_coefs, &
                         rttov_radiance, rttov_transmission, rttov_emissivity, &
                         rttov_chanprof
  USE rttov_CONST, ONLY: errorstatus_success, &
                         gas_id_watervapour, &
                         gas_unit_ppmvdry, gas_unit_specconc, &
                         gas_unit_ppmv, &
                         platform_name, inst_name
  USE parameters_m, ONLY: q_mr_to_ppmv, ch4_mr_to_ppmv, &
                          o3_mr_to_ppmv, co2_mr_to_ppmv, &
                          co_mr_to_ppmv, h2o_mr_to_ppmv
  USE parkind1, ONLY: jpim
  USE rttov_unix_env, ONLY: rttov_exit

  TYPE rttov_tlad_sp_t
    TYPE(State_t), POINTER :: X
    TYPE(rttov_utils_t) :: rttov_utils
    CHARACTER(LEN=1024) :: configFile
    INTEGER(i_kind)     :: locX(2)

  CONTAINS
    PROCEDURE, PUBLIC, PASS(this) :: initialize
    FINAL :: destructor_tlad_sp
    PROCEDURE, PUBLIC :: rttov_k_sp_simobs
    PROCEDURE, PUBLIC :: rttov_tl_sp
    PROCEDURE, PUBLIC :: rttov_ad_sp
    PROCEDURE, PUBLIC :: rttov_tl_sp_simobs
    PROCEDURE, PUBLIC :: rttov_ad_sp_simobs
  END TYPE rttov_tlad_sp_t

CONTAINS

  SUBROUTINE initialize(this, configFile, X, inst_name, platform_name)
    IMPLICIT NONE
    CLASS(rttov_tlad_sp_t) :: this
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

    CALL this%rttov_utils%read_coefs()
    CALL X%sg%mpddInfo_sg%barrier()
    PRINT *, "read_coefs of myrank = ", X%sg%mpddInfo_sg%myrank, ' is over'
    CALL X%sg%mpddInfo_sg%barrier()

  END SUBROUTINE initialize

  SUBROUTINE rttov_k_sp_simobs(this, X, ichan, locX, SatAngles)
    CLASS(rttov_tlad_sp_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    INTEGER(i_kind), INTENT(IN) :: ichan
    INTEGER(i_kind), INTENT(IN) :: locX(2)
    REAL(r_kind), INTENT(IN) :: SatAngles(4)
    REAL(r_kind) :: InLatLon(2)
    ! Local variables:
    INTEGER(i_kind) :: i, i_prof, it, nlevels, vLevel
    INTEGER(i_kind) :: nprofs, ninst, nsimobs
    INTEGER(i_kind) :: asw
    INTEGER(kind=jpim)  :: errstatus

    INCLUDE 'rttov_k.interface'
    INCLUDE 'rttov_init_emis_refl.interface'
    INCLUDE 'rttov_dealloc_coefs.interface'
    INCLUDE 'rttov_alloc_transmission.interface'
    INCLUDE 'rttov_alloc_rad.interface'

    vLevel = this%X%sg%vLevel
    InLatLon(1) = X%sg%cell_cntr(1, locX(1))
    InLatLon(2) = X%sg%cell_cntr(2, locX(1))

    CALL this%rttov_utils%rttov_GetMVars(X, locX)

    nlevels = vLevel
    nprofs = 1 ! sg%num_icell or sg%num_cell
    nchans = 1
    nsimobs = nprofs * nchans

    asw = 1 !allocate

    !CALL this%rttov_utils%alloc_profs(asw, nprofs, nchans, nlevels)
    !CALL this%rttov_utils%alloc_profs_k(asw, nprofs, nchans, nlevels)
    CALL this%rttov_utils%alloc_k(asw, nprofs, nchans, nlevels, ichan)

    CALL this%rttov_utils%rttov_MVar2RttovVar(nlevels, SatAngles, InLatLon)

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

    asw = 1 ! allocate
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

    IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%debug_dev) THEN
      PRINT *, 'Brightness Temperature (K) = '
      PRINT *, this%rttov_utils%rttov_in_calc%radiance%bt

      PRINT *, 'Radiance (mW/m2/sr/cm-1) = '
      PRINT *, this%rttov_utils%rttov_in_calc%radiance%total

      PRINT *, 'Radiance overcast (mW/m2/sr/cm-1) = '
      PRINT *, this%rttov_utils%rttov_in_calc%radiance%cloudy

      PRINT *, '-------------Print profile_k-----------'
      jchan = 1 ! print the 1st channel
      PRINT *, '================ Print from K ==================='
      PRINT *, "Level", "Pressure", "JAC(Temp) (K/K)", "JAC(WV) ", this%rttov_utils%rttov_in_calc%profs(1)%gas_units
      DO l = 1, this%rttov_utils%rttov_in_calc%profs_k(jchan)%nlevels
        PRINT *, 'Level = ', l
        PRINT *, this%rttov_utils%rttov_in_calc%profs(1)%p(l)
        PRINT *, this%rttov_utils%rttov_in_calc%profs_k(jchan)%t(l)
        PRINT *, this%rttov_utils%rttov_in_calc%profs_k(jchan)%q(l)
      END DO

      ! Deallocate structures for rttov_direct/rttov_k_sp
      ! 'DO NOT DEALLOCATION IN THIS SUBROUTINE'
      ! 'WE NEED THE CALCULATED TRAJECTORIES AND PARAMETERS'
    END IF

    !PRINT *, 'rttov_k_sp_simobs is successfully set'

  END SUBROUTINE rttov_k_sp_simobs

  SUBROUTINE rttov_tl_sp(this, X, dX, ichan, locX, SatAngles)
    CLASS(rttov_tlad_sp_t) :: this
    TYPE(State_t), INTENT(IN) :: X, dX
    INTEGER(i_kind), INTENT(IN) :: ichan
    INTEGER(i_kind), INTENT(IN) :: locX(2)
    REAL(r_kind), INTENT(IN) :: SatAngles(4)
    REAL(r_kind) :: InLatLon(2)
    ! Local variables:
    INTEGER(i_kind) :: i, i_prof, it, nlevels, vLevel
    INTEGER(i_kind) :: nprofs, ninst, nsimobs
    INTEGER(i_kind) :: asw
    INTEGER(kind=jpim)  :: errstatus
    REAL(r_kind), ALLOCATABLE :: pres(:)

    INCLUDE 'rttov_tl.interface'
    INCLUDE 'rttov_init_emis_refl.interface'
    INCLUDE 'rttov_dealloc_coefs.interface'
    INCLUDE 'rttov_alloc_transmission.interface'
    INCLUDE 'rttov_alloc_rad.interface'

    vLevel = this%X%sg%vLevel
    InLatLon(1) = X%sg%cell_cntr(1,locX(1))
    InLatLon(2) = X%sg%cell_cntr(2,locX(1)) 

    CALL this%rttov_utils%rttov_GetMVars(X, locX)
    pres = X%sg%FGPres(:, locX(1), locX(2))

    nlevels = vLevel
    nprofs = 1 ! sg%num_icell or sg%num_cell
    nchans = 1
    nsimobs = nprofs*nchans

    asw = 1 !allocate

    CALL this%rttov_utils%alloc_tl(asw, nprofs, nchans, nlevels, ichan)

    CALL this%rttov_utils%rttov_MVar2RttovVar(nlevels, SatAngles, InLatLon)

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

    CALL this%rttov_utils%rttov_TLVars(pres, dX, locX, nlevels)
    DEALLOCATE(pres)

    asw = 1 ! allocate
    CALL rttov_tl(errstatus, &! out   error flag
                  this%rttov_utils%rttov_in_calc%chanprof, &! in LOCAL channel and profile index structure
                  this%rttov_utils%opts, &! in    options structure
                  this%rttov_utils%rttov_in_calc%profs, &! in    profile array
                  this%rttov_utils%rttov_in_calc%profs_tl, &! in    profile array
                  this%rttov_utils%coefs, &! in    coefficients structure
                  this%rttov_utils%rttov_in_calc%transm, &! inout computed transmittances
                  this%rttov_utils%rttov_in_calc%transm_tl, &! inout computed transmittances
                  this%rttov_utils%rttov_in_calc%radiance, &! inout computed radiances
                  this%rttov_utils%rttov_in_calc%radiance_tl, &! inout computed radiances
                  calcemis=this%rttov_utils%rttov_in_calc%calcemis, &! in    flag for internal emissivity calcs
                  emissivity=this%rttov_utils%rttov_in_calc%emis, &!, &! inout input/output emissivities per channel
                  emissivity_tl=this%rttov_utils%rttov_in_calc%emis_tl)                    ! inout input/output emissivities per channel

    IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%debug_dev) THEN
      PRINT *, 'Brightness Temperature (K) = '
      PRINT *, this%rttov_utils%rttov_in_calc%radiance%bt, this%rttov_utils%rttov_in_calc%radiance_tl%bt

      PRINT *, 'Radiance (mW/m2/sr/cm-1) = '
      PRINT *, this%rttov_utils%rttov_in_calc%radiance%total, this%rttov_utils%rttov_in_calc%radiance_tl%total

      PRINT *, 'Radiance overcast (mW/m2/sr/cm-1) = '
      PRINT *, this%rttov_utils%rttov_in_calc%radiance%cloudy, this%rttov_utils%rttov_in_calc%radiance_tl%cloudy

      PRINT *, '-------------Print profile_tl-----------'
      jchan = 1 ! print the 1st channel
      PRINT *, '================ Print from K ==================='
      PRINT *, "Level", "Pressure", "JAC(Temp) (K/K)", "JAC(WV) ", this%rttov_utils%rttov_in_calc%profs(1)%gas_units
      DO l = 1, this%rttov_utils%rttov_in_calc%profs_tl(jchan)%nlevels
        PRINT *, 'Level = ', l
        PRINT *, this%rttov_utils%rttov_in_calc%profs(1)%p(l)
        PRINT *, this%rttov_utils%rttov_in_calc%profs_tl(jchan)%t(l)
        PRINT *, this%rttov_utils%rttov_in_calc%profs_tl(jchan)%q(l)
      END DO

    END IF

  END SUBROUTINE rttov_tl_sp

  SUBROUTINE rttov_ad_sp(this, norm_d_sp, X, ichan, locX, SatAngles)
    CLASS(rttov_tlad_sp_t) :: this
    REAL(r_kind), INTENT(IN) :: norm_d_sp
    TYPE(State_t), INTENT(IN) :: X
    INTEGER(i_kind), INTENT(IN) :: ichan
    INTEGER(i_kind), INTENT(IN) :: locX(2)
    REAL(r_kind), INTENT(IN) :: SatAngles(4)
    REAL(r_kind) :: InLatLon(2)
    ! Local variables:
    INTEGER(i_kind) :: i, i_prof, it, nlevels, vLevel
    INTEGER(i_kind) :: nprofs, ninst, nsimobs
    INTEGER(i_kind) :: asw
    INTEGER(kind=jpim)  :: errstatus

    INCLUDE 'rttov_ad.interface'
    INCLUDE 'rttov_init_emis_refl.interface'
    INCLUDE 'rttov_dealloc_coefs.interface'
    INCLUDE 'rttov_alloc_transmission.interface'
    INCLUDE 'rttov_alloc_rad.interface'

    vLevel = this%X%sg%vLevel
    InLatLon(1) = X%sg%cell_cntr(1,locX(1))
    InLatLon(2) = X%sg%cell_cntr(2,locX(1)) 

    CALL this%rttov_utils%rttov_GetMVars(X, locX)

    nlevels = vLevel
    nprofs = 1 ! sg%num_icell or sg%num_cell
    nchans = 1
    nsimobs = nprofs*nchans

    asw = 1 !allocate

    CALL this%rttov_utils%alloc_ad(asw, norm_d_sp, nprofs, nchans, nlevels, ichan)

    CALL this%rttov_utils%rttov_MVar2RttovVar(nlevels, SatAngles, InLatLon)

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

    asw = 1 ! allocate
    CALL rttov_ad(errstatus, &! out   error flag
                  this%rttov_utils%rttov_in_calc%chanprof, &! in LOCAL channel and profile index structure
                  this%rttov_utils%opts, &! in    options structure
                  this%rttov_utils%rttov_in_calc%profs, &! in    profile array
                  this%rttov_utils%rttov_in_calc%profs_ad, &! in    profile array
                  this%rttov_utils%coefs, &! in    coefficients structure
                  this%rttov_utils%rttov_in_calc%transm, &! inout computed transmittances
                  this%rttov_utils%rttov_in_calc%transm_ad, &! inout computed transmittances
                  this%rttov_utils%rttov_in_calc%radiance, &! inout computed radiances
                  this%rttov_utils%rttov_in_calc%radiance_ad, &! inout computed radiances
                  calcemis=this%rttov_utils%rttov_in_calc%calcemis, &! in    flag for internal emissivity calcs
                  emissivity=this%rttov_utils%rttov_in_calc%emis, &!, &! inout input/output emissivities per channel
                  emissivity_ad=this%rttov_utils%rttov_in_calc%emis_ad)                    ! inout input/output emissivities per channel

    IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%debug_dev) THEN
      PRINT *, 'Brightness Temperature (K) = '
      PRINT *, this%rttov_utils%rttov_in_calc%radiance%bt

      PRINT *, 'Radiance (mW/m2/sr/cm-1) = '
      PRINT *, this%rttov_utils%rttov_in_calc%radiance%total

      PRINT *, 'Radiance overcast (mW/m2/sr/cm-1) = '
      PRINT *, this%rttov_utils%rttov_in_calc%radiance%cloudy

      PRINT *, '-------------Print profile_ad-----------'
      jchan = 1 ! print the 1st channel
      PRINT *, '================ Print from K ==================='
      PRINT *, "Level", "Pressure", "JAC(Temp) (K/K)", "JAC(WV) ", this%rttov_utils%rttov_in_calc%profs(1)%gas_units
      DO l = 1, this%rttov_utils%rttov_in_calc%profs_ad(jchan)%nlevels
        PRINT *, 'Level = ', l
        PRINT *, this%rttov_utils%rttov_in_calc%profs(1)%p(l)
        PRINT *, this%rttov_utils%rttov_in_calc%profs_ad(jchan)%t(l)
        PRINT *, this%rttov_utils%rttov_in_calc%profs_ad(jchan)%q(l)
      END DO

    END IF

  END SUBROUTINE rttov_ad_sp

  SUBROUTINE rttov_tl_sp_simobs(this, X, ichan, locX, SatAngles, dx_sp, dhx_sp)

    CLASS(rttov_tlad_sp_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    TYPE(State_t), TARGET, INTENT(IN) :: dx_sp ! delta x
    INTEGER(i_kind), INTENT(IN) :: ichan
    INTEGER(i_kind), INTENT(IN) :: locX(2)
    REAL(r_kind), INTENT(IN) :: SatAngles(4)
    REAL(r_kind), INTENT(INOUT) :: dhx_sp  ! delta_BT

    ! Local variables:
    INTEGER(i_kind) :: i, i_prof, it, nlevels, vLevel
    INTEGER(i_kind) :: ivar, nch, nlevs
    INTEGER(i_kind) :: nprofs, ninst, nsimobs, nchans
    INTEGER(i_kind) :: asw
    INTEGER(kind=jpim)  :: errstatus
    INTEGER(i_kind) :: stride, rttov_p_start, rttov_p_end, layer_start, layer_end
    REAL(r_kind) :: dhx_sp_cld
    REAL(r_kind), POINTER :: temp(:), qvapor(:), qcloud(:), qice(:), pres(:)
    LOGICAL :: use_kmodel

    use_kmodel = this%rttov_utils%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%use_kmodel
    nprofs = 1 ! this%sg%num_icell
    IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%debug_dev) PRINT *, 'TL nprofs =', nprofs
    nchans = 1
    IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%debug_dev) PRINT *, 'rttov_tl_sp_simobs: nchans = ', nchans
    ! TODO : repalce profs%t etc
    nsimobs = nprofs * nchans
    IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%debug_dev) PRINT *, 'TL nsimobs =', nsimobs

    IF ( .NOT. use_kmodel ) THEN
      CALL this%rttov_tl_sp(X, dx_sp, ichan, locX, SatAngles)
      dhx_sp = this%rttov_utils%rttov_in_calc%radiance_tl%bt(1)
      ! print *, 'debug after calling rttov_tl_sp: dhx_sp = ', dhx_sp

      asw = 0 !deallocate
      ! CALL this%rttov_utils%alloc_profs(asw, nprofs, nchans, X%sg%vLevel)
      CALL this%rttov_utils%zero_tl()
      CALL this%rttov_utils%alloc_tl(asw, nprofs, nchans, X%sg%vLevel, ichan)
    ELSE
      CALL this%rttov_k_sp_simobs(X, ichan, locX, SatAngles)

      nch = 0

      DO jch = 1, nchans

        nch = nch + 1
        i_prof = this%rttov_utils%rttov_in_calc%chanprof(nch)%prof

        dhx_sp = 0.0

        temp => dx_sp%fields(dx_sp%getVarIdx('temp'))%data(:, locX(1), locX(2))
        qvapor => dx_sp%fields(dx_sp%getVarIdx('qvapor'))%data(:, locX(1), locX(2))
        pres => X%sg%FGPres(:, locX(1), locX(2))
        IF ( this%rttov_utils%opts%rt_ir%addclouds ) &
          qcloud => dx_sp%fields(dx_sp%getVarIdx('qcloud'))%data(:, locX(1), locX(2))
        IF ( this%rttov_utils%opts%rt_ir%addclouds ) &
          qice => dx_sp%fields(dx_sp%getVarIdx('qice'))%data(:, locX(1), locX(2))
        ! psfc => dx_sp%fields(dx_sp%getVarIdx('psfc'))%data(1, locX(1), locX(2)) )

        nlevs = SIZE(temp, 1)

        IF (pres(1) > pres(2) .OR. pres(2) > pres(3)) THEN
          stride = -1
          rttov_p_start = nlevs
          rttov_p_end = 1
          layer_start = nlevs - 1
          layer_end = 1
        ELSE
          stride = 1
          rttov_p_start = 1
          rttov_p_end = nlevs
          layer_start = 1
          layer_end = nlevs - 1
        END IF
        
        ! print *, 'tl: ', rttov_p_start, rttov_p_end, stride
        ! PRINT *, 'Check cloud shape: ', SIZE(this%rttov_utils%rttov_in_calc%profs_k(nch)%cloud)
        IF (this%rttov_utils%opts%rt_ir%addclouds) THEN
          dhx_sp_cld = &
          sum(this%rttov_utils%rttov_in_calc%profs_k(nch)%cloud(1,1:nlevs-1) * qcloud(layer_start:layer_end:stride) ) + &
          sum(this%rttov_utils%rttov_in_calc%profs_k(nch)%cloud(6,1:nlevs-1) * qice(layer_start:layer_end:stride) )
        ELSE
          dhx_sp_cld = 0.0D0
        END IF

        IF (this%rttov_utils%rttov_setup_input%rttov_in_gasvars%ppmv_gasunit) THEN
          dhx_sp = dhx_sp + dhx_sp_cld + &
          sum(this%rttov_utils%rttov_in_calc%profs_k(nch)%t * temp(rttov_p_start:rttov_p_end:stride)) + &
          sum(this%rttov_utils%rttov_in_calc%profs_k(nch)%q * qvapor(rttov_p_start:rttov_p_end:stride) * q_mr_to_ppmv) !+ &
        ! this%rttov_utils%rttov_in_calc%profs_k(nch)%s2m%p * psfc
        ELSE 
          dhx_sp = dhx_sp + dhx_sp_cld + &
          sum(this%rttov_utils%rttov_in_calc%profs_k(nch)%t * temp(rttov_p_start:rttov_p_end:stride)) + &
          sum(this%rttov_utils%rttov_in_calc%profs_k(nch)%q * qvapor(rttov_p_start:rttov_p_end:stride)) !+ &
        ! this%rttov_utils%rttov_in_calc%profs_k(nch)%s2m%p * psfc
        END IF
        ! print *, 'debug after calling rttov_k: dhx_sp = ', dhx_sp
        
      END DO

      NULLIFY(qvapor, temp)
      IF ( this%rttov_utils%opts%rt_ir%addclouds ) NULLIFY(qcloud, qice)

      asw = 0 !deallocate

      ! CALL this%rttov_utils%alloc_profs(asw, nprofs, nchans, X%sg%vLevel)
      CALL this%rttov_utils%zero_k()
      CALL this%rttov_utils%alloc_k(asw, nprofs, nchans, X%sg%vLevel, ichan)

    END IF

    CALL this%rttov_utils%rttov_dealloc_array()

  END SUBROUTINE rttov_tl_sp_simobs

  SUBROUTINE rttov_ad_sp_simobs(this, X, ichan, locX, SatAngles, norm_d_sp, temp, qvapor, qcloud, qice)

    CLASS(rttov_tlad_sp_t) :: this
    TYPE(State_t), INTENT(IN) :: X
    REAL(r_kind), INTENT(IN) :: norm_d_sp  !O^(-1).(hx-y); d is innovation
    INTEGER(i_kind), INTENT(IN) :: ichan
    INTEGER(i_kind), INTENT(IN) :: locX(2)
    REAL(r_kind), INTENT(IN) :: SatAngles(4)
    REAL(r_kind), INTENT(INOUT) :: temp(:), qvapor(:), qcloud(:), qice(:)

    ! Local variables:
    INTEGER(i_kind) :: i, i_prof, it, nlevs, vLevel
    INTEGER(i_kind) :: ivar, nch
    INTEGER(i_kind) :: nprofs, ninst, nsimobs
    INTEGER(i_kind) :: asw
    INTEGER(kind=jpim)  :: errstatus
    INTEGER(i_kind) :: stride, rttov_p_start, rttov_p_end, layer_start, layer_end
    LOGICAL :: use_kmodel

    use_kmodel = this%rttov_utils%rttov_setup_input%rttov_in_userdefs%rttov_insts(1)%use_kmodel

    nprofs = 1 !this%sg%num_icell
    IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%debug_dev) PRINT *, 'AD nprofs =', nprofs

    nchans = 1
    ! TODO : repalce profs%t etc
    nsimobs = nprofs * nchans
    IF (this%rttov_utils%rttov_setup_input%rttov_in_userdefs%debug_dev) PRINT *, 'AD nsimobs =', nsimobs

    IF ( .NOT. use_kmodel ) THEN
      CALL this%rttov_ad_sp(norm_d_sp, X, ichan, locX, SatAngles)

      nch = 0
      DO jch = 1, nchans
        nch = nch + 1
        ASSOCIATE(pres => X%sg%FGPres(:, locX(1), locX(2)))
          
          nlevs = SIZE(temp,1)
          ! PRINT *, 'Check pres: ', pres(1), pres(2)
          IF (pres(1) > pres(2) .OR. pres(2) > pres(3)) THEN
            stride = -1
            rttov_p_start = nlevs
            rttov_p_end = 1
            layer_start = nlevs - 1
            layer_end = 1
          ELSE
            stride = 1
            rttov_p_start = 1
            rttov_p_end = nlevs
            layer_start = 1
            layer_end = nlevs - 1
          END IF
          temp = 0.0D0; qvapor = 0.0D0; qcloud = 0.0D0; qice = 0.0D0
          temp = this%rttov_utils%rttov_in_calc%profs_ad(nch)%t(rttov_p_start:rttov_p_end:stride)
          IF (this%rttov_utils%rttov_setup_input%rttov_in_gasvars%ppmv_gasunit) THEN
            qvapor = this%rttov_utils%rttov_in_calc%profs_ad(nch)%q(rttov_p_start:rttov_p_end:stride) * q_mr_to_ppmv
          ELSE 
            qvapor = this%rttov_utils%rttov_in_calc%profs_ad(nch)%q(rttov_p_start:rttov_p_end:stride)
          END IF
          IF (this%rttov_utils%opts%rt_ir%addclouds) THEN
            qcloud(1:nlevs-1) = this%rttov_utils%rttov_in_calc%profs_ad(nch)%cloud(1, layer_start:layer_end:stride)
            qice(1:nlevs-1) = this%rttov_utils%rttov_in_calc%profs_ad(nch)%cloud(6, layer_start:layer_end:stride)
          END IF
          ! print *, 'debug after calling rttov_ad_sp: temp = ', maxval(temp), minval(temp)
          ! print *, 'debug after calling rttov_ad_sp: qvapor = ', maxval(qvapor), minval(qvapor)
          ! print *, 'debug after calling rttov_ad_sp: qcloud = ', maxval(qcloud), minval(qcloud)
          ! print *, 'debug after calling rttov_ad_sp: qice = ', maxval(qice), minval(qice)
        END ASSOCIATE
      END DO

      asw = 0 !deallocate
      ! CALL this%rttov_utils%alloc_profs(asw, nprofs, nchans, X%sg%vLevel)
      CALL this%rttov_utils%zero_ad()
      CALL this%rttov_utils%alloc_ad(asw, norm_d_sp, nprofs, nchans, nlevels, ichan)

    ELSE

      CALL this%rttov_k_sp_simobs(X, ichan, locX, SatAngles)

      nch = 0

      ! TODO : repalce profs%t etc

      DO jch = 1, nchans
        nch = nch + 1

        ASSOCIATE(pres => X%sg%FGPres(:, locX(1), locX(2)))
          
          nlevs = SIZE(temp,1)
          ! PRINT *, 'Check pres: ', pres(1), pres(2)
          IF (pres(1) > pres(2) .OR. pres(2) > pres(3)) THEN
            stride = -1
            rttov_p_start = nlevs
            rttov_p_end = 1
            layer_start = nlevs - 1
            layer_end = 1
          ELSE
            stride = 1
            rttov_p_start = 1
            rttov_p_end = nlevs
            layer_start = 1
            layer_end = nlevs - 1
          END IF

          temp = 0.0D0; qvapor = 0.0D0; qcloud = 0.0D0; qice = 0.0D0
          temp = norm_d_sp * this%rttov_utils%rttov_in_calc%profs_k(nch)%t(rttov_p_start:rttov_p_end:stride)

          IF (this%rttov_utils%rttov_setup_input%rttov_in_gasvars%ppmv_gasunit) THEN
            qvapor = norm_d_sp * this%rttov_utils%rttov_in_calc%profs_k(nch)%q(rttov_p_start:rttov_p_end:stride) * q_mr_to_ppmv
          ELSE 
            qvapor = norm_d_sp * this%rttov_utils%rttov_in_calc%profs_k(nch)%q(rttov_p_start:rttov_p_end:stride)
          END IF
          IF (this%rttov_utils%opts%rt_ir%addclouds) THEN
            qcloud(1:nlevs-1) = norm_d_sp * this%rttov_utils%rttov_in_calc%profs_k(nch)%cloud(1, layer_start:layer_end:stride)
            qice(1:nlevs-1) = norm_d_sp * this%rttov_utils%rttov_in_calc%profs_k(nch)%cloud(6, layer_start:layer_end:stride)
          END IF
          ! psfc = norm_d_sp * this%rttov_utils%rttov_in_calc%profs_k(nch)%s2m%p

        END ASSOCIATE

      END DO

      asw = 0
      ! CALL this%rttov_utils%alloc_profs(asw, nprofs, nchans, X%sg%vLevel)
      CALL this%rttov_utils%zero_k()
      CALL this%rttov_utils%alloc_k(asw, nprofs, nchans, X%sg%vLevel, ichan)
    END IF

    CALL this%rttov_utils%rttov_dealloc_array()

  END SUBROUTINE rttov_ad_sp_simobs

  IMPURE ELEMENTAL SUBROUTINE destructor_tlad_sp(this)
    IMPLICIT NONE
    TYPE(rttov_tlad_sp_t), INTENT(INOUT) :: this

    CALL this%rttov_utils%rttov_dealloc_types()

  END SUBROUTINE destructor_tlad_sp

END MODULE rttov_tlad_sp_m
