
MODULE cloud_rad

!     Cloud Radiation and Microphysics Parameters

!     Backscattering efficiencies
  REAL, PARAMETER :: bksct_eff_clwc = .063
  REAL, PARAMETER :: bksct_eff_cice = .14
  REAL, PARAMETER :: bksct_eff_rain = .063
  REAL, PARAMETER :: bksct_eff_snow = .14
  REAL, PARAMETER :: bksct_eff_graupel = .30

!     Scattering efficiencies
  REAL, PARAMETER :: q_clwc = 2.0
  REAL, PARAMETER :: q_cice = 2.0
  REAL, PARAMETER :: q_rain = 1.0
  REAL, PARAMETER :: q_snow = 1.0
  REAL, PARAMETER :: q_graupel = 1.0

!     Densities
  REAL, PARAMETER :: rholiq = 1E3 ! kilograms per cubic meter
  REAL, PARAMETER :: rhosnow = .07E3 ! kilograms per cubic meter
  REAL, PARAMETER :: rhograupel = .50E3 ! kilograms per cubic meter

!     Effective radii
  REAL, PARAMETER :: reff_clwc = .000007 ! m
  REAL, PARAMETER :: reff_cice = .000034 ! m
  REAL, PARAMETER :: reff_rain = .000750 ! m
  REAL, PARAMETER :: reff_snow = .004000 ! m
  REAL, PARAMETER :: reff_graupel = .010000 ! m

!     GHI related
  REAL, PARAMETER :: ghi_zen_toa = 1361.5 ! solar const, W/m**2 at 1AU
  REAL, PARAMETER :: zen_kt = 0.815       ! zenithal attenuation of GHI
  PRIVATE stretch2
  PUBLIC albedo_to_clouds, albedo_to_clouds2

CONTAINS
  SUBROUTINE stretch2(IL, IH, JL, JH, rArg)

    IMPLICIT NONE

    REAL A, B, IL, IH, JL, JH, rarg

    a = (jh - jl) / (ih - il)
    b = jl - il * a

    rarg = a * rarg + b

    RETURN
  END SUBROUTINE stretch2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE albedo_to_clouds(albedo & ! I
                              , cloud_rad_trans, cloud_od, cloud_opacity)   ! O

    USE mem_namelist, ONLY: r_missing_data

    clear_albedo = .2097063
    cloud_albedo = .4485300
!       cloud_albedo = .40

    IF (albedo .EQ. r_missing_data) THEN
      cloud_rad_trans = r_missing_data
      cloud_od = r_missing_data
      cloud_opacity = r_missing_data
      RETURN
    END IF

!       Uniform clouds are presently assumed in the grid box or cloud layer
    cloud_frac = 1.0

    arg = albedo

    CALL stretch2(clear_albedo, cloud_albedo, 0., 1., arg)

    cloud_albedo = MIN(arg, 0.99999)      ! Back Scattered

    cloud_rad_trans = 1.0 - cloud_albedo ! Fwd Scattered + Direct Transmission

    bksc_eff_od = -LOG(cloud_rad_trans)  ! Tau * Back Scat Efficiency

    cloud_od = bksc_eff_od / 0.10        ! Tau

    cloud_opacity = 1.0 - EXP(-cloud_od) ! 1 - Direct Transmission

    RETURN

  END SUBROUTINE albedo_to_clouds

  SUBROUTINE albedo_to_clouds2(albedo & ! I
                               , cloud_trans_l, cloud_trans_i & ! O
                               , cloud_od_l, cloud_od_i & ! O
                               , cloud_opac_l, cloud_opac_i)               ! O

    USE mem_namelist, ONLY: r_missing_data

    IF (albedo .EQ. r_missing_data) THEN
      cloud_trans_l = r_missing_data
      cloud_trans_i = r_missing_data
      cloud_od_l = r_missing_data
      cloud_od_i = r_missing_data
      cloud_opac_l = r_missing_data
      cloud_opac_i = r_missing_data
      RETURN
    END IF

    cloud_albedo = MIN(albedo, 0.99999)     ! Back Scattered

    cloud_trans_l = 1.0 - cloud_albedo     ! Fwd Scattered + Direct Transmission
    cloud_trans_i = 1.0 - cloud_albedo     ! Fwd Scattered + Direct Transmission

    cloud_od_l = cloud_albedo / (bksct_eff_clwc * (1.-cloud_albedo)) ! Tau
    cloud_od_i = cloud_albedo / (bksct_eff_cice * (1.-cloud_albedo)) ! Tau

    cloud_opac_l = 1.0 - EXP(-cloud_od_l) ! 1 - Direct Transmission
    cloud_opac_i = 1.0 - EXP(-cloud_od_i) ! 1 - Direct Transmission

    RETURN

  END SUBROUTINE albedo_to_clouds2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END MODULE cloud_rad
