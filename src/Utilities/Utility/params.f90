MODULE parameters_m

  !>
  !!=================================================================
  !!  This module defines a set of meteorological constant parameters
  !!  used in the models.
  !!
  !!  uthor Yuanfu Xie
  !!  \History: 2018-2
  !!=================================================================
  !
  USE kinds_m, ONLY: i_kind, r_kind, r_double

  REAL(r_double), PARAMETER :: g = 9.80665D0, f = 1.0D-4
  REAL(r_double), PARAMETER :: ae = 6371220.0D0        ! Radius of the Earth
  REAL(r_double), PARAMETER :: EarthRadius = 6371220.0D0  ! Radius of the Earth
  REAL(r_double), PARAMETER :: dry_air_gas_const = 287.058D0 ! Unit: J/(kg*K)
  REAL(r_double), PARAMETER :: spec_heat_const_vol = 717D0 ! Unit: J/(Kg*K) specific heat at constant volumes
  REAL(r_double), PARAMETER :: spec_heat_const_pres = 1004D0 ! Unit J/(Kg*K) specific heat at constant pressures
  REAL(r_double), PARAMETER :: surface_ref_pres = 1000.0D0 !1000.0 in Grapes model, 1013.0D0 in Grapes book ! Unit hPa surface reference pressure

  ! The Earth rotating rate: 2*pi/24 hours/3600 seconds
  REAL(r_double), PARAMETER :: Omega = 7.292D-05

  ! pi:
  REAL(r_double), PARAMETER :: pi = 3.141592653589793D0

  ! pi/180:
  REAL(r_double), PARAMETER :: degree2radian = 0.01745329252D0
  REAL(r_double), PARAMETER :: radian2degree = 57.2957795132D0
  REAL(r_double), PARAMETER :: machineEps = 1.1102230246251565D-016

  ! Coriolis parameters:
  REAL(r_double), PARAMETER :: coriolis(2) = (/1.0D-04, 1.0D-11/)

  ! For RTTOV usage:
  ! kg/kg to ppmv
  REAL(r_double), PARAMETER, PUBLIC :: q_mr_to_ppmv = 1.60771704E+6  ! kg/kg to ppmv
  REAL(r_double), PARAMETER, PUBLIC :: o3_mr_to_ppmv = 6.03504E+5
  REAL(r_double), PARAMETER, PUBLIC :: co2_mr_to_ppmv = 6.58114E+5
  REAL(r_double), PARAMETER, PUBLIC :: co_mr_to_ppmv = 1.0340699E+6
  REAL(r_double), PARAMETER, PUBLIC :: h2o_mr_to_ppmv = 6.58090E+5
  REAL(r_double), PARAMETER, PUBLIC :: ch4_mr_to_ppmv = 1.80548E+6
  REAL(r_double), PARAMETER, PUBLIC :: Pa2hPa = 0.01D0
  REAL(r_double), PARAMETER, PUBLIC :: hPa2Pa = 100.0D0
  REAL(r_double), PARAMETER, PUBLIC :: g_to_kg = 0.001D0
  REAL(r_double), PARAMETER, PUBLIC :: kg_to_g = 1000.0D0
  REAL(r_double), PARAMETER, PUBLIC :: KelvinsT0 = 273.15D0
  REAL(r_double), PARAMETER, PUBLIC :: km2m = 1.0E3
  REAL(r_double), PARAMETER, PUBLIC :: qv_limit = 1.0E-10, rhov_ctl_limit = 1.0E-10
  REAL(r_double), PARAMETER, PUBLIC :: ZERO = 0.0D0

  ! Air parcel constants
  REAL(r_double), PARAMETER :: r_d = 287.0D0       ! Gas constant
  REAL(r_double), PARAMETER :: cp = 7.0D0 * r_d / 2.0D0  ! Air Heat capacity at constant  pressure
  REAL(r_double), PARAMETER :: k_d = r_d / cp      ! Adiabatic index
  REAL(r_double), PARAMETER :: r_k_d = 1.0D0 / k_d     ! Reverse of Adiabatic index
  REAL(r_double), PARAMETER :: mwDry = 2.896D01    ! dry air molecular weight
  REAL(r_double), PARAMETER :: gamma_d = g / cp        ! dry lapse rate
  REAL(r_double), PARAMETER :: standard_Pres = 1.0D03      ! standard presure level
  REAL(r_double), PARAMETER :: TK = 2.7315D02   ! standard presure level

  ! Water vapor constants
  REAL(r_double), PARAMETER :: Lv = 2.5D06       ! Latent heat ratio
  REAL(r_double), PARAMETER :: mwWater = 1.802D01    ! water molecular weight
  REAL(r_double), PARAMETER :: epsw_r = mwWater / mwDry

  ! Standard error/out unit
  INTEGER(i_kind), PARAMETER, PUBLIC :: stderr = 12
  INTEGER(i_kind), PARAMETER, PUBLIC :: stdout = 13
  ! Set Missing Value
  REAL(r_double), PARAMETER :: missing = 99999999.0D0, invalid = 3.4D38
  REAL(r_double), PARAMETER :: rho_air = 1.29D0

CONTAINS

  SUBROUTINE machine
    DOUBLE PRECISION :: eps

    eps = 1.0D0

    DO WHILE (1.0D0 + eps .NE. 1.0D0)
      eps = eps / 2.0D0
    END DO
    PRINT *, 'Machine epsilon: ', eps
  END SUBROUTINE machine
END MODULE parameters_m
