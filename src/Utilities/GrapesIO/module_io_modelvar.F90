MODULE module_io_modelvar
  USE kinds_m, ONLY: i_kind, r_kind, r_single
  USE NameValueMap_m

  IMPLICIT NONE

  TYPE :: module_io_modelvar_t
    CHARACTER(LEN=1024), ALLOCATABLE  ::  modvarFns(:)
    CHARACTER(LEN=1024)               ::  modvarCtlFn
    INTEGER(i_kind)                   ::  num_files, vLevel
    REAL(r_kind), ALLOCATABLE         ::  lat2D(:, :), lon2D(:, :), lat1D(:), lon1D(:)
    REAL(r_kind), ALLOCATABLE       ::  zRHght_u(:, :, :), zs(:, :)
    REAL(r_kind), ALLOCATABLE         ::  zRHght_s(:, :, :), topo(:, :)
    REAL(r_single), DIMENSION(:, :, :, :), ALLOCATABLE :: pi, &
                                                          th, &
                                                          u, &
                                                          v, &
                                                          w, &
                                                          qv, &
                                                          zz
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE   :: temp, &
                                                          pres, &
                                                          uwnd, &
                                                          vwnd, &
                                                          wwnd, &
                                                          qvapor, &
                                                          height
    REAL(r_single), DIMENSION(:, :, :), ALLOCATABLE    :: rainc_s, &
                                                          rainnc_s, &
                                                          t2m_s, &
                                                          q2m_s, &
                                                          u10m_s, &
                                                          v10m_s, &
                                                          vis_s, &
                                                          psl_s, &
                                                          psfc_s
    REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE      :: rainc, &
                                                          rainnc, &
                                                          t2m, &
                                                          q2m, &
                                                          u10m, &
                                                          v10m, &
                                                          vis, &
                                                          psl, &
                                                          psfc
    INTEGER(i_kind)              :: idn, jdn, kdn, tdn
    INTEGER(i_kind)              :: ids, ide, jds, jde, kds, kde
    INTEGER(i_kind), ALLOCATABLE :: time_unix(:), time_gmt(:, :)

  CONTAINS
    PROCEDURE :: init_modelvar
    PROCEDURE :: read_modelvar
    PROCEDURE :: write_modelvar
    PROCEDURE :: u_c2a
    PROCEDURE :: v_c2a
    PROCEDURE :: compute_t_p
    PROCEDURE :: compute_zHght_s
    PROCEDURE :: writeBackToGrapesInput


  END TYPE module_io_modelvar_t

  SUBROUTINE init_modelvar(hgrid, nlFileName, giFileName)
    CLASS (module_io_modelvar_t) :: this

  END SUBROUTINE init_modelvar

  SUBROUTINE read_modelvar(this)
    CLASS (module_io_modelvar_t) :: this
  END SUBROUTINE read_modelvar

  SUBROUTINE u_c2a(this)
    CLASS (module_io_modelvar_t) :: this
  END SUBROUTINE u_c2a

  SUBROUTINE v_c2a(this)
    CLASS (module_io_modelvar_t) :: this
  END SUBROUTINE v_c2a

  SUBROUTINE compute_t_p(this)
    CLASS (module_io_modelvar_t) :: this
  END SUBROUTINE compute_t_p

  SUBROUTINE compute_zHght_s(this)
    CLASS (module_io_modelvar_t) :: this
  END SUBROUTINE compute_zHght_s

  SUBROUTINE writeBackToGrapesInput(this)
    CLASS (module_io_modelvar_t) :: this
  END SUBROUTINE writeBackToGrapesInput

END MODULE module_io_modelvar
