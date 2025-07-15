!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com.com), 2021/10/27, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!! This is a plotting routine to check the analytic solution in the unitTest
!
PROGRAM plot_sfc
  USE kinds_m, ONLY: i_kind, r_kind
  USE unitTest_sfc_m, ONLY: unitTest_sfc_t
  USE State_m, ONLY: State_t
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State2NC_m
  USE YAMLRead_m

  IMPLICIT NONE

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(State_t) :: X
  TYPE(unitTest_sfc_t) :: unitTest

  INTEGER(i_kind) :: n, i, j
  REAL(r_kind), ALLOCATABLE :: func(:), xyt(:, :)

  CHARACTER(LEN=1024) :: configFile, ncOutputFile

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/Application/App_3DVAR_SimpleCase.yaml"

  IF (yaml_get_var(configFile, 'IO', 'output_dir', ncOutputFile) /= 0) STOP
  ncOutputFile = TRIM(ncOutputFile)//"/"

  ! Initializer
  ! Auxtypes
  CALL mpddGlob%initialize()

  CALL mpddGlob%barrier

  ! Initialize geometry
  CALL geometry%initialize(configFile, mpddGlob)

  ASSOCIATE (sg => geometry%mg%sg(8))
    n = sg%num_cell
    CALL X%initialize(configFile, sg)       ! Initialize the state
    ! Assume just plot one time frame:
    ALLOCATE (func(n), xyt(3, n))

    BLOCK
      REAL(r_kind) :: minGlobLat, maxGlobLat, minGlobLon, maxGlobLon, swap

      swap = MINVAL(sg%cell_cntr(1, 1:sg%num_icell))
      CALL sg%mpddInfo_sg%AllReduceMinReal(swap, minGlobLat)

      swap = MAXVAL(sg%cell_cntr(1, :))
      CALL sg%mpddInfo_sg%AllReduceMaxReal(swap, maxGlobLat)

      swap = MINVAL(sg%cell_cntr(2, :))
      CALL sg%mpddInfo_sg%AllReduceMinReal(swap, minGlobLon)

      swap = MAXVAL(sg%cell_cntr(2, :))
      CALL sg%mpddInfo_sg%AllReduceMaxReal(swap, maxGlobLon)

      PRINT *, minGlobLat, maxGlobLat, minGlobLon, maxGlobLon

      xyt(1, :) = (sg%cell_cntr(1, :) - minGlobLat) / (maxGlobLat - minGlobLat)
      xyt(2, :) = (sg%cell_cntr(2, :) - minGlobLon) / (maxGlobLon - minGlobLon)
      xyt(3, :) = 0.5D0

    END BLOCK

    CALL unitTest%analytic(n, xyt, X%fields(1)%DATA)
    PRINT *, 'Max/min analytic function: ', MAXVAL(X%fields(1)%DATA), MINVAL(X%fields(1)%DATA)
  END ASSOCIATE

  CALL mpddGlob%barrier

  ! Output state file to NC for view.
  CALL Output_NC_State_SV(X, ncOutputFile, 'Plot', "temp")

  DEALLOCATE (func, xyt)

  ! Finalize
  CALL mpddGlob%finalize

END PROGRAM
