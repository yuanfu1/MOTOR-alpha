!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.BMatrix
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/4/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
PROGRAM Test_BMatrix
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State_m, ONLY: State_t
  USE BMatrix_m, ONLY: BMatrix_t
  USE kinds_m, ONLY: i_kind, r_kind

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(BMatrix_t) :: BMat
  TYPE(State_t) :: Xm
  TYPE(State_t) :: XX

  CHARACTER(LEN=1024) :: configFile
  REAL(r_kind) :: J

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/UnitTest/testBMat.yaml"

  ! Initializer
  CALL mpddGlob%initialize()                              ! Initialize the mpdd
  CALL geometry%initialize(configFile, mpddGlob)          ! Initialize the geometry

  ASSOCIATE (sg => geometry%mg%sg(8))
    CALL Xm%initialize(configFile, sg)
    CALL BMat%initialize(configFile, sg)

    ! XX = B^(- rac {1} {2})*Xm
    XX = BMat.SQRTINVMUL.Xm

    ! Try dot product.
    J = (XX.DOT.XX)
    PRINT *, 'Programe is done! J is ', J
    ! PRINT *, 'X is: ', Xm%fields(1)%data(1,:)
  END ASSOCIATE
  ! Destroy the structures
  CALL mpddGlob%finalize

END PROGRAM Test_BMatrix
