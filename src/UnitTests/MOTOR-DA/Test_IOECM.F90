PROGRAM Test_IOECM
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE singleGrid_m, ONLY: singleGrid_t
  USE State_m, ONLY: State_t
  USE IOECM_m, ONLY: IOECM_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE YAMLRead_m
  USE parameters_m, ONLY: degree2radian

  IMPLICIT NONE

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t) :: geometry
  TYPE(singleGrid_t) :: sg
  TYPE(State_t) :: Xm
  TYPE(IOECM_t) :: ioecm
  CHARACTER(LEN=1024) :: configFile, outputDir

  configFile = "/home/mgq/proj/MOTOR/input/hagongda/App_SurfaceAnalysis.yaml"
  outputDir = "/home/mgq/proj/MOTOR/input/hagongda/output"

  ! 初始化几何和网格
  PRINT *, 'Test_IOECM0: Initializing geometry and singleGrid objects...'
  CALL mpddGlob%initialize()
  CALL geometry%initialize(configFile, mpddGlob)

  ASSOCIATE (sg => geometry%mg%sg(7))

    ! 初始化状态变量
    CALL Xm%initialize(configFile, sg)
    PRINT *, 'Test_IOECM1: Initialized State_t object.'
    ! 初始化 IOECM 模块
    CALL ioecm%initialize(configFile, geometry)
    PRINT *, 'Test_IOECM2: Initialized IOECM_t object.'

    ! 读取背景数据并加载到状态变量中
    CALL ioecm%m_read_bcg_into_Xm(Xm, sg)

    ! 打印一些调试信息
    PRINT *, 'Test_IOECM5: Background data loaded successfully.'

  END ASSOCIATE
  ! 清理资源
  ! CALL ioecm%destructor()

END PROGRAM Test_IOECM
