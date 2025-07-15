!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/3/3, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
!! @author Zilong Qin
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
PROGRAM Test_IntepHP
  USE geometry_m, ONLY: geometry_t
  USE mpddGlob_m, ONLY: mpddGlob_t
  USE State_m, ONLY: State_t
  USE State2NC_m
  USE GrapesIO_m, ONLY: GrapesIO_t
  USE InterpHPNew_m, ONLY: InterpHPNew_t
  USE kinds_m, ONLY: r_single
  USE YAMLRead_m

  ! Define types
  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(geometry_t), TARGET :: geometry
  TYPE(State_t) :: Xm
  REAL(r_kind) :: t1, t2
  TYPE(GrapesIO_t)  :: GrapesIO
  INTEGER(i_kind) :: ifile

  CHARACTER(LEN=1024) :: configFile, outputDir

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/testInterpHP.yaml"

  ifile = yaml_get_var(configFile, 'IO', 'output_dir', outputDir)

  ! Initialize
  CALL mpddGlob%initialize()  ! Initialize the mpdd
  CALL geometry%initialize(configFile, mpddGlob)          ! Initialize the geometry

  ASSOCIATE (sg => geometry%mg%sg(10))

    ! Initialize the Xm
    CALL Xm%initialize(configFile, sg)

    CALL mpddGlob%barrier
    CALL CPU_TIME(t1)

    BLOCK
      USE parameters_m, ONLY: degree2radian

      REAL(r_kind), ALLOCATABLE :: cellCntrModel(:, :)
      INTEGER(i_kind) :: i, j, k, HorNumGridModel
      TYPE(InterpHPNew_t) :: InterpHPNew

      IF (sg%isBaseProc()) THEN

        GrapesIO = GrapesIO_t('/home/liuxuancheng/WORK-SPACE/SUB-SPACE/MOTOR/input/namelist.input', &
                              '/home/liuxuancheng/WORK-SPACE/SUB-SPACE/MOTOR/input/grapesinput-4dv-2021053006')

        HorNumGridModel = GrapesIO%idn * GrapesIO%jdn; 
        ALLOCATE (cellCntrModel(2, HorNumGridModel))

        cellCntrModel(1, :) = GrapesIO%lat * degree2radian
        cellCntrModel(2, :) = GrapesIO%lon * degree2radian
      END IF

      PRINT *, SIZE(cellCntrModel), SIZE(cellCntrModel(1, :))
      CALL sg%mpddInfo_sg%bCast(HorNumGridModel)

      IF (sg%isActiveProc() .AND. (.NOT. sg%isBaseProc())) THEN
        ALLOCATE (cellCntrModel(2, HorNumGridModel))
      END IF
      CALL sg%mpddInfo_sg%bCast(cellCntrModel)

      InterpHPNew = InterpHPNew_t(geometry%mpdd, sg, cellCntrModel, HorNumGridModel)

      BLOCK
        REAL(r_single) :: zRHght_s(GrapesIO%kde_s - GrapesIO%kds_s + 1, HorNumGridModel), tmp(HorNumGridModel)
        REAL(r_single) :: zRHght_s_temp(GrapesIO%kde_s - GrapesIO%kds_s + 1, HorNumGridModel)

        ! PRINT *, 'Here generate the model height at location of MOTOR-DA cell'
        IF (sg%isBaseProc()) THEN
          vLevel_s = GrapesIO%kde_s - GrapesIO%kds_s + 1
          vLevel_u = GrapesIO%kde_u - GrapesIO%kds_u + 1
          ! FORALL (i=1:vLevel_s) zRHght_s(i, :) = &
          !   RESHAPE(GrapesIO%zRHght_s(:, i + GrapesIO%kds_s - 1, :), (/HorNumGridModel/))

          DO i = 1, GrapesIO%idn
            DO j = 1, GrapesIO%jdn
              zRHght_s_temp(:, (j - 1) * GrapesIO%idn + i) = GrapesIO%zRHght_s(i, GrapesIO%kds_s:GrapesIO%kde_s, j)
            END DO
          END DO
          zRHght_s = zRHght_s_temp
          tmp = RESHAPE(GrapesIO%hgrid%ht(:, :), (/HorNumGridModel/))
        END IF

        PRINT *, SHAPE(zRHght_s)

        CALL InterpHPNew%interpN(sg%vLevel, zRHght_s, Xm%fields(1)%DATA(:, :, 1))
        CALL InterpHPNew%interp_singleLevel_s(tmp, Xm%fields(2)%DATA(1, :, 1))
      END BLOCK

    END BLOCK

    CALL mpddGlob%barrier
    CALL CPU_TIME(t2)
    IF (mpddGlob%isBaseProc()) PRINT *, 'Time cost:', t2 - t1, 's'

    ! Output Xm to NC files
    CALL Output_NC_State_AV(Xm, outputDir, 'Test_IntepHP')

  END ASSOCIATE

  ! Destroy the structures
  CALL mpddGlob%finalize

END PROGRAM Test_IntepHP
