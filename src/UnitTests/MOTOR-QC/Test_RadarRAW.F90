!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA.RadarRAW
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin, Yuanfu Xie
! VERSION           : V 0.1
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2022/1/4, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
PROGRAM Test_RadarRAW
  USE kinds_m, ONLY: i_kind, r_kind, r_double
  USE RadarRAW_m, ONLY: RadarRAW_t
  USE YAMLRead_m

  IMPLICIT NONE

  CHARACTER(LEN=1024) :: inputDir
  CHARACTER(LEN=1024) :: outputDir
  CHARACTER(len=1024) :: filenameInput, filenameOutput
  CHARACTER(len=1024) :: configFile
  CHARACTER(LEN=1024) :: staticDir
  INTEGER(i_kind) :: istatus

  TYPE(RadarRAW_t) :: RadarRAW
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", staticDir)
  configFile = TRIM(staticDir)//"/UnitTest"//"/testRadarRAW.yaml"

  ! Get the nc filenmae
  istatus = yaml_get_var(configFile, 'IO', 'input_dir', inputDir)
  filenameInput = TRIM(inputDir)//"/Obs"//"/Z_RADR_I_Z9200_20220217060001_O_DOR_SAD_CAP_FMT.nc"
  filenameInput = TRIM(inputDir)//"/Obs"//"/ZSE01_20220526225928.nc"
  filenameInput = "/public/home/simi/optest/3DAnalysis.alpha/task/221116/221116_2335/input/obs/radar/Z_RADR_I_Z9660_20221116232401_O_DOR_SAD_CAP_FMT.nc"
  filenameInput = "/public/home/simi/optest/Case0908/2022090809/input/obs/radar/Z_RADR_I_Z9755_20220908085401_O_DOR_SAD_CAP_FMT.nc"
  filenameInput = "/Users/qzl/sources/MOTOR/input/220527_0000_3km/input/obs/radar/Z_RADR_I_Z9200_20220526230600_O_DOR_SAD_CAP_FMT.nc"

  istatus = yaml_get_var(configFile, 'IO', 'output_dir', outputDir)
  filenameOutput = TRIM(outputDir)//"/testRadar.nc"

  ! Constructer for RadarRAW
  RadarRAW = RadarRAW_t(configFile)

  ! Read the radar data
  CALL RadarRAW%ReadRadarData(filenameInput, 'vel')

  ! Output the raw data for debug
  CALL RadarRAW%outputNCForTest(filenameOutput)
END PROGRAM Test_RadarRAW
