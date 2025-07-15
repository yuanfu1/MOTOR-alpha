PROGRAM Namelist_Read_TEST

  USE NMLRead_m

  IMPLICIT NONE
  REAL                :: test_real4
  INTEGER*4           :: test_int4
  INTEGER*8           :: test_int8
  LOGICAL             :: test_bool
  CHARACTER(LEN=1024) :: solver
  CHARACTER(LEN=1024) :: NMLFilename
  CHARACTER(LEN=1024) :: fNMLFilename
  CHARACTER(LEN=10), ALLOCATABLE :: varList(:)
  CHARACTER(LEN=1024), ALLOCATABLE :: bakFilesNames(:)
  CHARACTER(LEN=1024) :: bakFilesDir
  INTEGER*4           :: bakFilesNum, i

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", NMLFilename)
  NMLFilename = TRIM(NMLFilename)//"/UnitTest/testNMLRead.nl"
  PRINT *, 'NMLFilename is: ', TRIM(NMLFilename)

  CALL namelist_read(TRIM(NMLFilename), "vLevel", test_int4)
  PRINT *, 'vLevel is: ', test_int4

  CALL namelist_read(TRIM(NMLFilename), 'solver', solver)
  PRINT *, 'solver is: ', TRIM(solver)

  !CALL namelist_read(TRIM(NMLFilename), 'NMLFileName', fNMLFilename)
  fNMLFilename = nml_get_string(TRIM(NMLFilename), 'NMLFileName')
  PRINT *, 'NMLFileName is: ', TRIM(fNMLFilename)

  CALL namelist_read(TRIM(NMLFilename), 'varList', varList)
  PRINT *, 'varList is: ', varList

  CALL GET_ENVIRONMENT_VARIABLE("SRC_DIR", bakFilesDir)
  bakFilesDir = TRIM(bakFilesDir)//"/MOTOR-DA/IO/IOModels/IOWRF/"
  PRINT *, 'bakFilesDir is: ', TRIM(bakFilesDir)

  CALL namelist_read(TRIM(NMLFilename), 'bakFileList', bakFilesNames)
  bakFilesNum = UBOUND(bakFilesNames, 1)
  DO i = 1, bakFilesNum
    bakFilesNames(i) = TRIM(bakFilesDir)//TRIM(bakFilesNames(i))
    PRINT *, 'bakFilesNames is: ', TRIM(bakFilesNames(i))
  END DO

  IF (test_int4 .EQ. 65) THEN
    PRINT *, 'Test passed!'
  END IF

END PROGRAM Namelist_Read_TEST
