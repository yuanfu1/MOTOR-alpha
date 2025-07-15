!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Xuancheng LIU
! VERSION           : V 0.0
! HISTORY           :
!   Created by Xuancheng Liu (lxclit@outlook.com), 2022/4/4, @SZSC, Shenzhen
!----------------------------------------------------------------------------------------

PROGRAM Test_yamlRead
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: OUTPUT_UNIT
  USE yaml!!, only: parse, error_length
  USE yaml_types!!, only: type_node, type_dictionary, type_error, real_kind, &
                !!        type_list, type_list_item, type_scalar

  CLASS(type_node), POINTER :: root
  CHARACTER(len=error_length) :: error

  CLASS(type_dictionary), POINTER :: dict, temp
  CLASS(type_list), POINTER :: list
  CLASS(type_list_item), POINTER :: item
  TYPE(type_error), POINTER :: io_err

  ! character(len=:), allocatable :: string
  CHARACTER(len=1024) :: configFile
  CHARACTER(len=10) :: string
  CHARACTER(len=1024) :: inputDir

  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//'/Test_yamlRead.yaml'

  root => parse(TRIM(configFile), error=error)
  IF (error /= '') THEN
    PRINT *, TRIM(error)
    STOP 1
  END IF

  SELECT TYPE (root)
  CLASS is (type_dictionary)

! !> INTEGER: icosTr_grid: glvl
!     dict => root%get_dictionary('icosTr_grid', required=.false., error=io_err)
!     print *, 'INTEGER: icosTr_grid: glvl: ', dict%get_integer('glvl', error=io_err)

! !> INTEGER from DICTIONARY: latlon_grid: num_grid: x
!     dict => root%get_dictionary('latlon_grid', required=.false., error=io_err)
!     temp => dict%get_dictionary('num_grid', required=.false., error=io_err)
!     print *, 'INTEGER from DICTIONARY: latlon_grid: num_grid: x: ', temp%get_integer('x', error=io_err)

! !> REAL from DICTIONARY: latlon_grid: domain_latlon: lat_s
!     temp => dict%get_dictionary('domain_latlon', required=.false., error=io_err)
!     print *, 'REAL from DICTIONARY: latlon_grid: domain_latlon: lat_s: ', temp%get_real('lat_s', error=io_err)

! !> STRING: poissonSovler: solver
!     dict => root%get_dictionary('poissonSovler', required=.false., error=io_err)
!     string = trim(dict%get_string('solver', error=io_err))
!     print *, 'STRING: poissonSovler: solver: ', string

! !> STRING: flog: LogFileName
!     dict => root%get_dictionary('flog', required=.false., error=io_err)
!     string = trim(dict%get_string('LogFileName', error=io_err))
!     print *, 'STRING: flog: LogFileName: ', string

! !> STRING from LIST: IO: ModelFileName
!     dict => root%get_dictionary('IO', required=.false., error=io_err)
!     list => dict%get_list('ModelFileName', required=.false., error=io_err)
!     item => list%first
!     do while (associated(item))
!       select type (element => item%node)
!       class is (type_scalar)
!         ! string = repeat(" ",128)    !! allocate string, then string can read
!         string = trim(element%string)
!         print *, 'STRING from LIST: IO: ModelFileName: ', string
!         item => item%next
!       end select
!     end do

    ! !> STRING from LIST: IO: ModelFileName
    ! dict => root%get_dictionary('analysis_para', required=.false., error=io_err)
    ! list => dict%get_list('start_time', required=.false., error=io_err)
    ! item => list%first
    ! do while (associated(item))
    !   select type (element => item%node)
    !   class is (type_scalar)
    !     ! string = repeat(" ",128)    !! allocate string, then string can read
    !     ! string = trim(element%string)
    !     print *, 'STRING from LIST: IO: ModelFileName: ', element%to_integer(0)
    !     item => item%next
    !   end select
    ! end do

  END SELECT

  CALL root%finalize()
  DEALLOCATE (root)

  BLOCK
    USE YAMLRead_m
    INTEGER, ALLOCATABLE :: array1(:)
    REAL(8), ALLOCATABLE :: array2(:)
    INTEGER*2 :: value1
    REAL(8) :: value2
    INTEGER*4, ALLOCATABLE :: array4_1(:)
    REAL(4), ALLOCATABLE :: array4_2(:)
    INTEGER*4 :: value4_1
    REAL(16) :: value4_2
    CHARACTER(len=10), ALLOCATABLE :: vars(:)
    CHARACTER(len=1024), ALLOCATABLE :: files(:)
    INTEGER :: i, ifile
    REAL(16), ALLOCATABLE :: array_PX(:)
    CLASS(*), ALLOCATABLE :: test(:)
    LOGICAL :: isOK

    PRINT *, 'yaml_get_var tests: '

    ifile = yaml_get_var(configFile, 'poissonSovler', 'nCycle', value4_1)
    PRINT *, 'Result is: ', value4_1, ', get status: ', ifile

    ifile = yaml_get_var(configFile, 'poissonSovler', 'nCycle', value4_2)
    PRINT *, 'Result is: ', value4_2, ', get status: ', ifile

    ifile = yaml_get_var(configFile, 'IO', 'ModelFileName', files)
    DO i = 1, SIZE(files)
      PRINT *, 'Result is: ', TRIM(files(i)), i, ', get status: ', ifile
    END DO

    ifile = yaml_get_var(configFile, 'poissonSovler', 'omegas', array_PX)
    DO i = 1, SIZE(array_PX)
      PRINT *, 'Result is: ', (array_PX(i)), i, ', get status: ', ifile
    END DO
    ! STOP

    ifile = yaml_get_var(configFile, 'IO', 'input_dir', inputDir)
    PRINT *, 'Result is: ', TRIM(inputDir), ifile

    ifile = yaml_get_var('/Users/qzl/sources/MOTOR/input/2024110300_T10p1/ObsSelection.yaml',&
     'G3', 'isDA', isOK, dictName2 = 'SYNOP', dictName3 = 'uwnd')
     print*, 'Result is: ', isOK, ifile

    PRINT *, "YAML tests succeed!!"

  END BLOCK

END PROGRAM

