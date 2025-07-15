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

MODULE YAMLRead_m
  USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: OUTPUT_UNIT
  USE yaml!!, only: parse, error_length
  USE yaml_types!!, only: type_node, type_dictionary, type_error, real_kind, &
                !!        type_list, type_list_item, type_scalar
  USE String_Functions
  IMPLICIT NONE

  CLASS(type_node), POINTER, PRIVATE :: root => NULL()
  LOGICAL, PRIVATE :: root_parse = .FALSE.
  CHARACTER(len=1024), PRIVATE :: YAMLName = ""

  INTERFACE yaml_get_var
    PROCEDURE readRealXArray
    PROCEDURE readReal8Array
    PROCEDURE readReal4Array
    PROCEDURE readInt8Array
    PROCEDURE readInt4Array
    PROCEDURE readInt2Array
    PROCEDURE readRealXVal
    PROCEDURE readReal8Val
    PROCEDURE readReal4Val
    PROCEDURE readInt8Val
    PROCEDURE readInt4Val
    PROCEDURE readInt2Val
    PROCEDURE readStringArray
    PROCEDURE readString
    PROCEDURE readLogicalArray
    PROCEDURE readLogical
  END INTERFACE

CONTAINS

  ! 辅助函数：确保YAML文件已解析
  SUBROUTINE ensure_parsed(fileName)
    CHARACTER(len=*), INTENT(IN) :: fileName
    CHARACTER(len=error_length) :: error

    IF (.NOT. ASSOCIATED(root) .OR. root_parse .OR. TRIM(YAMLName) .NE. TRIM(fileName)) THEN
      root => parse(TRIM(fileName), error=error)
      IF (error /= '') THEN
        PRINT *, TRIM(error)
        STOP 1
      END IF
      YAMLName = TRIM(fileName)
    END IF
  END SUBROUTINE ensure_parsed

  ! 处理环境变量替换
  SUBROUTINE process_env_vars(string)
    CHARACTER(len=*), INTENT(INOUT) :: string
    CHARACTER(len=1024) :: STATIC_DIR, OUTPUT_DIR, INPUT_DIR, LOG_DIR

    CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", STATIC_DIR)
    CALL GET_ENVIRONMENT_VARIABLE("OUTPUT_DIR", OUTPUT_DIR)
    CALL GET_ENVIRONMENT_VARIABLE("INPUT_DIR", INPUT_DIR)
    CALL GET_ENVIRONMENT_VARIABLE("LOG_DIR", LOG_DIR)

    string = TRIM(Replace_Text(string, '${STATIC_DIR}', TRIM(STATIC_DIR)))
    string = TRIM(Replace_Text(string, '${OUTPUT_DIR}', TRIM(OUTPUT_DIR)))
    string = TRIM(Replace_Text(string, '${INPUT_DIR}', TRIM(INPUT_DIR)))
    string = TRIM(Replace_Text(string, '${LOG_DIR}', TRIM(LOG_DIR)))
  END SUBROUTINE process_env_vars

  FUNCTION readRealXArray(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    REAL(16), ALLOCATABLE :: val(:)
    CLASS(*), ALLOCATABLE :: valTemp(:)

    fileStat = 1
    ALLOCATE (REAL(16) :: valTemp(0))
    CALL ensure_parsed(fileName)
    CALL getArrayFromYAML(fileName, dictName1, valueName, valTemp, fileStat, dictName2, dictName3)
    SELECT TYPE (valTemp)
    TYPE IS (REAL(16))
        val = valTemp
    END SELECT
    DEALLOCATE (valTemp)

  END FUNCTION readRealXArray

  FUNCTION readRealXVal(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    REAL(16) :: val

    fileStat = 1
    CALL ensure_parsed(fileName)
    CALL getVarFromYAML(fileName, dictName1, valueName, val, fileStat, dictName2, dictName3)
  END FUNCTION readRealXVal

  FUNCTION readReal8Array(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    REAL(8), ALLOCATABLE :: val(:)
    CLASS(*), ALLOCATABLE :: valTemp(:)

    fileStat = 1
    ALLOCATE (REAL(8) :: valTemp(0))
    CALL ensure_parsed(fileName)
    CALL getArrayFromYAML(fileName, dictName1, valueName, valTemp, fileStat, dictName2, dictName3)
    SELECT TYPE (valTemp)
    TYPE IS (REAL(8))
        val = valTemp
    END SELECT
    DEALLOCATE (valTemp)
  END FUNCTION readReal8Array

  FUNCTION readReal8Val(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    REAL(8) :: val

    fileStat = 1
    CALL ensure_parsed(fileName)
    CALL getVarFromYAML(fileName, dictName1, valueName, val, fileStat, dictName2, dictName3)
  END FUNCTION readReal8Val

  FUNCTION readReal4Array(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    REAL(4), ALLOCATABLE :: val(:)
    CLASS(*), ALLOCATABLE :: valTemp(:)

    fileStat = 1
    ALLOCATE (REAL(4) :: valTemp(0))
    CALL ensure_parsed(fileName)
    CALL getArrayFromYAML(fileName, dictName1, valueName, valTemp, fileStat, dictName2, dictName3)
    SELECT TYPE (valTemp)
    TYPE IS (REAL(4))
        val = valTemp
    END SELECT
    DEALLOCATE (valTemp)
  END FUNCTION readReal4Array

  FUNCTION readReal4Val(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    REAL(4) :: val

    fileStat = 1
    CALL ensure_parsed(fileName)
    CALL getVarFromYAML(fileName, dictName1, valueName, val, fileStat, dictName2, dictName3)
  END FUNCTION readReal4Val

  FUNCTION readInt8Array(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    INTEGER(8), ALLOCATABLE :: val(:)
    CLASS(*), ALLOCATABLE :: valTemp(:)

    fileStat = 1
    ALLOCATE (INTEGER(8) :: valTemp(0))
    CALL ensure_parsed(fileName)
    CALL getArrayFromYAML(fileName, dictName1, valueName, valTemp, fileStat, dictName2, dictName3)
    SELECT TYPE (valTemp)
    TYPE IS (INTEGER(8))
        val = valTemp
    END SELECT
    DEALLOCATE (valTemp)
  END FUNCTION readInt8Array

  FUNCTION readInt8Val(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    INTEGER(8) :: val

    fileStat = 1
    CALL ensure_parsed(fileName)
    CALL getVarFromYAML(fileName, dictName1, valueName, val, fileStat, dictName2, dictName3)
  END FUNCTION readInt8Val

  FUNCTION readInt4Array(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    INTEGER(4), ALLOCATABLE :: val(:)
    CLASS(*), ALLOCATABLE :: valTemp(:)

    fileStat = 1
    ALLOCATE (INTEGER(4) :: valTemp(0))
    CALL ensure_parsed(fileName)
    CALL getArrayFromYAML(fileName, dictName1, valueName, valTemp, fileStat, dictName2, dictName3)
    SELECT TYPE (valTemp)
    TYPE IS (INTEGER(4))
        val = valTemp
    END SELECT
    DEALLOCATE (valTemp)
  END FUNCTION readInt4Array

  FUNCTION readInt4Val(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    INTEGER(4) :: val

    fileStat = 1
    CALL ensure_parsed(fileName)
    CALL getVarFromYAML(fileName, dictName1, valueName, val, fileStat, dictName2, dictName3)
  END FUNCTION readInt4Val

  FUNCTION readInt2Array(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    INTEGER(2), ALLOCATABLE :: val(:)
    CLASS(*), ALLOCATABLE :: valTemp(:)

    fileStat = 1
    ALLOCATE (INTEGER(2) :: valTemp(0))
    CALL ensure_parsed(fileName)
    CALL getArrayFromYAML(fileName, dictName1, valueName, valTemp, fileStat, dictName2, dictName3)
    SELECT TYPE (valTemp)
    TYPE IS (INTEGER(2))
        val = valTemp
    END SELECT
    DEALLOCATE (valTemp)
  END FUNCTION readInt2Array

  FUNCTION readInt2Val(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    INTEGER(2) :: val

    fileStat = 1
    CALL ensure_parsed(fileName)
    CALL getVarFromYAML(fileName, dictName1, valueName, val, fileStat, dictName2, dictName3)
  END FUNCTION readInt2Val

  FUNCTION readString(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    CHARACTER(*) :: val

    fileStat = 1
    CALL ensure_parsed(fileName)
    CALL getVarFromYAML(fileName, dictName1, valueName, val, fileStat, dictName2, dictName3)
    CALL process_env_vars(val)
  END FUNCTION readString

  FUNCTION readLogicalArray(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
        CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    LOGICAL, ALLOCATABLE :: val(:)
    CLASS(*), ALLOCATABLE :: valTemp(:)

    fileStat = 1
    ALLOCATE (LOGICAL :: valTemp(0))
    CALL ensure_parsed(fileName)
    CALL getArrayFromYAML(fileName, dictName1, valueName, valTemp, fileStat, dictName2, dictName3)
    SELECT TYPE (valTemp)
    TYPE IS (LOGICAL)
        val = valTemp
    END SELECT
    DEALLOCATE (valTemp)
  END FUNCTION readLogicalArray

  FUNCTION readLogical(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
      CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3

    INTEGER :: fileStat
    LOGICAL :: val

    fileStat = 1
    CALL ensure_parsed(fileName)
    CALL getVarFromYAML(fileName, dictName1, valueName, val, fileStat, dictName2, dictName3)
  END FUNCTION readLogical


  FUNCTION readStringArray(fileName, dictName1, valueName, val, dictName2, dictName3) RESULT(fileStat)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3
    INTEGER :: fileStat
    CHARACTER(len=*), ALLOCATABLE :: val(:)
    CHARACTER(len=error_length) :: error
    CLASS(type_dictionary), POINTER :: dict
    CLASS(type_list), POINTER :: list
    CLASS(type_list_item), POINTER :: item
    TYPE(type_error), POINTER :: io_err
    INTEGER :: nArray = 0
    
    fileStat = 1
    CALL ensure_parsed(fileName)

    SELECT TYPE (root)
    CLASS is (type_dictionary)
      !> STRING from LIST: IO: ModelFileName
      dict => root%get_dictionary(TRIM(dictName1), required=.TRUE., error=io_err)
      IF(PRESENT(dictName2)) THEN
        IF (io_err%message /= '') THEN
          PRINT *, '!!->>Waring: key "', TRIM(valueName), '" reads fail. Value not change.'; 
          return
        ENDIF
        dict => dict%get_dictionary(TRIM(dictName2), required=.TRUE., error=io_err)
      END IF
      IF(PRESENT(dictName3)) THEN
        IF (io_err%message /= '') THEN
          PRINT *, '!!->>Waring: key "', TRIM(valueName), '" reads fail. Value not change.'; 
          return
        END IF
        dict => dict%get_dictionary(TRIM(dictName3), required=.TRUE., error=io_err)
      END IF
      
      IF (io_err%message /= '') THEN
        PRINT *, '!!->>Waring: key "', TRIM(valueName), '" reads fail. Value not change.'; 
        return
      ENDIF
      list => dict%get_list(TRIM(valueName), required=.TRUE., error=io_err)
      IF (.NOT. ASSOCIATED(io_err)) THEN
        ! 获取 list 的大小
        item => list%first
        nArray = 0
        DO WHILE (ASSOCIATED(item))
          item => item%next
          nArray = nArray + 1
        END DO
        IF(ALLOCATED(val)) DEALLOCATE(val)
        ALLOCATE (val(nArray))


        item => list%first
        nArray = 0
        DO WHILE (ASSOCIATED(item))
          nArray = nArray + 1
          ! 扩容 val 的大小
          SELECT TYPE (element => item%node)
          CLASS is (type_scalar)
            item => item%next
              val(nArray) = TRIM(element%string)
          END SELECT
        END DO
        fileStat = 0
      ELSE
        PRINT *, '"'//TRIM(valueName)//'" is not assigned.'
      END IF
    END SELECT

  END FUNCTION readStringArray

  SUBROUTINE getArrayFromYAML(fileName, dictName1, valueName, val, fileStat, dictName2, dictName3)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3
    CLASS(*), ALLOCATABLE :: val(:)
    CHARACTER(len=error_length) :: error
    INTEGER, INTENT(inout) :: fileStat

    CLASS(type_dictionary), POINTER :: dict
    CLASS(type_list), POINTER :: list
    CLASS(type_list_item), POINTER :: item
    TYPE(type_error), POINTER :: io_err
    CLASS(*), ALLOCATABLE :: tmp(:)
    INTEGER :: nArray = 0

    SELECT TYPE (root)
    CLASS is (type_dictionary)
      !> STRING from LIST: IO: ModelFileName
      dict => root%get_dictionary(TRIM(dictName1), required=.TRUE., error=io_err)
      IF(PRESENT(dictName2)) THEN
        IF (io_err%message /= '') THEN
          PRINT *, '!!->>Waring: key "', TRIM(valueName), '" reads fail. Value not change.'; 
          return
        ENDIF
        dict => dict%get_dictionary(TRIM(dictName2), required=.TRUE., error=io_err)
      END IF
      IF(PRESENT(dictName3)) THEN
        IF (io_err%message /= '') THEN
          PRINT *, '!!->>Waring: key "', TRIM(valueName), '" reads fail. Value not change.'; 
          return
        ENDIF
        dict => dict%get_dictionary(TRIM(dictName3), required=.TRUE., error=io_err)
      END IF
      
      IF (io_err%message /= '') THEN
        PRINT *, '!!->>Waring: key "', TRIM(valueName), '" reads fail. Value not change.'; 
        return
      ENDIF
      list => dict%get_list(TRIM(valueName), required=.TRUE., error=io_err)
      IF (.NOT. ASSOCIATED(io_err)) THEN
        ! 获取 list 的大小
        item => list%first
        nArray = 0
        DO WHILE (ASSOCIATED(item))
          item => item%next
          nArray = nArray + 1
        END DO
        
        ALLOCATE (tmp(nArray), MOLD=val)
        DEALLOCATE (val)
        ALLOCATE (val(nArray), MOLD=tmp)
        DEALLOCATE (tmp)

        item => list%first
        nArray = 0
        DO WHILE (ASSOCIATED(item))
          nArray = nArray + 1
          ! 扩容 val 的大小
          SELECT TYPE (element => item%node)
          CLASS is (type_scalar)
            item => item%next
            SELECT TYPE (val)
            TYPE is (REAL(16))
              val(nArray) = REAL(element%to_real(0.0D0), 16)
            TYPE is (REAL(8))
              val(nArray) = REAL(element%to_real(0.0D0), 8)
            TYPE is (REAL(4))
              val(nArray) = REAL(element%to_real(0.0D0), 4)
            TYPE is (INTEGER(8))
              val(nArray) = INT8(element%to_integer(0))
            TYPE is (INTEGER(4))
              val(nArray) = INT(element%to_integer(0), 4)
            TYPE is (INTEGER(2))
              val(nArray) = INT(element%to_integer(0), 2)
            TYPE is (CHARACTER(len=*))
              val(nArray) = TRIM(element%string)
            TYPE is (LOGICAL)
              val(nArray) = element%to_logical(.TRUE.)
            END SELECT
          END SELECT
        END DO

        fileStat = 0
      ELSE
        PRINT *, '"'//TRIM(valueName)//'" is not assigned.'
      END IF
    END SELECT

  END SUBROUTINE getArrayFromYAML

  SUBROUTINE getVarFromYAML(fileName, dictName1, valueName, VALUE, fileStat, dictName2, dictName3)
    CHARACTER(len=*), INTENT(IN) :: fileName, dictName1, valueName
    CHARACTER(len=*), INTENT(IN), OPTIONAL :: dictName2, dictName3
    CLASS(*) :: VALUE

    CHARACTER(len=error_length) :: error
    INTEGER, INTENT(inout) :: fileStat

    CLASS(type_dictionary), POINTER :: dict
    TYPE(type_error), POINTER :: io_err

    SELECT TYPE (root)
    CLASS is (type_dictionary)
      dict => root%get_dictionary(TRIM(dictName1), required=.TRUE., error=io_err)
      IF(PRESENT(dictName2)) THEN
        IF (io_err%message /= '') THEN
          PRINT *, '!!->>Waring: key "', TRIM(valueName), '" reads fail. Value not change.'; 
          return
        ENDIF
        dict => dict%get_dictionary(TRIM(dictName2), required=.TRUE., error=io_err)
      END IF
      IF(PRESENT(dictName3)) THEN
        IF (io_err%message /= '') THEN
          PRINT *, '!!->>Waring: key "', TRIM(valueName), '" reads fail. Value not change.'; 
          return
        ENDIF
        dict => dict%get_dictionary(TRIM(dictName3), required=.TRUE., error=io_err)
      END IF

      IF (io_err%message /= '') THEN
        PRINT *, '!!->>Waring: key "', TRIM(valueName), '" reads fail. Value not change.'
      ELSE
        SELECT TYPE (VALUE)
        TYPE IS (REAL(16))
          VALUE = dict%get_real(TRIM(valueName), error=io_err)
        TYPE IS (REAL(8))
          VALUE = dict%get_real(TRIM(valueName), error=io_err)
        TYPE IS (REAL(4))
          VALUE = dict%get_real(TRIM(valueName), error=io_err)
        TYPE IS (INTEGER(8))
          VALUE = dict%get_integer(TRIM(valueName), error=io_err)
        TYPE IS (INTEGER(4))
          VALUE = dict%get_integer(TRIM(valueName), error=io_err)
        TYPE IS (INTEGER(2))
          VALUE = dict%get_integer(TRIM(valueName), error=io_err)
        TYPE IS (CHARACTER(len=*))
          VALUE = dict%get_string(TRIM(valueName), error=io_err)
        TYPE IS (LOGICAL)
          VALUE = dict%get_logical(TRIM(valueName), error=io_err)
        END SELECT
        IF (.NOT. ASSOCIATED(io_err)) THEN
          fileStat = 0
        ELSE
          IF(io_err%message == '') fileStat = 0
        ENDIF
      END IF
    END SELECT
  END SUBROUTINE getVarFromYAML

END MODULE YAMLRead_m
