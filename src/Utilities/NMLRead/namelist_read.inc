    INTEGER            :: open_status, read_status
    character*1024     :: fileLine
    INTEGER            :: varNamePos, varBeginPos, varEndPos
    LOGICAL            :: hasPara = .false.
    OPEN (UNIT=1000, FILE=trim(nmlst_fn), ACTION="READ", IOSTAT=open_status)

    IF (open_status == 0) THEN
      DO
        read (1000, "(A)", IOSTAT=read_status) fileLine
        fileLine = AdjustL(trim(fileLine))
        varNamePos = INDEX(fileLine, varName)

        if (varNamePos == 1) then
          varBeginPos = INDEX(fileLine, "=")
          IF (TRIM(varName) .eq. TRIM(fileLine(1:varBeginPos-1))) then
            varEndPos = LEN(trim(fileLine))

            IF ((fileLine(varEndPos:varEndPos) .eq. ',') .or. (fileLine(varEndPos:varEndPos) .eq. '/')) varEndPos = varEndPos - 1
            rdVarStr = AdjustL(fileLine(varBeginPos + 1:varEndPos))
            hasPara = .true.
            EXIT
          end if
        end if

        IF (read_status /= 0) EXIT
      END DO
    ELSE
      PRINT*, 'Failed to open the NML config file. FileName is: ', nmlst_fn
    END IF

    close (UNIT=1000)

    IF(.not. hasPara) PRINT*, 'Cannot find ','varName',' !'
