!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Zilong Qin
! VERSION           : V 0.0
! HISTORY           :
!   Created by Zilong Qin (zilong.qin@gmail.com), 2021/1/26, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
MODULE Mock_m
  USE ObsSet_m, ONLY: ObsSet_t
  USE MPObs_m, ONLY: MPObs_t
  USE ObsField_m, ONLY: ObsField_t
  USE SingleGrid_m, ONLY: SingleGrid_t
  USE kinds_m, ONLY: i_kind, r_kind

CONTAINS

  SUBROUTINE Set_Mock_Single_Pt(configFile, Y, mpObs, VALUE, sg, varName)
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(ObsSet_t), INTENT(INOUT) :: Y
    TYPE(MPObs_t), TARGET, INTENT(IN) :: mpObs
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    REAL(r_kind) :: VALUE
    INTEGER(i_kind) :: numObs, i, j, numObsTotal = 1024, k
    CHARACTER(*) :: varName

    Y = ObsSet_t(configFile, mpObs)

    IF (ALLOCATED(Y%ObsFields)) RETURN

    ALLOCATE (Y%ObsFields(1))
    Y%ObsFields(1) = ObsField_t(configFile, mpObs)

    numObsTotal = 32
    numObs = numObsTotal / (sg%mpddInfo_sg%proc_layout(1) * sg%mpddInfo_sg%proc_layout(2))
    ALLOCATE (Y%ObsFields(1)%idx(numObs * sg%tSlots * sg%vLevel))
    ALLOCATE (Y%ObsFields(1)%values(numObs * sg%tSlots * sg%vLevel))

    PRINT *, 'sg%tSlots', sg%tSlots

    ! CALL RANDOM_SEED(PUT = 1)
    DO i = 1, numObs
    DO j = 1, sg%tSlots
    DO k = 1, sg%vLevel
      BLOCK
        REAL(r_kind) :: rnd
        CALL Y%ObsFields(1)%Set_Name(varName)
        CALL RANDOM_NUMBER(rnd)
        ! rnd = i*1.0D0/numObs

        Y%ObsFields(1)%idx(j + (i - 1) * sg%tSlots + (k - 1) * sg%tSlots * numObs)%tIdx = j
        Y%ObsFields(1)%idx(j + (i - 1) * sg%tSlots + (k - 1) * sg%tSlots * numObs)%hIdx = rnd * sg%dimCell(1) * sg%dimCell(2)

        IF (Y%ObsFields(1)%idx(j + (i - 1) * sg%tSlots + (k - 1) * sg%tSlots * numObs)%hIdx .EQ. 0) &
          Y%ObsFields(1)%idx(j + (i - 1) * sg%tSlots + (k - 1) * sg%tSlots * numObs)%hIdx = 1

        Y%ObsFields(1)%idx(j + (i - 1) * sg%tSlots + (k - 1) * sg%tSlots * numObs)%vIdx = k
        Y%ObsFields(1)%values(j + (i - 1) * sg%tSlots + (k - 1) * sg%tSlots * numObs) = VALUE
      END BLOCK
    END DO
    END DO
    END DO

  END SUBROUTINE

  SUBROUTINE Set_MP_ObsFields(configFile, Y, mpObs, VALUE, zenangle, azangle, sunzenangle, sunazangle, sg, ObsName, nchans)
    IMPLICIT NONE
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(ObsSet_t), INTENT(INOUT) :: Y
    TYPE(MPObs_t), TARGET, INTENT(IN) :: mpObs
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    CHARACTER(*), INTENT(IN) :: ObsName
    INTEGER(i_kind), INTENT(IN) :: nchans
    REAL(r_kind), INTENT(IN) :: VALUE, zenangle, azangle, sunzenangle, sunazangle
    INTEGER(i_kind) :: numChans, ichan
    INTEGER(i_kind) :: numObs, i, j, numObsTotal = 2, k

    Y = ObsSet_t(configFile, mpObs)

    IF (ALLOCATED(Y%ObsFields)) RETURN

    ALLOCATE (Y%ObsFields(nchans))

    DO ichan = 1, nchans
      Y%ObsFields(ichan) = ObsField_t(configFile, mpObs)

      CALL Y%ObsFields(ichan)%Set_Name(ObsName)

      numObsTotal = 40
      numObs = numObsTotal / (sg%mpddInfo_sg%proc_layout(1) * sg%mpddInfo_sg%proc_layout(2))
      numChans = nchans

      ! Note that sg%vLevel should be 1 for satellite radiance
      ALLOCATE (Y%ObsFields(ichan)%idx(numObs*sg%tSlots))
      ALLOCATE (Y%ObsFields(ichan)%values(numObs*sg%tSlots))
      ALLOCATE (Y%ObsFields(ichan)%errors(numObs*sg%tSlots))
      ALLOCATE (Y%obsFields(ichan)%ObsAttrSat%zenangle(numObs*sg%tSlots))
      ALLOCATE (Y%obsFields(ichan)%ObsAttrSat%azangle(numObs*sg%tSlots))
      ALLOCATE (Y%obsFields(ichan)%ObsAttrSat%sunzenangle(numObs*sg%tSlots))
      ALLOCATE (Y%obsFields(ichan)%ObsAttrSat%sunazangle(numObs*sg%tSlots))

      DO i = 1, numObs
        DO j = 1, sg%tSlots
          BLOCK
            REAL(r_kind) :: rnd
            CALL RANDOM_NUMBER(rnd)

            Y%ObsFields(ichan)%idx(j + (i - 1) * sg%tSlots)%tIdx = j
            Y%ObsFields(ichan)%idx(j + (i - 1) * sg%tSlots)%hIdx = rnd * sg%dimCell(1) * sg%dimCell(2)

              IF (Y%ObsFields(ichan)%idx(j + (i - 1)*sg%tSlots )%hIdx .eq. 0) &
                Y%ObsFields(ichan)%idx(j + (i - 1)*sg%tSlots )%hIdx = 1
              
              Y%ObsFields(ichan)%idx(j + (i - 1)*sg%tSlots )%vIdx = 1
              Y%ObsFields(ichan)%values(j + (i - 1)*sg%tSlots ) = value
              Y%ObsFields(ichan)%errors(j + (i - 1)*sg%tSlots ) = 1.0

            Y%ObsFields(ichan)%idx(j + (i - 1) * sg%tSlots)%vIdx = 1
            Y%ObsFields(ichan)%values(j + (i - 1) * sg%tSlots) = VALUE

            Y%obsFields(ichan)%ObsAttrSat%zenangle(j + (i - 1) * sg%tSlots) = zenangle
            Y%obsFields(ichan)%ObsAttrSat%azangle(j + (i - 1) * sg%tSlots) = azangle
            Y%obsFields(ichan)%ObsAttrSat%sunzenangle(j + (i - 1) * sg%tSlots) = sunzenangle
            Y%obsFields(ichan)%ObsAttrSat%sunazangle(j + (i - 1) * sg%tSlots) = sunazangle

          END BLOCK
        END DO
      END DO

    END DO
  END SUBROUTINE Set_MP_ObsFields

  SUBROUTINE Set_Mock_Single_Pt_For_ObsSet(configFile, Y, mpObs, VALUE, sg)
    CHARACTER(LEN=1024), INTENT(IN) :: configFile
    TYPE(ObsSet_t), INTENT(INOUT) :: Y
    TYPE(MPObs_t), TARGET, INTENT(IN) :: mpObs
    TYPE(SingleGrid_t), TARGET, INTENT(IN) :: sg
    REAL(r_kind) :: VALUE
    INTEGER(i_kind) :: numObs, i

    Y = ObsSet_t(configFile, mpObs)

    IF (ALLOCATED(Y%ObsFields)) RETURN

    ALLOCATE (Y%ObsFields(1))
    Y%ObsFields(1) = ObsField_t(configFile, mpObs)

    numObs = 1
    ALLOCATE (Y%ObsFields(1)%idx(numObs))
    ALLOCATE (Y%ObsFields(1)%values(numObs))

    CALL RANDOM_SEED()
    DO i = 1, numObs
      BLOCK
        REAL(r_kind) :: rnd
        CALL Y%ObsFields(1)%Set_Name("ua")
        CALL RANDOM_NUMBER(rnd)
        Y%ObsFields(1)%idx(i)%tIdx = 1
        Y%ObsFields(1)%idx(i)%hIdx = rnd * sg%dimCell(1) * sg%dimCell(2)
        Y%ObsFields(1)%idx(i)%vIdx = 5
        Y%ObsFields(1)%values(i) = VALUE
      END BLOCK
    END DO
    ! Y%ObsFields(1)%idx(2)%hIdx = 32768 + 10000
    ! Y%ObsFields(1)%idx(2)%vIdx = 5
    ! Y%ObsFields(1)%values(2) = value

  END SUBROUTINE

END MODULE Mock_m
