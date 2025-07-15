MODULE PostProcTools_m
  USE State_m, ONLY: State_t
  USE kinds_m, ONLY: i_kind
  USE conversions_m

CONTAINS

  SUBROUTINE add_uv_to_state(X)
    TYPE(state_t) :: X
    INTEGER(i_kind) :: i, j, k

    IF (X%getVarIdx('winds') == 0) RETURN
    IF (X%getVarIdx('windd') == 0) RETURN

    CALL X%addVar('uwnd')
    CALL X%addVar('vwnd')

    DO k = 1, X%sg%vLevel
    DO i = 1, X%sg%num_cell
      DO j = 1, X%sg%tSlots
        CALL wind_to_uv(X%Fields(X%getVarIdx('winds'))%DATA(k, i, j), X%Fields(X%getVarIdx('windd'))%DATA(k, i, j), &
                        X%Fields(X%getVarIdx('uwnd'))%DATA(k, i, j), X%Fields(X%getVarIdx('vwnd'))%DATA(k, i, j))
      END DO
    END DO
    END DO

  END SUBROUTINE

  SUBROUTINE add_sd_to_state(X)
    TYPE(state_t) :: X
    INTEGER(i_kind) :: i, j, k

    IF (X%getVarIdx('vwnd') /= 0 .AND. X%getVarIdx('uwnd') /= 0) THEN
      CALL X%addVar('winds')
      CALL X%addVar('windd')

      DO k = 1, X%sg%vLevel
        DO i = 1, X%sg%num_cell
          DO j = 1, X%sg%tSlots
            CALL uv_to_wind(X%Fields(X%getVarIdx('uwnd'))%DATA(k, i, j), X%Fields(X%getVarIdx('vwnd'))%DATA(k, i, j), &
                            X%Fields(X%getVarIdx('winds'))%DATA(k, i, j), X%Fields(X%getVarIdx('windd'))%DATA(k, i, j))
          END DO
        END DO
      END DO
    END IF

    IF (X%getVarIdx('vwnd2min') /= 0 .AND. X%getVarIdx('uwnd2min') /= 0) THEN
      CALL X%addVar('winds2min')
      CALL X%addVar('windd2min')

      DO k = 1, X%sg%vLevel
        DO i = 1, X%sg%num_cell
          DO j = 1, X%sg%tSlots
            CALL uv_to_wind(X%Fields(X%getVarIdx('uwnd2min'))%DATA(k, i, j), X%Fields(X%getVarIdx('vwnd2min'))%DATA(k, i, j), &
                            X%Fields(X%getVarIdx('winds2min'))%DATA(k, i, j), X%Fields(X%getVarIdx('windd2min'))%DATA(k, i, j))
          END DO
        END DO
      END DO
    END IF

    IF (X%getVarIdx('vwnd10min') /= 0 .AND. X%getVarIdx('uwnd10min') /= 0) THEN
      CALL X%addVar('winds10min')
      CALL X%addVar('windd10min')

      DO k = 1, X%sg%vLevel
        DO i = 1, X%sg%num_cell
          DO j = 1, X%sg%tSlots
            CALL uv_to_wind(X%Fields(X%getVarIdx('uwnd10min'))%DATA(k, i, j), X%Fields(X%getVarIdx('vwnd10min'))%DATA(k, i, j), &
                            X%Fields(X%getVarIdx('winds10min'))%DATA(k, i, j), X%Fields(X%getVarIdx('windd10min'))%DATA(k, i, j))
          END DO
        END DO
      END DO
    END IF

  END SUBROUTINE

  ! This is a temporary subroutine for Shenzhen
  SUBROUTINE add_10min_uv_to_state(X)
    TYPE(state_t) :: X
    INTEGER(i_kind) :: i, j, k

    IF (X%getVarIdx('uwnd') == 0) RETURN
    IF (X%getVarIdx('vwnd') == 0) RETURN

    CALL X%addVar('uwnd10min')
    CALL X%addVar('vwnd10min')

    X%Fields(X%getVarIdx('uwnd10min'))%DATA = X%Fields(X%getVarIdx('uwnd'))%DATA
    X%Fields(X%getVarIdx('vwnd10min'))%DATA = X%Fields(X%getVarIdx('vwnd'))%DATA

    X%Fields(X%getVarIdx('vwnd10min'))%DATA(:, :, X%sg%tSlots) = (X%Fields(X%getVarIdx('vwnd'))%DATA(:, :, X%sg%tSlots - 1) + &
                                                                  X%Fields(X%getVarIdx('vwnd'))%DATA(:, :, X%sg%tSlots - 2) + &
                                                                  X%Fields(X%getVarIdx('vwnd'))%DATA(:, :, X%sg%tSlots)) / 3.0D0

    X%Fields(X%getVarIdx('uwnd10min'))%DATA(:, :, X%sg%tSlots) = (X%Fields(X%getVarIdx('uwnd'))%DATA(:, :, X%sg%tSlots - 1) + &
                                                                  X%Fields(X%getVarIdx('uwnd'))%DATA(:, :, X%sg%tSlots - 2) + &
                                                                  X%Fields(X%getVarIdx('uwnd'))%DATA(:, :, X%sg%tSlots)) / 3.0D0

  END SUBROUTINE

  ! SUBROUTINE add_wind_pressure(X)
  !   USE parameters_m
  !   TYPE(state_t) :: X
  !   INTEGER(i_kind) :: i, j, k

  !   IF (X%getVarIdx('uwnd') == 0) return
  !   IF (X%getVarIdx('vwnd') == 0) return
  !   IF (X%getVarIdx('temp') == 0) return
  !   IF (X%getVarIdx('pres') == 0) return
  !   IF (X%getVarIdx('qvapor') == 0) return

  !   CALL X%addVar('windpres')

  !   ASSOCIATE (uwnd => X%Fields(X%getVarIdx('uwnd'))%data, vwnd => X%Fields(X%getVarIdx('vwnd'))%data, &
  !              temp => X%Fields(X%getVarIdx('temp'))%data, pres => X%Fields(X%getVarIdx('pres'))%data, &
  !              windpres => X%Fields(X%getVarIdx('windpres'))%data, qvapor => X%Fields(X%getVarIdx('qvapor'))%data)
  !   DO k = 1, X%sg%vLevel
  !     DO i = 1, X%sg%num_cell
  !       DO j = 1, X%sg%tSlots
  !         windpres(k, i, j) = (uwnd(k, i, j)**2 + vwnd(k, i, j)**2)*0.5*pres(k, i, j)/(temp(k, i, j)*(1 + 0.608*qvapor(k, i, j))) &
  !                             /dry_air_gas_const
  !       END DO
  !     END DO
  !   END DO
  !   END ASSOCIATE

  ! END SUBROUTINE

END MODULE PostProcTools_m
