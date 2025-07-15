!----------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Centers
!                     for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation (SIMI)
! AUTOHR(S)         : Jiongming Pang
! VERSION           : V 0.0
! HISTORY           : V 0.0 created WRFIO_t as a constructor function for reading wrfout
!
! V 0.0 Created by Jiongming Pang (pang.j.m@hotmail.com), 2021-10-31, @GBA-MWF, Shenzhen
!----------------------------------------------------------------------------------------

MODULE WRFIO_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE ncReadVar_m
  USE AdvanceTime_m
  USE Widgets_m
  USE parameters_m, ONLY: degree2radian

  TYPE WRFIO_t

    INTEGER(i_kind) :: num_files   ! the numbers of background files

    INTEGER(i_kind), DIMENSION(:), ALLOCATABLE :: time_unix
    INTEGER(i_kind), DIMENSION(:, :), ALLOCATABLE :: time_gmt

    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: lon2D, lat2D
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: lon2D_U, lat2D_U
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: lon2D_V, lat2D_V

    REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE :: t2m  !T2      ! X x Y x T
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: temp !T       ! X x Y x Z x T

    REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE :: u10m   !U10_a     ! X x Y x T
    REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE :: v10m   !V10_a     ! X x Y x T

    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: uwnd_s !U_c       ! X x Y x Z x T
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: vwnd_s !V_c       ! X x Y x Z x T
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: wwnd_s !W_c       ! X x Y x Z x T
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: uwnd   !U_a       ! X x Y x Z x T
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: vwnd   !V_a       ! X x Y x Z x T
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: wwnd   !W_a       ! X x Y x Z x T

    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: P       ! X x Y x Z x T
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: PB      ! X x Y x Z x T
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: PH      ! X x Y x Z x T
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: PHB     ! X x Y x Z x T
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: pres !Pres    ! X x Y x Z x T

    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: height_s !HEIGHT  ! X x Y x Z x T
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: height   !HEIGHT  ! X x Y x Z x T
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: topo   !HGT     ! X x Y
    REAL(r_kind), DIMENSION(:, :), ALLOCATABLE :: landmask !LANDMASK     ! X x Y

    REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE :: psfc !PSFC    ! X x Y x T
    REAL(r_kind), DIMENSION(:, :, :), ALLOCATABLE :: q2m  !Q2      ! X x Y x T
    REAL(r_kind), DIMENSION(:, :, :, :), ALLOCATABLE :: qvapor !QVAPOR  ! X x Y x Z x T

    INTEGER(i_kind) :: nx, ny, nz, nt

  END TYPE WRFIO_t

  INTERFACE WRFIO_t
    PROCEDURE :: constructor
  END INTERFACE

CONTAINS

  ! FUNCTION constructor(configFile) result(this)
  FUNCTION constructor(wrfFileNames, varNames, modeFlag) RESULT(this)
    TYPE(WRFIO_t) :: this

    TYPE(ncDimInfo_t), ALLOCATABLE  ::  dimInfo
    TYPE(ncVarInfo_t), ALLOCATABLE  ::  varInfo
    TYPE(ncData3D_t), ALLOCATABLE  ::  vardata3D
    TYPE(ncData4D_t), ALLOCATABLE  ::  vardata4D
    TYPE(ncStr_t), ALLOCATABLE  ::  varstr

    CHARACTER(LEN=1024), INTENT(IN) ::  wrfFileNames(:)
    CHARACTER(LEN=20), INTENT(IN) ::  varNames(:)
    CHARACTER(LEN=20), INTENT(IN), OPTIONAL ::  modeFlag

    CHARACTER(LEN=20), ALLOCATABLE  ::  dimName(:)

    CHARACTER(LEN=1024) :: wrfFileDir
    CHARACTER(LEN=1024) :: inFileName
    CHARACTER(LEN=10)   :: vartmp
    INTEGER(i_kind)     :: t, i

    CHARACTER(LEN=20), ALLOCATABLE  ::  usrNames(:)
    INTEGER(i_kind), ALLOCATABLE  ::  varIndx(:)

    ! GET the info of model domain
    this%num_files = UBOUND(wrfFileNames, 1)
    inFileName = TRIM(wrfFileNames(1))
    ALLOCATE (dimInfo, dimName(3))
    dimName(1) = "west_east"
    dimName(2) = "south_north"
    ! dimName(3) = "bottom_top_stag"
    dimName(3) = "bottom_top"
    ! PRINT *, "dimName:", dimName
    CALL dimInfo%ncFileDimsLen(TRIM(inFileName), dimName)
    this%nx = dimInfo%file_dimLen(1)
    this%ny = dimInfo%file_dimLen(2)
    this%nz = dimInfo%file_dimLen(3)
    PRINT *, "nx, ny, nz:", this%nx, this%ny, this%nz
    DEALLOCATE (dimInfo, dimName)

    !***** GET THE NECCESSARY VARIABLES, INCLUDING XLAT, XLONG, HGT, and LANDMASK. *****!

    ! XLONG
    ALLOCATE (this%lon2D(this%nx, this%ny))
    ALLOCATE (vardata3D)
    CALL vardata3D%ncGetData3d(TRIM(inFileName), TRIM("XLONG"))
    this%lon2D = vardata3D%ncVar(:, :, 1) * degree2radian
    PRINT *, 'min, max, medium of XLONG are: ', &
      MINVAL(this%lon2D / degree2radian), MAXVAL(this%lon2D / degree2radian), &
      this%lon2D(this%nx / 2, this%ny / 2) / degree2radian
    DEALLOCATE (vardata3D)
    PRINT *, "DONE READING LON..."

    ! XLAT
    ALLOCATE (this%lat2D(this%nx, this%ny))
    ALLOCATE (vardata3D)
    CALL vardata3D%ncGetData3d(TRIM(inFileName), TRIM("XLAT"))
    this%lat2D = vardata3D%ncVar(:, :, 1) * degree2radian
    PRINT *, 'min, max, medium of XLAT are: ', &
      MINVAL(this%lat2D) / degree2radian, MAXVAL(this%lat2D) / degree2radian, &
      this%lat2D(this%nx / 2, this%ny / 2) / degree2radian
    DEALLOCATE (vardata3D)
    PRINT *, "DONE READING LAT..."

    ! HGT
    ALLOCATE (this%topo(this%nx, this%ny))
    vartmp = "HGT"
    ALLOCATE (vardata3D)
    CALL vardata3D%ncGetData3d(TRIM(inFileName), TRIM(vartmp))
    this%topo(:, :) = vardata3D%ncVar(:, :, 1)
    PRINT *, 'max and min of ', TRIM(vartmp), ': ', MAXVAL(vardata3D%ncVar), MINVAL(vardata3D%ncVar)
    DEALLOCATE (vardata3D)

    ! LANDMASK
    ! ALLOCATE(this%landmask(this%nx, this%ny))
    ! vartmp = "LANDMASK"
    ! ALLOCATE(vardata3D)
    ! CALL vardata3D%ncGetData3d(TRIM(inFileName), TRIM(vartmp))
    ! this%landmask(:,:) = vardata3D%ncVar(:,:,1)
    ! PRINT *, 'max and min of ', TRIM(vartmp) ,': ', MAXVAL(vardata3D%ncVar), MINVAL(vardata3D%ncVar)
    ! DEALLOCATE(vardata3D)

    !! The staggered latitude and longitude were not neccesary in MOTOR-DA,
    !! Users can uncomment this part if you need.

    ! !!! BEGIN READING LONGITUDE_U AND LATITUDE_U !!!
    ! ! READING longitude_U and latitude_U from wrfout file
    ! ALLOCATE(this%lon2D_U(this%nx+1, this%ny), this%lat2D_U(this%nx+1, this%ny))

    ! ! get longitude_U
    ! ALLOCATE(vardata3D)
    ! CALL vardata3D%ncGetData3d(TRIM(inFileName), TRIM("XLONG_U"))
    ! this%lon2D_U = vardata3D%ncVar(:,:,1) * degree2radian
    ! DEALLOCATE(vardata3D)
    ! PRINT *, "DONE READING LON_U..."

    ! ! get latitude_U
    ! ALLOCATE(vardata3D)
    ! CALL vardata3D%ncGetData3d(TRIM(inFileName), TRIM("XLAT_U"))
    ! this%lat2D_U = vardata3D%ncVar(:,:,1) * degree2radian
    ! DEALLOCATE(vardata3D)
    ! PRINT *, "DONE READING LAT_U..."
    ! !!! END READING LONGITUDE_U AND LATITUDE_U !!!

    ! !!! BEGIN READING LONGITUDE_V AND LATITUDE_V !!!
    ! ! READING longitude_V and latitude_V from wrfout file
    ! ALLOCATE(this%lon2D_V(this%nx, this%ny+1), this%lat2D_V(this%nx, this%ny+1))

    ! ! get longitude_V
    ! ALLOCATE(vardata3D)
    ! CALL vardata3D%ncGetData3d(TRIM(inFileName), TRIM("XLONG_V"))
    ! this%lon2D_U = vardata3D%ncVar(:,:,1) * degree2radian
    ! DEALLOCATE(vardata3D)
    ! PRINT *, "DONE READING LON_V..."

    ! ! get latitude_V
    ! ALLOCATE(vardata3D)
    ! CALL vardata3D%ncGetData3d(TRIM(inFileName), TRIM("XLAT_V"))
    ! this%lat2D_V = vardata3D%ncVar(:,:,1) * degree2radian
    ! DEALLOCATE(vardata3D)
    ! PRINT *, "DONE READING LAT_V..."
    ! !!! END READING LONGITUDE_V AND LATITUDE_V !!!

    !!! BEGIN READING VARIABLES !!!
    ! allocate var
    this%nt = this%num_files
    PRINT *, "nx, ny, nz, nt, num_files:", this%nx, this%ny, this%nz, this%nt, this%num_files

    ! for t2m
    varIndx = getStrIndexs(varNames, (/"t2m"/))
    IF (MAXVAL(varIndx) > 0) ALLOCATE (this%t2m(this%nx, this%ny, this%nt))
    DEALLOCATE (varIndx)

    ! for temp
    varIndx = getStrIndexs(varNames, (/"temp"/))
    IF (MAXVAL(varIndx) > 0) THEN
      ALLOCATE (this%temp(this%nx, this%ny, this%nz, this%nt))
      IF (.NOT. ALLOCATED(this%P)) ALLOCATE (this%P(this%nx, this%ny, this%nz, this%nt))
      IF (.NOT. ALLOCATED(this%PB)) ALLOCATE (this%PB(this%nx, this%ny, this%nz, this%nt))
    END IF
    DEALLOCATE (varIndx)

    ! for u10m and v10m
    ALLOCATE (usrNames(2))
    usrNames(1) = "u10m"; usrNames(2) = "v10m"
    varIndx = getStrIndexs(varNames, usrNames)
    IF (MAXVAL(varIndx) > 0) ALLOCATE (this%u10m(this%nx, this%ny, this%nt), this%v10m(this%nx, this%ny, this%nt))
    DEALLOCATE (usrNames, varIndx)

    ! for uwnd, vwnd
    ALLOCATE (usrNames(4))
    usrNames(1) = "uwnd_s"; usrNames(2) = "vwnd_s"; usrNames(3) = "uwnd"; usrNames(4) = "vwnd"
    varIndx = getStrIndexs(varNames, usrNames)
    IF (MAXVAL(varIndx) > 0) &
      ALLOCATE (this%uwnd_s(this%nx + 1, this%ny, this%nz, this%nt), &
                this%vwnd_s(this%nx, this%ny + 1, this%nz, this%nt), &
                this%uwnd(this%nx, this%ny, this%nz, this%nt), &
                this%vwnd(this%nx, this%ny, this%nz, this%nt))
    DEALLOCATE (usrNames, varIndx)

    varIndx = getStrIndexs(varNames, (/"wwnd_s"/))
    IF (MAXVAL(varIndx) > 0) ALLOCATE (this%wwnd_s(this%nx, this%ny, this%nz + 1, this%nt))
    DEALLOCATE (varIndx)

    ! for pres
    ALLOCATE (usrNames(3))
    usrNames(1) = "pres"; usrNames(2) = "P"; usrNames(3) = "PB"
    varIndx = getStrIndexs(varNames, usrNames)
    IF (MAXVAL(varIndx) > 0) THEN
      ALLOCATE (this%pres(this%nx, this%ny, this%nz, this%nt))
      IF (.NOT. ALLOCATED(this%P)) ALLOCATE (this%P(this%nx, this%ny, this%nz, this%nt))
      IF (.NOT. ALLOCATED(this%PB)) ALLOCATE (this%PB(this%nx, this%ny, this%nz, this%nt))
    END IF
    DEALLOCATE (usrNames, varIndx)

    ! for psfc
    varIndx = getStrIndexs(varNames, (/"psfc"/))
    IF (MAXVAL(varIndx) > 0) ALLOCATE (this%psfc(this%nx, this%ny, this%nt))
    DEALLOCATE (varIndx)

    ! for q2m
    varIndx = getStrIndexs(varNames, (/"q2m"/))
    IF (MAXVAL(varIndx) > 0) ALLOCATE (this%q2m(this%nx, this%ny, this%nt))
    DEALLOCATE (varIndx)

    ! for qvapor
    varIndx = getStrIndexs(varNames, (/"qvapor"/))
    IF (MAXVAL(varIndx) > 0) ALLOCATE (this%qvapor(this%nx, this%ny, this%nz, this%nt))
    DEALLOCATE (varIndx)

    ! for height
    ALLOCATE (usrNames(4))
    usrNames(1) = "height"; usrNames(2) = "PH"; usrNames(3) = "PHB"; usrNames(4) = "height_s"
    varIndx = getStrIndexs(varNames, usrNames)
    IF (MAXVAL(varIndx) > 0) THEN
      ALLOCATE (this%height(this%nx, this%ny, this%nz, this%nt))
      ALLOCATE (this%height_s(this%nx, this%ny, this%nz + 1, this%nt))
      IF (.NOT. ALLOCATED(this%PH)) ALLOCATE (this%PH(this%nx, this%ny, this%nz + 1, this%nt))
      IF (.NOT. ALLOCATED(this%PHB)) ALLOCATE (this%PHB(this%nx, this%ny, this%nz + 1, this%nt))
    END IF
    DEALLOCATE (usrNames, varIndx)

    ! for time_lable
    ALLOCATE (this%time_gmt(this%num_files, 6), this%time_unix(this%num_files))

    ! PRINT*, 'Before T2'
    ! read variables from wrfout files
    DO t = 1, this%nt
      inFileName = TRIM(wrfFileNames(t))

      ! T2
      IF (ALLOCATED(this%t2m)) THEN
        vartmp = "T2"
        ALLOCATE (vardata3D)
        CALL vardata3D%ncGetData3d(TRIM(inFileName), TRIM(vartmp))
        this%t2m(:, :, t) = vardata3D%ncVar(:, :, 1)
        IF (t .EQ. 1) THEN
          PRINT *, 'max and min of ', TRIM(vartmp), ': ', MAXVAL(vardata3D%ncVar), MINVAL(vardata3D%ncVar)
        END IF
        DEALLOCATE (vardata3D)
      END IF

      ! PRINT*, 'Before T'
      ! T
      IF (ALLOCATED(this%temp)) THEN
        ! vartmp = "t"
        vartmp = "T"
        ALLOCATE (vardata4D)
        CALL vardata4D%ncGetData4d(TRIM(inFileName), TRIM(vartmp))
        this%temp(:, :, :, t) = vardata4D%ncVar(:, :, :, 1)
        ! IF (t .eq. 1) THEN
        !   PRINT *, 'max and min of ', TRIM(vartmp) ,': ', MAXVAL(vardata4D%ncVar), MINVAL(vardata4D%ncVar)
        ! END IF
        DEALLOCATE (vardata4D)
      END IF

      ! PRINT*, 'Before U10'
      ! U10
      IF (ALLOCATED(this%u10m)) THEN
        vartmp = "U10"
        ALLOCATE (vardata3D)
        CALL vardata3D%ncGetData3d(TRIM(inFileName), TRIM(vartmp))
        this%u10m(:, :, t) = vardata3D%ncVar(:, :, 1)
        IF (t .EQ. 1) THEN
          PRINT *, 'max and min of ', TRIM(vartmp), ': ', MAXVAL(vardata3D%ncVar), MINVAL(vardata3D%ncVar)
        END IF
        DEALLOCATE (vardata3D)
      END IF

      ! U
      IF (ALLOCATED(this%uwnd_s)) THEN
        vartmp = "U"
        ALLOCATE (vardata4D)
        CALL vardata4D%ncGetData4d(TRIM(inFileName), TRIM(vartmp))
        this%uwnd_s(:, :, :, t) = vardata4D%ncVar(:, :, :, 1)
        ! IF (t .eq. 1) THEN
        !   PRINT *, 'max and min of ', TRIM(vartmp) ,': ', MAXVAL(vardata4D%ncVar), MINVAL(vardata4D%ncVar)
        ! END IF
        DEALLOCATE (vardata4D)
      END IF

      ! V10
      IF (ALLOCATED(this%v10m)) THEN
        vartmp = "V10"
        ALLOCATE (vardata3D)
        CALL vardata3D%ncGetData3d(TRIM(inFileName), TRIM(vartmp))
        this%v10m(:, :, t) = vardata3D%ncVar(:, :, 1)
        IF (t .EQ. 1) THEN
          PRINT *, 'max and min of ', TRIM(vartmp), ': ', MAXVAL(vardata3D%ncVar), MINVAL(vardata3D%ncVar)
        END IF
        DEALLOCATE (vardata3D)
      END IF

      ! V
      IF (ALLOCATED(this%vwnd_s)) THEN
        vartmp = "V"
        ALLOCATE (vardata4D)
        CALL vardata4D%ncGetData4d(TRIM(inFileName), TRIM(vartmp))
        this%vwnd_s(:, :, :, t) = vardata4D%ncVar(:, :, :, 1)
        ! IF (t .eq. 1) THEN
        !   PRINT *, 'max and min of ', TRIM(vartmp) ,': ', MAXVAL(vardata4D%ncVar), MINVAL(vardata4D%ncVar)
        ! END IF
        DEALLOCATE (vardata4D)
      END IF

      ! W
      IF (ALLOCATED(this%wwnd_s)) THEN
        vartmp = "W"
        ALLOCATE (vardata4D)
        CALL vardata4D%ncGetData4d(TRIM(inFileName), TRIM(vartmp))
        this%wwnd_s(:, :, :, t) = vardata4D%ncVar(:, :, :, 1)
        ! IF (t .eq. 1) THEN
        !   PRINT *, 'max and min of ', TRIM(vartmp) ,': ', MAXVAL(vardata4D%ncVar), MINVAL(vardata4D%ncVar)
        ! END IF
        DEALLOCATE (vardata4D)
      END IF

      ! P
      IF (ALLOCATED(this%P)) THEN
        vartmp = "P"
        ALLOCATE (vardata4D)
        CALL vardata4D%ncGetData4d(TRIM(inFileName), TRIM(vartmp))
        this%P(:, :, :, t) = vardata4D%ncVar(:, :, :, 1)
        ! IF (t .eq. 1) THEN
        !   PRINT *, 'max and min of ', TRIM(vartmp) ,': ', MAXVAL(vardata4D%ncVar), MINVAL(vardata4D%ncVar)
        ! END IF
        DEALLOCATE (vardata4D)
      END IF

      ! PB
      IF (ALLOCATED(this%PB)) THEN
        vartmp = "PB"
        ALLOCATE (vardata4D)
        CALL vardata4D%ncGetData4d(TRIM(inFileName), TRIM(vartmp))
        this%PB(:, :, :, t) = vardata4D%ncVar(:, :, :, 1)
        ! IF (t .eq. 1) THEN
        !   PRINT *, 'max and min of ', TRIM(vartmp) ,': ', MAXVAL(vardata4D%ncVar), MINVAL(vardata4D%ncVar)
        ! END IF
        DEALLOCATE (vardata4D)
      END IF

      ! PH
      IF (ALLOCATED(this%PH)) THEN
        vartmp = "PH"
        ALLOCATE (vardata4D)
        CALL vardata4D%ncGetData4d(TRIM(inFileName), TRIM(vartmp))
        this%PH(:, :, :, t) = vardata4D%ncVar(:, :, :, 1)
        ! IF (t .eq. 1) THEN
        !   PRINT *, 'max and min of ', TRIM(vartmp) ,': ', MAXVAL(vardata4D%ncVar), MINVAL(vardata4D%ncVar)
        ! END IF
        DEALLOCATE (vardata4D)
      END IF

      ! PHB
      IF (ALLOCATED(this%PHB)) THEN
        vartmp = "PHB"
        ALLOCATE (vardata4D)
        CALL vardata4D%ncGetData4d(TRIM(inFileName), TRIM(vartmp))
        this%PHB(:, :, :, t) = vardata4D%ncVar(:, :, :, 1)
        ! IF (t .eq. 1) THEN
        !   PRINT *, 'max and min of ', TRIM(vartmp) ,': ', MAXVAL(vardata4D%ncVar), MINVAL(vardata4D%ncVar)
        ! END IF
        DEALLOCATE (vardata4D)
      END IF

      ! PSFC
      IF (ALLOCATED(this%psfc)) THEN
        vartmp = "PSFC"
        ALLOCATE (vardata3D)
        CALL vardata3D%ncGetData3d(TRIM(inFileName), TRIM(vartmp))
        this%psfc(:, :, t) = vardata3D%ncVar(:, :, 1)
        IF (t .EQ. 1) THEN
          PRINT *, 'max and min of ', TRIM(vartmp), ': ', MAXVAL(vardata3D%ncVar), MINVAL(vardata3D%ncVar)
        END IF
        DEALLOCATE (vardata3D)
      END IF

      ! Q2
      IF (ALLOCATED(this%q2m)) THEN
        vartmp = "Q2"
        ALLOCATE (vardata3D)
        CALL vardata3D%ncGetData3d(TRIM(inFileName), TRIM(vartmp))
        this%q2m(:, :, t) = vardata3D%ncVar(:, :, 1)
        IF (t .EQ. 1) THEN
          PRINT *, 'max and min of ', TRIM(vartmp), ': ', MAXVAL(vardata3D%ncVar), MINVAL(vardata3D%ncVar)
        END IF
        DEALLOCATE (vardata3D)
      END IF

      ! QVAPOR
      IF (ALLOCATED(this%qvapor)) THEN
        vartmp = "QVAPOR"
        ALLOCATE (vardata4D)
        CALL vardata4D%ncGetData4d(TRIM(inFileName), TRIM(vartmp))
        this%qvapor(:, :, :, t) = vardata4D%ncVar(:, :, :, 1)
        IF (t .EQ. 1) THEN
          PRINT *, 'max and min of ', TRIM(vartmp), ': ', MAXVAL(vardata4D%ncVar), MINVAL(vardata4D%ncVar)
        END IF
        DEALLOCATE (vardata4D)
      END IF

      ! Times str
      ALLOCATE (varstr)
      CALL varstr%ncGetStr(TRIM(inFileName), "Times")
      PRINT *, 'Time Label:', varstr%ncVar(1)
      READ (varstr%ncVar(1) (1:4), "(i4)") this%time_gmt(t, 1) ! year
      READ (varstr%ncVar(1) (6:7), "(i2)") this%time_gmt(t, 2) ! mon
      READ (varstr%ncVar(1) (9:10), "(i2)") this%time_gmt(t, 3) ! day
      READ (varstr%ncVar(1) (12:13), "(i2)") this%time_gmt(t, 4) ! hor
      READ (varstr%ncVar(1) (15:16), "(i2)") this%time_gmt(t, 5) ! min
      READ (varstr%ncVar(1) (18:19), "(i2)") this%time_gmt(t, 6) ! sec
      PRINT *, "GMT Time: ", this%time_gmt(t, :)
      CALL Time_GMT_to_Unix(this%time_gmt(t, :), this%time_unix(t))
      PRINT *, "Unix Time: ", this%time_unix(t)

      DEALLOCATE (varstr)

    END DO ! end loop of this%nt

    IF (ALLOCATED(this%height_s)) THEN
      this%height_s = (this%PH + this%PHB) / 9.81D0
      ! PRINT *, 'max and min of height_s: ', MAXVAL(this%height_s), MINVAL(this%height_s)
    END IF

    IF (ALLOCATED(this%temp)) THEN
      this%temp = (this%temp + 300.0D0) * (((this%P + this%PB) / 1.0D5)**0.2854D0)
      PRINT *, 'max and min of temp: ', MAXVAL(this%temp), MINVAL(this%temp)
    END IF

    ! IF ( ALLOCATED(this%temp) ) THEN
    !   this%temp = (this%temp + 273.15D0)
    !   PRINT *, 'max and min of temp: ', MAXVAL(this%temp), MINVAL(this%temp)
    ! END IF

    IF (ALLOCATED(this%pres)) THEN
      this%pres = this%P + this%PB
      PRINT *, 'max and min of pres: ', MAXVAL(this%pres), MINVAL(this%pres)
    END IF

    IF (ALLOCATED(this%height)) THEN
      DO i = 1, this%nz
        this%height(:, :, i, :) = (this%height_s(:, :, i, :) + this%height_s(:, :, i, :)) / 2.0D0
      END DO
      PRINT *, 'max and min of height: ', MAXVAL(this%height), MINVAL(this%height)
    END IF

    IF (ALLOCATED(this%uwnd)) THEN
      DO i = 1, this%nx
        this%uwnd(i, :, :, :) = (this%uwnd_s(i, :, :, :) + this%uwnd_s(i + 1, :, :, :)) / 2.0D0
      END DO
      PRINT *, 'max and min of uwnd: ', MAXVAL(this%uwnd), MINVAL(this%uwnd)
    END IF

    IF (ALLOCATED(this%uwnd)) THEN
      DO i = 1, this%ny
        this%vwnd(:, i, :, :) = (this%vwnd_s(:, i, :, :) + this%vwnd_s(:, i + 1, :, :)) / 2.0D0
      END DO
      PRINT *, 'max and min of vwnd: ', MAXVAL(this%vwnd), MINVAL(this%vwnd)
    END IF

    PRINT *, "DONE WRFIO_t"
    !!! END READING VARIABLES !!!

  END FUNCTION constructor

END MODULE WRFIO_m
