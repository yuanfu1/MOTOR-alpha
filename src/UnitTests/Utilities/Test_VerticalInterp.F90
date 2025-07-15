!!--------------------------------------------------------------------------------------------------
! PROJECT           : MOTOR-DA
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
!                     Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2025/7/3, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!> @brief
!!
!! @author Yuanfu Xie
!! @copyright (C) 2020 GBA-MWF, All rights reserved.
PROGRAM Test_VerticalInterp
  USE kinds_m, ONLY: r_kind, i_kind
  USE verticalInterp_m, ONLY: verticalInterp_t
  USE YAMLRead_m, ONLY: yaml_get_var
  USE mpddGlob_m, ONLY: mpddGlob_t

  IMPLICIT NONE

  TYPE(mpddGlob_t), TARGET :: mpddGlob
  TYPE(verticalInterp_t) :: vertInterp

  CHARACTER(LEN=1024) :: configFile
  INTEGER(i_kind) :: istatus, logOpt, nin, nout, ng, lvl = 8
  INTEGER(i_kind) :: i, l, ierr
  REAL(r_kind) :: top, z, dz
  REAL(r_kind), PARAMETER :: pi = 3.141592653589793238D0
  REAL(r_kind), PARAMETER :: off = 0.0D0 ! 100.0_r_kind
  REAL(r_kind), ALLOCATABLE :: fin(:), vin(:), fout(:), vout(:), ana(:), err(:)

  ! Get the configFile
  CALL GET_ENVIRONMENT_VARIABLE("STATIC_DIR", configFile)
  configFile = TRIM(configFile)//"/testVerticalInterp.yaml"

  ! Initialize the mpdd
  CALL mpddGlob%initialize()

  ! Initialize the vertical interpolation object
  ! CALL vertInterp%initialize(mpddGlob)

  istatus = 0 ! default pass number, yes
  ! Read in the configuration file
  ! ifile = yaml_get_var(configFile, 'IO', 'output_dir', vertInterp%outputDir)

  top = 30000.0D0

  ! Test case 1:
  nin = 5
  nout = 30
  logOpt = 1 ! flag for log interpolation

  ! Analytical vertical profile for comparison:
  ALLOCATE(ana(nout), fout(nout), vout(nout))
  DO i=1,nout
    z = REAL(i-1, r_kind) * top / REAL(nout-1, r_kind)
    vout(i) = z
    ana(i) = DSIN(pi * vout(i) / top)  ! Example vertical profile function
  END DO

  ! Multigrid vertical interpolation test
  ALLOCATE(err(lvl))
  err = 0.0_r_kind
  DO l=1,lvl

    nin = 2 * (nin-1) + 1

    ! vertical grid points:
    ng = 3 * (nin - 1) + 1

    PRINT*,'Test level', l, 'with vertical grid points:', ng, ' and output points:', nout

    ALLOCATE(fin(ng), vin(ng))

    ! 1 dimensional input and output arrays:
    ! Break the vertical into 2 sections, 1: from bottom to the middle, 2: from middle to the top
    ! The bottom section uses twice grid points as the top section
    dz = top / REAL(4*(nin-1)+1, r_kind)
    DO i=1, ng
      IF (i .LE. 2*(nin-1)+1) THEN
        z = REAL(i-1, r_kind) * dz
      ELSE
        z = 0.5D0 * top + REAL(i- 2*(nin-1) - 1) * 2.0D0 * dz
      END IF
      vin(i) = off + z
      fin(i) = DSIN(pi * vin(i) / top)   ! Example vertical profile function
    END DO

    ! Perform vertical interpolation test
    PRINT*,'Interpolating from', ng, 'to', nout, 'vertical levels at level', l, fin(1), vin(1)/top
    ! CALL vertInterp%verticalInterp(fin, vin, fout, vout, istatus, logOpt)
    CALL vertInterp%verticalInterp(fin, vin, fout, vout, istatus)
    IF (l .EQ. 8) THEN
      DO i=1, ng
        WRITE(10,2) i, fin(i), vin(i)/top
      END DO
      DO i=1, nout
        WRITE(20,2) i, fout(i)-ana(i), vout(i)/top
      END DO
    END IF
2   FORMAT('vlvl', I3, ': f: ', E12.6, ', z: ', E12.6)
    PRINT*,'Interpolation status:', istatus

    ! Calculate error:
    ierr = 0
    err(l) = 0.0_r_kind
    DO i=1, nout
      IF (ABS(fout(i) - ana(i)) .GT. err(l)) THEN
        ierr = ierr + 1
        err(l) = ABS(fout(i) - ana(i))
      END IF
    END DO
    WRITE(*,5) l, err(l), MAXVAL(ABS(fout - ana)), ierr, nout
5   FORMAT('Error at level', I2, ': ', E12.6, ' max error: ', E12.6, ' at ', I3, ' points out of ', I3)

    DEALLOCATE(fin, vin)

  END DO

  DO l=2,lvl
    IF (err(l-1)/err(l) .LT. 3.6) THEN
      istatus = -1  ! Fail if the error does not decrease significantly
      WRITE(*,1) l, err(l-1)/err(l), err(l-1), err(l), mpddGlob%myrank
1     FORMAT('Test level', I2, ' failed. Error reduction ratio:', E12.4,' with errs: ',2E12.4, ' pc: ',I1)
      EXIT
    ELSE
      IF (mpddGlob%isBaseProc()) THEN
        WRITE(*,3) l,  err(l-1)/err(l), err(l-1), err(l), mpddGlob%myrank
3       FORMAT('Test level', I2, ' passed with error reduction ratio:', E12.4, ' with errs: ',2E12.4, ' pc: ',I1)
      END IF
    END IF
  END DO

  IF (istatus .EQ. 0) THEN
    IF (mpddGlob%isBaseProc()) WRITE(*,*) 'Test passed'
  ELSE
    IF (mpddGlob%isBaseProc()) WRITE(*,*) 'Test failed with failure code:', istatus
  END IF

  DEALLOCATE(ana, fout, vout)

  CALL mpddGlob%finalize()  ! Finalize the mpdd

END PROGRAM Test_VerticalInterp