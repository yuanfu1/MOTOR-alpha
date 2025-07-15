!!--------------------------------------------------------------------------------------------------
! PROJECT         : MOTOR-DA.possionSolver_Kp.psMatrix
! AFFILIATION     : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring 
!                   Warning and Forecasting (GBA-MWF) Shenzhen Institute of Meteorological Innovation
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.0
! HISTORY           :
!   Created by Yuanfu Xie, on 2024-12-05
!!--------------------------------------------------------------------------------------------------
!> @brief
!! This module contains a solver of streamfunction and velocity potential from vorticity and divergence
!! with u and v as boundary conditions.
!! @copyright (C) 2024 GBA-MWF, All rights reserved.
!! @note
!! @warning
!! @attention
MODULE mixedPsiChiSolver_m
  USE kinds_m, ONLY: i_kind, r_kind
  USE geometry_m, ONLY: geometry_t
  USE gzm_m, ONLY: gzm_t

  TYPE :: mixedPsiChiSolver_t
    INTEGER(i_kind) :: mgStart,mgEnd
    INTEGER(i_kind) :: numVLevels,numCells
    TYPE(gzm_t), ALLOCATABLE :: gzm(:)
    CONTAINS
      PROCEDURE :: initialize_s
      PROCEDURE :: mixedSolver_s
      PROCEDURE :: rangeSpaceProjection_s
  END TYPE mixedPsiChiSolver_t

  CONTAINS
    SUBROUTINE initialize_s(this,mgStart,mgEnd,geo)
      CLASS(mixedPsiChiSolver_t) :: this
      INTEGER(i_kind), INTENT(IN) :: mgStart,mgEnd
      TYPE(geometry_t), INTENT(IN) :: geo

      ! Allocate memory and gzms:
      IF (mgStart .LT. geo%mg%mg_coarsest .OR. mgEnd .GT. geo%mg%mg_finest) THEN
        WRITE(*,1) mgStart,mgEnd,geo%mg%mg_coarsest,geo%mg%mg_finest
1       FORMAT('Multigrid levels: ',2I3,' requested are out of the geometry setting: ',2I3)
        STOP
      END IF

      ALLOCATE(this%gzm(mgStart:mgEnd))
      DO i=mgStart,mgEnd
        this%gzm(i) = gzm_t(geo%mg%sg(i))
      END DO
      this%numVLevels = geo%mg%sg(mgEnd)%vLevel
      this%numCells = geo%mg%sg(mgEnd)%num_cell

      this%mgStart = mgStart
      this%mgEnd = mgEnd
    END SUBROUTINE initialize_s

    ! Apply a multigrid method solving a mixed Poisson equations
    ! for Psi and Chi with u and v boundary conditions:
    SUBROUTINE mixedSolver_s(this, geo, vor, div, u, v, psi, chi, useFic)
      USE mpddGlob_m, ONLY: mpddGlob_t
      USE geometry_m, ONLY: geometry_t
      USE State_m, ONLY: State_t
      USE State2NC_m
      USE YAMLRead_m
      USE parameters_m, ONLY: EarthRadius

      IMPLICIT NONE
      CLASS(mixedPsiChiSolver_t) :: this
      TYPE(geometry_t), INTENT(IN) :: geo   ! Geometry is used to implement a multigrid solver
      REAL(r_kind), INTENT(IN) :: &
        vor(this%numVLevels,this%numCells), &
        div(this%numVLevels,this%numCells), &
        u(this%numVLevels,this%numCells), &
        v(this%numVLevels,this%numCells)
      
      REAL(r_kind), INTENT(INOUT) :: &
        psi(this%numVLevels,this%numCells), &
        chi(this%numVLevels,this%numCells)

      ! useFic: 0 - no use of Laplacian equations at fictitious points; 
      !         1 - use Laplacian
      !         2 - use velocity at the actual boundaries
      INTEGER(i_kind), INTENT(IN) :: useFic 

      ! Local variables:
      EXTERNAL :: ddot,dotprd,psolve,AdotX4PsiChi_s
      INTEGER(i_kind) :: glevel,vlevel,ncells,ngmres,mgmres,its,info,half,i,j,m
      REAL(r_kind) :: res,del,scaling(2)  ! Scaling for interior and BC
      REAL(r_kind), ALLOCATABLE :: rhs(:,:),ua(:),va(:)   ! Right hand sides for each multigrid levels
      REAL(r_kind), ALLOCATABLE :: xgmres(:,:),hgmres(:,:,:),vgmres(:,:,:)

      ! Arrays hold the boundary values:
      REAL(r_kind), ALLOCATABLE :: uwind(:,:),vwind(:,:)

      ! For plot: to delete when the debugging is done
      !========================================================
      CHARACTER(LEN=1024) :: yamlFile, ncOutputFile
      CHARACTER(LEN=20) :: task
      INTEGER(i_kind) :: istatus
      TYPE(State_t) :: X
      TYPE(mpddGlob_t) :: mpddGlob
      REAL(r_kind), ALLOCATABLE :: vor_num(:,:),div_num(:,:)
      !========================================================

      CALL getarg(1, yamlFile)

      its = 80000
      res = 1.0D-18 ! on input: |Ax-b|/|b|<res;
                  ! on exit:  residual reached
      del = 1.0D1 ! on input: if(del>0) then the x returned is the hookstep
                  ! on exit:  norm of next b predicted by hook
                  ! on exit:  number of its taken
      info = 1    ! on input: if(info==1) print* residuals
                  ! if(info==2) recalc hookstep with new del
                  ! on exit:  0 sucessful, 1 method breakdown, 2 max its

      ! Apply a multigrid solver:
      DO glevel=this%mgEnd,this%mgEnd

        ! Plot the solutions: to delete when the debugging is doneto delete when the debugging is done
        !========================================================
        istatus = yaml_get_var(TRIM(yamlFile), 'IO', 'output_dir', ncOutputFile)
        istatus = yaml_get_var(TRIM(yamlFile), 'RunMode', 'Task', task)
        CALL X%initialize(yamlFile,geo%mg%sg(glevel))
        !========================================================
        
        vlevel = geo%mg%sg(glevel)%vLevel
        ncells = geo%mg%sg(glevel)%num_cell
        half = vlevel*ncells
        ngmres = 2*half ! Number of variables to solve by GMRES
        mgmres = 100    ! GMRES dimension

        ALLOCATE(rhs(ngmres,this%mgStart:this%mgEnd), &
          xgmres(ngmres,this%mgStart:this%mgEnd), &
          hgmres(mgmres+1,mgmres,this%mgStart:this%mgEnd), &
          vgmres(ngmres,mgmres+1,this%mgStart:this%mgEnd))
        rhs = 0.0D0

        ALLOCATE(vor_num(geo%mg%sg(this%mgEnd)%vLevel,geo%mg%sg(this%mgEnd)%num_cell))
        ALLOCATE(div_num(geo%mg%sg(this%mgEnd)%vLevel,geo%mg%sg(this%mgEnd)%num_cell))
        ALLOCATE(ua(geo%mg%sg(this%mgEnd)%vLevel),va(geo%mg%sg(this%mgEnd)%vLevel))
        ALLOCATE(uwind(geo%mg%sg(this%mgEnd)%vLevel,geo%mg%sg(this%mgEnd)%num_cell), &
                 vwind(geo%mg%sg(this%mgEnd)%vLevel,geo%mg%sg(this%mgEnd)%num_cell))

        ! Testing the boundary BC calculation:
        CALL uvVelocityOnInterior(psi,chi,uwind,vwind,geo%mg%sg(this%mgEnd))
        PRINT*,'Got numerical wind',maxval(ABS(psi(1,:)))
        BLOCK
          INTEGER(i_kind) :: imx,jmx
          REAL(r_kind) :: amx
          imx = 0; jmx = 0
          amx = 0.0D0
          DO i=1,geo%mg%sg(glevel)%num_cell
            IF (geo%mg%sg(glevel)%cell_type(i) .LE. 2) then
              ! DO j=1,UBOUND(geo%mg%sg(glevel)%cell_adnb,1)
                ! IF (geo%mg%sg(glevel)%cell_adnb(j,i) .NE. 0) THEN
                  ! IF (geo%mg%sg(glevel)%cell_type(geo%mg%sg(glevel)%cell_adnb(j,i)) .EQ. 0) THEN
                    IF (ABS(v(1,i)-vwind(1,i)) .GT. amx) THEN
                      imx = i
                      jmx = j
                      amx = ABS(v(1,i)-vwind(1,i))
                    END IF
                  ! END IF
                ! END IF
              ! END DO
            END IF
          END DO
          WRITE(*,555) amx,imx,v(1,imx),vwind(1,imx)
555       FORMAT('V computed - anaV: ',E14.6,' at:',I6,' true/computed: ',2E14.6)
          PRINT*,'Calculated U component: ',MAXVAL(ABS(uwind(1,:)))
          PRINT*,'Calculated V 2: ',vwind(1,2),v(1,2)

          ! Check the velocity at corners:
          CALL this%gzm(glevel)%Laplacia(psi,vor_num,1)
          imx = 0
          amx = 0.0D0
          DO i=1,geo%mg%sg(glevel)%num_cell
            IF (ABS(vor_num(1,i)-vor(1,i)) .GT. amx) THEN
              imx = i
              amx = ABS(vor_num(1,i)-vor(1,i))
            END IF
          END DO
          WRITE(*,666) imx,amx,vor_num(1,imx),vor(1,imx)
  666     FORMAT('Max err in vor at: ',I6,' err: ',E14.6,' vors: ',2E14.6)
        END BLOCK

        ! Initial guess:
        xgmres = 1.0D0 !/DSQRT(DBLE(ngmres))
        ! Ideal guess:
        xgmres = 0.0D0
        PRINT*,'PSI Norm lvl 1: ',MAXVAL(ABS(psi(1,:))),MAXVAL(ABS(geo%mg%sg(glevel)%cell_cntr(2,:))), &
        MAXVAL(ABS(psi(1,:)))*MAXVAL(ABS(geo%mg%sg(glevel)%cell_cntr(2,:)))
        DO j=1,geo%mg%sg(glevel)%num_cell
          !IF (geo%mg%sg(glevel)%cell_type(j) .EQ. -1) THEN
            DO i=1,geo%mg%sg(glevel)%vLevel
              psi(i,j) = 1.0D0*psi(i,j) !+ geo%mg%sg(glevel)%cell_cntr(2,j)*MAXVAL(ABS(psi(1,:)))*1.0D2
            END DO
          !ELSE
          !  psi(:,j) = 0.0D0
          !END IF
        END DO

        ! Select a constant initial guess:
        psi = 1.0D0/DSQRT(dot_product(psi(1,:),psi(1,:)))
        xgmres(     1:  half,glevel) = RESHAPE(psi,(/ half /))

        X%fields(5)%DATA(:,:,1) = RESHAPE(xgmres(1:half,glevel), (/ vlevel, ncells /))
        PRINT*,'Initial guess at 1: ',xgmres(1,this%mgEnd)

        DO j=1,geo%mg%sg(glevel)%num_cell
          ! IF (geo%mg%sg(glevel)%cell_type(j) .EQ. 0) THEN ! Apply to interior points only
          IF (geo%mg%sg(glevel)%cell_type(j) .LE. 2) THEN
            DO i=1,geo%mg%sg(glevel)%vLevel
              rhs(i+geo%mg%sg(glevel)%vLevel*(j-1)     ,glevel) = vor(i,j)
              rhs(i+geo%mg%sg(glevel)%vLevel*(j-1)+half,glevel) = div(i,j)
            END DO
          END IF
        END DO

        ! Scale of all vertical levels:
        ! DO i=1,geo%mg%sg(glevel)%vLevel
        !   scaling = scaling + (dot_product(vor(i,:),vor(i,:))) ! Scale the boundary conditions
        ! END DO
        scaling = DOT_PRODUCT(rhs(:,glevel),rhs(:,glevel))
        scaling = DSQRT(scaling)

        scaling(2) = 0.0D-1 ! scaling(1)*0.0D0 ! Temporarily turn off the boundary conditions

        ! scaling = 1.0D0 ! no scaling
        PRINT*,'RHS before scaling: ',dot_product(rhs(:,glevel),rhs(:,glevel)),DSQRT(dot_product(rhs(:,glevel),rhs(:,glevel)))
        rhs = rhs/scaling(1)
        PRINT*,'Scaling: ',scaling(1),1.0D0/scaling(1), &
          UBOUND(geo%mg%sg(glevel)%edgeNorm2,1), &
          UBOUND(geo%mg%sg(glevel)%edgeTang2,2), &
          geo%mg%sg(glevel)%vLevel
        PRINT*,'Edge norm: ',geo%mg%sg(glevel)%edgeNorm2(:,1,20), &
        geo%mg%sg(glevel)%edgeTang2(:,1,20),maxval(geo%mg%sg(glevel)%edgeNorm2),minval(geo%mg%sg(glevel)%edgeNorm2)

        ! Replace the rhs value with U and V BC values:
        ! ASSOCIATE(sg => geo%mg%sg(glevel))
        !   DO i=1,sg%num_cell
        !     ! Apply BC to boundary only, find an interior point with a neighbor(s)
        !     IF (sg%cell_type(i) .EQ. 1) THEN
        !       DO j=1,UBOUND(sg%cell_adnb,1)
        !         ! Skip those neighbor cells outside the domain:
        !         IF (sg%cell_adnb(j,i) .NE. 0) THEN
                
        !           IF (sg%cell_type(sg%cell_adnb(j,i)) .EQ. 0) THEN
        !             rhs(sg%vLevel*(i-1)+1     :sg%vLevel*i,     glevel) = u(:,i)*scaling(2)
        !             rhs(sg%vLevel*(i-1)+1+half:sg%vLevel*i+half,glevel) = v(:,i)*scaling(2)
        !           END IF

        !         END IF
        !       END DO
        !     END IF
        !   END DO
        ! END ASSOCIATE

        ! Plot for debugging:
        X%fields(1)%DATA(:,:,1) = RESHAPE(rhs(1:half,glevel), (/ vlevel, ncells /))
        X%fields(2)%DATA(:,:,1) = psi
        psi = 0.0D0

        DO j=1,geo%mg%sg(glevel)%num_cell
          DO i=1,geo%mg%sg(glevel)%vLevel
            X%fields(3)%DATA(i,j,1) = vor(i,j)
          END DO
        END DO

        PRINT*,'Calling GMRES...',ngmres,this%gzm(glevel)%sg%gLevel,dot_product(rhs(:,glevel),rhs(:,glevel))
        ! Yuanfu Xie modified gmres to gmres_gzm to access gzm for AdotX4PsiChi_s:
        CALL gmresm_gzm(mgmres,ngmres,xgmres(:,glevel),rhs(:,glevel),AdotX4PsiChi_s,psolve,dotprd, &
          hgmres(:,:,glevel),vgmres(:,:,glevel),res,del,its,info, this%gzm(glevel),scaling)
        
        WRITE(*,1)
1       FORMAT('+----------------------------------------+')
        IF (info .EQ. 0) THEN
          WRITE(*,2)
2         FORMAT('|   GMRES successes:                    |')
        ELSE IF (info .EQ. 1) THEN
          WRITE(*,3) info
3         FORMAT('| GMRES breakdown, no solution output! ',I1,' |')
        ELSE
          WRITE(*,4)
4         FORMAT('| GMRES reaches max iterations!          |')
        END IF
        WRITE(*,1)

        IF (info .NE. 1000) THEN
          psi = RESHAPE(xgmres(1:half,glevel), (/ vlevel, ncells /))
          chi = RESHAPE(xgmres(half+1:2*half,glevel), (/ vlevel, ncells /))
        END IF
        psi = psi/scaling(1)
        chi = chi/scaling(1)
        PRINT*,'Norm stream: ',DSQRT(dot_product(psi(1,:),psi(1,:))),MAXVAL(ABS(psi(1,:))),geo%mg%sg(glevel)%num_cell,geo%mg%sg(glevel)%num_icell

        ! Calculate the numerical vorticity and divergence for verification:
        CALL this%gzm(glevel)%Laplacia(psi,vor_num,1)
        CALL this%gzm(glevel)%Laplacia(chi,div_num,1)
        X%fields(2)%DATA(:,:,1) = psi
        PRINT*,'vor_num corner 2: ',vor_num(1,51)

        ! Calculate the boundary values:
        CALL uvVelocityOnInterior(psi,chi,uwind,vwind,geo%mg%sg(glevel))

        ! Determining the general solution:
        BLOCK
          INTEGER(i_kind) :: interior,nintegral
          REAL(r_kind) :: a,b,c,umx,vmx,pmx,cmx,coef(2,2),dd,dll(2),matrix(2,2),r(2,2),cs,cs2,lon
          REAL(r_kind), ALLOCATABLE :: u1sum(:),ulsum(:),v1sum(:),vlsum(:),l1sum(:),l2sum(:)
          REAL(r_kind), ALLOCATABLE :: integral(:)
          ASSOCIATE(sg => geo%mg%sg(glevel))
            dll(1) = sg%cell_cntr(1,sg%cell_adnb(3,1))-sg%cell_cntr(1,1)
            dll(2) = sg%cell_cntr(2,sg%cell_adnb(2,1))-sg%cell_cntr(2,1)

            interior = 0
            umx = 0.0D0; vmx = 0.0D0; pmx = 0.0D0; cmx = 0.0D0
            ALLOCATE(u1sum(sg%vlevel),ulsum(sg%vlevel),v1sum(sg%vlevel),vlsum(sg%vlevel),l1sum(sg%vlevel),l2sum(sg%vlevel))
            ALLOCATE(integral(sg%num_cell))
            u1sum = 0.0D0; ulsum = 0.0D0; v1sum = 0.0D0; vlsum = 0.0D0; l1sum = 0.0D0; l2sum = 0.0D0
            integral = 0.0D0

            matrix = 0.0D0; r = 0.0D0
            DO i=1,sg%num_cell

              ! Integral of 1/cos(theta): THIS WORKS ONLY FOR Lat-lon grid
              nintegral = NINT((sg%cell_cntr(1,i)-sg%cell_cntr(1,1))/dll(1))
              integral(i) = 0.0D0
              DO j=1,nintegral
                ! Starting the interior cell south boundary:
                integral(i) = integral(i) + dll(1)/DCOS(sg%cell_cntr(1,1)+(DBLE(j)-0.5D0)*dll(1))

                ! Plotting the integral
                X%fields(1)%DATA(1,i,3) = integral(i)
              END DO

              ! For all interior point to find the least square solution determining the general solution:
              IF (sg%cell_type(i) .LE. 1) THEN
                interior = interior+1
                IF (umx .LT. ABS(uwind(1,i))) umx = ABS(uwind(1,i))
                IF (vmx .LT. ABS(vwind(1,i))) vmx = ABS(vwind(1,i))
                IF (pmx .LT. ABS(psi(1,i))) pmx = ABS(psi(1,i))
                IF (cmx .LT. ABS(chi(1,i))) cmx = ABS(chi(1,i))

                u1sum = u1sum + DCOS(sg%cell_cntr(1,i))*(u(:,i)-uwind(:,i))
                ulsum = ulsum + DCOS(sg%cell_cntr(1,i))*(u(:,i)-uwind(:,i))*sg%cell_cntr(2,i)
                v1sum = v1sum + DCOS(sg%cell_cntr(1,i))*(v(:,i)-vwind(:,i))
                vlsum = vlsum + DCOS(sg%cell_cntr(1,i))*(v(:,i)-vwind(:,i))*sg%cell_cntr(2,i)
                l1sum = l1sum + sg%cell_cntr(2,i)
                l2sum = l2sum + sg%cell_cntr(2,i)**2

                ! The corrent coefficient matrix:
                cs = DCOS(sg%cell_cntr(1,i))
                cs2 = cs**2
                lon = sg%cell_cntr(2,i)
                matrix(1,1) = matrix(1,1) + 1.0D0/cs2
                matrix(1,2) = matrix(1,2) + lon/cs2
                matrix(2,2) = matrix(2,2) + lon**2/cs2
                r(1,1) = r(1,1) +     (u(1,i)-uwind(1,i))/cs
                r(2,1) = r(2,1) - lon*(u(1,i)-uwind(1,i))/cs
                r(1,2) = r(1,2) +     (v(1,i)-vwind(1,i))/cs
                r(2,2) = r(2,2) + lon*(v(1,i)-vwind(1,i))/cs
              END IF
            END DO

            a = DBLE(interior); b = l1sum(1); c = l2sum(1); dd = a*c-b*b
            coef(1,1) = (c*u1sum(1)+b*ulsum(1))/dd
            coef(2,1) = (b*u1sum(1)+a*ulsum(1))/dd
            coef(1,2) = (c*v1sum(1)-b*vlsum(1))/dd
            coef(2,2) =(-b*v1sum(1)+a*vlsum(1))/dd

            WRITE(*,20) a,-b,u1sum(1),-b,c,ulsum(1)
    20      FORMAT('Coefficient matrix 1: ',2E14.6,' b1: ',E14.6)
            WRITE(*,21) a,b,v1sum(1),b,c,vlsum(1)
    21      FORMAT('Coefficient matrix 2: ',2E14.6,' b2: ',E14.6)
            WRITE(*,23) coef(1:2,1),coef(1:2,2)
    23      FORMAT('Solution to projection: ',4E14.6)
            v1sum = v1sum/DBLE(interior)
            PRINT*,'The mean constant: ',v1sum(1),interior,minval(sg%cell_cntr(2,:)),maxval(sg%cell_cntr(2,:))
            PRINT*,'Max abs u and v: ',umx,vmx,pmx,cmx

            ! Correct form:
            dd = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(1,2)
            coef(1,1) = (matrix(2,2)*r(1,1)+matrix(1,2)*r(2,1))/dd
            coef(1,2) = (matrix(1,2)*r(1,1)+matrix(1,1)*r(2,1))/dd
            coef(2,1) = (matrix(2,2)*r(1,2)-matrix(1,2)*r(2,2))/dd
            coef(2,2) =(-matrix(1,2)*r(1,2)+matrix(1,1)*r(2,2))/dd

            ! Add the general solution component to psi: consider v only
            ! DO i=1,sg%num_cell
            !   psi(:,i) = psi(:,i) + v1sum*EarthRadius*sg%cell_cntr(2,i)
            ! END DO
            ! Add the general solutions to both psi and chi:
            DO i=1,sg%num_cell
              psi(:,i) = psi(:,i) + EarthRadius*sg%cell_cntr(2,i)*(coef(2,1)+1.0D0*coef(1,2)*integral(i))
              chi(:,i) = chi(:,i) + EarthRadius*sg%cell_cntr(2,i)*(coef(1,1)+1.0D0*coef(2,2)*integral(i))
            END DO

            ! Check if the U component becomes smaller:
            u1sum = 0.0D0
            v1sum = 0.0D0
            DO i=1,sg%num_cell
              IF (sg%cell_type(i) .EQ. 0) &
                u1sum(1) = u1sum(1) + (uwind(1,i) + (coef(1,1)-coef(1,2))*sg%cell_cntr(2,i)/DCOS(sg%cell_cntr(1,i)))**2
                v1sum(1) = v1sum(1) + uwind(1,i)**2
            END DO
            PRINT*,'Modified U sum: ',u1sum(1),v1sum(1),maxval(ABS(uwind(1,:))),coef(1,:)

            DEALLOCATE(u1sum,ulsum,v1sum,vlsum,l1sum,l2sum,integral)
          END ASSOCIATE
        END BLOCK
        ! Recalculate the boundary values after adding the general solution:
        CALL uvVelocityOnInterior(psi,chi,uwind,vwind,geo%mg%sg(glevel))
        PRINT*,'Recalculated the wind...'

        ! Plot the residual of vorticity field:
        BLOCK
          INTEGER(i_kind) :: imx,inx
          REAL(r_kind) :: amx,anx,alow,high

          imx = 0; inx = 0
          amx = 0.0D0; anx = 0.0D0
          ! psi = RESHAPE(rhs(1:half,glevel), (/ vlevel, ncells /))
          ! psi = psi*scaling(1)
          ! chi = RESHAPE(rhs(half+1:2*half,glevel), (/ vlevel, ncells /))
          ! chi = chi*scaling(1)
          DO j=1,geo%mg%sg(glevel)%num_cell
            IF (geo%mg%sg(glevel)%cell_type(j) .LE. 2) THEN
              DO i=1,geo%mg%sg(glevel)%vLevel
                X%fields(4)%DATA(i,j,1) = vor_num(i,j)/scaling(1)
              END DO
              IF (ABS(vor(1,j)) .GT. amx) THEN
                imx = j; amx = ABS(vor(1,j))
              END IF
              IF (ABS(vor_num(1,j)-vor(1,j)) .GT. anx) THEN
                inx = j; anx = ABS(vor_num(1,j)-vor(1,j))
              END IF
            END IF
          END DO
          PRINT*,'MX residual in vor: ',inx,anx,vor_num(1,inx),vor(1,inx),geo%mg%sg(glevel)%num_cell

          ! For interior points, plot the residual:
          X%fields(5)%DATA(:,:,1) = vor_num - vor

          ! Clean the boundary values:
          DO j=1,geo%mg%sg(glevel)%num_cell
            IF (geo%mg%sg(glevel)%cell_type(j) .GT. 2) THEN
              DO i=1,geo%mg%sg(glevel)%vLevel
                X%fields(5)%DATA(i,j,1) = 0.0D0
              END DO
            END IF
          END DO

          ! For the boundary points, plot the residual:
          X%fields(5)%DATA(:,:,2) = 0.0D0
          ASSOCIATE(sg => geo%mg%sg(glevel))
            imx = 0; amx = 0.0D0; inx = 0
            alow = 1.0D10; high = -1.0D10
            
            DO i=1,sg%num_cell
              IF (sg%cell_type(i) .LE. 2) THEN
                X%fields(5)%DATA(:,i,2) = (vwind(:,i) - v(:,i))
                X%fields(1)%DATA(:,i,2) = uwind(:,i)
                X%fields(2)%DATA(:,i,2) = u(:,i)
                X%fields(3)%DATA(:,i,2) = vwind(:,i)
                X%fields(4)%DATA(:,i,2) = v(:,i)
                IF (ABS(X%fields(5)%DATA(1,i,2)) .GT. amx) THEN
                  imx = i
                  amx = ABS(X%fields(5)%DATA(1,i,2))
                END IF
              END IF
            END DO
            IF (imx .NE. 0) THEN
              WRITE(*,333) imx,amx,vwind(1,imx),v(1,imx)
333           FORMAT('Max diff in v: ',I6,' val: ',E14.6,' V: ',2E14.6)
            ELSE
              PRINT*,'No difference at the boundaries! ',scaling(2)
            END IF
          END ASSOCIATE
        END BLOCK

        CALL Output_NC_State_AV(X, TRIM(ncOutputFile), TRIM(task)//"_rhs", .TRUE., .TRUE.)

        ! WRITE(*,111) xgmres(67*68,glevel),psi(2,68),xgmres(1,glevel),psi(1,1)
111     FORMAT('XGMRES psi at 67*68: ',2E14.6,' at 1: ',2E14.6)
        ! WRITE(*,222) psi(1,63:67)
222     FORMAT('Psi at right corner: ',5E14.6)

        ! psi = RESHAPE(xgmres(1:half,glevel), (/ vlevel, ncells /))
        ! chi = RESHAPE(xgmres(half+1:2*half,glevel), (/ vlevel, ncells /))
        ! psi = psi/scaling(1)
        ! chi = chi/scaling(1)
        PRINT*,'Norm stream 2: ',DSQRT(dot_product(psi(1,:),psi(1,:))),MAXVAL(ABS(psi(1,:))),vlevel,ncells,half
      END DO

      ! Deallocate memory:
      DEALLOCATE(rhs,xgmres,hgmres,vgmres,vor_num,div_num,ua,va,uwind,vwind)

    END SUBROUTINE mixedSolver_s

    ! Project the residual in the null space out from the solution:
    SUBROUTINE rangeSpaceProjection_s(this,n,x,y)
      CLASS(mixedPsiChiSolver_t) :: this
      INTEGER(i_kind), INTENT(IN) :: n
      REAL(r_kind), INTENT(IN) :: x(n)
      REAL(r_kind), INTENT(OUT) :: y(n)

      ! Local variables:
      INTEGER(i_kind) :: i,j,k,m,nn,ioff(2)
      REAL(r_kind) :: s(n),z(n),d(n),xx,yy

      m = DSQRT(DBLE(n/2))
      ioff(1) = 0
      ioff(2) = m

      ! Project out the constant:
      d(1:n/2) = x(1:n/2)-SUM(x(1:n/2))/DBLE(n/2)
      d(n/2+1:n) = x(n/2+1:n)-SUM(x(n/2+1:n))/DBLE(n/2)
      y = d
      PRINT*,'MaxMin values: ',maxval(d),minval(d)
      PRINT*,'MaxMin values: ',maxval(x),minval(x)

      ! Project out the linear function:
      k = 0
      s = 0.0D0
      z = 0.0D0
      DO j=1,m
          DO i=1,m
              k = k+1
              s(k+ioff(1)) =  DBLE(i-1)/DBLE(m-1)-0.5D0
              z(k+ioff(1)) = -(DBLE(j-1)/DBLE(m-1)-0.5D0)
              s(k+ioff(2)) =  DBLE(j-1)/DBLE(m-1)-0.5D0
              z(k+ioff(2)) =  DBLE(i-1)/DBLE(m-1)-0.5D0
          END DO
      END DO
      d(1:n/2) = y(1:n/2) - DOT_PRODUCT(y(1:n/2),s(1:n/2))/DOT_PRODUCT(s(1:n/2),s(1:n/2))*s(1:n/2)
      y(1:n/2) = d(1:n/2) - DOT_PRODUCT(d(1:n/2),z(1:n/2))/DOT_PRODUCT(z(1:n/2),z(1:n/2))*z(1:n/2)
      d(n/2+1:n) = y(n/2+1:n) - DOT_PRODUCT(y(n/2+1:n),s(n/2+1:n))/DOT_PRODUCT(s(n/2+1:n),s(n/2+1:n))*s(n/2+1:n)
      y(n/2+1:n) = d(n/2+1:n) - DOT_PRODUCT(d(n/2+1:n),z(n/2+1:n))/DOT_PRODUCT(z(n/2+1:n),z(n/2+1:n))*z(n/2+1:n)

    END SUBROUTINE rangeSpaceProjection_s
END MODULE mixedPsiChiSolver_m