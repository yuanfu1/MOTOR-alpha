!
MODULE module_interpolation_m
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Description:
!   Module for Interpolation
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!  USE  namelist_define
!  use  obs_type_define
!  use  Meteoro_constants

  IMPLICIT NONE
!-----------------------------------------------------------------------------

  TYPE interpolation_t

  CONTAINS

    PROCEDURE, PUBLIC  :: Intp2DWeight
    PROCEDURE, PUBLIC  :: Intp2DWeight2
    PROCEDURE, PUBLIC  :: VerticalWeight
    PROCEDURE, PUBLIC  :: VerticalWeight_h
    PROCEDURE, PUBLIC  :: Gposition
    PROCEDURE, PUBLIC  :: UvConvert

  END TYPE interpolation_t

CONTAINS

  SUBROUTINE Intp2DWeight(lat_ob, lon_ob, lats, lonw, dlat, dlon, mi, mj, &
                          i, j, alpha, alpha1, beta, beta1)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Discription:
!    Transfer obs. to grid and calculate its nondimensional
!    distance
!
!  Method:
!    obs. coord at G
!    oy = (lat_ob-latS)/dlat
!    ox = (lon_ob-latW)/dlon
!        where:
!    Lon_ob    :  longitude degree of obs at G.
!    Lat_ob    :  latitude degree of obs at G.
!    latS      :  the start lat. degree of grid,  at the southwestern corner;
!    latW      :  the start lon. degree of grid.
!    dlat      :  mesh space in degree in y direction
!    dlon      :  mesh space in degree in x direction
!    mi        :  Total number of grid in y direction
!    mj        :  Total number of grid in x direction
!
!    i+1,j                 i+1,j+1
!    |-----------------------|-^-
!    |           |           | |
!    |           |           |alpha1
!    |           |           | |
!    |           |           | |
!    |-----------G-----------|-v-
!    |           |           | ^
!    |           |           | |
!    |           |           | |
!    |           |           |alpha
!    |           |           | |
!    |           |           | |
!    |-----------|-----------|-v-
!    |<---beta-->|<--beta1-->|
!    i,j                   i,j+1
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IMPLICIT NONE

!* Subroutine arguments:
    REAL, INTENT(in)  :: lat_ob, lon_ob
    REAL, INTENT(in)  :: dlat, dlon, lats, lonw
    INTEGER, INTENT(in)  :: mi, mj

    REAL, INTENT(out) :: alpha, alpha1, beta, beta1
    INTEGER, INTENT(out) :: i, j

!* Local parameters:
    REAL                                :: ox, oy
!- End of header ----------------------------------------------

    oy = (lat_ob - lats) / dlat + 1.0
    ox = (lon_ob - lonw) / dlon + 1.0

    i = INT(oy + 1.0E-20)
    j = INT(ox + 1.0E-20)

    alpha = oy - REAL(i)
!  if(i==0)alpha=1.0-alpha
    alpha1 = 1.0 - alpha

    beta = ox - REAL(j)
    beta1 = 1.0 - beta

! north pole
    IF (i >= mi) THEN
      alpha = alpha * 2
      IF (alpha > 1.0) THEN
        WRITE (0, *) 'lat_ob,lon_ob,lats,lonw,dlat,dlon,oy,i,ox,j,alpha,mi,=', lat_ob, lon_ob, lats, lonw, dlat, dlon, oy, i, ox, j, alpha, mi, &
          '  stop  alpha is too large '
      END IF
      alpha1 = 1.0 - alpha
    END IF
! south pole
    IF (i < 1) THEN
      alpha = alpha * 2 - 1
      IF (alpha > 1.0) THEN
        WRITE (0, *) 'lat_ob,lon_ob,lats,lonw,dlat,dlon,oy,i,ox,j,alpha,mi=', lat_ob, lon_ob, lats, lonw, dlat, dlon, oy, i, ox, j, alpha, mi, &
          'stop  alpha is too large '
      END IF
      IF (alpha < 0.0) THEN
        WRITE (0, *) 'lat_ob,lon_ob,lats,lonw,dlat,dlon,oy,i,ox,j,alpha=', lat_ob, lon_ob, lats, lonw, dlat, dlon, oy, i, ox, j, alpha, &
          'stop  alpha is too small'
      END IF
      alpha1 = 1.0 - alpha
    END IF

  END SUBROUTINE Intp2DWeight

  SUBROUTINE VerticalWeight(pre, ko, dp, dpm)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Discription:
!    Looking for the location of thr pre and calculate its weight in
!    vertical direction.
!
!  Method:
!   Bi-linear
!
! Current code owner:   RCNMP
!
! History:
! Version      Date       Comment
! --------   ---------    --------
!   1.0      2001.12.01   Coding           ZHANG HUA
!   1.0      2001.12.01   Debug & test     ZHANG HUA
!   1.8      2003.04.09   modified         Zhuang Shiyu
! Code Description:
!   Language:           Fortran 90
!   Software Standard:
!
! Parent module:
!   module_setupstructures
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IMPLICIT NONE

    REAL, INTENT(in)  :: pre     ! Presure of Obs.
    INTEGER, INTENT(out) :: ko      ! Level of Obs.
    REAL, INTENT(out) :: dp, dpm ! Weight

    REAL                   :: dp1, dp2, dpx
    INTEGER                :: m       ! loop counter
!---------------------------------------------------------------------------

    IF (pre >= p_mandat(1)) ko = 1       ! below the first modle level(1000hPa)

    DO m = 1, mkp - 1
      IF (pre < p_mandat(m) .AND. pre >= p_mandat(m + 1)) ko = m
    END DO

    IF (pre < p_mandat(mkp)) ko = mkp - 1  ! above the top stadard level(10hPa)

    dp1 = LOG(pre / p_mandat(ko + 1))
    dpx = LOG(p_mandat(ko) / p_mandat(ko + 1))

    dp2 = dpx - dp1

    dp = dp1 / dpx
    dpm = dp2 / dpx

  END SUBROUTINE VerticalWeight

  SUBROUTINE VerticalWeight_h(height, xbheight, ko, dp, dpm, mkp)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Discription:
!    Looking for the location of thr height and calculate its weight in
!    vertical direction.
!
!  Method:
!  linear
!
! Current code owner:   RCNMP
!
! History:
! Version      Date       Comment
! writted      20140404   taoshw
! --------   ---------    --------
! Code Description:
!   Language:           Fortran 90
!   Software Standard:
!
! Parent module:
!   module_setupstructures
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IMPLICIT NONE

    REAL, INTENT(in)  :: height     ! RPresure of Obs.
    INTEGER, INTENT(in) :: mkp     ! number of Levels of xb .
    INTEGER, INTENT(out) :: ko      ! Level NO. of Obs in xb.
    REAL, INTENT(out) :: dp, dpm ! Weight
    REAL, DIMENSION(mkp) :: xbheight

    REAL                   :: dp1, dp2, dpx
    INTEGER                :: m       ! loop counter
!---------------------------------------------------------------------------

    IF (height <= xbheight(1)) ko = 1       ! below the first modle level(1000hPa)

    DO m = 1, mkp - 1
      IF (height > xbheight(m) .AND. height <= xbheight(m + 1)) ko = m
    END DO

    IF (height > xbheight(mkp)) ko = mkp - 1  ! above the top stadard level(10hPa)

    dp1 = height - xbheight(ko + 1)
    dpx = xbheight(ko) - xbheight(ko + 1)

    dp2 = dpx - dp1

    dp = dp1 / dpx
    dpm = dp2 / dpx

  END SUBROUTINE VerticalWeight_h

  SUBROUTINE Gposition(iy, jx, i, j, isw, ise, inw, ine, jsw, jse, jnw, jne)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Description:
! Positioning for Global horizontal interpolation
!
!
! Current code owner: GRAPES/CNPR
!
! History:
! Version Date Comment
! -------- --------- --------
! 2.5 2006.08.15 Original Zhuang shiyu
!
! Code Description:
! Language: Fortran 90
! Software Standard:
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    IMPLICIT NONE
!*Subroutine arguments:
    INTEGER, INTENT(in) :: i, j, iy, jx
    INTEGER, INTENT(out) :: isw, ise, inw, ine, jsw, jse, jnw, jne
    INTEGER :: jmd
!! positioning for global horizontal interpolation
    jmd = jx / 2
    isw = i
    ise = i
    inw = i + 1
    ine = i + 1
    jsw = j
    jse = j + 1
    jnw = j
    jne = j + 1
    IF (j == jx) THEN
      jse = 1
      jne = 1
    END IF

    IF (i < 1) THEN
      isw = 1
      ise = 1
      inw = 1
      ine = 1
      jnw = j + 1 + jmd
      jne = j + jmd
      IF (jnw > jx) jnw = jnw - jx
      IF (jne > jx) jne = jne - jx
    END IF

    IF (i >= iy) THEN
      isw = iy
      ise = iy
      inw = iy
      ine = iy
      jnw = j + 1 + jmd
      jne = j + jmd
      IF (jnw > jx) jnw = jnw - jx
      IF (jne > jx) jne = jne - jx
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE Gposition

  SUBROUTINE UvConvert(lat, lon, U, V, MTYPE)
!--------------------------------------------------------------------
! Purpose
! Projection of the wind field between two different coordinate
!
! MTYPE=1 mapping U, V from Latlon into Cartesian coordinate
! MTYPE=2 mappong U, V from Cartesian into Latlon coordiante
!
! 2003/12/03 Provider : Zhu Zhongshen
! 2003/1204/ modified and adjoint coding : Zhuang Shiyu
!
!--------------------------------------------------------------------
!
    REAL, INTENT(in) :: lat, lon
    REAL, INTENT(inout) :: U, V
    INTEGER, INTENT(in) :: MTYPE
    REAL :: SX, CX
    REAL :: WX, WY
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    SX = SIN(lon * pi180)
    CX = COS(lon * pi180)
    IF (lat < 0.) CX = -CX
    SELECT CASE (MTYPE)
    CASE (1)
      WX = U * CX - V * SX
      WY = U * SX + V * CX
      U = WX
      V = WY
    CASE (2)
      WX = U * CX + V * SX
      WY = -U * SX + V * CX
      U = WX
      V = WY
    END SELECT
    RETURN
  END SUBROUTINE UvConvert

!END MODULE module_interpolation

  SUBROUTINE Intp2DWeight2(lat_ob, lon_ob, slat, wlon, dlat, dlon, iy, jx, &
                           a1, a2, a3, a4, isw, ise, inw, ine, jsw, jse, jnw, jne, id)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!  Discription:
!    Transfer obs. x to grid j and calculate its nondimensional
!    distance to grid j and j+1
!
!  Method:
!    obs. coord at G
!    oy = (lat_ob-latS)/dlat
!    ox = (lon_ob-latW)/dlon
!        where:
!    Lon_ob    :  longitude degree of obs at G.
!    Lat_ob    :  latitude degree of obs at G.
!    latS      :  the start lat. degree of grid,  at the southwestern corner;
!    latW      :  the start lon. degree of grid.
!    dlat      :  mesh space in degree in y direction
!    dlon      :  mesh space in degree in x direction
!    iy        :  Total number of grid in y direction
!    jx        :  Total number of grid in x direction
!    id        :  Option for u(pi,th,hum) and v. 1 for u/mss/hum, others for v
!
!   Define the position of Obs. in the interpolation.
!   Bilinear interpolation degrades to tringle interpolation at poles.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    IMPLICIT NONE

!* Subroutine arguments:
    REAL, INTENT(in)  :: lat_ob, lon_ob
    REAL, INTENT(in)  :: dlat, dlon, slat, wlon
    INTEGER, INTENT(in)  :: iy, jx
    REAL, INTENT(out) :: a1, a2, a3, a4
    INTEGER, INTENT(out) :: isw, ise, inw, ine, jsw, jse, jnw, jne
    INTEGER, INTENT(in)  :: id

!* Local parameters:
    INTEGER :: i, j
    REAL    :: ox, oy
    REAL    :: alpha, beta, alpha1, beta1, atmp      ! Weight

!-------------------------------------------------------------------------------
! [1.0] Calculation of the distance from the obs. to the most sounthwestern grid
!-------------------------------------------------------------------------------
    oy = (lat_ob - slat) / dlat + 1.0
    ox = (lon_ob - wlon) / dlon + 1.0

    i = INT(oy + 1.0E-8)
    j = INT(ox + 1.0E-8)

    beta = ox - REAL(j)
    beta1 = 1.0 - beta

    alpha = oy - REAL(i)

    IF (id .EQ. 0) THEN                     ! 0 for v, 1 for u, 2 for mss/hum;
      IF (i >= iy - 1) alpha = alpha / 2.0          ! for north polar region
      IF (i < 1) alpha = (1.0 - alpha) / 2.0      ! for south polar region
    ELSE
      IF (i < 1) alpha = 1.0 - alpha              ! for south polar region
    END IF
    alpha1 = 1.0 - alpha

!------------------------------------------------------------------------------
! [2.0] calculate the interpolation weight coefficients.
!------------------------------------------------------------------------------
    a1 = alpha1 * beta1; a2 = alpha1 * beta
    a3 = alpha * beta; a4 = alpha * beta1

    SELECT CASE (id)
    CASE (0)      ! for v
      IF (i == 0 .OR. i == iy - 1) THEN
        atmp = a3
        a3 = -a4
        a4 = -atmp
      END IF

    CASE (1)     ! for u
      IF (i == 0 .OR. i == iy) THEN
        atmp = a3
        a3 = -a4
        a4 = -atmp
      END IF

    CASE (2)    ! for others
      IF (i == 0 .OR. i == iy) THEN
        atmp = a3
        a3 = a4
        a4 = atmp
      END IF
    END SELECT

!------------------------------------------------------------------------------
! [3.0] define the location of four grids
!------------------------------------------------------------------------------
    isw = i; jsw = j
    ise = i; jse = j + 1
    ine = i + 1; jne = j + 1
    inw = i + 1; jnw = j

!expand east and west boundary

    IF (j == jx) THEN
      jse = jse - jx; jne = jne - jx
    END IF

    IF (j > jx) THEN
      jsw = jsw - jx; jnw = jnw - jx
      jse = jse - jx; jne = jne - jx
    END IF

    IF (j < 1) THEN
      jsw = jsw + jx; jnw = jnw + jx
    END IF

!!  cross the poles
    IF (id .EQ. 0) THEN            ! 0 for v, 1 for u, 2 for mss/hum;
      IF (i >= iy - 1) THEN
        ine = isw; jne = jsw + jx / 2
        inw = ise; jnw = jse + jx / 2
      END IF

    ELSE
      IF (i >= iy) THEN
        ine = isw; jne = jsw + jx / 2
        inw = ise; jnw = jse + jx / 2
      END IF
    END IF

    IF (i < 1) THEN
      inw = 1; ine = 1
      isw = 1; ise = 1
      jne = jsw + jx / 2; jnw = jse + jx / 2
    END IF
    IF (jne > jx) jne = jne - jx
    IF (jnw > jx) jnw = jnw - jx

  END SUBROUTINE Intp2DWeight2

END MODULE module_interpolation_m
