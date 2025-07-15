&icosTr_grid
  glvl = 4
  slvl = 2
/
&latlon_grid
  num_grid = 1025,1025
  domain_latlon = 21.0, 25.0, 111.5, 115.0
/
&analysis_para
  start_time = 2022,  2, 11, 13, 36,  0
  end_time = 2022,  2, 11, 15, 12,  0
/
&geometry
mpi_layout_g14 = 2, 2
mpi_layout_g58 = 4, 4
mpi_layout_g9L = 4, 4
time_steps_g14 = 9
time_steps_g58 = 17
time_steps_g9L = 17
/
&poissonSovler
  solver = 'FMV'
  nCycle = 3
  nIterPre = 3
  nIterPost = 6
  nRelax = 2
  max_band = 10
  omegas = 1.951627, 0.569688
/
&modelState
  vLevel = 60
  varList = temp, uwnd, vwnd
/
&IO
  bk_model = WRF
  NMLFileName = '/namlist.input'
  ModelFileName = '/grapesinput'
  bakFileList = wrfout_d01_2022-02-11_13:00:00, wrfout_d01_2022-02-11_14:00:00, wrfout_d01_2022-02-11_15:00:00, wrfout_d01_2022-02-11_16:00:00
  obsFileList_Surface = 20220211_1335, 20220211_1340, 20220211_1345, 20220211_1350, 20220211_1355, 20220211_1400, 20220211_1405, 20220211_1410, 20220211_1415, 20220211_1420, 20220211_1425, 20220211_1430, 20220211_1435, 20220211_1440, 20220211_1445, 20220211_1450, 20220211_1455, 20220211_1500, 20220211_1505, 20220211_1510
  obsFileList_Sound = rec_RTEMP_20220211120000_g.dat, rec_RTEMP_20220211150000_g.dat
  obsFileList_Vwpw = vwpw22021112.dat
  mgStart = 3
  mgEnd = 7
/ 
&flog
  LogFileName = '/log/test'
/
&BMat
  ScaleParaX = 30000000.0
  ScaleParaY = 30000000.0
  ScaleParaZ = 3000.0
  ScaleParaT = 1.0
/
&Minimization
  MaxOptStep=100
/
C poissonSovler
C	solver:	     	multigrid scheme, VVV: w cycle
C	nCycle:		    number of multigrid cycles
C	nIterPre:		  number of iterations in pre-cycle
C	nIterPost:		number of iterations in post-cycle
C	nRelax:		    number of over relaxations
C	omegas:	    	relaxations coefficients
C max_band:     max band of poisson coefficent matrix 
C icosTr_grid
C glvl: icos triangle grid G level;
C slvl: Icos triangle grid start level;
C ordr: order of accuracy;
C latlon_grid
C domain_latlon specifies the (lat0,lat1) X (lon0,lon1)
C geometry
C mpi_layout:  layout of mpi
