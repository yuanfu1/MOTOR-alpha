&icosTr_grid
  glvl = 4
  slvl = 2
/

&latlon_grid
  num_grid = 1025,1025
  domain_latlon = 18, 27, 107, 119
/

&analysis_para
  start_time = 2022, 1, 25, 12, 0, 0
  end_time = 2022, 1, 25, 12, 10, 0
/

&geometry
  mpi_layout_g14 = 1, 1
  mpi_layout_g58 = 1, 1
  mpi_layout_g9L = 1, 1
  time_steps_g14 = 2
  time_steps_g58 = 2
  time_steps_g9L = 2
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
  varList = temp
/

&obs_thinning
  obs_qcThreshold = 5.0D2, 5.0D2, 5.0D2
  obs_missing = 999999999.0D0
  obs_invalid = 3.4D38
  obs_radius = 5.0D0, 5.0D0, 5.0D2, 5.0D1
  thinning_threshold = 0.25D0
  interpolation = 1
/

&IO
  bk_model = GRAPES
  ensModel = WRF
  NMLFileName = '/namlist.input'
  ModelFileName = '/grapesinput'
  bakFileList = wrfout_d01_2022-01-25_12:00:00, wrfout_d01_2022-01-26_00:00:00
  ObsFlag_Surface = 1
  ObsFlag_Sound = 1
  ObsFlag_Vwpw = 1
  obsVarList_Surface = temperature
  obsFileList_Surface = 20220125_1200, 20220125_1205
  obsFileList_Sound = rec_RTEMP_20220113000000_g.dat, rec_RTEMP_20220113030000_g.dat
  obsFileList_Vwpw = vwpw22011300.dat, vwpw22011306.dat
  mgStart = 6
  mgEnd = 6
  varWRFList = temp, height
/ 

&flog
  LogFileName = '/log/test'
/

&BMat
  ScaleParaX = 30000000.0
  ScaleParaY = 30000000.0
  ScaleParaZ = 3000.0
  ScaleParaT = 1.0
  NumEns = 6
  ensForm = state
  ensBKDiag = 1
  ensBEC = 1
  ensPath = '/EnsData_EnLoc'
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

C obs_thinning:
C obs_qcThreshold matches the model state variables, specify what range of values are acceptable
C obs_radius specifies the influence radius in horizontal, height, land, and time.
C thinning_threshold specifies those points where no thinning data if Gaussian is smaller than this.
C interpolation: 1 using unstructured grid interpolation scheme; 2 using lat-lon grid interpolation
