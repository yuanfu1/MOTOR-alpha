&icosTr_grid
  glvl = 4
  slvl = 2
/

&latlon_grid
  num_grid = 1025,1025
  domain_latlon = 18, 27, 107, 119
/

&analysis_para
  start_time = 2022, 6, 1, 0, 0, 0
  end_time = 2022, 6, 1, 1, 0, 0
  datetime = 2022060100

/

&geometry
mpi_layout_g14 = 2, 2
mpi_layout_g58 = 2, 2
mpi_layout_g9L = 2, 2
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
  vLevel = 10
  varList = uwnd,vwnd,temp,pres
/

&ObsFlag_satob
  logic_satob_qc = .true. 
  dqcuv_satob=3.0
/

&IO
  bk_model = GRAPES
  NMLFileName = '/namlist.input'
  ModelFileName = '/grapesinput'
  bakFileList = wrfout_d01_2022-01-13_00:00:00, wrfout_d01_2022-01-13_01:00:00
  ObsFlag_Surface = 1
  ObsFlag_satob = 1
  ObsFlag_Vwpw = 1
  obsFileList_Surface = 20220113_0000, 20220113_0005
  obsVarList_Surface = temperature
  obsFileList_Satob = FileList_Satob.txt
  obsFileList_Vwpw = vwpw22011300.dat, vwpw22011306.dat
  obsType = vwpw
  mgStart = 3
  mgEnd = 6
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

