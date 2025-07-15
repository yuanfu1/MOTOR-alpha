&icosTr_grid
  glvl = 4
  slvl = 2
/

&latlon_grid
  num_grid = 1025,1025
  domain_latlon = -5, 5, -5, 5
/

&analysis_para
  start_time = 2020, 1, 1, 6, 0, 0
  end_time = 2020, 1, 1, 7, 0, 0
/

&geometry
mpi_layout_g14 = 1, 1
mpi_layout_g58 = 2, 2
mpi_layout_g9L = 2, 2
time_steps_g14 = 1
time_steps_g58 = 1
time_steps_g9L = 1
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
  vLevel = 66
  varList = ua, va, T, p, q
/

&IO
  bk_model = GRAPES
  NMLFileName = '/namlist.input'
  ModelFileName = '/grapesinput'
/ 

&flog
  LogFileName = '/log/test'
/

C poissonSovler
C	solver:	     	multigrid scheme, VVV: w cycle
C	nCycle:		    number of multigrid cycles
C	nIterPre:		  number of iterations in pre-cycle
C	nIterPost:		number of iterations in post-cycle
C	nRelax:		    number of over relaxations
C	omegas:	    	relaxations coefficients
C   max_band      max band of poisson coefficent matrix 

C icosTr_grid
C glvl: icos triangle grid G level;
C slvl: Icos triangle grid start level;
C ordr: order of accuracy;

C latlon_grid
C domain_latlon specifies the (lat0,lat1) X (lon0,lon1)

C geometry
C mpi_layout:  layout of mpi

