MODULE mem_grid

  REAL, ALLOCATABLE, DIMENSION(:, :) :: lat, lon, topo, ldf
  INTEGER, ALLOCATABLE, DIMENSION(:, :) :: IstartIend
  INTEGER, ALLOCATABLE, DIMENSION(:)   :: recvcounts, displs
  INTEGER                              :: nPEs, rank

END MODULE
