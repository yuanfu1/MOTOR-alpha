# poissonSolver_Kp

poissonSolver_Kp is a multigrid solver for solving stream and potential functions with vorticity and divergence and Dirichlet boundary condition.

To use the poissonSolver_Kp, please use:

```fortran
! Also link the lib `poissonSolver` to the CMakeLists.txt
USE poissonSolver_m, only poissonSolver_t
...

CALL ps%initialize(configFile, geometry)   ! Initialize the possionSolver
...

CALL ps%PoissonSol(sgFinest%gLevel, vLevel, rights, solutn)   ! PS solver
...

CALL ps%destroy                            ! Destroy the mpdd
```      

`poissonSolver_Kp` depends on the class `geometry` and the a `.nl` config file:

```fortran
&poissonSovler
  solver = 'FMV'       
  nCycle = 3
  nIterPre = 3
  nIterPost = 6
  nRelax = 2
  max_band = 10
  omegas = 1.951627, 0.569688
/

C poissonSovler
C solver:          multigrid scheme, VVV: w cycle
C nCycle:          number of multigrid cycles
C nIterPre:        number of iterations in pre-cycle
C nIterPost:       number of iterations in post-cycle
C nRelax:          number of over relaxations
C omegas:          relaxations coefficients
C max_band         max band of poisson coefficent matrix 
```

## Main function: PoissonSol

```fortran
  SUBROUTINE PoissonSol(this, gLevel, vLevel, rights, solutn)
    INTEGER(i_kind), INTENT(IN) :: vLevel,                       
                                   gLevel
    REAL(r_kind), INTENT(INOUT) :: rights(vlevel, this%geometry%mg%sg(gLevel)%num_cell)
    REAL(r_kind), INTENT(OUT) :: solutn(vlevel, this%geometry%mg%sg(gLevel)%num_cell)

```

Here the `rights` and `solutn` are both distributive variables operated by `geometry`. The size of them are `sg%numiCell`. `solutn` is the initial value at the beginning, and the result finally. The boundary condition of the possion solver is included in the `rights`. 
