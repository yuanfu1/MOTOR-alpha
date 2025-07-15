# Numerical Prediction Model development plan
###Jan 1, 2023:
We break down our model development plan at this moment into two effort, regional and global.

1. regional plan: extend the Z-grid model for terrain following full 3-D model and adopt the NCAR CCPP (Common Community Physical Package) to build a regional model. 

	#### requirement
	All development must strictly follow the OO design so that all will be easily switch to the global model effort.

	Time line:
	
	2023-03: Extend the Z-grid model to full 3-D terrain following system;
	2023-06: Adopt CCPP to the system;
	2023-09: Testing;
	2023-12: Initial parallel testing with CMA-GD.


2. global plan: Continue the Galewsky test case for a  shallow water equation. testing different initialization to stablize the 144 hour run.

	Time line:
	
	2023-03: complete various initialization scheme to get a stable Galewsky run.
	2023-06: merge with the regional full 3-D model development.
	2023-09: Test hydrostatic model run
	2023-12: build a global hydrostatic model with CCPP.

###Mar 15, 2023:
#### Status
1. Regional Z-grid model has been ported to MOTOR-PS under alpha branch.