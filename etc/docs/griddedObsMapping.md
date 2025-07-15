This document lays out the strategy for projecting a gridded observation data structure at a finer grid to a coarser one. 

1. we use the gridded obs indices to mask a local grid function with 1 if obs is available; otherwise 0.

   a. Construct a local grid mask function, but there is a question regarding the c_tc_idx, waiting for Zilong's reply.

2. Map all local mask grid to a global grid.

3. Construct a global gridded obs grid on the coarser grid, with missing value at no obs grid points. 

4. Construct a local gridded observation data structure at coarser level using the g_t_sp_idx at coarser level.

Note that the current MOTOR has the same values of sp_t_g_idx_toCoarser and sp_t_g_idx_toFiner and sp_t_g_idx, we may consider removing those redundants. 