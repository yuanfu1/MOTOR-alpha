# Global Z-grid plotting
On my new laptop as Mar. 2023, the Z-grid model output plotting program is under
/Users/xiey/developments/models/plottings/build

1. GZM model outputs icosahedral grid value in the following form in binary: 
  
  WRITE(10) gzm%grid%num_cell
  
  WRITE(10) gzm%grid%cell_cntr
  
  WRITE(10) r(1,:)-aj(1,:)
2. Under current directory, set a build directoy to compile the src directory.
3. Symbolic link the output file to model.dat
3. 