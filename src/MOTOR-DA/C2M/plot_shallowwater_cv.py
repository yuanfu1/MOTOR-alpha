import xarray as xr
import matplotlib.pyplot as plt

# Load forward model data
forward_ds = xr.open_dataset('output_forward.nc')

# Plot vorticity at the final time step
vorticity = forward_ds['vorticity'].isel(time=-1)
vorticity.plot()
plt.title("Forward Model Vorticity")
plt.show()

# Load adjoint model data
adjoint_ds = xr.open_dataset('output_adjoint.nc')

# Plot adjoint vorticity at the final time step
adj_vorticity = adjoint_ds['adj_vorticity'].isel(time=-1)
adj_vorticity.plot()
plt.title("Adjoint Model Vorticity")
plt.show()
