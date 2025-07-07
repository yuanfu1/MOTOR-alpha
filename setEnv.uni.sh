source ./pathEnv.sh
export HDF5_USE_FILE_LOCKING=FALSE
git config --global core.filemode false
git config core.hooksPath .githooks

# Set environment for some compilers
# source /public/home/simi/opt/intel/oneapi/setvars.sh 

# Set environment variables for each library and include path
export LIBRARY_HDF5="SearchByCMake"
export INCLUDE_HDF5="SearchByCMake"

export LIBRARY_NETCDF_FORTRAN="SearchByCMake"
export INCLUDE_NETCDF_FORTRAN="SearchByCMake"

export LIBRARY_NETCDF_C="SearchByCMake"
export INCLUDE_NETCDF_C="SearchByCMake"

export LIBRARY_RTTOV="$(pwd)/external/RTTOV/RTTOV4MOTOR/lib"
export INCLUDE_RTTOV="$(pwd)/external/RTTOV/RTTOV4MOTOR/include"
export MOD_RTTOV="$(pwd)/external/RTTOV/RTTOV4MOTOR/mod"

export LIBRARY_MPI="SearchByCMake"
export INCLUDE_MPI="SearchByCMake"

export LIBRARY_Z="SearchByCMake"

# Select the Slint path based on the platform and architecture
# PATH_SLINT="${ETC_DIR}/Slint/Linux/gnu"
# PATH_SLINT="${ETC_DIR}/Slint/Linux/ifort"
# PATH_SLINT="${ETC_DIR}/Slint/Linux/szhpc"
PATH_SLINT="${ETC_DIR}/Slint/Darwin/arm64/gnu"
# PATH_SLINT="${ETC_DIR}/Slint/Darwin/x86_64/gnu"

export LIBRARY_SLINT="${PATH_SLINT}/libs"
export INCLUDE_SLINT="${PATH_SLINT}/mods"

