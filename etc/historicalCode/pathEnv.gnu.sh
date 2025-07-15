# This shell set the essential path env to the predefined folders under this directory.

echo $(pwd)
source ./pathEnv.sh

# Add the reference to Slint lib. 同时判断uname -m 为x86-64
if [[ $(uname -a) =~ "Darwin" ]]; then
    if [[ $(uname -m) =~ "x86_64" ]]; then
        export INCLUDE_SLINT="${ETC_DIR}/Slint/Darwin/x86_64/gnu/mods"
        export LIBRARY_PATH_SLINT="${ETC_DIR}/Slint/Darwin/x86_64/gnu/libs"
        echo "The architecture is x86_64."
    elif [[ $(uname -m) =~ "arm64" ]]; then
        export INCLUDE_SLINT="${ETC_DIR}/Slint/Darwin/arm64/gnu/mods"
        export LIBRARY_PATH_SLINT="${ETC_DIR}/Slint/Darwin/arm64/gnu/libs"
        echo "The architecture is arm64."
    fi
elif [[ $(uname -a) =~ "Linux" ]]; then
    export INCLUDE_SLINT="${ETC_DIR}/Slint/Linux/gfortran/mods"
    export LIBRARY_PATH_SLINT="${ETC_DIR}/Slint/Linux/gfortran/libs"
    echo "The platform is Linux."
else
    echo $a" is not supported currently."
fi

# This is conserved for compilation options on sz hpc.
if [[ $1 =~ "szhpc" ]]; then
    export INCLUDE_SLINT="${ETC_DIR}/Slint/Linux/szhpc/mods"
    export LIBRARY_PATH_SLINT="${ETC_DIR}/Slint/Linux/szhpc/libs"
fi

export INCLUDE_RTTOV="$(pwd)/external/RTTOV/RTTOV4MOTOR/include"
export MOD_RTTOV="$(pwd)/external/RTTOV/RTTOV4MOTOR/mod"
export LIBRARY_PATH_RTTOV="$(pwd)/external/RTTOV/RTTOV4MOTOR/lib"

export INCLUDE_HDF5="/usr/include/hdf5/serial"
export LIBRARY_PATH_HDF5="/usr/lib/x86_64-linux-gnu/hdf5/serial"

# Jilong's configuration.
# export INCLUDE_HDF5="/home/jilongchen/apps/hdf5/include"
# export LIBRARY_PATH_HDF5="/home/jilongchen/apps/hdf5/lib"


# git configurations
git config --global core.filemode false
