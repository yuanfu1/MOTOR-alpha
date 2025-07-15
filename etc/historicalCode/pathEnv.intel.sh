# This shell set the essential path env to the predefined folders under this directory.

echo $(pwd)
source ./pathEnv.sh

# Add the reference to Slint lib.
if [[ $(uname -a) =~ "Darwin" ]]; then
    export INCLUDE_SLINT="${ETC_DIR}/Slint/Darwin/ifort/mods"
    export LIBRARY_PATH_SLINT="${ETC_DIR}/Slint/Darwin/ifort/libs"
elif [[ $(uname -a) =~ "Linux" ]]; then
    export INCLUDE_SLINT="${ETC_DIR}/Slint/Linux/ifort/mods"
    export LIBRARY_PATH_SLINT="${ETC_DIR}/Slint/Linux/ifort/libs"
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

# git configurations
git config --global core.filemode false
