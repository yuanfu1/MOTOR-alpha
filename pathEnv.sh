echo $(pwd)
export CTEMPLATES="$(pwd)/ctemplates"
export COMMON_LIBS="$(pwd)/libs"
export COMMON_MODS="$(pwd)/mods"
export STATIC_DIR="$(pwd)/static"
export SRC_DIR="$(pwd)/src"
export GRID_DIR="$(pwd)/grid"
export BUILD_DIR="$(pwd)/build"
export BIN_DIR="$(pwd)/build/Debug"
export ETC_DIR="$(pwd)/etc"
export LOG_DIR="$(pwd)/log"
export CMAKETMP_DIR="$(pwd)/ctemplates"
export OUTPUT_DIR="$(pwd)/output"
export INPUT_DIR="$(pwd)/input"
export RTTOV_INCLUDE_DIR="$(pwd)/external/RTTOV/RTTOV4MOTOR/include"
export RTTOV_MOD_DIR="$(pwd)/external/RTTOV/RTTOV4MOTOR/mod"

# Clean the code if given the “clean” argument.
if [[ $1 =~ "clean" ]]; then
    rm -rf ${COMMON_MODS}
    # rm -rf ${COMMON_LIBS}
    rm -rf ${BUILD_DIR}
fi

# Build the essential directory for compilation and running.
if [ ! -d ${CTEMPLATES} ]; then
    mkdir ${CTEMPLATES}
    echo "./ctemplates has been built."
else
    echo "./ctemplates has been found!"
fi

if [ ! -d ${COMMON_LIBS} ]; then
    mkdir ${COMMON_LIBS}
    echo "./libs has been built."
else
    echo "./libs has been found!"
fi

if [ ! -d ${COMMON_MODS} ]; then
    mkdir ${COMMON_MODS}
    echo "./mods has been built."
else
    echo "./mods has been found!"
fi

if [ ! -d ${STATIC_DIR} ]; then
    mkdir ${STATIC_DIR}
    echo "./static has been built."
else
    echo "./static has been found!"
fi

if [ ! -d ${SRC_DIR} ]; then
    mkdir ${SRC_DIR}
    echo "./src has been built."
else
    echo "./src has been found!"
fi

if [ ! -d ${GRID_DIR} ]; then
    mkdir ${GRID_DIR}
    echo "./grid has been built."
else
    echo "./grid has been found!"
fi

if [ ! -d ${BUILD_DIR} ]; then
    mkdir ${BUILD_DIR}
    echo "${BUILD_DIR} has been built."
else
    echo "${BUILD_DIR} has been found!"
fi

if [ ! -d ${BIN_DIR} ]; then
    mkdir ${BIN_DIR}
    echo "${BIN_DIR} has been built."
else
    echo "${BIN_DIR} has been found!"
fi

if [ ! -d ${ETC_DIR} ]; then
    mkdir ${ETC_DIR}
    echo "${ETC_DIR} has been built."
else
    echo "${ETC_DIR} has been found!"
fi

if [ ! -d ${LOG_DIR} ]; then
    mkdir ${LOG_DIR}
    echo "${LOG_DIR} has been built."
else
    echo "${LOG_DIR} has been found!"
fi

if [ ! -d ${OUTPUT_DIR} ]; then
    mkdir ${OUTPUT_DIR}
    echo "${OUTPUT_DIR} has been built."
else
    echo "${OUTPUT_DIR} has been found!"
fi

if [ ! -d ${INPUT_DIR} ]; then
    mkdir ${INPUT_DIR}
    echo "${INPUT_DIR} has been built."
else
    echo "${INPUT_DIR} has been found!"
fi
