# Function: Compile the external libraries for MOTOR.

# 判断命令行输入
if [ $# -eq 0 ]; then
    echo "Usage: ./compileExt.gnu.sh [compiler]"
    echo "compiler: gnu-XX or intel"
    exit
fi

# 读取命令行输入
compiler=$1
echo "Compiler: $compiler"
if [ $compiler == "gnu-13" ]; then
    echo "gfortran-13 has some bugs in the running of MOTOR, please use gfortran-12 instead."
    exit
elif [ $compiler == "gnu-12" ]; then
    export FC=gfortran-12
    export CC=gcc-12
    export CXX=g++-12
elif [ $compiler == "gnu-11" ]; then
    export FC=gfortran-11
    export CC=gcc-11
    export CXX=g++-11
elif [ $compiler == "gnu-10" ]; then
    export FC=gfortran-10
    export CC=gcc-10
    export CXX=g++-10
elif [ $compiler == "gnu-9" ]; then
    export FC=gfortran
    export CC=gcc
    export CXX=g++
elif [ $compiler == "intel" ]; then
    export FC=ifort
    export CC=icc
    export CXX=icpc
else
    echo "Invalid compiler: $compiler"
    exit
fi

MYPATH=$(pwd)

# Compile LAPACK
cd external/lapack
rm -rf build
mkdir build
cd build
cmake -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX ..
make -j4
cp ./lib/lib* ../../../libs
cd ../../../

# Compile fortran-yaml-cpp
cd external/fortran-yaml-cpp
rm -rf build
mkdir build
cd build
cmake -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX ..
make -j4
cp ./lib* ../../../libs
cp ./*.mod ../../../mods
cp ./yaml-cpp/lib* ../../../libs
cd ../../../

# Compile BLAS
cd external/blas-3.8
cp make_gnu.inc make.inc
make clean
make
cp ./lib*.a ../../libs
rm -f *.o
cd ../../

# Compile EISPACK
cd external/eispack
cp makefile.gnu makefile
make clean
make
cp ./lib*.a ../../libs
cd ../../

# Compile LINPACK
cd external/linpack
cp makefile.gnu makefile
make clean
make
cp ./lib*.a ../../libs
cd ../../

# Compile TAPENADE
cd $MYPATH
cd external/tapenade
mkdir build
cd build
cmake -DCMAKE_Fortran_COMPILER=$FC -DCMAKE_C_COMPILER=$CC -DCMAKE_CXX_COMPILER=$CXX ..
make clean
make
cp *.a ${MYPATH}/libs
cd $MYPATH

# Compile RTTOV
# Take ~10 min on Shiyan supercomputer.
# if [[ $(uname -a) =~ "Linux" ]]; then
#  source ./setEnv.sh # Enable the MOTOR environments on Shiyan
# fi
export motor_path=$(pwd)
cd external/RTTOV/RTTOV4MOTOR/src
../build/rttov_compile.gnu.sh
cd $MYPATH


MYPATH=$(pwd)
cd external/fortran_stdlib
export PATH=$MYPATH/external/fortran_stdlib:$PATH
mkdir build
cd build
cmake -DCMAKE_Fortran_COMPILER=$FC -DBUILD_TESTING=FALSE -DCMAKE_INSTALL_PREFIX=. ..
make clean
make
make install
cp lib/*.a $MYPATH/libs
find include -name "*.mod" -exec cp {} $MYPATH/mods \;
