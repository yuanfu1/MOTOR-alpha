# This script is built for finding the SLINT library automatically.

find_path(SLINT_Fortran_INCLUDE_DIR slint.mod $ENV{INCLUDE_SLINT})
find_library(SLINT_LIBRARY NAMES Slint PATHS $ENV{LIBRARY_SLINT})

if (SLINT_Fortran_INCLUDE_DIR AND SLINT_LIBRARY)
    set(SLINT_FOUND TRUE)
endif (SLINT_Fortran_INCLUDE_DIR AND SLINT_LIBRARY)
