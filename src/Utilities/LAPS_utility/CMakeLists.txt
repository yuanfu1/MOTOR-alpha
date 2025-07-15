cmake_minimum_required(VERSION 3.5)

# Change to your project:
# build my project:
# Note: this cmake builds the OO design of Multiscale Data Assimilation System:
PROJECT(laps_utility)

#add_subdirectory(./src/grid)
add_subdirectory(./src/lib)
add_subdirectory(./src/lib/modules)
add_subdirectory(./src/lib/mthermo)
