cmake_minimum_required(VERSION 3.5)
PROJECT(TapenadeAux)

# Use Fortran compiler:
ENABLE_LANGUAGE( Fortran )
ENABLE_LANGUAGE( C )
ENABLE_LANGUAGE( CXX )

# Generate lib file:
ADD_LIBRARY(TapenadeAux STATIC adBuffer.F adStack.c)