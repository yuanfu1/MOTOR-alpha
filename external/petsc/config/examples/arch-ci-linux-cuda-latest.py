#!/usr/bin/env python3

import os
petsc_hash_pkgs=os.path.join(os.getenv('HOME'),'petsc-hash-pkgs')

if __name__ == '__main__':
  import sys
  import os
  sys.path.insert(0, os.path.abspath('config'))
  import configure
  configure_options = [
    '--package-prefix-hash='+petsc_hash_pkgs,
    '--with-make-test-np=2',
    'COPTFLAGS=-g -O',
    'FOPTFLAGS=-g -O',
    'CXXOPTFLAGS=-g -O',
    '--with-precision=double',
    '--with-clanguage=c',
    '--with-mpi-dir=/software/mpich-430p2-cuda129',
    '--with-cuda-dir=/usr/local/cuda-12.9',
    '--download-hypre=1',
    '--with-strict-petscerrorcode',
  ]

  configure.petsc_configure(configure_options)
