#!/bin/bash 
#
mkdir -p build
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j


cd build
mpirun -np 4 ./lu_mpi

