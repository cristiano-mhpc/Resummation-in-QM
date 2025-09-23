#!/bin/bash 
#
mkdir -p build
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j


#run the program
mpirun -np 2 ./build/auto
