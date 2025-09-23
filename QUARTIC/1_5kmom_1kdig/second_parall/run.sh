#!/bin/bash 

mkdir -p build
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j

mpirun -np 5 ./build/second 

