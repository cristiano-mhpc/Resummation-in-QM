#!/bin/bash 

mkdir -p build
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j


# Run the executable 
mpirun -np 1 ./build/result

