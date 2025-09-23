#!/bin/bash
rm -r build 
mkdir -p build
cmake -S . -B build -DCMAKE_CXX_COMPILER=mpicxx
cmake --build build -j


# run the binary 
mpirun -np 1 ./build/moments 
