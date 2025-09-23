#!/bin/bash
rm -rf build
mkdir -p build
cmake -S . -B build
cmake --build build -j

# 
./build/constant 
