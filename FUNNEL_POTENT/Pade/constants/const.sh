#!/bin/bash
#prep
#Author:Christian
rm constant constant.o 

g++ -I /usr/include/eigen3 -c constant.cpp

g++ constant.o -o constant -O3 -lmpfr -lgmp 

./constant


