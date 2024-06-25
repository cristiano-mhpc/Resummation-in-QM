#!/bin/bash

rm fourth

mpicxx -I/usr/local/include/boost/ -o fourth fourth.cpp -O3 -L/usr/local/lib -lboost_mpi -lboost_serialization -lboost_thread -lboost_system -lboost_chrono -pthread -lmpfr -lgmp

mpirun --oversubscribe -np 4 ./fourth
