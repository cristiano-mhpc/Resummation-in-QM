#!/bin/bash

rm result

mpicxx -I/usr/local/include/boost/ -o result result.cpp -L/usr/local/lib -lboost_mpi -lboost_serialization -lboost_thread -lboost_system -lboost_chrono -pthread -lmpfr -lgmp

mpirun --oversubscribe -np 1 ./result
