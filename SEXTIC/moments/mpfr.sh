#!/bin/bash

rm moments

mpicxx -I/usr/local/include/boost/ -o moments moments.cpp -O3 -L/usr/local/lib -lboost_mpi -lboost_serialization -lboost_thread -lboost_system -lboost_chrono -pthread -lmpfr -lgmp

mpirun -np 1 ./moments 
