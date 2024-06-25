#!/bin/bash

rm Pade

mpicxx -I/usr/local/include/boost/ -o Pade pade.cpp -L/usr/local/lib -O3 -lboost_mpi -lboost_serialization -lboost_thread -lboost_system -lboost_chrono -pthread -lmpfr -lgmp

mpirun --oversubscribe -np 1 ./Pade

