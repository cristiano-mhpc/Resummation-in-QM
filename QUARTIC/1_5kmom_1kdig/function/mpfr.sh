#!/bin/bash


rm  function

mpicxx -I/usr/local/include/boost/ -o function function.cpp -L/usr/local/lib -lboost_mpi -lboost_serialization -lboost_thread -lboost_system -lboost_chrono -pthread -lmpfr -lgmp

mpirun --oversubscribe -np 1 ./function
