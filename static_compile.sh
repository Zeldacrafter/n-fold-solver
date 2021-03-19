#!/usr/bin/bash
g++ static_main.cpp -std=c++17 -I /usr/include/eigen3 -I third-party/hopscotch-map/include/ -ldl -lbacktrace -Ofast -g -DN_NFOLD=$1 -DR_NFOLD=$2 -DS_NFOLD=$3 -DT_NFOLD=$3
