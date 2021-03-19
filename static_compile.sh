#!/usr/bin/bash
local_opts=-I third-party/hopscotch-map/include -I /usr/include/eigen3 -ldl -backtrace -DUSING_BOOST

if [ "$#" -eq 4 ]; then
    echo "Compiling for n = $1, r = $2, s = $3 and t = $4"
    g++ static_main.cpp -o static_a.out -std=c++17 $local_opts -Ofast -g -DN_NFOLD=$1 -DR_NFOLD=$2 -DS_NFOLD=$3 -DT_NFOLD=$4
elif [ "$#" -eq 1 ]; then
    echo "Looking for params in file $1."
    line=$(head -n 1 $1)
    vals=($line)
    echo "Found n = ${vals[0]}, r = ${vals[1]}, s = ${vals[2]}, t = ${vals[3]}. Starting compilation"
    g++ static_main.cpp -o static_a.out -std=c++17 -I /usr/include/eigen3 -I third-party/hopscotch-map/include/ -ldl -lbacktrace -Ofast -g -DN_NFOLD=${vals[0]} -DR_NFOLD=${vals[1]} -DS_NFOLD=${vals[2]} -DT_NFOLD=${vals[3]}
else
    echo "Please provide 1 param (file path) or 4 params (values for n, r, s and t)"
fi
