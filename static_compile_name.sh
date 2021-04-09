#!/usr/bin/bash
local_opts="-I $(dirname $0)/third-party/hopscotch-map/include -I /usr/include/eigen3"
warn_flags="-Wall -Wextra -Wpedantic -Wshadow"
compile_flags="-std=c++17 -Ofast -g"

if [ "$#" -eq 5 ]; then
    echo "Compiling for n = $1, r = $2, s = $3 and t = $4"
    g++ $(dirname $0)/src/main.cpp -o $5.out $compile_flags $warn_flags $local_opts -DN_NFOLD=$1 -DR_NFOLD=$2 -DS_NFOLD=$3 -DT_NFOLD=$4
elif [ "$#" -eq 2 ]; then
    echo "Looking for params in file $1."
    line=$(head -n 1 $1)
    vals=($line)
    echo "Found n = ${vals[0]}, r = ${vals[1]}, s = ${vals[2]}, t = ${vals[3]}. Starting compilation"
    g++ $(dirname $0)/src/main.cpp -o $2.out $compile_flags $warn_flags $local_opts -DN_NFOLD=${vals[0]} -DR_NFOLD=${vals[1]} -DS_NFOLD=${vals[2]} -DT_NFOLD=${vals[3]}
else
    echo "Please provide 2 params (file path and exec name) or 5 params (values for n, r, s, t and exec name)"
    echo "For example \"./compile.sh input.in name\" prodoces the executable \"name.out\""
fi
