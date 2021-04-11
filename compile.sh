#!/usr/bin/bash
# Include path to local libraries here if they are not in $PATH
# Needed are eigen3 and boost.
local_opts="-I /usr/include/eigen3"

library_includes="-I $(dirname $0)/third-party/hopscotch-map/include"
warn_flags="-Wall -Wextra -Wpedantic -Wshadow"
compile_flags="-std=c++17 -Ofast -g"

if [ "$#" -eq 1 ]; then
    echo "Looking for params in file $1"
    line=$(head -n 1 $1)
    vals=($line)
    echo "Found n = ${vals[0]}, r = ${vals[1]}, s = ${vals[2]}, t = ${vals[3]}. Starting compilation for executable ${vals[0]}-${vals[1]}-${vals[2]}-${vals[3]}.out"
    g++ $(dirname $0)/src/main.cpp -o ${vals[0]}-${vals[1]}-${vals[2]}-${vals[3]}.out $compile_flags $warn_flags $local_opts $library_includes\
      -DN_NFOLD=${vals[0]} -DR_NFOLD=${vals[1]} -DS_NFOLD=${vals[2]} -DT_NFOLD=${vals[3]}
    echo "Finished compilation with status code $?."
elif [ "$#" -eq 2 ]; then
    echo "Looking for params in file $1."
    line=$(head -n 1 $1)
    vals=($line)
    echo "Found n = ${vals[0]}, r = ${vals[1]}, s = ${vals[2]}, t = ${vals[3]}. Starting compilation"
    g++ $(dirname $0)/src/main.cpp -o $2 $compile_flags $warn_flags $local_opts $library_includes -DN_NFOLD=${vals[0]} -DR_NFOLD=${vals[1]} -DS_NFOLD=${vals[2]} -DT_NFOLD=${vals[3]}
    echo "Finished compilation with status code $?."
elif [ "$#" -eq 4 ]; then
    echo "Compiling for n = $1, r = $2, s = $3 and t = $4. Producing executable $1-$2-$3-$4.out}"
    g++ $(dirname $0)/src/main.cpp -o $1-$2-$3-$4.out $compile_flags $warn_flags $local_opts $library_includes -DN_NFOLD=$1 -DR_NFOLD=$2 -DS_NFOLD=$3 -DT_NFOLD=$4
    echo "Finished compilation with status code $?."
elif [ "$#" -eq 5 ]; then
    echo "Compiling for n = $1, r = $2, s = $3 and t = $4"
    g++ $(dirname $0)/src/main.cpp -o $5 $compile_flags $warn_flags $local_opts $library_includes -DN_NFOLD=$1 -DR_NFOLD=$2 -DS_NFOLD=$3 -DT_NFOLD=$4
    echo "Finished compilation with status code $?."
else
    echo "Usage like one of the following:"
    echo "    ./compile.sh input.in                  # Compilation for input 'input.in'"
    echo "    ./compile.sh input.in exec_name.out    # Produce 'exec_name.out' as executable"
    echo "    ./compile.sh 1 2 3 4                   # Compilation for n = 1, r = 2, s = 3, t = 4"
    echo "    ./compile.sh 1 2 3 4 exec_name.out     # Produce 'exec_name.out' as executable"
fi
