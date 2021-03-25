#!/usr/bin/bash
local_opts="-I $(dirname $0)/third-party/hopscotch-map/include -I /usr/include/eigen3"
compile_flags="-std=c++17 -Ofast -g" # -fsanitize=address" # -fsanitize=undefined"

if [ "$#" -eq 1 ]; then
    g++ $(dirname $0)/src/main.cpp -o $1.out $compile_flags $local_opts
else
    echo "Please provide exactly one parameter: The name of the resulting executable"
fi
