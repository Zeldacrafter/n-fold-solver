#!/bin/bash

# No arguments provided. Build everything.
if [ $# -eq 0 ] ; then
    enable_src=true
    enable_test=true
    enable_bench=true
fi

# Handle command line arguments
while test $# -gt 0
do
    case "$1" in
        --clean) rm -rf build nFold_* ;;
        -c)      rm -rf build nFold_* ;;
        --src)   enable_src=true   ;;
        -s)      enable_src=true   ;;
        --test)  enable_test=true  ;;
        -t)      enable_test=true  ;;
        --bench) enable_bench=true ;;
        -b)      enable_bench=true ;;
        --debug) enable_debug=true ;;
        -d)      enable_debug=true ;;
    esac
    shift
done

# Construct defines for build command
args=""
if [ "$enable_src" = true ] ; then
    args="${args} -DBUILD_SRC=ON"
fi
if [ "$enable_test" = true ] ; then
    args="${args} -DBUILD_TEST=ON"
fi
if [ "$enable_bench" = true ] ; then
    args="${args} -DBUILD_BENCH=ON"
fi
if [ "$enable_debug" = true ] ; then
    args="${args} -DMY_DEBUG=ON"
fi

# Build with cmake and ninja
cmake -Wno-deprecated -GNinja -B build . $args
ninja -C build

# Create links to finished builds
if [ "$enable_src" = true ] ; then
    ln -sf build/src/nFold_main nFold_main
fi
if [ "$enable_test" = true ] ; then
    ln -sf build/test/nFold_tests nFold_tests
fi
if [ "$enable_bench" = true ] ; then
    ln -sf build/bench/nFold_bench nFold_bench
fi
