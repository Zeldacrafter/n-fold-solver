#!/bin/bash

: '
# No arguments provided. Build everything.
if [ $# -eq 0 ] ; then
    enable_test=true
    enable_bench=true
fi
'
if [ $# -eq 0 ] ; then
    echo "No arguments provided."
    echo "Possible:"
    echo "   --test"
    echo "   --bench"
    echo "   --clean"
    echo "   --debug"
    exit -1
fi

# Handle command line arguments
while test $# -gt 0
do
    case "$1" in
        --clean)  echo "Cleaning"
                  rm -rf nFold_*
                  rm -rf build       ;;
        -c)       echo "Cleaning"
                  rm -rf nFold_*
                  rm -rf build       ;;
        --test)   enable_test=true   ;;
        -t)       enable_test=true   ;;
        --bench)  enable_bench=true  ;;
        -b)       enable_bench=true  ;;
        --debug)  enable_debug=true  ;;
        -d)       enable_debug=true  ;;
    esac
    shift
done

# Construct defines for build command
args=""
if [ "$enable_test" = true ] ; then
    args="${args} -DBUILD_TEST=ON"
fi
if [ "$enable_bench" = true ] ; then
    args="${args} -DBUILD_BENCH=ON"
fi
if [ "$enable_debug" = true ] ; then
    args="${args} -DMY_DEBUG=ON"
fi

if [ -z "$args" ] ; then
    echo "No module to build provided"
    exit 1
fi

# Build with cmake and ninja
echo cmake -Wno-deprecated -GNinja -B build . $args
cmake -Wno-deprecated -GNinja -B build . $args
echo ninja -C build
ninja -C build

ln -sf build/nFold_main nFold_main

# Create links to finished builds
if [ "$enable_test" = true ] ; then
    ln -sf build/test/nFold_test nFold_test
fi
if [ "$enable_bench" = true ] ; then
    ln -sf build/bench/nFold_bench nFold_bench
fi
