#!/bin/bash
cmake -Wno-deprecated -GNinja -B build .
ninja -C build
