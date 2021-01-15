#!/bin/bash
cmake -Wno-deprecated -G "Ninja" -S . -B build
ninja -C build
