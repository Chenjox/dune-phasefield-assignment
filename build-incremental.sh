#!/bin/bash

clang-format -i src/biw4-67-phasefield.cc
# okay, wir gehen ins Verzeichnis über uns
cd build-cmake
# führen dunecontrol aus mit den argumenten hier
make
cd ..