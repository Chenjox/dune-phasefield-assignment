#!/bin/bash

# okay, wir gehen ins Verzeichnis über uns
cd build-cmake
# führen dunecontrol aus mit den argumenten hier
make build_tests
make test
cd ..