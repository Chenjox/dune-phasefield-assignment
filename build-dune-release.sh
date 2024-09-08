#!/bin/bash

# okay, wir gehen ins Verzeichnis über uns
cd ..
# führen dunecontrol aus mit den argumenten hier
dunecontrol --opts=biw4-07/release.opts all
cd biw4-07


[ -f "compile_commands.json" ] && rm "compile_commands.json"

[ -f "build-cmake/compile_commands.json" ] && mv "build-cmake/compile_commands.json" "compile_commands.json"