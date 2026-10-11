# Oracle application for the electron code-to-code comparison
# (docs/validation.md, "Electron oracles: Nebula and Geant4 MicroElec").
#
# Our own code (MIT, like the rest of this tree). It links against a Geant4
# installation the user builds outside this tree; nothing of Geant4 is copied
# here. Written from the Geant4 Book For Application Developers (release
# 11.4): "How to Define the main() Program", "How to Specify Materials in the
# Detector", "How to Make an Executable Program" and "Physics List" (the
# G4EmParameters::AddMicroElec(region) interface).
#
# Build (outside this tree):
#   cmake -S validation/oracles/geant4_microelec -B /path/to/build \
#         -DGeant4_DIR=/path/to/geant4-install/lib/cmake/Geant4
#   cmake --build /path/to/build
# then point GEANT4_MICROELEC_BIN at /path/to/build/lindhard_microelec.
cmake_minimum_required(VERSION 3.16)
project(lindhard_microelec CXX)

find_package(Geant4 REQUIRED)
include(${Geant4_USE_FILE})

add_executable(lindhard_microelec lindhard_microelec.cc)
target_link_libraries(lindhard_microelec ${Geant4_LIBRARIES})
