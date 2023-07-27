#!/bin/bash

mkdir -p "Build"

echo "building ChordLengthDistributions for CPU"
echo "--------------------------------------------------"
echo "compiling auxiliary.cpp"
g++ -w -fexceptions -O3 -std=c++11 -fopenmp  -c $PWD"/auxiliary.cpp" -o $PWD"/Build/auxiliary.o"
echo "compiling hdcommunication.cpp"
g++ -w -fexceptions -O3 -std=c++11 -fopenmp  -c $PWD"/hdcommunication.cpp" -o $PWD"/Build/hdcommunication.o"
echo "compiling chordlengthdistribution.cpp"
g++ -w -fexceptions -O3 -std=c++11 -fopenmp  -c $PWD"/chordlengthdistribution.cpp" -o $PWD"/Build/chordlengthdistribution.o"
echo "compiling neldermead.cpp"
g++ -w -fexceptions -O3 -std=c++11 -fopenmp  -c $PWD"/neldermead.cpp" -o $PWD"/Build/neldermead.o"
echo "compiling main.cpp"
g++ -w -fexceptions -O3 -std=c++11 -fopenmp  -c $PWD"/main.cpp" -o $PWD"/Build/main.o"
echo "linking ChordLengthDistributions"
g++  -o $PWD/ChordLengthDistributions $PWD/Build/main.o  $PWD/Build/auxiliary.o $PWD/Build/hdcommunication.o $PWD/Build/chordlengthdistribution.o  $PWD/Build/neldermead.o  -ltiff -lgomp
echo "--------------------------------------------------"

################################################################################################################################################################
################################################################################################################################################################
