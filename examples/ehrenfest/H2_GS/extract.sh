#!/bin/bash

out=H2-ehrenfest.out
grep temperature $out | awk '{print $3}' | nl > TEMPERATURE.dat
grep const H2-ehrenfest.out | awk '{print $6}' | nl >ENERGY.dat
awk '/ATOMIC_POSITIONS/ {getline; H1x=$2; H1y=$3; H1z=$4; getline; H2x=$2; H2y=$3; H2z=$4; print H1x*10, H2x*10}' $out | nl >HYDROGEN.dat

