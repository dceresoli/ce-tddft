#!/bin/bash

pw="mpirun -np 8 $HOME/Codes/qe-6.8/bin/pw.x"
tddft="mpirun -np 8 ../../bin/tddft.x"

molecule=C6H6

$pw <$molecule-scf.in > $molecule-scf.out

for edir in x y z
do
  $tddft <$molecule-tddft_$edir.in > $molecule-tddft_$edir.out
  grep ^DIP $molecule-tddft_$edir.out >dip_$edir.dat
done



