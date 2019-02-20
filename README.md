# CE-TDDFT
This is a real time TDDFT code for Quantum-Espresso.

### News
* Feb 20, 2019: **(experimental)** electron wavepacket
* Feb 18, 2019: save time dependent charge density
* Feb 17, 2019: Ehrenfest dynamics implemented
* Nov 26, 2018: ce-tddft is now compatible with QE-6.3
* Jul 17, 2017: ce-tddft is now compatible with QE-6.1


## Features
* Optical absorption spectrum of molecules
* Ehrenfest dynamics, norm-conserving only


## Authors and contributors
D. Ceresoli, X. Qian, A. Genova, A. Krisihtal, M. Pavanello


## Build instructions:
### From Quantum-Espresso distribution
not available yet!


### Stand-alone 
1. Configure and compile QE as usual, then:
2. ```git clone https://github.com/dceresoli/ce-tddft.git```
3. ```cd ce-tddft```
4. ```./configure --with-qe-source="QE folder containing make.inc"```
5. ```make```

### Source releases
Official source releases are found [here](https://github.com/dceresoli/ce-tddft/releases).

### TODO
- Ehrenfest dynamics with USPP and PAW pseudos
- arbitrary pulses
- improve performance of NL-terms in circular dichroism
- restart/checkpointing
- in-memory evolution for molecules
- projection over AO's
- integrators: taylor-n, cranck-nichoson-2
- spin-orbit
- periodic crystals

