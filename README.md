# CE-TDDFT
This is a real time TDDFT code for Quantum-Espresso.

## Features
* NMR shielding tensors, EFG tensors
* EPR g-tensor, hyperfine couplings
* Mössbauer
* Norm-conserving (g-tensor only), USPP and PAW
* LDA and GGA functionals
* isolated and periodic systems
* integration with Quantum-Environment (solvent effects)


## Authors and contributors
D. Ceresoli, X. Qian


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


### News
* Jul 17, 2017: ce-tddft is now compatible with QE-6.1

### TODO
- improve performance of NL-terms in circular dichroism
- restart/checkpointing
- in-memory evolution for molecules
- projection over AO's
- integrators: taylor-n, cranck-nichoson-2
- spin-orbit
- periodic crystals

