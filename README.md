# ARACNE-multicore
Python wrapper to run a multicore version of [ARACNe2](http://califano.c2b2.columbia.edu/aracne) with a [SingularityCE](https://github.com/sylabs/singularity/) container companion. 

## Pre-requisites
- singularity-ce >= 3.8.0

## Installation
1. Install [SingularityCE](https://github.com/sylabs/singularity/)
2. Clone this repository
   ```
   git clone https://github.com/ddiannae/ARACNE-multicore.git
   cd ARACNE-multicore/
   ```
4. Build the cointainer
   ```
   sudo singularity build aracne_multicore.sif aracne_multicore.def
   ```
6. Run the container with an appropiate script
   ```
   ./aracne_multicore.sif examples/aracne_matrix.py
   ```
   

