# Variational Cluster Approximation

An equilibrium Variational Cluster Approximation code.

The code is based on:  
* SciFortran [https://github.com/aamaricci/SciFortran]  
* DMFT_Tools [https://github.com/aamaricci/DMFTtools]

--

### REQUIRES

SciFortran
DMFTtools

### INSTALLATION
```
mkdir build
cd build
cmake .. -DBUILD_TYPE=(TESTING/DEBUG) -DEXE=(FILE IN ../drivers) 
make
```
### HOW IT WORKS

The user has to provide the hopping matrices for the cluster and the lattice (in the cluster basis). Optionally, a bath can be fed to the program (type "normal" of DMFT-ED). VCA_SOLVE returns the value of the Potthoff potential for the input the user provides. Minimization procedure or plots of Ω must be implemented by the user. Periodization must be implemented by the user (examples are provided in drivers).




