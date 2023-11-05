# Cluster-DMFT Lanczos Exact Diagonalization

Lanczos based solver for the **Cluster** Dynamical Mean-Field Theory, as used for the numerical calculations in [https://arxiv.org/abs/2301.05588]

This folder contains the CDMFT code proper, relying on the following libraries (also contained in this repo):  

* SciFortran 

* DMFT_Tools

The code structure is as follow:  

* The set of modules compile into a top layer named `CDMFT_ED.f90`  
* The actual implementation of the DMFT equations is case by case performed in a driver program, usually placed in the directory `drivers`. 
* In the driver code the user must includes the `CDMFT_ED` module and call the necessary procedures to solve the DMFT equations.

The code relies on the cmake build environment. To build a driver, assuming SciFortran and DMFTtools are built and accessible, follow the steps

```
mkdir build
cd build
cmake .. -DEXE=NAME_OF_DRIVER -DBUILD_TYPE=TESTING
make
```
