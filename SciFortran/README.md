# SciFortran

This is a unitary collection of fortran modules and procedures for scientific calculations. The library aims to provide a simple and generic environment for any scientific or mathematic computations. The project is largely inspired by *SciPy* for Python and tries to closely follow its guidelines and naming convention. 

There are large areas that are still not covered.  
Anyone is welcome to contribute or to test the software. 

### Dependencies

* gfortran > 4.9x **OR** ifort  > 13.0
* cmake > 3.0.0    
* lapack  ( https://github.com/aamaricci/Lapack )   
* blas  ( https://github.com/aamaricci/Blas )   
* MPI ( https://github.com/open-mpi/ompi )  [optional, recommended]


If libraries are not available in your system, please use the provided links to install them. All libraries listed can be installed using `CMake` 



### Installation

Installation is now available using CMake. Experimental support for Intel MKL is provided but this still not universal and may end up in wrong linking directives. 
And from the repository directory (`cd scifor`) make a standard out-of-source CMake compilation:

`mkdir build`  
`cd build`  
`cmake ..`      (*)  
`make`     
`make install`   

(*) *In some cases CMake fails to find the MPI Fortran compiler, even if it is effectively installed and loaded into the system. An easy fix is to setup and export the `FC=mpif90` environment variable before invoking `cmake`.* 

Please follow the instructions on the screen to complete installation on your environment.  
The library can be loaded using one of the following, automatically generated, files :  

* pkg-config file in `~/.pkg-config.d/scifor.pc`  
* environment module file `~/.modules.d/scifor/<PLAT>/<VERSION>`  
* homebrew `bash` script `<PREFIX>/bin/configvars.sh`


The `CMake` compilation can be controlled using the following additional variables, default values between `< >`:   

* `-DPREFIX=prefix directory <~/opt/scifor/PLAT/VERSION>` 

* `-DUSE_MPI=<yes>/no`  

* `-DVERBOSE=yes/<no> `  

* `-DBUILD_TYPE=<RELEASE>/TESTING/DEBUG`  

### UNINSTALL

`Cmake` does not officially provide uninstall procedure in the generate Makefile. Yet, it is possible to construct one such uninstall mode in a simple way. SCIFOR provides a way to uninstall the files generated inside any out-of-source compilation calling: 
`make uninstall`  

