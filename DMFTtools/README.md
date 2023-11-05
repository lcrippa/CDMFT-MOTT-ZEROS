# DMFTtools
A collection of fortran modules and routines for DMFT and derivatives based on SciFortran

### Dependencies

* SciFortran (also available in this repo)

If libraries are not available in your system, please use the provided links to install them. All libraries listed can be installed using `CMake` 

### Installation

Installation is now available using CMake. Experimental support for Intel MKL is provided but this still not universal and may end up in wrong linking directives. 
To build this library, follow these steps:

`mkdir build`  
`cd build`  
`cmake ..` 
`make`     
`make install`   
