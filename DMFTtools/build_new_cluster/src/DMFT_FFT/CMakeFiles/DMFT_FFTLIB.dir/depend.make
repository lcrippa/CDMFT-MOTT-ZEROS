# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Note that incremental build could trigger a call to cmake_copy_f90_mod on each re-build

src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTAUX.f90.o: src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/dmft_fftgf.mod.stamp
src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTAUX.f90.o.provides.build: src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/dmft_fftaux.mod.stamp
src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/dmft_fftaux.mod.stamp: src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTAUX.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/dmft_fftaux.mod src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/dmft_fftaux.mod.stamp GNU
src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTAUX.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTAUX.f90.o.provides.build
src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/build: src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTAUX.f90.o.provides.build

src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTGF.f90.o: /home/lcrippa/opt/scifor/gnu/for_slurm/4.7.3/include/sf_arrays.mod
src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTGF.f90.o: /home/lcrippa/opt/scifor/gnu/for_slurm/4.7.3/include/sf_constants.mod
src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTGF.f90.o: /home/lcrippa/opt/scifor/gnu/for_slurm/4.7.3/include/sf_fft_fftpack.mod
src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTGF.f90.o: /home/lcrippa/opt/scifor/gnu/for_slurm/4.7.3/include/sf_interpolate.mod
src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTGF.f90.o: /home/lcrippa/opt/scifor/gnu/for_slurm/4.7.3/include/sf_linalg.mod
src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTGF.f90.o: /home/lcrippa/opt/scifor/gnu/for_slurm/4.7.3/include/sf_special.mod
src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTGF.f90.o.provides.build: src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/dmft_fftgf.mod.stamp
src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/dmft_fftgf.mod.stamp: src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTGF.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/dmft_fftgf.mod src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/dmft_fftgf.mod.stamp GNU
src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTGF.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTGF.f90.o.provides.build
src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/build: src/DMFT_FFT/CMakeFiles/DMFT_FFTLIB.dir/DMFT_FFTGF.f90.o.provides.build
