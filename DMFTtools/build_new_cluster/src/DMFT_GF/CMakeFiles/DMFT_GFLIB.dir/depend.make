# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.15

# Note that incremental build could trigger a call to cmake_copy_f90_mod on each re-build

src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/DMFT_GF.f90.o: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gf_common.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/DMFT_GF.f90.o: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gk_matsubara.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/DMFT_GF.f90.o: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gk_realaxis.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/DMFT_GF.f90.o: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gloc_matsubara.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/DMFT_GF.f90.o: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gloc_realaxis.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/DMFT_GF.f90.o.provides.build: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/dmft_gf.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/dmft_gf.mod.stamp: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/DMFT_GF.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/dmft_gf.mod src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/dmft_gf.mod.stamp GNU
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/DMFT_GF.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/DMFT_GF.f90.o.provides.build
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/build: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/DMFT_GF.f90.o.provides.build

src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GF_COMMON.f90.o: /home/lcrippa/opt/scifor/gnu/for_slurm/4.7.3/include/sf_arrays.mod
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GF_COMMON.f90.o: /home/lcrippa/opt/scifor/gnu/for_slurm/4.7.3/include/sf_constants.mod
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GF_COMMON.f90.o: /home/lcrippa/opt/scifor/gnu/for_slurm/4.7.3/include/sf_iotools.mod
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GF_COMMON.f90.o: /home/lcrippa/opt/scifor/gnu/for_slurm/4.7.3/include/sf_linalg.mod
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GF_COMMON.f90.o: /home/lcrippa/opt/scifor/gnu/for_slurm/4.7.3/include/sf_misc.mod
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GF_COMMON.f90.o: /home/lcrippa/opt/scifor/gnu/for_slurm/4.7.3/include/sf_mpi.mod
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GF_COMMON.f90.o: /home/lcrippa/opt/scifor/gnu/for_slurm/4.7.3/include/sf_timer.mod
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GF_COMMON.f90.o.provides.build: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gf_common.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gf_common.mod.stamp: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GF_COMMON.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/gf_common.mod src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gf_common.mod.stamp GNU
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GF_COMMON.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GF_COMMON.f90.o.provides.build
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/build: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GF_COMMON.f90.o.provides.build

src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GK_MATSUBARA.f90.o: include/dmft_ctrl_vars.mod
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GK_MATSUBARA.f90.o: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gf_common.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GK_MATSUBARA.f90.o.provides.build: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gk_matsubara.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gk_matsubara.mod.stamp: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GK_MATSUBARA.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/gk_matsubara.mod src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gk_matsubara.mod.stamp GNU
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GK_MATSUBARA.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GK_MATSUBARA.f90.o.provides.build
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/build: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GK_MATSUBARA.f90.o.provides.build

src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GK_REALAXIS.f90.o: include/dmft_ctrl_vars.mod
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GK_REALAXIS.f90.o: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gf_common.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GK_REALAXIS.f90.o.provides.build: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gk_realaxis.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gk_realaxis.mod.stamp: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GK_REALAXIS.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/gk_realaxis.mod src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gk_realaxis.mod.stamp GNU
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GK_REALAXIS.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GK_REALAXIS.f90.o.provides.build
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/build: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GK_REALAXIS.f90.o.provides.build

src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GLOC_MATSUBARA.f90.o: include/dmft_ctrl_vars.mod
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GLOC_MATSUBARA.f90.o: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gf_common.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GLOC_MATSUBARA.f90.o.provides.build: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gloc_matsubara.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gloc_matsubara.mod.stamp: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GLOC_MATSUBARA.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/gloc_matsubara.mod src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gloc_matsubara.mod.stamp GNU
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GLOC_MATSUBARA.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GLOC_MATSUBARA.f90.o.provides.build
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/build: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GLOC_MATSUBARA.f90.o.provides.build

src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GLOC_REALAXIS.f90.o: include/dmft_ctrl_vars.mod
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GLOC_REALAXIS.f90.o: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gf_common.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GLOC_REALAXIS.f90.o.provides.build: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gloc_realaxis.mod.stamp
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gloc_realaxis.mod.stamp: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GLOC_REALAXIS.f90.o
	$(CMAKE_COMMAND) -E cmake_copy_f90_mod include/gloc_realaxis.mod src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/gloc_realaxis.mod.stamp GNU
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GLOC_REALAXIS.f90.o.provides.build:
	$(CMAKE_COMMAND) -E touch src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GLOC_REALAXIS.f90.o.provides.build
src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/build: src/DMFT_GF/CMakeFiles/DMFT_GFLIB.dir/GLOC_REALAXIS.f90.o.provides.build