# Install script for directory: /home/lcrippa/tools/DMFTtools

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RELEASE")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/include/")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/include" TYPE DIRECTORY FILES "/home/lcrippa/tools/DMFTtools/build_new_cluster/include/")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  execute_process(COMMAND "/opt2/ohpc/pub/utils/cmake/3.15.4/bin/cmake" -E remove /home/lcrippa/tools/DMFTtools/src/dmft_tools_version_current_sha1.inc)
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/lib/libdmft_tools.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/lib" TYPE STATIC_LIBRARY FILES "/home/lcrippa/tools/DMFTtools/build_new_cluster/libdmft_tools.a")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/etc/")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/etc" TYPE DIRECTORY FILES "/home/lcrippa/tools/DMFTtools/build_new_cluster/etc/")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/bin/dmft_tools_config_user.sh")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/bin" TYPE FILE PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_WRITE WORLD_READ WORLD_EXECUTE SETUID FILES "/home/lcrippa/tools/DMFTtools/build_new_cluster/etc/dmft_tools_config_user.sh")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/bin/dmft_tools_config_global.sh")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/bin" TYPE FILE PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_WRITE WORLD_READ WORLD_EXECUTE SETUID FILES "/home/lcrippa/tools/DMFTtools/build_new_cluster/etc/dmft_tools_config_global.sh")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/version")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0" TYPE FILE PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_WRITE WORLD_READ WORLD_EXECUTE SETUID FILES "/home/lcrippa/tools/DMFTtools/build_new_cluster/version")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/lcrippa/.pkgconfig.d/dmft_tools.pc")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/lcrippa/.pkgconfig.d" TYPE FILE PERMISSIONS OWNER_WRITE OWNER_READ OWNER_EXECUTE GROUP_WRITE GROUP_READ GROUP_EXECUTE WORLD_WRITE WORLD_READ WORLD_EXECUTE SETUID FILES "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/etc/dmft_tools.pc")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/lcrippa/.modules.d/")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/lcrippa/.modules.d" TYPE DIRECTORY FILES "/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/etc/modules/")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  MESSAGE(
"
ADD LIBRARY TO YOUR SYSTEM: 
Pick ONE method below  [add it in your bash profile to make it permanent, e.g. ~/.bashrc]:
[33mMethod 1: use the provided dmft_tools environment module[m:
   $ module use $HOME/.modules.d
   $ module load dmft_tools/

[33mMethod 2: source the config script[m:
   $ source /home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/bin/dmft_tools_config_user.sh

[33mMethod 3: use pkg-config with the provided dmft_tools.pc[m:
   $ export PKG_CONFIG_PATH=/home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/etc/:$PKG_CONFIG_PATH
   $ pkg-config --cflags --libs dmft_tools

[33mMethod ADMIN: Add this line to the system shell configuration file, e.g. /etc/bash.bashrc[m
   $ source /home/lcrippa/opt/dmft_tools/gnu/for_slurm/2.3.0/bindmft_tools_config_global.sh
")

endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/lcrippa/tools/DMFTtools/build_new_cluster/src/DMFT_CTRL_VARS/cmake_install.cmake")
  include("/home/lcrippa/tools/DMFTtools/build_new_cluster/src/DMFT_GF/cmake_install.cmake")
  include("/home/lcrippa/tools/DMFTtools/build_new_cluster/src/DMFT_GFIO/cmake_install.cmake")
  include("/home/lcrippa/tools/DMFTtools/build_new_cluster/src/DMFT_WEISS_FIELD/cmake_install.cmake")
  include("/home/lcrippa/tools/DMFTtools/build_new_cluster/src/DMFT_EKIN/cmake_install.cmake")
  include("/home/lcrippa/tools/DMFTtools/build_new_cluster/src/DMFT_TIGHT_BINDING/cmake_install.cmake")
  include("/home/lcrippa/tools/DMFTtools/build_new_cluster/src/DMFT_FFT/cmake_install.cmake")
  include("/home/lcrippa/tools/DMFTtools/build_new_cluster/src/DMFT_CONVERGENCE/cmake_install.cmake")
  include("/home/lcrippa/tools/DMFTtools/build_new_cluster/src/DMFT_MISC/cmake_install.cmake")
  include("/home/lcrippa/tools/DMFTtools/build_new_cluster/src/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/lcrippa/tools/DMFTtools/build_new_cluster/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
