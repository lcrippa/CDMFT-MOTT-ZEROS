
if(NOT EXISTS "/home/lcrippa/tools/DMFTtools/build_new_cluster/install_manifest.txt")
  MESSAGE(FATAL_ERROR "Cannot find install manifest: /home/lcrippa/tools/DMFTtools/build_new_cluster/install_manifest.txt")
endif()

FILE(READ "/home/lcrippa/tools/DMFTtools/build_new_cluster/install_manifest.txt" files)
STRING(REGEX REPLACE "\n" ";" files "${files}")
FOREACH(file ${files})
  MESSAGE(STATUS "Uninstalling $ENV{DESTDIR}${file}")
  IF(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
    EXEC_PROGRAM(
      "/opt2/ohpc/pub/utils/cmake/3.15.4/bin/cmake" ARGS "-E rm -rf \"$ENV{DESTDIR}${file}\""
      OUTPUT_VARIABLE rm_out
      RETURN_VALUE rm_retval
      )
    IF(NOT "${rm_retval}" STREQUAL 0)
      MESSAGE(FATAL_ERROR "Problem when removing $ENV{DESTDIR}${file}")
    ENDIF()
  ELSE(IS_SYMLINK "$ENV{DESTDIR}${file}" OR EXISTS "$ENV{DESTDIR}${file}")
    MESSAGE(STATUS "File $ENV{DESTDIR}${file} does not exist.")
  ENDIF()
ENDFOREACH()
