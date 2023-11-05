file(REMOVE_RECURSE
  "libdmft_tools.a"
  "libdmft_tools.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang Fortran)
  include(CMakeFiles/dmft_tools.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
