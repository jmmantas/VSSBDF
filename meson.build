project('VSSBDF', 'cpp')
project_sources=['src/IVP_Solver.cpp', 'src/VSBDF/SBDF_Solver.cpp']
openblasdep=dependency('openblas')
LISdir_include='/usr/local/include'
LIS_lib='/usr/local/lib/liblis.a'



inc_dir=['include', 'include/IVPs', LISdir_include]

#LIS_inc = include_directories(LISdir/include)
#LIS_lib = static_library(...)
#LIS_dep = declare_dependency(link_with : my_lib,
#include_directories : my_inc)
  

executable('IVP_Solver', project_sources,
dependencies:[openblasdep],include_directories : inc_dir,objects:LIS_lib)
