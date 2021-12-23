src_files = [ 'workpackage2.cpp', 'Options.cpp' ]
deps = [ btllib_dep, threads_dep, openmp_dep, tcmalloc_dep, sdsl_dep, zlib_dep ]
inc_dirs = [ 'Common' ]
link_args = [ '-ltcmalloc' ]

executable('workpackage2',
    src_files,
    dependencies : deps,
    include_directories : inc_dirs,
    link_args : link_args
)