src_files = [ 'workpackage2.cpp', 'Options.cpp', 'read_hashing.cpp', 'opt.cpp', 'spaced_seeds.cpp' ]
deps = [ btllib_dep, threads_dep, openmp_dep, tcmalloc_dep, divsufsort_dep, divsufsort64_dep, sdsl_dep ]
inc_dirs = [ 'Common' ]
link_args = [ '-ltcmalloc', '-lsdsl', '-ldivsufsort', '-ldivsufsort64' ]

executable('workpackage2',
    src_files,
    dependencies : deps,
    include_directories : inc_dirs,
    link_args : link_args,
    install : true
)