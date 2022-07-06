src_files = [ 'goldrush_path.cpp', 'Options.cpp', 'read_hashing.cpp', 'opt.cpp', 'spaced_seeds.cpp', 'calc_phred_average.cpp' ]
deps = [ btllib_dep, threads_dep, openmp_dep, tcmalloc_dep, divsufsort_dep, divsufsort64_dep, sdsl_dep, btllib_dep ]
inc_dirs = [ 'Common' ]
link_args = [ '-ltcmalloc', '-lsdsl', '-ldivsufsort', '-ldivsufsort64', '-lbtllib' ]

executable('goldrush-path',
    src_files,
    dependencies : deps,
    include_directories : inc_dirs,
    link_args : link_args,
    install : true
)
