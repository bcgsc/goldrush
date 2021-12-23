project('goldrush', 'cpp',
        version : '0.0.1',
        license : 'GPL3',
        default_options : [ 'cpp_std=c++20', 'warning_level=2', 'werror=true' ])

btllib_dep = declare_dependency(include_directories : 'external/btllib/include')

threads_dep = dependency('threads')
openmp_dep = dependency('openmp')
tcmalloc_dep = meson.get_compiler('cpp').find_library('tcmalloc')
sdsl_dep = dependency('sdsl-lite')
zlib_dep = dependency('zlib')

subdir('workpackage2')