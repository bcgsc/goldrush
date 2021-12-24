project('goldrush', 'cpp',
        version : '0.0.1',
        license : 'MIT',
        default_options : [ 'cpp_std=c++17', 'warning_level=3', 'werror=true' ])

btllib_dep = declare_dependency(include_directories : 'external/btllib/include')

threads_dep = dependency('threads')
openmp_dep = dependency('openmp')

compiler = meson.get_compiler('cpp')
tcmalloc_dep = compiler.find_library('tcmalloc')
divsufsort_dep = compiler.find_library('divsufsort')
divsufsort64_dep = compiler.find_library('divsufsort64')
sdsl_dep = compiler.find_library('sdsl')
zlib_dep = compiler.find_library('z')

subdir('workpackage2')