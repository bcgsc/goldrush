project('goldrush', 'cpp',
        version : '1.2.0',
        license : 'GPL-3',
        default_options : [ 'cpp_std=c++17', 'warning_level=3', 'werror=true' ])

threads_dep = dependency('threads')
openmp_dep = dependency('openmp')

compiler = meson.get_compiler('cpp')
tcmalloc_dep = compiler.find_library('tcmalloc')
divsufsort_dep = compiler.find_library('divsufsort')
divsufsort64_dep = compiler.find_library('divsufsort64')
sdsl_dep = compiler.find_library('sdsl')
btllib_dep = compiler.find_library('btllib')

subdir('goldrush_path')
subdir('bin')
subproject('goldpolish')

rootpath = meson.source_root()

clang_format = find_program('clang-format', required : false)
if clang_format.found()
  run_target('format', command: join_paths(rootpath, 'scripts/format'))
endif
