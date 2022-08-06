project('goldrush', 'cpp',
        version : '0.9.3',
        license : 'GPL-3',
        default_options : [ 'cpp_std=c++17', 'warning_level=3', 'werror=true' ])


threads_dep = dependency('threads')
openmp_dep = dependency('openmp')

compiler = meson.get_compiler('cpp')
divsufsort_dep = compiler.find_library('divsufsort')
divsufsort64_dep = compiler.find_library('divsufsort64')
btllib_dep = compiler.find_library('btllib')
sdsl_dep = compiler.find_library('sdsl')

subdir('goldrush_path')
subdir('bin')
subproject('goldrush-edit')

rootpath = meson.source_root()

clang_format = find_program('clang-format', required : false)
if clang_format.found()
  run_target('format', command: join_paths(rootpath, 'scripts/format'))
endif
