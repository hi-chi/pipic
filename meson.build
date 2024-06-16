project(
    'pipic',
    'cpp',
    version: '1.1',
    default_options: [
        'cpp_std=c++11',
    ],
)

#
# Define directory variables
#
dir = 'pipic'
extensions_dir = 'src/extensions'


#
# Decide main module
#
py = import('python').find_installation(pure: false)
cpp = meson.get_compiler('cpp')


#
# Compiler flags
#
_global_cpp_args = cpp.get_supported_arguments(
    '-Wno-unused-but-set-variable',
    '-Wno-unused-functions',
    '-Wno-conversion',
    '-Wno-misleading-indentation',
    '-O3',
)
add_project_arguments(_global_cpp_args, language : 'cpp')

# macOS
if host_machine.system() == 'darwin'
    if cpp.has_link_argument('-Wl,-ld_classic')
        # New linker introduced in macOS 14 not working yet
        add_project_link_arguments('-Wl,-ld_classic', language : 'cpp')
    endif
    # if cpp.has_link_argument('-Wl,-dead_strip')
    #     # Allow linker to strip unused symbols
    #     add_project_link_arguments('-Wl,-dead_strip', language : 'cpp')
    # endif
endif


#
# Install python package files (found under pipic/)
#
install_subdir('pipic', install_dir: py.get_install_dir() / dir, strip_directory: true)


#
# Install pipic module
#
pybind11_dep = dependency('pybind11')
fftw3_dep = dependency('fftw3')

py.extension_module('_pipic',
    'src/pipic.cpp',
    subdir: dir,
    install: true,
    dependencies : [pybind11_dep, fftw3_dep],
    link_args: ['-fopenmp']
)

#
# Install extensions
#
interface_dep = declare_dependency(include_directories: include_directories('src'))
extension_deps = [pybind11_dep, fftw3_dep, interface_dep]

extensions = [
    'downsampler_gonoskov2022',
    'landau_lifshitz',
    'qed_gonoskov2015',
    'qed_volokitin2023',
    'x_converter_c',
    'x_reflector_c',

]
foreach e : extensions
    subdir(extensions_dir / e)
endforeach
