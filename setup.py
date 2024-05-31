#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from glob import glob
import sys
import os

import setuptools
from pybind11.setup_helpers import Pybind11Extension
from setuptools import find_packages, setup
from setuptools.command.build_ext import build_ext

# The first argument defines the name of the binary.
# It is customary to name such modules with a leading underscore.
pipic_cpp_module = Pybind11Extension(
    'pipic._pipic',
    sorted(glob('src/*.cpp')),
    language='c++')

# Build one module for each subfolder in extensions_source_directory.
# The modules are named with a leading underscore and placed in extensions_directory
extensions_source_directory = 'src/extensions/'
extensions_directory = 'pipic.extensions.'
extension_modules = []
for name in os.listdir(extensions_source_directory):
    if os.path.isdir(extensions_source_directory + name):
        cpp_module = Pybind11Extension(
            extensions_directory + '_' + name,
            sorted(glob(extensions_source_directory + name + '/*.cpp')),
            language='c++',
            include_dirs=['src/', extensions_source_directory + name + '/'])
        extension_modules.append(cpp_module)

# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689
def has_flag(compiler, flagname):
    """
    Returns a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile
    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        print("\033[33mNOTE: The following is just a compiler flag test.\033[0m")
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            print('\033[33mCompiler flag test failed! Flag will be ignored.\033[0m')
            return False
    print("\033[33mCompiler flag test completed.\033[0m")
    return True


def cpp_flag(compiler):
    """
    Returns the -std=c++[11/14] compiler flag.

    The c++14 is prefered over c++11 (when it is available).
    """
    if has_flag(compiler, '-std=c++14'):
        return '-std=c++14'
    elif has_flag(compiler, '-std=c++11'):
        return '-std=c++11'
    else:
        raise RuntimeError('Unsupported compiler -- at least C++11 support '
                           'is needed!')


class BuildExt(build_ext):
    """A custom build extension for adding compiler-specific options."""
    c_opts = {
        'msvc': ['/EHsc'],
        'unix': [],
    }

    if sys.platform == 'darwin':
        c_opts['unix'] += ['-mmacosx-version-min=10.7']

    def build_extensions(self):
        ct = self.compiler.compiler_type
        opts = self.c_opts.get(ct, [])
        if sys.platform == 'darwin':
            if has_flag(self.compiler, '-Xclang=-fopenmp'):
                # If using clang compiler on macOS
                opts.append('-stdlib=libc++')
                opts.append('-lomp')
                opts.append('-Xclang')  # Must come before the next -fopenmp (see ct=='unix' below)
            else:
                opts.append('-stdlib=libstdc++')
                # Fixes issue with xcode 15 linker 
                opts.append('-ld_classic')
        if ct == 'unix':
            opts.append('-fopenmp')
            opts.append("-DVERSION_INFO='{}'"
                        .format(self.distribution.get_version()))
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append('-fopenmp')
            opts.append("/DVERSION_INFO=\\'{}\\'"
                        .format(self.distribution.get_version()))
        opts.append('-O3')
        opts.append('-fPIC')
        opts.append('-lfftw3')

        for ext in self.extensions:
            ext.extra_compile_args = opts
            ext.extra_link_args = opts
        build_ext.build_extensions(self)


if sys.version_info < (3, 8, 0, 'final', 0):
    raise SystemExit('Python 3.8 or later is required!')


if __name__ == '__main__':
    setup(
        ext_modules=[pipic_cpp_module, *extension_modules],
        packages=find_packages(exclude=('tests',)),
        cmdclass={'build_ext': BuildExt}
    )
