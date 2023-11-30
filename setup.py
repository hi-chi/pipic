#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from glob import glob
import sys

import setuptools
from pybind11.setup_helpers import Pybind11Extension
from setuptools import find_packages, setup
from setuptools.command.build_ext import build_ext

# The first argument defines the name of the binary.
# It is customary to name such modules with a leading underscore.
pipic_cpp_module = Pybind11Extension(
    '_pipic',
    sorted(glob('src/*.cpp')),
    language='c++')


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
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            print('Fails!\n')
            return False
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
            if has_flag(self.compiler, '-stdlib=libc++'):
                pass# opts.append('-stdlib=libc++')
        if ct == 'unix':
            opts.append("-DVERSION_INFO='{}'"
                        .format(self.distribution.get_version()))
            opts.append(cpp_flag(self.compiler))
            if has_flag(self.compiler, '-fvisibility=hidden'):
                opts.append('-fvisibility=hidden')
        elif ct == 'msvc':
            opts.append("/DVERSION_INFO=\\'{}\\'"
                        .format(self.distribution.get_version()))
        opts.append('-O3')
        opts.append('-fPIC')
        opts.append('-fopenmp')
        opts.append('-lfftw3')
        for ext in self.extensions:
            ext.extra_compile_args = opts
            ext.extra_link_args = opts
        build_ext.build_extensions(self)


if sys.version_info < (3, 8, 0, 'final', 0):
    raise SystemExit('Python 3.8 or later is required!')


if __name__ == '__main__':
    setup(
        ext_modules=[pipic_cpp_module],
        package_data={"pipic": ["__init__.pyi"]},
        packages=find_packages(),
        cmdclass={'build_ext': BuildExt},
        zip_safe=False,
    )
