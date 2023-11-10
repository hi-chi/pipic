#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from glob import glob
import re
import sys

import setuptools
from pybind11.setup_helpers import Pybind11Extension
from setuptools import find_packages, setup
from setuptools.command.build_ext import build_ext

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
        #opts.append('-shared')
        opts.append('-fPIC')
        opts.append('-fopenmp')
        opts.append('-lfftw3')
        for ext in self.extensions:
            ext.extra_compile_args = opts
            ext.extra_link_args = opts
        build_ext.build_extensions(self)


if sys.version_info < (3, 8, 0, 'final', 0):
    raise SystemExit('Python 3.8 or later is required!')

with open('README.md', encoding='utf-8') as fd:
    long_description = fd.read()

with open('pipic/__init__.py', encoding='utf-8') as fd:
    try:
        lines = ''
        for item in fd.readlines():
            item = item
            lines += item + '\n'
    except Exception as exc:
        raise Exception('Caught exception {}'.format(exc))


version = re.search("__version__ = '(.*)'", lines).group(1)
maintainer = re.search("__maintainer__ = '(.*)'", lines).group(1)
description = re.search("__description__ = '(.*)'", lines).group(1)
url = re.search("__url__ = '(.*)'", lines).group(1)
license = re.search("__license__ = '(.*)'", lines).group(1)

classifiers = [
    'Development Status :: 5 - Production/Stable',
    'Operating System :: OS Independent',
    'Programming Language :: Python :: 3 :: Only',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.8',
    'Programming Language :: Python :: 3.9',
    'Programming Language :: Python :: 3.10',
    'Programming Language :: Python :: 3.11',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: {}'.format(license),
    'Topic :: Scientific/Engineering :: Physics']


if __name__ == '__main__':

    setup(
        name='pipic',
        version=version,
        author='pipic developer group',
        description=description,
        long_description=long_description,
        ext_modules=[pipic_cpp_module],
        install_requires=['matplotlib',
                          'numpy',
                          'numba'],
        packages=find_packages(),
        python_requires='>=3.8',
        cmdclass={'build_ext': BuildExt},
        zip_safe=False,
        classifiers=classifiers,
        license=license,
        url=url,
    )
