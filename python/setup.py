#!/usr/bin/env python
#
# setup.py: this file is part of the EZmock python package.
#
# EZmock: Effective Zel'dovich approximation mock generator.
#
# Github repository:
#       https://github.com/cheng-zhao/EZmock
#
# Copyright (c) 2023 Cheng Zhao <zhaocheng03@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize
import os, shutil

package_name = 'EZmock'
package_version = '1.0.0'
lib_path = os.path.abspath('../build')  # build directory of libEZmock

# Copy libEZmock.so before building
class CustomBuildExt(build_ext):
  def run(self):
    lib_src = os.path.join(lib_path, 'lib', 'libEZmock.so')
    lib_dst = os.path.join(package_name)
    shutil.copy(lib_src, lib_dst)
    super().run()

ext_modules = [
  Extension(
    f'{package_name}.{package_name}',
    sources = [os.path.join(package_name, 'EZmock.pyx')],
    include_dirs = [os.path.join(lib_path, 'include')],
    libraries = ['EZmock'],
    library_dirs = [os.path.join(lib_path, 'lib')],
    runtime_library_dirs = ['$ORIGIN'],
    extra_compile_args = ['-std=c99','-O3','-Wall','-DSINGLE_PREC']
  ),
]

setup(
  name = package_name,
  version = package_version,
  packages = [package_name],
  package_data = {package_name: ['libEZmock.so']},
  ext_modules = cythonize(
   ext_modules,
   compiler_directives = {'language_level': '3'}
  ),
  cmdclass = {'build_ext': CustomBuildExt}
)
