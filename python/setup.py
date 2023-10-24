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
from Cython.Build import cythonize

setup(
    ext_modules = cythonize(
      Extension('EZmock', ['EZmock.pyx'],
        include_dirs=['../include/'],
        libraries=['EZmock'],
        library_dirs=['../lib/'],
        extra_compile_args=['-std=c99','-O3','-Wall','-DSINGLE_PREC']),
      compiler_directives = {'language_level': '3'}
    )
)
