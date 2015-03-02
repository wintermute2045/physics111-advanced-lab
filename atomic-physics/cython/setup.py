#!/usr/bin/env python2.7

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules=[Extension('gsl_wlsq',
                       ['gsl_wlsq.pyx'],
                       libraries=['gsl', 'gslcblas'])]

setup(
    cmdclass = {"build_ext": build_ext},
    ext_modules = ext_modules
)
