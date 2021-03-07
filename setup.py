#!/usr/bin/env python

from os import listdir
from numpy.distutils.core import setup, Extension

sourcefiles = ['levmpy/fortranlib/levmpyfortran_edited.pyf']
sourcefiles += ['levmpy/fortranlib/' + f for f in listdir('levmpy/fortranlib/') if f.endswith('.f')]

fortranlib = Extension(
    name='levmpy.fortranlib.LEVMpyFortran',
    extra_f77_compile_args=['-Wall', '-Wno-tabs', '-O3', '-ffast-math'],
    sources=sourcefiles,
)

setup(
    name='levmpy',
    description='Python module and wrapper for the LEVM complex nonlinear least squares program',
    version='1.1',
    author='Jeremy Smith',
    author_email='j.smith.03@cantab.net',
    url='https://github.com/jzmnd/LEVMpy',
    download_url='https://github.com/jzmnd/LEVMpy',
    packages=['levmpy', 'levmpy.fortranlib'],
    ext_modules=[fortranlib]
)
