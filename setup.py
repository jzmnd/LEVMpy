#!/usr/bin/env python

from numpy.distutils.core import setup, Extension

fortranlib = Extension(name = 'fortranlib',
    extra_compile_args = ['-Wall -Wno-tabs -O3 -ffast-math'],
    sources = ['levmpy/fortranlib/levmpyfortran_edited.pyf',
               'levmpy/fortranlib/*.f']
    )

setup(name = 'levmpy',
    description = 'Python module and wrapper for the LEVM complex nonlinear least squares program'
    version = '1.1',
    author = 'Jeremy Smith',
    author_email = 'j.smith.03@cantab.net',
    url = 'http://jeremysmithscientist.com',
    py_modules = ['levmpy'],
    ext_modules = [fortranlib]
    )
