#! /usr/bin/env python
"""
levmpy.py
Python module and wrapper for the LEVM complex nonlinear least squares program

Fortran source code by J R McDonald

Created by Jeremy Smith on 2017-04-16
"""

import os
import sys

import LEVMpyFortran as lv
import numpy as np

__author__ = "Jeremy Smith"
__version__ = "1.0"


class Experiment():
    """Class for single LEVM input file"""
    def __init__(self, infile):

        # Read infile by line and set parameters and data


    def fit(self):

        # Set common variables

        # Run MAINCLC if MAXFEV > 3

        # If MAXFEV = 0 no fit calc new data
        # If MAXFEV = 1 no fit convert
        # If MAXFEV = 2 no fit calc new data without parameters with NFREE = 3

        # Write outputs
