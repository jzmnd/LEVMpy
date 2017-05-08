#! /usr/bin/env python
"""
outputfunc.py
Output functions for levmpy.py

Fortran source code by J R McDonald
Based on LEVM version 8.06

Created by Jeremy Smith on 2017-04-16
"""

import os
import sys

__author__ = "Jeremy Smith"
__version__ = "1.1"


def print_fit_info(expt):
    """Print fit information"""
    print "\nLEVM : COMPLEX NONLINEAR LEAST SQUARES IMMITTANCE DATA FITTING PROGRAM"
    print "       VERSION 8.06 - 2/05\n"
    print "{:s}\n".format(expt.alpha)
    print "  DATA ENTERED IN {:s}{:s} FORMAT TO BE USED IN {:s}{:s} FIT".format(expt.dinp, expt.pinp, expt.dfit, expt.pfit)
    print "  CIRCUIT MODEL : {:s}\n".format(expt.fun)
    print "  *****  FIT OF {:s} DATA  *****\n".format(expt.dattyp)

    print "  # OF DATA POINTS = {:d}   WEIGHT: IRCH = {:d}   # OF FREE PARAMETERS = {:d}".format(expt.md, expt.irch, expt.nfrei)
    print "  PRINTS EVERY {:d} ITERATIONS   MAX # ITERATIONS = {:d}   MAIN WT USES: {:s}".format(expt.nprint, expt.maxfev, expt.qq)
    print "  CELL CAPACITANCE = {:e}".format(expt.celcap)

    if expt.ixi == 0:
        if expt.irch == 0:
            print "  WEIGHTS (STANDARD DEVIATIONS) ARE READ IN"
        if expt.irch == 1:
            print "  UNIT WEIGHTING ASSIGNED TO EACH POINT"
        if expt.irch == 2:
            print "  WEIGHTS INVOLVE THE MAGNITUDES OF DATA OR FUNCTION VALUES RAISED TO THE POWER {:e}".format(expt.parameters[31])
        if expt.irch == 3:
            print "  MODULUS WEIGHTING: RESULTS RAISED TO THE POWER {:e}".format(expt.parameters[31])
        if expt.irch == 4:
            print "  #1174 SPINOLO FRA WEIGHTING"
        if expt.irch == 5:
            print "  #1250 & 1286 ORAZEM-AGARWAL FRA WEIGHTING"
        if expt.irch == 6:
            print "  #1250 & PAR 273 ORAZEM-AGARWAL FRA WEIGHTING"
    else:
        print "  WEIGHTS USE XI OR U**2 AND XI: BOTH MAY BE FREE PARAMETERS"

    print "\n  INITIAL PARAMETER GUESSES AND FIXED (0) OR FREE (1 OR 2) STATUS"
    for i in range(16):
        j = i + 16
        print "     P({:2d}) = {: e}   {:d}     P({:2d}) = {: e}   {:d}".format(i + 1, expt.parameters[i], expt.nfree[i],
                                                                                j + 1, expt.parameters[j], expt.nfree[j])
    if expt.n > 32:
        print "\n  THE FOLLOWING PARAMETERS ARE ALWAYS FIXED"
        for i in range(32, 36):
            j = i + 4
            print "     P({:2d}) = {: e}         P({:2d}) = {: e}".format(i + 1, expt.parameters[i], j + 1, expt.parameters[j])
    return

def print_outputs(expt):
    """Write output data"""
    print "\n  PARAMETER ESTIMATES"
    for i, n in enumerate(expt.ns1):
        print "     P({:2d}) = {: e}".format(n, expt.x[i])
    print "  NFEV =", expt.result[1]
    print "\n  OBSERVED VARIABLES:"
    print "  M  T   MEASURED      ESTIMATED     UNCERTANTY    RESIDUALS     RESID/MODEL"
    for i in range(expt.md):
        j = i + expt.md
        print "{:3d}  R  {: e} {: e} {: e} {: e} {: e}".format(i + 1, expt.y[i], expt.outputvals[i],
                                                               expt.r[i], expt.res[i], expt.resmod[i])
        print "{:3d}  I  {: e} {: e} {: e} {: e} {: e}".format(i + 1, expt.y[j], expt.outputvals[j],
                                                               expt.r[j], expt.res[j], expt.resmod[j])
    return
