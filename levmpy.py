#! /usr/bin/env python
"""
levmpy.py
Python module and wrapper for the LEVM complex nonlinear least squares program

Fortran source code by J R McDonald
Based on LEVM version 8.06

Created by Jeremy Smith on 2017-04-16
"""

import os
import sys
import LEVMpyFortran as lv
import numpy as np

__author__ = "Jeremy Smith"
__version__ = "1.0"


EPV = 8.85418782e-14

def dfloat(val):
    """Converts Fortran D string format to float"""
    return float(val.replace('D', 'E'))


class Experiment():
    """Class for single LEVM input file"""
    def __init__(self, infile, path=os.getcwd()):
        self.path = path
        self.infile = infile

        _readinfile()

    def _readinfile(self):
        """Reads LEVM input file"""
        with open(os.path.join(self.path, self.infile), 'r') as f:
            # Line 1 - description of run
            self.alpha = f.readline().strip()
            
            line2 = f.readline()
            line3 = f.readline()

            # Line 2 - data and fit types
            self.iopt = int(line2[:4])
            self.iacy = 0
            if self.iopt < 0:
                self.iacy = abs(self.iopt)
                self.iopt = 0

            self.dinp = line2[6]
            self.dfit = line2[7]
            self.pinp = line2[8]
            self.pfit = line2[9]
            self.ipf = ['R', 'P', 'D', 'L'].index(self.pfit)

            self.freeq = line2[10]
            self.neg = line2[11]
            self.fun = line2[12]
            self.celcap = dfloat(line2[14:23])
            self.dattyp = line2[23]
            self.ipar = int(line2[25:34])
            self.iorig = 0
            if self.ipar < 0:
                self.iorig = 1

            self.roe = dfloat(line2[35:44])
            self.ifp = int(line2[45:50])
            self.jfp = 1
            if self.ifp < 0:
                self.ifp = abs(self.ifp)
                self.jfp = 0
                
            self.ire = int(line2[51:56])

            # Line 3 - data, weighting, and fitting specifications
            self.md = int(line3[:5])
            self.imd = 1
            if self.md < 0:
                self.md = abs(self.md)
                self.imd = -1

            self.np = int(line3[6:10])
            self.inp = 0
            if self.np < 0:
                self.np = abs(self.np)
                self.inp = 1

            self.maxfev = int(line3[11:15])
            self.nprint = int(line3[16:20])
            self.irch = int(line3[21:25])
            if self.irch < 0:
                self.irch = abs(self.irch)
                self.iwt = 1
                self.qq = 'FUNC'
            else:
                self.iwt = 0
                if self.irch == 1:
                    self.qq = 'UNIT'
                else:
                    self.qq = 'DATA'

            self.mode = int(line3[26:30])
            self.icp = int(line3[31:35])
            self.iprint = int(line3[36:40])
            self.ipl = 0
            if self.iprint < 0:
                self.iprint = abs(self.iprint)
                self.ipl = self.iprint                

            self.igacc = int(line3[41:45])
            self.atemp = dfloat(line3[46:55])

            # Default CELCAP
            if self.celcap == 0:
                self.celcap = 1.0
            if self.imd == -1:
                self.celcap = EPV

            # Lines 4 and 5 - Set model parameters, and nfree/x/ns arrays
            plist = []
            for i in range(int(self.np / 5) + 1):
                line = f.readline().strip().split()
                for p in line:
                    plist.append(dfloat(p))

            self.parameters = np.array(plist, dtype=float)

            self.irg = 0
            if self.parameters[30] < 0:
                self.irg = int(round(abs(self.parameters[30])))
                self.parameters[30] = 0

            self.nfree = np.array(list(f.readline()), dtype=int)

            self.ixi = 0
            if (self.nfree[31] > 0) and (self.irch > 1):
                self.ixi = 1

            self.ns = np.where(self.nfree > 0)[0]   # note 0-start index
            self.x = self.parameters[self.ns]
            self.nfrei = len(self.x)
            self.nfix = self.np - self.nfrei

            # Read input data
            flist = []
            y1list = []
            y2list = []
            for i in range(self.md):
                line = f.readline().strip().split()
                flist.append(dfloat(line[1]))
                y1list.append(dfloat(line[2]))
                if self.neg == 'N':
                    y2list.append(-dfloat(line[3]))
                else:
                    y2list.append(dfloat(line[3]))

            self.freq = np.array(flist, dtype=float)
            self.y = np.array(y1list + y2list, dtype=float)

            if self.dattyp == 'C':
                self.jy = 1
                self.ky = 2 * self.md
            elif self.dattyp == 'R':
                self.jy = 1
                self.ky = self.md
                self.iopt = 0
            elif self.dattyp == 'I':
                self.jy = 1 + self.md
                self.ky = 2 * self.md
                self.iopt = 0

            if self.roe < 0:
                self.y = abs(self.roe) * self.y
                self.roe = 0

            # Perform data transformations if required
            if self.freeq == 'F':
                self.freq = 2 * np.pi * self.freq


            # Read fixed weighting parameters if they exist (IRCH=0)


        return


    def fit(self):

        # Set COMMON block variables

        # Run MAINCLC if MAXFEV > 3

        # If MAXFEV = 0 no fit calc new data
        # If MAXFEV = 1 no fit convert
        # If MAXFEV = 2 no fit calc new data without parameters with NFREE = 3

        # Write outputs

        return
