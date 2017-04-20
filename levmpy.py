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
import levmpyfortran.LEVMpyFortran as lv
import numpy as np

__author__ = "Jeremy Smith"
__version__ = "1.0"


EPV = 8.85418782e-14
NPAFR = len(lv.cm2.ns)
NTOT = len(lv.cm2.nfree)
NPTS = len(lv.cm1.freq)
NPT2 = len(lv.cm2.y)


def dfloat(val):
    """Converts Fortran D string format to float"""
    return float(val.replace('D', 'E'))


def convert_to_rect(mod, theta, t='P'):
    """Converts polar to rectangular data"""
    if t == 'L':
        mod = np.exp(mod)
    if t == 'D':
        theta = np.pi * theta / 180
    rx = mod * np.cos(theta)
    ry = mod * np.sin(theta)
    return rx, ry


def convert_to_polar(rx, ry, t='P'):
    """Converts rectangular to polar data"""
    mod = np.sqrt(rx**2 + ry**2)
    theta = np.arctan(ry / rx)
    if t == 'L':
        mod = np.log(mod)
    if t == 'D':
        theta = 180 * theta / np.pi
    return mod, theta


def convert_to_inverse(y1, y2):
    """Converts to inverse of input data"""
    norm = y1**2 + y2**2
    y1 = y1 / norm
    y2 = -y2 / norm
    return y1, y2


def convert_m_to_z(y1, y2, fq, clcp=1.0):
    """Converts modulus to impedance"""
    r1 = y2 / (fq * clcp)
    r2 = -y1 / (fq * clcp)
    return r1, r2


def convert_z_to_m(y1, y2, fq, clcp=1.0):
    """Converts impedance to modulus"""
    r1 = -y2 * fq * clcp
    r2 = y1 * fq * clcp
    return r1, r2


class Experiment():
    """Class for single LEVM input file"""
    def __init__(self, infile, path=os.getcwd()):
        self.path = path
        self.infile = infile

        self._readinfile()

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
            if self.iacy == 0:
                self.ftol = 1e-30
                self.xtol = 1e-48
            else:
                self.ftol = 10**(self.iacy)
                self.xtol = self.ftol

            self.dinp = line2[6]
            self.dfit = line2[7]
            self.pinp = line2[8]
            self.pfit = line2[9]
            self.ipf = ['R', 'P', 'D', 'L'].index(self.pfit)

            self.freqtyp = line2[10]
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
            self.iopr = 0
            if self.ire < -10:
                self.iopr = 1

            # Line 3 - data, weighting, and fitting specifications
            self.md = int(line3[:5])
            self.imd = 1
            if self.md < 0:
                self.md = abs(self.md)
                self.imd = -1

            self.n = int(line3[6:10])
            self.inp = 0
            self.izr = 0
            if self.n < 0:
                self.n = abs(self.n)
                self.inp = 1
                if (self.dattyp == 'C') and (self.dinp == 'Z') and (self.pinp == 'R') and (self.dfit == 'Z') and (self.pfit == 'R'):
                    self.izr = 1

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

            # Set default CELCAP
            if self.celcap == 0:
                self.celcap = 1.0
            if self.imd == -1:
                self.celcap = EPV

            # Line 4 - set model parameters
            plist = []
            for i in range(int(self.n / 5)):
                line = f.readline().strip().split()
                for p in line:
                    plist.append(dfloat(p))

            self.parameters = np.array(plist, dtype=float)

            self.irg = 0
            if self.parameters[30] < 0:
                self.irg = int(round(abs(self.parameters[30])))
                self.parameters[30] = 0

            # Line 5 - set nfree/x/ns arrays
            nflist = list(f.readline().strip())
            self.nfree = np.array(nflist, dtype=int)

            self.ixi = 0
            if (self.nfree[31] > 0) and (self.irch > 1):
                self.ixi = 1

            self.ns = np.where(self.nfree > 0)[0]   # 0-start index
            self.ns1 = self.ns + 1                  # 1-start index
            self.x = self.parameters[self.ns]
            self.nfrei = len(self.x)
            self.nfixi = self.n - self.nfrei

            self.ine = 0
            if (self.fun == 'R') or (self.fun == 'K'):
                if self.nfree[28] > 0:
                    self.ine += 1
                if self.nfree[29] > 0:
                    self.ine += 1

            self.ins = 0
            if (self.nfree[30] == 0) and (self.nfree[31] != 0):
                self.ins = 1
            if (self.nfree[30] != 0) and (self.nfree[31] != 0):
                self.ins = 2

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
            self.y1 = np.array(y1list, dtype=float)
            self.y2 = np.array(y2list, dtype=float)

            # Set start and end index
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

            # Scale by A/L if ROE<0
            if self.roe < 0:
                self.y1 = abs(self.roe) * self.y1
                self.y2 = abs(self.roe) * self.y2
                self.roe = 0

            # Perform data transformations if required
            if self.freqtyp == 'F':
                self.omega = 2 * np.pi * self.freq
            else:
                self.omega = self.freq

            if not ((self.dinp == self.dfit) and (self.pinp == self.pfit)):
                if self.pinp != 'R':
                    self.y1, self.y2 = convert_to_rect(self.y1, self.y2, t=self.pinp)
                if (self.dinp == 'Y') or (self.dinp == 'E'):
                    self.y1, self.y2 = convert_to_inverse(self.y1, self.y2)
                if (self.dinp == 'E') or (self.dinp == 'M'):
                    self.y1, self.y2 = convert_m_to_z(self.y1, self.y2, self.omega, clcp=self.celcap)

                if (self.dfit == 'E') or (self.dfit == 'M'):
                    self.y1, self.y2 = convert_z_to_m(self.y1, self.y2, self.omega, clcp=self.celcap)
                if (self.dfit == 'Y') or (self.dfit == 'E'):
                    self.y1, self.y2 = convert_to_inverse(self.y1, self.y2)
                if self.pfit != 'R':
                    self.y1, self.y2 = convert_to_polar(self.y1, self.y2, t=self.pfit)

            self.y = np.hstack((self.y1, self.y2))

            # Read fixed weighting parameters if they exist (IRCH=0)


        return


    def fit(self):
        """Performs the CNLS fit"""
        # Set all COMMON block variables
        # Resize arrays
        self.omega.resize(NPTS)
        self.y.resize(NPT2)
        self.parameters.resize(NTOT)
        self.ns1.resize(NPAFR)
        self.nfree.resize(NPAFR)
        # CM1
        lv.cm1.freq = self.omega
        lv.cm1.m = self.md
        lv.cm1.dattyp = self.dattyp
        # CM2
        lv.cm2.y = self.y
        #lv.cm2.r = 0
        lv.cm2.p = self.parameters
        lv.cm2.drss = 1.0e-8
        lv.cm2.roe = self.roe
        lv.cm2.ns = self.ns1
        lv.cm2.nfree = self.nfree
        lv.cm2.np = self.n
        lv.cm2.irch = self.irch
        lv.cm2.ixi = self.ixi
        lv.cm2.dattyq = self.dattyp
        # CM3
        lv.cm3.celcap = self.celcap
        lv.cm3.fun = self.fun
        lv.cm3.dfit = self.dfit
        lv.cm3.pfit = self.pfit
        # CM10
        lv.cm10.epsg = 10**(-abs(self.igacc))
        lv.cm10.izr = self.izr
        # CM11
        lv.cm11.mqy = self.nfrei - self.ins
        # CM12
        lv.cm12.clcap = self.celcap
        lv.cm12.atemp = self.atemp
        lv.cm12.maxfev = self.maxfev
        lv.cm12.mde = self.mode
        # CM16
        lv.cm16.iop = self.iopt
        lv.cm16.iorig = self.iorig
        lv.cm16.jy = self.jy
        lv.cm16.iprint = self.iprint
        lv.cm16.ldfjac = self.ky
        lv.cm16.mode = self.mode
        lv.cm16.ifp = self.ifp
        lv.cm16.ire = self.ire
        lv.cm16.jfp = self.jfp
        #lv.cm16.nph = self.icf / 2
        lv.cm16.ine = self.ine
        # CM18
        lv.cm18.ipar = self.ipar
        lv.cm18.iopr = self.iopr
        # CM34
        lv.cm34.mda = self.md
        lv.cm34.iwt = self.iwt
        lv.cm34.ixw = 1 * ((self.ixi == 1) or (self.iwt == 1))
        lv.cm34.ipl = self.ipl
        # CM35
        lv.cm35.ipf = self.ipf
        lv.cm35.nprint = self.nprint
        # CM39
        lv.cm39.p8 = self.parameters[7]
        lv.cm39.p18 = self.parameters[17]
        # CM78
        lv.cm78.dattyc = self.dattyp
        # CM79
        lv.cm79.ytt = self.y

        # Run MAINCLC if MAXFEV > 3
        if self.maxfev > 3:
            self.outputvals = np.zeros(self.ky)
            lv.mainclc(self.ky, self.ftol, 0.0, self.xtol, self.x, self.maxfev, self.nprint, 2, self.parameters, self.nfrei, self.outputvals)

        # If MAXFEV = 0 no fit calc new data
        # If MAXFEV = 1 no fit convert
        # If MAXFEV = 2 no fit calc new data without parameters with NFREE = 3

        # Write outputs

        return
