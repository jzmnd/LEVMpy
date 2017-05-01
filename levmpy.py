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


def resized(arr, s):
    """Returns resized array padded with zeros"""
    tmparr = np.copy(arr)
    tmparr.resize(s)
    return tmparr


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
        self.fitted = False

        # Read input file
        self._readinfile()

        # Set flags
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
        self.ipf = ['R', 'P', 'D', 'L'].index(self.pfit)
        self.iorig = 0
        if self.ipar < 0:
            self.iorig = 1
        self.jfp = 1
        if self.ifp < 0:
            self.ifp = abs(self.ifp)
            self.jfp = 0
        self.iopr = 0
        if self.ire < -10:
            self.iopr = 1
        self.imd = 1
        if self.md < 0:
            self.md = abs(self.md)
            self.imd = -1
        self.inp = 0
        self.izr = 0
        if self.n < 0:
            self.n = abs(self.n)
            self.inp = 1
            if (self.dattyp == 'C') and (self.dinp == 'Z') and (self.pinp == 'R') and (self.dfit == 'Z') and (self.pfit == 'R'):
                self.izr = 1
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
        self.ipl = 0
        if self.iprint < 0:
            self.iprint = abs(self.iprint)
            self.ipl = self.iprint

        # Set default CELCAP
        if self.celcap == 0:
            self.celcap = 1.0
        if self.imd == -1:
            self.celcap = EPV

        # Set parameter flags
        self.irg = 0
        if self.parameters[30] < 0:
            self.irg = int(round(abs(self.parameters[30])))
            self.parameters[30] = 0

        # Set NS and X parameter arrays
        self.ns = np.where(self.nfree > 0)[0]   # 0-start index
        self.ns1 = self.ns + 1                  # 1-start index
        self.x = self.parameters[self.ns]
        self.nfrei = len(self.x)
        self.nfixi = self.n - self.nfrei

        # Set free parameter flags
        self.ixi = 0
        if (self.nfree[31] > 0) and (self.irch > 1):
            self.ixi = 1

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

        # Set start and end indices for data
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

        # Set NELEM for O-circuit
        self.nelem = 0
        if self.fun == 'O':
            if (self.parameters[34] != 4) and (self.parameters[40] == 4):
                self.nelem = self.parameters[19]
            else:
                self.nelem = self.parameters[9]
        self.inde = 2
        if (self.fun == 'O') and (self.nelem in [7, 10, 32, 36]):
            self.inde = 0
        elif self.ire < 0:
            self.inde = 1
        else:
            self.inde = 2
        if (self.nelem == 7) and (self.parameters[6] != 1):
            self.parameters[9] = 6
            self.nelem = 6
            self.inde = 1

        # Combine y and r into single array
        self.y = np.hstack((self.y1, self.y2))
        self.r = np.hstack((self.r1, self.r2))

        # Set pex values as "exact" inputs
        if self.iorig == 1:
            self.pex = resized(self.parameters, NTOT)
        else:
            self.pex = np.zeros(NTOT)

    def _readinfile(self):
        """Reads LEVM input file"""
        with open(os.path.join(self.path, self.infile), 'r') as f:
            # Line 1 - description of run
            self.alpha = f.readline().strip()
            
            line2 = f.readline()
            line3 = f.readline()

            # Line 2 - data and fit types
            self.iopt = int(line2[:4])
            self.dinp = line2[6]
            self.dfit = line2[7]
            self.pinp = line2[8]
            self.pfit = line2[9]
            self.freqtyp = line2[10]
            self.neg = line2[11]
            self.fun = line2[12]
            self.celcap = dfloat(line2[14:23])
            self.dattyp = line2[23]
            self.ipar = int(line2[25:34])
            self.roe = dfloat(line2[35:44])
            self.ifp = int(line2[45:50])
            self.ire = int(line2[51:56])

            # Line 3 - data, weighting, and fitting specifications
            self.md = int(line3[:5])
            self.n = int(line3[6:10])
            self.maxfev = int(line3[11:15])
            self.nprint = int(line3[16:20])
            self.irch = int(line3[21:25])
            self.mode = int(line3[26:30])
            self.icp = int(line3[31:35])
            self.iprint = int(line3[36:40])
            self.igacc = int(line3[41:45])
            self.atemp = dfloat(line3[46:55])

            # Line 4 - set model parameters
            plist = []
            for i in range(int(self.n / 5)):
                line = f.readline().strip().split()
                for p in line:
                    plist.append(dfloat(p))

            self.parameters = np.array(plist, dtype=float)

            # Line 5 - set nfree
            nflist = list(f.readline().strip())
            self.nfree = np.array(nflist, dtype=int)

            # Line 6 - read input data
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

            # Line 7 - read fixed weighting parameters if they exist (IRCH=0)
            if self.irch == 0:
                r1list = []
                r2list = []
                for i in range(self.md):
                    line = f.readline().strip().split()
                    r1list.append(dfloat(line[1]))
                    r2list.append(dfloat(line[2]))
                self.r1 = np.array(r1list, dtype=float)
                self.r2 = np.array(r2list, dtype=float)
            else:
                self.r1 = np.zeros(self.md) 
                self.r2 = np.zeros(self.md)
        return

    def fit(self):
        """Performs the CNLS fit"""

        # Define empty arrays for data output
        self.outputvals = np.zeros(NPT2)

        # Print fit information
        print "\nLEVM : COMPLEX NONLINEAR LEAST SQUARES IMMITTANCE DATA FITTING PROGRAM"
        print "       VERSION 8.06 - 2/05\n"
        print "{:s}\n".format(self.alpha)
        print "  DATA ENTERED IN {:s}{:s} FORMAT TO BE USED IN {:s}{:s} FIT".format(self.dinp, self.pinp, self.dfit, self.pfit,)
        print "  CIRCUIT MODEL : {:s}\n".format(self.fun)
        print "  *****  FIT OF {:s} DATA  *****\n".format(self.dattyp)

        print "  # OF DATA POINTS = {:d}   WEIGHT: IRCH = {:d}   # OF FREE PARAMETERS = {:d}".format(self.md, self.irch, self.nfrei)
        print "  PRINTS EVERY {:d} ITERATIONS   MAX # ITERATIONS = {:d}   MAIN WT USES: {:s}".format(self.nprint, self.maxfev, self.qq)
        print "  CELL CAPACITANCE = {:e}".format(self.celcap)

        if self.ixi == 0:
            if self.irch == 0:
                print "  WEIGHTS (STANDARD DEVIATIONS) ARE READ IN"
            if self.irch == 1:
                print "  UNIT WEIGHTING ASSIGNED TO EACH POINT"
            if self.irch == 2:
                print "  WEIGHTS INVOLVE THE MAGNITUDES OF DATA OR FUNCTION VALUES RAISED TO THE POWER {:e}".format(self.parameters[31])
            if self.irch == 3:
                print "  MODULUS WEIGHTING: RESULTS RAISED TO THE POWER {:e}".format(self.parameters[31])
            if self.irch == 4:
                print "  #1174 SPINOLO FRA WEIGHTING"
            if self.irch == 5:
                print "  #1250 & 1286 ORAZEM-AGARWAL FRA WEIGHTING"
            if self.irch == 6:
                print "  #1250 & PAR 273 ORAZEM-AGARWAL FRA WEIGHTING"
        else:
            print "  WEIGHTS USE XI OR U**2 AND XI: BOTH MAY BE FREE PARAMETERS"

        print "\n  INITIAL PARAMETER GUESSES AND FIXED (0) OR FREE (1 OR 2) STATUS"
        for i in range(16):
            j = i + 16
            print "     P({:2d}) = {:e}   {:d}     P({:2d}) = {:e}   {:d}".format(i + 1, self.parameters[i], self.nfree[i],
                                                                                  j + 1, self.parameters[j], self.nfree[j])
        if self.n > 32:
            print "\n  THE FOLLOWING PARAMETERS ARE ALWAYS FIXED"
            for i in range(32, 36):
                j = i + 4
                print "     P({:2d}) = {:e}   {:d}     P({:2d}) = {:e}   {:d}".format(i + 1, self.parameters[i], self.nfree[i],
                                                                                      j + 1, self.parameters[j], self.nfree[j])
        # Set the levmpy common variables
        self._setcommon()

        # Run MAINCLC if MAXFEV > 3
        if self.maxfev > 3:
            self.result = lv.mainclc(self.ky, self.ftol, 0.0, self.xtol, self.x, self.maxfev, self.nprint, 2, self.pex, self.nfrei, self.outputvals)
            self.fitted = True

        # If MAXFEV = 0 no fit calc new data
        # If MAXFEV = 1 no fit convert
        # If MAXFEV = 2 no fit calc new data without parameters with NFREE = 3

        # Write outputs
        print "\nFITTED PARAMETERS:"
        print self.result[0]
        print "NFEV:", self.result[1]
        print "OUTPUT VALUES:"
        print self.result[2]

        return

    def _setcommon(self):
        """Set all COMMON block variables"""
        # CM1
        lv.cm1.freq = resized(self.omega, NPTS)
        lv.cm1.m = self.md
        lv.cm1.dattyp = self.dattyp
        # CM2
        lv.cm2.y = resized(self.y, NPT2)
        lv.cm2.r = resized(self.r, NPT2)
        lv.cm2.p = resized(self.parameters, NTOT)
        lv.cm2.drss = 1.0e-8
        lv.cm2.roe = self.roe
        lv.cm2.ns = resized(self.ns1, NPAFR)
        lv.cm2.nfree = resized(self.nfree, NTOT)
        lv.cm2.n = self.n
        lv.cm2.icnt = 0
        lv.cm2.mn = 0
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
        #####lv.cm12.wf = 
        #####lv.cm12.icf = 
        lv.cm12.maxfev = self.maxfev
        lv.cm12.mde = self.mode
        lv.cm12.jcdx = 0
        # CM13
        lv.cm13.icav = 0
        lv.cm13.nelem = self.nelem
        # CM16
        lv.cm16.iop = self.iopt
        lv.cm16.iorig = self.iorig
        lv.cm16.nyc = 0
        lv.cm16.jy = self.jy
        lv.cm16.iprint = self.iprint
        lv.cm16.ldfjac = self.ky
        lv.cm16.mode = self.mode
        lv.cm16.ifp = self.ifp
        lv.cm16.ire = self.ire
        lv.cm16.jfp = self.jfp
        #####lv.cm16.nph = self.icf / 2
        lv.cm16.ine = self.ine
        # CM18
        lv.cm18.ipar = self.ipar
        lv.cm18.iopr = self.iopr
        # CM34
        lv.cm34.mda = self.md
        lv.cm34.iwt = self.iwt
        lv.cm34.ixw = 1 * ((self.ixi == 1) or (self.iwt == 1))
        lv.cm34.infp = 0
        lv.cm34.ipl = self.ipl
        # CM35
        lv.cm35.jit = 0
        lv.cm35.ipf = self.ipf
        lv.cm35.nprint = self.nprint
        # CM36
        #####lv.cm36.shw = 
        lv.cm36.isw = 0
        # CM39
        lv.cm39.p8 = self.parameters[7]
        lv.cm39.p18 = self.parameters[17]
        # CM47
        lv.cm47.icount = 0
        # CM78
        lv.cm78.dattyc = self.dattyp
        # CM79
        lv.cm79.ytt = resized(self.y, NPT2)
        return
