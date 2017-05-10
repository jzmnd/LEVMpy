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
import levmpyfortran.LEVMpyFortran as _lv
import numpy as np
from utils import *

__author__ = "Jeremy Smith"
__version__ = "1.1"


EPV = 8.85418782e-14
NPAFR = len(_lv.cm2.ns)
NTOT = len(_lv.cm2.nfree)
NPTS = len(_lv.cm1.freq)
NPT2 = len(_lv.cm2.y)


class Experiment():
    """Class for single LEVM input file"""
    def __init__(self, infile, outsuffix="OUTFILE", path=os.getcwd()):
        self.path = path
        self.infile = infile
        self.outfile = infile + '_' + outsuffix
        self.fitted = False
        # Read input file
        self._readinfile()
        # Set flags for MAINCLC
        self._setflags()
        # Set defaults
        self._defaults()
        # Set parameter flags, arrays and indices
        self._setparams()
        # Combine y into single array pre-transform
        self.ytt = np.hstack((self.y1, self.y2))
        # Perform data transformations if required
        self._transform()
        # Set circuit flags
        self._setcircuitsp()
        # Combine y and r into single array
        self.y = np.hstack((self.y1, self.y2))
        self.r = np.hstack((self.r1, self.r2))
        # Define empty arrays for data output
        self.outputvals = np.zeros(NPT2)

    def _readinfile(self):
        """Reads LEVM input file"""
        with open(os.path.join(self.path, self.infile), 'r') as f:
            # Line 1 - description of run
            self.alpha = f.readline().strip()

            # Line 2 - data and fit types
            line2 = f.readline()
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
            line3 = f.readline()
            self.md = int(line3[:5])
            self.n = int(line3[6:10])
            self.maxfev = int(line3[11:15])
            self.nprint = int(line3[16:20])
            self.irch = int(line3[21:25])
            self.mode = int(line3[26:30])
            self.icp = int(line3[31:35])
            self.iprint = int(line3[36:40])
            self.igacc = int(line3[41:45])
            self.atemp = dfloat(line3[45:55])

            # Line 4 - set model parameters
            plist = []
            for i in range(int(np.ceil(self.n / 5.0))):
                line = f.readline().strip().split()
                for p in line:
                    plist.append(dfloat(p))

            self.parameters = np.array(plist, dtype=float)
            self.parameters.resize(NTOT)

            # Line 5 - set nfree
            nflist = list(f.readline().strip())
            self.nfree = np.array(nflist, dtype=int)

            # Line 6 - read input data
            flist = []
            y1list = []
            y2list = []
            for i in range(abs(self.md)):
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
                for i in range(abs(self.md)):
                    line = f.readline().strip().split()
                    r1list.append(dfloat(line[1]))
                    r2list.append(dfloat(line[2]))
                self.r1 = np.array(r1list, dtype=float)
                self.r2 = np.array(r2list, dtype=float)
            else:
                self.r1 = np.zeros(abs(self.md)) 
                self.r2 = np.zeros(abs(self.md))
        return

    def _setflags(self):
        """Set flags for MAINCLC"""
        self.iacy = 0 
        if self.iopt < 0:
            self.iacy = abs(self.iopt)
            self.iopt = 0
        if self.iacy == 0:
            self.ftol = 1e-15
            self.xtol = 1e-28
        else:
            self.ftol = 10**(self.iacy)
            self.xtol = self.ftol

        self.ipf = ['R', 'P', 'D', 'L'].index(self.pfit)

        self.iorig = 1 if self.ipar < 0 else 0
        self.jfp = 1 if self.ifp > 0 else 0
        self.ifp = abs(self.ifp)
        self.iopr = 1 if self.ire < -10 else 0
        self.imd = 1 if self.md > 0 else -1
        self.md = abs(self.md)

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
            self.qq = 'UNIT' if self.irch == 1 else 'DATA'

        self.ipl = self.iprint if self.iprint < 0 else 0
        self.iprint = abs(self.iprint)

        self.irg = 0
        if self.parameters[30] < 0:
            self.irg = int(round(abs(self.parameters[30])))
            self.parameters[30] = 0
        return

    def _defaults(self):
        """Set default values"""
        if self.celcap == 0:
            self.celcap = 1.0
        if self.imd == -1:
            self.celcap = EPV
        self.clcap = 1.0 if self.celcap < 0 else self.celcap
        return

    def _setparams(self):
        """Set parameter flags, arrays and indices"""
        # Set NS and X parameter arrays
        self.ns = np.where(self.nfree > 0)[0]   # 0-start index
        self.ns1 = self.ns + 1                  # 1-start index
        self.x = self.parameters[self.ns]
        self.nfrei = len(self.x)
        self.nfixi = self.n - self.nfrei

        # Set free parameter flags
        self.ixi = 1 if (self.nfree[31] > 0) and (self.irch > 1) else 0
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

        # Set pex values as "exact" inputs
        if self.iorig == 1:
            self.pex = resized(self.parameters, NTOT)
        else:
            self.pex = np.zeros(NTOT)
        return

    def _transform(self):
        """Perform data transformations"""
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
        return

    def _setcircuitsp(self):
        """Set circuit specific flags"""
        # K,R,S-circuit flags
        self.icf = 0
        self.wf = 0.0
        if (self.fun == 'K') or (self.fun == 'R') or (self.fun == 'S'):
            if (self.fun == 'K'):
                jia = 22
                jib = 18
            else:
                jia = 28
                jib = 12
            ip34a = int(round(abs(self.parameters[33])))
            ip34 = np.sign(self.parameters[33]) * ip34a
            ip35a = int(round(abs(self.parameters[34])))
            ip38a = 1 if ip34 > 0 else int(round(abs(self.parameters[37])))
            ip40a = int(round(abs(self.parameters[39])))
            ip40s = np.sign(self.parameters[39])
            if self.ire > -10:
                self.ire = -10
            if ip35a > 0:
                self.icf = 2 * ip35a
            else:
                self.icf = np.count_nonzero(self.parameters[:jia])
                if self.n > 40:
                    self.icf += np.count_nonzero(self.parameters[40:])

            if ip38a == 1:
                if int(self.parameters[35]) == 0:
                    if self.freq[self.md - 1] > self.freq[0]:
                        taumin = 1.0 / self.freq[self.md - 1]
                        taumax = 1.0 / self.freq[0]
                    else:
                        taumin = 1.0 / self.freq[0]
                        taumax = 1.0 / self.freq[self.md - 1]
                    frat = taumax / taumin
                else:
                    taumin = self.parameters[35]
                    frat = self.parameters[38] / taumin
                self.wf = (1.0/(self.icf / 2 - 1)) * np.log(frat)

            if not ((ip38a != 1) and (ip40a <= 1) and (ip34 <= 0)):
                if 0 < ip34 < 8:
                    self.maxfev = 1
                    self.irch = 2
                    if (ip34 < 3) or (ip34 == 5):
                        self.parameters[37] = 1
                        if ip40a != 4: self.parameters[39] = 3 * ip40s
                    elif (ip34 == 3) or (ip34 == 4):
                        self.parameters[37] = 2
                        if ip40a != 4: self.parameters[39] = 2 * ip40s
                    elif ip34 == 6:
                        self.parameters[37] = 2
                    else:
                        self.parameters[37] = 2
                        self.parameters[39] = 2 * ip40s
                #if self.fun != 'S':
                    # RESORT(P,NFREE,TAUM,WF,X,NTOT,NS,JIA,JIB)

        # O-circuit NELEM and flags
        self.nelem = 0
        if self.fun == 'O':
            if (int(self.parameters[34]) != 4) and (int(self.parameters[39]) == 4):
                self.nelem = int(self.parameters[19])
            else:
                self.nelem = int(self.parameters[9])
        if (self.fun == 'O') and (self.nelem in [7, 10, 32, 36]):
            self.inde = 0
        elif self.ire < 0:
            self.inde = 1
        else:
            self.inde = 2
        if (self.nelem == 7) and (int(self.parameters[6]) != 1):
            self.parameters[9] = 6
            self.nelem = 6
            self.inde = 1
        return

    def fit(self):
        """Performs the CNLS fit"""
        # Reset outputs
        self.outputvals = np.zeros(NPT2)
        self.nfev = 2
        self.fitted = False
        # Set the levmpy common variables
        self._setcommon()
        # Print fit information
        self._print_fit_info(sys.stdout.write)

        # Run MAINCLC if MAXFEV > 3
        if self.maxfev > 3:
            self.result = _lv.mainclc(self.ky, self.ftol, 0.0, self.xtol, self.x,
                                      self.maxfev, self.nprint, self.nfev, self.pex, self.nfrei, self.outputvals)
            self.fitted = True

        # If MAXFEV = 0 no fit calc new data
        # If MAXFEV = 1 no fit convert
        # If MAXFEV = 2 no fit calc new data without parameters with NFREE = 3

        # Resdiduals and errors
        self.r = _lv.cm2.r[self.jy - 1:self.ky]
        self.res = self.y - self.outputvals[self.jy - 1:self.ky]
        self.resmod = self.res / self.outputvals[self.jy - 1:self.ky]
        # Relative standard deviation of parameters
        self.rxsd = _lv.cm20.rxsd[:self.nfrei]
        # Fit information
        self.info = int(_lv.cm34.infp)
        self.fnorm = float(_lv.cm4.fnorm)
        self.fqq = float(_lv.cm11.fqq)
        self.ndf = int(_lv.cm11.ndf)
        # Write outputs
        self._print_outputs(sys.stdout.write)
        return

    def _setcommon(self):
        """Set all COMMON block variables"""
        # CM1
        _lv.cm1.freq = resized(self.omega, NPTS)
        _lv.cm1.m = self.md
        _lv.cm1.dattyp = self.dattyp
        # CM2
        _lv.cm2.y = resized(self.y, NPT2)
        _lv.cm2.r = resized(self.r, NPT2)
        _lv.cm2.fj = np.zeros(NPT2)
        _lv.cm2.p = self.parameters
        _lv.cm2.drss = 1.0e-8
        _lv.cm2.roe = self.roe
        _lv.cm2.ns = resized(self.ns1, NPAFR)
        _lv.cm2.nfree = resized(self.nfree, NTOT)
        _lv.cm2.n = self.n
        _lv.cm2.icnt = 0
        _lv.cm2.mn = 0
        _lv.cm2.irch = self.irch
        _lv.cm2.ixi = self.ixi
        _lv.cm2.dattyq = self.dattyp
        # CM3
        _lv.cm3.celcap = self.celcap
        _lv.cm3.fun = self.fun
        _lv.cm3.dfit = self.dfit
        _lv.cm3.pfit = self.pfit
        # CM4
        _lv.cm4.fnorm = 0.0
        # CM5
        _lv.cm5.tomega = 0.0
        _lv.cm5.phicom = 0.0
        _lv.cm5.iv = 0
        # CM9
        _lv.cm9.tomegax = 0.0
        _lv.cm9.phix = 0.0
        _lv.cm9.ichg = 0
        _lv.cm9.iwtx = 0
        # CM10
        _lv.cm10.epsg = 10**(-abs(self.igacc))
        _lv.cm10.izr = self.izr
        # CM11
        _lv.cm11.mqy = self.nfrei - self.ins
        _lv.cm11.ispr = 0
        _lv.cm11.icx = 0
        _lv.cm11.ndf = 0
        _lv.cm11.fqq = 0.0
        # CM12
        _lv.cm12.clcap = self.clcap
        _lv.cm12.atemp = self.atemp
        _lv.cm12.wf = self.wf
        _lv.cm12.maxfev = self.maxfev
        _lv.cm12.icf = self.icf
        _lv.cm12.mde = self.mode
        _lv.cm12.jcdx = 0
        # CM13
        _lv.cm13.rx = 0.0
        _lv.cm13.tx = 0.0
        _lv.cm13.ux = 0.0
        _lv.cm13.phi = 0.0
        _lv.cm13.xxm1 = 0.0
        _lv.cm13.xx1 = 0.0
        _lv.cm13.xx2 = 0.0
        _lv.cm13.xx3 = 0.0
        _lv.cm13.rn = 0.0
        _lv.cm13.ain = 0.0
        _lv.cm13.icav = 0
        _lv.cm13.nelem = self.nelem
        _lv.cm13.nch = 0
        # CM14
        _lv.cm14.fjacc = np.zeros((NPT2, NPAFR))
        # CM16
        _lv.cm16.iop = self.iopt
        _lv.cm16.iorig = self.iorig
        _lv.cm16.nyc = 0
        _lv.cm16.jy = self.jy
        _lv.cm16.iprint = self.iprint
        _lv.cm16.ldfjac = self.ky
        _lv.cm16.mode = self.mode
        _lv.cm16.ifp = self.ifp
        _lv.cm16.ire = self.ire
        _lv.cm16.istp = 0
        _lv.cm16.jfp = self.jfp
        _lv.cm16.nph = self.icf / 2
        _lv.cm16.ine = self.ine
        # CM18
        _lv.cm18.sdwc = 0.0
        _lv.cm18.sdrc = 0.0
        _lv.cm18.diag = np.zeros(NPAFR)
        _lv.cm18.ipar = self.ipar
        _lv.cm18.iopr = self.iopr
        # CM20
        _lv.cm20.rxsd = np.zeros(NPAFR)
        # CM29
        _lv.cm29.gam = 0.0
        _lv.cm29.sn1 = 0.0
        _lv.cm29.cs1 = 0.0
        _lv.cm29.pii = 0.0
        _lv.cm29.cdn = 0.0
        _lv.cm29.qpi = 0.0
        _lv.cm29.mdex = 0
        # CM34
        _lv.cm34.mda = self.md
        _lv.cm34.iwt = self.iwt
        _lv.cm34.ixw = 1 * ((self.ixi == 1) or (self.iwt == 1))
        _lv.cm34.infp = 0
        _lv.cm34.ipl = self.ipl
        # CM35
        _lv.cm35.jit = 0
        _lv.cm35.ipf = self.ipf
        _lv.cm35.nprint = self.nprint
        # CM36
        _lv.cm36.shw = np.zeros(NPT2)
        _lv.cm36.isw = 0
        # CM39
        _lv.cm39.p8 = self.parameters[7]
        _lv.cm39.p18 = self.parameters[17]
        # CM47
        _lv.cm47.icount = 0
        # CM55
        _lv.cm55.px1 = 0.0
        _lv.cm55.px41 = 0.0
        _lv.cm55.px45 = 0.0
        # CM73
        _lv.cm73.p39 = self.parameters[38]
        # CM78
        _lv.cm78.dattyc = self.dattyp
        # CM79
        _lv.cm79.ytt = resized(self.ytt, NPT2)
        # CM80
        _lv.cm80.omegat = 0.0
        _lv.cm80.psi = 0.0
        _lv.cm80.gama = 0.0
        _lv.cm80.deli = 0.0
        return

    def _print_fit_info(self, fp):
        """Prints fit information"""
        fp("\nLEVM : COMPLEX NONLINEAR LEAST SQUARES IMMITTANCE DATA FITTING PROGRAM\n")
        fp("       VERSION 8.06 - 2/05\n\n")
        fp("{:s}\n".format(self.alpha))
        fp("  DATA ENTERED IN {:s}{:s} FORMAT TO BE USED IN {:s}{:s} FIT\n".format(self.dinp, self.pinp, self.dfit, self.pfit))
        fp("  CIRCUIT MODEL : {:s}\n\n".format(self.fun))
        fp("  *****  FIT OF {:s} DATA  *****\n\n".format(self.dattyp))

        fp("  # OF DATA POINTS = {:d}   WEIGHT: IRCH = {:d}   # OF FREE PARAMETERS = {:d}\n".format(self.md, self.irch, self.nfrei))
        fp("  PRINTS EVERY {:d} ITERATIONS   MAX # ITERATIONS = {:d}   MAIN WT USES: {:s}\n".format(self.nprint, self.maxfev, self.qq))
        fp("  CELL CAPACITANCE = {:e}\n".format(self.celcap))

        if self.ixi == 0:
            if self.irch == 0:
                fp("  WEIGHTS (STANDARD DEVIATIONS) ARE READ IN\n")
            if self.irch == 1:
                fp("  UNIT WEIGHTING ASSIGNED TO EACH POINT\n")
            if self.irch == 2:
                fp("  WEIGHTS INVOLVE THE MAGNITUDES OF DATA OR FUNCTION VALUES RAISED TO THE POWER {:e}\n".format(self.parameters[31]))
            if self.irch == 3:
                fp("  MODULUS WEIGHTING: RESULTS RAISED TO THE POWER {:e}\n".format(self.parameters[31]))
            if self.irch == 4:
                fp("  #1174 SPINOLO FRA WEIGHTING\n")
            if self.irch == 5:
                fp("  #1250 & 1286 ORAZEM-AGARWAL FRA WEIGHTING\n")
            if self.irch == 6:
                fp("  #1250 & PAR 273 ORAZEM-AGARWAL FRA WEIGHTING\n")
        else:
            fp("  WEIGHTS USE XI OR U**2 AND XI: BOTH MAY BE FREE PARAMETERS\n")

        fp("\n  INITIAL PARAMETER GUESSES AND FIXED (0) OR FREE (1 OR 2) STATUS\n")
        for i in range(16):
            j = i + 16
            fp("     P({:2d}) = {: e}   {:d}     P({:2d}) = {: e}   {:d}\n".format(i + 1, self.parameters[i], self.nfree[i],
                                                                                       j + 1, self.parameters[j], self.nfree[j]))
        if self.n > 32:
            fp("\n  THE FOLLOWING PARAMETERS ARE ALWAYS FIXED\n")
            for i in range(32, 36):
                j = i + 4
                fp("     P({:2d}) = {: e}         P({:2d}) = {: e}\n".format(i + 1, self.parameters[i], j + 1, self.parameters[j]))
        fp("\n")
        return

    def _print_outputs(self, fp):
        """Prints output data"""
        if self.info < 0:
            fp("  U & XI OR XI CONVERGENCE PROBLEM\n\n")
        if self.info == 0:
            fp("  IMPROPER INPUT PARAMETERS\n\n")
        if self.info == 1:
            fp("  ALGORITHM TERMINATES HAVING ESTIMATED THAT THE RELATIVE ERROR\n")
            fp("  IN THE SUM OF SQUARES IS AT MOST FTOL\n\n")
        if self.info == 2:
            fp("  ALGORITHM TERMINATES HAVING ESTIMATED THAT THE RELATIVE ERROR\n")
            fp("  BETWEEN THE VECTOR OF FREE PARAMETERS AND THE DESIRED SOLUTION IS AT MOST XTOL\n\n")
        if self.info == 3:
            fp("  ALGORITHM TERMINATES HAVING ESTIMATED THAT THE RELATIVE ERROR\n")
            fp("  IN THE SUM OF SQUARES IS LESS THAN FTOL AND THE L2 NORM OF\n")
            fp("  THE VECTOR OF FREE PARAMETERS AND THE DESIRED SOLUTION IS LESS THAN XTOL\n\n")
        if self.info == 4:
            fp("  ALGORITHM TERMINATES HAVING ESTIMATED THAT FVEC IS ORTHOGONAL\n")
            fp("  TO THE COLUMNS OF THE JACOBIAN TO PRECISION GTOL\n\n")
        if self.info == 5:
            fp("  TERMINATION: NUMBER OF CALLS TO FCN WITH IFLAG=1 EXCEEDS MAXFEV\n\n")
        if self.info == 6:
            fp("  TERMINATION: FTOL IS TOO SMALL\n")
            fp("  NO FURTHER REDUCTION IN THE SUM OF SQUARES IS POSSIBLE\n\n")
        if self.info == 7:
            fp("  TERMINATION: XTOL IS TOO SMALL\n")
            fp("  NO FURTHER IMPROVEMENT IN THE FREE PARAMETERS IS POSSIBLE\n\n")

        fp("  FINAL L2 NORM OF THE RESIDUALS = {: e}\n".format(self.fnorm))
        fp("  NUMBER OF FUNCTION EVALUATIONS = {: d}\n".format(self.result[1]))
        fp("  FIT QUALITY FACTOR             = {: e}\n".format(self.fqq))
        fp("  NO. DEGREES OF FREEDOM         = {: d}\n".format(self.ndf))

        fp("\n  PARAMETER ESTIMATES         STD DEV         REL STD DEV\n")
        for i, n in enumerate(self.ns1):
            fp("     P({:2d}) = {: e}    {:e}    {:e}\n".format(n, self.x[i], self.rxsd[i] * self.x[i], self.rxsd[i]))
        fp("\n  PDAV = {:e}\n".format(np.sum(self.rxsd) / (self.nfrei - self.ins)))
        fp("  PDRMS = {:e}\n".format(np.sqrt(np.sum(self.rxsd**2) / (self.nfrei - self.ins))))

        fp("\n  OBSERVED VARIABLES\n")
        fp("  M  T   MEASURED      ESTIMATED     UNCERTANTY    RESIDUALS     RESID/MODEL\n")
        for i in range(self.md):
            j = i + self.md
            fp("{:3d}  R  {: e} {: e} {: e} {: e} {: e}\n".format(i + 1, self.y[i], self.outputvals[i],
                                                                      self.r[i], self.res[i], self.resmod[i]))
            fp("{:3d}  I  {: e} {: e} {: e} {: e} {: e}\n".format(i + 1, self.y[j], self.outputvals[j],
                                                                      self.r[j], self.res[j], self.resmod[j]))
        return

    def writetofile(self):
        """Writes current data to output file"""
        if self.fitted:
            with open(os.path.join(self.path, self.outfile), 'w') as f:
                self._print_fit_info(f.write)
                self._print_outputs(f.write)
        return
