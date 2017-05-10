#! /usr/bin/env python
"""
utils.py
Utility functions for levmpy

Fortran source code by J R McDonald
Based on LEVM version 8.06

Created by Jeremy Smith on 2017-04-16
"""

import os
import sys
import numpy as np

__author__ = "Jeremy Smith"
__version__ = "1.1"


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
