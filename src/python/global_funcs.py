# -*- coding:utf-8 -*-

import gmpy2
from gmpy2 import mpz, mpq, mpfr#, mpc
from gmpy2 import RoundDown, RoundUp

from matplotlib import pyplot as plt
import numpy as np

class MemGlobal(object):
    def __init__(self, num_library=gmpy2):
        self.libmath = num_library
        #if num_library == np:
        #    mem.libmath_random = np.random.uniform

    def change_to_numpy(self):
        self.libmath = np
        return self

    def change_to_gmpy2(self):
        self.libmath = gmpy2
        return self



mem = MemGlobal()

mpfr_type = type(gmpy2.mpfr(0))


# Some constants
def calc_pi():
    if mem.libmath == np:
        pi = np.pi
    elif mem.libmath == gmpy2:
        pi = gmpy2.const_pi()
    else:
        raise TypeError("Constants not implemented. Error in 'mem.libmath'")

    return pi

def calc_inf():
    if mem.libmath==np:
        inf = np.inf
    elif mem.libmath==gmpy2:
        inf = gmpy2.inf(1)
    else:
        raise TypeError("Constants not implemented. Error in 'mem.libmath'")

    return inf

# Access to functions
def exp(a):
    try:
        return a.exp()
    except:
        return mem.libmath.exp(a)

def log(a):
    try:
        return a.log()
    except:
        return mem.libmath.log(a)

def sqrt(a):
    try:
        return a.sqrt()
    except:
        return mem.libmath.sqrt(a)

def sin(a):
    try:
        return a.sin()
    except:
        return mem.libmath.sin(a)

def cos(a):
    try:
        return a.cos()
    except:
        return mem.libmath.cos(a)

def tan(a):
    try:
        return a.tan()
    except:
        return mem.libmath.tan(a)
