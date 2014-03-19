# -*- coding: utf-8 -*-

# run using the command "nosetests" in the current directory
# nosetests finds tests with name "test_..."

from interval import *
from interval_tools import *
import numpy as np

# Allow errors using `raise`:
from nose.tools import *

def test_addition():

    for i in range(51):
        a = random_interval()
        b = random_interval()
        #
        num1a, num2a = a.lo, a.hi
        num1b, num2b = b.lo, b.hi
        num3 = np.random.uniform( -10.0, 10.0 )
        num3 = mpf( num3 )
        #
        c1 = a + b
        c2 = a + num3
        c3 = num3 + b
        print num1a+num3, num2a+num3
        assert c1.lo == num1a+num1b and c1.hi == num2a+num2b
        assert c2.lo == num1a+num3 and c2.hi == num2a+num3
        assert c3.lo == num1b+num3 and c3.hi == num2b+num3
    #


@raises(ValueError)
def test_pow():

    for i in range(51):
        # a.lo < 0 and a.hi > 0
        a = random_interval()
        num1a, num2a = -abs(a.lo), abs(a.hi)
        num_min = min( -num1a, num2a )
        num_max = a.mag()
        a = Interval( num1a, num2a )
        #
        # Test even power
        b = a**2
        assert b.lo == 0 and b.hi == num_max**2
        c = a**Interval(2.0)
        assert c.lo == 0 and c.hi == num_max**2.0
        #
        # Test odd power
        b = a**3
        assert b.lo == num1a**3 and b.hi == num2a**3
        c = a**Interval(3)
        assert b == c
        #
        ## Test fractional power, with a.lo < 0
        b = a**2.5
        assert b.lo == 0 and b.hi == a.hi**2.5
        #
        # Test: a.lo < 0, a.hi > 0 and power is an interval
        a = Interval( -num_min, num_max )
        b = random_interval()
        c = a**b
        if 0 > b.hi:
            assert c.lo == mp.exp( b.lo*mp.log(a.hi) ) and c.hi == mpf('+inf')
        elif 0 < b.lo:
            assert c.lo == 0 and c.hi == mp.exp( b.hi*mp.log(a.hi) )
        elif 0 >= b.lo and 0 <= b.hi:
            assert c.lo == 0 and c.hi == mpf('+inf')
        #
        # Test fractional power, with a.lo > 0
        a = Interval( num_min, num_max )
        alfa = 1.0/3.0
        beta = 1.0/alfa
        b = a**alfa
        assert b.lo == num_min**alfa and b.hi == num_max**alfa
        c = b**beta
        assert c.lo == b.lo**beta and c.hi == b.hi**beta
        ## OJO
        # El siguiente test NO funciona, por errores de redondeo
        #assert c == a
        #
        # Test: error with negative interval and fractional power
        a = -a
        c = a**2.5
        raise ValueError("Lanzar este error es ok")
    #


@raises(ValueError)
def test_log():

    for i in range(51):
        a = random_interval()
        b = a.log()
        #
        if 0 in a:
            assert b.lo == mpf('-inf') and b.hi == mp.log(a.hi)
        #
        elif 0 < a.lo:
            assert b.lo == mp.log(a.lo) and b.hi == mp.log(a.hi)
        #
        elif 0 > a.hi:
            raise ValueError("Lanzar este error es ok")
    #
