# -*- coding: utf-8 -*-

# run using the command "nosetests" in the current directory
# nosetests finds tests with name "test_..."

from interval import *
from interval_tools import *
import numpy as np

# Allow errors using `raise`:
from nose.tools import *

a = Interval(1.1,0.1)
b = Interval(0.9,2.0)
c = Interval(0.25,4.0)

def test_instances():
    assert isinstance(a, Interval)
    assert isinstance(Interval(mpfr("0.1")), Interval)

def test_equalities():
    assert Interval(0.25)-1.0/4 == Interval(0)
    assert not(a == b)
    assert a != b
    assert a == Interval(a.lo,a.hi)
    assert Interval(0.1) == Interval(0.0999999999999999999,0.1)

def test_add_sub_mul_div():
    assert a+b == Interval(mpfr('0.99999999999999989'),mpfr('3.1'))
    assert -a == Interval(-1.1,-0.0999999999999999999)
    assert a-b == Interval(mpfr("-1.9000000000000001e+00"), mpfr("2.0000000000000018e-01"))
    assert a*b == Interval(a.lo*b.lo, a.hi*b.hi)
    assert 10*a == Interval(mpfr("9.9999999999999989e-01"), mpfr("1.1000000000000002e+01"))
    assert b/a == Interval(mpfr("8.181818181818179e-01"), mpfr("2.0000000000000004e+01"))
    assert a/c == Interval(mpfr("0.024999999999999998"),mpfr("4.4"))
    assert c/4.0 == Interval(0.0625,1.0)
    assert Interval.reciprocal(Interval(0)) == Interval(mpfr('inf'), mpfr('inf'))
    assert Interval.reciprocal(Interval(0,1)) == Interval(1,mpfr('inf'))
    assert Interval.reciprocal(Interval(1,mpfr('inf'))) == Interval(0,1)
    assert Interval.reciprocal(c) == c
    assert 1/b == Interval.reciprocal(b)
    assert 0.1 in Interval(0.1)
    assert not(Interval.strictly_contains(Interval(0.1), 0.1))

def test_pow():
    assert Interval(-3,2)**2 == Interval(mpfr('0.0'), mpfr('9.0'))
    assert Interval(-3,2)**3 == Interval(mpfr('-27.0'), mpfr('8.0'))
    assert Interval(-3,4)**0.5 == Interval(mpfr('0.0'), mpfr('2.0'))
    assert Interval(1,27)**(1.0/3) == Interval(1,3)
    assert Interval(-3,2)**Interval(2) == Interval(mpfr('0.0'), mpfr('9.0'))
    assert Interval(-3,4)**Interval(0.5) == Interval(mpfr('0.0'), mpfr('2.0'))
    assert Interval(0.1,0.7)**(1.0/3) == Interval(mpfr('0.46415888336127786'), mpfr('0.88790400174260076'))


# @raises(ValueError)
# def test_log():

#     for i in range(51):
#         a = random_interval()
#         b = a.log()
#         #
#         if 0 in a:
#             assert b.lo == mpf('-inf') and b.hi == mp.log(a.hi)
#         #
#         elif 0 < a.lo:
#             assert b.lo == mp.log(a.lo) and b.hi == mp.log(a.hi)
#         #
#         elif 0 > a.hi:
#             raise ValueError("Lanzar este error es ok")
#     #
