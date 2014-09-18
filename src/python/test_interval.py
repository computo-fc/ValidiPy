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
inf = calc_inf()

def test_instances():
    assert isinstance(a, Interval)
    assert isinstance(Interval(mpfr("0.1")), Interval)

def test_equalities():
    assert Interval(0.25)-1.0/4 == Interval(0)
    assert not(a == b)
    assert a != b
    assert a == Interval(a.lo,a.hi)
    assert Interval(0.1) == Interval(0.0999999999999999999,0.1)
    assert Interval(0).reciprocal() == Interval(inf,inf)
    assert Interval(0,1).reciprocal() == Interval(1,inf)
    assert Interval(1,inf).reciprocal() == Interval(0,1)
    assert c.reciprocal() == c
    assert 1/b == b.reciprocal()
    assert a.intersection( b.hull(a)) == a

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

def test_mig_mag():
    assert 0.1 in Interval(0.1)
    assert (Interval(-2,2)).mig() == 0.0
    assert (-b).mag() == b.hi
    assert a.diam() == a.hi - a.lo

def test_pow():
    assert Interval(-3,2)**2 == Interval(mpfr('0.0'), mpfr('9.0'))
    assert Interval(-3,2)**3 == Interval(mpfr('-27.0'), mpfr('8.0'))
    assert Interval(-3,4)**0.5 == Interval(mpfr('0.0'), mpfr('2.0'))
    assert Interval(1,27)**(1.0/3) == Interval(1,3)
    assert Interval(-3,2)**Interval(2) == Interval(mpfr('0.0'), mpfr('9.0'))
    assert Interval(-3,4)**Interval(0.5) == Interval(-3,4)**0.5 == Interval(0, 2)
    assert Interval(0.1,0.7)**(1.0/3) == Interval(mpfr('0.46415888336127786'), mpfr('0.88790400174260076'))

def test_exp_log():
    assert exp(Interval(0.5)) == Interval(mpfr('1.648721270700128'), mpfr('1.6487212707001282'))
    assert exp(Interval(0.1)) == Interval(mpfr('1.1051709180756475'), mpfr('1.1051709180756477'))
    assert log(Interval(0.5)) == Interval(mpfr('-6.931471805599454e-01'), mpfr('-6.9314718055994529e-01'))
    assert log(Interval(0.1)) == Interval(mpfr('-2.3025850929940459e+00'), mpfr('-2.3025850929940455e+00'))

def test_sin_cos_tan():
    assert sin(Interval(0.5)) == Interval(mpfr('0.47942553860420295'), mpfr('0.47942553860420301'))
    assert sin(Interval(0.5,1.67)) == Interval(mpfr('0.47942553860420295'), mpfr('1.0'))
    assert sin(Interval(2.1, 5.6)) == Interval(mpfr('-1.0'), mpfr('0.86320936664887404'))
    assert sin(Interval(0.5,8.5)) == Interval(-1.0, 1.0)
    # The following does not correspond to the result of julia
    #?assert sin(Interval(1.67,3.2)) == Interval(mpfr('-0.058374143427580093'), mpfr('0.99508334981018021'))
    assert sin(Interval(1.67,3.2)) == Interval(mpfr('-0.058374143427580086'), mpfr('0.9950833498101801'))

    assert cos(Interval(0.5)) == Interval(mpfr('0.87758256189037265'), mpfr('0.87758256189037276'))
    assert cos(Interval(0.5,1.67)) == Interval(mpfr('-0.099041036598728246'), mpfr('0.87758256189037276'))
    assert cos(Interval(2.1, 5.6)) == Interval(mpfr('-1.0'), mpfr('0.77556587851025016'))
    assert cos(Interval(0.5,8.5)) == Interval(mpfr('-1.0'), mpfr('1.0'))
    assert cos(Interval(1.67,3.2)) == Interval(mpfr('-1.0'), mpfr('-0.099041036598728011'))

    assert tan(Interval(0.5)) == Interval(mpfr('0.54630248984379048'), mpfr('0.5463024898437906'))
    assert tan(Interval(0.5,1.67)) == Interval(mpfr('-inf'), mpfr('inf'))
    assert tan(Interval(1.67,3.2)) == Interval(mpfr('-10.047182299210307'), mpfr('0.058473854459578652'))

