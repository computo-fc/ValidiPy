# -*- coding:utf-8 -*-

"""The x in the argument of functions is actually an object which is itself a derivative pair"""


import numpy as math  # NB
import numpy as np
from global_funcs import *


class AutoDiff:
    """Multidim AutoDiff"""
    def __init__(self, valor, deriv=0.0):
        self.valor = valor
        self.deriv = deriv

    def __add__(self, otro):
        if not isinstance(otro, AutoDiff):
            otro = AutoDiff(otro)

        return AutoDiff(self.valor+otro.valor, self.deriv+otro.deriv)

    def __radd__(self, otro):
        return self+otro

    def __sub__(self, otro):
        if not isinstance(otro, AutoDiff):
            otro = AutoDiff(otro)

        return AutoDiff(self.valor-otro.valor, self.deriv-otro.deriv)

    def __neg__(self, otro):
        return AutoDiff(-self.valor, -self.otro)

    def __rsub__(self, otro):
        return (self - otro).neg()

    def __pow__(self, n):
        if not isinstance(n, int):
            raise ValueError("Powers not implemented for non-integers")

        return AutoDiff(self.valor**n, n*self.valor**(n-1) * self.deriv)

    def __repr__(self):
        return "AutoDiff({}, {})".format(self.valor, self.deriv)

    def __mul__(self, otro):
        if not isinstance(otro, AutoDiff):
            otro = AutoDiff(otro)

        return AutoDiff(self.valor*otro.valor, 
                    self.valor*otro.deriv + self.deriv*otro.valor)

    def __rmul__(self, otro):
        return self*otro

    def __div__(self, otro):
        if not isinstance(otro, AutoDiff):
            otro = AutoDiff(otro)

        ratio = self.valor / otro.valor

        return AutoDiff(ratio, 
                    (self.deriv - ratio * otro.deriv) / otro.valor)

    def __rdiv__(self, otro):
        if not isinstance(otro, AutoDiff):
            otro = AutoDiff(otro)

        return (otro / self)

    def sin(self):
        return AutoDiff(sin(self.valor), cos(self.valor) * self.deriv)

    def cos(self):
        return AutoDiff(cos(self.valor), -sin(self.valor) * self.deriv)



    def exp(self):
        value = math.exp(self.valor)
        return AutoDiff(value, value*self.deriv) 

    # def __getitem__(self, n):
    #     if n==0:
    #         return self.valor
    #     elif n==1:
    #         return self.deriv

    #     raise ValueError("Item number must be 0 or 1; was given %d" % n)

    def tuple(self):
        return (self.valor, self.deriv)


def differentiate(f, a):
    """Calculates the AutoDiff (value and derivative) of the funtion f(x) in the point x=a
    Returns the result as a pair (value, derivative) """

    x = AutoDiff(a, np.ones_like(a))

    result = f(x)
    return result.tuple()

# def sin(x):
#     #print "Trying sin of ", x
#     try:
#         return x.sin()
#     except:
#         return math.sin(x)

# def cos(x):
#     try:
#         return x.cos()
#     except:
#         return math.cos(x)


# def exp(x):
#     try:
#         return x.exp()
#     except:
#         return math.exp(x)

    



