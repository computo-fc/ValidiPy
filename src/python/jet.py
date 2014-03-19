# -*- coding:utf-8 -*-

"""The x in the argument of functions is actually an object which is itself a derivative pair"""


import numpy as math  # NB
import numpy as np


class Jet:
    """Multidim jet"""
    def __init__(self, valor, deriv=0.0):
        self.valor = valor
        self.deriv = deriv

    def __add__(self, otro):
        if not isinstance(otro, Jet):
            otro = Jet(otro)

        return Jet(self.valor+otro.valor, self.deriv+otro.deriv)

    def __radd__(self, otro):
        return self+otro

    def __sub__(self, otro):
        if not isinstance(otro, Jet):
            otro = Jet(otro)

        return Jet(self.valor-otro.valor, self.deriv-otro.deriv)

    def __neg__(self, otro):
        return Jet(-self.valor, -self.otro)

    def __rsub__(self, otro):
        return (self - otro).neg()

    def __pow__(self, n):
        if not isinstance(n, int):
            raise ValueError("Powers not implemented for non-integers")

        return Jet(self.valor**n, n*self.valor**(n-1) * self.deriv)

    def __repr__(self):
        return "Jet({}, {})".format(self.valor, self.deriv)

    def __mul__(self, otro):
        if not isinstance(otro, Jet):
            otro = Jet(otro)

        return Jet(self.valor*otro.valor, 
                    self.valor*otro.deriv + self.deriv*otro.valor)

    def __rmul__(self, otro):
        return self*otro

    def __div__(self, otro):
        if not isinstance(otro, Jet):
            otro = Jet(otro)

        ratio = self.valor / otro.valor

        return Jet(ratio, 
                    (self.deriv - ratio * otro.deriv) / otro.valor)

    def __rdiv__(self, otro):
        if not isinstance(otro, Jet):
            otro = Jet(otro)

        return (otro / self)

    def sin(self):
        return Jet(sin(self.valor), cos(self.valor) * self.deriv)

    def cos(self):
        return Jet(cos(self.valor), -sin(self.valor) * self.deriv)



    def exp(self):
        value = math.exp(self.valor)
        return Jet(value, value*self.deriv) 

    # def __getitem__(self, n):
    #     if n==0:
    #         return self.valor
    #     elif n==1:
    #         return self.deriv

    #     raise ValueError("Item number must be 0 or 1; was given %d" % n)

    def tuple(self):
        return (self.valor, self.deriv)


def derivar(f, a):
    """Calcular el jet (valor y derivada) de la funcion f(x) en el punto x=a"""

    x = Jet(a, np.ones_like(a))

    result = f(x)
    return result.tuple()

def sin(x):
    #print "Trying sin of ", x
    try:
        return x.sin()
    except:
        return math.sin(x)

def cos(x):
    try:
        return x.cos()
    except:
        return math.cos(x)


def exp(x):
    try:
        return x.exp()
    except:
        return math.exp(x)

    



