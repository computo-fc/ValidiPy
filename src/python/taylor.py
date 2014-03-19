import numpy as np
from math import factorial

from global_funcs import *


class Taylor:
    """Taylor expansion"""

    def __init__(self, coeffs, order=None):

        try:
            _ = coeffs[0]   # checar si es una lista
        except:
            coeffs = [coeffs]


        if order is None:
            order = len(coeffs) - 1


        self.order = order
        self.coeffs = [0] * (self.order+1)   # lista de todos ceros
        # no utilizar array de numpy todavia para poder meter intervalos etc

        self.coeffs[:len(coeffs)] = coeffs

        self.coeffs = np.array(self.coeffs)

    def __repr__(self):
        return "Taylor({})".format(self.coeffs)

    def __add__(self, other):
        if not isinstance(other, Taylor):
            other = Taylor(other, order=self.order+1)

        return Taylor(self.coeffs + other.coeffs)

    def __radd__(self, other):
        return self + other

    def __sub__(self, other):
        if not isinstance(other, Taylor):
            other = Taylor(other, order=self.order+1)

        return Taylor(self.coeffs - other.coeffs)

    def __rsub__(self, other):
        return Taylor( - (other-self).coeffs )

    def __mul__(self, other):

        if not isinstance(other, Taylor):
            other = Taylor(other)

        order = max(self.order, other.order)

        # ensure both of same order:
        a = Taylor(self.coeffs, order)
        b = Taylor(other.coeffs, order)

        coeffs = [0] * (order + 1)

        for k in range(order+1):
            for i in range(k+1):
                coeffs[k] += a.coeffs[i]*b.coeffs[k-i]

        return Taylor(coeffs)

    def __rmul__(self, other):
        return self*other

    def deriv(self, n):
        return self.coeffs[n] * factorial(n)

    def exp(self):
        exp_coeffs = self.coeffs.copy()

        exp_coeffs[0] = exp(self.coeffs[0])

        for k in range(1, len(self.coeffs)):
            suma = 0.0

            for i in range(1, k+1):
                suma += i * self.coeffs[i] * exp_coeffs[k-i]

            exp_coeffs[k] = suma / k


            # exp_coeffs[k] = np.sum(np.arange(1, k+1) * self.coeffs[1:k+1] * exp_coeffs[k:0:-1])
            # exp_coeffs[k] /= k
            
            
        return Taylor(exp_coeffs)

    