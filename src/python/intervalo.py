# -*- coding: utf-8 -*-

from global_funcs import *

ctx = gmpy2.get_context()

class Intervalo(object):
    """
    Class 'Intervalo': defines a one dimensional interval and sets up the basic operations
    with intervals by overloading the operators +, -, *, / and **, and implements also
    some elementary functions over intervals. Default floats are `mpf`, which is charged
    internally, in order to use some elementary functions with extended precision and
    implement (in the future) directed rounding.
    """

    def __init__(self, a, b=None):
        """ 
        We define an interval by fixing its lower (lo) and upper (hi) limits,
        including directed rounding, i.e., the lower limit is rounded down, 
        and the upper one is rounded up. 

        Extended precision is implemented using gmpy2

        Attributes
        ----------
        a : lower limit; it must be provided
        b : upper limit; if it is not provided, b = a.

        Returns
        -------
        An interval: self.lo and self.hi

        Irrespective of the magnitude of a and b, the returned value satisfies 
        self.lo <= self.hi_half

        """
        if b is None:
            # single argument makes thin interval
            b = a

        elif (b < a):
            # limits wrong way round; not right approach for extended IA
            a, b = b, a

        # This implements extended precision and directed rounding
        if not isinstance(a, mpfr_type) and mem.libmath==gmpy2:
            #with gmpy2.local_context(gmpy2.get_context(), round=RoundDown) as ctx:
            #    a = mpfr( str(a) )
            ctx.round = RoundDown
            a = mpfr( str(a) )
        
        if not isinstance(b, mpfr_type) and mem.libmath==gmpy2:
            #with gmpy2.local_context(gmpy2.get_context(), round=RoundUp) as ctx:
            #    b = mpfr( str(b) )
            ctx.round = RoundUp
            b = mpfr( str(b) )

        self.lo, self.hi = a, b
        

    # Lo siguiente sirve para dar información bonita del objeto `Intervalo`
    #
    def __repr__(self):
        return "Intervalo[{}, {}]".format(repr(self.lo),repr(self.hi))
    
    def __str__(self):
        return "[{},{}]".format(repr(self.lo),repr(self.hi))

    # Representaciones especiales para el IPython Notebook:
    def _repr_html_(self):
        reprn = "Intervalo[{}, {}]".format(repr(self.lo),repr(self.hi))
        reprn = reprn.replace("inf", r"&infin;")
        return reprn

    def _repr_latex_(self):
        return "$[{}, {}]$".format(repr(self.lo),repr(self.hi))


    # Aquí vienen las operaciones aritméticas y varias funciones
    #
    def __add__(self, otro):
        """
        Suma de intervalos
        """
        otro = self.make_interval(otro)

        ctx.round = RoundDown
        lower = self.lo + otro.lo
        
        ctx.round = RoundUp
        upper = self.hi + otro.hi
        
        return Intervalo( lower, upper )
        
    def __radd__(self, otro):
        return self + otro


    def __sub__(self, otro):
        """
        Resta de intervalos
        """
        otro = self.make_interval(otro)

        ctx.round = RoundDown
        lower = self.lo - otro.hi
        
        ctx.round = RoundUp
        upper = self.hi - otro.lo
        
        return Intervalo( lower, upper )
                
    def __rsub__(self, otro):
        return -(self - otro)


    def __pos__(self):
        return self

    def __neg__(self):
        """
        El negativo de un intervalo
        """
        return Intervalo( -self.hi, -self.lo )
        

    def __mul__(self, otro):
        """
        Se implementa la multiplicación usando `multFast`
        """

        otro = self.make_interval(otro)
        return self._mult3(otro)
        
    def __rmul__(self, otro):
        return self * otro


    def _mult1(self,otro):
        """ Algorítmo de la multiplicación ingenuo """
        S = [ self.lo*otro.lo, self.lo*otro.hi, 
              self.hi*otro.lo, self.hi*otro.hi ]
        return Intervalo( min(S), max(S) )

    def _mult2(self,otro):
        """
        Algorítmo de la multiplicación que distingue los nueve casos posibles
        """
        if (self.lo >= 0.0 and otro.lo >= 0.0):
            return Intervalo( self.lo*otro.lo, self.hi*otro.hi )
        elif (self.hi < 0.0 and otro.hi < 0.0):
            return Intervalo( self.hi*otro.hi, self.lo*otro.lo )
        elif (self.lo >= 0.0 and otro.hi < 0.0):
            return Intervalo( self.lo*otro.hi, self.hi*otro.lo )
        elif (self.hi < 0.0 and otro.lo >= 0.0):
            return Intervalo( self.hi*otro.lo, self.lo*otro.hi )
        elif (self.lo >= 0.0 and otro.lo*otro.hi < 0.0):
            return Intervalo( self.hi*otro.lo, self.hi*otro.hi )
        elif (self.hi < 0.0 and otro.lo*otro.hi < 0.0):
            return Intervalo( self.lo*otro.hi, self.lo*otro.lo )
        elif (otro.lo >= 0.0 and self.lo*self.hi < 0.0):
            return Intervalo( self.lo*otro.hi, self.hi*otro.hi )
        elif (otro.hi < 0.0 and self.lo*self.hi < 0.0):
            return Intervalo( self.hi*otro.lo, self.lo*otro.lo )

        else:
            #(self.lo*self.hi < 0.0 and otro.lo*otro.hi < 0.0):
            S1 = [ self.lo*otro.lo, self.hi*otro.hi ]
            S2 = [ self.hi*otro.lo, self.lo*otro.hi ]
            return Intervalo( min(S2), max(S1) )

    def _mult3(self,otro):
        """ Algorítmo de la multiplicación ingenuo; incorpora el redonde """

        ctx.round = RoundDown
        S_lower = [ self.lo*otro.lo, self.lo*otro.hi, self.hi*otro.lo, self.hi*otro.hi ]

        ctx.round = RoundUp
        S_upper = [ self.lo*otro.lo, self.lo*otro.hi, self.hi*otro.lo, self.hi*otro.hi ]

        S1, S2 = min( S_lower ), max( S_upper )
        return Intervalo( S1, S2 )


    def __div__(self, otro):
        """
        División de intervalos: producto del primero por el recíproco del segundo
        """
        otro = self.make_interval(otro)

        try:
            return self * otro.reciprocal()
        except ZeroDivisionError:
            print ("To divide by an interval containining 0,"
                   " we need to implement extended intervals!" )

    def __rdiv__(self, otro):
        # Esto se encarga de cosas tipo numero/intervalo; self es el intervalo
        return (self / otro).reciprocal()


    def reciprocal(self):
        """
        Esto define el recíproco de un intervalo
        """
        inf = calc_inf()

        ctx.round = RoundDown
        try:
            lower = 1.0/self.hi
        except:
            lower = -inf
            
        ctx.round = RoundUp
        try:
            upper = 1.0/self.lo
        except:
            upper = inf
            
        if self.strictly_contains( 0.0 ):
            lower = -inf
            upper = inf
            print "\nInterval {} in denominator contains 0.".format(self)
            
        return Intervalo( lower, upper )


    def __contains__(self, x):
        """
        Esto verifica si el intervalo contiene (o no) un número real;
        implementa al operador `in`
        """

        if isinstance(x, Intervalo):
            return self.lo <= x.lo and x.hi <= self.hi

        return self.lo <= x <= self.hi


    def strictly_contains(self, x):

        if isinstance(x, Intervalo):
            return self.lo < x.lo and x.hi < self.hi

        return self.lo < x < self.hi
        

    # pow, rpow, abs, sin, cos, ...
    def exp(self):
        """
        Exponencial de un intervalo: 'self.exp()'
        """
        ctx.round = RoundDown
        lower = exp( self.lo )
        
        ctx.round = RoundUp
        upper = exp( self.hi )

        return Intervalo( lower, upper )


    def log(self):
        """
        Logaritmo de un intervalo: 'self.log()'

        NOTA: Si el intervalo contiene al 0, pero no es estrictamente negativo,
        se calcula el logaritmo de la intersección del intervalo con el dominio
        natural del logaritmo, i.e., [0,+inf].
        """
        
        inf = calc_inf()
        domainNatural = Intervalo(0, inf)
        if 0 in self:
            intervalRestricted = self.intersection( domainNatural )

            print ("\nWARNING:\n"
                   "Interval {} contains 0 or negative numbers.\n".format(self) )
            
            print ("Restricting to the intersection with the natural domain of log(x),\n"
                   "i.e. {}\n".format(intervalRestricted) )
                        
            ctx.round = RoundDown
            lower = log( intervalRestricted.lo )
            
            ctx.round = RoundUp
            upper = log( intervalRestricted.hi )

            return Intervalo( lower, upper )

        elif 0 > self.hi:
            print ( "Interval {} < 0\nlog(x) cannot be computed "
                    "for negative numbers.".format(self) )

        else:
            ctx.round = RoundDown
            lower = log( self.lo )
            
            ctx.round = RoundUp
            upper = log( self.hi )

            return Intervalo( lower, upper )

    def sqrt(self):
        """An implementation of the square root"""
        inf = calc_inf()
        domainNatural = Intervalo(0, inf)

        if 0 in self:

            intervalRestricted = self.intersection( domainNatural )

            print ("\nWARNING:\n"
                   "Interval {} contains 0 or negative numbers.\n".format(self) )
            
            print ("Restricting to the intersection with the natural domain of sqrt(x),\n"
                   "i.e. {}\n".format(intervalRestricted) )
            
            ctx.round = RoundDown
            lower = sqrt( intervalRestricted.lo )
            
            ctx.round = RoundUp
            upper = sqrt( intervalRestricted.hi )

            return Intervalo( lower, upper )

        elif 0 > self.hi:
            print ("Interval {} < 0\nsqrt(x) cannot be computed "
                   "for negative numbers.".format(self) )
            pass

        else:
            ctx.round = RoundDown
            lower = sqrt( self.lo )
            
            ctx.round = RoundUp
            upper = sqrt( self.hi )

            return Intervalo( lower, upper )



    def __pow__(self, exponent):
        """
        Se calcula la potencia de un intervalo; operador '**'
        """
        inf = calc_inf()

        if isinstance( exponent, Intervalo ):
            # here exponent is an interval
            tmp_dummy = exponent * self.log()
            return tmp_dummy.exp()

        # exponent is a number (int, float, mpfr, ...)
        if exponent == int(exponent):
            # exponent is an integer
            if exponent >= 0:
                if exponent%2 == 0:
                    # even exponent
                    ctx.round = RoundDown
                    lower = (self.mig())**exponent
                    
                    ctx.round = RoundUp
                    upper = (self.mag())**exponent
                    return Intervalo( lower, upper )

                else:
                    # odd exponent
                    ctx.round = RoundDown
                    lower = self.lo**exponent
                    
                    ctx.round = RoundUp
                    upper = self.hi**exponent
                    return Intervalo( lower, upper )

            else:
                # exponent < 0
                result = self**(-exponent)
                return result.reciprocal()

        else:
            # exponent is a generic float
            domainNatural = Intervalo(0, inf)
            if exponent >= 0:
                if 0 in self:
                    intervalRestricted = self.intersection( domainNatural )

                    print "\nWARNING: Interval {} contains 0.\n".format(self)

                    print ("Restricting to the intersection to the natural domain of **, "
                           "i.e. {}\n".format(intervalRestricted) )

                    ctx.round = RoundDown
                    lower = mpfr( '0.0' )**exponent
                    
                    ctx.round = RoundUp
                    upper = elf.hi**exponent
                    
                    return Intervalo( lower, upper )

                elif 0 > self:
                    print ("Negative intervals can not be raised to a fractional power")

                else:
                    ctx.round = RoundDown
                    lower = min( self.lo**exponent, self.hi**exponent )
                    
                    ctx.round = RoundUp
                    upper = max( self.lo**exponent, self.hi**exponent )
                    
                    return Intervalo( lower, upper )

            else:
                result = self**(-exponent)
                return result.reciprocal()

    def __rpow__(self,exponent):
        return Intervalo(exponent)**self


    def sin(self):
        """
        Se calcula el seno de un intervalo
        """
        pi = calc_pi()
        pi_half = 0.5 * pi
        dospi = 2.0 * pi
        xlow, xhig = self.lo, self.hi
        whole_range = Intervalo(-1.0,1.0)

        # Check the specific case:
        if xhig > xlow + dospi:
            return whole_range
        
        # within 1 full period of sin(x); 20 cases
        # some abreviations
        lo_mod2pi = xlow % dospi
        hi_mod2pi = xhig % dospi
        lo_quarter = mem.libmath.floor( lo_mod2pi / pi_half )
        hi_quarter = mem.libmath.floor( hi_mod2pi / pi_half )
        
        if lo_quarter == hi_quarter:
            if lo_mod2pi <= hi_mod2pi:
                
                ctx.round = RoundDown
                sin_xlo = sin( xlow )
                
                ctx.round = RoundUp
                sin_xhi = sin( xhig )

                return Intervalo( sin_xlo, sin_xhi )
                
            else:
                return whole_range

        else:
            if (( lo_quarter == 3 and hi_quarter==0 ) or 
                ( lo_quarter == 1 and hi_quarter==2 )):
                
                ctx.round = RoundDown
                sin_xlo = sin( xlow )
                
                ctx.round = RoundUp
                sin_xhi = sin( xhig )
                
                return Intervalo( sin_xlo, sin_xhi )
            
            elif (( lo_quarter == 0 or lo_quarter==3 ) and 
                  ( hi_quarter==1 or hi_quarter==2 )):
                
                ctx.round = RoundDown
                sin_xlo = sin( xlow )
                sin_xhi = sin( xhig )
                min_sin = min( sin_xlo, sin_xhi )
                
                return Intervalo( min_sin, 1.0 )
            
            elif (( lo_quarter == 1 or lo_quarter==2 ) and 
                  ( hi_quarter==3 or hi_quarter==0 )):
                
                ctx.round = RoundUp
                sin_xlo = sin( xlow )
                sin_xhi = sin( xhig )
                max_sin = max( sin_xlo, sin_xhi )
              
                return Intervalo( -1.0, max_sin )
            
            elif (( lo_quarter == 0 and hi_quarter==3 ) or 
                  ( lo_quarter == 2 and hi_quarter==1 )):
                return whole_range
            
            else:
                # This should be never reached!
                raise NotImplementedError( ('SOMETHING WENT WRONG.\n'
                    'This should have never been reached') )


    def cos(self):
        """
        Se calcula el coseno de un intervalo
        """
        pi = calc_pi()
        pi_half = 0.5 * pi
        dospi = 2.0 * pi
        xlow, xhig = self.lo, self.hi
        whole_range = Intervalo(-1.0,1.0)

        # Check the specific case:
        if xhig > xlow + dospi:
            return whole_range
        
        # within 1 full period of sin(x); 20 cases
        # some abreviations
        lo_mod2pi = xlow % dospi
        hi_mod2pi = xhig % dospi
        lo_quarter = mem.libmath.floor( lo_mod2pi / pi_half )
        hi_quarter = mem.libmath.floor( hi_mod2pi / pi_half )
        
        if lo_quarter == hi_quarter:
            if lo_mod2pi <= hi_mod2pi:

                ctx.round = RoundDown
                cos_xhi = cos( xhig )
                
                ctx.round = RoundUp
                cos_xlo = cos( xlow )
                
                return Intervalo( cos_xhi, cos_xlo )
            else:
                return whole_range

        else:
            if (( lo_quarter == 2 and hi_quarter==3 ) or 
                ( lo_quarter == 0 and hi_quarter==1 )):
                
                ctx.round = RoundDown
                cos_xhi = cos( xhig )
                
                ctx.round = RoundUp
                cos_xlo = cos( xlow )
                
                return Intervalo( cos_xhi, cos_xlo )
            
            elif (( lo_quarter == 2 or lo_quarter==3 ) and 
                  ( hi_quarter==0 or hi_quarter==1 )):
               
                ctx.round = RoundDown
                cos_xlo = cos( xlow )
                cos_xhi = cos( xhig )
                min_cos = min( cos_xlo, cos_xhi )
               
                return Intervalo( min_cos, 1.0 )
            
            elif (( lo_quarter == 0 or lo_quarter==1 ) and 
                  ( hi_quarter==2 or hi_quarter==3 )):
            
                ctx.round = RoundUp
                cos_xlo = cos( xlow )
                cos_xhi = cos( xhig )
                max_cos = max( cos_xlo, cos_xhi )
            
                return Intervalo( -1.0, max_cos )
            
            elif (( lo_quarter == 3 and hi_quarter==2 ) or 
                  ( lo_quarter == 1 and hi_quarter==0 )):
                return whole_range
            
            else:
                raise NotImplementedError( ('SOMETHING WENT WRONG.\n'
                    'This should have never been reached') )


    def tan(self):
        """
        Se calcula la tangente de un intervalo
        """
        inf = calc_inf()
        pi = calc_pi()
        pi_half = 0.5 * pi
        xlow, xhig = self.lo, self.hi
        whole_range = Intervalo(-inf,inf)

        # Check the specific case:
        if xhig > xlow + pi:
            return whole_range

        # within 1 full period of tan(x) --> 4 cases
        # some abreviations
        lo_modpi = xlow % pi
        hi_modpi = xhig % pi
        lo_half = mem.libmath.floor( lo_modpi / pi_half )
        hi_half = mem.libmath.floor( hi_modpi / pi_half )
       
        ctx.round = RoundDown
        tan_xlo = tan( xlow )
        
        ctx.round = RoundUp
        tan_xhi = tan( xhig )

        if lo_half > hi_half or ( lo_half == hi_half and lo_modpi < hi_modpi):
            return Intervalo( tan_xlo, tan_xhi )

        else:

            disjoint_interval2 = Intervalo( tan_xlo, inf )
            disjoint_interval1 = Intervalo( -inf, tan_xhi )

            print ("\n The resulting interval is disjoint:\n {}U{}\n"
                   "\n We consider the hull of the disjoint subintervals:\n "
                   "{}\n".format(disjoint_interval1, disjoint_interval2, whole_range))

            return whole_range


    # Las relaciones que sirven para checar el orden parcial
    def __eq__(self, otro):
        """
        Aquí se checa la igualdad de dos intervalos; operador '=='
        """
        try:
            return self.lo == otro.lo and self.hi == otro.hi
        except:
            return self == Intervalo(otro)

    def __ne__(self, otro):
        """
        Aquí se checa la NO igualdad de dos intervalos; operador '!='
        """
        return not self == otro

    def __le__(self, otro):
        """
        Se checa el ordenamiento de los intervalos (ver Tucker); 
        operador '<='
        """
        try:
            return self.lo <= otro.lo and self.hi <= otro.hi
        except:
            return self <= Intervalo(otro)

    def __rle__(self, otro):
        return self >= otro

    def __ge__(self, otro):
        """ 
        Se checa el ordenamiento de los intervalos (ver Tucker); 
        operador '>='
        """
        try:
            return self.lo >= otro.lo and self.hi >= otro.hi
        except:
            return self >= Intervalo(otro)

    def __rge__(self, otro):
        return self <= otro

    def __lt__(self, otro):
        """
        Se checa el ordenamiento de los intervalos (ver Tucker); 
        operador '<'
        """
        try:
            return self.lo < otro.lo and self.hi < otro.hi
        except:
            return self < Intervalo(otro)

    def __rlt__(self, otro):
        return self > otro

    def __gt__(self, otro):
        """
        Se checa el ordenamiento de los intervalos (ver Tucker); 
        operador '>'
        """
        try:
            return self.lo > otro.lo and self.hi > otro.hi
        except:
            return self > Intervalo(otro)

    def __rgt__(self, otro):
        return self < otro

    # Las operaciones entre intervalos vistas como conjuntos 
    def _is_empty_intersection(self, otro):
        """Verifica si la intersección de los intervalos es vacía"""
        return self.hi < otro.lo or otro.hi < self.lo

    def intersection(self, otro):
        """
        Intersección de intervalos
        """
        otro = self.make_interval(otro)

        if self._is_empty_intersection(otro):
            print ("Intersection is empty:\n" 
                   "Intervals {} and {} are disjoint".format(self,otro) )
            return

        else:
            return Intervalo( max(self.lo,otro.lo), min(self.hi,otro.hi) )

    def hull(self, otro):
        """Envoltura/casco de dos intervalos"""
        lower = min(self.lo,otro.lo)
        upper = max(self.hi,otro.hi)

        return Intervalo( lower, upper )

    def union(self, otro):
        """Unión de intervalos"""
        otro = self.make_interval(otro)

        if self._is_empty_intersection(otro):
            print ("Union yields no connected interval:\n" 
                   "Intervals {} and {} are disjoint".format(self,otro) )
        else:
            return self.hull(otro)

    # Algunas funciones escalares de intervalos (ver Tucker)
    def diam(self):
        return self.hi - self.lo

    def rad(self):
        return 0.5*self.diam()

    def mid(self):
        return 0.5*( self.lo + self.hi )

    def mag(self):
        """Distancia máxima (magnitude) al origen"""
        return max( abs(self.lo), abs(self.hi) )

    def mig(self):
        """Distancia mínima (mignitude) al origen"""
        if 0 in self:
            return 0
        else:
            return min( abs(self.lo), abs(self.hi) )

    def __abs__(self):
        """
        Esto define la función abs() de un intervalo, cuyo resultado
        es un intervalo (ver Tucker).

        NOTA: La función que regresa la máxima distancia al origen es `self.mag()`
        (magnitud) y la que regresa la mínima distancia es `self.mig()` ('mignitud').
        """
        return Intervalo( self.mig(), self.mag() )

    def abs(self):
        return abs(self)

    def dist(self, otro):
        """
        Esto define la distancia de Hausdorff entre dos intervalos; ver Tucker.
        """
        return max( abs(self.lo-otro.lo), abs(self.hi-otro.hi) )

    def make_interval(self, a):
        if isinstance(a, Intervalo):
            return a

        return Intervalo(a)


# Funciones extras

#def make_mpfr(a, rounding=None):
#    """This creates a mpfr-number with specified rounding"""
#
#    if rounding == "RoundToNearest" or rounding==0 or rounding is None:
#        with gmpy2.local_context(gmpy2.get_context(), round=0) as ctx:
#            a = mpfr( str(a) )
#    elif rounding == "RoundUp" or rounding==2:
#        with gmpy2.local_context(gmpy2.get_context(), round=2) as ctx:
#            a = mpfr( str(a) )
#    elif rounding == "RoundDown" or rounding==3:
#        with gmpy2.local_context(gmpy2.get_context(), round=3) as ctx:
#            a = mpfr( str(a) )
#    elif rounding == "RoundToZero" or rounding==1:
#        with gmpy2.local_context(gmpy2.get_context(), round=1) as ctx:
#            a = mpfr( str(a) )
#    elif rounding == "RoundAwayZero" or rounding==4:
#        with gmpy2.local_context(gmpy2.get_context(), round=4) as ctx:
#            a = mpfr( str(a) )
#    else:
#        raise ValueError("No such rounding mode; setting it to RoundToNearest")
#
#    return a


def random_interval( infimum=-10.0, supremum=10.0 ):
    num1a = np.random.uniform( infimum, supremum )
    num2a = np.random.uniform( infimum, supremum )

    return Intervalo( num1a, num2a )

def split_interval( x, num_divisions=1 ):
    """
    Divide un intervalo en n=num_divisions intervalos iguales
    """
    num_divisions = int(num_divisions)
    if num_divisions < 1:
        num_divisions = 1

    edge_points = np.linspace(x.lo, x.hi, num_divisions+1)
    splited_intervals = [Intervalo(a, b) for (a,b) in 
                         zip(edge_points[0:num_divisions+1],edge_points[1:num_divisions+2]) ]

    return splited_intervals

def range_interval_f( fun, subdivided_interval ):
    """
    Evalua la función f(x) extendida sobre intervalos, en una lista de subintervalos
    y regresa el hull de todos ellos, es decir, una cota del rango de la función
    """
    if not isinstance( subdivided_interval, list ):
        subdivided_interval = [ subdivided_interval ]

    range_fun = [ fun(i) for i in subdivided_interval ]
    range_tot = range_fun[0]

    for i in range_fun[1:]:
        range_tot = range_tot.hull(i)

    return range_tot

def plot_interval_f( fun, x, pow2=0, num_points=101 ):
    """
    This plots the interval extension of a function `fun` over the interval `x`,
    which is diveded in num=1,2,4,...,2**pow2 uniform subintervals.
    """
    num_intervals = [ 2**p for p in range(pow2+1) ]
    plt.figure()
    plt.subplot(1, 1, 1)

    for num in num_intervals:
        fact_alfa = num*1.0/num_intervals[-1]   # for plotting

        # Se divide los subintervaloe en 2**num subintervalos iguales
        subdivided_intervals = split_interval( x, num )
        # Se calculan las extensiones de la función sobre el intervalo, usando los subintervalos
        rango_total = range_interval_f( fun, subdivided_intervals )
        print "Rango_tot (N={}) = {}".format(num,rango_total)
        
        # Hago el dibujo
        for x1 in subdivided_intervals:
            low = float(x1.lo)
            high = float(x1.hi)
            Ffun = fun(x1)
            xa1 = np.array([low, low, high, high])
            ya1 = np.array([float(Ffun.lo), float(Ffun.hi), float(Ffun.hi), float(Ffun.lo) ])
            plt.fill( xa1, ya1, 'b', alpha=fact_alfa )
    
    low = float(x.lo)
    high = float(x.hi)
    xx = np.linspace(low,high,num_points)
    yy = fun(xx)
    plt.plot( xx, yy, 'red')
    return 

