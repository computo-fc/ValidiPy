# -*- coding: utf-8 -*-

from global_funcs import *
# global_funcs defines sin, cos etc. in a sensible way so that they work on any object or default to a number

ctx = gmpy2.get_context()

class Interval(object):
    """
    Represents a one-dimensional interval and sets up the basic operations
    with intervals by overloading the operators +, -, *, / and **.
    The basic elementary functions over intervals are implemented. 
    Default floats are `mpfr`, which is loaded internally, in order to use some 
    elementary functions with extended precision. Directed rounding is implemented.

    Many ideas are taken from the book *Validated Numerics* by Warwick Tucker.
    """

    def __init__(self, a, b=None):
        """ 
        We define an interval by fixing its lower (lo) and upper (hi) limits,
        including *directed rounding*: the lower limit is rounded down, the upper limit is rounded up.

        Extended precision is implemented using gmpy2

        Attributes
        ----------
        a : lower limit; must be given
        b : upper limit; if not given, b = a.

        Returns
        -------
        An interval: self.lo and self.hi

        Irrespective of the magnitude of a and b, the returned value satisfies 
        self.lo <= self.hi_half

        """
        if b is None:       # single argument

            try a[0]:       # if a is a tuple, list etc.
                a, b = a    # then unpack a

            except:
                b = a       # else create a thin interval

        elif (b < a):       # limits "wrong way round"
            
            a, b = b, a

            # (Note that this is *not* the right approach for so-called "extended IA"
            # in which inverted intervals represent *excluded* intervals)



        # If using arbitrary precision, convert the types if necessary:

        if not isinstance(a, mpfr_type) and mem.libmath==gmpy2:
            # maybe more elegant to use contexts, but the code is perhaps less clear
            # with gmpy2.local_context(gmpy2.get_context(), round=RoundDown) as ctx:
            #    a = mpfr( str(a) )
            ctx.round = RoundDown
            a = mpfr( str(a) )   # convert to mpfr
        
        if not isinstance(b, mpfr_type) and mem.libmath==gmpy2:
            #with gmpy2.local_context(gmpy2.get_context(), round=RoundUp) as ctx:
            #    b = mpfr( str(b) )
            ctx.round = RoundUp
            b = mpfr( str(b) )

        self.lo, self.hi = a, b
        

    # Formatting functions:
    #
    def __repr__(self):
        return "Interval[{}, {}]".format(repr(self.lo),repr(self.hi))
    
    def __str__(self):
        return "[{},{}]".format(repr(self.lo),repr(self.hi))

    # Special formatting functions for the IPython Notebook:
    #
    def _repr_html_(self):
        reprn = "Interval[{}, {}]".format(repr(self.lo),repr(self.hi))
        reprn = reprn.replace("inf", r"&infin;")
        return reprn

    def _repr_latex_(self):
        return "$[{}, {}]$".format(repr(self.lo),repr(self.hi))


    # Arithmetic operations etc.
    #
    def __add__(self, other):
        """Sum of intervals"""

        other = self.make_interval(other)

        ctx.round = RoundDown
        lower = self.lo + other.lo
        
        ctx.round = RoundUp
        upper = self.hi + other.hi
        
        return Interval( lower, upper )
        
    def __radd__(self, other):
        return self + other


    def __sub__(self, other):
        """Difference of intervals"""

        other = self.make_interval(other)

        ctx.round = RoundDown
        lower = self.lo - other.hi
        
        ctx.round = RoundUp
        upper = self.hi - other.lo
        
        return Interval( lower, upper )
                
    def __rsub__(self, other):
        return -(self - other)


    def __pos__(self):
        return self

    def __neg__(self):
        """
        El negativo de un Interval
        """
        return Interval( -self.hi, -self.lo )
        

    def __mul__(self, other):
        """
        Se implementa la multiplicación usando `multFast`
        """

        other = self.make_interval(other)
        return self._mult3(other)
        
    def __rmul__(self, other):
        return self * other


    def _mult1(self,other):
        """ Algorítmo de la multiplicación ingenuo """
        S = [ self.lo*other.lo, self.lo*other.hi, 
              self.hi*other.lo, self.hi*other.hi ]
        return Interval( min(S), max(S) )

    def _mult2(self,other):
        """
        Algorítmo de la multiplicación que distingue los nueve casos posibles
        """
        if (self.lo >= 0.0 and other.lo >= 0.0):
            return Interval( self.lo*other.lo, self.hi*other.hi )
        elif (self.hi < 0.0 and other.hi < 0.0):
            return Interval( self.hi*other.hi, self.lo*other.lo )
        elif (self.lo >= 0.0 and other.hi < 0.0):
            return Interval( self.lo*other.hi, self.hi*other.lo )
        elif (self.hi < 0.0 and other.lo >= 0.0):
            return Interval( self.hi*other.lo, self.lo*other.hi )
        elif (self.lo >= 0.0 and other.lo*other.hi < 0.0):
            return Interval( self.hi*other.lo, self.hi*other.hi )
        elif (self.hi < 0.0 and other.lo*other.hi < 0.0):
            return Interval( self.lo*other.hi, self.lo*other.lo )
        elif (other.lo >= 0.0 and self.lo*self.hi < 0.0):
            return Interval( self.lo*other.hi, self.hi*other.hi )
        elif (other.hi < 0.0 and self.lo*self.hi < 0.0):
            return Interval( self.hi*other.lo, self.lo*other.lo )

        else:
            #(self.lo*self.hi < 0.0 and other.lo*other.hi < 0.0):
            S1 = [ self.lo*other.lo, self.hi*other.hi ]
            S2 = [ self.hi*other.lo, self.lo*other.hi ]
            return Interval( min(S2), max(S1) )

    def _mult3(self,other):
        """ Algorítmo de la multiplicación ingenuo; incorpora el redonde """

        ctx.round = RoundDown
        S_lower = [ self.lo*other.lo, self.lo*other.hi, self.hi*other.lo, self.hi*other.hi ]

        ctx.round = RoundUp
        S_upper = [ self.lo*other.lo, self.lo*other.hi, self.hi*other.lo, self.hi*other.hi ]

        S1, S2 = min( S_lower ), max( S_upper )
        return Interval( S1, S2 )


    def __div__(self, other):
        """
        División de Intervals: producto del primero por el recíproco del segundo
        """
        other = self.make_interval(other)

        try:
            return self * other.reciprocal()
        except ZeroDivisionError:
            print ("To divide by an interval containining 0,"
                   " we need to implement extended intervals!" )

    def __rdiv__(self, other):
        # Esto se encarga de cosas tipo numero/Interval; self es el Interval
        return (self / other).reciprocal()


    def reciprocal(self):
        """
        Esto define el recíproco de un Interval
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
            
        return Interval( lower, upper )


    def __contains__(self, x):
        """
        Esto verifica si el Interval contiene (o no) un número real;
        implementa al operador `in`
        """

        if isinstance(x, Interval):
            return self.lo <= x.lo and x.hi <= self.hi

        return self.lo <= x <= self.hi


    def strictly_contains(self, x):

        if isinstance(x, Interval):
            return self.lo < x.lo and x.hi < self.hi

        return self.lo < x < self.hi
        

    # pow, rpow, abs, sin, cos, ...
    def exp(self):
        """
        Exponencial de un Interval: 'self.exp()'
        """
        ctx.round = RoundDown
        lower = exp( self.lo )
        
        ctx.round = RoundUp
        upper = exp( self.hi )

        return Interval( lower, upper )


    def log(self):
        """
        Logaritmo de un Interval: 'self.log()'

        NOTA: Si el Interval contiene al 0, pero no es estrictamente negativo,
        se calcula el logaritmo de la intersección del Interval con el dominio
        natural del logaritmo, i.e., [0,+inf].
        """
        
        inf = calc_inf()
        domainNatural = Interval(0, inf)
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

            return Interval( lower, upper )

        elif 0 > self.hi:
            print ( "Interval {} < 0\nlog(x) cannot be computed "
                    "for negative numbers.".format(self) )

        else:
            ctx.round = RoundDown
            lower = log( self.lo )
            
            ctx.round = RoundUp
            upper = log( self.hi )

            return Interval( lower, upper )

    def sqrt(self):
        """An implementation of the square root"""
        inf = calc_inf()
        domainNatural = Interval(0, inf)

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

            return Interval( lower, upper )

        elif 0 > self.hi:
            print ("Interval {} < 0\nsqrt(x) cannot be computed "
                   "for negative numbers.".format(self) )
            pass

        else:
            ctx.round = RoundDown
            lower = sqrt( self.lo )
            
            ctx.round = RoundUp
            upper = sqrt( self.hi )

            return Interval( lower, upper )



    def __pow__(self, exponent):
        """
        Se calcula la potencia de un Interval; operador '**'
        """
        inf = calc_inf()

        if isinstance( exponent, Interval ):
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
                    return Interval( lower, upper )

                else:
                    # odd exponent
                    ctx.round = RoundDown
                    lower = self.lo**exponent
                    
                    ctx.round = RoundUp
                    upper = self.hi**exponent
                    return Interval( lower, upper )

            else:
                # exponent < 0
                result = self**(-exponent)
                return result.reciprocal()

        else:
            # exponent is a generic float
            domainNatural = Interval(0, inf)
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
                    
                    return Interval( lower, upper )

                elif 0 > self:
                    print ("Negative intervals can not be raised to a fractional power")

                else:
                    ctx.round = RoundDown
                    lower = min( self.lo**exponent, self.hi**exponent )
                    
                    ctx.round = RoundUp
                    upper = max( self.lo**exponent, self.hi**exponent )
                    
                    return Interval( lower, upper )

            else:
                result = self**(-exponent)
                return result.reciprocal()

    def __rpow__(self,exponent):
        return Interval(exponent)**self


    def sin(self):
        """
        Se calcula el seno de un Interval
        """
        pi = calc_pi()
        pi_half = 0.5 * pi
        dospi = 2.0 * pi
        xlow, xhig = self.lo, self.hi
        whole_range = Interval(-1.0,1.0)

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

                return Interval( sin_xlo, sin_xhi )
                
            else:
                return whole_range

        else:
            if (( lo_quarter == 3 and hi_quarter==0 ) or 
                ( lo_quarter == 1 and hi_quarter==2 )):
                
                ctx.round = RoundDown
                sin_xlo = sin( xlow )
                
                ctx.round = RoundUp
                sin_xhi = sin( xhig )
                
                return Interval( sin_xlo, sin_xhi )
            
            elif (( lo_quarter == 0 or lo_quarter==3 ) and 
                  ( hi_quarter==1 or hi_quarter==2 )):
                
                ctx.round = RoundDown
                sin_xlo = sin( xlow )
                sin_xhi = sin( xhig )
                min_sin = min( sin_xlo, sin_xhi )
                
                return Interval( min_sin, 1.0 )
            
            elif (( lo_quarter == 1 or lo_quarter==2 ) and 
                  ( hi_quarter==3 or hi_quarter==0 )):
                
                ctx.round = RoundUp
                sin_xlo = sin( xlow )
                sin_xhi = sin( xhig )
                max_sin = max( sin_xlo, sin_xhi )
              
                return Interval( -1.0, max_sin )
            
            elif (( lo_quarter == 0 and hi_quarter==3 ) or 
                  ( lo_quarter == 2 and hi_quarter==1 )):
                return whole_range
            
            else:
                # This should be never reached!
                raise NotImplementedError( ('SOMETHING WENT WRONG.\n'
                    'This should have never been reached') )


    def cos(self):
        """
        Se calcula el coseno de un Interval
        """
        pi = calc_pi()
        pi_half = 0.5 * pi
        dospi = 2.0 * pi
        xlow, xhig = self.lo, self.hi
        whole_range = Interval(-1.0,1.0)

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
                
                return Interval( cos_xhi, cos_xlo )
            else:
                return whole_range

        else:
            if (( lo_quarter == 2 and hi_quarter==3 ) or 
                ( lo_quarter == 0 and hi_quarter==1 )):
                
                ctx.round = RoundDown
                cos_xhi = cos( xhig )
                
                ctx.round = RoundUp
                cos_xlo = cos( xlow )
                
                return Interval( cos_xhi, cos_xlo )
            
            elif (( lo_quarter == 2 or lo_quarter==3 ) and 
                  ( hi_quarter==0 or hi_quarter==1 )):
               
                ctx.round = RoundDown
                cos_xlo = cos( xlow )
                cos_xhi = cos( xhig )
                min_cos = min( cos_xlo, cos_xhi )
               
                return Interval( min_cos, 1.0 )
            
            elif (( lo_quarter == 0 or lo_quarter==1 ) and 
                  ( hi_quarter==2 or hi_quarter==3 )):
            
                ctx.round = RoundUp
                cos_xlo = cos( xlow )
                cos_xhi = cos( xhig )
                max_cos = max( cos_xlo, cos_xhi )
            
                return Interval( -1.0, max_cos )
            
            elif (( lo_quarter == 3 and hi_quarter==2 ) or 
                  ( lo_quarter == 1 and hi_quarter==0 )):
                return whole_range
            
            else:
                raise NotImplementedError( ('SOMETHING WENT WRONG.\n'
                    'This should have never been reached') )


    def tan(self):
        """
        Se calcula la tangente de un Interval
        """
        inf = calc_inf()
        pi = calc_pi()
        pi_half = 0.5 * pi
        xlow, xhig = self.lo, self.hi
        whole_range = Interval(-inf,inf)

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
            return Interval( tan_xlo, tan_xhi )

        else:

            disjoint_interval2 = Interval( tan_xlo, inf )
            disjoint_interval1 = Interval( -inf, tan_xhi )

            print ("\n The resulting interval is disjoint:\n {}U{}\n"
                   "\n We consider the hull of the disjoint subintervals:\n "
                   "{}\n".format(disjoint_interval1, disjoint_interval2, whole_range))

            return whole_range


    # Las relaciones que sirven para checar el orden parcial
    def __eq__(self, other):
        """
        Aquí se checa la igualdad de dos Intervals; operador '=='
        """
        try:
            return self.lo == other.lo and self.hi == other.hi
        except:
            return self == Interval(other)

    def __ne__(self, other):
        """
        Aquí se checa la NO igualdad de dos Intervals; operador '!='
        """
        return not self == other

    def __le__(self, other):
        """
        Se checa el ordenamiento de los Intervals (ver Tucker); 
        operador '<='
        """
        try:
            return self.lo <= other.lo and self.hi <= other.hi
        except:
            return self <= Interval(other)

    def __rle__(self, other):
        return self >= other

    def __ge__(self, other):
        """ 
        Se checa el ordenamiento de los Intervals (ver Tucker); 
        operador '>='
        """
        try:
            return self.lo >= other.lo and self.hi >= other.hi
        except:
            return self >= Interval(other)

    def __rge__(self, other):
        return self <= other

    def __lt__(self, other):
        """
        Se checa el ordenamiento de los Intervals (ver Tucker); 
        operador '<'
        """
        try:
            return self.lo < other.lo and self.hi < other.hi
        except:
            return self < Interval(other)

    def __rlt__(self, other):
        return self > other

    def __gt__(self, other):
        """
        Se checa el ordenamiento de los Intervals (ver Tucker); 
        operador '>'
        """
        try:
            return self.lo > other.lo and self.hi > other.hi
        except:
            return self > Interval(other)

    def __rgt__(self, other):
        return self < other

    # Las operaciones entre Intervals vistas como conjuntos 
    def _is_empty_intersection(self, other):
        """Verifica si la intersección de los Intervals es vacía"""
        return self.hi < other.lo or other.hi < self.lo

    def intersection(self, other):
        """
        Intersección de Intervals
        """
        other = self.make_interval(other)

        if self._is_empty_intersection(other):
            print ("Intersection is empty:\n" 
                   "Intervals {} and {} are disjoint".format(self,other) )
            return

        else:
            return Interval( max(self.lo,other.lo), min(self.hi,other.hi) )

    def hull(self, other):
        """Envoltura/casco de dos Intervals"""
        lower = min(self.lo,other.lo)
        upper = max(self.hi,other.hi)

        return Interval( lower, upper )

    def union(self, other):
        """Unión de Intervals"""
        other = self.make_interval(other)

        if self._is_empty_intersection(other):
            print ("Union yields no connected interval:\n" 
                   "Intervals {} and {} are disjoint".format(self,other) )
        else:
            return self.hull(other)

    # Algunas funciones escalares de Intervals (ver Tucker)
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
        Esto define la función abs() de un Interval, cuyo resultado
        es un Interval (ver Tucker).

        NOTA: La función que regresa la máxima distancia al origen es `self.mag()`
        (magnitud) y la que regresa la mínima distancia es `self.mig()` ('mignitud').
        """
        return Interval( self.mig(), self.mag() )

    def abs(self):
        return abs(self)

    def dist(self, other):
        """
        Esto define la distancia de Hausdorff entre dos Intervals; ver Tucker.
        """
        return max( abs(self.lo-other.lo), abs(self.hi-other.hi) )


    # convert a number, tuple etc. to an interval if necessary:
    def make_interval(self, a):
        if isinstance(a, Interval):
            return a

        return Interval(a)


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

    return Interval( num1a, num2a )

def split_interval( x, num_divisions=1 ):
    """
    Divide un Interval en n=num_divisions Intervals iguales
    """
    num_divisions = int(num_divisions)
    if num_divisions < 1:
        num_divisions = 1

    edge_points = np.linspace(x.lo, x.hi, num_divisions+1)
    splited_intervals = [Interval(a, b) for (a,b) in 
                         zip(edge_points[0:num_divisions+1],edge_points[1:num_divisions+2]) ]

    return splited_intervals

def range_interval_f( fun, subdivided_interval ):
    """
    Evalua la función f(x) extendida sobre Intervals, en una lista de subIntervals
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

        # Se divide los subIntervale en 2**num subIntervals iguales
        subdivided_intervals = split_interval( x, num )
        # Se calculan las extensiones de la función sobre el Interval, usando los subIntervals
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
