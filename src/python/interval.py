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

        if isinstance(a, Interval):
            a, b = a.lo, a.hi

        elif b is None:       # single argument

            try:
              a, b = a      # try to unpack a (check if it's a tuple etc.)

            except:
                b = a       # else create a thin interval

        elif (b < a):         # limits "wrong way round"

            a, b = b, a


        # If using arbitrary precision, convert the types as necessary:

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

        ctx.round = RoundToNearest

        self.lo, self.hi = a, b


    # Formatting functions:
    #
    def __repr__(self):
        return "Interval({}, {})".format(repr(self.lo), repr(self.hi))
        #return "Interval[{}, {}]".format(repr(self.lo), repr(self.hi))


    def __str__(self):
        return "[{},{}]".format(self.lo, self.hi)

    # Special formatting functions for the IPython Notebook:
    #
    def _repr_html_(self):
        reprn = "Interval[{}, {}]".format(self.lo, self.hi)
        reprn = reprn.replace("inf", r"&infin;")
        return reprn

    def _repr_latex_(self):
        return "$[{}, {}]$".format(self.lo, self.hi)


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
        other = self.make_interval(other)
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
        other = self.make_interval(other)
        return other - self


    def __pos__(self):
        return self

    def __neg__(self):
        """
        Unary negative
        """
        return Interval( -self.hi, -self.lo )


    def __mul__(self, other):
        """
        Multiplication: use the naive version with directed rounding
        """

        other = self.make_interval(other)
        return self._mult(other)

    def __rmul__(self, other):
        other = self.make_interval(other)
        return self * other

    # def _mult1(self,other):
    #     """ Naive multiplication; does not yet have directed rounding """

    #     S = [ self.lo*other.lo, self.lo*other.hi,
    #           self.hi*other.lo, self.hi*other.hi ]
    #     return Interval( min(S), max(S) )

    # def _mult2(self,other):
    #     """
    #     Multiplication algorithm considering explicitly the nine possible cases
    #     (No directed rounding)
    #     """
    #     if (self.lo >= 0.0 and other.lo >= 0.0):
    #         return Interval( self.lo*other.lo, self.hi*other.hi )
    #     elif (self.hi < 0.0 and other.hi < 0.0):
    #         return Interval( self.hi*other.hi, self.lo*other.lo )
    #     elif (self.lo >= 0.0 and other.hi < 0.0):
    #         return Interval( self.lo*other.hi, self.hi*other.lo )
    #     elif (self.hi < 0.0 and other.lo >= 0.0):
    #         return Interval( self.hi*other.lo, self.lo*other.hi )
    #     elif (self.lo >= 0.0 and other.lo*other.hi < 0.0):
    #         return Interval( self.hi*other.lo, self.hi*other.hi )
    #     elif (self.hi < 0.0 and other.lo*other.hi < 0.0):
    #         return Interval( self.lo*other.hi, self.lo*other.lo )
    #     elif (other.lo >= 0.0 and self.lo*self.hi < 0.0):
    #         return Interval( self.lo*other.hi, self.hi*other.hi )
    #     elif (other.hi < 0.0 and self.lo*self.hi < 0.0):
    #         return Interval( self.hi*other.lo, self.lo*other.lo )

    #     else:
    #         #(self.lo*self.hi < 0.0 and other.lo*other.hi < 0.0):
    #         S1 = [ self.lo*other.lo, self.hi*other.hi ]
    #         S2 = [ self.hi*other.lo, self.lo*other.hi ]
    #         return Interval( min(S2), max(S1) )

    def _mult(self,other):
        """ Naive multiplication with direct rounding """

        ctx.round = RoundDown
        S_lower = [ self.lo*other.lo, self.lo*other.hi, self.hi*other.lo, self.hi*other.hi ]
        S1 = min(S_lower)

        ctx.round = RoundUp
        S_upper = [ self.lo*other.lo, self.lo*other.hi, self.hi*other.lo, self.hi*other.hi ]
        S2 = max(S_upper)

        return Interval( S1, S2 )


    def __div__(self, other):
        """
        Division: product of the first by the reciprocal of the second
        """
        other = self.make_interval(other)

        try:
            return self * other.reciprocal()
        except ZeroDivisionError:
            print ("To divide by an interval containining 0,"
                   " we need to implement extended intervals!" )

    def __rdiv__(self, other):
        other = self.make_interval(other)
        return other / self


    def reciprocal(self):
        """
        Reciprocal of an interval
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
        Implements the `in` operator to check if the number x lies in the interval
        Also works if x is an interval to check if subset
        """

        if isinstance(x, Interval):
            return self.lo <= x.lo and x.hi <= self.hi

        return self.lo <= x <= self.hi


    def strictly_contains(self, x):
        """ Version of `in` with strict inequalities """

        if isinstance(x, Interval):
            return self.lo < x.lo and x.hi < self.hi

        return self.lo < x < self.hi


    # Elementary mathematical functions: pow, rpow, abs, sin, cos, exp, etc.
    #
    def __pow__(self, exponent):
        """
        Power -- implements the `**` operator
        """

        if 0 > self:
            if not exponent == int(exponent):  # noninteger exponent not allowed
                print ("Negative intervals cannot be raised to a fractional power")
                return  # Should throw error

        inf = calc_inf()
        naturalDomain = Interval(0, inf)
        intervalRestricted = self.intersection( naturalDomain )

        if isinstance( exponent, Interval ):   # exponent is an interval
            if exponent.diam() == 0:
                return self**(exponent.lo)

            ctx.round = RoundDown
            lolo = intervalRestricted.lo**(exponent.lo)
            lohi = intervalRestricted.lo**(exponent.hi)
            lower = min( lolo, lohi)

            ctx.round = RoundUp
            hilo = intervalRestricted.hi**(exponent.lo)
            hihi = intervalRestricted.hi**(exponent.hi)
            upper = max( hilo, hihi)

            return( Interval( lower, upper ) )

        if exponent == int(exponent):    # exponent is an integer
            if exponent >= 0:
                if exponent%2 == 0:     # even exponent
                    ctx.round = RoundDown
                    lower = (self.mig())**exponent

                    ctx.round = RoundUp
                    upper = (self.mag())**exponent
                    return Interval( lower, upper )

                else:                # odd exponent
                    ctx.round = RoundDown
                    lower = self.lo**exponent

                    ctx.round = RoundUp
                    upper = self.hi**exponent
                    return Interval( lower, upper )

            else:   # exponent < 0
                result = self**(-exponent)
                return result.reciprocal()

        else:  # exponent is a generic float

            if exponent >= 0:
                if 0 in self:
                    print "\nWARNING: Interval {} contains 0.\n".format(self)

                    print ("Restricting to the intersection to the natural domain of **, "
                           "i.e. {}\n".format(intervalRestricted) )

                    ctx.round = RoundDown
                    lower = mpfr( '0.0' )**exponent

                    ctx.round = RoundUp
                    upper = self.hi**exponent

                    return Interval( lower, upper )

                elif 0 > self:
                    print ("Negative intervals cannot be raised to a fractional power")

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
        exponent = self.make_interval(exponent)
        return exponent**self


    def exp(self):
        """
        Exponential
        """
        ctx.round = RoundDown
        lower = exp( self.lo )

        ctx.round = RoundUp
        upper = exp( self.hi )

        return Interval( lower, upper )


    def log(self):
        """
        Logarithm

        Note: If the Interval contains 0, but is not strictly negative,
        we calculate the logarithm of the intersection of the Interval with the natural domain of
        log, i.e. [0, +inf].
        """

        inf = calc_inf()
        naturalDomain = Interval(0, inf)

        if 0 in self:
            intervalRestricted = self.intersection( naturalDomain )

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
        naturalDomain = Interval(0, inf)

        if 0 in self:

            intervalRestricted = self.intersection( naturalDomain )

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



    def sin(self):
        """
        Sine
        """

        pi = calc_pi()
        half_pi = 0.5 * pi
        two_pi = 2.0 * pi
        xlow, xhig = self.lo, self.hi
        whole_range = Interval(-1.0,1.0)

        # Check the specific case:
        if xhig > xlow + two_pi:
            return whole_range

        # within 1 full period of sin(x); 20 cases
        # some abreviations
        lo_mod2pi = xlow % two_pi
        hi_mod2pi = xhig % two_pi
        lo_quarter = mem.libmath.floor( lo_mod2pi / half_pi )
        hi_quarter = mem.libmath.floor( hi_mod2pi / half_pi )

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
        Cosine
        """
        pi = calc_pi()
        half_pi = 0.5 * pi
        two_pi = 2.0 * pi
        xlow, xhig = self.lo, self.hi
        whole_range = Interval(-1.0,1.0)

        # Check the specific case:
        if xhig > xlow + two_pi:
            return whole_range

        # within 1 full period of sin(x); 20 cases
        # some abreviations
        lo_mod2pi = xlow % two_pi
        hi_mod2pi = xhig % two_pi
        lo_quarter = mem.libmath.floor( lo_mod2pi / half_pi )
        hi_quarter = mem.libmath.floor( hi_mod2pi / half_pi )

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
        Tangent
        """
        inf = calc_inf()
        pi = calc_pi()
        half_pi = 0.5 * pi
        xlow, xhig = self.lo, self.hi
        whole_range = Interval(-inf,inf)

        # Check the specific case:
        if xhig > xlow + pi:
            return whole_range

        # within 1 full period of tan(x) --> 4 cases
        # some abreviations
        lo_modpi = xlow % pi
        hi_modpi = xhig % pi
        lo_half = mem.libmath.floor( lo_modpi / half_pi )
        hi_half = mem.libmath.floor( hi_modpi / half_pi )

        ctx.round = RoundDown
        tan_xlo = tan( xlow )

        ctx.round = RoundUp
        tan_xhi = tan( xhig )

        if lo_half > hi_half or ( lo_half == hi_half and lo_modpi <= hi_modpi):
            return Interval( tan_xlo, tan_xhi )

        else:

            disjoint_interval2 = Interval( tan_xlo, inf )
            disjoint_interval1 = Interval( -inf, tan_xhi )

            print ("\n The resulting interval is disjoint:\n {}U{}\n"
                   "\n We consider the hull of the disjoint subintervals:\n "
                   "{}\n".format(disjoint_interval1, disjoint_interval2, whole_range))

            return whole_range


    # Partial order relations
    #

    def __eq__(self, other):
        """ Equality: `==` operator """
        other = self.make_interval(other)
        return self.lo == other.lo and self.hi == other.hi

    def __ne__(self, other):
        """ Inequality; `!=` operator """
        return not self == other

    def __le__(self, other):
        """ Inequality `<=` """
        other = self.make_interval(other)
        return self.lo <= other.lo and self.hi <= other.hi

    def __rle__(self, other):
        return self >= other

    def __ge__(self, other):
        """ Inequality `>=` """
        other = self.make_interval(other)
        return self.lo >= other.lo and self.hi >= other.hi

    def __rge__(self, other):
        return self <= other

    def __lt__(self, other):
        """ Inequality `<` """
        other = self.make_interval(other)
        return self.lo < other.lo and self.hi < other.hi

    def __rlt__(self, other):
        return self > other

    def __gt__(self, other):
        """ Inequality `>` """
        other = self.make_interval(other)
        return self.lo > other.lo and self.hi > other.hi

    def __rgt__(self, other):
        return self < other


    # Operations between intervals viewed as sets
    #
    def _is_empty_intersection(self, other):
        """ Check if the intersection between intervals is empty """
        return self.hi < other.lo or other.hi < self.lo

    def intersection(self, other):
        """ Intersection """

        other = self.make_interval(other)

        if self._is_empty_intersection(other):
            print ("Intersection is empty:\n"
                   "Intervals {} and {} are disjoint".format(self,other) )
            return

        return Interval( max(self.lo,other.lo), min(self.hi,other.hi) )

    def hull(self, other):
        """Hull of two Intervals"""
        lower = min(self.lo,other.lo)
        upper = max(self.hi,other.hi)

        return Interval( lower, upper )

    def union(self, other):
        """Union of Intervals"""
        other = self.make_interval(other)

        if self._is_empty_intersection(other):
            print ("Union yields no connected interval:\n"
                   "Intervals {} and {} are disjoint".format(self,other) )
        else:
            return self.hull(other)

    # Scalar functions of Intervals:
    #
    def diam(self):
        return self.hi - self.lo

    def rad(self):
        return 0.5*self.diam()

    def mid(self):
        return 0.5*( self.lo + self.hi )

    def mag(self):
        """Maximum distance to the origin (magnitude)."""
        return max( abs(self.lo), abs(self.hi) )

    def mig(self):
        """Minimum distance to the origin (mignitude)."""
        if 0 in self:
            return 0
        else:
            return min( abs(self.lo), abs(self.hi) )

    def __abs__(self):
        """
        Defines the function abs() for an Interval, whose result is an interval.
        """
        return Interval( self.mig(), self.mag() )

    def dist(self, other):
        """
        Hausdorff distance between two intervals.
        """
        return max( abs(self.lo-other.lo), abs(self.hi-other.hi) )


    # convert a number, tuple etc. to an interval if necessary:
    # with the current implementation, this could be replaced by Interval(a),
    # which, however, unnecessarily constructs a new interval instead of just reusing the current one,
    # in the case that a is an interval:

    def make_interval(self, a):
        if isinstance(a, Interval):
            return a

        return Interval(a)



# def make_mpfr(a, rounding=None):
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
