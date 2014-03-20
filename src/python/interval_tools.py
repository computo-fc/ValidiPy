# Miscellaneous functions for working with intervals

import numpy as np
from interval import *

from matplotlib import pyplot as plt


def random_interval( infimum=-10.0, supremum=10.0 ):
    num1a = np.random.uniform( infimum, supremum )
    num2a = np.random.uniform( infimum, supremum )

    return Interval( num1a, num2a )

def split_interval( x, num_divisions=1 ):
    """
    Divide an Interval into n=num_divisions Intervals of equal lengths
    """
    num_divisions = int(num_divisions)
    if num_divisions < 1:
        num_divisions = 1

    edge_points = np.linspace(x.lo, x.hi, num_divisions+1)
    intervals = [Interval(a, b) for (a,b) in
                         zip(edge_points[0:num_divisions+1], edge_points[1:num_divisions+2]) ]

    return intervals

def range_interval_f( fun, subdivided_interval ):
    """
    Evaluate the function f(x) extended over Intervals, in a list of subIntervals,
    and return the hull of them all, i.e. a bound on the range of the function
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
    Plots the interval extension of a function `fun` over the interval `x`,
    which is divided into num=1, 2, 4, ..., 2**pow2 uniform subintervals.
    The original function is plotted at num_points points
    """
    num_intervals = [ 2**p for p in range(pow2+1) ]

    # plt.figure()
    # plt.subplot(1, 1, 1)

    for num in num_intervals:
        fact_alfa = num * 1.0/num_intervals[-1]   # for plotting

        subdivided_intervals = split_interval( x, num )
        total_range = range_interval_f( fun, subdivided_intervals )
        print "Total range (N={}) = {}".format(num, total_range)

        # Draw each subinterval and its range:
        for x1 in subdivided_intervals:
            low = float(x1.lo)
            high = float(x1.hi)
            Ffun = fun(x1)

            # draw a block of the correct size:
            xa1 = np.array([low, low, high, high])
            ya1 = np.array([float(Ffun.lo), float(Ffun.hi), float(Ffun.hi), float(Ffun.lo) ])
            plt.fill( xa1, ya1, 'b', alpha=fact_alfa )

    # plot the graph of the function:
    low = float(x.lo)
    high = float(x.hi)
    xx = np.linspace(low,high,num_points)
    yy = [fun(x) for x in xx]
    plt.plot( xx, yy, 'red')
