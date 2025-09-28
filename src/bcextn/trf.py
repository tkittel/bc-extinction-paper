import numpy as np
import math

# Here we are implementing the transformation (x,y) -> (xprime,yprime) and its
# inverse. The point is that yprime(xprime) curves are more appropriate for
# fitting with Legendre polynomials than y(x) curves.

def _primeimpl_yshift( xp ):
    # We want to shift y values by some function f(xprime) so they have a
    # limiting value of 0.0 at both xprime=-1 (x=0) and xprime=+1 (x=inf). This
    # means we have to subtract 1 at xprime=-1, and 0 at
    # xprime=+1. Additionally, we do not want to change the slope at these two
    # limits, so we want the derivative of f to vanish at xprime=+-1. This gives
    # us four conditions, so by picking to describe f with a third order
    # polynomial for, this determines all coefficients.

    ##
    ##sage: var('C0','C1','C2','C3')
    ##sage: f=C0+C1*x+C2*x**2+C3*x**3
    ##sage: solve([f(x=-1)==-1,f(x=1)==0,f.diff(x,1)(x=-1)==0,f.diff(x,1)(x=1)==0],(C0,C1,C2,C3))
    ##[[C0 == (-1/2), C1 == (3/4), C2 == 0, C3 == (-1/4)]]

    #Using these coefficients, we get the following expression:
    return ( 0.25*xp*(3.0 - xp*xp) -0.5 )

def _primeimpl_yscale( xp ):
    # Having already shifted y-values according to _primeimpl_yshift, we know
    # that the y-curve must now pass through 0 at xprime = +-1 (and also with
    # unchanged slope).
    #
    # Now, we ultimately need to enforce the fitted curves to also pass through
    # the points (-1,0) and (+1,0). However, we also want to be able to perform
    # the fit with the standard (and rather awesome) Legendre polynomials, which
    # are not constrained to go through these points. So as a first idea, we
    # might want to fit with (1+xp)*(1-xp)*LEGENDRE_POLYNOMIALS. However, it is
    # better to simply divide out this factor of (1+xp)*(1-xp) first, because
    # then we can rely on completely standard Legendre polynomial
    # routines. Additionally, the factors of (1+-xp) could appear in any power,
    # and it turns out that we get vastly superior result if we use a power of 2
    # for the (1+xp) factor and a power of 1 for the factor of (1-xp). This was
    # simply seen by trying out various powers, but is most likely related to
    # the way the y(x) curves behavour near x=0.
    return  ( ( 1 + xp )**2 ) * ( 1 - xp )

#def y_unprime( xp, yp ):
#    y = yp.copy() if hasattr(yp,'copy') else yp
#    y *= _primeimpl_yscale(xp)
#    y -= _primeimpl_yshift(xp)
#    return y

def xy_unprime( xp, yp ):
    x = ( (xp+1)/(1-xp) )**2
    if yp is None:
        return xp
    y = yp.copy()
    y *= _primeimpl_yscale(xp)
    y -= _primeimpl_yshift(xp)
    return x, y

def x_to_xprime(x):
    xpow = np.sqrt(x)
    return ( xpow-1.0 ) / (xpow+1.0)

@np.vectorize
def xy_prime( x, y ):
    if math.isinf(x) and x>0:
        xp = 1.0
    else:
        xp = x_to_xprime(x)
    if y is None:
        return xp
    yp = y.copy() if hasattr(y,'copy') else y
    yp += _primeimpl_yshift(xp)
    yp /= _primeimpl_yscale(xp)
    return xp, yp

def transform_to_xprime(xvals):
    return xy_prime( xvals, None )

def transform_to_yprime(yvals,xvals):
    xp, yp = xy_prime(yvals, xvals)
    return yp
