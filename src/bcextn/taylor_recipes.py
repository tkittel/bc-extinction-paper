#Taylor expansions, for user recipes (FIXME: Unclear if we will use this ultimately):

import math

def taylor_lowx_scndgauss( x, sintheta ):
    u = 2*x
    s_sqrt = math.sqrt(sintheta)
    s2 = sintheta*sintheta
    s3, s4 = sintheta*s2, s2*s2
    s5 = s2*s3
    sqrt2, sqrt3 = 1.4142135623730951, 1.7320508075688772
    sqrt5, sqrt7 = 2.23606797749979, 2.6457513110645907
    c0 = 1.0
    c1 = sqrt2*(-3./8.)
    c2 = sqrt3*(2./15.)+sqrt3*(1./15.)*s2*s_sqrt
    c3 = (-1./12.)+(-1./6.)*s3
    c4 = sqrt5*(2./175.)+sqrt5*(13./175.)*s3*s_sqrt
    c5 = sqrt3*sqrt2*(-1./360.)+sqrt3*sqrt2*(-43./720.)*s4
    c6 = sqrt7*(4./6615.)+sqrt7*(311./6615.)*s4*s_sqrt
    c7 = sqrt2*(-1./4200.)+sqrt2*(-157./2100.)*s5
    c8 = (2./31185.)+(2833./31185.)*s5*s_sqrt
    return c0+u*(c1+u*(c2+u*(c3+u*(c4+u*(c5+u*(c6+u*(c7+u*c8)))))))
