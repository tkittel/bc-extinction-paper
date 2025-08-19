import numpy as np
import math
kDeg2Rad = np.pi/180

def safediv( a, b ):
    return a/b if abs(b)>1e-99 else 1e99

def safesqrt( x ):
    return math.sqrt(max( 0.0, x ))

@np.vectorize
def yp_bc1974_curve( x, A, B ):
    k = 1 + B*x
    k2 = 1.0+2.0*x+safediv(A*x*x,k)
    return safediv( 1.0, safesqrt( k2 ) )

@np.vectorize
def yp_proposed_curve( x, A, B, C ):
    k = 1.0 + B*x + C*x/(1.0+x)
    k2 = 1.0+2.0*x+safediv( A*x*x-0.1*x, k )
    return safediv( 1.0, safesqrt( k2 ) )

class ProposedCurve:

    def __init__(self):
        from .data import load_json_data
        data = load_json_data('fitted_curves_ABC.json')
        def sinth(theta):
            return math.sin( float(theta)*kDeg2Rad )
        pA0, pA1, pA2, pA3, pA4, pA5, pA6 = data['A']
        pB0, pB1, pB2, pB3, pB4, pB5, pB6 = data['B']
        pC0, pC1, pC2, pC3, pC4, pC5, pC6 = data['C']
        def fA( theta ):
            u = sinth(theta)
            return pA0 + pA1*u + pA2*u**2 + pA3*u**3 +pA4*u**4 +pA5*u**5 + pA6*u**6
        def fB( theta ):
            u = sinth(theta)
            return pB0 + pB1*u + pB2*u**2 + pB3*u**3 +pB4*u**4 +pB5*u**5 + pB6*u**6
        def fC( theta ):
            u = sinth(theta)
            return pC0 + pC1*u + pC2*u**2 + pC3*u**3 +pC4*u**4 +pC5*u**5 + pC6*u**6
        self._A = np.vectorize(fA)
        self._B = np.vectorize(fB)
        self._C = np.vectorize(fC)

    def A( self, theta ):
        return self._A(theta)

    def B( self, theta ):
        return self._B(theta)

    def C( self, theta ):
        return self._C(theta)

    def __call__( self, x, theta ):
        A, B, C = self._A(theta), self._B(theta), self._C(theta)
        return yp_proposed_curve( x, A, B, C )

class UpdatedClassicCurve:

    def __init__(self):
        from .data import load_json_data
        data = load_json_data('fitted_curves_ABC.json')
        def cos2th(theta):
            return math.cos( 2*float(theta)*kDeg2Rad )
        pA0, pA1 = data['origA']
        pB0, pB1 = data['origB']
        def fA( theta ):
            u = cos2th(theta)
            return pA0 + pA1*u
        def fB( theta ):
            u = cos2th(theta)
            return pB0 + pB1*u
        self._A = np.vectorize(fA)
        self._B = np.vectorize(fB)

    def A( self, theta ):
        return self._A(theta)

    def B( self, theta ):
        return self._B(theta)

    def __call__( self, x, theta ):
        A, B = self._A(theta), self._B(theta)
        return yp_bc1974_curve( x, A, B )

class ClassicCurve:
    def __init__(self):
        def cos2th(theta):
            return math.cos( 2*float(theta)*kDeg2Rad )
        def fA( theta ):
            u = cos2th(theta)
            return 0.2 + 0.45*u
        def fB( theta ):
            u = cos2th(theta)
            return 0.22 -0.12*(0.5-u)**2
        self._A = np.vectorize(fA)
        self._B = np.vectorize(fB)

    def A( self, theta ):
        return self._A(theta)

    def B( self, theta ):
        return self._B(theta)

    def __call__( self, x, theta ):
        A, B = self._A(theta), self._B(theta)
        return yp_bc1974_curve( x, A, B )
