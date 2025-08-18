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
    k = 1 + B*x
    k2 = 1.0+2.0*x+safediv(A*x*x+C*x,k)
    return safediv( 1.0, safesqrt( k2 ) )

class ProposedCurve:

    def __init__(self):
        from .data import load_json_data
        data = load_json_data('fitted_curves_ABC.json')
        def sinth(theta):
            return math.sin( float(theta)*kDeg2Rad )
        p0, p1, p2, p3, p4, p5, p6 = data['A']
        def f( theta, p0, p1, p2,p3,p4,p5,p6 ):
            u = sinth(theta)
            return p0 + p1*u + p2*u**2 + p3*u**3 +p4*u**4 +p5*u**5 + p6*u**6
        self._A = np.vectorize(f)
        p0, p1, p2, p3, p4, p5, p6 = data['B']
        def f( theta, p0, p1, p2,p3,p4,p5,p6 ):
            u = sinth(theta)
            return p0 + p1*u + p2*u**2 + p3*u**3 +p4*u**4 +p5*u**5 + p6*u**6
        self._B = np.vectorize(f)
        p0, p1, p2, p3, p4, p5, p6 = data['C']
        def f( theta, p0, p1, p2,p3,p4,p5,p6 ):
            u = sinth(theta)
            return p0 + p1*u + p2*u**2 + p3*u**3 +p4*u**4 +p5*u**5 + p6*u**6
        self._C = np.vectorize(f)

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
        p0, p1 = data['origA']
        def f( theta, p0, p1 ):
            u = cos2th(theta)
            return p0 + p1*u
        self._A = np.vectorize(f)
        p0, p1 = data['origB']
        def f( theta, p0, p1 ):
            u = cos2th(theta)
            return p0 + p1*(0.5-u)**2
        self._B = np.vectorize(f)

    def A( self, theta ):
        return self._A(theta)

    def B( self, theta ):
        return self._B(theta)

    def __call__( self, x, theta ):
        A, B = self._A(theta), self._B(theta)
        return yp_bc1974_curve( x, A, B )
