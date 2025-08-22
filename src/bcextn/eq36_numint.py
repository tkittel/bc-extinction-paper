from .mpmath import mp, mpf

def numint_eq36( theta_degree, x, eps ):
    assert mpf(x) > 1e-20, "should clearly have used a taylor expansion!"
    eps = mpf(eps) * 0.1 #safety
    yp_0a = FiniteIntegralEq36( theta_degree, x, eps )
    while True:
        yp_ainf = TailIntegral( x, yp_0a.a )
        v = yp_0a.value + yp_ainf.value
        e = yp_0a.error + yp_ainf.error
        assert 0 < v < 1.0
        assert 0 < e < 1.0
        if e < eps * v:
            assert v-e > 0
            assert v+e < 1
            return v, e
        yp_0a.grow_a()

def integrand_fct_eq36( theta_degree, x ):
    from .eval_fofeta import F_of_Eta
    from .eval_phipi import PhiPi
    from .eval_phi0 import Phi0
    f = F_of_Eta()
    phi0 = Phi0()
    phipi = PhiPi()
    x = mpf(x)
    sinth = mp.sin(mpf(theta_degree) * mp.pi/mpf(180) )
    sqrtsinth = mp.sqrt( sinth )
    sinth_to_3div2 = sinth * sqrtsinth
    if not sinth:
        phi = phi0
    elif sinth == mpf(1):
        phi = phipi
    else:
        def phi( sr ):
            #BC 1974 eq. 34 ( sr = sigma*r )
            assert sr >= 0.0
            u = sr * sqrtsinth
            return phi0( sr ) + sinth_to_3div2 * ( phipi(u) - phi0(u) )
    def integrand( eta ):
        fval = f(eta)
        return fval * phi( x * fval )
    return integrand

class FiniteIntegralEq36:

    #integrates eq. 36 from 0 to a finite value, a.
    def __init__( self, theta_degree, x, eps ):
        self.__g = integrand_fct_eq36( theta_degree, x )
        self.__eps = eps
        assert 0 < eps < 1e-4
        #We always integrate over boundaries in steps of 1, since integrand
        #contains trigonometric terms with arguments eta or 2*eta.
        self.__step = max( mpf(1000), mp.ceil(mpf(x)) )
        self.__a = mpf(0)
        self.__val, self.__err = mpf(0), mpf(0)
        self.grow_a()

    def _get_raw_integrand( self ):
        #For unit testing!
        return self.__g

    #Results: integral over [0,a] gives an estimated value with a relative error:
    @property
    def a( self ):
        return mpf( self.__a )

    @property
    def value( self ):
        return mpf( self.__val )

    @property
    def relerror( self ):
        return mpf( self.__err / self.__val )

    @property
    def error( self ):
        return mpf( self.__err )

    #If a is deemed too low, it can be grown by calling the following function:
    def grow_a( self ):
        if self.__a > 1e7:
            raise RuntimeError('Growing "a" beyond 1e7!')
        eta_low = self.__a
        self.__a += self.__step
        v, e = self.__integrate( eta_low, self.__a )
        self.__val += v
        self.__err += e #conservative, assuming fully correlated errors

    def __integrate( self, eta_low, eta_high ):
        g = self.__g
        assert eta_high == int(eta_high)
        assert eta_low == int(eta_low)
        assert mpf(0) <= eta_low < eta_high
        n_intervals = eta_high - eta_low
        assert n_intervals == int(n_intervals)
        bounds = mp.linspace( eta_low, eta_high, int(n_intervals)+1 )
        maxdegree_min = 3
        maxdegree_max = 6
        ok = False
        for maxdegree in range(maxdegree_min,maxdegree_max+1):
            print("mp.quad( .., len(bounds)=%i, maxdegree=%i )"%( len(bounds),
                                                                  maxdegree) )
            res, error = mp.quad( g,
                                  bounds,
                                  method='gauss-legendre',
                                  maxdegree=maxdegree,
                                  error=True )
            assert res>0
            assert 0<error<1
            ok = error < 0.1 * self.__eps * res
            if ok:
                break
        if not ok:
            raise RuntimeError("WARNING: Problems integrating!")
        k = ( mpf(6)/ (mpf(4)*mp.pi) )
        return k * res, k * error

class TailIntegral:
    def __init__(self,x,a):
        #NB: Does not depend on theta!
        a = mpf(a)
        x = mpf(x)
        assert a >= x
        assert a >= 100
        a1 = mpf(1) / mpf(a)
        a2 = a1 * a1
        a3 = a2 * a1
        v = a1 - mpf('1/2') * x * a3
        e = mpf('1/2') * a2 + mpf(2) * a3
        k = ( mpf(6)/ (mpf(4)*mp.pi) )
        self.__val = k * v
        self.__err = k * e

    @property
    def a( self ):
        return mpf( self.__a )

    @property
    def value( self ):
        return mpf( self.__val )

    @property
    def relerror( self ):
        return mpf( self.__err / self.__val )

    @property
    def error( self ):
        return mpf( self.__err )
