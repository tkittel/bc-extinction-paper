from .mpmath import mp,mpf

def _calc_taylor_coeff_phi0(i):
    assert isinstance(i,int)
    assert 0 <= i <= 500
    sgn = -1 if i%2 else 1
    i = mpf(i)
    return mpf(sgn*3)*(i+2)*(mpf(4)**i)/mp.factorial(i+3)

def _calc_taylor_coeff_phipi(i):
    assert isinstance(i,int)
    assert 0 <= i <= 500
    sgn = -1 if i%2 else 1
    i = mpf(i)
    return mpf(sgn*3)*(mpf(2)**i)/(i+3)

def _calc_taylor_coeff_phi(i, sqrt_sinth):
    if sqrt_sinth <= 0.0 or i<2:#FIXME: add "or i<2:"?
        return _calc_taylor_coeff_phi0(i)
    if sqrt_sinth >= 1.0:
        return _calc_taylor_coeff_phipi(i)
    cn_0  = _calc_taylor_coeff_phi0(i)
    cn_pi = _calc_taylor_coeff_phipi(i)
    return cn_0 + (mpf(sqrt_sinth)**(mpf(i)+3))*(cn_pi-cn_0)

class TailIntegral:

    @staticmethod
    def estimate_lowest_a( x ):
        return 2.0+4.0*float(x)

    def __init__(self,theta_degree, x, a, order):
        if not a >= TailIntegral.estimate_lowest_a(x):
            raise RuntimeError('%g not at least %g (x=%g)'%(a,TailIntegral.estimate_lowest_a(x),x))
        assert isinstance(order,int)
        assert 2 <= order <= 1000
        sinth = mp.sin(mpf(theta_degree) * mp.pi/mpf(180) )
        sqrt_sinth = mp.sqrt(sinth)
        _cache_cn_phi = {}
        def cphi_n(n_raw):
            n = int(n_raw)
            assert n_raw==n
            if n not in _cache_cn_phi:
                _cache_cn_phi[n] = _calc_taylor_coeff_phi(n,sqrt_sinth)
            return _cache_cn_phi[n]
        #cn_phi = [ _calc_taylor_coeff_phi(i,sqrt_sinth)
        #           for i in range(order+3) ]
        b = mpf('4/3')*mpf(x)
        #assert a>b
        inva = mpf(1)/mpf(a)

        def d_alt( k, n):
            return (-1)**k*mp.binomial(n+k-1,k)
        def d(k, n):
            v=(-1)**k*mp.factorial(n+k-1)/( mp.factorial(k) * mp.factorial(n-1) )
            assert d_alt(k,n)==v
            return v
        def D( k,n ):
            return d(k,n+1)/(2*(k+n)+1)
        def em( m_raw ):
            m = int(m_raw)
            assert m_raw == m
            res=mpf(0)
            for j in range(m+1):
                res += cphi_n(j)*D(m-j,j)*(b**j)
            return res
        terms = [em(m)*(inva**(2*m+1)) for m in range(order+1)]
        assert order >= 4
        norm = mpf(2) / mp.pi #= factor of 6*c_M/4pi from eq. 2 in paper
        self.__val = norm * sum(terms)
        self.__err = norm * sum(abs(e) for e in terms[-3:])#FIXME: Too conservative?

    @property
    def value( self ):
        return mpf( self.__val )

    #@property
    #def relerror( self ):
    #    return mpf( self.__err / self.__val )

    @property
    def error( self ):
        return mpf( self.__err )

def tailint_auto( theta_degree, x ):
    a = TailIntegral.estimate_lowest_a(x)
    #print("INITIAL a=%g"%a)
    #a = 40
    #NOTE: order=40 implies terms up to (1/a)**(2*40+1)=(1/a)**(2*40+1)
    order = 30
    eps=1e-20#at least 1e-8, but with loooots of safety
    while True:
        assert 2<a<5000
        tf=TailIntegral( theta_degree=theta_degree,
                         x=x, a=a, order=order )
        if ( not (0.0<tf.value<1.0)
             or not (0.0<tf.error<1.0)
             or not (0.0<tf.value+tf.error<1.0)
             or not (0.0<tf.value-tf.error<1.0)
             or tf.error > eps*min(tf.value,1.0-tf.value) ):
            a *= 2
            print("WARNING: Increasing a to a=%g"%a)
            raise RuntimeError("Had to increase a!!!")
            continue
        return dict( theta_degree = theta_degree,
                     x=x,
                     a=a,
                     tailint_val = tf.value,
                     tailint_err = tf.error )
