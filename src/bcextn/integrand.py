def create_integrand_fct_eq36( f_of_eta_fct, theta_degree, x, xscale = 1 ):
    #Integrand: f(eta)*phi(theta,xscale*x*f(eta))
    from .mpmath import mp, mpf
    from .eval_phipi import PhiPi
    from .eval_phi0 import Phi0
    f = f_of_eta_fct
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
    xscale = mpf(xscale)
    def integrand( eta ):
        fval = f(eta)
        return fval * phi( xscale * x * fval )
    return integrand
