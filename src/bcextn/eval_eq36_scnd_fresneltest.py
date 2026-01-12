from .mpmath import mp, mpf
from .integrand import create_integrand_fct_eq36

def evalref( theta_degree, x,
             nmin, nmax,
             eps = 1e-20,
             print_ref = False,
             maxdegree0 = 5, #4 is typically faster than 3, since 3 often
                             #triggers retry below
            ):
    assert x > 0.0
    assert 0.0 <= theta_degree <= 90.0
    assert x<=1000
    assert nmin>=0 and nmin==int(nmin)
    assert nmax>=0 and nmax==int(nmax)
    assert nmax >= nmin
    nmin, nmax = int(nmin), int(nmax)

    from .eval_fofeta_scnd_fresnel import F_of_Eta
    #If all individual terms have the desired relative error, we should
    #be good to go also for the overall relative error. However, we do
    #add a safety factor of 1e-6:
    eps = mpf(eps) * 1e-6

    g = create_integrand_fct_eq36( f_of_eta_fct = F_of_Eta(),
                                   theta_degree = theta_degree,
                                   x = x )

    #Zeroes are at multiples of ((4/3)*pi), so we define the n'th contribution
    #as the integral over [ n * ((4/3)*pi),  (n+1) * ((4/3)*pi) ]
    def contribn( n ):
        assert n>=0 and n==int(n)
        n=int(n)
        bounds = [ mpf(n*4)*mp.pi/(mpf(3)),
                   mpf((n+1)*4)*mp.pi/(mpf(3)) ]
        if n>0:
            assert g(bounds[0])<1e-30
        assert g(bounds[1])<1e-30
        maxdegree = maxdegree0
        while True:
            assert maxdegree < 15
            val, err = mp.quad( g,
                                bounds,
                                method='gauss-legendre',
                                maxdegree=maxdegree,
                                error=True )
            if err < eps * val:
                return val, err
            maxdegree += 1
            print("WARNING: Increase maxdegree to %i"%maxdegree)

    k_norm = mpf('3/2') / mp.pi

    res = []
    for n in range( nmin, nmax+1 ):
        val, err = contribn(n)
        val *= k_norm
        err *= k_norm
        res.append( (n, val, err ) )
    return res
