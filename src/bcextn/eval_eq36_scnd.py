from .mpmath import mp, mpf
from .integrand import create_integrand_fct_eq36
from .eq36_lowx_taylor import (taylor_lowx_eq36_scnd_gauss,
                               taylor_lowx_eq36_scnd_lorentz,
                               taylor_lowx_eq36_scnd_fresnel)

eps_val = 1e-20

def eval_eq36_scnd_gauss( theta_degree, x ):
    t = taylor_lowx_eq36_scnd_gauss(theta_degree,x,eps=eps_val)
    if t is not None:
        return t
    return numint_eq36_scnd_gauss( theta_degree, x, eps=eps_val )

def eval_eq36_scnd_lorentz( theta_degree, x ):
    t = taylor_lowx_eq36_scnd_lorentz(theta_degree,x,eps=eps_val)
    if t is not None:
        return t
    return numint_eq36_scnd_lorentz( theta_degree, x, eps=eps_val )

def eval_eq36_scnd_fresnel( theta_degree, x ):
    t = taylor_lowx_eq36_scnd_fresnel(theta_degree,x,eps=eps_val)
    if t is not None:
        return t
    return numint_eq36_scnd_fresnel( theta_degree, x, eps=eps_val )

def numint_eq36_scnd_lorentz_old( theta_degree, x, eps ):

    from .eq36_tail_scnd_lorentz import TailIntegral
    from .eval_fofeta_scnd_lorentz import F_of_Eta

    eps = mpf(eps) * 0.1 #safety
    assert x > 0
    g = create_integrand_fct_eq36( f_of_eta_fct = F_of_Eta(),
                                   theta_degree = theta_degree,
                                   x = x,
                                   xscale = mpf('4/3') )

    #Very careful bounds (this is most likely overkill):
    def step_to_at_least( boundslist, aa, nn ):
        step = aa / nn
        while boundslist[-1]<aa:
            boundslist.append( boundslist[-1] + step )

    bounds = [0]
    n=1
    step_to_at_least( bounds, 0.001, 5*n )
    step_to_at_least( bounds, 0.01, 5*n )
    step_to_at_least( bounds, 0.1, 5*n )
    step_to_at_least( bounds, 1, 10*n )
    step_to_at_least( bounds, 10, 10*n )
    step_to_at_least( bounds, 2*max(x,10), 10*n )
    k_norm = mpf(2)/mp.pi
    def do_quad( thebounds, maxdegree = 8 ):
        #print("do_quad ( nbounds = %i, maxdegree=%i)"%(len(thebounds),maxdegree))
        assert thebounds[-1] < 1e6
        assert maxdegree < 15
        res, err = mp.quad( g,
                            thebounds,
                            method='gauss-legendre',
                            maxdegree=maxdegree,
                            error=True )
        res *= k_norm
        err *= k_norm
        if err < res * eps:
            return res, err
        return do_quad( thebounds, maxdegree + 1 )

    val_quad, err_quad = do_quad( bounds )
    while True:
        tail = TailIntegral( theta_degree, x, bounds[-1] )
        #Adding errors conservatively
        val, err = val_quad  + tail.value, err_quad + tail.error
        if err < val * eps:
            assert err < eps * abs(val)
            assert val-err > 0.0
            assert val+err < 1.0
            return val, err
        #Update quad result:
        new_bounds = [ bounds[-1] ]
        step_to_at_least( new_bounds, bounds[-1]*20, 10*n )
        bounds = new_bounds
        val2_quad, err2_quad = do_quad( bounds )
        val_quad += val2_quad
        err_quad += err2_quad
        del new_bounds, val2_quad, err2_quad

def numint_eq36_scnd_lorentz( theta_degree, x, eps ):

    from .eq36_tail_scnd_lorentz_NEW import tailint_auto
    from .eval_fofeta_scnd_lorentz import F_of_Eta

    eps = mpf(eps) * 0.1 #safety
    assert x > 0
    g = create_integrand_fct_eq36( f_of_eta_fct = F_of_Eta(),
                                   theta_degree = theta_degree,
                                   x = x,
                                   xscale = mpf('4/3') )

    #Very careful bounds (this is most likely overkill):
    def step_to_at_least( boundslist, aa, nn ):
        step = aa / nn
        while boundslist[-1]<aa:
            boundslist.append( boundslist[-1] + step )

    bounds = [0]
    n=1
    step_to_at_least( bounds, 0.001, 5*n )
    step_to_at_least( bounds, 0.01, 5*n )
    step_to_at_least( bounds, 0.1, 5*n )
    step_to_at_least( bounds, 1, 10*n )
    step_to_at_least( bounds, 10, 10*n )
    step_to_at_least( bounds, 2*max(x,10), 10*n )
    k_norm = mpf(2)/mp.pi
    def do_quad( thebounds, maxdegree = 8 ):
        #print("do_quad ( nbounds = %i, maxdegree=%i)"%(len(thebounds),maxdegree))
        assert thebounds[-1] < 1e6
        assert maxdegree < 15
        res, err = mp.quad( g,
                            thebounds,
                            method='gauss-legendre',
                            maxdegree=maxdegree,
                            error=True )
        res *= k_norm
        err *= k_norm
        if err < res * eps:
            return res, err
        return do_quad( thebounds, maxdegree + 1 )

    ti = tailint_auto( theta_degree=theta_degree, x=x )
    while ti['a']<bounds[-1]:
        bounds = bounds[:-1]
    bounds.append(ti['a'])

    val_quad, err_quad = do_quad( bounds )

    val, err = val_quad  + ti['tailint_val'], err_quad + ti['tailint_err']
    assert err < val * eps
    assert err < eps * abs(val)
    assert err < eps * (1.0-abs(val))
    assert err < eps * (1.0-val)
    assert val-err > 0.0
    assert val+err < 1.0
    return val, err

def numint_eq36_scnd_gauss( theta_degree, x, eps ):
    from .eval_fofeta_scnd_gauss import F_of_Eta
    eps = mpf(eps) * 0.1 #safety
    g = create_integrand_fct_eq36( f_of_eta_fct = F_of_Eta(),
                                   theta_degree = theta_degree,
                                   x = x )

    bounds = [ mpf(0), mpf(0.01), mpf(0.1), mpf(0.5), mpf(1),
               mpf(2), mpf(4), mpf(6), mpf(8), mpf(10), mpf(20),
               mpf(30), mpf(50), mpf(70), mpf(100), mpf(1000), mpf('inf') ]
    val, err = mp.quad( g,
                        bounds,
                        method='gauss-legendre',
                        maxdegree=6,
                        error=True )
    k_norm = mpf(6) / (mpf(4)*mp.pi)
    val *= k_norm
    err *= k_norm
    assert 0 < val < 1
    assert err < eps * abs(val)
    return val, err

def numint_eq36_scnd_fresnel( theta_degree, x, eps, print_ref = False ):
    if x == 0:
        #Avoid issues, but obviously prefer Taylor expansion in this case:
        return mpf(1), mpf(0)

    assert x<=1000,"does not seem numerically stable at higher x"

    #Very tricky, but we try with Richardsons extrapolation.
    from .eval_fofeta_scnd_fresnel import F_of_Eta
    eps = mpf(eps) * 0.01 #safety
    g = create_integrand_fct_eq36( f_of_eta_fct = F_of_Eta(),
                                   theta_degree = theta_degree,
                                   x = x )

    tol = eps * 1e-15

    #Zeroes are at multiples of ((4/3)*pi), so we define the n'th contribution
    #as the integral over [ n * ((4/3)*pi),  (n+1) * ((4/3)*pi) ]
    def contribn( n ):
        assert n>=0 and n==int(n)
        n=int(n)
        maxdegree = 8# 4 is fastest, since 3 triggers retry below
        bounds = [ mpf(n*4)*mp.pi/(mpf(3)),
                   mpf((n+1)*4)*mp.pi/(mpf(3)) ]
        if n>0:
            assert g(bounds[0])<1e-30, g(bounds[0])
        assert g(bounds[1])<1e-30
        while True:
            assert maxdegree < 15
            val, err = mp.quad( g,
                                bounds,
                                method='gauss-legendre',
                                maxdegree=maxdegree,
                                error=True )
            if err < tol * val * 1e-20: #individual terms better have a lot
                                        #smaller error than final result.
                return val
            maxdegree += 1
    res = mp.nsum( contribn,
                   [0,mp.inf],#n=0,1,2,...
                   tol = tol,
                   method='richardson',#easier to describe in paper than a multi-method
                   maxterms = 10000*mp.dps,#default is 10*dps
                   steps=[400,1]
                  )
    k_norm=mpf('3/2')/mp.pi
    if print_ref:
        #Brute-force testing:
        olddps = mp.dps
        for dps in 15, 20, 25, 40:
            mp.dps = dps
            res2 = k_norm*mp.quadosc( g, [0,mp.inf], omega=mpf('3/4') )
            print(f'QUADOSC (dps={dps}:',res2)
        mp.dps = olddps

    return ( res*k_norm, eps*k_norm )

def test_fresnel():

    #Reference values found by setting print_ref=True:
    ref_releps = 1e-18
    for (x,theta,ref) in [ (1000.0, 55, 0.01633068365545878076598065008548644904862),
                           (1000.0, 0, 0.01363204187306092208666920248561948425878),
                           (1000.0, 90, 0.01709396590396128947302559276288498214923),
                           #(1500.0, 55, 0.0133311338221847897148609546026948984359),#last few digits uncertain
                           (3.5,0,0.246707496779635973634137),
                           ( 5.0, 90, 0.2587106717313447081308667 ),
                           (0,45,1.0),
                           (100.0, 70, 0.05354615515435589189056355 ),
                           (0.0001, 45, 0.9999000106488811257694316),
                           (1.0, 30, 0.5097359403462856004969412  ),
                       ]:
        print_ref = (ref is None)
        print( f'( x,theta ) = ( {x:g}, {theta:g} )' )
        val, err = numint_eq36_scnd_fresnel(theta_degree=theta,x=x,eps=1e-10,
                                            print_ref = print_ref)
        print( float(val), float(err) )
        assert 0 < val - err
        assert val + err <= 1
        toterr = abs(ref_releps*ref) + abs(err)
        print("   val: %g +- %g"%(float(val), float(err)))
        print("   ref: %g +- %g"%(float(ref), float(ref_releps*ref)))
        print("   error from ref: {float(ref_releps*ref):g}")
        print("   error from nsum: {float(err):g}")
        toterr = abs(ref_releps*ref) + abs(err)
        print("   toterr: {float(toterr):g}")
        print("   ys_fresnel ndiff: ",abs(val-ref) / (toterr or 1e-200) )
        assert abs(val-ref) < toterr

def test():
    test_fresnel()

    from .eq36_lowx_taylor import (taylor_lowx_eq36_scnd_gauss,
                                   taylor_lowx_eq36_scnd_lorentz,
                                   taylor_lowx_eq36_scnd_fresnel)

    for name, numintfct, taylorfct in [ ( 'lorentz',
                                          numint_eq36_scnd_lorentz,
                                          taylor_lowx_eq36_scnd_lorentz ),
                                        ( 'gauss',
                                          numint_eq36_scnd_gauss,
                                          taylor_lowx_eq36_scnd_gauss ),
                                        ( 'fresnel',
                                          numint_eq36_scnd_fresnel,
                                          taylor_lowx_eq36_scnd_fresnel ),

                                       ]:
        for (x,sintheta) in [ (10.0, 0.5 ),
                              (0.0001, 0.0),
                              (0.1, 0.05),
                              (1.0, 0.05),
                              (5.0, 0.05),
                              (30.0, 0.05),
                              (0.1, 0.5),
                              (1.0, 0.5),
                              (5.0, 0.5),
                              (30.0, 0.5),
                              (0.1, 0.9),
                              (1.0, 0.9),
                              (5.0, 0.9),
                              (30.0, 0.9),
                              (0.00001, 1.0),
                              (1.0, 1.0),
                              (30.0, 1.0),

                              (0.4, 0.0),
                              (0.4, 0.01),
                              (0.4, 0.2),
                              (0.4, 0.4),
                              (0.4, 0.6),
                              (0.4, 0.8),
                              (0.4, 0.9999),
                              (0.4, 1.0),

                              (0.2, 0.0),
                              (0.2, 0.01),
                              (0.2, 0.2),
                              (0.2, 0.4),
                              (0.2, 0.6),
                              (0.2, 0.8),
                              (0.2, 0.9999),
                              (0.2, 1.0),
                             ]:
            theta_degree = mp.asin(sintheta)*180/mp.pi
            res_taylor = taylorfct(theta_degree, x = x, eps=1e-50)
            res,err = numintfct( theta_degree = theta_degree,
                                 x = x,
                                 eps = 1e-50 )
            print(f"y_scnd_{name}(x={x:4g}, sintheta={sintheta:4g}) = {float(res):g} +- {float(err):3g}")
            if res_taylor is not None:
                v,e = res_taylor
                reldiff = 2*abs(v-res)/abs(v+res)
                print(f"y_scnd_{name}(x={x:4g}, sintheta={sintheta:4g}) = {float(v):g} +- {float(e):3g} (TAYLOR reldiff: %.3g)"%float(reldiff))
                assert reldiff < e + err
