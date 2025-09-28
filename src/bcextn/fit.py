import numpy as np
from inspect import signature
import scipy.optimize
import matplotlib.pyplot as plt
import math

class FitFunc:
    def __init__( self, rawfct, p0 = None ):
        s = signature(rawfct)
        self.__name = rawfct.__name__
        self.__fct = np.vectorize(rawfct)
        self.__nparams = len(s.parameters)-1
        if p0 is None:
            p0 = [ 1.0, ]*self.__nparams
        else:
            assert len(p0) == self.__nparams
        self.__p0 = p0
    @property
    def nparams( self ):
        return self.__nparams
    @property
    def p0( self ):
        return self.__p0
    @property
    def name( self ):
        return self.__name
    def __call__(self, x, *p):
        return self.__fct(x,*p)
    def fit( self, xvals, yvals, p0 = None ):
        p0 = self.__p0 if p0 is None else p0
        assert len(p0)==self.__nparams
        if self.nparams:
            pres,_ = scipy.optimize.curve_fit( self.__fct, xvals, yvals, p0 = p0 )
        else:
            pres = []
        pres = [float(e) for e in pres]

        class FitResult:
            def __init__(self,name,fct,p):
                self.__n = name
                self.__f = fct
                self.__p = p
            @property
            def name( self ):
                return self.__n
            @property
            def parameters( self ):
                return self.__p
            def __call__(self, x):
                return self.__f(x,*(self.__p))
        return FitResult(self.__name,self.__fct,pres)

def fitfunction( func ):
    return FitFunc(func)

def _print_taylor_code( poly_coefficients, varname='xp', mode='py',
                       resvarname = 'yp'):
    np = len(poly_coefficients)
    assert np>1
    assert mode in ['cpp','py']
    def fmtcoeff(val):
        return '%.14g'%val
    s = 'const double ' if mode=='cpp' else ''
    s += '%s = '%resvarname
    for i,c in enumerate(poly_coefficients):
        if i+1==np:
            if c < 0.0:
                s += '-%s*%s'%(varname,fmtcoeff(abs(c)))
            else:
                s += '+%s*%s'%(varname,fmtcoeff(c))
        elif i:
            s += '+%s*(%s'%(varname,fmtcoeff(c))
        else:
            s += fmtcoeff(c)
    s+=')'*(np-2)
    if mode=='cpp':
        s+=';'
    print(s)
    return s


########################################################
###### Transform + Legendre fit ########################
########################################################

def legfit_weight_fct(y):
    return ( 1.0/np.minimum(abs(y),abs(1.0-y)) )

def main_legfit( args ):
    from .load import load_xscan090
    do_lux = False
    while 'lux' in args:
        args.remove('lux')
        do_lux = True
    include90 = False
    while 'include90' in args:
        args.remove('include90')
        include90 = True
    from . import new_recipes
    taylor_cutoff = new_recipes.recipe_taylor_cutoff(do_lux)

    #We only need to fit down to the Taylor cutoff (and go a bit below for less
    #abrupt change at the threshold):
    minx = 0.9*taylor_cutoff

    assert len(args)==0
    modes=['primary',
           'scndfresnel',
           'scndlorentz',
           'scndgauss',
           #'curve::sabineprimary',
           #'curve::sabinescndrec',
           #'curve::sabinescndtri',
           ]
    thvals=[0]
    if include90:
        #Just for testing
        thvals += [90]
    datasets = []
    for m in modes:
        data = load_xscan090(m)
        xvals = data['xvals'].copy()
        mask = xvals >= minx
        xvals = xvals[mask]
        for th, y in data['theta_2_ypvals'].items():
            if thvals and not any( abs(float(th)-e)<1e-6 for e in thvals ):
                continue
            assert int(float(th))==float(th)
            datasets.append( ( xvals,
                               y.copy()[mask],
                               f'{m} ({float(th):g} deg)',
                               f'polycoeff_xprime_to_y-{m}_at_th{int(float(th))}' ) )
    do_legfit( datasets,
               do_lux=do_lux )

def try_legfit( x, y, order, nchop  ):
    from .trf import xy_prime
    from .poly import raw_yprime_poly_to_y_poly

    xp, yp = xy_prime( x, y )
    fit_weights = legfit_weight_fct(y)
    legendre_coefs = np.polynomial.legendre.legfit( xp, yp,
                                                    order,
                                                    w = fit_weights )
    poly_coeffs_yp = np.polynomial.legendre.leg2poly(legendre_coefs)
    poly_coeffs = raw_yprime_poly_to_y_poly( poly_coeffs_yp )

    #Chop values (and remove highest orders if necessary):
    def do_chop( iorder, val ):
        return float((f'%.{nchop}g')%val)
    poly_coeffs = [do_chop(i,val) for i,val in enumerate(poly_coeffs)]

    while not poly_coeffs[-1]:
        poly_coeffs = poly_coeffs[:-1]

    ypfit = np.polynomial.polynomial.polyval(xp, poly_coeffs_yp)
    yfit = np.polynomial.polynomial.polyval(xp, poly_coeffs)
    diff_yp = abs(ypfit-yp)
    diff_y = ( abs(yfit-y)/np.minimum(abs(y),abs(1-y)) )
    maxdiff = diff_y.max()
    return dict( maxdiff = maxdiff,
                 diff_y = diff_y,
                 diff_yp = diff_yp,
                 poly_coeffs = poly_coeffs,
                 poly_coeffs_yp = poly_coeffs_yp,
                 xp = xp,
                 yp = yp,
                 ypfit = ypfit,
                 yfit = yfit )

def cost_legfit_recipe( order, nchop ):
    #Estimate "cost" of recipe, when written as c0+xp*(c1+xp*(...)). Meaning
    #that we want to find the balance between having as few polynomial orders as
    #possible, but also as few number of digits needed in each coefficient. To
    #not get near ulp, we will allow at most 14 digits.
    assert 2 <= nchop <= 14
    assert 2 <= order <= 100
    string_length_approx = (10+nchop)*order
    speed_cost = order
    cost = float(string_length_approx)*float(speed_cost)**3
    return cost

def do_legfit( datasets, do_lux ):
    from .trf import xy_prime, xy_unprime
    from .new_recipes import recipe_target_prec_each_main_term
    target_prec = recipe_target_prec_each_main_term(do_lux)
    target_prec *= 0.9 #safety

    output_filename='legendre_coefficients'
    if do_lux:
        output_filename += '_lux'
    output_filename += '.json'

    results = {}
    fitresults = []

    order_nchop_tries = []

    #Try various combinations of order and nchop, and try lower "cost"
    #combinations first.
    for order in range( 4, 60 ):
        for nchop in list(range(4,14+1)):
            order_nchop_tries.append( ( cost_legfit_recipe(order,nchop),
                                        order, nchop ) )
    order_nchop_tries.sort()

    for x, y, lbl, resvarname in datasets:
        ok = False
        for _, order, nchop in order_nchop_tries:
            print(f"    Trying legendre {order}-order fit for {resvarname} (nchop={nchop})")
            lf = try_legfit( x, y, order, nchop)
            if lf['maxdiff'] < target_prec:
                ok = True
                break
        if not ok:
            raise SystemExit('Failed to find proper fit')
        fitresults.append( ( lf, order, nchop) )
        results[resvarname] = [ float(e) for e in lf['poly_coeffs'] ]
        _print_taylor_code(lf['poly_coeffs'],resvarname=resvarname)
        print(f"DONE at order={order}")
    print()
    print()
    print()

    highest_maxdiff = 0.0
    for (x, y, lbl, resvarname),(lf,order,nchop) in zip(datasets,fitresults):
        highest_maxdiff = max(highest_maxdiff,lf['maxdiff'])
        print("Worst precision of %s : %g   (order=%i, nchop=%i)"%(resvarname,lf['maxdiff'],order,nchop))
    print("Worst precision overall: %g"%highest_maxdiff)

    print()
    for (x, y, lbl, resvarname),(lf,order,nchop) in zip(datasets,fitresults):
        _print_taylor_code(lf['poly_coeffs'],resvarname=resvarname)

    from .json import save_json
    save_json(
        output_filename,
        results,
        force = True
    )

    fig, axs = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    ax, axdiff = axs

    data_colors = []

    for (x, y, lbl, resvarname), (lf,order,nchop) in zip(datasets,fitresults):
        xp, yp = xy_prime( x, y )
        col = ax.plot( xp, yp, label=lbl )[0].get_color()
        ax.plot(xp,lf['ypfit'],ls='-.',lw=5,color=col)
        axdiff.plot(xp,lf['diff_yp'],color=col)
        data_colors.append( col )

    ax.grid()
    axdiff.grid()
    axdiff.semilogy()
    axdiff.set_ylabel('yprimefit-yprime')
    ax.legend()
    ax.set_xlabel('xprime')
    ax.set_ylabel('yprime')
    plt.show()

    fig, axs = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    ax, axdiff = axs
    for (x, y, lbl,resvarname),(lf,_,_),col in zip(datasets,fitresults,data_colors):
        xp,yp = xy_prime( x, y )
        _, yfit = xy_unprime( xp, lf['ypfit'] )
        ax.plot(x,y,label=lbl,color=col)
        ax.plot(x,yfit,ls='-.',lw=5,color=col)
        axdiff.plot(x,lf['diff_y'],color=col)

    axdiff.plot( [x[0],x[-1]],[target_prec,]*2, color='black',lw=4, ls=':' )
    ax.grid()
    axdiff.grid()
    ax.semilogx()
    axdiff.semilogy()
    ax.legend()
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show()

def main_deltalegfit( args ):
    #legendre fits directly to y90(x)-y0(x)
    from .load import load_xscan090
    do_lux = False
    while 'lux' in args:
        args.remove('lux')
        do_lux = True
    from . import new_recipes
    #We only need to fit down to the Taylor cutoff (and go a bit below for less
    #abrupt change at the threshold):
    fit_xrange =   [ new_recipes.recipe_taylor_cutoff(do_lux)*0.9,
                     1e3 ]

    assert len(args)==0
    modes=['primary',
           'scndfresnel',
           'scndlorentz',
           'scndgauss',
           #'curve::sabineprimary',
           #'curve::sabinescndrec',
           #'curve::sabinescndtri',
           ]
    #thvals=[0,90]#Due to final interpolation, only these matters
    datasets = []
    for m in modes:
        data = load_xscan090(m)
        xvals = data['xvals'].copy()
        mask = np.logical_and( xvals >= fit_xrange[0]*0.999999, xvals <= fit_xrange[1]*1.000001 )
        xvals = xvals[mask]
        assert set(data['theta_2_ypvals'].keys()) == set(['0.0','90.0'])
        y0 = data['theta_2_ypvals']['0.0'][mask]
        y90 = data['theta_2_ypvals']['90.0'][mask]
        ydelta = (y90-y0)
        datasets.append( ( xvals, ydelta, y0,y90,
                           f'{m} (ydelta)',
                           f'ydelta_{m}' ) )
    do_deltalegfit( datasets, do_lux=do_lux )

def deltalegfit_weight_fct(ydelta,y0,y90):
    if False:
        #m = np.minimum(abs(y90),abs(1.0-y90))
        #m = np.minimum(abs(y90),abs(1.0-y90))
        m = abs(ydelta)
    else:
        m0 = np.minimum(abs(y0),abs(1.0-y0))
        m90 = np.minimum(abs(y90),abs(1.0-y90))
        #m = m0
        m = np.minimum(m0,m90)
    return ( 1.0/m )**1.0 #Fixme: revisit this, make it easy to describe in paper

def try_deltalegfit( x, ydelta, y0, y90, order, nchop  ):
    from .trf import x_to_xprime

    xp = x_to_xprime(x)
    fit_weights = deltalegfit_weight_fct(ydelta,y0,y90)
    legendre_coefs = np.polynomial.legendre.legfit( xp, ydelta,
                                                    order,
                                                    w = fit_weights )
    poly_coeffs = np.polynomial.legendre.leg2poly(legendre_coefs)

    #Chop values (and remove highest orders if necessary):
    def do_chop( iorder, val ):
        return float((f'%.{nchop}g')%val)
    poly_coeffs = [do_chop(i,val) for i,val in enumerate(poly_coeffs)]

    while not poly_coeffs[-1]:
        poly_coeffs = poly_coeffs[:-1]

    ydeltafit = np.polynomial.polynomial.polyval(xp, poly_coeffs)
    #_, yfit = xy_unprime( xp, ypfit )
    diff_ydelta = abs(ydeltafit-ydelta)
    m0 = np.minimum(abs(y0),abs(1.0-y0))
    m90 = np.minimum(abs(y90),abs(1.0-y90))
    #diff_y = ( abs(diff_ydelta)/np.minimum(m0,m90) )
    reldiff_ydelta = ( abs(diff_ydelta)/ydelta )
    reldiff_ydelta_to_y0 = ( abs(diff_ydelta)/y0 )
    reldiff_ydelta_to_smallest = ( abs(diff_ydelta)/np.minimum(m0,m90) )

    return dict( maxdiff = reldiff_ydelta_to_smallest.max(),
                 reldiff_ydelta = reldiff_ydelta,
                 reldiff_ydelta_to_smallest = reldiff_ydelta_to_smallest,
                 reldiff_ydelta_to_y0 = reldiff_ydelta_to_y0,
                 diff_ydelta = diff_ydelta,
                 poly_coeffs = poly_coeffs,
                 xp = xp,
                 ydeltafit = ydeltafit,
                 ydelta = ydelta )


def do_deltalegfit( datasets, do_lux ):
    from .trf import x_to_xprime

    from .new_recipes import recipe_target_prec_each_main_term
    target_prec = recipe_target_prec_each_main_term(do_lux)
    target_prec *= 0.9 #safety

    output_filename='legendre_coefficients_ydelta'
    if do_lux:
        output_filename += '_lux'
    output_filename += '.json'

    results = {}
    fitresults = []

    order_nchop_tries = []
    if do_lux:
        #Just use the same always, we anyway need a mega table (and tries shows
        #that all fits needs nchop in the range 12-14 when do_lux).
        nchop_vals = [8,9,10,11,12,13,14]#list(range(14,16+1))
        order_ini = 2
    else:
        nchop_vals = list(range(3,14+1))
        order_ini = 2

    for order in range( order_ini, 80 ):
        for nchop in nchop_vals:
            order_nchop_tries.append( ( cost_legfit_recipe(order,nchop),
                                        order, nchop ) )
    #order_nchop_tries.sort()

    for x, ydelta, y0, y90, lbl, resvarname in datasets:
        ok = False
        for _, order, nchop in order_nchop_tries:
            print(f"    Trying legendre {order}-order fit for {resvarname} (nchop={nchop})")
            lf = try_deltalegfit( x, ydelta, y0, y90, order, nchop)
            print(lf['maxdiff'])
            if lf['maxdiff'] < target_prec:
                ok = True
                break
        if not ok:
            raise SystemExit('Failed to find proper fit')
        fitresults.append( ( lf, order, nchop) )
        results[resvarname] = [ float(e) for e in lf['poly_coeffs'] ]
        _print_taylor_code(lf['poly_coeffs'],resvarname=resvarname)
        print(f"DONE at order={order}")
    print()
    print()
    print()

    highest_maxdiff = 0.0
    for (x, ydelta, y0, y90, lbl, resvarname),(lf,order,nchop) in zip(datasets,fitresults):
        highest_maxdiff = max(highest_maxdiff,lf['maxdiff'])
        print("Worst precision of %s : %g   (order=%i, nchop=%i)"%(resvarname,lf['maxdiff'],order,nchop))
    print("Worst precision overall: %g"%highest_maxdiff)

    print()
    for (x, ydelta, y0, y90, lbl, resvarname),(lf,order,nchop) in zip(datasets,fitresults):
        _print_taylor_code(lf['poly_coeffs'],resvarname=resvarname)

    from .json import save_json
    save_json(
        output_filename,
        results,
        force = True
    )

    fig, axs = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    ax, axdiff = axs

    data_colors = []

    for (x, ydelta, y0, y90, lbl, resvarname), (lf,order,nchop) in zip(datasets,fitresults):
        xp = x_to_xprime(x)
        col = ax.plot( xp, ydelta, label=lbl )[0].get_color()
        ax.plot(xp,lf['ydeltafit'],ls='-.',lw=5,color=col)
        axdiff.plot(xp,lf['diff_ydelta'],color=col)
        data_colors.append( col )

    ax.grid()
    axdiff.grid()
    axdiff.semilogy()
    axdiff.set_ylabel('ydeltafit-ydelta')
    ax.legend()
    ax.set_xlabel('xprime')
    ax.set_ylabel('ydelta')
    plt.show()

    fig, axs = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    ax, axdiff = axs
    for (x, ydelta, y0, y90, lbl,resvarname),(lf,_,_),col in zip(datasets,fitresults,data_colors):
        ydeltafit = lf['ydeltafit']
        ax.plot(x,ydelta,label=lbl,color=col)
        ax.plot(x,ydeltafit,ls='-.',lw=5,color=col)
        axdiff.plot(x,lf['reldiff_ydelta_to_smallest'],color=col)

    axdiff.plot( [x[0],x[-1]],[target_prec,]*2, color='black',lw=4, ls=':' )
    ax.grid()
    axdiff.grid()
    ax.semilogx()
    axdiff.semilogy()
    ax.legend()
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    plt.show()

def main_fitglobal( args ):
    mode = args[0] if args else 'missing'
    modes=['primary','scndfresnel','scndlorentz','scndgauss']
    if mode not in modes:
        raise SystemExit('Please provide mode as one of: %s'%(' '.join(modes)))
    if len(args)!=1:
        raise SystemExit('Please only provide a single (mode) argument')
    do_fitglobal(mode)

def do_fitglobal( mode ):
    assert mode in ['primary','scndfresnel','scndlorentz','scndgauss']
    output_filename='global_refitted_classic_curves'
    if mode != 'primary':
        output_filename += f'_{mode}'
    output_filename += '.json'

    from .load import load_thetascan
    from . import curves
    if mode == 'primary':
        FitClassic = curves.ClassicCurve_Primary
    elif mode == 'scndgauss':
        FitClassic = curves.ClassicCurve_ScndGauss
    elif mode == 'scndlorentz':
        FitClassic = curves.ClassicCurve_ScndLorentz
    elif mode == 'scndfresnel':
        FitClassic = curves.ClassicCurve_ScndFresnel
    else:
        assert False

    def flatten_data( thedata ):
        #Flatten data into list of (x, theta, y):
        xvals = thedata['xvals']
        flattened_data = []
        for theta_degree_str, yvals in thedata['theta_2_ypvals'].items():
            theta = float(theta_degree_str)
            for x, y in zip( xvals, yvals ):
                flattened_data.append( (x, theta, y) )
            #    if len(flattened_data)==10:
            #        break#FIXME
            #if len(flattened_data)==10:
            #    break#FIXME
        return np.asarray(flattened_data)

    data = load_thetascan(mode)
    flatdata = flatten_data(data)

    def do_fit_by_fct( flatdata, fitter ):

        #@np.vectorize
        def fitfct( xth, *params ):
            print("FITFCT called with params:",params)
            xarr, tharr  = xth
            return np.vectorize( lambda xval, thval
                                 : fitter( xval, thval, params ) )(xarr, tharr)
        #yvals = flatdata[:,2]
        #errors = 1/(0.0001+abs(yvals-0.5))

        def estimate_ysigma():
            #Make sure the fit won't ignore the tail where y <<1:
            yy = flatdata[:,2]
            ysigma = abs(np.minimum(yy,1.0-yy))
            if mode=='scndlorentz':
                ysigma = ysigma**1.5 #Oddly sensitive here!!
            #if mode=='scndgauss':
            #    ysigma[flatdata[:,0] < 0.5] /=10000
            return ysigma


        ysigma = estimate_ysigma()

        fit_flatdata, fit_ysigma = flatdata, ysigma

        minx = None
        #Fixme: fitter.taylor_cutoff obsolete in this context where we are only
        #fitting the classical recipes??
        if hasattr(fitter,'taylor_cutoff'):
            minx = fitter.taylor_cutoff()
        if minx is not None and minx>0:
            mask = flatdata[:,0] >= minx
            fit_flatdata = fit_flatdata[mask]
            fit_ysigma = fit_ysigma[mask]

        if False:#FIXME
            mask = np.logical_and(flatdata[:,0] >= 0.05,flatdata[:,0] <= 5.0)
            fit_flatdata = fit_flatdata[mask]
            fit_ysigma = fit_ysigma[mask]

        fitparams,fitcov = scipy.optimize.curve_fit( fitfct,
                                                     fit_flatdata[:,:2].T,
                                                     fit_flatdata[:,2],
                                                     sigma=fit_ysigma,
                                                     p0 = fitter.p0(),
                                                     nan_policy='raise')
        corrmat = np.corrcoef(fitcov)
        corrlist = [ (abs(corrmat[i, j]),i,j)
                     for i in range(len(corrmat))
                     for j in range(len(corrmat))
                     if i<j ]
        corrlist.sort()
        for _,i,j in corrlist[::-1][0:10]:
            print( 'Correlation between parameters %i and %i : %.4g'%(i,j,corrmat[i,j]))
        print(fitparams)
        print(fitter.p0())
        return fitparams

    fitpars_classic = do_fit_by_fct( flatdata, FitClassic() )

    from .json import save_json
    def tostdlist_or_none( arr ):
        return [float(e) for e in arr] if arr is not None else None
    save_json(
        output_filename,
        dict( classic = tostdlist_or_none(fitpars_classic) ),
        force = True
    )

def main_fittest( args ):
    #Try to fit new candidates.
    print("WARNING: Not really relevant after this was superseded by the Legendre fit.")
    #Fixme: remove this code?

    mode = args[0] if args else 'missing'
    modes=['primary','scndfresnel','scndlorentz','scndgauss']
    if len(args)>1 or mode not in modes:
        raise SystemExit('Please provide mode as one of: %s'%(' '.join(modes)))

    from .load import load_xscan as dataload
    from .curves import safediv, safesqrt
    from .fit import fitfunction

    val2 = 2.0 if mode!='scndgauss' else 3/math.sqrt(2)

    fcts = []
    @fitfunction
    def f0( x, A, B ):
        k = 1 + B*x
        k2 = 1.0+val2*x+safediv(A*x*x,k)
        return safediv( 1.0, safesqrt( k2 ) )

    fcts = []#fixme
    @fitfunction
    def f1( x, A, B, C ):
        k = 1 + B*x
        k2 = 1.0+val2*x+safediv(A*x*x+C*x,k)
        return safediv( 1.0, safesqrt( k2 ) )

    @fitfunction
    def f2( x, A, B, C, D ):
        b = D*(x/(1.0+x))
        k = 1.0 + B*x +b
        k2 = 1.0+val2*x+safediv(A*x*x+C*x,k)
        return safediv( 1.0, safesqrt( k2 ) )

    @fitfunction
    def f3( x, A, B, C ):
        k = 1.0 + B*x + C*x/(1.0+x)
        k2 = 1.0+val2*x+safediv( A*x*x-0.1*x, k )
        return safediv( 1.0, safesqrt( k2 ) )

    @fitfunction
    def f4( x, A, B, C ):
        b = C*(x/(1.0+x))
        k = 1.0 + B*x + b
        k2 = 1.0+val2*x+safediv(A*x*x,k)
        return safediv( 1.0, safesqrt( k2 ) )


    fcts.append((f0,'orig','green'))
    fcts.append( (f1,f1.name,'blue') )
    fcts.append( (f2,f2.name,'red') )
    fcts.append( (f3,f3.name,'orange') )
    fcts.append( (f4,f4.name,'brown') )

    data = dataload(mode)
    print( data.keys() )
    x = data['xvals']
    ycurves = [ ( f'$\\theta={th:g}\\degree$',
                  ls,
                  lookup_float_key(data['theta_2_ypvals'],th) )
                for th,ls in [
                        (0,':'),
                        (20,'-'),
                        (45,'-'),
                        (70,'--'),
                        (90,'--'),
                ]]

    fig, axs = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    ax, axdiff = axs

    for theta_lbl, ls, y  in ycurves:
        print(theta_lbl)
        ax.plot(x,y,ls=ls,color='black',label=f'reference ({theta_lbl})')
        for f,fname,col in fcts:
            fit = f.fit(x,y)
            print(fname,fit.parameters)
            ax.plot(x,fit(x),ls=ls,color=col,label=f'fit ({fname}, {theta_lbl})')
            axdiff.plot(x,abs((fit(x)-y)/y),ls=ls,color=col,alpha=0.5,lw=3)
    ax.grid()
    ax.semilogx()
    axdiff.semilogy()
    ax.legend()
    ax.set_xlabel('x')
    axdiff.grid()
    plt.show()

def lookup_float_key( adict, key, eps=1e-12 ):
    vals = [ v for k,v in adict.items()
             if abs(float(key)-float(k))<eps ]
    assert len(vals)==1
    return vals[0]
