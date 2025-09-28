import matplotlib.pyplot as plt
from .plotutils import th2color, array
import pandas as pd
import numpy as np
import scipy
import math

kDeg2Rad = np.pi/180
std_modemap = { 'yp':'theta_2_ypvals',
                'yperr':'theta_2_ypvals_maxerr',
                't':'theta_2_calctime'
                }

def std_plots( data,
               mode,
               do_map = False,
               do_cmpwithstd = False,
               finalise_plot = True,
               legend = True ):

    from . import trf
    if mode == 'all':
        for m in std_modemap.keys():
            std_plots(data,m,do_map=do_map)
        return

    data_key = std_modemap[mode]

    xvals_orig = data['xvals']
    xvals = xvals_orig.copy()

    mask = xvals >= 100
    mask = xvals >= -1 # all
    xvals = xvals[mask]

    yvals_std = np.sqrt(xvals_orig*2.0 + 1.0)**(-1.0)
    if do_cmpwithstd:
        yprime_std = trf.transform_to_yprime( yvals_std, xvals_orig )
    def trf_cmpy(y):
        return ( y/yprime_std ) if do_cmpwithstd else y

    def trf_y( y):
        return trf_cmpy(y)
    if do_map:
        xvals_orig = xvals.copy()
        xvals = trf.transform_to_xprime(xvals)
        def trf_y( y):
            return trf_cmpy(trf.transform_to_yprime( y, xvals_orig ))

    for thstr in data['th_keys']:
        thval = float(thstr)
        color = th2color(thstr)
        plt.plot( xvals,
                  trf_y( data[data_key][thstr][mask] ),
                  label = f'$\\theta={thval:g}\\degree$',
                  color = color )
    if finalise_plot:
        plt.ylabel(mode)
        if do_map and mode=='yp':
            plt.ylabel('yprime')
        plt.xlabel('xprime' if do_map else 'x')
        if not do_map:
            plt.semilogx()
        if legend:
            plt.legend()
        plt.grid()
        plt.show()

def main_xscan090( args ):
    return main_xscan(args,is_xscan090=True)

def main_xscan( args, is_xscan090 = False ):
    plotmodes = [e for e in std_modemap.keys()]+['all']
    modes=['primary',
           'scndfresnel',
           'scndlorentz',
           'scndgauss',
           'curve::sabineprimary',
           'curve::sabinescndrec',
           'curve::sabinescndtri',
           'all']
    flagnames = ['map','cmpwithstd']
    def usage():
        print('Please provide data-mode (e.g. "primary") followed by plot-mode (e.g. "all")')
        print('  data-modes: %s'%(' '.join(modes)))
        print('  plot-modes: %s'%(' '.join(plotmodes)))
        print('You can also provide the following additional keywords: %s'%(' '.join(flagnames)))
        raise SystemExit(1)
    flags = set()
    for fn in flagnames:
        while fn in args:
            args.remove(fn)
            flags.add(fn)
    if len(args)!=2:
        usage()
    mode, plotmode = args
    if mode not in modes:
        usage()
    if plotmode not in plotmodes:
        usage()
    from .load import load_xscan, load_xscan090
    loadfct = load_xscan090 if is_xscan090 else load_xscan
    common_args = dict( mode = plotmode,
                        do_map = ( 'map' in flags),
                        do_cmpwithstd = ( 'cmpwithstd' in flags) )
    if mode != 'all':
        std_plots( loadfct(mode), **common_args )
    else:
        actual_modes = [e for e in modes if e!='all']
        for i,m in enumerate(actual_modes):
            std_plots( loadfct(m), **common_args,
                       legend = True,
                       finalise_plot = (len(actual_modes)==i+1 ) )

def _tableN_as_dataframe( thedata ):
    flat = []
    xvals = thedata['xvals']
    for th, ypvals in thedata['theta_2_ypvals'].items():
        sinth = float('%g'%(np.sin(kDeg2Rad*float(th))))
        for x, yp in zip( xvals, ypvals ):
            flat.append( ( x, sinth, yp ) )
    return pd.DataFrame(flat, columns=['x', 'sinth', 'value'])

def _tableN_df_as_table( df ):
    return df.pivot(index='x', columns='sinth', values='value')


def plot_table( args, tablename, load_origptfct, load_allptsfct, origtable_fct ):
    modes = ['all','plot','diff','updated']
    mode = args[0] if (len(args)==1 and args[0] in modes) else None
    if not mode:
        raise SystemExit(f'Please provide {tablename} sub-mode, one of: %s'%(' '.join(modes)))

    data = load_origptfct()
    data_all = load_allptsfct()

    if mode in ('all','plot'):
        std_plots(data_all,'all')
        if mode=='plot':
            return

    orig_table = origtable_fct()
    df_orig = pd.DataFrame(orig_table['x_sinth_yp_list'],
                           columns=['x', 'sinth', 'value'])

    as_table = _tableN_df_as_table
    if mode in ('all','diff'):
        df = _tableN_as_dataframe(data)
        from .print_table1_diff import write_table1_diff_heatmap as ft
        ft( as_table(df_orig), as_table(df), f'{tablename}_diff.html',
            do_print = True )
        ft( as_table(df_orig), as_table(df), f'{tablename}_diff.tex' )

    if mode in ('all','updated'):
        df = _tableN_as_dataframe(data_all)
        from .print_table1_diff import write_updated_table1 as ft
        ft( as_table(df), f'{tablename}_updated.html',do_print = True )
        ft( as_table(df), f'{tablename}_updated.tex' )


def main_bc1974tables( args ):
    assert not args
    from .bc1974_tables import table1, table3, table4
    def doit( tablename, tablefct ):
        t0 = tablefct()
        df = pd.DataFrame(t0['x_sinth_yp_list'],
                          columns=['x', 'sinth', 'value'])
        t = _tableN_df_as_table(df)
        #Follow style of BC1974 Table 1, multiply by 1e4 and show now decimals:
        t *= 10000
        t.columns = ['%4g'%e for e in t.columns]
        t.index = ['%5g'%e for e in t.index]
        t.columns.name = 'sinth'
        t.index.name = '    x'
        print(t.style.format("{:04.0f}").to_string())

    doit('Table 1 (BC1974)', table1)
    doit('Table 3 (BC1974)', table3)
    doit('Table 3 (BC1974)', table4)

def main_table1( args ):
    from .load import ( load_table1scan_origpts,
                        load_table1scan_allpts )
    from .bc1974_tables import table1 as orig_table
    plot_table( args,
                'table1',
                load_table1scan_origpts,
                load_table1scan_allpts,
                orig_table )

def main_table3( args ):
    from .load import ( load_table3scan_origpts,
                        load_table3scan_allpts )
    from .bc1974_tables import table3 as orig_table
    plot_table( args,
                'table3',
                load_table3scan_origpts,
                load_table3scan_allpts,
                orig_table )

def main_table4( args ):
    from .load import ( load_table4scan_origpts,
                        load_table4scan_allpts )
    from .bc1974_tables import table4 as orig_table
    plot_table( args,
                'table4',
                load_table4scan_origpts,
                load_table4scan_allpts,
                orig_table )

def lookup_float_key( adict, key, eps=1e-12 ):
    vals = [ v for k,v in adict.items()
             if abs(float(key)-float(k))<eps ]
    assert len(vals)==1
    return vals[0]

def main_cmprecipes( args ):
    do_strict = True
    while 'nostrict' in args:
        args.remove('nostrict')
        do_strict = False
    do_lux = False
    while 'lux' in args:
        args.remove('lux')
        do_lux = True

    do_xscan090 = False
    while 'xscan090' in args:
        args.remove('xscan090')
        do_xscan090 = True

    mode = args[0] if args else 'missing'
    modes=['primary','scndfresnel','scndlorentz','scndgauss']
    if mode not in modes:
        raise SystemExit('Please provide mode as one of: %s'%(' '.join(modes)))
    args = args[1:]
    assert len(args) in ( 0,1 )
    do_reldiff = True
    if 'abs' in args:
        do_reldiff=False
        args.remove('abs')
    do_vs_old = False
    if 'vsold' in args:
        do_vs_old = True
        args.remove('vsold')
    assert not args
    do_cmprecipes( do_reldiff = do_reldiff,
                   do_vs_old = do_vs_old,
                   mode = mode,
                   do_strict = do_strict,
                   do_lux = do_lux,
                   do_xscan090 = do_xscan090 )

def do_cmprecipes( do_reldiff, do_vs_old, mode, do_strict, do_lux, do_xscan090  ):
    assert mode in ['primary','scndfresnel','scndlorentz','scndgauss']
    from . import curves
    ( ClassicCurve,
      UpdatedClassicCurve,
      ProposedCurve,
      ProposedLuxCurve ) = curves.mode_curves(mode)
    if do_lux:
      ProposedCurve = ProposedLuxCurve
    del ProposedLuxCurve

    from .load import load_xscan090, load_xscan

    if do_vs_old:
        assert do_reldiff, "vsold only for reldiff"
    data = (load_xscan090 if do_xscan090 else load_xscan)(mode=mode)
    #print( data.keys() )
    xvals = data['xvals']

    def fleq( a, b, eps=1e-12 ):
        return abs(float(a)-float(b))<eps

    th2yp = dict( (th,ypvals)
                  for th,ypvals in data['theta_2_ypvals'].items() )
    worst_reldiff_classic = [1e-99]
    worst_reldiff_updatedclassic = [1e-99]
    worst_reldiff_proposed = [1e-99]
    if do_reldiff:
        def ytrf( yp, ypref, worst_reldiff_obj = None ):
            if do_strict:
                #Relative difference to smallest of y and 1-y:
                v = ( abs(yp-ypref)/np.minimum(abs(ypref),abs(1.0-ypref)) )
            else:
                #Relative difference to absolute y value
                v = abs(yp/ypref-1.0)
            v = np.clip(v,0.0,1.0)
            if worst_reldiff_obj is not None:
                worst_reldiff_obj[0] = max(v.max(),worst_reldiff_obj[0])
            return v
    else:
        def ytrf( yp, ypref ):
            return np.clip(yp,0.0,1.0)

    seen_lbls = set()
    for th,yref in th2yp.items():
        th = float(th)
        color = th2color(th)
        yp_classic = ClassicCurve()(xvals,th)
        yp_updatedclassic = UpdatedClassicCurve()(xvals,th) if UpdatedClassicCurve else None
        yp_proposed = ProposedCurve()(xvals,th) if ProposedCurve else None
        if do_vs_old:
            yref = yp_classic

        common = dict( alpha = 0.4,
                       lw = 2 )
        if abs(th-0.0)<1e-7 or abs(th-90.0)<1e-7:
            common['lw'] = 5
            common['alpha'] = 0.8

        def lbltrf(curve_type, split_theta = None):
            lbl = curve_type
            if split_theta is not None:
                if th<split_theta:
                    lbl += f' ($\\theta<{split_theta:g}\\degree$)'
                else:
                    lbl += f' ($\\theta\\geq{split_theta:g}\\degree$)'
            if lbl in seen_lbls:
                return
            seen_lbls.add(lbl)
            return lbl
            if not any( fleq(e,th) for e in [90] ):
                return None
            thstr = f'$\\theta={th:g}\\degree$'
            return f'{thstr} ({curve_type or "reference"})'
        if not do_reldiff:
            #color='black'
            plt.plot( xvals, yref,
                      **common,
                      color=color,
                      label = lbltrf('') )
        #color='red'
        color = 'blue'#th2color(th*0.3,'Reds')
        split_th = 60
        if th > split_th:
            color='green'
        if not do_vs_old:
            plt.plot( xvals, ytrf(yp_classic,yref,worst_reldiff_classic),
                      **common,
                      color=color,
                      ls = '-.',
                      label = lbltrf('BC1974',split_th) )
        #color='green'
        color = 'red'#th2color(th*0.3,'Greens')
        if do_vs_old and th > split_th:
            color='green'
        if yp_updatedclassic is not None:
            plt.plot( xvals, ytrf(yp_updatedclassic,yref,worst_reldiff_updatedclassic),
                      **common,
                      color=color,
                      ls = '--',
                      label = lbltrf('updated BC1974',
                                     split_theta = split_th if do_vs_old else None) )


        #color='blue'
        color = 'black'#th2color(th*0.2,'Blues')
        if do_vs_old and th > split_th:
            color='gray'
        #color = th2color(th)
        if yp_proposed is not None:
            #bla = ytrf(yp_proposed,yref,worst_reldiff_proposed)
            #print( 'BLA',th,xvals[bla>1e-5],yref[bla>1e-5])
            plt.plot( xvals, ytrf(yp_proposed,yref,worst_reldiff_proposed),
                      **common,
                      color=color,
                      ls = '-',
                      label = lbltrf('new form',
                                     split_theta = split_th if do_vs_old else None) )
    yvalname = 'yp' if mode=='primary' else 'ys'

    if do_reldiff:
        plt.semilogy()
        plt.ylim(1e-10,1.0)
        if do_strict:
            plt.ylabel(f'max({yvalname}/{yvalname}ref-1,(1-{yvalname})/(1-{yvalname}ref)-1)')
        else:
            plt.ylabel(f'{yvalname}/{yvalname}ref-1')
    else:
        plt.ylim(0.0,1.0)
        plt.ylabel(yvalname)
    plt.title(mode)
    plt.xlim(xvals[0],xvals[-1])
    plt.grid()
    plt.semilogx()
    plt.legend()
    plt.xlabel('x')
    print("Worst reldiff (CLASSIC)    : %g"%worst_reldiff_classic[0])
    print("Worst reldiff (UpdCLASSIC) : %g"%worst_reldiff_updatedclassic[0])
    print("Worst reldiff (Proposed)   : %g"%worst_reldiff_proposed[0])

    for cn,curve in [ ('classic',ClassicCurve),
                      ('updclassic',UpdatedClassicCurve),
                      ('proposed',ProposedCurve) ]:
        if curve is None:
            continue
        print('%20s at x=1e6 and theta=0,45,90: %g   %g   %g'
              %(cn, curve()(1e6,0), curve()(1e6,45), curve()(1e6,90) ) )



    plt.show()

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

def main_par2latex( args ):
    assert len(args)==0
    from .curves import load_fitted_curve_parameters
    data = load_fitted_curve_parameters('primary')
    raise SystemExit("FIXME: This needs update for the new proposed curves")

    print(r'\begin{align}')
    pA0, pA1 = data['origA']
    pB0, pB1 = data['origB']
    sA0 = '%g'%pA0
    sB0 = '%g'%pB0
    sA1 = '%g'%pA1
    sB1 = '%g'%pB1
    if not sA1.startswith('-'):
        sA1 = '+' + sA1
    if not sB1.startswith('-'):
        sB1 = '+' + sB1
    print(r'  A(\theta) &= %s%s\cos(2\theta)\nonumber\\'%(sA0,sA1))
    print(r'  B(\theta) &= %s%s(0.5-\cos(2\theta))^2'%(sB0,sB1))
    print(r'\labeqn{ypfitfctparsupdated}')
    print(r'\end{align}')
    print()
    print(r'\begin{align}')
    for pn in 'ABC':
        pars = data[pn]
        s = r'  %s(\theta) &= '%pn
        for i,p in enumerate(pars):
            f = '%.5g'%p
            if i==0:
                s += f
                continue
            if i==4:
                print(s+r'\nonumber\\')
                s=r'             &\hspace*{1em} '
            if not f.startswith('-'):
                f = f'+{f}'
            if i==1:
                s += r'%s\sin(\theta)'%f
            else:
                s += r'%s\sin^%i(\theta)'%(f,i)
        if pn!='C':
            s += r'\nonumber\\'
        print(s)
    print(r'\labeqn{ypfitfctluxpars}')
    print(r'\end{align}')


@np.vectorize
def find_taylor_maxx( theta, eps ):
    from .eq36_lowx_taylor import taylor_lowx_eq36
    def is_ok( x ):
        return taylor_lowx_eq36( theta, x, eps ) is not None

    xmin, xmax = 1e-20, 100.0
    assert is_ok(xmin) and not is_ok(xmax)

    target_precision = 1e-7
    while True:
        xmid = (xmax+xmin)*0.5
        if is_ok(xmid):
            xmin = xmid
        else:
            xmax = xmid
        if (xmax-xmin) < target_precision:
            return 0.5*(xmax+xmin)

def main_taylorreach( args ):
    print("WARNING: This is only for the primary mode right now (and does"
          " not show reach of taylor user recipes)")
    n = 1000
    split_x = 5
    theta_vals = np.concat( ( np.linspace(0,split_x,n//2,endpoint=False),
                              np.linspace(split_x,90,n//2+1) ) )
    for eps in [1e-3, 1e-8, 1e-16]:
        print("Producing data for eps=%g"%eps)
        plt.plot( theta_vals,
                  find_taylor_maxx( theta_vals, eps ),
                  label = 'eps=%g'%eps)
    plt.grid()
    plt.xlabel(r'$\theta$ ($\degree$)')
    plt.ylabel('x')
    plt.title('Maximal x value where eq. 36 can be evaluated via Taylor expansion')
    plt.legend()
    plt.show()

def main_spread( args ):
    mode = args[0] if args else 'missing'
    modes=['primary','scndfresnel','scndlorentz','scndgauss']
    if mode not in modes:
        raise SystemExit('Please provide mode as one of: %s'%(' '.join(modes)))
    if len(args)!=1:
        raise SystemExit('Please only provide a single (mode) argument')
    do_spread(mode)

def do_spread( mode ):
    #Fixme: just plot all modes in same plot.
    from .curves import mode_curves
    _,_,_,ProposedLuxCurve = mode_curves(mode)
    from .load import load_xscan as dataload
    data = dataload(mode)
    xvals = data['xvals']
    def key_theta( theta ):
        keys = [ k for k in data['th_keys']
                 if abs(float(k)-float(theta))<1e-12 ]
        assert len(keys)==1
        return keys[0]
    yp0 = data['theta_2_ypvals'][key_theta(0)]
    #yp45 = data['theta_2_ypvals'][key_theta(0)]
    yp90 = data['theta_2_ypvals'][key_theta(90)]
    def reldiff(yp90,yp0):
        return (yp90-yp0)/(0.5*(yp90+yp0))
    plt.plot( xvals, reldiff(yp90,yp0), label='numerical integration' )
    plt.plot( xvals,
              reldiff(ProposedLuxCurve()(xvals,90),ProposedLuxCurve()(xvals,0)),
              label='Via new parameterisation' )
    #plt.semilogx()
    plt.loglog()
    plt.grid()
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('Max relative spread in yp values for different theta')
    plt.show()

def find_nearest_idx( arr, val):
    assert isinstance(arr, np.ndarray)
    idx = (np.abs(arr - val)).argmin()
    assert idx < len(arr)
    return idx

def plot_breakdown( mode, curve, ref ):
    assert curve in ('classic','updatedclassic','proposed','proposedlux')
    assert ref in ('numint','proposedlux')
    assert curve != ref

    print(f"Performing reldiff 2D plot of {curve} vs {ref} for {mode}")
    from .curves import mode_curves

    ( ClassicCurve,
      UpdatedClassicCurve,
      ProposedCurve,
      ProposedLuxCurve ) = mode_curves(mode)

    if ref=='proposedlux' and ProposedLuxCurve is None:
        print("Skipping plot breakdown vs. proposed lux curve since it is not yet ready")
        return
    if curve=='updatedclassic' and UpdatedClassicCurve is None:
        print("Skipping plot breakdown of updated classic curve since it is not yet ready")
        return
    if curve=='proposed' and ProposedCurve is None:
        print("Skipping plot breakdown of proposed curve since it is not yet ready")
        return
    if curve=='proposedlux' and ProposedLuxCurve is None:
        print("Skipping plot breakdown of proposed curve since it is not yet ready")
        return

    from .load import load_xscan as dataload

    if ref=='numint':
        refdata = dataload(mode)
        x = refdata['xvals']
        th = array(sorted(float(e) for e in refdata['th_keys']))
        _last_lookup = [None,None]
        def lookup_ref(xval,thval):
            if _last_lookup[0] != thval:
                for _th,_yp in refdata['theta_2_ypvals'].items():
                    if abs(float(_th)-thval)<1e-12:
                        _last_lookup[0] = thval
                        _last_lookup[1] = _yp
            assert _last_lookup[0] == thval
            return _last_lookup[1][find_nearest_idx(x,xval)]
    else:
        assert ref == 'proposedlux'
        x = np.geomspace(0.1, 100.0, 400)  # Adjust the range and resolution as needed
        th = np.linspace(0.1, 100.0, 400)
        ref_curve = ProposedLuxCurve()
        def lookup_ref(xval, thval):
            return ref_curve(xval,thval)

    X, TH = np.meshgrid(x, th)
    if curve=='classic':
        curvefct = ClassicCurve()
    elif curve=='updatedclassic':
        curvefct = UpdatedClassicCurve()
    elif curve=='proposed':
        curvefct = ProposedCurve()
    else:
        assert curve=='proposedlux'
        curvefct = ProposedLuxCurve()

    @np.vectorize
    def z_calc( xval, thval ):
        if hasattr(curve,'breakdown') and curvefct.breakdown(xval, thval):
            return 1.0
        y = curvefct(xval, thval)
        yref = lookup_ref(xval,thval)
        return max(0.0,min(1.0,abs((y-yref)/yref)))

    map_top = ( 25, 30 )
    @np.vectorize
    def val_fix( val ):
        val = 100.0*val#percentage
        a,newb = map_top
        b = 100
        if val > a:
            rel = (val-a)/(b-a)
            val = a + rel*(newb-a)
        return max(0.0,min(float(newb),val))

    #Custom colour map:
    import matplotlib.colors as mcol

    colors = [(0,0,0.1), (0,0,0.3), (0,0,0.35),   (0,0,0.7),    (0,0.7,0), (0,0.3,0),     "yellow",  "orange", "red"]
    positions = [0,      2.0,       2.1,      5-0.5,       5 + 0.5,         10-0.5,       10+0.5,     25, 100]

    if True:
        colors = ["darkblue", "blue", "green", "darkgreen", "yellow", "red", "black"]
        _eps = 1.0
        positions = [0,        5-_eps,       5 + _eps,         10-_eps,       10+_eps,     25, 100]



    positions = [ val_fix(e/100.0)/val_fix(100/100.0) for e in positions ]
    positions[-1] = 1.0
    cmap = mcol.LinearSegmentedColormap.from_list("custombreakdowncmap",
                                                  list(zip(positions,
                                                           colors)))


    Z = z_calc(X, TH)
    Z = val_fix(Z)
    #Z = np.clip(Z,0.0,20.0)

    plt.pcolormesh(x, th, Z,
                   cmap=cmap, vmin=0.0, vmax=map_top[-1],
                   edgecolors='none', rasterized=True,

                   )#, shading='auto')#, cmap='binary')
    plt.semilogx()
    cbar = plt.colorbar(label='Precision (%)')

    assert map_top[-1]==30
    yticks = [0,2,5,10,15, 20, 25, map_top[-1]]
    cbar.ax.set_yticks(yticks)
    yticklbls = [ str(e) for e in yticks ]
    yticklbls[-1] = str(100)+'+'
    cbar.ax.set_yticklabels(yticklbls)

    #plt.title('Relative model error')
    plt.grid()
    plt.ylim(0.0,90.0)
    plt.xlim(x[0],x[-1])
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\theta$')
    plt.savefig('%s_vs_%s.pdf'%(curve,ref),format='pdf',dpi=300)
    plt.show()

def main_breakdown( args ):
    mode = args[0] if args else 'missing'
    modes=['primary','scndfresnel','scndlorentz','scndgauss']
    if len(args)>1 or mode not in modes:
        raise SystemExit('Please provide mode as one of: %s'%(' '.join(modes)))

    plot_breakdown( mode, 'proposedlux', ref='numint')
    plot_breakdown( mode, 'classic', ref='proposedlux')
    #plot_breakdown( mode, 'classic', ref='numint')
    #plot_breakdown( mode, 'updatedclassic', ref='numint')
    plot_breakdown( mode, 'updatedclassic', ref='proposedlux')
    plot_breakdown( mode, 'proposed', ref='proposedlux')

def main_highx( args ):
    mode = args[0] if args else 'missing'
    modes=['primary','scndfresnel','scndlorentz','scndgauss']
    if mode not in modes:
        raise SystemExit('Please provide mode as one of: %s'%(' '.join(modes)))
    if len(args)!=1:
        raise SystemExit('Please only provide a single (mode) argument')
    do_highx(mode)

def do_highx( mode ):
    assert mode in ['primary','scndfresnel','scndlorentz','scndgauss']

    from .load import load_xscan as load
    data = load(mode)
    xvals = data['xvals']
    mask = xvals>=100.0
    xvals = xvals[mask]

    def fpow(x,norm,power):
        return norm * (x**power)

    for thstr in data['th_keys']:
        thval = float(thstr)
        color = th2color(thstr)
        yvals = data['theta_2_ypvals'][thstr][mask]
        plt.plot( xvals,
                  yvals,
                  label = f'$\\theta={thval:g}\\degree$',
                  color = color )
        do_fit = any(abs(thval-e)<1e-3 for e in [0,45,90])
        if do_fit:
            res,_ = scipy.optimize.curve_fit( fpow,
                                              xvals, yvals,
                                              p0 = [1,-0.5], nan_policy='raise')
            print(res[1])
            plt.plot(xvals,fpow(xvals,*res),ls=':',lw=3,color=color)

    plt.ylabel(mode)
    plt.xlabel('x')
    plt.legend()
    plt.grid()
    plt.loglog()
    plt.show()

def main_printcpp( args ):
    raise SystemExit("FIXME: Code needs to be updated for new legendre fits")
    mode = args[0] if args else 'missing'
    modes=['primary','scndfresnel','scndlorentz','scndgauss']
    if len(args)>1 or mode not in modes:
        raise SystemExit('Please provide mode as one of: %s'%(' '.join(modes)))
    from .curves import load_fitted_curve_parameters
    d = load_fitted_curve_parameters(mode)
    print(d)
    dp = d['proposed']
    npars = 7
    varnames = 'ABC' if len(dp)==21 else 'AB'
    assert len(dp)==len(varnames)*npars
    for i,v in enumerate(varnames):
        dv = dp[i:i+npars]
        assert len(dv) == npars
        for i,pi in enumerate(dv):
            print(f'constexpr double c{v}{i} = %g;'%pi)
    print('const double s = sintheta;')
    for v in varnames:
        c = f'c{v}'
        s = f'{c}0'
        for i in range(1,npars):
            s+=f'+s*({c}{i}'
        s+=')'*(npars-1)
        print(f'const double {v} = {s};')
    print()
    for i,v in enumerate(varnames):
        pars = dp[i:i+npars]
        assert npars == len(pars)
        s = '%g'%pars[0]
        for pv in pars[1:-1]:
            s+='+s*(%g'%pv
        pv ='%g'%pars[-1]
        if not pv.startswith('-'):
            pv=f'+{pv}'
        s+='%s*s'%pv
        s+=')'*(npars-2)
        print(f'const double {v} = {s};')

def main_investigatefcts( args ):
    print("FIXME: This mode might be outdated")
    mode = args[0] if args else 'missing'
    modes=['primary','scndfresnel','scndlorentz','scndgauss']
    if mode not in modes:
        raise SystemExit('Please provide mode as one of: %s'%(' '.join(modes)))
    if len(args)!=1:
        raise SystemExit('Please only provide a single (mode) argument')
    from .load import load_xscan as load
    data = load(mode)
    is_lorentz = (mode=='scndlorentz')
    is_gauss = (mode=='scndgauss')
    xvals = data['xvals']
    fitmask = xvals > 0.1 #focus on area where taylor expansion is not enough

    from .curves import safediv, safesqrt

    #TEST FIT FCT:
    fitnpars = 5
    p0 = [ 0.0, 1.0, 0.0, 0.0, 0.0 ]
    def fitfct( x, A, B, C, D, E ):
        powval = 2.87 if is_gauss else 2.0
        def step(v):
            return v*v / (1+v*v)
        k = 1.0 + (2+E)*x  + ( (1+A)*x**powval + D*x) / (1+B*x+C*step(x))
        return safediv( 1.0, safesqrt( k ))

    fitnpars_lorentz = 3
    def fitfct_lorentz( x, A, B, C ):
        k = ( 1.0
              + 2*x
              + ( A*x**(2+C) ) / (1+B*x)
              )
        return safediv( 1.0, safesqrt( k ))


    if is_lorentz and False:
        fitfct = fitfct_lorentz
        fitnpars = fitnpars_lorentz


    ############### TEMP TRY:
    fitnpars = 3
    def fitfct( x, A, B, C ):
        #powval = 2.87 if is_gauss else 2.0
        #def step(v):
        #    return v*v / (1+v*v)
        factor1 = safediv( 1.0, safesqrt( 1.0 + A * x ))
        factor2 = safediv( ( 1.0 ), ( 1.0 + C*x**(B-0.5)) )
        return factor1*factor2

    dummyx = np.geomspace(1e-6,1e6,1000)
    plt.plot( dummyx,
              np.vectorize( lambda x : fitfct(x,2.0,0.9,3.0) )(dummyx) )
    plt.semilogx()
    plt.grid()
    plt.show()
    ################################################################

    fig, axs = plt.subplots(2, 1, sharex=True)
    fig.subplots_adjust(hspace=0)
    ax, axdiff = axs

    worst_rd = -1.0
    for thstr in data['th_keys']:
        thval = float(thstr)
        if not any(abs(thval-e)<1e-3 for e in [0,30,60,90]):
            continue
        color = th2color(thstr)
        yvals = data['theta_2_ypvals'][thstr]
        ax.plot( xvals,
                 yvals,
                 label = f'$\\theta={thval:g}\\degree$',
                 color = color )
        #Fit:
        xfit, yfitdata = xvals[fitmask], yvals[fitmask]
        ysigma = abs(np.minimum(yfitdata,1.0-yfitdata))
        if is_lorentz:
            ysigma = ysigma**1.5 #Oddly sensitive here!!

        p0=[0.0,]*fitnpars

        def objective( params ):
            abs(np.minimum(yfitdata,1.0-yfitdata))
            yfit=np.vectorize(lambda xval : fitfct(xval,*params))(xfit)
            rd = abs(yfit-yfitdata)/(np.minimum(yfitdata,1.0-yfitdata))
            return rd.max()


        if False:
            print("Annealing...")
            anneal = scipy.optimize.dual_annealing(objective, [[-10,10],]*fitnpars)
            print(anneal)
            p0 = [e for e in anneal.x]



        def do_fit(p0):
            return scipy.optimize.curve_fit( fitfct, xfit, yfitdata,
                                             sigma = ysigma,
                                             p0 = p0,
                                             nan_policy='raise',
                                             maxfev = 10000,
                                            )[0]
        print("Fit1...")
        res = do_fit(p0)
        def fmt(a):
            return [float('%g'%e) for e in a]
        print("Fit2...")
        res = do_fit(res)
        print("Fit3...")
        res = do_fit(res)
        print("Final result:",fmt(res))
        yfit = fitfct(xfit,*res)
        ax.plot(xfit,yfit,ls=':',lw=3,color=color)
        rd = abs(yfit-yfitdata)/(np.minimum(yfitdata,1.0-yfitdata))
        worst_rd = max(worst_rd,rd.max())
        axdiff.plot(xfit,rd,ls=':',lw=3,color=color)

    print("WORST RELDIFF PERFORMANCE: %.3g"%worst_rd)
    ax.grid()
    ax.loglog()
    axdiff.semilogy()
    ax.legend()
    ax.set_xlabel('x')
    axdiff.grid()
    ax.set_xlim(0.05)
    plt.show()

#Fixme: to utils.py:
def print_taylor_code( poly_coefficients, varname='xp', mode='py',
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

def legfit_weight_fct(x,y,xp,yp):
    return ( 1.0/np.minimum(abs(y),abs(1.0-y)) )**0.3 #Fixme: revisit this, make it easy to describe in paper

def main_legfit( args ):
    from .load import load_xscan090
    do_lux = False
    while 'lux' in args:
        args.remove('lux')
        do_lux = True
    assert len(args)==0
    modes=['primary',
           'scndfresnel',
           'scndlorentz',
           'scndgauss',
           #'curve::sabineprimary',
           #'curve::sabinescndrec',
           #'curve::sabinescndtri',
           ]
    thvals=[0,90]#Due to final interpolation, only these matters
    datasets = []
    for m in modes:
        data = load_xscan090(m)
        for th, y in data['theta_2_ypvals'].items():
            if thvals and not any( abs(float(th)-e)<1e-6 for e in thvals ):
                continue
            assert int(float(th))==float(th)
            datasets.append( ( data['xvals'].copy(),
                               y.copy(),
                               f'{m} ({float(th):g} deg)',
                               f'yprime_{m}_{int(float(th))}' ) )
    do_legfit( datasets,
               do_lux=do_lux )

def try_legfit( x, y, order, nchop  ):
    from .trf import xy_prime, xy_unprime

    xp, yp = xy_prime( x, y )
    fit_weights = legfit_weight_fct(x,y,xp,yp)
    legendre_coefs = np.polynomial.legendre.legfit( xp, yp,
                                                    order,
                                                    w = fit_weights )
    poly_coeffs = np.polynomial.legendre.leg2poly(legendre_coefs)

    #Chop values (and remove highest orders if necessary):
    def do_chop( iorder, val ):
        return float((f'%.{nchop}g')%val)
    poly_coeffs = [do_chop(i,val) for i,val in enumerate(poly_coeffs)]

    while not poly_coeffs[-1]:
        poly_coeffs = poly_coeffs[:-1]

    ypfit = np.polynomial.polynomial.polyval(xp, poly_coeffs)
    _, yfit = xy_unprime( xp, ypfit )
    diff_yp = abs(ypfit-yp)
    diff_y = ( abs(yfit-y)/np.minimum(abs(y),abs(1-y)) )
    maxdiff = diff_y.max()
    return dict( maxdiff = maxdiff,
                 diff_y = diff_y,
                 diff_yp = diff_yp,
                 poly_coeffs = poly_coeffs,
                 xp = xp,
                 yp = yp,
                 ypfit = ypfit,
                 yfit = yfit )


def cost_legfit_recipe( order, nchop ):
    #Estimate "length" of recipe, when written as c0+xp*(c1+xp*(...)).
    #Penalise slightly for order, since higher order also means more multiplications
    str_length_approx = (10+nchop)*order
    return str_length_approx + order*3

def do_legfit( datasets, do_lux ):
    from .trf import xy_prime, xy_unprime

    target_prec = 1e-7 if do_lux else 1e-3
    target_prec *= 0.9 #safety (but do NOT decrease lux lvl further, since ref
                       #data can not handle it).

    output_filename='legendre_coefficients'
    if do_lux:
        output_filename += '_lux'
    output_filename += '.json'

    results = {}
    fitresults = []

    order_nchop_tries = []
    order_ini, nchop_ini = (7, 8) if do_lux else (4,3)
    nchop_vals = list(range(nchop_ini,16))
    for order in range( order_ini, 60 ):
        for nchop in nchop_vals:
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
        fitresults.append( lf )
        results[resvarname] = [ float(e) for e in lf['poly_coeffs'] ]
        print_taylor_code(lf['poly_coeffs'],resvarname=resvarname)
        print(f"DONE at order={order}")
    print()
    print()
    print()

    highest_maxdiff = 0.0
    for (x, y, lbl, resvarname),lf in zip(datasets,fitresults):
        highest_maxdiff = max(highest_maxdiff,lf['maxdiff'])
        print("Worst precision of %s : %g"%(resvarname,lf['maxdiff']))
    print("Worst precision overall: %g"%highest_maxdiff)

    print()
    for (x, y, lbl, resvarname),lf in zip(datasets,fitresults):
        print_taylor_code(lf['poly_coeffs'],resvarname=resvarname)

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

    for (x, y, lbl, resvarname), lf in zip(datasets,fitresults):
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
    for (x, y, lbl,resvarname),lf,col in zip(datasets,fitresults,data_colors):
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
