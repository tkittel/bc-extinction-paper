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
    norm_is_min_y_1minusy = True
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
        ft(
            as_table(df_orig), as_table(df), f'{tablename}_diff.html',
            do_print = True,
            norm_is_min_y_1minusy = norm_is_min_y_1minusy
           )
        ft(
            as_table(df_orig), as_table(df), f'{tablename}_diff.tex',
            norm_is_min_y_1minusy = norm_is_min_y_1minusy
        )

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

def main_recipevstaylor( args ):
    from .eq36_lowx_taylor import mode_taylorfct
    from .curves import mode_curves
    from .mpmath import mpf

    do_oneminusy = True
    do_showrelprec = True

    taylor_eps = 1e-20

    assert not args
    xvals = np.geomspace( 1e-8, 1.8e-1, 500 )
    thetavals = np.asarray( [0,
                             #15, 30, 45, 60, 75, 90
                             ], dtype=float )
    modes = ['primary','scndfresnel','scndlorentz','scndgauss']

    curve_lbls_col = [ ( 'BC1974', 'blue' ),
                       ( 'updated BC1974', 'red' ),
                       ( 'Proposed', 'black' ),
                       ( 'Proposed (Lux)', 'green' ) ]

    def mapy( y ):
        return y
    if do_oneminusy:
        def mapy( y ):
            return mpf(1)-y

    for imode,mode in enumerate(modes):
        for itheta,theta in enumerate(thetavals):
            def taylor( x ):
                return mode_taylorfct(mode)( theta, x, eps = taylor_eps )[0]
            yref = np.vectorize(lambda xx : mapy(taylor(xx)))(xvals)
            curveobjs = [c() for c in mode_curves(mode)]
            if not do_showrelprec:
                plt.plot( xvals, yref, color='gray', label='ref', ls=':',lw=4 )
            for curvefct, (lbl,col) in zip(curveobjs,curve_lbls_col):
                print (curvefct.__class__)
                if itheta > 0 or imode > 0:
                    lbl = None
                y = np.vectorize(lambda xx : mapy(curvefct(xx,theta)))(xvals)
                yprec = abs(y-yref) / (1e-300 + np.minimum(yref,1.0-yref) )
                print("itheta,lbl = ",itheta,lbl)
                if do_showrelprec:
                    plt.plot( xvals, yprec, color=col, label=lbl )
                else:
                    plt.plot( xvals, y, color=col, label=lbl )

    plt.semilogx()
    plt.semilogy()
    plt.grid()
    plt.legend()
    plt.xlabel('x')
    plt.show()

def main_taylororders( args ):
    from .new_recipes import recipe_target_prec_each_main_term

    ignore_c0 = None
    while '1minusy' in args:
        args.remove('1minusy')
        ignore_c0 = True
    while 'raw' in args:
        args.remove('raw')
        ignore_c0 = False
    if args or ignore_c0 is None:
        raise SystemExit('Please specify either "raw" or "1minusy" (and no other args).')

    #ignore_c0 means that we do all comparisions with ignore_c0 = True, to
    #leave out the leading c0=1 in the Taylor expansions.

    from .eq36_lowx_taylor import ( mode_taylorfct,
                                    create_taylorfct_lowx_eq36_limited_order )

    modes = ['primary','scndfresnel','scndlorentz','scndgauss']
    xvals = np.geomspace( 1e-5 if ignore_c0 else 1e-17,
                          0.2,
                          2000//10 )#FIXME
    thetavals = np.asarray( [0, 15, 30, 45, 60, 75, 90], dtype=float )
    orders = list(range(1,8+1)) + [12, 16 ]#, 20]

    min_eps = 1e-17

    mode2yreflist = {}
    for mode in modes:
        print(f"Preparing reference values for mode {mode}")
        yref_list = []
        raw_highresfct = mode_taylorfct(mode)
        for th in thetavals:
            @np.vectorize
            def highres_fct( x ):
                eps = 1e-20 if (mode!='scndfresnel' or x<0.01) else 1e-4
                eps = 1e-4
                v = raw_highresfct(x=x,theta_degree=th,
                                   eps=eps,
                                   ignore_c0=ignore_c0)
                if v is None:
                    raise SystemExit(f'Failed at mode={mode}, theta={th} and x={x}')
                return float(v[0])
            yref_list.append( highres_fct( xvals ) )
        mode2yreflist[mode] = yref_list


    plt.plot( [ xvals[0], xvals[-1] ],
              [ recipe_target_prec_each_main_term(lux=True), ]*2,
              color='black',ls=':',lw=3)
    plt.plot( [ xvals[0], xvals[-1] ],
              [ recipe_target_prec_each_main_term(lux=False), ]*2,
              color='black',ls=':',lw=3)
    from .new_recipes import recipe_taylor_cutoff as tc
    plt.plot( [ tc(lux=False), ]*2, [ min_eps,1.0, ], color='black',ls=':',lw=3)
    plt.plot( [ tc(lux=True), ]*2, [ min_eps,1.0, ], color='black',ls=':',lw=3)

    for i,order in enumerate(orders):
        #Find worst precision for any theta and any of the four modes:
        rd_worst = None
        for mode in modes:
            yref_list = mode2yreflist[mode]
            for th, yref in zip(thetavals,yref_list):
                sinth = float(math.sin(th * np.pi/180 ))
                f = create_taylorfct_lowx_eq36_limited_order( mode = mode,
                                                              sinth=sinth,
                                                              order=order,
                                                              ignore_c0=ignore_c0 )
                y = f(xvals)
                rd = np.clip(abs(yref-y)/np.maximum(1e-300,np.minimum(abs(yref),abs(1.0-yref))),min_eps*1e-5,1.0)
                if rd_worst is None:
                    rd_worst = rd
                else:
                    rd_worst = np.maximum(rd_worst,rd)
        plt.plot( xvals, rd_worst, label=f'Taylor order {order}', lw=3, alpha=1.0 )

    plt.ylabel('Worst relative precision of 1-y(x) for any model and theta_bragg')
    plt.semilogx()
    plt.semilogy()
    plt.ylim(min_eps,1.0)
    plt.xlim(xvals[0],0.2)#xvals[-1])#fixme?
    plt.grid()
    plt.legend()
    plt.xlabel('x')
    plt.show()

def main_par2latex( args ):
    assert len(args)==0
    from .curves import mode_curves
    for mode in ['primary','scndfresnel','scndlorentz','scndgauss']:
        if not hasattr(mode_curves(mode)[1](),'params2latex'):
            print("WARNING: Skipping %s for now"%mode)
            continue
        p2l = mode_curves(mode)[1]().params2latex()
        print()
        print(r'\begin{align}')
        for i,(k,v) in enumerate(p2l):
            is_last = (i+1==len(p2l))
            print(r'  %s &= %s%s\\'%(k,v,('' if is_last else r'\nonumber')))
        print(r'  \labeqn{refitted1974%s}'%mode)
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
    from .curves import mode_curves
    from .load import load_xscan090 as dataload
    show_ref_numint = False
    while 'ref' in args:
        args.remove('ref')
        show_ref_numint = True
    assert not args
    modes=['primary','scndfresnel','scndlorentz','scndgauss']
    for mode in modes:
        def reldiff(_yp90,_yp0):
            return (_yp90-_yp0)/(0.5*(_yp90+_yp0))
        if show_ref_numint:
            data = dataload(mode)
            xvals = data['xvals']
            def key_theta( theta ):
                keys = [ k for k in data['th_keys']
                         if abs(float(k)-float(theta))<1e-12 ]
                assert len(keys)==1
                return keys[0]
            yp0 = data['theta_2_ypvals'][key_theta(0)]
            yp90 = data['theta_2_ypvals'][key_theta(90)]
            plt.plot( xvals, reldiff(yp90,yp0), label=f'{mode} reference integration' )

        _,_,_,ProposedLuxCurve = mode_curves(mode)
        xvals_recipe = np.geomspace(1e-5,1e5,10000)
        plt.plot( xvals_recipe,
                  reldiff( ProposedLuxCurve()(xvals_recipe,90),
                           ProposedLuxCurve()(xvals_recipe,0)),
                  label=f'{mode} (BC2025 Lux recipe)' )
    plt.loglog()
    plt.grid()
    plt.legend()
    plt.xlabel('x')
    plt.ylabel('Max relative spread in y values for different theta')
    plt.show()

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
        def find_nearest_idx( arr, val):
            assert isinstance(arr, np.ndarray)
            idx = (np.abs(arr - val)).argmin()
            assert idx < len(arr)
            return idx
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
        #Fixme: higher res for paper:
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
        return max(0.0,min(1.0,float(abs((y-yref)/yref))))#FIXME: DIFFERENT PRECISION DEFINITION

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

    print('x',x.dtype)
    print('th',th.dtype)
    print('Z',Z.dtype)

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
    mask = xvals>=870# FIXME??
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
