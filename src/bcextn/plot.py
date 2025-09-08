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

def std_plots( data, mode ):
    if mode == 'all':
        for m in std_modemap.keys():
            std_plots(data,m)
        return
    data_key = std_modemap[mode]

    for thstr in data['th_keys']:
        thval = float(thstr)
        color = th2color(thstr)
        plt.plot( data['xvals'],
                  data[data_key][thstr],
                  label = f'$\\theta={thval:g}\\degree$',
                  color = color )
    plt.ylabel(mode)
    plt.xlabel('x')
    plt.legend()
    plt.grid()
    plt.show()

def main_xscan( args ):
    from .load import load_xscan
    assert len(args)==1 and (args[0]=='all' or args[0] in std_modemap.keys())
    std_plots( load_xscan(), args[0] )

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
        #df = _tableN_as_dataframe( t )
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

def main_fit( args ):
    from .load import load_thetascan
    from .curves import yp_proposed_curve as f
    from .curves import yp_bc1974_curve as orig_f
    assert len(args)==0
    data = load_thetascan()
    xvals = data['xvals']

    orig_optimalA = []
    orig_optimalB = []
    optimalA = []
    optimalB = []
    optimalC = []
    thvals = []

    for theta_degree_str, ypvals in data['theta_2_ypvals'].items():
        thvals.append(float(theta_degree_str))
        (A,B,C),_ = scipy.optimize.curve_fit(f,xvals, ypvals,
                                             p0 = [1.0,1.0,1.0] )
        optimalA.append(A)
        optimalB.append(B)
        optimalC.append(C)
        (A,B),_ = scipy.optimize.curve_fit(orig_f,xvals, ypvals,
                                           p0 = [1.0,1.0] )
        orig_optimalA.append(A)
        orig_optimalB.append(B)

    thvals = array(thvals)
    optimalA = array(optimalA)
    optimalB = array(optimalB)
    optimalC = array(optimalC)
    orig_optimalA = array(orig_optimalA)
    orig_optimalB = array(orig_optimalB)

    #Let us fit and visualise proposed and original forms:
    fctA_p0 = [ 1.0, ] * 7
    @np.vectorize
    def fctA( theta, p0, p1, p2,p3,p4,p5,p6 ):
        u = math.sin( theta*kDeg2Rad )
        return p0 + p1*u + p2*u**2 + p3*u**3 +p4*u**4 +p5*u**5 + p6*u**6
    res_fctA,_ = scipy.optimize.curve_fit( fctA, thvals, optimalA, p0 = fctA_p0 )
    print('fctA(theta) pars:',res_fctA)

    fctB_p0 = [ 1.0, ]*7
    @np.vectorize
    def fctB( theta, p0, p1, p2, p3, p4, p5, p6 ):
        u = math.sin( theta*kDeg2Rad )
        return p0 + p1*u + p2*u**2 + p3*u**3 + p4*u**4 + p5*u**5 + p6*u**6
    res_fctB,_ = scipy.optimize.curve_fit( fctB, thvals, optimalB, p0 = fctB_p0 )
    print('fctB(theta) pars:',res_fctB)

    fctC_p0 = [ 1.0, ] * 7
    @np.vectorize
    def fctC( theta, p0, p1, p2, p3, p4, p5, p6 ):
        u = math.sin( theta*kDeg2Rad )
        return p0 + p1*u + p2*u**2 + p3*u**3 + p4*u**4 + p5*u**5 + p6*u**6
    res_fctC,_ = scipy.optimize.curve_fit( fctC, thvals, optimalC, p0 = fctC_p0 )
    print('fctC(theta) pars:',res_fctC)

    orig_fctA_p0 = [ 1.0, ] * 2
    @np.vectorize
    def orig_fctA( theta, p0, p1 ):
        u = math.cos( 2 * theta*kDeg2Rad )
        return p0 + p1*u
    orig_res_fctA,_ = scipy.optimize.curve_fit( orig_fctA, thvals,
                                                orig_optimalA, p0 = orig_fctA_p0 )
    print('orig fctA(theta) pars:',orig_res_fctA)

    orig_fctB_p0 = [ 1.0, ] * 2
    @np.vectorize
    def orig_fctB( theta, p0, p1 ):
        u = math.cos( 2 * theta*kDeg2Rad )
        return p0 + p1*(0.5-u)**2
    orig_res_fctB,_ = scipy.optimize.curve_fit( orig_fctB, thvals,
                                                orig_optimalB, p0 = orig_fctB_p0 )
    print('orig fctB(theta) pars:',orig_res_fctB)


    plots = [ ( thvals, 'theta [degree]' ),
              ( np.cos(thvals*(2*np.pi/180.)), 'cos(2*theta)' ),
              ( np.sin(thvals*(np.pi/180.)), 'sin(theta)' ),
             ]
    for xvals, xlabel in plots:
        plt.plot( xvals, optimalA, label = 'A (optimal at each theta)',
                  color='blue')
        plt.plot( xvals, fctA(thvals,*res_fctA),
                  label = 'proposed A(theta)',
                  color='blue',linestyle='--',lw=3)
        plt.plot( xvals, optimalB, label = 'B (optimal at each theta)',
                  color='green')
        plt.plot( xvals, fctB(thvals,*res_fctB),
                  label = 'proposed B(theta)',
                  color='green',linestyle='--',lw=3)
        plt.plot( xvals, optimalC, label = 'C (optimal at each theta)',
                  color='red')
        plt.plot( xvals, fctC(thvals,*res_fctC),
                  label = 'proposed C(theta)',
                  color='red',linestyle='--',lw=3)
        plt.grid()
        plt.xlabel(xlabel)
        plt.legend()
        plt.show()
    for xvals, xlabel in plots:
        plt.plot( xvals, orig_optimalA, label = 'A (optimal at each theta)',
                  color='blue')
        plt.plot( xvals, orig_fctA(thvals,*orig_res_fctA),
                  label = 'proposed A(theta)',
                  color='blue',linestyle='--',lw=3 )

        plt.plot( xvals, orig_fctA(thvals,0.20,0.45),
                  label = 'original 1974 A(theta)',
                  color='blue',linestyle=':',lw=3, alpha=0.5 )
        plt.plot( xvals, orig_optimalB, label = 'B (optimal at each theta)',
                  color='green')
        plt.plot( xvals, orig_fctB(thvals,*orig_res_fctB),
                  label = 'proposed B(theta)',
                  color='green',linestyle='--',lw=3 )
        plt.plot( xvals, orig_fctB(thvals,0.22,-0.12),
                  label = 'original 1974 B(theta)',
                  color='green',linestyle=':',lw=3, alpha=0.5 )

        plt.grid()
        plt.xlabel(xlabel)
        plt.legend()
        plt.show()

    results = { 'A' : [ p for p in res_fctA ],
                'B' : [ p for p in res_fctB ],
                'C' : [ p for p in res_fctC ],
                'origA' : [ p for p in orig_res_fctA ],
                'origB' : [ p for p in orig_res_fctB ]
               }

    from .json import save_json
    save_json(
        'fitted_curves_ABC.json',
        results,
        force = True
    )

def lookup_float_key( adict, key, eps=1e-12 ):
    vals = [ v for k,v in adict.items()
             if abs(float(key)-float(k))<eps ]
    assert len(vals)==1
    return vals[0]

def main_cmprecipes( args ):
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
                   do_vs_old = do_vs_old )

def do_cmprecipes( do_reldiff, do_vs_old ):
    from .curves import ClassicCurve, UpdatedClassicCurve, ProposedCurve
    from .load import load_xscan as dataload

    if do_vs_old:
        assert do_reldiff, "vsold only for reldiff"
    data = dataload()
    print( data.keys() )
    xvals = data['xvals']

    def fleq( a, b, eps=1e-12 ):
        return abs(float(a)-float(b))<eps

    th2yp = dict( (th,ypvals)
                  for th,ypvals in data['theta_2_ypvals'].items()
                  #if ( int(float(th)) == float(th) and int(float(th))==80 ) #int(float(th))%30==0
                  #if ( int(float(th)) == float(th) and int(float(th))%5==0 )
                  )
    if do_reldiff:
        def ytrf( yp, ypref ):
            return np.clip(abs(yp/ypref-1.0),0.0,1.0)
    else:
        def ytrf( yp, ypref ):
            return np.clip(yp,0.0,1.0)

    seen_lbls = set()
    for th,yref in th2yp.items():
        th = float(th)
        #if th > 61:
        #    continue
        color = th2color(th)
        yp_classic = ClassicCurve()(xvals,th)
        yp_updatedclassic = UpdatedClassicCurve()(xvals,th)
        yp_proposed = ProposedCurve()(xvals,th)
        if do_vs_old:
            yref = yp_classic

        common = dict( alpha = 0.5,
                       lw = 2 )

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
            plt.plot( xvals, ytrf(yp_classic,yref),
                      **common,
                      color=color,
                      ls = '-.',
                      label = lbltrf('BC1974',split_th) )
        #color='green'
        color = 'red'#th2color(th*0.3,'Greens')
        if do_vs_old and th > split_th:
            color='green'
        plt.plot( xvals, ytrf(yp_updatedclassic,yref),
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
        plt.plot( xvals, ytrf(yp_proposed,yref),
                  **common,
                  color=color,
                  ls = '-',
                  label = lbltrf('new form',
                                 split_theta = split_th if do_vs_old else None) )
    if do_reldiff:
        plt.semilogy()
        plt.ylim(1e-6,1.0)
        plt.ylabel('yp/ypref-1')
    else:
        plt.ylim(0.0,1.0)
        plt.ylabel('yp')
    plt.xlim(xvals[0],xvals[-1])
    plt.grid()
    plt.semilogx()
    plt.legend()
    plt.xlabel('x')
    plt.show()

def main_fittest( args ):
    #Try to fit new candidates.
    from .load import load_xscan as dataload
    from .curves import safediv, safesqrt
    from .fit import fitfunction
    assert len(args)==0

    fcts = []
    @fitfunction
    def f0( x, A, B ):
        k = 1 + B*x
        k2 = 1.0+2.0*x+safediv(A*x*x,k)
        return safediv( 1.0, safesqrt( k2 ) )

    fcts = []#fixme
    @fitfunction
    def f1( x, A, B, C ):
        k = 1 + B*x
        k2 = 1.0+2.0*x+safediv(A*x*x+C*x,k)
        return safediv( 1.0, safesqrt( k2 ) )

    @fitfunction
    def f2( x, A, B, C, D ):
        b = D*(x/(1.0+x))
        k = 1.0 + B*x +b
        k2 = 1.0+2.0*x+safediv(A*x*x+C*x,k)
        return safediv( 1.0, safesqrt( k2 ) )

    @fitfunction
    def f3( x, A, B, C ):
        k = 1.0 + B*x + C*x/(1.0+x)
        k2 = 1.0+2.0*x+safediv( A*x*x-0.1*x, k )
        return safediv( 1.0, safesqrt( k2 ) )
    fcts.append((f0,'orig','green'))
    fcts.append( (f1,f1.name,'blue') )
    fcts.append( (f2,f2.name,'red') )
    fcts.append( (f3,f3.name,'orange') )

    data = dataload()
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
    data = load_fitted_curve_parameters()

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
    assert not args
    from .curves import ProposedCurve
    from .load import load_xscan as dataload
    data = dataload()
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
              reldiff(ProposedCurve()(xvals,90),ProposedCurve()(xvals,0)),
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

def plot_breakdown( curve, ref ):
    from .curves import ProposedCurve, ClassicCurve, UpdatedClassicCurve
    from .load import load_xscan as dataload

    assert curve in ('classic','updatedclassic','proposed')
    assert ref in ('numint','proposed')
    assert curve != ref

    if ref=='numint':
        refdata = dataload()
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
        x = np.geomspace(0.1, 100.0, 400)  # Adjust the range and resolution as needed
        th = np.linspace(0.1, 100.0, 400)
        ref_curve = ProposedCurve()
        def lookup_ref(xval, thval):
            return ref_curve(xval,thval)

    X, TH = np.meshgrid(x, th)
    if curve=='classic':
        curvefct = ClassicCurve()
    elif curve=='updatedclassic':
        curvefct = UpdatedClassicCurve()
    else:
        assert curve=='proposed'
        curvefct = ProposedCurve()

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
    assert not args
    plot_breakdown( 'classic', ref='proposed')
    plot_breakdown( 'classic', ref='numint')
    plot_breakdown( 'updatedclassic', ref='numint')
    plot_breakdown( 'updatedclassic', ref='proposed')
    plot_breakdown( 'proposed', ref='numint')


def main_printcpp( args ):
    assert not args
    from .curves import load_fitted_curve_parameters
    d = load_fitted_curve_parameters()
    print(d)
    npars = 7
    for v in 'ABC':
        assert len(d[v]) == npars
        for i,pi in enumerate(d[v]):
            print(f'constexpr double c{v}{i} = %g;'%pi)
    print('const double s = sintheta;')
    for v in 'ABC':
        c = f'c{v}'
        s = f'{c}0'
        for i in range(1,npars):
            s+=f'+s*({c}{i}'
        s+=')'*(npars-1)
        print(f'const double {v} = {s};')
    print()
    for v in 'ABC':
        pars = d[v]
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
