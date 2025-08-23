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


def main_table1( args ):
    from .load import load_table1scan
    from .bc1974_tables import table1 as orig_table
    assert len(args)==0 or (len(args)==1 and args[0]=='all')
    do_all = ( args and args[0] == 'all' )

    data = load_table1scan()

    if do_all:
        std_plots(data,'all')

    def as_table( df ):
        return df.pivot(index='x', columns='sinth', values='value')

    def print_table( df ):
        print( as_table(df) )


    orig_table = orig_table()
    df_orig = pd.DataFrame(orig_table['x_sinth_yp_list'],
                           columns=['x', 'sinth', 'value'])
    def new_as_dataframe():
        flat = []
        xvals = data['xvals']
        for th, ypvals in data['theta_2_ypvals'].items():
            sinth = float('%g'%(np.sin(kDeg2Rad*float(th))))
            for x, yp in zip( xvals, ypvals ):
                flat.append( ( x, sinth, yp ) )
        return pd.DataFrame(flat, columns=['x', 'sinth', 'value'])

    new_df = new_as_dataframe()

    print_table( df_orig )
    print_table( new_df )

    from .print_table1_diff import format_table1_heatmap as ft
    ft( as_table(df_orig), as_table(new_df), 'table1_diff.html',
        do_print = True )
    ft( as_table(df_orig), as_table(new_df), 'table1_diff.tex' )

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
