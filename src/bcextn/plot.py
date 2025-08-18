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
        #print(data[data_key])
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
        plt.plot( xvals, orig_optimalB, label = 'B (optimal at each theta)',
                  color='green')
        plt.plot( xvals, orig_fctB(thvals,*orig_res_fctB),
                  label = 'proposed B(theta)',
                  color='green',linestyle='--',lw=3 )
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
