from .eval_eq36 import BCYPCalc
import numpy as np
import pathlib
import sys
import time

def worker( args ):
    t = time.time()
    theta_deg, xval = args
    print(f"Starting worker at theta={theta_deg:g} x={xval:g}")
    calc = BCYPCalc( theta_deg )
    val, max_err = calc.calc_yp(xval)
    t = time.time() - t

    return ( float(theta_deg),
             float(xval),
             float(val),
             float(max_err),
             float(t) )

def find_output_fn( quick, table1_points, thetascan ):
    assert not (table1_points and thetascan)
    bn = 'bcdata'
    if table1_points:
        bn += '_table1pts'
    if thetascan:
        bn += '_thetascan'
    if quick:
        bn += '_quick'

    for i in range(1,10000000):
        p = pathlib.Path('.').joinpath(
            '%s%s.json'%(bn,'' if i==1 else '_%i'%i)
        )
        if not p.exists():
            assert p.parent.is_dir()
            return p
    assert False

def multiproc_run_worklist( worklist ):
    import multiprocessing
    import tqdm
    import random
    random.shuffle(worklist)#more reliable progress bar

    with multiprocessing.Pool() as pool:
        results = list(tqdm.tqdm( pool.imap(worker, worklist),
                                 total=len(worklist)))
    return sorted(results)

def main():
    quick = '--quick' in sys.argv[1:]
    table1_points = '--table1' in sys.argv[1:]
    thetascan = '--thetascan' in sys.argv[1:]
    assert not (table1_points and thetascan)

    x_range = ( 1e-3, 1e3 )
    nx = 10 if quick else 100
    nth = 3 if quick else 18

    #Hit nicer values:
    nx += 1
    nth += 1

    theta_vals = np.linspace(0.0, 90.0, nth)
    xvals = np.geomspace( *x_range, nx )

    if table1_points:
        from .bc1974_tables import table1
        t = table1()
        xvals = t['xvals']
        theta_vals = np.asin( t['sinthvals'] ) * ( 180 / np.pi )
        theta_vals.sort()
        if quick:
            xmax = xvals[-1]
            xvals = xvals[::7]
            xvals[-1] = xmax
            theta_max = theta_vals[-1]
            theta_vals = theta_vals[::7]
            theta_vals[-1] = theta_max
    if thetascan:
        #Remember tiny xvals are cheap
        xvals = np.asarray( [1e-10, 1e-5, 0.001, 0.01, 0.1, 0.5, 1.0, 1.5,
                             2.0, 3.0, 5.0, 10.0, 30.0, 100 ], dtype='float' )
        if quick:
            xmax = xvals[-1]
            xvals = xvals[::4]
            xvals[-1] = xmax
        theta_vals = np.linspace( 0.0, 90.0,
                                  (9 if quick else 180) +1 )

    def find_outname():
        return find_output_fn( quick, table1_points, thetascan )

    outfile = find_outname()
    print(f"Target file: {outfile}")


    worklist = []
    for th in theta_vals:
        for x in xvals:
            worklist.append( (th, x) )

    results = multiproc_run_worklist( worklist )

    r_vals = {}
    r_err = {}
    r_time = {}
    for th, x, yp, yp_maxerr, t in results:
        if th not in r_vals:
            r_vals[th] = []
            r_err[th] = []
            r_time[th] = []
        r_vals[th].append( float(yp) )
        r_err[th].append( float(yp_maxerr) )
        r_time[th].append( float(t) )

    final = { 'xvals' : [float(x) for x in xvals],
              'theta_2_ypvals' : r_vals,
              'theta_2_ypvals_maxerr' : r_err,
              'theta_2_calctime' : r_time,
             }

    if outfile.exists():
        print("WARNING: {outfile} found to exist now!")
        outfile = find_outname()
        print("WARNING: Writing to {outfile} instead!")

    from .json import save_json
    save_json( outfile,final )

