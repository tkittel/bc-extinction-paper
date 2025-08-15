from .eval_eq36 import BCYPCalc
import numpy
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

def find_output_fn( quick ):
    for i in range(1,10000000):
        p = pathlib.Path('.').joinpath(
            '%s%s.json'%('bcdata_quick' if quick else 'bcdata',
                         '' if i==1 else '_%i'%i)
        )
        if not p.exists():
            assert p.parent.is_dir()
            return p

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

    x_range = ( 1e-3, 1e3 )
    nx = 10 if quick else 100
    nth = 3 if quick else 18

    #Hit nicer values:
    nx += 1
    nth += 1

    theta_vals = numpy.linspace(0.0, 90.0, nth)
    xvals = numpy.geomspace( *x_range, nx )

    outfile = find_output_fn( quick )
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
        outfile = find_output_fn( quick )
        print("WARNING: Writing to {outfile} instead!")

    from .json import save_json
    save_json( outfile,final )

