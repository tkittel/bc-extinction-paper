from .eval_eq36 import eval_eq36
from .eval_eq36_scnd import ( eval_eq36_scnd_gauss,
                              eval_eq36_scnd_lorentz,
                              eval_eq36_scnd_fresnel )
import numpy as np
import pathlib
import sys
import time

def worker( args ):
    t = time.time()
    fcttype,theta_deg, xval = args
    assert fcttype in ('primary','scndgauss','scndlorentz','scndfresnel')
    if fcttype=='primary':
        fct = eval_eq36
    elif fcttype=='scndlorentz':
        fct = eval_eq36_scnd_lorentz
    elif fcttype=='scndfresnel':
        fct = eval_eq36_scnd_fresnel
    else:
        assert fcttype == 'scndgauss'
        fct = eval_eq36_scnd_gauss
    #print(f"Starting worker at theta={theta_deg:g} x={xval:g}")
    val, max_err = fct( theta_deg, xval )
    t = time.time() - t
    #print("WORKER DONE!!")
    return ( float(theta_deg),
             float(xval),
             float(val),
             float(max_err),
             float(t) )

def worker2( args ):
    fcttype,theta_deg, xval = args
    res = worker( args )
    return fcttype, res

def find_output_fn( quick, table1_points, thetascan, xscan090, xscannear45,
                    is_scnd_gauss, is_scnd_lorentz, is_scnd_fresnel,
                    is_user_test_data, is_highx ):
    assert sum(int(bool(e)) for e in [is_scnd_gauss,is_scnd_lorentz,is_scnd_fresnel]) in (0,1)
    assert not (table1_points and thetascan)
    bn = 'bcdata'
    if is_scnd_gauss:
        bn += '_scndgauss'
    if is_scnd_lorentz:
        bn += '_scndlorentz'
    if is_scnd_fresnel:
        bn += '_scndfresnel'
    if table1_points:
        bn += '_table1pts'
    if thetascan:
        bn += '_thetascan'
    if xscan090:
        bn += '_xscan090'
    if xscannear45:
        bn += '_xscannear45'
    if is_user_test_data:
        #Override
        bn = 'bc2025_reference_x_sintheta_yp'
    if is_highx:
        #Override
        bn = 'bcdata_special_highx'
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

def multiproc_run_worklist( worklist, workfct = worker ):
    import multiprocessing
    import tqdm
    import random
    random.shuffle(worklist)#more reliable progress bar

    with multiprocessing.Pool() as pool:
        results = list(tqdm.tqdm( pool.imap_unordered(workfct, worklist,chunksize=1),
                                  total=len(worklist) ))
    return sorted(results)

def main():
    args = sys.argv[1:]
    if '--all' not in args:
        doit( args )
        return
    is_quick = '--quick' in args
    import multiprocessing
    subjobs = []
    for m in [None,'--scnd-gauss','--scnd-lorentz','--scnd-fresnel']:
        for t in [None,'--table1','--thetascan','--xscan090','--xscannear45']:
            j = []
            if m is not None:
                j.append(m)
            if t is not None:
                j.append(t)
            if is_quick:
                j.append('--quick')
            subjobs.append(j)
    j = ['--usertestdata']
    if is_quick:
        j.append('--quick')
    subjobs.append(j)

    j = ['--highx']
    if is_quick:
        j.append('--quick')
    subjobs.append(j)

    for i,jobargs in enumerate(subjobs):
        print()
        print( '-'*80 )
        print( '-'*80 )
        print( '-'*80 )
        print( '------ LAUNCHING SUBJOB %i of %i: args=%s'%(i+1,
                                                            len(subjobs),
                                                            jobargs) )
        print( '-'*80 )
        print( '-'*80 )
        print( '-'*80 )
        print()

        #Run in separate processes for safety (but not concurrently, since there
        #is already concurrency at a lower level):
        p = multiprocessing.Process(target=doit, args=(jobargs,))
        p.start()
        p.join()
        #doit( jobargs )

def doit( args ):
    is_user_test_data = '--usertestdata' in args
    is_highx = '--highx' in args

    quick = '--quick' in args
    table1_points = '--table1' in args
    thetascan = '--thetascan' in args
    xscan090 = '--xscan090' in args
    xscannear45 = '--xscannear45' in args
    xscan = not( xscannear45 or xscan090 or thetascan or table1_points or is_user_test_data or is_highx )
    is_scnd_gauss = '--scnd-gauss' in args
    is_scnd_lorentz = '--scnd-lorentz' in args
    is_scnd_fresnel = '--scnd-fresnel' in args

    assert sum(int(bool(e)) for e in (xscan,table1_points,thetascan,xscan090,xscannear45,is_user_test_data,is_highx))==1
    assert not (is_scnd_gauss and is_scnd_lorentz)
    assert not (is_scnd_gauss and is_scnd_fresnel)
    assert not (is_scnd_fresnel and is_scnd_lorentz)
    fcttype = 'primary'
    if is_scnd_gauss:
        fcttype = 'scndgauss'
    if is_scnd_lorentz:
        fcttype = 'scndlorentz'
    if is_scnd_fresnel:
        fcttype = 'scndfresnel'

    if is_user_test_data and fcttype!='primary':
        raise SystemExit('Do not put a mode keyword with --usertestdata')

    if is_highx and fcttype!='primary':
        raise SystemExit('Do not put a mode keyword with --highx')

    x_range = ( 1e-3, 1e3 )
    nx = 3 if quick else 100
    nth = 3 if quick else 18

    #Hit nicer values:
    nx += 1
    nth += 1

    theta_vals = np.linspace(0.0, 90.0, nth)
    xvals = np.geomspace( *x_range, nx )

    if table1_points:
        from .bc1974_tables import table1
        t = table1()
        xvals = [e for e in t['xvals']]
        xvals += [ 0.001, 0.01, 100., 1000. ] #Extra values for new table
        xvals = np.asarray(xvals,dtype=float)
        xvals.sort()
        sinth_vals = [ e for e in t['sinthvals'] ]
        sinth_vals += [ 0.0, 0.99, 1.0 ] #Extra values for new table
        sinth_vals = np.asarray(sinth_vals, dtype=float)
        theta_vals = np.asin( sinth_vals ) * ( 180 / np.pi )
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
    if xscan090:
        #For fitting at theta=0 or theta=90
        x_range = ( 1e-3, 1e3 )
        nx = 10 if quick else 1000
        xvals = np.geomspace( *x_range, nx )
        theta_vals = np.asarray([0.0, 90.0],dtype=float)

    if xscannear45:
        #For providing xscan points also at (45-eps,45+eps), since the classic
        #scndlorentz recipe has a breakdown at 45<th<45+eps.
        theta_vals = np.asarray([44.99, 45.01],dtype=float)

    def find_outname():
        return find_output_fn( quick, table1_points, thetascan, xscan090, xscannear45,
                               is_scnd_gauss, is_scnd_lorentz, is_scnd_fresnel,
                               is_user_test_data = is_user_test_data,
                               is_highx = is_highx )

    outfile = find_outname()
    print(f"Target file: {outfile}")


    worklist = []
    if is_highx:
        del theta_vals
        del xvals
        xvals = [ 999.0, 1000.0 ]
        theta_vals = [ 0, 90 ]
        if quick:
            #No change!
            pass
        fcttypes = ['primary','scndgauss','scndlorentz','scndfresnel']

        worklist = []

        for fcttype in fcttypes:
            for thdeg in theta_vals:
                for x in xvals:
                    worklist.append( (fcttype, thdeg, x) )
        results = multiproc_run_worklist( worklist, workfct = worker2 )
        data = {}
        for fcttype, (th, x, yp, yp_maxerr, t) in results:
            assert int(round(th)) == th
            th = int(round(th))
            assert th in [0,90]
            if fcttype not in data:
                data[fcttype] = []
            assert yp_maxerr < 1e-9
            assert yp+yp_maxerr <= 1.0
            assert yp-yp_maxerr >= 0.0
            data[fcttype].append( ( th, float(x), float(yp) ) )

        assert set(data.keys()) == set(fcttypes)
        for f in fcttypes:
            data[f] = sorted(data[f][:])
        from .json import save_json
        save_json( outfile,data )
        raise SystemExit

    if is_user_test_data:
        del theta_vals
        del xvals
        xvals_extrapolation_threshold = 1000
        xvals = [ 1e-20, 1e-13,1e-12,1e-4,1e-3,
                  0.009,0.01,0.011,#around lux taylor threshold
                  0.09,0.1,0.11,#around std taylor threshold
                  0.5, 1.5, 4.0, 10.0, 30.0,#normal values
                  999,xvals_extrapolation_threshold
                 ]
        xvals_extrapolated = [ 1000.1, 1e4, 1e12, 1e99, 1e200 ]
        sinthvals = [ 0.0, 0.37, 0.7071, 0.93, 1.0 ]
        if quick:
            #xvals = xvals[::3]
            #sinthvals = sinthvals[::3]
            xvals = [1.0]
        fcttypes = ['primary','scndgauss','scndlorentz','scndfresnel']

        worklist = []

        for fcttype in fcttypes:
            for sinth in sinthvals:
                for x in xvals:
                    thdeg = float(np.asin(sinth)*180.0/np.pi)
                    if sinth==0.0:
                        thdeg = 0
                    if sinth==1.0:
                        thdeg = 90
                    worklist.append( (fcttype, thdeg, x) )
        results = multiproc_run_worklist( worklist, workfct = worker2 )
        data = {}
        for fcttype, (th, x, yp, yp_maxerr, t) in results:
            if fcttype not in data:
                data[fcttype] = []
            assert yp_maxerr < 1e-9
            assert yp+yp_maxerr <= 1.0
            assert yp-yp_maxerr >= 0.0
            sinth_approx = float(np.sin(float(th)*np.pi/180.))
            sinth = [ v for v in sinthvals if abs(v-sinth_approx)<1e-9 ]
            assert len(sinth)==1
            sinth = sinth[0]
            data[fcttype].append( ( float(x), sinth, float(yp) ) )
            if abs(x-xvals_extrapolation_threshold)<1e-9:
                #Fire off extrapolation values:
                from .new_recipes import recipe_highx_pow
                epower = recipe_highx_pow(fcttype)
                for xe in xvals_extrapolated:
                    yp_e = yp*( xe/xvals_extrapolation_threshold) ** (-epower)
                    data[fcttype].append( ( float(xe), sinth, float(yp_e) ) )

        assert set(data.keys()) == set(fcttypes)
        for f in fcttypes:
            data[f] = sorted(data[f][:])
        from .json import save_json
        save_json( outfile,data )
        raise SystemExit


    for th in theta_vals:
        for x in xvals:
            worklist.append( (fcttype, th, x) )

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

