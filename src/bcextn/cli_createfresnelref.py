from .eval_eq36_scnd_fresneltest import evalref
import pathlib
import sys
import time
from .mpmath import mpf, mpf_pack_to_str, mpf_unpack_from_str
from .json import save_json, load_json

def worker( args ):
    t = time.time()
    x, th, nmin, nmax = args
    result = evalref( theta_degree = th, x=x, nmin=nmin, nmax=nmax,
                      eps=1e-10, maxdegree0 = 4 )
    tot_val = mpf(0)
    tot_err = mpf(0)
    assert len(result) == (nmax+1)-nmin
    for n, val, err in result:
        tot_val += val
        tot_err += err
    t = time.time() - t
    return ( x, th, nmin, nmax, t,
             mpf_pack_to_str(tot_val),
             mpf_pack_to_str(tot_err) )

def multiproc_run_worklist( worklist ):
    import multiprocessing
    import tqdm
    import random
    random.shuffle(worklist)#more reliable progress bar

    with multiprocessing.Pool() as pool:
        results = list(tqdm.tqdm( pool.imap_unordered(worker, worklist,chunksize=1),
                                  total=len(worklist) ))
    tot_val = mpf(0)
    tot_err = mpf(0)
    tot_time = 0.0
    for x, th, nmin, nmax, t, mpfstr_sum_val, mpfstr_sum_err in results:
        tot_val += mpf_unpack_from_str( mpfstr_sum_val )
        tot_err += mpf_unpack_from_str( mpfstr_sum_err )
        tot_time += t
    return tot_time, tot_val, tot_err

def main():
    args = sys.argv[1:]

    if args and args[0] == '--merge':
        return main_merge(args[1:])

    quick = False
    while '--quick' in args:
        quick = True
        args.remove('--quick')
    assert len(args)==4, "please specify: x theta jobidx njobs [--quick]"

    xval = float(args[0])
    xval_str = '%.6g'%xval
    assert xval == float(xval_str)

    thval = float(args[1])
    assert int(thval)==thval
    thval = int(thval)
    thval_str = '%i'%thval
    assert thval == int(thval_str)

    job_idx = int(args[2])
    njobs = int(args[3])

    assert 0 <= job_idx < njobs
    assert 1 <= njobs < 10000
    assert 1e-4 <= xval <= 1000.0
    assert 0 <= thval <= 90
    contribn_step = int(1e3) #contrib n vals per job
    contribn_max = int(1e9)

    if quick:
        contribn_step = int(1e2)
        contribn_max = int(1e5)

    assert njobs+1 < contribn_max

    mergekey = '_mergeablesubjob%iof%i'%(job_idx+1,njobs) if njobs>1 else ''
    outfile =  pathlib.Path('.').joinpath( 'fresnel_contribn_sum_x%s_th%i%s%s.json'
                                           %( xval_str.replace('.','d'),
                                              thval,
                                              ('_quick' if quick else ''),
                                              mergekey
                                             )
                                          )
    print(f"Target file: {outfile}")

    if outfile.exists():
        raise RuntimeError(f"file already exists: {outfile}")
    assert outfile.parent.is_dir()

    worklist = []


    #special work unit for n=0 for simplicity:
    worklist.append( (xval, thval, 0, 0) )
    nmin_used = 0
    nmin = 1
    nmax = contribn_step
    nmax_used = None
    while nmin <= contribn_max:
        worklist.append( (xval, thval, nmin, nmax) )
        nmax_used = nmax
        nmin += contribn_step
        nmax += contribn_step

    if njobs > 1:
        worklist = worklist[ job_idx::njobs ]#start at job_idx, take strides of size njobs

    print("NWORKERS:",len(worklist))
    worklist_print = worklist[:]
    if len(worklist_print)<30:
        for e in worklist_print:
            print('    ',e)
    else:
        for e in worklist_print[:14]:
            print('    ',e)
        print('    ','...(snip)...')
        for e in worklist_print[-14:]:
            print('    ',e)
    tot_time, tot_val, tot_err = multiproc_run_worklist( worklist )

    data = dict( setup = dict( nmin = nmin_used,
                               nmax = nmax_used,
                               xval_str = xval_str,
                               thval_str = thval_str ),
                 results = dict( time = tot_time,
                                 sum_contribn_val = mpf_pack_to_str(tot_val),
                                 sum_contribn_err = mpf_pack_to_str(tot_err) ) )
    if njobs==1:
        save_json( outfile, data )
    else:
        d = dict( data = data,
                  subjobinfo = dict( jobidx = job_idx,
                                     njobs = njobs ) )
        save_json( outfile, d )

def merge( infiles, outfile ):
    jobidxs = []
    d_cfg = None
    assert len(infiles) >= 1
    assert not  outfile.is_file()
    sum_contribn_err = mpf(0)
    sum_contribn_val = mpf(0)
    sum_time = 0.0

    for i,e in enumerate(infiles):
        d = load_json(e)
        assert 'data' in d, f'file not meant for merging: {e}'
        assert 'subjobinfo' in d, f'file not meant for merging: {e}'
        d_data = d['data']
        d_setup = d_data['setup']
        d_results = d_data['results']
        d_subjobinfo = d['subjobinfo']

        d_cfg = ( d_setup['nmin'],
                  d_setup['nmax'],
                  d_setup['xval_str'],
                  d_setup['thval_str'],
                  d_subjobinfo['njobs'] )
        assert d_subjobinfo['jobidx'] not in jobidxs, "same jobidx in multiple files"
        jobidxs.append(d_subjobinfo['jobidx'])
        assert d_subjobinfo['jobidx'] < d_subjobinfo['njobs']
        if i==0:
            cfg = d_cfg
        else:
            assert cfg == d_cfg, "files not at same (x,theta,njobs) point"

        sum_contribn_err += mpf_unpack_from_str(d_results['sum_contribn_err'])
        sum_contribn_val += mpf_unpack_from_str(d_results['sum_contribn_val'])
        sum_time += d_results['time']

    nmin, nmax, xval_str, thval_str, njobs = cfg
    assert len(jobidxs) == njobs, "missing some job files"

    data = dict( setup = dict( nmin = nmin,
                               nmax = nmax,
                               xval_str = xval_str,
                               thval_str = thval_str ),
                 results = dict( time = sum_time,
                                 sum_contribn_val = mpf_pack_to_str(sum_contribn_val),
                                 sum_contribn_err = mpf_pack_to_str(sum_contribn_err) ) )

    save_json( outfile, data )

def main_merge( args ):
    def badusage():
        print("Usage:")
        print("")
        print("$0 --merge JSONFILE0 .. JSONFILEN OUTFILENAME")
        raise SystemExit(1)
    if len(args)<2:
        return badusage()
    infiles = [pathlib.Path(a) for a in args[:-1]]
    outfile = pathlib.Path(args[-1])
    assert all(e.is_file() for e in infiles)
    assert not outfile.is_file()
    assert outfile.parent.is_dir()
    if not infiles:
        return badusage()
    merge(infiles,outfile)
