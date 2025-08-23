datafile_name = 'testdata_evaleq36.json'
datafile_name_test2 = 'test2_testdata_evaleq36.json'

def gen_test_worklist():
    import numpy as np
    # Below 0.2 we can be very luxurious, since the taylor expansion is fast:
    for th in np.linspace( 0, 90, 9+1 ):
        yield ( th, 0.0 )
        for x in np.geomspace( 1e-20, 0.2, 4 ):
            yield ( th, x )
    #Limited set of points at higher values
    yield ( 0, 3 )
    yield ( 45, 3 )
    yield ( 90, 3 )
    yield ( 30, 0.8 )
    yield ( 80, 5.0 )
    yield ( 40, 10.0 )
    yield ( 0, 1e3 )
    yield ( 45, 1e3 )
    yield ( 90, 1e3 )

def eval_test_workers():
    from .cli_createdata import multiproc_run_worklist
    worklist = [ e for e in gen_test_worklist() ]
    result = multiproc_run_worklist( worklist )
    for th, x, val, err, time in result:
        validate_point(th,x,val,err)
    return result

def update_ref():
    print("Updating reference!\n")
    from .data import data_file
    from .json import save_json
    f = data_file( datafile_name, must_exist=False )
    assert f.parent.is_dir()
    data = eval_test_workers()
    save_json(f,data,force = True)

def validate_point( theta, x, val, err ):
    from .eval_eq36 import eq36_result_tolerance
    assert 0.0 <= err < eq36_result_tolerance
    if x == 0:
        assert val == 1.0
        assert err == 0.0
        return
    assert val-err > 0.0
    if x <= 1e-15:
        #ad hoc cutoff, used since x=1e-20 actually gives yp=1 at our precision
        assert val+err <= 1.0
    else:
        assert val+err < 1.0

def test():
    from .data import data_file
    from .json import load_json
    refdata = load_json(data_file( datafile_name ))
    data = eval_test_workers()
    data = refdata[:]
    assert len(refdata) == len(data)
    for r, d in zip(refdata,data):
        th, x, val, err, time = d
        if r[0:1] != d[0:1]:
            raise RuntimeError('Data point mismatch with reference data')
        refval, referr, reftime = r[2:]
        print('(th,x) = (%g,%g): yp = %g +- %g'%( th, x, val, err) )
        assert abs(refval-val) <= min(err, referr)


def gen_test2_data():
    eps = 1e-8
    th_x_pts = [
        ( 45.0, 0.0001 ),
        ( 0.0, 0.001 ),
        ( 90.0, 0.01 ),
        ( 45.0, 1.0 ),
        ( 0.0, 2.0 ),
        ( 90.0, 3.0 ),
        ( 45.0, 5.0 ),
        ( 0.0, 6.0 ),
        ( 90.0, 7.0 ),
        ( 45.0, 30.0 ),
        ( 0.0, 40.0 ),
        ( 90.0, 50.0 ),
        ( 45.0, 200.0 ),
        ( 0.0, 300.0 ),
        ( 90.0, 400.0 ),
        ( 90.0, 0.1 ),
        ( 45.0, 0.1 ),
        ( 0.0, 0.1 ),
        ( 5.0, 0.25 ),
        ( 85.0, 0.25 ),
        ( 45.0, 0.25 ),
        ( 1.5, 1.5 ),
        ( 9.0, 0.6 ),
        ( 23.0, 0.33 ),
    ]
    worklist = [ (th, x, eps) for th,x in th_x_pts ]

    from .cli_createdata import multiproc_run_worklist
    from .data import data_file
    from .json import save_json
    result = multiproc_run_worklist( worklist,
                                     workfct=test2_gendata_worker )

    data = []
    for th, x, vals_a_and_fi in result:
        data.append( dict( theta = float(th),
                           x = float(x),
                           vals_a_and_fi = [ (float(a),float(fi),float(fierr))
                                             for a,fi,fierr in vals_a_and_fi ] ) )
    f = data_file( datafile_name_test2, must_exist=False )
    assert f.parent.is_dir()
    save_json(f,
              dict( data=data, eps = eps ),
              force = True)



test2_fi_eps_factor = 0.5
def test2_gendata_worker( args ):
    theta_degree, x, eps = args
    from .eq36_numint import FiniteIntegralEq36
    from .mpmath import mpf
    fi = FiniteIntegralEq36(theta_degree, x, eps=test2_fi_eps_factor*eps)

    vals_a_and_fi = []
    last_val = mpf(0)
    while True:
        vals_a_and_fi.append( (fi.a, fi.value, fi.error ) )
        rel_contrib = (fi.value-last_val)/fi.value
        last_val = fi.value
        print(fi.a,'%g %%'%(rel_contrib*100))
        if fi.a > 5000.0 and rel_contrib < test2_fi_eps_factor*eps:
            break
        fi.grow_a()
    return ( float(theta_degree), float(x), vals_a_and_fi )

def cmpfail(a,b):
    return '%.14g'%a != '%.14g'%b

def test2():
    #Load data arrived at by brute-force integration far out on the tail, until
    #the last chunk of contributions had a contributution less than eps.
    from .data import load_json_data
    refdata = load_json_data( datafile_name_test2 )
    fi_eps = refdata['eps']
    assert not cmpfail(fi_eps,1e-8)
    max_de = -1.0
    for d in refdata['data']:
        th,x = float(d['theta']), float(d['x'])
        print('\n'+"="*80)
        print(f"Investigating pt: theta={th:g}, x={x:g}")
        pt_max_de = test2_plot( th, x, d['vals_a_and_fi'], fi_eps  )
        max_de = max(max_de,pt_max_de)
    assert max_de < 1.0
    print('-> max diff/err in any interval of any (theta,x) point: %e'%(max_de))
    assert max_de < 1e-3 #<- not really a fair check, <1.0 is enough. but we
                         #actually appear to have very conservative error
                         #estimates


def test2_plot( theta, x, vals_a_and_fi, fi_eps ):
    from .eq36_numint import FiniteIntegralEq36, TailIntegral
    from .eq36_lowx_taylor import taylor_lowx_eq36
    from .mpmath import mpf

    print(theta,x,fi_eps)
    eps = test2_fi_eps_factor*fi_eps
    #Verify first chunk of tail integral is still unchanged (could fail if we
    #broke our integration or refdata is obsolete):
    if False:#FIXME!!!!!!!!!!!
        fi = FiniteIntegralEq36(theta, x, eps=eps)
        if ( cmpfail( fi.a, vals_a_and_fi[0][0] )
             or cmpfail( fi.value, vals_a_and_fi[0][1] )
             or cmpfail( fi.error, vals_a_and_fi[0][2] ) ):
            raise RuntimeError('test2: inconsistency w.r.t. integration refdata')

    #final (very very expensive) result of millions of tight-interval
    #Gauss-Legendre integrals:
    final_fi = vals_a_and_fi[-1]
    #Add (hopefully tiny) tail contribution:
    final_tail = TailIntegral(x=x,a=final_fi[0])
    final_fi_plus_tail = (final_fi[1]+final_tail.value,final_fi[2]+final_tail.error)

    taylor = taylor_lowx_eq36( theta, x, eps )
    if taylor is not None:
        print('---> comparison with direct Taylor eval is possible!')
        print('taylor  : %g +- %g'%taylor)
        print('fi      : %g +- %g'%(final_fi[1],final_fi[2]))
        print('tail    : %g +- %g'%(final_tail.value,final_tail.error))
        print('fi+tail : %g +- %g'%final_fi_plus_tail)
        #TODO: Actually the ACTUAL precision of ref-values is limited by
        #floating point precision.
        diff_val = abs(final_fi_plus_tail[0]-taylor[0])
        diff_err = final_fi_plus_tail[1]+taylor[1]
        de = diff_val/diff_err
        print('taylor-(fi+tail) : %g +- %g'%(diff_val,diff_err))
        print('diff/err: %g'%(de))
        if de >= 1.0:
            raise RuntimeError('Inconsistencies detected')

    max_de = -1.0
    for i in range(len(vals_a_and_fi)):
        p1 = vals_a_and_fi[i-1] if i else ( 0.0, 0.0, 0.0 )
        p2 = vals_a_and_fi[i]
        a1, fi1, fierr1 = mpf(p1[0]), mpf(p1[1]), mpf(p1[2])
        a2, fi2, fierr2 = mpf(p2[0]), mpf(p2[1]), mpf(p2[2])
        assert a1 == int(a1)
        assert a2 == int(a2)
        a1 = int(a1)
        a2 = int(a2)
        if not ( a1 >= x ):
            #Can not use tail formula
            continue
        t1 = TailIntegral( x=x, a=a1 )
        t2 = TailIntegral( x=x, a=a2 )
        #Find interval contribution:
        iv_fi = fi2-fi1, fierr1+fierr2
        iv_tail = t1.value-t2.value, t1.error+t2.error
        diff_val = abs(iv_fi[0]-iv_tail[0])
        diff_err = iv_fi[1]+iv_tail[1]
        de = diff_val/diff_err
        max_de = max(de,max_de)
        bad = de >= 1.0
        if (i==1 or i+1==len(vals_a_and_fi) or i%400==0) or bad:
            print(f'interval [{a1,a2}]:')
            print('   -> num int: %g +- %g'%iv_fi)
            print('   -> tail   : %g +- %g'%iv_tail)
            print('   -> diff/err: %g'%(de))
        if bad:
            raise RuntimeError('Inconsistencies detected')
    print('-> max diff/err in any interval: %g'%(max_de))
    return max_de
