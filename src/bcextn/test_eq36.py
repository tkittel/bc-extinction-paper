datafile_name = 'testdata_evaleq36.json'

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
