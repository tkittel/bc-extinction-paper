
from .json import load_json
from .mpmath import mp, mpf, mpf_unpack_from_str
from .curves import mode_curves

def tail_limit_estimate(nmax):
    #if we summed contribn from 1 to nmax, that covered the integral over
    #[0,nmax*4pi/3]. The upper limit of the integral over [nmax*4pi/3,inf] is
    #2/(pi^2*(nmax+1)):
    k = mpf(2) / mp.pi**2
    return k / ( mpf(nmax) + 1 )

def _load_and_embellish(datafile):
    d = load_json(datafile)

    assert len(d)==2
    setup = d['setup']
    results = d['results']
    def fixval(v):
        if isinstance(v,str) and v.startswith('mpf('):
            return mpf_unpack_from_str(v)
        return v
    results = dict( (k,fixval(v)) for k,v in results.items() )
    tail_lim = tail_limit_estimate( setup['nmax'] )

    results['tail_contrib_upper_limit'] = tail_lim
    val = results['sum_contribn_val'] + tail_lim/2
    err = results['sum_contribn_err'] + tail_lim/2
    results['val'] = val
    results['err'] = err

    _,_,_,curve_luxrecipe = mode_curves( 'scndfresnel',
                                         use_actual_proposed_pyrecipe = True )
    setup['xval'] = float(setup['xval_str'])
    setup['thval'] = float(setup['thval_str'])

    lux_recipe_val = float(curve_luxrecipe()( x = setup['xval'],
                                              theta = setup['thval'] ))
    results['lux_recipe_val'] = lux_recipe_val
    return dict( results = results,
                 setup = setup )

_cache = [None]
def load_refpts():
    res = _cache[0]
    if res is None:
        from .data import datadir
        ddir = datadir()
        res = {}
        for f in sorted(ddir.glob('extnfresnelref_x*_th*.json')):
            d = _load_and_embellish(f)
            key = ( d['setup']['xval_str'], d['setup']['thval_str'] )
            assert key not in res
            res[key] = d
        res = dict(sorted(res.items()))
        _cache[0] = res
    import copy
    return copy.deepcopy( res )

_cache_sd = [None]
def load_stepdata():
    res = _cache_sd[0]
    if res is None:
        from .data import load_json_data
        _raw = load_json_data('fresnelref_stepdata.json')
        res={}
        for d in _raw:
            key = (d['x_str'],d['th_str'])
            pts = []
            for n, mpfstr_contribn, mpfstr_contribn_err in d['pts']:
                pts.append( ( n,
                              mpf_unpack_from_str(mpfstr_contribn),
                              mpf_unpack_from_str(mpfstr_contribn_err) ) )
            assert key not in res
            res[key] = dict( setup = dict( x_str = d['x_str'],
                                           th_str = d['th_str'],
                                           x = float(d['x_str']),
                                           th = float(d['th_str']) ),
                             pts = sorted(pts) )
        res = dict(sorted(res.items()))
        _cache_sd[0] = res
    import copy
    return copy.deepcopy( res )
