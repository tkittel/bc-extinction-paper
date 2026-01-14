
# Load data previously created with the createdata script

def load_xscan( mode = 'primary', split45 = False ):
    if mode.startswith('curve::'):
        return load_curves( mode[len('curve::'):], 'xscan' )
    assert mode in ['primary','scndfresnel','scndlorentz','scndgauss']
    if mode == 'primary':
        fn ='bcdata.json'
    else:
        fn = f'bcdata_{mode}.json'
    fn_split45 = None
    if split45:
        fn_split45 = fn.replace('.json','_xscannear45.json')
    return _loadcache(fn, fn_split45=fn_split45)

def load_xscan090( mode = 'primary' ):
    if mode.startswith('curve::'):
        return load_curves( mode[len('curve::'):], 'xscan090' )
    assert mode in ['primary','scndfresnel','scndlorentz','scndgauss']
    if mode == 'primary':
        return _loadcache('bcdata_xscan090.json')
    else:
        return _loadcache(f'bcdata_{mode}_xscan090.json')

def load_thetascan( mode = 'primary' ):
    if mode.startswith('curve::'):
        return load_curves( mode[len('curve::'):], 'thetascan' )
    assert mode in ['primary','scndfresnel','scndlorentz','scndgauss']
    if mode == 'primary':
        return _loadcache('bcdata_thetascan.json')
    else:
        return _loadcache(f'bcdata_{mode}_thetascan.json')

def _impl_load_tablescan_origpts( tablename ):
    assert tablename!='tableF', "not fresnel table in BC1974 paper"
    assert tablename in ('table1','table3','table4')
    cachekey = f'load_{tablename}scan_origpts_cache'
    if cachekey in _cache:
        return _cache[cachekey]

    from .bc1974_tables import table1, table3, table4
    from .plotutils import array
    import numpy as np
    import copy

    if tablename == 'table1':
        orig_table_data = table1()
        d_all = load_table1scan_allpts()
    elif tablename == 'table3':
        orig_table_data = table3()
        d_all = load_table3scan_allpts()
    elif tablename == 'table4':
        orig_table_data = table4()
        d_all = load_table4scan_allpts()
    else:
        assert False
    d_all = copy.deepcopy(d_all)
    orig_xvals = orig_table_data['xvals']
    orig_th_vals = [ np.asin(float(e))*180/np.pi
                     for e in orig_table_data['sinthvals'] ]
    def x_ok( xval ):
        return any(abs(float(xval)-float(e))<1e-13 for e in orig_xvals)
    def th_ok( thval ):
        return any(abs(float(thval)-float(e))<1e-13 for e in orig_th_vals)
    mask = np.array( [x_ok(e) for e in d_all['xvals']],dtype = bool )
    d_filtered = {}
    d_filtered['xvals'] = array([ e for e in d_all['xvals'] if x_ok(e) ])
    for k in ['theta_2_calctime', 'theta_2_ypvals', 'theta_2_ypvals_maxerr']:
        d_filtered[k] = dict( (k,v[mask]) for k,v in d_all[k].items() if th_ok(k) )
    d_filtered = _prepdata( d_filtered )
    assert set(d_all.keys()) == set(d_filtered.keys())
    _cache[cachekey] = d_filtered
    return d_filtered

def load_table1scan_origpts():
    return _impl_load_tablescan_origpts('table1')

def load_table3scan_origpts():
    return _impl_load_tablescan_origpts('table3')

def load_table4scan_origpts():
    return _impl_load_tablescan_origpts('table4')

def load_table1scan_allpts():
    return _loadcache('bcdata_table1pts.json')

def load_table3scan_allpts():
    return _loadcache('bcdata_scndgauss_table1pts.json')

def load_table4scan_allpts():
    return _loadcache('bcdata_scndlorentz_table1pts.json')

def load_tableFscan_allpts():
    return _loadcache('bcdata_scndfresnel_table1pts.json')

_cache={}
def _loadcache( fn, fn_split45 = None ):
    from .data import load_json_data as _load
    key = (fn,fn_split45)
    d = _cache.get(key)
    if not d:
        rawdata = _load(fn)
        if fn_split45 is not None:
            rawdata45 = _load(fn_split45)
            rawdata = _apply_split45( rawdata, rawdata45 )
        d = _prepdata(rawdata)
        _cache[key] = d
    return d

def _prepdata( data ):
    from .plotutils import array
    d = data
    ll = sorted( (float(e),e) for e in d['theta_2_ypvals'].keys() )
    d[ 'th_keys'] = [ thstr for th,thstr in ll ]
    d['xvals'] = array(d['xvals'])
    for k in d.keys():
        if k.startswith('theta_2_'):
            for thstr,v in d[k].items():
                d[k][thstr] = array( d[k][thstr] )
    return d


_cache_curvedata = {}
def load_curves( curvekey, mode ):
    assert mode in ('xscan','thetascan', 'xscan090')
    if curvekey in _cache_curvedata:
        return _cache_curvedata[curvekey]
    from . import curves
    if curvekey == 'sabineprimary':
        fct = curves.SabineCurve_Primary()
    elif curvekey == 'sabinescndrec':
        fct = curves.SabineCurve_ScndRec()
    elif curvekey == 'sabinescndtri':
        fct = curves.SabineCurve_ScndTri()
    else:
        assert False, "bad curvekey"
    loadfct = dict(xscan=load_xscan,
                   xscan090=load_xscan090,
                   thetascan=load_thetascan,)[mode]
    refdata = loadfct('primary')
    res = _prepdata( _load_curves_impl(fct,refdata) )
    _cache_curvedata[curvekey] = res
    return res

def _load_curves_impl( curvefct, refdata ):
    xvals = refdata['xvals'].copy()
    th2y, th2ye, th2t = {}, {}, {}
    for thkey, ypvals in refdata['theta_2_ypvals'].items():
        th = float(thkey)
        th2y[thkey] = [ curvefct(x,th) for x in xvals ]
        th2ye[thkey] = [ 0.0, ] * len(xvals)
        th2t[thkey] = [ 0.0, ] * len(xvals)
    return dict( xvals = xvals,
                 theta_2_ypvals = th2y,
                 theta_2_ypvals_maxerr = th2ye,
                 theta_2_calctime = th2t )

def load_legendre_lux( mode ):
    return load_legendre( mode, is_lux = True )

_cache_leg = {}
def load_legendre( mode, is_lux = False ):
    fn1 = 'legendre_coefficients%s.json'%('_lux' if is_lux else '')
    fn2 = 'legendre_coefficients_ydelta%s.json'%('_lux' if is_lux else '')
    key = ('leg',fn1, fn2, mode, is_lux )
    if key in _cache_leg:
        return _cache_leg[key]

    #Open files:
    key_filedata = ('filedata',fn1,fn2)
    data = _cache_leg.get( key_filedata )
    if not data:
        from .data import load_json_data
        data1 = load_json_data(fn1)
        data2 = load_json_data(fn2)
        _cache_leg[ key_filedata] = (data1,data2)
    else:
        data1, data2 = data

    #Extract relevant polynomial coefficients:
    p0 = data1[f'polycoeff_xprime_to_y-{mode}_at_th0'][:]
    pdelta = data2[f'ydelta_{mode}'][:]

    result = dict(p0=p0, pdelta=pdelta)
    _cache_leg[key] = result
    return result

def _apply_split45( rawdata, rawdata45 ):
    import copy
    assert all( abs(a-b)/(abs(a)+abs(b))<1e-14
                for a,b in zip(rawdata['xvals'],rawdata45['xvals']))
    assert set(rawdata.keys()) == set(rawdata45.keys())
    assert set(rawdata.keys()) == set(['theta_2_calctime', 'theta_2_ypvals', 'theta_2_ypvals_maxerr', 'xvals'])
    res = { 'xvals' : [e for e in rawdata['xvals']] }
    for k in ['theta_2_calctime', 'theta_2_ypvals', 'theta_2_ypvals_maxerr']:
        d = dict( (th, yvals) for th,yvals in rawdata[k].items()
                  if abs(float(th)-45)>1e-9 )
        for th,yvals in rawdata45[k].items():
            d[th] = yvals
        res[k] = dict( sorted(d.items()) )
    return copy.deepcopy(res)

_cache_highx = {}
def load_highx( mode ):
    if mode in _cache_highx:
        return _cache_highx[mode]
    if '<raw>' not in _cache_highx:
        from .data import load_json_data
        _cache_highx['<raw>'] = load_json_data('bcdata_special_highx.json')
    _cache_highx[mode] = _cache_highx['<raw>'][mode]
    return _cache_highx[mode]

def load_fresnelref():
    from .fresnelref import load
    return load()
