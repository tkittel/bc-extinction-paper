
# Load data previously created with the createdata script

def load_xscan( mode = 'primary' ):
    if mode.startswith('curve::'):
        return load_curves( mode[len('curve::'):], 'xscan' )
    assert mode in ['primary','scndfresnel','scndlorentz','scndgauss']
    if mode == 'primary':
        return _loadcache('bcdata.json')
    else:
        return _loadcache(f'bcdata_{mode}.json')

def load_thetascan( mode = 'primary' ):
    if mode.startswith('curve::'):
        return load_curves( mode[len('curve::'):], 'thetascan' )
    assert mode in ['primary','scndfresnel','scndlorentz','scndgauss']
    if mode == 'primary':
        return _loadcache('bcdata_thetascan.json')
    else:
        return _loadcache(f'bcdata_{mode}_thetascan.json')

def _impl_load_tablescan_origpts( tablename ):
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

_cache={}
def _loadcache( fn ):
    from .data import load_json_data as _load
    d = _cache.get(fn)
    if not d:
        d = _prepdata(_load(fn))
        _cache[fn] = d
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
    assert mode in 'xscan','thetascan'
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
    refdata = ( load_xscan('primary')
                if mode == 'xscan'
                else load_thetascan('primary') )
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
