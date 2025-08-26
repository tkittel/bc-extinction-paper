
# Load data previously created with the createdata script

def load_xscan():
    return _loadcache('bcdata.json')

def load_thetascan():
    return _loadcache('bcdata_thetascan.json')

def load_table1scan_origpts():
    cachekey = 'load_table1scan_origpts_cache'
    if cachekey in _cache:
        return _cache[cachekey]

    from .bc1974_tables import table1
    from .plotutils import array
    import numpy as np
    import copy

    orig_table_data = table1()
    orig_xvals = orig_table_data['xvals']
    orig_th_vals = [ np.asin(float(e))*180/np.pi
                     for e in orig_table_data['sinthvals'] ]
    def x_ok( xval ):
        return any(abs(float(xval)-float(e))<1e-13 for e in orig_xvals)
    def th_ok( thval ):
        return any(abs(float(thval)-float(e))<1e-13 for e in orig_th_vals)
    d_all = copy.deepcopy(load_table1scan_allpts())
    mask = np.array( [x_ok(e) for e in d_all['xvals']],dtype = bool )
    d_filtered = {}
    d_filtered['xvals'] = array([ e for e in d_all['xvals'] if x_ok(e) ])
    for k in ['theta_2_calctime', 'theta_2_ypvals', 'theta_2_ypvals_maxerr']:
        d_filtered[k] = dict( (k,v[mask]) for k,v in d_all[k].items() if th_ok(k) )
    d_filtered = _prepdata( d_filtered )
    assert set(d_all.keys()) == set(d_filtered.keys())
    _cache[cachekey] = d_filtered
    return d_filtered

def load_table1scan_allpts():
    return _loadcache('bcdata_table1pts.json')

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
