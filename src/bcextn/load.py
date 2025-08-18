
# Load data previously created with the createdata script

def load_xscan():
    return _loadcache('bcdata.json')

def load_thetascan():
    return _loadcache('bcdata_thetascan.json')

def load_table1scan():
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
