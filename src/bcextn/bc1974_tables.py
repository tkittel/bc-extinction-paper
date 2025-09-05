
def _do_load(fn, yp_or_ys):
    assert yp_or_ys in ('yp','ys')
    from .data import load_json_data
    import numpy as np
    def array( x ):
        return np.asarray( x, dtype=float )
    all_values = []
    raw = load_json_data(fn)
    assert set(raw.keys()) == set(['xvals','data'])
    xvals = array( raw['xvals'] )
    sinthvals = []
    x_sinth_to_yps = {}
    for e in raw['data']:
        #print( set(e.keys()))
        #print(set(['sinth',f'{yp_or_ys}*1e4']))
        assert set(e.keys()) == set(['sinth',f'{yp_or_ys}*1e4'])
        sinth = float(e['sinth'])
        sinthvals.append( sinth )
        ypsvals = array( [ yps*1e-4 for yps in e[f'{yp_or_ys}*1e4'] ] )
        assert len(ypsvals) == len(xvals)
        for x, yps in zip(xvals,ypsvals):
            all_values.append( ( x, sinth, yps ) )
            key = ( float(x), sinth )
            assert key not in x_sinth_to_yps
            x_sinth_to_yps[key] = float(yps)

    #FIXME: 'yp' keys to yy?
    return { 'xvals' : xvals,
             'sinthvals' : array(sinthvals),
             'x_sinth_to_yp' : x_sinth_to_yps,
             'x_sinth_yp_list' : array( sorted(all_values) ) }

_cache_table1 = [None]
def table1():
    if _cache_table1[0] is None:
        _cache_table1[0] = _do_load('bc1974_table1.json','yp')
    return _cache_table1[0]

_cache_table3 = [None]
def table3():
    if _cache_table3[0] is None:
        _cache_table3[0] = _do_load('bc1974_table3.json','ys')
    return _cache_table3[0]

_cache_table4 = [None]
def table4():
    if _cache_table4[0] is None:
        _cache_table4[0] = _do_load('bc1974_table4.json','ys')
    return _cache_table4[0]


