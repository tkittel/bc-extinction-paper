
_cache_table1 = [None]
def table1():
    if _cache_table1[0] is None:
        from .data import load_json_data
        import numpy as np
        def array( x ):
            return np.asarray( x, dtype=float )
        all_values = []
        raw = load_json_data('bc1974_table1.json')
        assert set(raw.keys()) == set(['xvals','data'])
        xvals = array( raw['xvals'] )
        #sinth2ypvals = {}
        sinthvals = []
        x_sinth_to_yp = {}
        for e in raw['data']:
            assert set(e.keys()) == set(['sinth','yp*1e4'])
            sinth = float(e['sinth'])
            sinthvals.append( sinth )
            ypvals = array( [ yp*1e-4 for yp in e['yp*1e4'] ] )
            assert len(ypvals) == len(xvals)
            #sinth2ypvals[sinth] = ypvals
            for x, yp in zip(xvals,ypvals):
                all_values.append( ( x, sinth, yp ) )
                key = ( float(x), sinth )
                assert key not in x_sinth_to_yp
                x_sinth_to_yp[key] = float(yp)

        _cache_table1[0] = { 'xvals' : xvals,
                             'sinthvals' : array(sinthvals),
                             'x_sinth_to_yp' : x_sinth_to_yp,
                             #'sinth_2_yp' : sinth2ypvals,
                             'x_sinth_yp_list' : array( sorted(all_values) ) }

    return _cache_table1[0]


