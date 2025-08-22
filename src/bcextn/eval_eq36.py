eq36_result_tolerance = 1e-7

def eval_eq36( theta_degree, x, eps = None ):
    from .eq36_lowx_taylor import taylor_lowx_eq36
    from .eq36_numint import numint_eq36
    if eps is None:
        eps = eq36_result_tolerance
    t = taylor_lowx_eq36( theta_degree, x, eps )
    if t is not None:
        return t
    return numint_eq36( theta_degree, x, eps )

def benchmark():
    from .mpmath import mpf
    import time
    res = []
    for theta_degree in 0, 30, 60, 90:
        for x in [ 0, 1e-20, 1e-3, 0.2, 0.5, 1.0, 1.01, 2.0, 3, 30, 100, 1e5, 1e9 ]:
            print()
            print('-'*70)
            print(f' theta={theta_degree:g},  x={x:g}')
            t0 = time.time()
            v, err = eval_eq36( theta_degree=theta_degree, x = x )
            t0 = time.time() - t0
            print('   RESULT: %g'%v)
            print('   1-RESULT: %g'%(mpf(1)-v))
            print('           +- %e'%err)
            print(f'   TIME: {t0:.2g}s')
            res.append( (t0, theta_degree, x, v, err ) )
