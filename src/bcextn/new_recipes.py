

def recipe_taylor_cutoff( lux ):
    return 0.01 if lux else 0.1

def recipe_highx_cutoff():
    return 1000.0

def recipe_highx_pow( mode ):
    return 0.933 if mode == 'scndgauss' else 0.5

def recipe_taylor_order():
    return 5

def recipe_target_prec( lux ):
    return 1e-6 if lux else 1e-3

def recipe_target_prec_each_main_term( lux ):
    #terms y0 or ydelta must both share the precision, so we add a factor of
    #0.5.
    return recipe_target_prec( lux ) * 0.5
