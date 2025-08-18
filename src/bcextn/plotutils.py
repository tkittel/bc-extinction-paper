import matplotlib
import numpy as np

def array( x ):
    return np.asarray( x, dtype=float )

_cmap = matplotlib.colormaps['Spectral']
def th2color( theta_degree ):
    return _cmap( 1.0 - float( theta_degree ) / 90.0 )

