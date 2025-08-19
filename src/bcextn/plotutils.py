import matplotlib
import numpy as np

def array( x ):
    return np.asarray( x, dtype=float )

def th2color( theta_degree, cmap = 'Spectral' ):
    cmap = matplotlib.colormaps[cmap]
    return cmap( 1.0 - float( theta_degree ) / 90.0 )

