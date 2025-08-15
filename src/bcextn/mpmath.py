import mpmath
import mpmath.libmp
assert mpmath.libmp.BACKEND=='gmpy', "even slower without gmpy"
mp = mpmath.mp
mpf = mp.mpf
#We aim for 100 digits of precision in f(eta) and phi*(sr) functions and terms
#of e.g. taylor expansions, so pick dps slighter above 100:
mp.dps = 105
