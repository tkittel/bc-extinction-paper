import mpmath
import mpmath.libmp
assert mpmath.libmp.BACKEND=='gmpy', "even slower without gmpy"
mp = mpmath.mp
mpf = mp.mpf
#We aim for 100 digits of precision in f(eta) and phi*(sr) functions and terms
#of e.g. taylor expansions, so pick dps slighter above 100:
#mp.dps = 105
#Update: Increase even more for added safety:
mp.dps = 150


def mpf_pack_to_str( v ):
    return repr(mpf(v))
def mpf_unpack_from_str( v ):
    assert v.startswith("mpf('")
    assert v.endswith("')")
    return mpf(v[5:-2])
