import numpy as np
import math
from .mpmath import mp, mpf
from . import trf

kDeg2Rad = np.pi/180

def mode_curves( mode ):
    assert mode in ['primary','scndfresnel','scndlorentz','scndgauss']
    if mode=='primary':
        return ( ClassicCurve_Primary,
                 UpdatedClassicCurve_Primary,
                 ProposedCurve_Primary,
                 ProposedLuxCurve_Primary)
    if mode=='scndgauss':
        return ( ClassicCurve_ScndGauss,
                 UpdatedClassicCurve_ScndGauss,
                 ProposedCurve_ScndGauss,
                 ProposedLuxCurve_ScndGauss )
    if mode=='scndlorentz':
        return ( ClassicCurve_ScndLorentz,
                 UpdatedClassicCurve_ScndLorentz,
                 ProposedCurve_ScndLorentz,
                 ProposedLuxCurve_ScndLorentz )
    assert mode=='scndfresnel'
    return ( ClassicCurve_ScndFresnel,
             UpdatedClassicCurve_ScndFresnel,
             ProposedCurve_ScndFresnel,
             ProposedLuxCurve_ScndFresnel)

#####################################################
############# ACTUAL CURVE OBJECTS ##################
#####################################################

class ClassicCurve_Primary:

    def __init__(self):
        self._defaultparams = self.bc1974_params()

    def __call__( self, x, theta, params = None ):
        if params is None:
            params = self._defaultparams
        A, B = self.calcThetaPars( theta, params )
        #Using same recipe form as primary:
        return yp_bc1974_curve( x, A, B )

    def bc1974_params( self ):
        return [0.2,0.45,0.22,-0.12]

    def p0( self ):
        return [0.56,0.54,0.6,-0.2]#self.bc1974_params()[:]

    def bestfit_params( self ):
        if not hasattr(self,'_fitparams'):
            setattr( self, '_fitparams',
                     load_fitted_curve_parameters('primary')['classic'][:])
        assert len(self._fitparams)==4
        return self._fitparams

    def calcThetaPars( self, theta, params ):
        assert len(params)==4
        p = params
        u = math.cos( 2*float(theta)*kDeg2Rad )
        A = p[0]+p[1]*u
        B = p[2]+p[3]*((0.5-u)**2)
        return A, B

class ClassicCurve_ScndGauss:

    def __init__(self):
        self._defaultparams = self.bc1974_params()

    def __call__( self, x, theta, params = None ):
        if params is None:
            params = self._defaultparams
        A, B = self.calcThetaPars( theta, params )
        return ys_bc1974_curve_scndgauss( x, A, B )

    def bc1974_params( self ):
        return [ 0.58, 0.48, 0.24, 0.02, -0.025 ]

    def p0( self ):
        return [0.94422166,0.56110532, 0.18150111, 0.03482274, 0.00949109]
        #return self.bc1974_params()

    def bestfit_params( self ):
        if not hasattr(self,'_fitparams'):
            setattr( self, '_fitparams',
                     load_fitted_curve_parameters('scndgauss')['classic'][:])
        assert len(self._fitparams)==5
        return self._fitparams

    def calcThetaPars( self, theta, params ):
        assert len(params)==5
        p = params
        u = math.cos( 2*float(theta)*kDeg2Rad )#fixme : mpf?
        A = p[0]   + p[1]*u  + p[2]*u**2
        B = p[3]   + p[4]*u
        return A, B

class ClassicCurve_ScndLorentz:

    def __init__(self):
        self._defaultparams = self.bc1974_params()

    def __call__( self, x, theta, params = None ):
        if params is None:
            params = self._defaultparams
        A, B = self.calcThetaPars( theta, params )
        #Using same recipe form as primary:
        return yp_bc1974_curve( x, A, B )

    def bc1974_params( self ):
        return [0.025,0.285,-0.45,0.15,-0.2,0.75]

    def p0( self ):
        #Manually found to give reasonable performance:
        return [0.0859,0.27381,-0.9,0.38,-0.14,0.56]

    def bestfit_params( self ):
        if not hasattr(self,'_fitparams'):
            setattr( self, '_fitparams',
                     load_fitted_curve_parameters('scndlorentz')['classic'][:])
        assert len(self._fitparams)==6
        return self._fitparams

    def calcThetaPars( self, theta, params ):
        assert len(params)==6
        p = params
        u = math.cos( 2*float(theta)*kDeg2Rad )#fixme : mpf?
        A = p[0]+p[1]*u
        B = ( (p[2]*u)
              if u <= 0
              else (p[3] +p[4]*((p[5]-u)**2) ) )
        return A, B

class ClassicCurve_ScndFresnel:

    def __init__(self):
        self._defaultparams = self.bc1974_params()

    def __call__( self, x, theta, params = None ):
        if params is None:
            params = self._defaultparams
        A, B = self.calcThetaPars( theta, params )
        #Using same recipe form as primary:
        return yp_bc1974_curve( x, A, B )

    def bc1974_params( self ):
        return [0.48,0.6,0.2,-0.06,0.2]

    def p0( self ):
        return self.bc1974_params()[:]

    def bestfit_params( self ):
        if not hasattr(self,'_fitparams'):
            setattr( self, '_fitparams',
                     load_fitted_curve_parameters('scndfresnel')['classic'][:])
        assert len(self._fitparams)==5
        return self._fitparams

    def calcThetaPars( self, theta, params ):
        assert len(params)==5
        p = params
        u = math.cos( 2*float(theta)*kDeg2Rad )#fixme : mpf?
        A = p[0]+p[1]*u
        B = p[2]+p[3]*((p[4]-u)**2)
        return A, B

#Refitted classic curves:

class UpdatedClassicCurve_Primary( ClassicCurve_Primary ):
    def __init__(self):
        super().__init__()
        self._defaultparams = self.bestfit_params()
class UpdatedClassicCurve_ScndGauss( ClassicCurve_ScndGauss ):
    def __init__(self):
        super().__init__()
        self._defaultparams = self.bestfit_params()
class UpdatedClassicCurve_ScndLorentz( ClassicCurve_ScndLorentz ):
    #FIXME: Need to follow up on the issues near theta=45 seen in ./bin/plot breakdown scndlorentz
    def __init__(self):
        super().__init__()
        self._defaultparams = self.bestfit_params()
class UpdatedClassicCurve_ScndFresnel( ClassicCurve_ScndFresnel ):
    def __init__(self):
        super().__init__()
        self._defaultparams = self.bestfit_params()


#Proposed curves:

_cache_stdlegfit = {}
class ProposedCurve_StdLegFit:

    def __do_init(self):
        assert not hasattr(self,'_fcts')
        datafilename, key = self._poly_coeff_src()
        if datafilename not in _cache_stdlegfit:
            from .data import load_json_data
            _cache_stdlegfit[datafilename] = load_json_data(datafilename)
        data = _cache_stdlegfit[datafilename]
        p0_fct  = np.polynomial.polynomial.Polynomial( data[f'yprime_{key}_0'] )
        p90_fct = np.polynomial.polynomial.Polynomial( data[f'yprime_{key}_90'] )
        setattr( self,'_fcts',( p0_fct, p90_fct ) )

    def __call__( self, x, theta ):
        if not hasattr(self,'_fcts'):
            self.__do_init()
        fct_yp_theta0, fct_yp_theta90 = self._fcts

        sinth =  math.sin( float(theta)*kDeg2Rad )
        sqrt_sinth = math.sqrt(sinth)
        #need y0 at x=x
        xp = trf.transform_to_xprime( x )
        yp0 = fct_yp_theta0( xp )
        y0 = trf.y_unprime( xp, yp0 )
        #need y0 and y90 at x=x*sqrt(sintheta):=u_x
        u_x = x * sqrt_sinth
        u_xp = trf.transform_to_xprime( u_x )
        u_yp0 = fct_yp_theta0( u_xp )
        u_yp90 = fct_yp_theta90( u_xp )
        u_y0 = trf.y_unprime( u_xp, u_yp0 )
        u_y90 = trf.y_unprime( u_xp, u_yp90 )
        #Now combine as, taking care to avoid numerical issues by using fsum:
        #y(x) = y0(x) + sinth^(3/2)*(y90(u_x)-y0(u_x))
        s32 = sinth*sqrt_sinth
        return np_fsum3(y0,s32*u_y90,-s32*u_y0)

class ProposedCurve_Primary(ProposedCurve_StdLegFit):
    def _poly_coeff_src( self ):
        return 'legendre_coefficients.json', 'primary'
class ProposedLuxCurve_Primary(ProposedCurve_StdLegFit):
    def _poly_coeff_src( self ):
        return 'legendre_coefficients_lux.json', 'primary'
class ProposedCurve_ScndGauss(ProposedCurve_StdLegFit):
    def _poly_coeff_src( self ):
        return 'legendre_coefficients.json', 'scndgauss'
class ProposedLuxCurve_ScndGauss(ProposedCurve_StdLegFit):
    def _poly_coeff_src( self ):
        return 'legendre_coefficients_lux.json', 'scndgauss'
class ProposedCurve_ScndLorentz(ProposedCurve_StdLegFit):
    def _poly_coeff_src( self ):
        return 'legendre_coefficients.json', 'scndlorentz'
class ProposedLuxCurve_ScndLorentz(ProposedCurve_StdLegFit):
    def _poly_coeff_src( self ):
        return 'legendre_coefficients_lux.json', 'scndlorentz'
class ProposedCurve_ScndFresnel(ProposedCurve_StdLegFit):
    def _poly_coeff_src( self ):
        return 'legendre_coefficients.json', 'scndfresnel'
class ProposedLuxCurve_ScndFresnel(ProposedCurve_StdLegFit):
    def _poly_coeff_src( self ):
        return 'legendre_coefficients_lux.json', 'scndfresnel'

########################################################
############# Sabine model wrappers ####################
########################################################

class SabineCurve_Primary:

    def __call__( self, x, theta ):
        # Redoing Sabine's calculations leading to eq. 6.4.5.3 and 6.4.5.4, I (TK)
        # found the following exact form:
        #
        # El = exp(-y) * exp(-x) * ( I0(x)+I1(x) )
        #
        # Where I0(x) and I1(x) are the modified bessel functions of the first kind,
        # for nu=0 and nu=1 respectively. In mpmath these are given as besseli(0,x)
        # and besseli(1,x) (in sagemath the function name is bessel_I).
        x = mpf(x)
        x *= mpf('3/2')#FIXME: Different x definition. CHECK!!!!
        s =  mp.sin(mp.radians(mpf(theta)))
        ssq = s*s
        csq =  max(0.0,1.0 - ssq)
        El = mp.exp( -x ) * ( mp.besseli(0,x) + mp.besseli(1,x) )
        Eb = mpf(1)/mp.sqrt(mpf(1)+mpf(x))
        return El*csq+Eb*ssq

class SabineCurve_ScndRec:

    def __call__( self, x, theta ):
        x = mpf(x)
        x *= mpf('3/2')#FIXME: Different x definition. CHECK!!!!
        s =  mp.sin(mp.radians(mpf(theta)))
        ssq = s*s
        csq =  max(0.0,1.0 - ssq)

        one = mpf(1)
        invx = one/x

        El = mpf('1/2') * invx * ( one - mp.exp(-2*x) )
        Eb = one / ( one + x )
        return float( El*csq+Eb*ssq )

class SabineCurve_ScndTri:

    def __call__( self, x, theta ):
        x = mpf(x)
        x *= mpf('3/2')#FIXME: Different x definition. CHECK!!!!
        s =  mp.sin(mp.radians(mpf(theta)))
        ssq = s*s
        csq =  max(0.0,1.0 - ssq)

        one = mpf(1)
        invx = one/x

        El = invx * ( one - mpf('1/2')*invx * ( one-mp.exp(-mpf(2)*x) ) )
        Eb = ( mpf(2)*invx**2 ) * ( x - mp.log(one+x) )
        return float( El*csq+Eb*ssq )

#####################################################
############# STANDALONE CURVE FUNCTIONS ############
#####################################################

@np.vectorize
def yp_bc1974_curve( x, A, B ):
    #Note, this function also used by secondary Lorentz + Fresnel
    k = 1 + B*x
    k2 = 1.0+2.0*x+safediv(A*x*x,k)
    return safediv( 1.0, safesqrt( k2 ) )

@np.vectorize
def yp_bc1974_curve_breakdown( x, A, B ):
    #Same, but only testing for strict numerical breakdown (zero division, sqrt
    #of negative, yp outside [0,1]):
    k = 1 + B*x
    if abs(k) <= 1e-99:
        return True
    k2 = 1.0+2.0*x+safediv(A*x*x,k)
    if k2 <= 1e-99:
        return True
    res = safediv( 1.0, safesqrt( k2 ) )
    return not ( 0.0 <= res <= 1.0 )

@np.vectorize
def ys_bc1974_curve_scndgauss( x, A, B ):
    k = 1 + B*x
    #Notice: 2.12 which is mentioned in the paper and comes from
    #        3/sqrt(2) ~= 2.1213...
    k2 = 1.0+2.12*x+safediv(A*x*x,k) # NOTICE 2.12 here
    return safediv( 1.0, safesqrt( k2 ) )

#####################################################
############# PLUMBING AND UTILS ####################
#####################################################

@np.vectorize
def safediv( a, b ):
    return a/b if abs(b)>1e-99 else 1e99

@np.vectorize
def safesqrt( x ):
    return math.sqrt(max( 0.0, x ))

_cache_fcp = {}
def load_fitted_curve_parameters(mode):
    assert mode in ['primary','scndfresnel','scndlorentz','scndgauss']
    key = mode
    if key in _cache_fcp:
        return _cache_fcp[key]

    #Clip values for appearing nicer in the paper. But also make sure that we
    #actually use the clipped values in performance plots!
    def fmt( key, x ):
        if 'orig' in key:
            #FIXME: return float('%.3g'%x)
            return float('%.5g'%x)
        else:
            #FIXME: return float('%.4g'%x)
            return float('%.7g'%x)
    from .data import load_json_data
    fn = 'global_refitted_classic_curves.json'
    if mode != 'primary':
        fn = f'global_refitted_classic_curves_{mode}.json'
    data = load_json_data(fn)
    for k in data.keys():
        data[k] = [ fmt(k,e) for e in data[k] ] if data[k] is not None else None
    _cache_fcp[key] = data
    return data

@np.vectorize
def np_fsum3( a, b, c ):
    return math.fsum([a,b,c])
