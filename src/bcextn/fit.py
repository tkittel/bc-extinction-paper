import numpy as np
from inspect import signature
import scipy.optimize

class FitFunc:
    def __init__( self, rawfct, p0 = None ):
        s = signature(rawfct)
        self.__name = rawfct.__name__
        self.__fct = np.vectorize(rawfct)
        self.__nparams = len(s.parameters)-1
        if p0 is None:
            p0 = [ 1.0, ]*self.__nparams
        else:
            assert len(p0) == self.__nparams
        self.__p0 = p0
    @property
    def nparams( self ):
        return self.__nparams
    @property
    def p0( self ):
        return self.__p0
    @property
    def name( self ):
        return self.__name
    def __call__(self, x, *p):
        return self.__fct(x,*p)
    def fit( self, xvals, yvals, p0 = None ):
        p0 = self.__p0 if p0 is None else p0
        assert len(p0)==self.__nparams
        if self.nparams:
            pres,_ = scipy.optimize.curve_fit( self.__fct, xvals, yvals, p0 = p0 )
        else:
            pres = []
        pres = [float(e) for e in pres]

        class FitResult:
            def __init__(self,name,fct,p):
                self.__n = name
                self.__f = fct
                self.__p = p
            @property
            def name( self ):
                return self.__n
            @property
            def parameters( self ):
                return self.__p
            def __call__(self, x):
                return self.__f(x,*(self.__p))
        return FitResult(self.__name,self.__fct,pres)

def fitfunction( func ):
    return FitFunc(func)
