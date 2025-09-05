from .mpmath import mp, mpf

class F_of_Eta:

    # Evaluate f(eta) for secondary extinction (BC, gauss) given by BC1974
    # eq. 40a, normalised so f(eta=0) = 1:
    #
    # f = exp( - k * eta^2 ), k = 9/(16pi)
    #

    def __init__( self ):
        k = -mpf('9/16') / mp.pi
        expfct = mp.exp
        self.__f = lambda eta : expfct ( k * eta * eta )

    def __call__( self, eta ):
        return self.__f( eta )
