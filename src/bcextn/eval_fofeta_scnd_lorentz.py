from .mpmath import mpf

class F_of_Eta:

    # Evaluate f(eta) for secondary extinction (BC, gauss) given by BC1974
    # eq. 40a, normalised so f(eta=0) = 1:
    #
    # f = 1 / ( 1 + eta^2 )
    #

    def __init__( self ):
        one = mpf(1)
        self.__f = lambda eta : one / ( one + eta*eta )

    def __call__( self, eta ):
        return self.__f( eta )
