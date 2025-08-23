from .mpmath import mp, mpf

def taylor_lowx_eq36( theta_degree, x, eps ):
    #Returns None if x is not low enough!
    sinth = mp.sin(mpf(theta_degree) * mp.pi/mpf(180) )
    x = mpf(x)
    eps = mpf(eps)*0.1 #safety
    taylor_fct = _taylor_lowx_eq36_create_fct( sinth )
    if x > mpf(10):
        #x is clearly too large
        return None
    if mpf(x)==0:
        return mpf(1), mpf(0)
    _, max_err = taylor_fct( x, only_error = True )
    if max_err > eps:
        #x is too large
        return None
    return taylor_fct( x )

def _taylor_lowx_eq36_create_fct( sinth ):
    # Taylor expansion to be used for small x. How small depends on the value of
    # sinth, so the code also returns an error bound. Code validated and
    # produced by the script bin/sagemath_taylor_eq36_lowx:

    m = mpf
    sinth = m(sinth)
    def s( power ):
        return sinth**m(power)
    c0 = m(1)
    c1 = m("-33/35")
    c2 = m("4738/5775")+m("2369/5775")*s("5/2")
    c3 = m("-23431877/39492724")+m("-23431877/19746362")*s(3)
    c4 = m("15818303/43402776")+m("19420368/8197883")*s("7/2")
    c5 = m("-34704149/178624650")+m("-34802273/8331614")*s(4)
    c6 = m("2645475/28932104")+m("22646387/3185480")*s("9/2")
    c7 = m("-7076575/183745624")+m("-39600909/3274684")*s(5)
    c8 = m("3778143/257403272")+m("73661039/3542884")*s("11/2")
    c9 = m("-889661/174179924")+m("-123768155/3419406")*s(6)
    c10 = m("741823/453654554")+m("247568569/3883970")*s("13/2")
    c11 = m("-1288490/2658208971")+m("-271100680/2391297")*s(7)
    c12 = m("215491/1610887604")+m("508831179/2502014")*s("15/2")
    c13 = m("-332827/9637359104")+m("-347802339/946354")*s(8)
    c14 = m("257465/30741493436")+m("66669808/99737")*s("17/2")
    c15 = m("-59163/30895090268")+m("-646003435/528329")*s(9)
    c16 = m("3659/8834425590")+m("4480847453/1993368")*s("19/2")
    c17 = m("-5521/64965502756")+m("-6183925831/1489698")*s(10)
    c18 = m("2581/155607420587")+m("4730431045/614594")*s("21/2")
    c19 = m("-1111/359936825592")+m("-1370036037/95651")*s(11)
    c20 = m("1243/2264668850872")+m("3757982407/140522")*s("23/2")
    c21 = m("-185/1979788177454")+m("-6280699297/125406")*s(12)
    c22 = m("69/4521826134880")+m("3357387539/35697")*s("25/2")
    c23 = m("-4/1670791357815")+m("-14064592969/79428")*s(13)
    c24 = m("17/47034812190742")+m("13816539855/41347")*s("27/2")
    c25 = m("-4/76072947608839")+m("-13446500609/21277")*s(14)
    c26 = m("1/135483517891918")+m("30124753055/25154")*s("29/2")
    c27 = m("-1/999042046860106")+m("-76908109889/33824")*s(15)
    c28 = m("1/7616623882587748")+m("2914716557/674")*s("31/2")
    c29 = m("-1/59972934662820664")+m("-11261686010/1367")*s(16)
    c30 = m("1/487218728395132032")+m("48175741043/3065")*s("33/2")
    def calc_taylor( x ):
        return c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*(c8+x*(c9+x*(c10+x*(c11+x*(c12+x*(c13+x*(c14+x*(c15+x*(c16+x*(c17+x*(c18+x*(c19+x*(c20+x*(c21+x*(c22+x*(c23+x*(c24+x*(c25+x*(c26+x*(c27+x*(c28+x*(c29+x*c30)))))))))))))))))))))))))))))
    def calc_taylor_errorbound( x ):
        # very conservative error bound, based on absolute size of
        # contributions of last three orders
        return (x**28)*(abs(c28)+x*(abs(c29)+x*abs(c30)))

    #Wrap up and return:
    def taylor( x, only_error = False ):
        v = None if only_error else calc_taylor(x)
        return v, calc_taylor_errorbound(x)
    return taylor

