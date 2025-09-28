from .mpmath import mpf

def poly_mult(a,b):
    assert len(a) and len(b)
    degree_a = len(a)-1
    degree_b = len(b)-1
    degree_new = degree_a + degree_b
    p = [ mpf(0), ]*(degree_new +1)
    for ia in range(len(a)):
        for ib in range(len(b)):
            #Add: a[ia]*b[ib] )*x^( ia+ib )
            p[ia+ib] += mpf(a[ia])*mpf(b[ib])
    while p[-1]==mpf(0):
        p = p[:-1]
    return p

def poly_add(a,b):
    p = [ mpf(0), ] * max(len(a),len(b))
    for i in range(len(p)):
        if i < len(a):
            p[i] += mpf(a[i])
        if i < len(b):
            p[i] += mpf(b[i])
    while p[-1]==mpf(0):
        p = p[:-1]
    return p

def raw_yprime_poly_to_y_poly( polycoeffs ):
    #If the coefficients define a given polynomial P which defines yp=P(xp),
    #then we return coefficients for a polynomial given by:
    #    P(xp)*(1+xp)**2*(1-xp)-(0.25*xp*(3.0-xp*xp)-0.5)
    #I.e. exactly what is needed to go directly from xp to y, rather than
    #from xp to yp as the original polynomial does.

    return [ float(e) for e in
             poly_add( poly_mult( polycoeffs, [1,1,-1,-1] ),
                       ['1/2','-3/4','0','1/4'] )
            ]
