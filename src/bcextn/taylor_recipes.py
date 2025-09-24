def taylor_y0_coeffs( mode ):
    from .mpmath import mp, mpf
    #5th order taylor coefficients for y(x,theta=0)
    if mode == 'primary':
        return [ mpf(1),
                 mpf("-33/35"),
                 mpf("4738/5775"),
                 mpf("-54120476/91216125"),
                 mpf("350769113251/962451740250"),
                 mpf("-276478446363113/1423053644512500") ]
    if mode == 'scndgauss':
        return [ mpf(1),
                 mp.sqrt(2)*mpf("-3/4"),
                 mp.sqrt(3)*mpf("8/15"),
                 mpf("-2/3"),
                 mp.sqrt(5)*mpf("32/175"),
                 mp.sqrt(3)*mp.sqrt(2)*mpf("-4/45") ]
    if mode == 'scndlorentz':
        return [ mpf(1),
                 mpf("-1"),
                 mpf("16/15"),
                 mpf("-80/81"),
                 mpf("64/81"),
                 mpf("-224/405") ]
    if mode == 'scndfresnel':
        return [ mpf(1),
                 mpf("-1"),
                 mpf("22/25"),
                 mpf("-604/945"),
                 mpf("15619/39690"),
                 mpf("-655177/3118500") ]
    assert False

def taylor_ydelta_coeffs( mode ):
    from .mpmath import mp, mpf
    #5th order taylor coefficients for y(x,theta=pi)-y(x,theta=0):
    if mode == 'primary':
        return [ mpf(0), mpf(0),
                 mpf("2369/5775"),
                 mpf("-108240952/91216125"),
                 mpf("350769113251/148069498500"),
                 mpf("-11888573193613859/2846107289025000") ]
    if mode == 'scndgauss':
        return [ mpf(0), mpf(0),
                 mp.sqrt(3)*mpf("4/15"),
                 mpf("-4/3"),
                 mp.sqrt(5)*mpf("208/175"),
                 mp.sqrt(3)*mp.sqrt(2)*mpf("-86/45") ]
    if mode == 'scndlorentz':
        return [ mpf(0), mpf(0),
                 mpf("8/15"),
                 mpf("-160/81"),
                 mpf("416/81"),
                 mpf("-4816/405") ]
    if mode == 'scndfresnel':
        return [ mpf(0), mpf(0),
                 mpf("11/25"),
                 mpf("-1208/945"),
                 mpf("203047/79380"),
                 mpf("-28172611/6237000") ]
    assert False
