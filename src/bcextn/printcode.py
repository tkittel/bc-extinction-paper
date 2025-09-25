import math

languages = ['py','c','cpp']

def doit( lang, do_print = True ):
    from .data import load_json_data
    from .mpmath import mp
    from .taylor_recipes import taylor_ydelta_coeffs, taylor_y0_coeffs
    from .load import load_legendre
    from .new_recipes import recipe_target_prec

    assert lang in languages

    o = []
    eps_lux = recipe_target_prec(lux=True)
    eps_std = recipe_target_prec(lux=False)

    modes = ['primary','scndgauss','scndlorentz','scndfresnel']

    if lang=='c':
        o.append('#include <math.h>')
    elif lang=='cpp':
        o.append('#include <cmath>')
    else:
        assert lang=='py'

    for lux  in [False,True]:
        eps = recipe_target_prec(lux=lux)
        ndigits = -math.log10(eps)
        assert abs(round(ndigits)-ndigits)<1e-7
        ndigits = round(ndigits)
        for e in textbox(lang,
                         'BC2025 %s Recipes'%('Luxury' if lux else 'Standard'),
                         f'Precision guarantee for x<1000: {ndigits} significant digits',
                         'Reference: T. Kittelmann 2025 (publication in preparation)'):
            o.append(e)
        for mode in modes:
            luxpostfix = '_lux' if lux else ''
            fctname = f'bc2025_y_{mode}{luxpostfix}'
            ns = 'std::' if lang=='cpp' else ''
            taylorthr = '0.01' if lux else '0.1'

            ndigits_taylor = 14 if lux else 5#FIXME: Justify this choice! AND MAKE SURE WE USE THE SAME IN curves.py!!!!

            def fmt_taylor_ydelta():
                p = taylor_ydelta_coeffs( mode )
                assert len(p)==6 and p[0]==0 and p[1]==0
                #lowest order is x^2
                return 'u*u*(%s)'%(fmt_horner([mp.nstr(e,n=ndigits_taylor) for e in p[2:]],varname='u'))
            def fmt_taylor_y0():
                p = taylor_y0_coeffs( mode )
                assert len(p)==6 and p[0]==1
                return '%s'%(fmt_horner([mp.nstr(e,n=ndigits_taylor) for e in p]))
            def fmt_yp0():
                p = load_legendre(mode,is_lux=lux)['p0']
                assert 5 < len(p)< (50 if lux else 15)
                return '%s'%(fmt_horner(['%.15g'%e for e in p],varname='xp'))
            def fmt_ydelta():
                p = load_legendre(mode,is_lux=lux)['pdelta']
                assert 5 < len(p)< (50 if lux else 15)
                return '%s'%(fmt_horner(['%.15g'%e for e in p],varname='up'))

            if lang=='py':
                o.append(f'def {fctname}( x, sintheta ):' )
                o.append( '    if x < %s:'%taylorthr)
                o += _wrap_on_plus(f'        y0 = ( {fmt_taylor_y0()} )',indent=' '*15)
                o.append( '    else:')
                o.append( '        if x > 1e3:')
                if mode == 'scndgauss':
                    o.append(f'            return {fctname}(1e3,sintheta)*(x*1e-3)**-0.93')
                else:
                    o.append(f'            return {fctname}(1e3,sintheta)*(x*1e-3)**-0.5')
                o.append( '        xs = x**0.5')
                o.append( '        xp = (xs-1.0)/(xs+1.0)')
                o += _wrap_on_plus(f'        yp0 = ( {fmt_yp0()} )',indent=' '*16)

                o.append( '        y0 = yp0*(1.0+xp)**2*(1-xp)-(0.25*xp*(3.0-xp*xp)-0.5)')
                o.append( '    s = sintheta**0.5')
                o.append( '    u = x*s')
                o.append( '    if u < %s:'%taylorthr)
                o += _wrap_on_plus(f'        ydelta = ( {fmt_taylor_ydelta()} )',indent=' '*19)
                o.append( '    else:')
                o.append( '        us = u**0.5')
                o.append( '        up = (us-1.0)/(us+1.0)')
                o += _wrap_on_plus(f'        ydelta = ( {fmt_ydelta()} )',indent=' '*19)
                o.append( '    return y0 + sintheta*s*ydelta')
            else:
                assert lang in ('c','cpp')
                o.append(f'double {fctname}( double x, double sintheta )')
                o.append( '{')
                o.append( '  double y0;')
                o.append( '  if ( x < %s ) {'%taylorthr)
                o += _wrap_on_plus(f'    y0 = {fmt_taylor_y0()};',indent=' '*6,rjust=False)
                #o.append(f'    y0 = 1.0 + y0minus1_taylor5_{mode}(x);')
                o.append( '  } else {')
                o.append( '    if ( x > 1e3 )')
                if mode == 'scndgauss':
                    o.append(f'      return {fctname}(1e3,sintheta)*{ns}pow(x*1e-3,-0.93);')
                else:
                    o.append(f'      return {fctname}(1e3,sintheta)*{ns}sqrt(1e3/x);')


                o.append(f'    const double xs = {ns}sqrt(x);')
                o.append( '    const double xp = (xs-1.0)/(xs+1.0);')
                o += _wrap_on_plus(f'    const double yp0 = {fmt_yp0()};',indent=' '*6,rjust=False)
                o.append( '    const double xpp1 = 1.0+xp;')
                o.append( '    y0 = yp0*xpp1*xpp1*(1.0-xp)-(0.25*xp*(3.0-xp*xp)-0.5);')
                o.append( '  }')
                o.append(f'  const double s = {ns}sqrt(sintheta);')
                o.append( '  const double u = x*s;')
                o.append( '  double ydelta;')
                o.append( '  if ( u < %s ) {'%taylorthr)
                o += _wrap_on_plus(f'    ydelta = {fmt_taylor_ydelta()};',indent=' '*6,rjust=False)
                o.append( '  } else {')
                o.append(f'    const double us = {ns}sqrt(u);')
                o.append( '    const double up = (us-1.0)/(us+1.0);')
                o += _wrap_on_plus(f'    ydelta = {fmt_ydelta()};',indent=' '*6,rjust=False)
                o.append( '  }')
                o.append( '  return y0 + sintheta*s*ydelta;')
            if lang!='py':
                o.append('}')
            o.append('')

    refdata = load_json_data('bc2025_reference_x_sintheta_yp.json')
    assert set(refdata.keys()) == set(modes)

    for e in textbox(lang,
                     'Test code for BC2025 Recipes',
                     'Test function to call: bc2015_test_implementation()'):
        o.append(e)

    if lang=='cpp':
        o.append( '#include <iostream>')
        o.append( '#include <sstream>')
        o.append( '#include <stdexcept>')

    if lang=='c':
        o.append( '#include "stdio.h"')
        o.append( '#include "stdlib.h"')
        o.append('')
        for mode in modes:
            fctname = f'bc2025_y_{mode}'
            fctname_lux = f'{fctname}_lux'
            o.append(f'void bc2015_test_implementation_{mode}( double x, double sinth, double refval)')
            o.append( '{')
            o.append(f'  const double val = {fctname}(x,sinth);')
            o.append(f'  const double val_lux = {fctname}_lux(x,sinth);')
            o.append(f'  if (!(fabs(val-refval) <= {eps_std}*fmin(refval,1.0-refval))) '+'{')
            o.append(f'    printf("bc2015_test_implementation: Failure in {fctname}(%g,%g)\\n",x,sinth);')
            o.append( '    exit(1);')
            o.append( '  }')
            o.append(f'  if (!(fabs(val_lux-refval) <= {eps_lux}*fmin(refval,1.0-refval))) '+'{')
            o.append(f'    printf("bc2015_test_implementation: Failure in {fctname_lux}(%g,%g)\\n",x,sinth);')
            o.append( '    exit(1);')
            o.append( '  }')
            o.append( '}')

    if lang=='py':
        o.append('def bc2015_test_implementation():')
    else:
        assert lang in ('c','cpp')
        o.append('void bc2015_test_implementation()')
        o.append('{')

    if lang=='cpp':
        for imode,mode in enumerate(modes):
            fctname = f'bc2025_y_{mode}'
            fctname_lux = f'{fctname}_lux'
            if imode:
                o.append('')
            o.append(f'  auto test_{mode} = []( double x, double sinth, double refval)')
            o.append( '  {')
            o.append(f'    const double val = {fctname}(x,sinth);')
            o.append(f'    const double val_lux = {fctname}_lux(x,sinth);')
            o.append(f'    if (!({ns}fabs(val-refval) <= {eps_std}*{ns}fmin(refval,1.0-refval))) '+'{')
            o.append( '      std::ostringstream ss;;')
            o.append(f'      ss << "bc2015_test_implementation: Failure in {fctname}("<<x<<","<<sinth<<")\\n";')
            o.append( '      throw std::runtime_error(ss.str());')
            o.append( '    }')
            o.append(f'    if (!({ns}fabs(val_lux-refval) <= {eps_lux}*{ns}fmin(refval,1.0-refval))) '+'{')
            o.append( '      std::ostringstream ss;;')
            o.append(f'      ss << "bc2015_test_implementation: Failure in {fctname_lux}("<<x<<","<<sinth<<")\\n";')
            o.append( '      throw std::runtime_error(ss.str());')
            o.append( '    }')
            o.append( '  };')
        for mode in modes:
            o.append('')
            o.append(f'  std::cout<<"Testing {fctname} + {fctname_lux}"<<std::endl;;')
            for x,sinth,y in refdata[mode]:
                o.append( '  test_%s(%.14g,%.14g,%.18g);'%(mode,x,sinth,y))
        o.append('  std::cout<<"All tests completed"<<std::endl;')

    if lang=='c':
        for mode in modes:
            fctname = f'bc2025_y_{mode}'
            fctname_lux = f'{fctname}_lux'
            o.append('')
            o.append(f'  printf("Testing {fctname} + {fctname_lux}\\n");')
            for x,sinth,y in refdata[mode]:
                o.append(f'  bc2015_test_implementation_{mode}(%.14g,%.14g,%.18g);'%(x,sinth,y))
        o.append('  printf("All tests completed\\n");')
        o.append('}')

    for mode in modes:
        fctname = f'bc2025_y_{mode}'
        fctname_lux = f'{fctname}_lux'
        if lang=='py':
            o.append( '    def t(x,sinth,refval):')
            o.append(f'        val = {fctname}(x,sinth)')
            o.append(f'        val_lux = {fctname}_lux(x,sinth)')
            o.append(f'        if not ( abs(val-refval) <= {eps_std}*min(refval,1.0-refval) ):')
            o.append(f'            raise RuntimeError("bc2015_test_implementation: Failure in {fctname}(%g,%g)"%(x,sinth))')
            o.append(f'        if not (abs(val_lux-refval) <= {eps_lux}*min(refval,1.0-refval)):')
            o.append(f'            raise RuntimeError("bc2015_test_implementation: Failure in {fctname_lux}(%g,%g)"%(x,sinth))')
            o.append(f'    print("Testing {fctname} + {fctname_lux}")')
            for x,sinth,y in refdata[mode]:
                o.append( '    t(%.14g,%.14g,%.18g)'%(x,sinth,y))

    if lang=='py':
        o.append('    print("All tests completed")')

    if lang=='cpp':
        o.append('}')

    for e in textbox(lang,'Hook for executing file as test application'):
        o.append(e)

    o.append('')

    if lang=='py':
        o.append('if __name__== "__main__":')
        o.append('    bc2015_test_implementation()')
    if lang in ('c','cpp'):
        o.append('#ifndef BC2025_NO_MAIN')
    if lang=='c':
        o.append('int main( int argc, char ** argv )')
        o.append('{')
        o.append('  (void)argc;')
        o.append('  (void)argv;')
        o.append('  bc2015_test_implementation();')
        o.append('  return 0;')
        o.append('}')
    if lang=='cpp':
        o.append('int main()')
        o.append('{')
        o.append('  bc2015_test_implementation();')
        o.append('  return 0;')
        o.append('}')
    if lang in ('c','cpp'):
        o.append('#endif')
    if print:
        print('\n'.join(o))
    return o

def fmt_horner( polycoeffs, varname = 'x' ):
    #https://en.wikipedia.org/wiki/Horner%27s_method
    def c(i):
        e = polycoeffs[i]
        e = e if isinstance(e,str) else '%.19g'%e
        return e+'.0' if e.isdigit() or (e[0]=='-' and e[1:].isdigit()) else e
    n,j = len(polycoeffs),f'+{varname}*'
    assert n>=1
    if n==1:
        return c(0)
    if n==2:
        return c(0)+j+c(1)
    return c(0)+j + ''.join('(%s%s'%(c(i),j) for i in range(1,n-1))+c(n-1)+')'*(n-2)

def _wrap_on_plus( s, width = 80, indent='    ',rjust=False):
    if len(s)<=width:
        return [s]
    tmp = ''
    out = []
    for i,p in enumerate(s.split('+')):
        if i:
            p = '+'+p
        if len(tmp)>len(indent) and len(tmp)+len(p)>width:
            #flush tmp:
            if rjust and out:
                tmp = tmp.rjust(width)
            out += [tmp]
            tmp = indent
        tmp += p
    if tmp:
        if rjust and out:
            tmp = tmp.rjust(width)
        out+=[tmp]
    return out

def textbox(lang,*msg_lines):
    yield ''
    if lang in ('py','cpp'):
        c = '#' if lang=='py' else '/'
        yield c*80
        yield c*80
        for msg in msg_lines:
            yield c*3+'  '+msg.center(70)+'  '+c*3
        yield c*80
        yield c*80
    elif lang=='c':
        hr = '/' + '*'*78 + '/'
        yield hr
        yield hr
        for msg in msg_lines:
            yield '/*   '+msg.center(70)+'   */'
        yield hr
        yield hr
    else:
        assert False
    yield ''
