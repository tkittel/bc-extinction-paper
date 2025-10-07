import math

languages = ['py','c','cpp','latex']

def doit( lang, outfile = None, include_main = True, include_testfcts = True ):
    from .data import load_json_data
    from .taylor_recipes import taylor_ydelta_coeffs, taylor_y0_coeffs
    from .load import load_legendre
    from .new_recipes import recipe_target_prec, recipe_highx_pow, recipe_taylor_cutoff
    import pathlib

    if lang=='latex':
        include_main = False
        include_testfcts = False

    if include_main:
        assert include_testfcts, "must include test functions if including main"

    assert lang in languages

    o_shared = [ [] ]
    o = o_shared[0]

    flush_if_latex_mode = False
    def flush_if_latex(tag):
        pass
    if outfile and 'TAGHERE' in outfile:
        if lang!='latex':
            raise SystemExit('TAGHERE multi output only supported for lang=latex')
        flush_if_latex_mode = True
        def flush_if_latex(tag):
            assert len(o_shared[0])>10
            c = '\n'.join(o_shared[0]) + '\n'
            o_shared[0].clear()
            pathlib.Path(outfile.replace('TAGHERE',tag)).write_text(c)

    eps_lux = recipe_target_prec(lux=True)
    eps_std = recipe_target_prec(lux=False)

    modes = ['primary','scndgauss','scndlorentz','scndfresnel']

    if lang=='c':
        o.append('#include <math.h>')
    elif lang=='cpp':
        o.append('#include <cmath>')

    for lux  in [False,True]:
        eps = recipe_target_prec(lux=lux)
        ndigits_recipe = -math.log10(eps)
        assert abs(round(ndigits_recipe)-ndigits_recipe)<1e-7
        ndigits_recipe = round(ndigits_recipe)
        for e in textbox(lang,
                         'BC2025 %s Recipes'%('Luxury' if lux else 'Standard'),
                         f'Precision guarantee for x<1000: Error less than 1e-{ndigits_recipe}*min(y,1-y)',
                         'Reference: T. Kittelmann 2025 (publication in preparation)'):
            o.append(e)
        for mode in modes:
            fctname = 'bc2025_y_%s%s'%(mode,'_lux' if lux else '')
            ns = 'std::' if lang=='cpp' else ''
            taylorthr = recipe_taylor_cutoff(lux)
            fpmode = 'luxrecipe' if lux else 'stdrecipe'
            pydelta_taylor = taylor_ydelta_coeffs( mode, fpmode = fpmode )
            pydelta = load_legendre(mode,is_lux=lux)['pdelta']
            py0_taylor = taylor_y0_coeffs( mode, fpmode = fpmode )
            py0 = load_legendre(mode,is_lux=lux)['p0']

            def add_polydef( resvarname, varname, polycoef ):
                kw = dict(
                    polycoeffs = polycoef,
                    resvarname = resvarname,
                    varname = varname,
                    fmtfct = lambda ci, i : '%.15g'%ci,
                    lang = lang)

                l1 = fmt_horner(**kw)
                if len(l1)==1:
                    return l1
                #needs multiline parens:
                return fmt_horner(**kw, enclose_in_parens = True)

            highx_pow = recipe_highx_pow(mode)
            if lang=='py':
                o.append(f'def {fctname}( x, sintheta ):' )
                o.append( '    if x < %s:'%taylorthr)
                o += add_polydef('y0','x',py0_taylor)
                o.append( '    else:')
                o.append( '        if x > 1e3:')
                if mode == 'scndgauss':
                    assert 0.92 < highx_pow < 0.94
                    o.append(f'            return {fctname}(1e3,sintheta)*(x*1e-3)**-{highx_pow:.13g}')
                else:
                    assert highx_pow == 0.5
                    o.append(f'            return {fctname}(1e3,sintheta)*(x*1e-3)**-0.5')
                o.append( '        xs = x**0.5')
                o.append( '        xp = (xs-1.0)/(xs+1.0)')
                o += add_polydef('y0','xp',py0)
                o.append( '    s = sintheta**0.5')
                o.append( '    u = x*s')
                o.append( '    if u < %s:'%taylorthr)
                o += add_polydef('ydelta','u',pydelta_taylor)
                o.append( '    else:')
                o.append( '        us = u**0.5')
                o.append( '        up = (us-1.0)/(us+1.0)')
                o += add_polydef('ydelta','up',pydelta)
                o.append( '    return y0 + sintheta*s*ydelta')
            elif lang in ('c','cpp'):
                o.append(f'double {fctname}( double x, double sintheta )')
                o.append( '{')
                o.append( '  double y0;')
                o.append( '  if ( x < %s ) {'%taylorthr)
                o += add_polydef('y0','x',py0_taylor)
                o.append( '  } else {')
                o.append( '    if ( x > 1e3 )')
                if mode == 'scndgauss':
                    assert 0.92 < highx_pow < 0.94
                    o.append(f'      return {fctname}(1e3,sintheta)*{ns}pow(x*1e-3,-{highx_pow:.13g});')
                else:
                    assert highx_pow == 0.5
                    o.append(f'      return {fctname}(1e3,sintheta)*{ns}sqrt(1e3/x);')
                o.append(f'    const double xs = {ns}sqrt(x);')
                o.append( '    const double xp = (xs-1.0)/(xs+1.0);')
                o += add_polydef('y0','xp',py0)
                o.append( '  }')
                o.append(f'  const double s = {ns}sqrt(sintheta);')
                o.append( '  const double u = x*s;')
                o.append( '  double ydelta;')
                o.append( '  if ( u < %s ) {'%taylorthr)
                o += add_polydef('ydelta','u',pydelta_taylor)
                o.append( '  } else {')
                o.append(f'    const double us = {ns}sqrt(u);')
                o.append( '    const double up = (us-1.0)/(us+1.0);')
                #o += do_wrap_assign('ydelta',fmt_ydelta())
                o += add_polydef('ydelta','up',pydelta)
                o.append( '  }')
                o.append( '  return y0 + sintheta*s*ydelta;')
                o.append('}')
            else:
                assert lang=='latex'
                ml = mode2letter(mode).upper()
                o.append(r'\begin{algorithm}')
                o.append(r'    \caption{Algorithm for calculating $y_%s(x,\sin\theta)$.'%ml)
                o.append(r'             When $x\le1000$ the result has an absolute error less than $10^{-%i}\cdot\min(y_\text{true},1-y_\text{true})$.'%ndigits_recipe)
                o.append(r"             All variables must be double precision floating point types or better, as per the \texttt{binary64} format\protect\cite{IEEE754}.}")
                o.append(r'    \begin{algorithmic}[1]')
                o.append(r'        \small')
                fctname = 'Y%s%s'%(ml,'Lux' if lux else '')
                o.append(r'        \Procedure{%s}{$x, \sin\theta$}'%fctname)
                o.append(r'            \If{$x < %s$}'%taylorthr)
                o += fmt_horner_for_latex(py0_taylor,'x',resvarname='y_0')
                o.append(r'            \Else')
                o.append(r'                \If{$x>1000$}')
                if mode == 'scndgauss':
                    assert 0.92 < highx_pow < 0.94
                    tmppow = r'(0.001{\cdot}x)^{-%.13g}'%highx_pow
                else:
                    assert highx_pow == 0.5
                    tmppow = r'\sqrt{1000/x}'
                o.append(r'                    \State \textbf{return} $\text{%s}(1000,\sin\theta)\cdot%s$'%(fctname,tmppow))
                o.append(r'                \EndIf')
                o.append(r'                \State $x_s \gets \sqrt{x} $')
                o.append(r'                \State $x%s \gets (x_s-1)/(x_s+1) $'%"'")
                o += fmt_horner_for_latex(py0,"x'",resvarname='y_0')
                o.append(r'            \EndIf')
                o.append(r'                \State $s \gets \sqrt{\sin\theta} $')
                o.append(r'                \State $u \gets x{\cdot}s$')
                o.append(r'            \If{$u < %s$}'%taylorthr)
                o += fmt_horner_for_latex(pydelta_taylor,'u',resvarname=r'y^\Delta')
                o.append(r'            \Else')
                o.append(r'                \State $u_s \gets \sqrt{u} $')
                o.append(r'                \State $u%s \gets (u_s-1)/(u_s+1) $'%"'")
                o += fmt_horner_for_latex(pydelta,"u'",resvarname=r'y^\Delta')
                o.append(r'            \EndIf')
                o.append(r'            \State \textbf{return} $y_0 + \sin\theta{\cdot}s{\cdot}y^\Delta$')
                o.append(r'        \EndProcedure')
                o.append(r'    \end{algorithmic}')
                o.append(r'    \label{alg:yfct_%s}'%(('%s_lux' if lux else '%s')%mode))
                o.append(r'\end{algorithm}')
                flush_if_latex('%s%s'%(mode, '_lux' if lux else ''))
            if not flush_if_latex_mode:
                o.append('')

    if include_testfcts:

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
                o.append(f'void bc2015_test_impl_{mode}( double x, double sinth, double refval)')
                o.append( '{')
                o.append(f'  const double val = {fctname}(x,sinth);')
                o.append(f'  const double val_lux = {fctname}_lux(x,sinth);')
                o.append(f'  if (!(fabs(val-refval) <= {eps_std}*fmin(refval,1.0-refval))) '+'{')
                o.append( '    printf("bc2015_test_implementation: "')
                o.append(f'           "Failure in {fctname}(%g,%g)\\n",x,sinth);')
                o.append( '    exit(1);')
                o.append( '  }')
                o.append(f'  if (!(fabs(val_lux-refval) <= {eps_lux}*fmin(refval,1.0-refval))) '+'{')
                o.append( '    printf("bc2015_test_implementation: "')
                o.append(f'           "Failure in {fctname_lux}(%g,%g)\\n",x,sinth);')
                o.append( '    exit(1);')
                o.append( '  }')
                o.append( '}')
                o.append( '')

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
                o.append( '      std::ostringstream ss;')
                o.append( '      ss << "bc2015_test_implementation: "')
                o.append(f'         << "Failure in {fctname}("<<x<<","<<sinth<<")\\n";')
                o.append( '      throw std::runtime_error(ss.str());')
                o.append( '    }')
                o.append(f'    if (!({ns}fabs(val_lux-refval) <= {eps_lux}*{ns}fmin(refval,1.0-refval))) '+'{')
                o.append( '      std::ostringstream ss;')
                o.append( '      ss << "bc2015_test_implementation: "')
                o.append(f'         << "Failure in {fctname_lux}("<<x<<","<<sinth<<")\\n";')
                o.append( '      throw std::runtime_error(ss.str());')
                o.append( '    }')
                o.append( '  };')
            for mode in modes:
                fctname = f'bc2025_y_{mode}'
                fctname_lux = f'{fctname}_lux'
                o.append('')
                o.append(f'  std::cout << "Testing {fctname} + {fctname_lux}"')
                o.append( '            << std::endl;')
                for x,sinth,y in refdata[mode]:
                    o.append( '  test_%s(%.14g,%.14g,%.18g);'%(mode,x,sinth,y))
            o.append('  std::cout<<"All tests completed"<<std::endl;')

        if lang=='c':
            for imode,mode in enumerate(modes):
                fctname = f'bc2025_y_{mode}'
                fctname_lux = f'{fctname}_lux'
                if imode:
                    o.append('')
                o.append(f'  printf("Testing {fctname} + {fctname_lux}\\n");')
                for x,sinth,y in refdata[mode]:
                    o.append(f'  bc2015_test_impl_{mode}(%.14g,%.14g,%.18g);'%(x,sinth,y))
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
                o.append( '            raise RuntimeError(')
                o.append( '                "bc2015_test_implementation: "')
                o.append(f'                "Failure in {fctname}(%g,%g)"%(x,sinth)')
                o.append( '            )')
                o.append(f'        if not (abs(val_lux-refval) <= {eps_lux}*min(refval,1.0-refval)):')
                o.append( '            raise RuntimeError(')
                o.append( '                "bc2015_test_implementation: "')
                o.append(f'                "Failure in {fctname_lux}(%g,%g)"%(x,sinth)')
                o.append( '            )')
                o.append(f'    print("Testing {fctname} + {fctname_lux}")')
                for x,sinth,y in refdata[mode]:
                    o.append( '    t(%.14g,%.14g,%.18g)'%(x,sinth,y))

        if lang=='py':
            o.append('    print("All tests completed")')

        if lang=='cpp':
            o.append('}')

    if include_main:
        for e in textbox(lang,'Hook for executing file as test application'):
            o.append(e)
        #o.append('')
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

    content = '\n'.join(o) + '\n'
    if outfile:
        if not flush_if_latex_mode:
            pathlib.Path(outfile).write_text(content)
    else:
        print(content,end='')
    return o

def _wrap_on_plus( s, width = 80, indent='    '):
    if len(s)<=width:
        return [s]
    tmp = ''
    out = []
    for i,p in enumerate(s.split('+')):
        if i:
            p = '+'+p
        if len(tmp)>len(indent) and len(tmp)+len(p)>width:
            #flush tmp:
            out += [tmp]
            tmp = indent
        tmp += p
    if tmp:
        out+=[tmp]
    return out

def fmt_horner( polycoeffs,
                resvarname,
                varname,
                fmtfct,
                lang,
                enclose_in_parens = False ):
    #https://en.wikipedia.org/wiki/Horner%27s_method
    aggressive = False
    def fmtfct_19g( ci, i ):
        return '%.19g'%ci
    if fmtfct is None:
        fmtfct = fmtfct_19g

    n = len(polycoeffs)
    assert n>=3
    nleadzeroes = 0
    while nleadzeroes < n and polycoeffs[nleadzeroes]==0:
        nleadzeroes += 1
    polycoeffs = polycoeffs[nleadzeroes:]
    n -= nleadzeroes
    current_sign = [1]
    def c(i):
        e = polycoeffs[i]*current_sign[0]
        e = e if isinstance(e,str) else fmtfct(e,i)
        res = (e+'.0' if e.isdigit() or (e[0]=='-' and e[1:].isdigit()) else e)
        if aggressive:
            if res.endswith('.0'):
                res = res[:-1]
            elif res.startswith('0.'):
                res = res[1:]
        return res
    o = []
    def c_is_neg(i):
        return polycoeffs[i]*current_sign[0] < 0.0
    def flip_signs():
        current_sign[0] *= -1

    i0 = 1
    nfinal_paren = n - 2
    if nleadzeroes:
        s = (f'{varname}*')*nleadzeroes + '(%s'%c(0)
        nfinal_paren = n + 1 - nleadzeroes
    else:
        s = c(0)
    if enclose_in_parens:
        s = '( %s'%s

    for i in range(i0,n):
        ci_is_neg = c_is_neg(i)
        if ci_is_neg:
            flip_signs()
        pm = '-' if ci_is_neg else '+'
        ci = c(i)
        assert not ci.startswith('-') and not ci.startswith('+')
        if i+1==n:
            tmp = '%s%s*%s%s'%(pm,ci,varname,')'*nfinal_paren)
            if enclose_in_parens:
                tmp += ' )'
            if lang in ('c','cpp'):
                tmp += ';'
        else:
            tmp = '%s%s*(%s'%(pm,varname,ci)
        w = 61 if o else 71-len(resvarname)-(2 if enclose_in_parens else 0)
        if len(s+tmp)>w:
            o.append(s)
            s = tmp
        else:
            s += tmp
    if s:
        o.append(s)
    tabwidth = 2 if lang in ('c','cpp') else 4
    o[0] = ' '*(2*tabwidth)+resvarname + ' = '+o[0]
    for i in range(1,len(o)):
        o[i] = ' '*(5 + 2*tabwidth + len(resvarname)) + o[i]
    return o

def fmt_horner_for_latex( polycoeffs, varname, resvarname ):
    n = len(polycoeffs)
    assert n>=3
    current_sign = [1]
    def c(i):
        e=polycoeffs[i]*current_sign[0]
        if float(round(e))==float(e):
            return str(int(round(e)))
        e = e if isinstance(e,str) else '%.15g'%e
        assert 'e' not in e.lower(), 'exponential notation not handled for latex yet'
        return e

    #While composing, using 'x' as placeholder for the variable and '*' for
    #multiplication, to be replaced at the end.
    nfinal_paren = n - 2
    if c(0)=='0' and c(1)=='0':
        s = 'x*x*(%s'%c(2)
        i0 = 3
    else:
        s = c(0)
        i0 = 1

    o = []
    def c_is_neg(i):
        return polycoeffs[i]*current_sign[0] < 0.0
    def flip_signs():
        current_sign[0] *= -1
    for i in range(i0,n):
        ci_is_neg = c_is_neg(i)
        if ci_is_neg:
            flip_signs()
        pm = '-' if ci_is_neg else '+'
        ci = c(i)
        assert not ci.startswith('-') and not ci.startswith('+')
        if i+1==n:
            tmp = '%s%s*x%s'%(pm,ci,')'*nfinal_paren)
        else:
            tmp = '%sx*(%s'%(pm,ci)
        w = 60
        if o:
            w += 2
        if len(s+tmp)>w:
            o.append(s)
            s = tmp
        else:
            s += tmp
    if s:
        o.append(s)
    o[0] = r'                \State ${%s} \gets '%resvarname + o[0] + '$'
    for i in range(1,len(o)):
        o[i] = r'                  \hspace{6em}$' + o[i]+'$'
    for i in range(0,len(o)-1):
        o[i] += r'\\'

    for i in range(len(o)):
        o[i] = o[i].replace('*',r'{\cdot}').replace('x','{%s}'%varname)
    return o

def textbox(lang,*msg_lines):
    if lang!='latex':
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
    elif lang=='latex':
        pass
    else:
        assert False
    if lang!='latex':
        yield ''

def mode2letter(mode):
    return { 'primary':'P',
             'scndgauss':'G',
             'scndlorentz':'L',
             'scndfresnel':'F'}[mode]
