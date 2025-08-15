
################################################################################
##                                                                            ##
##  This file is part of NCrystal (see https://mctools.github.io/ncrystal/)   ##
##                                                                            ##
##  Copyright 2015-2025 NCrystal developers                                   ##
##                                                                            ##
##  Licensed under the Apache License, Version 2.0 (the "License");           ##
##  you may not use this file except in compliance with the License.          ##
##  You may obtain a copy of the License at                                   ##
##                                                                            ##
##      http://www.apache.org/licenses/LICENSE-2.0                            ##
##                                                                            ##
##  Unless required by applicable law or agreed to in writing, software       ##
##  distributed under the License is distributed on an "AS IS" BASIS,         ##
##  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  ##
##  See the License for the specific language governing permissions and       ##
##  limitations under the License.                                            ##
##                                                                            ##
################################################################################

#import NCrystalDev as NC
from .mpmath import mp, mpf
import numpy as np

eq36_result_tolerance = 1e-7

class BCYPCalc:

    def __init__( self, theta_deg ):
        from .eval_fofeta import F_of_Eta
        from .eval_phipi import PhiPi
        from .eval_phi0 import Phi0
        self._f_of_eta = F_of_Eta()
        self._phi0 = Phi0()
        self._phipi = PhiPi()

        self._theta = mpf(theta_deg)
        self._sinth = mp.sin(self._theta * mp.pi/mpf(180) )
        self._sqrtsinth = mp.sqrt( self._sinth )
        self._sinth_to_3div2 = self._sinth * self._sqrtsinth
        self._taylor_lowx_fct = _taylor_lowx_eq36_create_fct( self._sinth )

    @property
    def theta_deg( self):
        return self._theta

    @property
    def sintheta( self):
        return self._sinth

    def phi( self, sigma_r ):
        #BC 1974 eq. 34 ( sr = sigma*r )
        sr = mpf(sigma_r)
        sinth = self._sinth
        assert sr >= 0.0
        if not sinth:
            return self._phi0( sr )
        if sinth == mpf(1):
            return self._phipi( sr )
        u = sr * self._sqrtsinth
        return self._phi0( sr ) + self._sinth_to_3div2 * ( self._phipi(u) - self._phi0(u) )

    def integrand_fct_eq36( self, x ):
        x = mpf(x)
        phi = self.phi
        f = self._f_of_eta
        def integrand(eta):
            fval = f(eta)
            return fval * phi( x * fval )
        return integrand

    def _integral_tail( self, x, a ):
        # Integration from a to infinity. Returns (estimate,max_error).
        # This tail estimate is found by noticing that:
        #
        #  (eta^2-eta)/eta^4 < f(eta) < (eta^2+eta+1)/eta^4
        #

        # Using these two simpler functions as lower and upper bound, it is
        # possible to perform the integration over eq. 36 from [a,inf]
        # analytically. For instance the lower value can be found with (here g
        # is the integrand, xx is our "x" and x is our "eta"):

        # sage: assume(x>0,xx>0,u>0,sinth>=-1,sinth<=1)
        # sage: p0=(3/(64*x**3))*(8*(x^2)+4*x*exp(-4*x)-(1-exp(-4*x)))
        # sage: ppi=(3/(4*x**3))*(x^2-x+log(1+2*x)/2)
        # sage: phi=p0(x=x)+sinth**(3/2)*(ppi(x=x*sqrt(sinth))-p0(x=x*sqrt(sinth)))
        # sage: fmin=(x^2-x)/x^4
        # sage: gmin=fmin(x=x)*phi(x=xx*fmin(x=x))
        # sage: amin=(gmin(x=1/u).taylor(u,0,20)(u=1/x)).integrate(x,a,infinity)
        # sage: fmax=(x^2+x+1)/x^4
        # sage: gmax=fmax(x=x)*phi(x=xx*fmax(x=x))
        # sage: amax=(gmax(x=1/u).taylor(u,0,20)(u=1/x)).integrate(x,a,infinity)
        # sage:

        # sage: ((amax+amin)/2)(a=1/t).taylor(t,0,6).collect(t)
        #    6       5 ⎛    2 ⎛     5/2    ⎞        ⎞    3
        #   t ⋅xx   t ⋅⎝4⋅xx ⋅⎝sinth    + 2⎠ - 15⋅xx⎠   t ⋅(3⋅xx - 1)
        # - ───── + ───────────────────────────────── - ───────────── + t
        #     4                    25                         6
        # sage: ((amax-amin)/2)(a=1/t).taylor(t,0,6).collect(t)
        #  6 ⎛    2 ⎛     5/2    ⎞       ⎞      5         4       3    2
        # t ⋅⎝8⋅xx ⋅⎝sinth    + 2⎠ - 5⋅xx⎠   3⋅t ⋅xx   3⋅t ⋅xx   t    t
        # ──────────────────────────────── - ─────── - ─────── + ── + ──
        #                20                    10         4      6    2
        #
        # Notice that we helped sagemath along by temporarily switching variable
        # to u=1/x (1/eta), and performing a Taylor expansion in u (which is
        # small "near infinity"). The final .taylor(t,0,6) is just used to drop
        # higher orders of t=1/a.

        #Needs lower bound on a to ensure convergence:
        if a < x or a < 10.0:
            return None, 1e99

        x = mpf(x)
        x2 = x*x
        t = mpf(1)/mpf(a)
        t2 = t*t
        t3 = t2*t
        t4, t5, t6 = t2*t2, t2*t3, t3*t3

        s52 = self._sinth_to_3div2 * self._sinth
        self._sinth = mp.sin(self._theta * mp.pi/mpf(180) )
        self._sqrtsinth = mp.sqrt( self._sinth )
        self._sinth_to_3div2 = self._sinth * self._sqrtsinth
        val = ( t - t3*(mpf(3)*x - mpf(1))/mpf(6)
                + t5 * ( mpf(4)*x2*(s52 + mpf(2)) - mpf(15)*x ) / mpf(25)
                - t6*x / mpf(4)
               )
        max_err = ( t2/mpf(2) + t3 / mpf(6) - mpf('3/4')*t4*x - mpf('3/10')*t5*x
                    + (mpf(8)*(s52 + mpf(2))*x2 - mpf(5)*x)*t6/mpf(20)
                   )
        #Finally, the missing factor from outside the integral in eq. 36 + a
        #factor of 2 since we are only integrating one of the tails.
        k = ( mpf(6)/ (mpf(4)*mp.pi) )
        val *= k
        max_err *= k
        #Sanity checks:
        if not ( (1e-50 < val < 0.1) and  (1e-50 < max_err < 1e-6) ):
            return val, 1e99
        return val, max_err

    def _calc_yp_taylorx( self, x ):
        # If x is small, we can do a taylor expansion in x and then evaluate the
        # resulting integral analytically, yielding the following formula:
        if x > 1.0:
            #Not applicable
            return None
        if mpf(x)==0:
            return mpf(1), mpf(0)
        _, max_err = self._taylor_lowx_fct( x, only_error = True )
        if max_err > eq36_result_tolerance*0.1:
            return None
        return self._taylor_lowx_fct( x )

    def calc_yp( self, x ):
        #WARNING: This gets very slow if x>1e9!!
        assert x<1e11, "gets very slow at such high x"

        #We integrate the function via the mpmath quad method, with
        #Gauss-Legendre (since we are smooth and have no end-point
        #singularities). One complication is that our f(eta) function contains
        #sin(2*eta) and sin^2(eta)~=cos(2eta) terms, and therefore our integrand
        #contains slight oscilations of period eta=2. Thus, we make sure to
        #integrate in steps of at most 10-20 times that period, with extra care
        #near eta=0, where most of the contribution comes from.

        #Handle very small x values:
        x = mpf(x)
        eps = eq36_result_tolerance
        taylorx = self._calc_yp_taylorx(x)
        if taylorx is not None:
            print(f"  Using taylorx at theta={float(self._theta):g} x={float(x):g}")
            return taylorx

        tail_a = ( mpf(10.0)**(-mp.log(eps)/(mpf(2)*mp.log(10))) )
        tail_a *= mpf(4)
        if tail_a < x:
            tail_a = x
        if tail_a < 10:
            tail_a = mpf(10.0)
        #print("tail_a=",tail_a)

        while True:
            tail_contrib, tail_maxerr = self._integral_tail( x, tail_a )
            if tail_maxerr > 0.1*eps:
                tail_a *= mpf('1.1')
            else:
                break
        #print("tail_a=",tail_a)
        #print('tail_error: %e'%tail_maxerr)
        assert eps > 10**(-(mp.dps+10)), f"eps={eps} requires higher mp.dps"

        bounds = [ mpf(float(e)) for e in np.linspace(0,tail_a,1+round(tail_a/2)) ]
        bounds[0] = mpf(0)
        bounds[-1] = tail_a
        i_fct = self.integrand_fct_eq36( x )

        def calc( bounds ):
            #print("CALC x=",x,"bounds=",len(bounds))
            maxdegree_min = 3
            maxdegree_max = 6
            ok = False
            for maxdegree in range(maxdegree_min,maxdegree_max+1):
                print("mp.quad( .., len(bounds)=%i, maxdegree=%i )"%( len(bounds),
                                                                      maxdegree) )
                res, error = mp.quad( i_fct,
                                      bounds,
                                      method='gauss-legendre',
                                      maxdegree=maxdegree,
                                      error=True )
                error /= res
                #print('error: %e'%error)
                ok = error < 0.1 * eps
                if ok:
                    break
            if not ok:
                print("WARNING: Problems integrating!")
            k = ( mpf(6)/ (mpf(4)*mp.pi) )
            return k * res, k * error


        result, max_error = calc( bounds )
        #We assume the worst-case of error correlation and just add them
        #linearly:
        print( "contrib (quad) = %g"%result)
        print( "contrib (tail) = %g"%tail_contrib)
        print( "max_error (quad) = %g"%max_error)
        print( "max_error (tail) = %g"%tail_maxerr)
        result += tail_contrib
        max_error += tail_maxerr

        assert mpf(-eps) < result < mpf(1)+mpf(eps)
        result = ( mpf(1)
                   if result > mpf(1)
                   else ( mpf(0) if result<mpf(0) else result ) )
        return result, max_error

def _taylor_lowx_eq36_create_fct( sinth ):
    #found by taylor expanding integrand, then integrating.
    sinth = mpf(sinth)
    def s( power ):
        return sinth**mpf(power)
    m = mpf
    c0 = m('2/3')
    c1 = m('-22/35')
    c2 = m('4738/17325')*(s('5/2') + m(2))
    c3 = m('-108240952/273648375')*(m(2)*s(3) + m(1))
    c4 = m('350769113251/2887355220750')*(m(13)*s('7/2') + m(2))
    c5 = m('-276478446363113/4269160933537500')*(m(43)*s(4) + m(2))
    c6 = m('19564394058822418571/1283790057226395468750')*(m(311)*s('9/2') + m(4))
    c7 = m('-36703165412564741196614548/1429514849663909766179296875')*(m(314)*s(5) + m(1))
    c8 = m('36928840127061822463096205687/7547838406225443565426687500000')*(m(2833)*s('11/2') + m(2))
    c9 = m('-33748450139852875393232090287374607/19822052940875302766351098083187500000')*(m(14173)*s(6) + m(2))
    c10 = m('19603807095172946325281028797499429174403/71931091630720835149483462686690932812500000')*(m(155921)*s('13/2') + m(4))
    c11 = m('-4404234687448218327421971987673239988236198193/27258363254423962203857012446476293758552734375000')*(m(467773)*s(7) + m(2))
    c12 = m('271729330015704263744140826096102038543392468179939/12187759378318041980588547405068480465324098593750000000')*(m(6081071)*s('15/2') + m(4))
    c13 = m('-1905969038362018082597633645314669909788695669993864257013/331136140946503929461065911291835084567987451682796875000000000')*(m(42567521)*s(8) + m(4))
    c14 = m('5903507038519501143784298379985297238620823232177510866209599031/8458592330951096425880657546389540873931381509862333763621093750000000')*(m('638512867')*s('17/2') + m(8))
    c15 = m('-86204560112829310767154769236518867366235607759705804898584404601/67524407134503654639728921504399745124869018714237960860804443359375000')*(m('638512874')*s(9) + m(1))
    c16 = m('107725131117286353966210744394213813735305091374757266317751547055561765579/780286680810597625462526165143453528771265344415719867227838046691750000000000000')*(m('10854718873')*s('19/2') + m(2))
    c17 = m('-192093005779411254048810194953839432177316626944299322096503088051203039754862607/6781064316086202228435193119109777788134142587932894136953163381366167781250000000000000')*(m('97692469873')*s(10) + m(2))
    c18 = m('670908254986834689460447914327944885431110936072775648170485974671797992908741805619503811/242692684251902441619302255459530234086751841320264246970676199247886507506710760156250000000000000')*(m('1856156927621')*s('21/2') + m(4))
    c19 = m('-3883055495278357133884255336777672572411560830118646314653682162493090494680949274801539986561033/3774045009634729664276442067293310929462816696122544614468334665095836417237230378582105273437500000000000')*(m('9280784638123')*s(11) + m(2))
    c20 = m('25236193726733596737393854392950088667786986729280271069226782191667243067562027856461321896295162131811/275872671830273162630771835836562395791082781009510122381022311711940114970643056390074204808828125000000000000000')*(m('194896477400621')*s('23/2') + m(4))
    c21 = m('-201539919426668053017335610062723121822123353849294299602444827316408533110530087827827897752046913031400914391/12940746478895665105389357764765738641481545149741055523207015122863233421080567004847024102993317419531250000000000000000')*(m('2143861251406871')*s(12) + m(4))
    c22 = m('2889264471964207176908817976433611491392188244874063997073974100578921872573819298482261808158845535960359869600036957/2272130713027470353147989040071527678708496061561741098228359005595002961935649407634367578489055804555904326171875000000000000000')*(m('49308808782358117')*s('25/2') + m(8))
    c23 = m('-6087701393757930630438121662583163938754484330691513123432393618842768864489972489705926413199627387240106928665975838337/7628459158236821035495130118194298208068428714574385285921933845437983305261585901166084878100783044721008227769775390625000000000000')*(m('147926426347074373')*s(13) + m(2))
    c24 = m('499714821155154724192055635481466587305910177223126680494392952981041348051319844000045227087413383516333063208052480567932829696577/8295526857163382584936304267877645974181447795211823937137090698271612731156144960790018233051685943880693964918409216511406250000000000000000000')*(m('3698160658676859371')*s('27/2') + m(4))
    c25 = m('-30023508819070147056421028397853996081209759326484733842284085459179872766314764749647383390174507097964831907699529399401040176584232425159/3425965220139958272877856949243609353797163574636677235577169445046645650935249997372084024271271329149012123350329907433249839355468750000000000000000000')*(m('48076088562799171871')*s(14) + m(4))
    c26 = m('97471086673754123351303360366393885205969956967007807067901921269028883083266758799395886751768733570630030921350239812197561580387661065687651/158468708583698875818502204753167876115932492377912998880911711441711158075700594801049854397643123104734254633446211075985404140178222656250000000000000000000')*(m('1298054391195577640617')*s('29/2') + m(8))
    c27 = m('-6597568037817919166909905608458713969799886602032341262598901949811589786577119561876883656116963358894695241296915364544115030236750944993898430133781/39547487260802617011427204661827127923176792536630306878360074164592051629668380962681407357507933532578557974553238891413631813211470798178100585937500000000000000000')*(m('9086380738369043484371')*s(15) + m(4))
    c28 = m('144578244898200343405729653171148907555881332292329108434502237934874802816999830381166211362755697489005682966070041466655033861384804153899043456104587567203/13214377355931040974872333709252554531838812381948098847646416812252544738231233782890677695857973106556142855147845469535912617826362158691359452880859375000000000000000000000')*(m('263505041412702261046867')*s('31/2') + m(8))
    c29 = m('-41361643146937523165339505029738300623538449472235010288124044622522408489165973547935084445790262235643862309268797849534137010560133850184044650006346053727313128301/29766949463978269387093459230005906222756578592645196959733628244747684548850866774265477421541889757936595766250000013451584483143930400307372679588876333984375000000000000000000000000')*(m('3952575621190533915703117')*s(16) + m(8))
    c30 = m('21686398646550445719100867314692340193996176456320187545075676437386683908382214774538875230473256608873436679296823238958988250858727723100113089038658629933938891023263919/253584469729013452034764343370165855255957535961315100584095566277198747118134134469044358649112465970823085593511892413701113001814624834193863216415479160912738037109375000000000000000000000')*(m('122529844256906551386796859')*s('33/2') + m(16))
    def taylor( xx, only_error = False ):
        x = xx
        v = None
        if not only_error:
            v = c0+x*(c1+x*(c2+x*(c3+x*(c4+x*(c5+x*(c6+x*(c7+x*(c8+x*(c9+x*(c10+x*(c11+x*(c12+x*(c13+x*(c14+x*(c15+x*(c16+x*(c17+x*(c18+x*(c19+x*(c20+x*(c21+x*(c22+x*(c23+x*(c24+x*(c25+x*(c26+x*(c27+x*(c28+x*(c29+x*c30)))))))))))))))))))))))))))))
            v *= m('3/2')
        #We should in principle use the 31st term to estimate the error, but
        #conservatively we use the 30th:
        max_err = m('3/2') * c30 * x**30
        return v, max_err
    return taylor


def benchmark():
    import time
    res = []
    for theta_deg in 0, 30, 60, 90:
        for x in [ 0, 1e-20, 1e-3, 0.2, 0.5, 1.0, 1.01, 2.0, 3, 30, 100, 1e5, 1e9 ]:
            print()
            print('-'*70)
            print(f' theta={theta_deg:g},  x={x:g}')
            t0 = time.time()
            calc = BCYPCalc( theta_deg=theta_deg )
            v, err = calc.calc_yp( x=x )
            t0 = time.time() - t0
            print('   RESULT: %g'%v)
            print('   1-RESULT: %g'%(mpf(1)-v))
            print('           +- %e'%err)
            print(f'   TIME: {t0:.2g}s')
            res.append( (t0, theta_deg, x, v, err ) )
