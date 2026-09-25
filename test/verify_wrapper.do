version 17
clear all
set more off
set seed 55
adopath ++ "/Users/nfukuda/lillardhaz-stata/pkg"

* ================================================================
* End-to-end test of the lillardhaz.ado front-end command (not the
* raw lillardhaz_d0 evaluator): builds the equation bracket, calls
* ml model/ml maximize, and exercises lillardhaz_p (predict) --
* pgompertz x pgompertz, correlated, with x1/x2 covariates.
* ================================================================
quietly {
    set obs 25000
    gen double x1 = rnormal()
    gen double x2 = rbinomial(1, 0.5)

    scalar b0=-1.3
    scalar b1=0.3
    scalar alpha0=-0.4
    scalar alpha1=0.35
    scalar alpha2=-0.1
    scalar lev2_0=-0.9
    scalar lev2_1=-0.25
    scalar gamma0=-0.1
    scalar gamma1=0.25
    scalar gamma2=-0.2
    scalar rho_true = 0.35

    gen double L1a = b0 + b1*x1 + alpha0
    gen double L1b_off = alpha1*5
    gen double L2a = lev2_0 + lev2_1*x2 + gamma0
    gen double L2b_off = gamma1*4

    gen double e1 = rnormal()
    gen double e2 = rho_true*e1 + sqrt(1-rho_true^2)*rnormal()
    gen double U1 = normal(e1)
    gen double U2 = normal(e2)
    gen double Htarget1 = -ln(1-U1)
    gen double Htarget2 = -ln(1-U2)

    gen double seg1a_true = exp(L1a)/alpha1*(exp(alpha1*5)-1)
    gen double t1 = cond(Htarget1 < seg1a_true, ln(1+Htarget1*alpha1/exp(L1a))/alpha1, ///
        5 + ln(1+(Htarget1-seg1a_true)*alpha2/exp(L1a+L1b_off))/alpha2)

    gen double seg2a_true = exp(L2a)/gamma1*(exp(gamma1*4)-1)
    gen double t2 = cond(Htarget2 < seg2a_true, ln(1+Htarget2*gamma1/exp(L2a))/gamma1, ///
        4 + ln(1+(Htarget2-seg2a_true)*gamma2/exp(L2a+L2b_off))/gamma2)

    gen double c1_ = -ln(runiform())/0.03
    gen double c2_ = -ln(runiform())/0.03
    gen double time1 = min(t1,c1_)
    gen byte event1 = (t1<=c1_)
    gen double time2 = min(t2,c2_)
    gen byte event2 = (t2<=c2_)
    replace time1 = 1e-4 if time1<1e-4
    replace time2 = 1e-4 if time2<1e-4
}

display as result _n "TRUE: level1=" b0+alpha0 " b1=" b1 " slope1a=" alpha1 " slope1b=" alpha2
display as result "      level2=" lev2_0+gamma0 " lev2_1=" lev2_1 " slope2a=" gamma1 " slope2b=" gamma2 " rho=" rho_true

lillardhaz time1 event1 time2 event2, eq1(pgompertz) eq2(pgompertz) ///
    x1(x1) x2(x2) nodes1(5) nodes2(4) technique(bfgs) iterate(250)

display as result _n "e(rho) = " e(rho) "  e(eq1type)=" e(eq1type) "  e(eq2type)=" e(eq2type)

predict double s1hat, surv1
predict double s2hat, surv2
predict double f1hat, dens1
predict double f2hat, dens2
predict double xb1hat, xb1
predict double xb2hat, xb2

summarize s1hat s2hat f1hat f2hat xb1hat xb2hat

display as result _n "===== probit x lognormal, correlated, via wrapper ====="
quietly {
    drop _all
    set obs 20000
    gen double x1 = rnormal()
    gen double x2 = rbinomial(1, 0.5)

    scalar g0 = 0.3
    scalar g1 = -0.5
    scalar b20 = 0.7
    scalar b21 = 0.35
    scalar lnsig2 = -0.15
    scalar rho_true2 = -0.4

    gen double w1idx = g0 + g1*x1
    gen double mu2 = b20 + b21*x2
    gen double sig2 = exp(lnsig2)

    gen double e1 = rnormal()
    gen double u2 = rho_true2*e1 + sqrt(1-rho_true2^2)*rnormal()

    gen byte y1 = (w1idx + e1 > 0)
    gen double t2 = exp(mu2 + sig2*u2)
    gen double c2 = -ln(runiform())/0.03
    gen double time2  = min(t2, c2)
    gen byte   event2 = (t2 <= c2)
}
display as result _n "TRUE: g0=" g0 " g1=" g1 " b20=" b20 " b21=" b21 " ln_sigma2=" lnsig2 " rho=" rho_true2

lillardhaz y1 time2 event2, eq1(probit) eq2(lognormal) x1(x1) x2(x2) technique(bfgs)
display as result _n "e(rho) = " e(rho)
predict double pr1hat, pr1
predict double s2hat_b, surv2
summarize pr1hat s2hat_b

display as result _n "===== lognormal x lognormal, NOCORR, via wrapper ====="
quietly {
    drop _all
    set obs 20000
    gen double x1 = rnormal()
    gen double x2 = rbinomial(1, 0.5)

    scalar b10=1.0
    scalar b11=-0.3
    scalar lnsig1 = -0.2
    scalar b20c=0.6
    scalar b21c=0.4
    scalar lnsig2c = 0.1

    gen double mu1 = b10 + b11*x1
    gen double mu2 = b20c + b21c*x2
    gen double t1 = exp(mu1 + exp(lnsig1)*rnormal())
    gen double t2 = exp(mu2 + exp(lnsig2c)*rnormal())
    gen double c1 = -ln(runiform())/0.03
    gen double c2 = -ln(runiform())/0.03
    gen double time1  = min(t1, c1)
    gen byte   event1 = (t1 <= c1)
    gen double time2  = min(t2, c2)
    gen byte   event2 = (t2 <= c2)
}
lillardhaz time1 event1 time2 event2, eq1(lognormal) eq2(lognormal) x1(x1) x2(x2) nocorr technique(bfgs)
display as result "e(corr) = " e(corr)
