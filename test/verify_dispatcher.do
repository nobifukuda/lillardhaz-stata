version 17
clear all
set more off
set seed 55
adopath ++ "/Users/nfukuda/lillardhaz-stata/pkg"

* ================================================================
* Full lillardhaz_d0 dispatcher test: pgompertz x pgompertz, correlated,
* WITH x1/x2 covariates this time (fixing the DGP oversight from the
* subroutine-only test), arbitrary 2-segment nodes on each side.
* ================================================================
global LH_eq1type "pgompertz"
global LH_eq2type "pgompertz"
global LH_nodes1 "5"
global LH_nodes2 "4"
global LH_corr = 1

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
    gen double L1b_off = alpha1*5   /* added to L1a per obs below */
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

display as result _n "TRUE: level1(b0+alpha0)=" b0+alpha0 " b1=" b1 " slope1a=" alpha1 " slope1b=" alpha2
display as result "      level2(lev2_0+gamma0)=" lev2_0+gamma0 " lev2_1=" lev2_1 " slope2a=" gamma1 " slope2b=" gamma2 " rho=" rho_true

capture noisily ml model d0 lillardhaz_d0 (eq1: time1 event1 time2 event2 x1 = x1) (s1a:) (s1b:) (eq2: x2) (s2a:) (s2b:) (atanh_rho:), technique(bfgs)
capture noisily ml maximize, iterate(250)
display as result _n "first attempt rc = " _rc

display as result _n "===== retry with informed starting values ====="
ml model d0 lillardhaz_d0 (eq1: time1 event1 time2 event2 x1 = x1) (s1a:) (s1b:) (eq2: x2) (s2a:) (s2b:) (atanh_rho:), technique(bfgs)
ml init eq1:_cons=-1.7 eq1:x1=0.3 s1a:_cons=0.35 s1b:_cons=-0.1 eq2:_cons=-1.0 eq2:x2=-0.25 s2a:_cons=0.25 s2b:_cons=-0.2 atanh_rho:_cons=0.36
ml maximize, iterate(250)
