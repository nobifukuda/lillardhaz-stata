version 17
clear all
set more off
set seed 31
adopath ++ "/Users/nfukuda/lillardhaz-stata/pkg"

* ================================================================
* Full lillardhaz_d0 dispatcher test: probit x piecewise Gompertz,
* correlated. Same DGP as verify_probit_pgomp.do (standalone lf
* program) but now routed through the actual dispatcher (which
* calls _lillardhaz_pgomp for eq2). 3-depvar layout (y1=binary,
* y2=time2, y3=event2) -- the case already correct before the
* y2time fix; this run confirms no regression.
* ================================================================
global LH_eq1type "probit"
global LH_eq2type "pgompertz"
global LH_nodes2 "5"
global LH_corr = 1

quietly {
    set obs 20000
    gen double x1 = rnormal()
    scalar g0=0.2
    scalar g1=-0.4
    scalar c0=-1.1
    scalar gamma1=0.3
    scalar gamma2=-0.12
    scalar rho_true=0.3
    scalar node1 = 5

    local L2a = c0
    local L2b = `L2a'+gamma1*node1

    gen double e1 = rnormal()
    gen double e2 = rho_true*e1 + sqrt(1-rho_true^2)*rnormal()
    gen byte y1 = (g0+g1*x1+e1 > 0)
    gen double U2 = normal(e2)
    gen double Htarget2 = -ln(1-U2)
    scalar seg2a_true = exp(`L2a')/gamma1*(exp(gamma1*node1)-1)
    gen double t2 = cond(Htarget2 < seg2a_true, ln(1+Htarget2*gamma1/exp(`L2a'))/gamma1, ///
        node1 + ln(1+(Htarget2-seg2a_true)*gamma2/exp(`L2b'))/gamma2)
    gen double c2 = -ln(runiform())/0.03
    gen double time2 = min(t2,c2)
    gen byte event2 = (t2<=c2)
    replace time2 = 1e-4 if time2 < 1e-4
}

display as result _n "TRUE: g0=" g0 " g1=" g1 " level2=" c0 " slope2a=" gamma1 " slope2b=" gamma2 " rho=" rho_true

ml model d0 lillardhaz_d0 (eq1: y1 time2 event2 x1 = x1) (level2:) (slope2a:) (slope2b:) (atanh_rho:), technique(bfgs)
ml maximize, iterate(250)
display as result _n "rho estimate = " tanh(_b[atanh_rho:_cons])
