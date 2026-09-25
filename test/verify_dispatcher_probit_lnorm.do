version 17
clear all
set more off
set seed 7
adopath ++ "/Users/nfukuda/lillardhaz-stata/pkg"

* ================================================================
* Full lillardhaz_d0 dispatcher test: probit x lognormal hazard,
* correlated. Same DGP as verify_probit_lnorm.do (standalone lf
* program) but now routed through the actual dispatcher. This is
* the 3-depvar layout (y1=binary, y2=time2, y3=event2), which is
* the case that was already correct before the y2time fix -- this
* run confirms the fix did not regress it.
* ================================================================
global LH_eq1type "probit"
global LH_eq2type "lognormal"
global LH_corr = 1

quietly {
    set obs 20000
    gen double x1 = rnormal()
    gen double x2 = rbinomial(1, 0.5)

    scalar g0 = 0.3
    scalar g1 = -0.5
    scalar b20 = 0.7
    scalar b21 = 0.35
    scalar lnsig2 = -0.15
    scalar rho_true = -0.4

    gen double w1idx = g0 + g1*x1
    gen double mu2 = b20 + b21*x2
    gen double sig2 = exp(lnsig2)

    gen double e1 = rnormal()
    gen double u2 = rho_true*e1 + sqrt(1-rho_true^2)*rnormal()

    gen byte y1 = (w1idx + e1 > 0)
    gen double t2 = exp(mu2 + sig2*u2)
    gen double c2 = -ln(runiform())/0.03
    gen double time2  = min(t2, c2)
    gen byte   event2 = (t2 <= c2)
}

display as result _n "TRUE: g0=" g0 " g1=" g1 " b20=" b20 " b21=" b21 " ln_sigma2=" lnsig2 " rho=" rho_true

ml model d0 lillardhaz_d0 (eq1: y1 time2 event2 = x1) (eq2: x2) (ln_sigma2:) (atanh_rho:), technique(bfgs)
ml maximize, iterate(200)

display as result _n "atanh_rho -> rho = " tanh(_b[atanh_rho:_cons])
