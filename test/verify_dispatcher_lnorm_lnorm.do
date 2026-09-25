version 17
clear all
set more off
set seed 42
adopath ++ "/Users/nfukuda/lillardhaz-stata/pkg"

* ================================================================
* Full lillardhaz_d0 dispatcher test: lognormal x lognormal,
* correlated. Same DGP as verify_lnorm_lnorm.do (which used a
* standalone lf program) but now routed through the actual
* dispatcher used by the lillardhaz.ado front end.
* ================================================================
global LH_eq1type "lognormal"
global LH_eq2type "lognormal"
global LH_corr = 1

quietly {
    set obs 20000
    gen double x1 = rnormal()
    gen double x2 = rbinomial(1, 0.5)

    scalar b10=1.0
    scalar b11=-0.3
    scalar lnsig1 = -0.2
    scalar b20=0.6
    scalar b21=0.4
    scalar lnsig2 = 0.1
    scalar rho_true = 0.5

    gen double mu1 = b10 + b11*x1
    gen double mu2 = b20 + b21*x2
    gen double sig1 = exp(lnsig1)
    gen double sig2 = exp(lnsig2)

    gen double u1 = rnormal()
    gen double u2 = rho_true*u1 + sqrt(1-rho_true^2)*rnormal()

    gen double t1 = exp(mu1 + sig1*u1)
    gen double t2 = exp(mu2 + sig2*u2)

    gen double c1 = -ln(runiform())/0.03
    gen double c2 = -ln(runiform())/0.03

    gen double time1  = min(t1, c1)
    gen byte   event1 = (t1 <= c1)
    gen double time2  = min(t2, c2)
    gen byte   event2 = (t2 <= c2)
}

display as result _n "TRUE: b10=" b10 " b11=" b11 " ln_sigma1=" lnsig1 " b20=" b20 " b21=" b21 " ln_sigma2=" lnsig2 " rho=" rho_true

ml model d0 lillardhaz_d0 (eq1: time1 event1 time2 event2 = x1) (ln_sigma1:) (eq2: x2) (ln_sigma2:) (atanh_rho:), technique(bfgs)
ml maximize, iterate(200)

display as result _n "atanh_rho -> rho = " tanh(_b[atanh_rho:_cons])
