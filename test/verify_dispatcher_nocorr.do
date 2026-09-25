version 17
clear all
set more off
set seed 88
adopath ++ "/Users/nfukuda/lillardhaz-stata/pkg"

* ================================================================
* Confirms the existing lillardhaz_d0 evaluator, run with
* $LH_corr=0 (rho fixed at 0, no atanh_rho equation in the ml
* model bracket), correctly reduces to the product of the two
* independent marginal likelihoods -- i.e. no separate univariate
* evaluator programs are needed for the "nocorr" case; Lillard's
* own identification result is that independence collapses the
* joint likelihood to two separate fits, and the algebra confirms
* the shared bivariate formulas factor correctly when rho=0.
* lognormal x lognormal, TRUE rho = 0.
* ================================================================
global LH_eq1type "lognormal"
global LH_eq2type "lognormal"
global LH_corr = 0

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

    gen double mu1 = b10 + b11*x1
    gen double mu2 = b20 + b21*x2
    gen double sig1 = exp(lnsig1)
    gen double sig2 = exp(lnsig2)

    gen double t1 = exp(mu1 + sig1*rnormal())
    gen double t2 = exp(mu2 + sig2*rnormal())

    gen double c1 = -ln(runiform())/0.03
    gen double c2 = -ln(runiform())/0.03

    gen double time1  = min(t1, c1)
    gen byte   event1 = (t1 <= c1)
    gen double time2  = min(t2, c2)
    gen byte   event2 = (t2 <= c2)
}

display as result _n "TRUE: b10=" b10 " b11=" b11 " ln_sigma1=" lnsig1 " b20=" b20 " b21=" b21 " ln_sigma2=" lnsig2 " (rho=0, fixed)"

ml model d0 lillardhaz_d0 (eq1: time1 event1 time2 event2 = x1) (ln_sigma1:) (eq2: x2) (ln_sigma2:), technique(bfgs)
ml maximize, iterate(200)

display as result _n "===== compare to two separate univariate streg fits ====="
stset time1, failure(event1)
streg x1, dist(lognormal) time
stset time2, failure(event2)
streg x2, dist(lognormal) time
