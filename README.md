# lillardhaz (Stata)

A Stata 17+ command for **simultaneous-equations hazard/probit models with
a Gaussian-copula correlation**, in the style of Lillard (1993): a pair of
processes — each a binary probit outcome or a continuous-time hazard
duration — linked through a single correlation parameter between their
underlying error terms, estimated jointly by maximum likelihood.

Four model families are supported via `eq1()`/`eq2()`:

| eq1() | eq2() | typical use |
|---|---|---|
| `probit` | `lognormal` | a binary decision correlated with a log-normal AFT duration |
| `lognormal` | `lognormal` | two correlated log-normal AFT durations |
| `probit` | `pgompertz` | a binary decision correlated with a flexible piecewise-Gompertz duration |
| `pgompertz` | `pgompertz` | two correlated piecewise-Gompertz durations |

`pgompertz` is a piecewise-linear-in-time log-hazard (piecewise Gompertz)
with an arbitrary, user-specified number of segments (`nodes1()`/`nodes2()`).
Every combination can also be fit independently with `nocorr` (rho fixed
at 0), which — as derived in the manual — reduces algebraically to two
separate univariate fits, giving a natural nested baseline for a
likelihood-ratio test of correlation.

See [`lillardhaz.sthlp`](lillardhaz.sthlp) (`help lillardhaz` once
installed) for the full syntax reference, and
[`docs/manual.html`](https://htmlpreview.github.io/?https://github.com/nobifukuda/lillardhaz-stata/blob/main/docs/manual.html)
for the model derivation (the joint likelihood for every eq1()×eq2()
combination, the piecewise-Gompertz closed-form hazard, the independence
reduction, and simulation-based verification results).

## Author

**Nobutaka Fukuda**, Tohoku University — <nobutaka.fukuda@tohoku.ac.jp>

## Installation

```stata
net install lillardhaz, from("https://raw.githubusercontent.com/nobifukuda/lillardhaz-stata/main/") replace
help lillardhaz
```

(Once accepted to SSC — see below — installation will simply be `ssc install lillardhaz`.)

## Usage

```stata
* Probit & log-normal hazard, correlated
lillardhaz y1 time2 event2, eq1(probit) eq2(lognormal) x1(x1) x2(x2)

* Piecewise Gompertz & piecewise Gompertz, correlated, 2 interior nodes each side
lillardhaz time1 event1 time2 event2, eq1(pgompertz) eq2(pgompertz) ///
    x1(x1) x2(x2) nodes1(5) nodes2(4)

* Same, but independent (nested test of rho=0)
lillardhaz time1 event1 time2 event2, eq1(pgompertz) eq2(pgompertz) ///
    x1(x1) x2(x2) nodes1(5) nodes2(4) nocorr

* Post-estimation
predict S2_hat            // fitted eq2 survival (default)
predict f2_hat, dens2     // fitted eq2 density
predict xb1_hat, xb1      // eq1 linear index
```

## Verification

Every eq1()×eq2() combination, correlated and independent, was verified by
simulating 20,000–25,000 observations under known parameters and confirming
`lillardhaz` recovers them; see [`test/`](test) for the scripts and their
logs, and the manual's verification section for a summary table. A
depvar-indexing bug caught during this process (equation 2's duration read
from the wrong `$ML_y#` slot when equation 1 is a hazard rather than a
probit) is documented there as well.

## License

MIT — see [LICENSE](LICENSE).

## References

Lillard, L. A. 1993. Simultaneous equations for hazards: Marriage duration
and fertility timing. *Journal of Econometrics* 56(1–2): 189–217.

Waite, L. J., & Lillard, L. A. 1991. Children and marital disruption.
*American Journal of Sociology* 96(4): 930–953.
