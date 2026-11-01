* Spot check of the fasttrop covariate results against the Stata package
* gorgeousfish/TROP. This file only calls the Stata command. It does not
* include code from that package.
*
* Usage: start Stata in the fasttrop source folder, then type
*   do inst/validation/stata/trop_covariates.do
* The file writes inst/validation/stata/stata-covariates.csv. Then run
*   Rscript inst/validation/stata/compare_stata.R
*
* The panels are tests/testthat/fixtures/XC1.csv, XC2.csv, and XC3.csv. The
* tuning parameters are fixed. They are the same as in
* tests/testthat/fixtures/lambdas.csv. Stata uses 1e10 for lambda_nn = Inf.
* Each fit is made without covariates and with covariates. Thus, the compare
* script can separate a difference in the covariate step from a difference in
* the estimator without covariates.

version 16
clear all
set more off

capture which trop
if _rc {
    net install trop, from("https://raw.githubusercontent.com/gorgeousfish/TROP/main") replace
}
trop, version

local fixdir "tests/testthat/fixtures"
local outdir "inst/validation/stata"

tempname post
tempfile results
postfile `post' str8 panel str16 config str8 method byte covariates ///
    double(att gamma1 gamma2) byte converged int rc using `results', replace

local configs  "twfe weighted_twfe mc trop_a trop_b"
local l_twfe          "0 0 1e10"
local l_weighted_twfe "0.5 0.5 1e10"
local l_mc            "0 0 0.1"
local l_trop_a        "0.3 0.8 0.05"
local l_trop_b        "1 0.2 1"

foreach p in XC1 XC2 XC3 {
    foreach c of local configs {
        foreach m in local global {
            foreach x in 0 1 {
                import delimited using "`fixdir'/`p'.csv", clear
                local copt ""
                if `x' local copt "covariates(x1 x2)"
                capture noisily trop y d, panelvar(unit) timevar(time) method(`m') ///
                    fixedlambda(`l_`c'') bootstrap(0) tol(1e-10) maxiter(100000) `copt'
                local rc = _rc
                local att = .
                local g1 = .
                local g2 = .
                local conv = .
                if `rc' == 0 {
                    local att = e(att)
                    local conv = e(converged)
                    if `x' {
                        tempname G
                        matrix `G' = e(gamma)
                        local g1 = `G'[1, 1]
                        if colsof(`G') >= 2 local g2 = `G'[1, 2]
                    }
                }
                post `post' ("`p'") ("`c'") ("`m'") (`x') (`att') (`g1') (`g2') ///
                    (`conv') (`rc')
            }
        }
    }
}
postclose `post'

use `results', clear
export delimited using "`outdir'/stata-covariates.csv", replace
display as text "Done. The results are in `outdir'/stata-covariates.csv"
