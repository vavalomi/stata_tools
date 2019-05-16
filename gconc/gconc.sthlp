
{smcl}
{* *! version 1.0  11jan2010}{...}
{cmd:help gconc}
{hline}

{title:Title}

{p2colset 5 14 22 2}{...}
{p2col :{cmd:gconc} {hline 2}}Generalized measures of concentraction{p_end}
{p2colreset}{...}


{title:Syntax}

{p 8 16 2}
{opt gconc} {it:incomevar} [{it:outcomevar}] {ifin} {weight}
  [{cmd:, nu(}{it:#}{cmd:) keepnegatives nozeros over(}{varlist}{cmd:)}
   {cmd:vce(}{help vce_option}{cmd:) robust svy} {cmdab:cl:uster(}{varlist}}{cmd:)}
   {cmd:subpop(}{it:subpopspec}{cmd:)}
  ]



{title:Description}

{pstd}
{cmd:gconc} computes generalized measures of inequality and concentration including
Gini, generalized Gini (S-Gini) and concentration indices.
Probability weights (pweights) and importance weights (iweights) are allowed.
If only {it:incomevar} is specified, the S-Gini coefficient
is computed for this variable. If both {it:incomevar} and {it:outcomevar}
are specified, then the coefficient of generalized concentration
of {it:outcomevar} with respect to ranking on {it:incomevar} is computed.

{title:Options}

{synoptset}
{synopthdr}
{synoptline}
{synopt:{cmd:nu(}{it:#}{cmd:)}}The shape parameter for generalized concentration coefficient.
                  The default value is 2 corresponding to the ranks of {it:incomevar}
                  entering linearly.

{synopt:{cmd:keepnegatives}}Specifies whether the negative values of {it:incomevar} are allowed.

{synopt:{cmd:nozeros}}Specifies whether the zero values of {it:incomevar} are allowed.
                     By default, they are used in estimation.

{synopt:{cmd:over(}{varlist}{cmd:)}}Requests estimation for separate subgroups in the sample.
                  Note that generalized measures of concentration are not decomposable
                  into subgroups in the strict technical sense: the total concentration is
                  not equal to the sum of subgroup concentrations plus the between group
                  concentration.

{synopt:{cmd:vce(}{help vce_option}{cmd:)}}Specifies the method to compute standard errors. See {help vce_option}.

{synopt:{cmd:robust}}Requests computation of the standard errors robust to heteroskedasticity in {it:incomvar}.

{synopt:{cmdab:cl:uster(}{varlist}{cmd:)}}Requests computation of the standard errors robust to
                  intraclass correlation due to {it:varlist}.

{synopt:{cmd:svy}}Requests estimation that respects complex sampling designs. The resampling
                  variance estimators {cmd:svy, vce(jackknife): gconc ...} and
                  {cmd:svy, vce(brr): gconc ...} will work without this option.
                  The linearized standard errors, {cmd:svy, vce(linearized)},
                  are not directly supported, but can be obtained with this option.

{synopt:{cmd:subpop(}{it:subpopspec}{cmd:)}}For complex survey samples, {cmd:subpop} option
                  should be specified to subset the data, rather than {cmd:if} condition.
                  See {manlink SVY subpopulation estimation}.


{title:Notes and examples}

{pstd}To compute Gini coefficient and its standard error:{p_end}

{phang2}{cmd:. sysuse nlsw88 }{p_end}
{phang2}{cmd:. gconc wage}{p_end}

{pstd}To compute Gini coefficient with the bootstrap standard error:{p_end}

{phang2}{cmd:. sysuse nlsw88 }{p_end}
{phang2}{cmd:. bootstrap, reps(200): gconc wage}{p_end}

{pstd}To compute concentration coefficient for complex survey data:{p_end}

{phang2}{cmd:. webuse nhanes2}{p_end}
{phang2}{cmd:. gconc height weight, svy}{p_end}

{pstd}To compute concentration coefficient with the jackknife standard errors:{p_end}

{phang2}{cmd:. webuse nhanes2}{p_end}
{phang2}{cmd:. svy, vce(jackknife): gconc height weight}{p_end}

{pstd}
Extensions of the Gini index that allow for different sensitivities in different
parts of distribution were given by Yitzhaki (1983).
Computations by {cmd:gconc} utilize the representations of the Gini
method family of concentration measures given by Yitzhaki (1991).
Sandstrom, Wretman and Walden (1988) report that linearization variance estimator
of the Gini coefficient was biased by a factor of 10 in their simulations, and
recommend using jackknife estimator. Jackknife is also the preferred estimatior
in Yitzhaki (1991). However results in Barrett and Donald (2009) showed quite accurate
performance of the linearization estimator.


{title:References}

{phang}
Barrett, G. F., and Donald, S. G. (2009).
Statistical Inference with Generalized Gini Indices of Inequality, Poverty, and Welfare.
{it:Journal of Business and Economic Statistics}, {cmd:27} (1), 1-17,
{browse "http://dx.doi.org/10.1198/jbes.2009.0001":doi:10.1198/jbes.2009.0001}.

{phang}
Sandstrom, A., Wretman, J. H., and Walden, B. (1988). Variance Estimators
of the Gini Coefficient: Probability Sampling.
{it:Journal of Business and Economic Statistics},
{cmd:6} (1), 113--119, {browse "http://www.jstor.org/stable/1391424"}.

{phang}Yitzhaki, S. (1983). On an extension of the Gini inequality index.
{it:International Economic Review}, {cmd:24} (3), 617--628,
{browse "http://dx.doi.org/10.2307/2648789":doi:10.2307/2648789}.

{phang}Yitzhaki, S. (1991).
Calculating jackknife variance estimators for parameters of the gini method.
{it:Journal of Business and Economic Statistics}, {cmd:9} (2), 235--239,
{browse "http://www.jstor.org/stable/1391792"}.


{title:Also see}

{psee}
Help:  {help inequality}, {help glcurve} (if installed).


{title:Authors}

{phang}
Stanislav Kolenikov (skolenik at gmail dot com)

{phang}
Zurab Sajaya (zsajaia at worldbank dot org)
