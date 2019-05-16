{smcl}
{* *! version 1.0  2Sep2009}{...}
{cmd:help switch_probit}{right: ({browse "http://www.stata-journal.com/article.html?article=st0233":SJ11-3: st0233})}
{hline}

{title:Title}

{p2col 5 22 24 2 :{hi: switch_probit} {hline 2}}Maximum likelihood estimation of endogenous switching probit model{p_end}


{title:Syntax for switch_probit}

{phang}Basic model

{p 8 17 2}
{cmd:switch_probit}
{it:{help depvar:depvar}}
[{varlist}]
{ifin}
{weight}{cmd:,}
{cmdab:sel:ect:(}{it:depvar_s} {it:varlist_s}{cmd:)}
[{it:{help switch_probit##options:options}}]


{phang}Fully specified model

{p 8 17 2}
{cmd:switch_probit}
{cmd:(}{it:{help depvar:depvar1}} [{cmd:=}] {varlist:1}{cmd:)}
{cmd:(}{it:depvar0} [{cmd:=}] {it:varlist0}{cmd:)}
{ifin}
{weight}{cmd:,}
{cmdab:sel:ect:(}{it:depvar_s} [{cmd:=}] {it:varlist_s}{cmd:)}
[{it:{help switch_probit##options:options}}]


{synoptset 26 tabbed}{...}
{marker options}{...}
{synopthdr}
{synoptline}
{p2coldent:* {opt select()}}specify selection equation: dependent and independent variables{p_end}
{synopt :{opt nocon:stant}}suppress constant term{p_end}
{synopt :{opth offset_s(varname)}}include {it:varname} in selection equation with coefficient constrained to 1{p_end}
{synopt :{opth offset1(varname)}}include {it:varname} in regime 1 with coefficient constrained to 1{p_end}
{synopt :{opth offset0(varname)}}include {it:varname} in regime 0 with coefficient constrained to 1{p_end}
{synopt :{cmdab:const:raints(}{it:{help estimation options##constraints():constraints}}{cmd:)}}apply specified linear constraints{p_end}
{synopt :{opt col:linear}}keep collinear variables{p_end}
{synopt :{opt r:obust}}robust estimator of variance{p_end}
{synopt :{opth cl:uster(varname)}}adjust standard errors for intragroup correlation{p_end}
{synopt :{opt l:evel(#)}}set confidence level; default is {cmd:level(95)}{p_end}
{synopt :{opt noskip}}perform likelihood-ratio test{p_end}
{synopt :{it:{help switch_probit##maximize_options:maximize_options}}}control the maximization process; seldom used{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
* {opt select()} is required.{p_end}
{p 4 6 2}
{opt pweight}s, {opt fweight}s, and {opt iweight}s are allowed; see {help weight}.{p_end}


{title:Syntax for predict}

{p 8 16 2}
{cmd:predict}
{dtype}
{it:{help newvar}}
{ifin}
[{cmd:,} {it:statistic} {opt nooff:set}]

{p 8 16 2}
{cmd:predict}
{dtype}
{c -(}{it:stub}{cmd:*}{c |}{it:newvar_sel} {it:newvar_eq1} {it:newvar_eq0} {it:newvar_athrho1} {it:newvar_athrho0}{c )-}
{ifin}
[{cmd:,}
{opt sc:ores}]

{synoptset 11}{...}
{synopthdr :statistic}
{synoptline}
{synopt :{opt p11}}Pr(depvar_s=1, depvar1=1); the default{p_end}
{synopt :{opt p10}}Pr(depvar_s=1, depvar1=0){p_end}
{synopt :{opt p01}}Pr(depvar_s=0, depvar0=1){p_end}
{synopt :{opt p00}}Pr(depvar_s=0, depvar0=0){p_end}
{synopt :{opt psel}}Pr(depvar_s=1){p_end}
{synopt :{opt pcond1}}Pr(depvar1=1 | depvar_s=1){p_end}
{synopt :{opt pcond0}}Pr(depvar0=1 | depvar_s=0){p_end}
{synopt :{opt zb}}linear prediction for selection equation{p_end}
{synopt :{opt xb1}}linear prediction for regime 1{p_end}
{synopt :{opt xb0}}linear prediction for regime 0{p_end}
{synopt :{opt stdpsel}}standard error of the linear prediction for selection equation{p_end}
{synopt :{opt stdp1}}standard error of the linear prediction for regime 1{p_end}
{synopt :{opt stdp0}}standard error of the linear prediction for regime 0{p_end}
{synopt :{opt tt}}treatment effect on the treated{p_end}
{synopt :{opt tu}}treatment effect on the untreated{p_end}
{synopt :{opt te}}treatment effect{p_end}
{synopt :{opt mte}}marginal treatment effect{p_end}
{synoptline}
{p2colreset}{...}


{title:Description}

{pstd}
{opt switch_probit} fits the full maximum likelihood model of binary choice
with binary endogenous regressors.  It is implemented using the {cmd:d2}
evaluator to calculate the overall log likelihood together with its first and
second derivatives.


{title:Options for switch_probit}

{phang}
{opt select(...)} specifies variables in the selection equation. This option is
an integral part of the {cmd:switch_probit} estimation and is required. Both
instruments and exogenous variables must be specified in {it:varlist_s}. If
there are no instrumental variables in the model, the model will be identified
by nonlinearities.

{phang}
{opt noconstant} suppresses the constant terms.

{phang}
{opth offset_s(varname)}, {opt offset1(varname)}, and
{opt offset0(varname)} include variables in each equation with coefficients
constrained to {cmd:1}.

{pmore}
For more information, see 
{helpb estimation options:[R] estimation options}.

{phang}
{opt constraints(constraints)}, {opt collinear}; see
{helpb estimation options:[R] estimation options}.

{phang}
{opt robust} specifies that the Huber/White/sandwich estimator of the
variance is to be used in place of the conventional maximum likelihood variance
estimator.  {opt robust} combined with {opt cluster()} further allows
observations that are not independent within cluster (although they must be
independent between clusters).  If you specify {cmd:pweight}s, then
{cmd:robust} is implied.  See {findalias frrobust}.

{phang}
{opth cluster(varname)} specifies that the observations are
independent across groups (clusters) but not necessarily within groups.
{it:varname} specifies to which group each observation belongs; for example,
{cmd:cluster(personid)} refers to data with repeated observations on
individuals.  Specifying {opt cluster()} affects the estimated standard errors
and variance-covariance matrix of the estimators (VCE) but not the estimated
coefficients.  {cmd:cluster()} can be used with {helpb pweight}s to produce
estimates for unstratified cluster-sampled data.  Specifying {cmd:cluster()}
implies {cmd:robust}.

{phang}
{opt level(#)}; see {helpb estimation options:[R] estimation options}.

{phang}
{opt noskip} specifies that a full maximum likelihood model with only a
constant for the regression equation be fit.  This model is not displayed
but is used as the base model to compute a likelihood-ratio test for the model
test statistic displayed in the estimation header.  By default, the overall
model test statistic is an asymptotically equivalent Wald test that all the
parameters in the regression equation are zero (except the constant).  For
many models, this option can substantially increase estimation time.

{marker maximize_options}{...}
{phang}
{it:maximize_options} control the maximization process; see
{helpb maximize:[R] maximize}.  With the possible exception of {cmd:iterate(0)}
and {cmd:trace}, you should specify these options only if the model is unstable.
The maximization uses the {cmd:difficult} option by default. This option need
not be specified.


{title:Options for predict}

{phang}
{opt p11}, the default, calculates the probability of being treated and having
a positive outcome [Pr({it:depvar_s}=1, {it:depvar1}=1)].

{phang}
{opt p10} calculates the probability of being treated and having a zero
outcome [Pr({it:depvar_s}=1, {it:depvar1}=0)].

{phang}
{opt p01} calculates the probability of not being treated and having a
positive outcome [Pr({it:depvar_s}=0, {it:depvar0}=1)].

{phang}
{opt p00} calculates the probability of not being treated and having a zero
outcome [Pr({it:depvar_s}=0, {it:depvar0}=0)].

{phang}
{opt psel} calculates the probability of being treated (that is, the
probability of being in regime 1).

{phang}
{opt pcond1} calculates the probability of a positive outcome conditional on
being treated [Pr({it:depvar1}=1 | {it:depvar_s}=1)].

{phang}
{opt pcond0} calculates the probability of a positive outcome conditional on
not being treated [Pr({it:depvar0}=1 | {it:depvar_s}=0)].

{phang}
{opt zb} calculates the probit linear prediction for the selection equation.

{phang}
{opt xb1} calculates the probit linear prediction based on the coefficients of
the outcome equation in regime 1.

{phang}
{opt xb0} calculates the probit linear prediction based on the coefficients of
the outcome equation in regime 0.

{phang}
{opt stdpsel} calculates the standard error of the linear prediction of
the selection equation.

{phang}
{opt stdp1} calculates the standard error of the linear prediction of regime 1.

{phang}
{opt stdp0} calculates the standard error of the linear prediction of regime 0.

{phang}
{opt tt} calculates the treatment effect on the treated.

{phang}
{opt tu} calculates the treatment effect on the untreated.

{phang}
{opt te} calculates the treatment effect.

{phang}
{opt mte} calculates the marginal treatment effect.

{phang}
{opt nooffset} is relevant only if you specified {opt offset_s(varname)},
{opt offset1(varname)}, or {opt offset2(varname)} for {cmd:switch_probit}. It
modifies the calculations made by {opt predict} so that they ignore the offset
variables; the linear predictions are treated as xb1, xb0, and zb rather than
xb1 + offset1, xb0 + offset0, and zb + offset_s.

{phang}
{opt scores} calculates equation-level score variables.

{pmore}
The first new variable will contain the derivative of the log likelihood with
respect to the selection equation.

{pmore}
The second new variable will contain the derivative of the log likelihood with
respect to the regression equation for regime 1.

{pmore}
The third new variable will contain the derivative of the log likelihood with
respect to the regression equation for regime 0.

{pmore}
The fourth new variable will contain the derivative of the log likelihood with
respect to the fourth equation ({hi:athrho1}).

{pmore}
The fifth new variable will contain the derivative of the log likelihood with
respect to the fifth equation ({hi:athrho0}).


{title:Examples}

{pstd}Estimation{p_end}
{phang2}{cmd:. switch_probit y x1 x2 x3, select(d x1 x2 x3 z1)}{p_end}
{phang2}{cmd:. switch_probit (y x1 x2 x3) (y x1 x2), select(d x1 x2 x3 z1)}

{pstd}Prediction{p_end}
{phang2}{cmd:. switch_probit y x1 x2 x3 x4, select(d= x1 x2 x3 z1 z2)}{p_end}
{phang2}{cmd:. predict at, at}{p_end}
{phang2}{cmd:. predict scr*, scores}

{pstd}Example from the {it:Stata Journal}{p_end}
{phang2}{cmd:. use switch_probit_example}{p_end}
{phang2}{cmd:. switch_probit works age age2 wedu_2-wedu_5 hhsize hhsize2 reg_*,}{break}
{cmd:select(migrant age age2 wedu_2-wedu_5 hhsize hhsize2 reg_* pmigrants)}{break}
{it:({stata `"do http://www.adeptanalytics.org/download/ado/switch_probit/switch_probit_example.do"':click to run})}


{title:Saved results}

{pstd}
{cmd:switch_probit} saves the following in {cmd:e()} (* indicates saved
parameters specific for {cmd:switch_probit}):

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:e(N)}}number of observations{p_end}
{synopt:{cmd:e(k)}}number of parameters{p_end}
{synopt:{cmd:e(k_eq)}}number of equations in {cmd:e(b)}{p_end}
{synopt:{cmd:e(k_eq_model)}}number of equations in overall model test{p_end}
{synopt:{cmd:e(k_aux)}}number of auxiliary parameters{p_end}
{synopt:{cmd:e(k_dv)}}number of dependent variables{p_end}
{synopt:{cmd:e(df_m)}}model degrees of freedom{p_end}
{synopt:{cmd:e(ll)}}log likelihood{p_end}
{synopt:{cmd:e(ll_0)}}log likelihood, constant-only model ({cmd:noskip} only){p_end}
{synopt:{cmd:e(N_clust)}}number of clusters{p_end}
{synopt:{cmd:e(chi2_c)}}value of the chi-squared test on the equality of both correlation coefficients to 0{p_end}
{synopt:{cmd:e(p_c)}*}probability of rejecting the chi-squared test{p_end}
{synopt:{cmd:e(p)}}significance of comparison test{p_end}
{synopt:{cmd:e(rho1)}*}estimated coefficient of correlation between the error terms of the selection equation and the outcome equations in regime 1{p_end}
{synopt:{cmd:e(rho0)}*}estimated coefficient of correlation between the error terms of the selection equation and the outcome equations in regime 0{p_end}
{synopt:{cmd:e(rank)}}rank of {cmd:e(V)}{p_end}
{synopt:{cmd:e(rank0)}}rank of {cmd:e(V)} for constant-only model{p_end}
{synopt:{cmd:e(ic)}}number of iterations{p_end}
{synopt:{cmd:e(rc)}}return code{p_end}
{synopt:{cmd:e(converged)}}{cmd:1} if converged, {cmd:0} otherwise{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:e(cmd)}}{cmd:switch_probit}{p_end}
{synopt:{cmd:e(depvar)}}name of dependent variable{p_end}
{synopt:{cmd:e(wtype)}}weight type{p_end}
{synopt:{cmd:e(wexp)}}weight expression{p_end}
{synopt:{cmd:e(title)}}title in estimation output{p_end}
{synopt:{cmd:e(clustvar)}}name of cluster variable{p_end}
{synopt:{cmd:e(offset1)}}offset for selection equation{p_end}
{synopt:{cmd:e(offset2)}}offset for equation 1{p_end}
{synopt:{cmd:e(offset3)}}offset for equation 0{p_end}
{synopt:{cmd:e(chi2type)}}{cmd:Wald} or {cmd:LR}; type of model chi-squared
	test{p_end}
{synopt:{cmd:e(chi2_ct)}}{cmd:Wald} or {cmd:LR}; type of model chi-squared test
	corresponding to {cmd:e(chi2_c)}{p_end}
{synopt:{cmd:e(vce)}}{it:vcetype} specified in {cmd:vce()}{p_end}
{synopt:{cmd:e(vcetype)}}title used to label Std. Err.{p_end}
{synopt:{cmd:e(opt)}}type of optimization{p_end}
{synopt:{cmd:e(ml_method)}}type of {cmd:ml} method{p_end}
{synopt:{cmd:e(user)}}name of likelihood-evaluator program{p_end}
{synopt:{cmd:e(technique)}}maximization technique{p_end}
{synopt:{cmd:e(crittype)}}optimization criterion{p_end}
{synopt:{cmd:e(properties)}}{cmd:b V}{p_end}
{synopt:{cmd:e(predict)}}program used to implement {cmd:predict}{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Matrices}{p_end}
{synopt:{cmd:e(b)}}coefficient vector{p_end}
{synopt:{cmd:e(ilog)}}iteration log (up to 20 iterations){p_end}
{synopt:{cmd:e(gradient)}}gradient vector{p_end}
{synopt:{cmd:e(V)}}variance-covariance matrix of the estimators{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Functions}{p_end}
{synopt:{cmd:e(sample)}}marks estimation sample{p_end}
{p2colreset}{...}


{title:Authors}

    Michael Lokshin
    The World Bank
    Washington, DC
    mlokshin@worldbank.org

    Zurab Sajaia
    The World Bank
    Washington, DC


{title:Also see}

{p 4 14 2}Article:  {it:Stata Journal}, volume 11, number 3: {browse "http://www.stata-journal.com/article.html?article=st0233":st0233}
{p_end}

{p 5 14 2}
Online:  {helpb biprobit:[R] biprobit}, {helpb heckprob:[R] heckprob},
{helpb ml:[R] ml}
{p_end}
