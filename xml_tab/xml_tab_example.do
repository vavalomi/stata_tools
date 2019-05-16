#delim ;

capture program drop table1;
program define table1;
 	quietly {; // quietly
	use xml_tab_example, clear;

	local dstr ulog_wage mother_smoked father_smoked
			   age married divorced widowed relig_* educat_* speaks_*
			   lhhsize s_* urban;

	tabstat `dstr' if ufilter, by(smokes) s(mean sd) save notot;
	matrix SMOKERS  = r(Stat2);
	matrix NSMOKERS = r(Stat1);
	matrix TAB = SMOKERS',NSMOKERS';

	xml_tab TAB, format(S2103 N2203 N2123) sheet(Table1) noi
				 rblank(age         "Marital status" S2123,
					 	widowed     "Religion"       S2123,
						relig_5     "Education"      S2123,
						hours_pweek "Subjective health assesment" S2123
				 )
				 replace
				 ;

	}; // end quietly
end; // program table1

/*******************************************************************/
/*******************************************************************/
capture program drop table2;
program define table2;
	quietly {;
	use http://siteresources.worldbank.org/INTPOVRES/Resources/smoking, clear;
	keep if (ufilter);
	local dstr $ind $edu $hh $dst;
	reg ulog_wage smokes `dstr' if (ufilter);                           		estimates store ols;
	reg smoke mother_sm father_sm `dstr' if (ufilter);                  		estimates store _2sls_1;
	ivreg ulog_wage `dstr' (smokes= mother_sm father_sm) if (ufilter), first;  estimates store _2sls_2;

	noi xml_tab ols _2sls_1 _2sls_2,
		rblank(age3         "Marital Status"     S2123,
			   widowed      "Religion"           S2123,
			   relig_6      "Education"          S2123,
			   lhours_pweek "HH Characteristics" S2123)
		drop(dist_* mother_e* father_e*)
		style(s1)
		sheet(Table 2)
		title(2SLS estimation) append
		;
	}; // end quietly
end; // program table2

preserve;
clear;
set more off;
set matsize 400;
set mem 10m;

global ind    age age2 age3 married divorced widowed relig_2-relig_5;
global edu    educat_2-educat_6 speaks_* lhours_pweek mother_edu_2-mother_edu_6 father_edu_2-father_edu_6;
global hh     lhhsize* urban s_*;
global dst    dist_2-dist_15 dist_17-dist_36;

// Descriptive statistics
table1;
// OLS and 2SLS
table2;

restore;
display _n;
display as text  "type {help xml_tab:help xml_tab} for more examples"
