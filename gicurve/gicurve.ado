*! version 2.15 05Sep2008 M. Lokshin

#delimit ;

program define gicurve, rclass;
	version 8.0;

        syntax using [if] [in] [fweight aweight iweight],
                 var1(varname numeric)                   // expenditure in the first year
                 var2(string)                            // expenditure in the second year
                 [OUTputfile(string)                     // dta file with saved results
                 HCindex(numlist min=1 max=1 >=0 <=100)  // P0 to calculate the PPGrowth at
                 PLine(numlist min=1 max=1)              // pline to calculate the PPGrowth at
                 np(int 100)                             // number of percentiles to use
                 YPeriod(int 1)                          // number of years between two periods
                 ratio                                   // ratio instead of percentage growth
                 noGraph                                 // supresses graph output
                 ginmean                                 // outputs growth in means line
                 gatmedian                               // outputs growth in median line
                 mgrpp                                   // outputs mean growth at a percentile line
                 meangr                                  // output mean growth rate line
                 bands(int 100)                          // number of bands for gicurve smoothing
                 knots(int 1)                            // number of nots for gicurve smoothing
                 ci(numlist min=1 max=2 integer >0)      // generates confidence intervals
                 ci_u(numlist min=1 max=1)               // cuts the graph of ci on specified number
                 ci_l(numlist min=1 max=1)               // cuts the graph of ci on specified number
                 minmax                                  // delets from the graph max and min values
                 addplot(string) *]                      // add extra plot to the graph
		 ;

	local main_file $S_FN;  // name of the file that is currently opened

	tempfile ci_file _bs1;
	tempvar welf1 welf2 pctl pr_growth hcnt_pl;

   	quietly {;

        tempfile _tf1;
             summarize `var1'      `if' `in' [`weight' `exp'], meanonly;
             local _tot1=r(mean);  // absolute value in year 1
             local n_obs1 = r(N);  // number of observation in the first file
             pctile `welf1'=`var1' `if' `in' [`weight' `exp'], nq(`np') genp(`pctl');
             _pctile `var1'        `if' `in' [`weight' `exp'], p(50); local w1=r(r1);
             keep `welf1' `pctl'; drop if missing(`welf1');
             sort `pctl'; save `_tf1', replace;

        use `using', clear;
             summarize `var2'      `if' `in' [`weight' `exp'], meanonly;
             local _tot2=r(mean); // absolute value in year 2
             local n_obs2 = r(N); // number of observation in the second file
             pctile `welf2'=`var2' `if' `in' [`weight' `exp'], nq(`np') genp(`pctl');
             _pctile `var2'        `if' `in' [`weight' `exp'], p(50); local w2=r(r1);
             keep `welf2' `pctl'; drop if missing(`welf2');
             sort `pctl';

	merge `pctl' using `_tf1'; drop _merge;
	order `pctl' `welf1' `welf2';

        if  missing("`ratio'") {;
            generate `pr_growth' =(`welf2'-`welf1')/`welf1';    // growth of percentiles
            local    tot_growth  =(`_tot2'-`_tot1')/`_tot1';    // growth of the mean
            local    med_growth  =(`w2'-`w1')/`w1';             // growth at the median
            sum      `pr_growth', meanonly; local mean_of_g=r(sum)/(`np'-1);

            replace  `pr_growth' = ((`pr_growth' +1)^(1/`yperiod')-1)*100; // adjust for annual growth
            local     mean_of_g  = ((`mean_of_g' +1)^(1/`yperiod')-1)*100;
            local     med_growth = ((`med_growth'+1)^(1/`yperiod')-1)*100;
            local     tot_growth = ((`tot_growth'+1)^(1/`yperiod')-1)*100;
        }; // endif
        else {;      // case when request ratio instead of percentage growth
            generate `pr_growth' = `welf2'/`welf1';                     // growth of percentiles
            local    tot_growth  = `_tot2'/`_tot1';                     // growth of the mean
            local    med_growth  = `w2'/`w1';                           // growth at the median
	    summarize pr_growth'; local mean_of_g=r(sum)/(`np'-1);
        }; // endelse

        // integration if poverty line is specified
        if !missing("`pline'") {;
                summarize `pctl' if `welf1'<=`pline', meanonly;
                local hcount_pl = r(max); // determine headcount based on pl
                summarize `pr_growth' if `welf1'<=`pline', meanonly;
		local  integral = r(sum)/((`hcount_pl'/100)*`np');
		};

        // integration if headcount is specified
        if !missing("`hcindex'") {;
        	summarize `pr_growth' if `pctl'<=`hcindex', meanonly;
			local  hcount_pl = `hcindex';
			local  integral   = r(sum)/((`hcount_pl'/100)*`np');
		};

        // neither headcount no pline is specified *
        if (missing("`hcindex'") & missing("`pline'")) {;
			forvalues i=10(5)30 {;
            	summarize `pr_growth' if (`pctl'<=`i'), meanonly;
				local intgrl_`i' = r(sum)/((`i'/100)*`np');
			}; /* end forvalues */
   		};

	tempvar intgrl1 gr_in_mean gr_in_median mean_of_growth;

	generate `intgrl1'        = sum(`pr_growth')/_n;
	generate `gr_in_mean'     =`tot_growth';
	generate `gr_in_median'	  =`med_growth';
	generate `mean_of_growth' =`mean_of_g';

    label variable `pr_growth'      "Growth incidence curve";
    label variable `intgrl1'        "Mean growth rate for poorest p%";
	label variable `gr_in_mean' 	"Growth rate in mean";
	label variable `gr_in_median' 	"Growth rate at median";
    label variable `pctl'           "Percentiles";
	label variable `mean_of_growth' "Mean of growth rates";

        tempfile ff1;
        save `ff1';

        if ~missing("`ci'") {; // Bootstraping for standard errors
           local nreps  : word 1 of `ci';
           local c_level: word 2 of `ci';

           use "`main_file'", clear;

           set seed 123456789;

           forvalues i = 1/`nreps' {;
                preserve;

                bsample `if';
                         pctile `welf1'=`var1' `if' [`weight' `exp'], nq(`np') genp(`pctl');
                         keep `welf1' `pctl'; drop if missing(`welf1');
                         sort `pctl'; save `_bs1', replace;

                use `using', clear;
                bsample `if';
                        pctile `welf2'=`var2' `if' [`weight' `exp'], nq(`np') genp(`pctl');
                        keep `welf2' `pctl'; drop if missing(`welf2'); sort `pctl';

                merge `pctl' using `_bs1';   drop _merge;
                order `pctl' `welf1' `welf2';

                if missing("`ratio'")
                generate `pr_growth' = ((((`welf2'-`welf1')/`welf1') +1)^(1/`yperiod')-1)*100;  // adjust for anual growth
                else generate `pr_growth' = `welf2'/`welf1';                                    // growth of percentiles        */
                rename `pr_growth' pg`i'; keep `pctl' pg`i'; sort `pctl';

                if (`i'==1) save `ci_file', replace;
                else {; merge  `pctl'  using `ci_file'; drop _merge; sort `pctl'; save `ci_file', replace;};
                restore;
                noisily display "." _continue;
           }; // end forvalues

           use `ci_file', clear;
           aorder   pg*;
           egen     pg_sd   =   rsd(pg1-pg`nreps');
           egen     pg_mean = rmean(pg1-pg`nreps');

           if missing("`c_level'") local c_level=95;

           local c_l       =(100-(100-`c_level')/2)/100;
           local d_freedom =((`n_obs1'+`n_obs2')/2)-1;

           generate pg_ci_u = pg_mean-invttail(`d_freedom',`c_l')*pg_sd;
           generate pg_ci_l = pg_mean+invttail(`d_freedom',`c_l')*pg_sd;
           keep     `pctl' pg_ci_l pg_ci_u;
           sort     `pctl';
           save     `ci_file', replace;

           use `ff1', clear;

           sort `pctl'; merge `pctl' using `ci_file'; drop _merge;
           label variable pg_ci_u "Upper `c_level'% confidence bound";
           label variable pg_ci_l "Lower `c_level'% confidence bound";

           if (!missing("`ci_u'")) replace pg_ci_u=`ci_u' if (pg_ci_u >`ci_u');
           if (!missing("`ci_l'")) replace pg_ci_l=`ci_l' if (pg_ci_l <`ci_l');

	   local keepci="pg_ci_u pg_ci_l";

        }; // if1

        if !missing("`outputfile'") {;
           preserve;
           rename `pr_growth'      pr_growth;
           rename `intgrl1'        intgrl1;
           rename `gr_in_mean'     gr_in_mean;
           rename `gr_in_median'   gr_in_median;
           rename `pctl'           pctl;
           rename `mean_of_growth' mean_of_growth;
	   sort pctl;
           keep  pctl pr_growth gr_in_mean gr_in_median mean_of_growth intgrl1 `keepci';
           order pctl pr_growth gr_in_mean gr_in_median mean_of_growth intgrl1 `keepci';
           save `outputfile', replace;
           restore;
        };

        };   // end quietly

        if missing("`graph'") {;  // nograth is not specified         if1
                if !missing("`mgrpp'")     local gstring                  `intgrl1';
                if !missing("`ginmean'")   local gstring `gstring'     `gr_in_mean';
                if !missing("`gatmedian'") local gstring `gstring'   `gr_in_median';
                if !missing("`meangr'")    local gstring `gstring' `mean_of_growth';

                if !missing("`minmax'") {; // if2
                   summarize `pr_growth', meanonly;
                   quietly drop if (`pr_growth' >= r(max)) | (`pr_growth'<=r(min));
                }; // if2

                local g_mspline mspline `pr_growth' `pctl', bands(`bands') n(`knots');
                local g_rarea   rarea    pg_ci_u pg_ci_l `pctl', sort bcolor(gs14);
                local g_line    line    `gstring'       `pctl';

                if !missing("`ci'")      local g_main `g_rarea' || `g_mspline';
                else                     local g_main `g_mspline';
                if !missing("`gstring'") local g_main `g_main' || `g_line',; // important to have a comma here
                if !missing("`addplot'") local g_main `g_main' || `addplot';

                twoway `g_main' `options';

        };// end if1

	display _newline(2);
    display in green "   Growth rate in mean"          _col(32) "= " %5.2f `tot_growth';
    display in green "   Growth rate at median"        _col(32) "= " %5.2f `med_growth';
	display in green "   Mean percentile growth rate"  _col(32) "= " %5.2f `mean_of_g' _newline;

        if !missing("`pline'") disp in green "   Poverty line"   _col(32) "= " `pline' _newline;

	display in green "   Corresponding percentile {c |} Rate of pro-poor growth ";
	display in green "{hline 28}{c +}{hline 25}";
        if (!missing("`pline'") | !missing("`hcindex'")) {;
			display _col(12) as res %5.2f `hcount_pl' _col(29) as text "{c |}" _col(37) as res %5.2f `integral';
		};
	else {;
		display _col(12) as res 10  _col(29) as text "{c |}" _col(37) as res %5.2f `intgrl_10' ;
		display _col(12) as res 15  _col(29) as text "{c |}" _col(37) as res %5.2f `intgrl_15' ;
		display _col(12) as res 20  _col(29) as text "{c |}" _col(37) as res %5.2f `intgrl_20' ;
		display _col(12) as res 25  _col(29) as text "{c |}" _col(37) as res %5.2f `intgrl_25' ;
		display _col(12) as res 30  _col(29) as text "{c |}" _col(37) as res %5.2f `intgrl_30' ;
		};
	display in green "{hline 28}{c BT}{hline 25}";

	use "`main_file'", clear;

        // returned values
        return scalar gr_inmean   = `tot_growth';
        return scalar gr_inmedian = `med_growth';
        return scalar mean_pct_gr = `mean_of_g' ;

        if (missing("`hcindex'") & missing("`pline'")) {;
           return scalar gr_in_10 = `intgrl_10';
           return scalar gr_in_15 = `intgrl_15';
           return scalar gr_in_20 = `intgrl_20';
           return scalar gr_in_25 = `intgrl_25';
           return scalar gr_in_30 = `intgrl_30';
           };
        else return scalar gr_in_sp  = `integral'; // returns specified growth rate

        return scalar mean_var1   = `_tot1';
        return scalar mean_var2   = `_tot2';
        return scalar med_var1    = `w1';
        return scalar med_var2    = `w2';
        return scalar N1          = `n_obs1';
        return scalar N2          = `n_obs2';

        return local  title  "Growth Incidence Curve";
        return local  cmd    "gicurve";
        return local  var1   `var1';
        return local  var2   `var2';
end;
