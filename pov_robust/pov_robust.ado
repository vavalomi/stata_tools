*! version 3.05  27Feb2008 M. Lokshin
#delimit ;

program define pov_robust, rclass;
version 8.0;

syntax [using] [if] [aw fw], var1(varname numeric) [var2(string) automax pi pd ps *];

if missing(`"`using'"')!=missing("`var2'") {;
	display as error "using and var2 options must be used together";
	exit 198;
};

preserve;

tempfile file1;
tempvar c_year c_welf;

quietly {;
	if ~missing("`if'") keep `if'; // restrict observation for user-specified if conditions
	generate byte `c_year'=1;
	if ~missing("`exp'") {;
		tempvar c_wgt;
		generate `c_wgt' `exp';
	};
	rename `var1' `c_welf';
	keep `c_welf' `c_year' `c_wgt';
	save `file1';

	if ~missing(`"`using'"') {;
		use `using', clear;
		if ~missing("`if'") keep `if'; // restrict observation for user-specified if conditions
		generate byte `c_year'=2;
		rename `var2' `c_welf';
		if ~missing("`exp'") {;
			generate `c_wgt' `exp';
		};
		keep `c_welf' `c_year' `c_wgt';
		append using `file1';
		local years = 2;
		local lab " in year \`y'";
	};
	else local years = 1;

	if ~missing("`exp'") local exp "=`c_wgt'";

	drop if `c_welf'==.;
	label variable `c_welf' "Welfare indicator";

	forvalues y = 1/`years' {;
		tempvar step`y' pi`y';
		cumul `c_welf' if `c_year'==`y' [`weight' `exp'], gen(`pi`y''); // poverty incidence curve
		label variable `pi`y'' "Poverty incidence curve`lab'";

		sort `c_year' `c_welf';
		generate `step`y''=(`c_welf'-`c_welf'[_n-1]) if `c_year'==`y';
	 	 replace `step`y''=. if `step`y''<0;	// for the first obs in year 2

		if ~missing("`pd'`ps'") {;
			tempvar pd`y';
			generate `pd`y''=sum(`step`y''*`pi`y'') if `c_year'==`y'; // poverty deficit curve
			label variable `pd`y'' "Poverty depth curve`lab'";
		};
		if ~missing("`ps'") {;
			tempvar ps`y';
			generate `ps`y''=sum(`step`y''*`pd`y'') if `c_year'==`y'; // poverty severity curve
			label variable `ps1' "Poverty severity curve`lab'";
		};
	};
	                      local p_line `pi1' `pi2';
	if (!missing("`ps'")) local p_line `ps1' `ps2';
	if (!missing("`pd'")) local p_line `pd1' `pd2';
}; // end quietly

 // Everything is prepared. Now - rescale
if ~missing("`automax'") {;
	local FLAG=1;
	while `FLAG' {;
    	// Determine if we can divide ALL the variables in p_line by 1000
		foreach pline of local p_line {;
			summarize `pline',meanonly;
			if `r(max)'<1000 local FLAG=0; // if the max is less than 1000 then switch off the flag
		};
		// If FLAG is set, then the max of each variable in p_line is more than 1000
		if `FLAG' foreach pline of local p_line {;
					quietly replace `pline'=`pline'/1000;
				  };
	};

	summarize `c_welf',meanonly;
	roundnice `r(max)';	     // make a nice maximum
	round_det `r(roundnice)';  // round to 000s, mlns or blns
	local max_c_welf=`r(number)'; local max_c_welf_u="`r(units)'";
	if "`max_c_welf_u'"=="'000" replace `c_welf'=`c_welf'/1000;
	if "`max_c_welf_u'"=="mln"  replace `c_welf'=`c_welf'/1000000;
	if "`max_c_welf_u'"=="bln"  replace `c_welf'=`c_welf'/1000000000;
	local stepnice=`max_c_welf'/5;

	local lbl : variable label `c_welf';

	twoway line `p_line' `c_welf',
		xtitle("`lbl', `max_c_welf_u'", size(3) margin(0 0 0 2))
		xlabel(0(`stepnice')`max_c_welf', labsize(2.9))
		`options';
};
else twoway line `p_line' `c_welf', `options';

restore;

end;

program define round_det, rclass;
	args number;
   	local units "units";
	if `number'>1000 {; local number=round(`number'/1000, 0.1); local units="'000";	};
	if `number'>1000 {; local number=round(`number'/1000, 0.1); local units="mln";	};
	if `number'>1000 {; local number=round(`number'/1000, 0.1); local units="bln";	};
	return local number `number'; return local units "`units'";
end;

#delimit cr

program roundnice, rclass
   local roundwhat `1'
   local round_method="ceil"
   if "`2'"=="0" local round_method="round"
   if "`2'"=="-1" local round_method="floor"

   local l1 50 100 500 1000 2000 5000 10000 20000 100000 200000 500000 1000000 5000000   /* If X is less than this number */
   local l2  1  10  50  100  200  500  1000  2000   5000  10000  20000   50000  100000   /* round it to multiples of this number */

   capture confirm number `roundwhat'
   if _rc {
     window stopbox stop "Error, `roundwhat' is not a number or missing"
     exit
   }

   local roundnice = `roundwhat'

   local indx=1
   foreach t of local l1 {
     if `roundwhat'<`t' {
       local r=word("`l2'",`indx')
       local roundnice = `=`round_method'(`roundwhat'/`r')*`r''
       continue, break
     }
     local indx=`indx'+1
   }
   return scalar roundnice=`roundnice'
end



































