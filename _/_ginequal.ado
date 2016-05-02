// Last modified: 2-May-2016
# delimit ;
program define _ginequal;
	syntax newvarname =/exp [if] [in] [, BY(varlist) Weights(varname) INDex(string)];

	local inc : subinstr local exp "(" "";
	local inc : subinstr local inc ")" "";
	local inc = trim("`inc'");
	confirm variable `inc';

	tempvar badinc touse nk nk0 fik meanyk meanyk0 lnm pyk i tmptmp N;

	quietly {; // quietly
		marksample touse , novarlist;
		count if `inc' <= 0 & `touse';
		local ct = r(N);
		if `ct' > 0 {;
			noisily display;
			noisily display in blue "Warning: `inc' has `ct' values <= 0." _c;
			noisily display in blue " Not used in calculations";
			};
		generate `badinc' = 0;
		replace  `badinc' =. if (`inc' <= 0 | `inc'==.);
		markout `touse' `badinc';

		markout `touse' `by', strok;

		if ("`by'" == "") {;
			tempvar by;
			generate byte `by' = 1;
		};
		else {;
			unab by : `by';
			local bystr "by(`by')";
		};

		generate `i' = .;
		generate `tmptmp' = .;

		gsort -`touse' + `by' `inc';

		if ("`weights'" == "")  {;
			tempvar wi;
			generate byte `wi' = 1;
			bys `by': replace `i' = _n;
		};
		else {;
			tempname wi;
			local wi `weights';
			bys `by': replace `tmptmp' = sum(`wi');
			replace `i' = ((2*`tmptmp')-`wi'+1)/2;
		};
		if ("`index'" == "") {;
			local index "gini";
		};

		count if (`wi' < 0);
		if (_result(1) > 0) {;
			noisily display as error "`wi' contains negative values";
			exit 402;
		};

		egen     `nk'  = sum(`wi') if `touse', by(`by');
		generate `fik' = `wi'/`nk' if `touse';
		egen `meanyk'  = sum(`fik'*`inc') if `touse', by(`by');

		// relative mean deviation
	    if ("`index'" == "rmd") {;
	       egen `varlist' = sum(`fik'*abs(`inc'-`meanyk')/(2*`meanyk')) if `touse', by(`by');
	    };
		// coefficient of variation
		else if ("`index'" == "cov") {;
        	egen    `varlist' = sum(`wi'*(`inc'-`meanyk')^2) if `touse', by(`by');
			egen    `N'       = count(`inc') if `touse', by(`by');
			replace `varlist' =  sqrt(`varlist'/`nk'*`N'/(`N'-1))/`meanyk' if `touse';
		};
        // standard deviation of logs
		else if ("`index'" == "sdl") {;
			egen `lnm' = sum(`fik'*ln(`inc')) if `touse', by(`by');
        	egen    `varlist' = sum(`wi'*(ln(`inc')-`lnm')^2) if `touse', by(`by');
        	egen    `N'       = count(`inc') if `touse', by(`by');
        	replace `varlist' =  sqrt(`varlist'/`nk'*`N'/(`N'-1)) if `touse';
		};
		// gini
	    else if ("`index'" == "gini") {;
	   	  sort `by' `inc';
			by `by': generate `pyk' = (2*sum(`wi') - `wi' + 1)/(2 * `nk' ) if `touse';
			egen `varlist' = sum(`fik'*(2/`meanyk')*`pyk'*(`inc'-`meanyk')) if `touse', by(`by');
	    };
		// mehran
		else if ("`index'" == "mehran") {;
        	egen `varlist' = sum(3*`wi'*`i'*(2*`nk'+1 -`i')*(`inc' - `meanyk')/(`nk'^3*`meanyk'))  if `touse', by(`by');
		};
		// piesch
		else if ("`index'" == "piesch") {;
        	egen    `varlist' = sum(`wi'*`i'*(`i'-1)*(`inc'-`meanyk')) if `touse', by(`by');
        	replace `varlist' = 3*`varlist'/(2*`nk'^3*`meanyk') if `touse';
		};
		// kakwani
		else if ("`index'" == "kakwani") {;
			egen    `varlist' = sum(`wi'*((`inc'^2+`meanyk'^2)^0.5)) if `touse', by(`by');
        	replace `varlist' = (1/(2-2^0.5))*((`varlist'/(`nk'*`meanyk')-2^0.5))  if `touse';
		};
		else {;
			egen     `nk0' = sum(`wi') if `touse' & `inc' > 0, by(`by');
			egen `meanyk0' = sum(`fik'*`inc') if `touse' & `inc' > 0, by(`by');
			// theil
			if ("`index'" == "theil") {;
        		egen `varlist' = sum(`wi'*((`inc'/`meanyk0')*(log(`inc'/`meanyk0')))/`nk0') if `touse', by(`by');
			};
			// mean log deviation
			else if ("`index'" == "mld") {;
        		egen `varlist' = sum(`wi'*(log(`meanyk0'/`inc'))/`nk0') if `touse', by(`by');
			};
			// GE -1
			else if ("`index'" == "entropy") {;
				egen    `varlist' = sum(`wi'*(`meanyk0'/`inc')) if `touse', by(`by');
        		replace `varlist' = ((`varlist'/`nk0')-1)/2 if `touse';
			};
    		// GE 2
			else if ("`index'" == "half") {;
        		egen    `varlist' = sum(`wi'*((`inc'/`meanyk0')^2)) if `touse', by(`by');
        		replace `varlist' = ((`varlist'/`nk0')-1)/2 if `touse';
			};
		};
		label variable `varlist' "(`index') `inc' `bystr'";
	}; // end quietly
end; // program _ginequal
