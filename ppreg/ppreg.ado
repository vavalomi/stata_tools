// Last saved 8Dec2008
#delimit ;

program define ppreg;
	version 9;
	if replay() {;
		if ("`e(cmd)'" ~= "ppreg") error 301;
		ppreg_replay `0';
	};
	else ppreg_est `0';
end; // end program ppreg

program define ppreg_replay;
	syntax [,Level(cilevel)];
	_coef_table_header;
		display;
		ereturn display, level(`level');
end; // end program ppreg_reply

program define ppreg_est, eclass;
	syntax [varlist(default=none)] [if] [in] [aweight pweight iweight fweight], Time(varname) Cohort(varname);
  	quietly {;
		tabulate `time';   scalar T=r(r);
		tabulate `cohort'; scalar C=r(r);

   		marksample touse;
		markout `touse' `varlist' `time' `cohort', strok;
		gettoken y x : varlist;
		preserve;
		keep if `touse';

		tempvar w e tw;
		if (~missing("`weight'")) {;
			generate double `w'`exp';
			local ww "[`weight'=`w']";
		};
		else generate `w'=1;

			egen double `tw' = total(`w'), by(`time' `cohort');
			summarize `w'; local ttw = r(sum);
			foreach var of varlist `varlist' {;
				tempvar dvar tm`var';
				egen double `tm`var''   = total(`var'*`w'/`ttw');
				egen double `dvar' = total(`var'*`w'/`tw'), by(`time' `cohort');
			  	replace     `dvar' = (`var' - `dvar') / sqrt((`tw' - 1)*`tw');
				local dvarlist  "`dvarlist' `dvar'";
				local tmvarlist "`tmvarlist' `tm`var''";
			};
			tempname S Syy Sxy Sxx M Mxy Mxx invOMEGA b V V1 Xb EX;
			matrix accum `S'   = `dvarlist' `ww', noconstant;
  			matrix 		 `S'   = `S' / (C*T);
			matrix       `Sxy' = `S'[2..., 1];
			matrix       `Sxx' = `S'[2..., 2...];
			collapse (mean) `varlist' (min) `tmvarlist' (count) `w'=`y' `ww', by(`time' `cohort');

			if (T*C > r(r)) {; noisily display as error "not every cohort-time cell has observations"; exit 2000; };

			egen double `tw' = total(`w'), by(`cohort');
			foreach var of varlist `varlist' {;
				tempvar mean;
		  		egen double `mean' = total(`var'*`w'/`tw'), by(`cohort');
				replace      `var' = (`var' - `mean');
			};

			matrix accum `M'   = `varlist' `ww', noconstant;
			matrix       `M'   = `M' / (C*(T-1));
			matrix       `Mxy' = `M'[2..., 1];
			matrix       `Mxx' = `M'[2..., 2...];
	  		matrix   `invOMEGA' = inv(`Mxx'-`Sxx');
			matrix          `b' = (`invOMEGA' * (`Mxy'-`Sxy'))';
   			matrix rownames `b' = `y';

			matrix colnames `b' = `x';
			matrix score double `Xb' = `b';
			generate     double  `e' = `y' - `Xb';

			matrix accum `EX' = `e' `x' `ww', noconstant;
			matrix       `EX' = `EX' / (C*(T-1));
			matrix 	     `V'  = `invOMEGA'*(`Mxx'*`EX'[1,1]                              +`EX'[2...,1]*`EX'[1,2...])*`invOMEGA'/(C*(T-1));
			matrix       `V1' = `invOMEGA'*(`Mxx'*(`M'[1,1]-2*`Mxy''*`b''+`b'*`Mxx'*`b'')+`EX'[2...,1]*`EX'[1,2...])*`invOMEGA'/(C*(T-1));
	   	restore;
		count if `touse';

			matrix `V' =`V'+`V1'/`r(N)'*C*T;

		 	matrix colnames `V' = `x';
  		 	matrix rownames `V' = `x';

		ereturn post `b' `V', esample(`touse') depname("`y'") obs(`r(N)');

		ereturn local wtype  "`weight'";
		ereturn local wexp   "`exp'";
		ereturn local cmd    "ppreg";
		ereturn local depvar "`y'";
		ereturn local title "Pseudo panel regression";

   	}; // quietly
	ppreg_replay;
end; // program ppreg_est

