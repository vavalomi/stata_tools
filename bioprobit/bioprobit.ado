*! version 1.11  17Aug2007
#delimit ;

program define bioprobit;
	version 9;
	if replay() {;
		if ("`e(cmd)'" ~= "bioprobit") error 301;
		bioprobit_replay `0';
	}; // end if
	else bioprobit_mx `0';
end; // end program bioprobit

program define bioprobit_replay;
	syntax [,Level(cilevel)];
	local rho diparm(athrho  , tanh label("rho"));
	_coef_table_header;
	display;
	_coef_table, level(`level') `rho' notest;
   	display	in green e(chi2_ct) " test of indep. eqns. :"
		_col(38) "chi2(" in yellow "1" in green ") = " in yellow %8.2f e(chi2_c)
		_col(59) in green  "Prob > chi2 = " in yellow %6.4f e(p_c);
		if strpos("`e(title)'", "Simultaneous")	display	in smcl in green "{hline 78}";

end; // end program bioprobit_reply

program define bioprobit_mx, eclass;

	gettoken first : 0, match(paren);
	if missing("`paren'") {;	                 // syntax 1, bivariate ordered probit models
 	  	syntax varlist [if] [in] [pweight iweight fweight]
		   			  		   	[, offset1(varname) offset2(varname) COLlinear
								   Robust CLuster(varname) Level(cilevel) end *];
		gettoken y1 varlist : varlist;
		gettoken y2 varlist : varlist;
		local x1 `varlist';
		local x2 `varlist';
		local pref "B";
	};
	else {;                                      // syntax 2, seemingly unrelated bivariate ordered probit model
 	  	syntax anything(id="equation id" equalok) [if] [in] [pweight iweight fweight]
		   			     	   	[, offset1(varname) offset2(varname) COLlinear
								   Robust CLuster(varname) Level(cilevel) end *];
		gettoken eq1 eq2_ :  anything,  match(parns) bind;
		tokenize "`eq1'", parse("=");
		if ("`2'"!="=") {;
			tokenize "`eq1'";
			local y1 `1';
			macro shift;
			local x1 `*';
		};
		else {;
			local y1 `1';
			local x1 `3';
		};
		gettoken eq2 : eq2_ , match(parns) bind;
		if ("`parns'"!="(") local eq2 "`eq2_'";
		tokenize "`eq2'", parse("=");
		if ("`2'"!="=") {;
			tokenize "`eq2'";
			local y2 `1';
			macro shift;
			local x2 `*';
		};
		else {;
			local y2 `1';
			local x2 `3';
		};
		local pref "Seemingly unrelated b";
	}; // if
	global END =~missing("`end'");

	quietly {;

	// define standard dep vars
	tempvar sy1 sy2;
	egen `sy1' = group(`y1');
	egen `sy2' = group(`y2');

	marksample touse;
	markout `touse' `sy1' `sy2' `x1' `x2' `cluster' `offset1' `offset2', strok;

	mlopts mlopts, `options';

	if (~missing("`offset1'")) local offo1 "offset(`offset1')";
    if (~missing("`offset2'")) local offo2 "offset(`offset2')";
   	if (~missing("`weight'"))  local weight "[`weight'`exp']";
	if (~missing("`cluster'")) local clopt cluster(`cluster');
	if "`weight'" == "pweight" | (~missing("`cluster'")) local robust "robust";

    // Remove collinear variables
	noisily _rmdcoll `sy1' `x1' `weight' if `touse', `collinear';
	local x1 "`r(varlist)'";
noi tabulate `sy2' if `touse';

	noisily _rmdcoll `sy2' `x2' `weight' if `touse', `collinear';
	local x2 "`r(varlist)'";

	summarize `sy1'; global NC1  = r(max);
 	summarize `sy2'; global NC2  = r(max);    // number of categiries
	if ($NC1==1) {;
		noisily display as error "`y1' does not vary";
		exit 2000;
	};
	if ($NC2==1) {;
		noisily display as error "`y2' does not vary";
		exit 2000;
	};
    local NC1_1   = $NC1-1;
	local NC2_1   = $NC2-1;

	// define initial values
	tempname Ib_op1 Ib_op2 Ic_op Ic_op1 Ic_op2 I_rho TMP ll_0;
   	scalar `ll_0' = 0;

   	oprobit `sy1' `x1' `weight' if `touse', `offo1';
		matrix 			`TMP' = e(b);
		matrix       `Ib_op1' = `TMP'[1,"`sy1':"];
		matrix coleq `Ib_op1' = `y1';
		matrix 		 `Ic_op1' = `TMP'[1,"cut1:_cons".."cut`NC1_1':_cons"];
	scalar `ll_0'=`ll_0'+e(ll);
	if $END {;
		tempvar xb1;
		matrix score `xb1' = `Ib_op1';
	};
	oprobit `sy2' `xb1' `x2' `weight' if `touse', `offo2';
		matrix 			`TMP' = e(b);
		matrix 		 `Ib_op2' = `TMP'[1,"`sy2':"];
		matrix coleq `Ib_op2' = `y2';
		matrix 		 `Ic_op2' = `TMP'[1,"cut1:_cons".."cut`NC2_1':_cons"];
	scalar `ll_0'=`ll_0'+e(ll);


	matrix `Ic_op'  = `Ic_op1'[1,1];
	forvalues i = 2/`NC1_1' {; matrix `Ic_op' =`Ic_op',sqrt(`Ic_op1'[1,`i']-`Ic_op1'[1,`i'-1]); };
	matrix `Ic_op'  = `Ic_op', `Ic_op2'[1,1];
	forvalues i = 2/`NC2_1' {; matrix `Ic_op' =`Ic_op',sqrt(`Ic_op2'[1,`i']-`Ic_op2'[1,`i'-1]);	};

	correlate `sy1' `sy2' if `touse';
	matrix          `I_rho' = .5;//atanh(r(rho));
	matrix colnames `I_rho' = athrho:_cons;

	local cuts;
	forvalues k = 1/`NC1_1' {; local cuts "`cuts' /cut1`k'"; local ceqs "`ceqs' cut1`k'"; };
	forvalues k = 1/`NC2_1' {; local cuts "`cuts' /cut2`k'"; local ceqs "`ceqs' cut2`k'"; };
	matrix coleq    `Ic_op' = `ceqs';
	matrix colnames `Ic_op' = _cons;
	}; // end quietly

 	// maximization
	if $END {;
		matrix `I_rho'  =`I_rho',`Ib_op2'[1,1];
		matrix colnames `I_rho' = athrho:_cons gamma:_cons;

		matrix `Ib_op2' =`Ib_op2'[1,2...];
		local gamma "/gamma";
		local title "title(Simultaneous bivariate ordered probit regression)";
	};
	else local title "title(`pref'ivariate ordered probit regression)";

	ml model d2 bioprobit_d2
	     (`y1' : `sy1' = `x1', noconstant)
	     (`y2' : `sy2' = `x2', noconstant)
		 /athrho
		 `gamma'
		 `cuts'
		 `weight'
		 if `touse'
      	 ,
		 `title'
		 difficult
		 collinear
		 missing
	     search(on)
		 init(`Ib_op1' `Ib_op2' `I_rho' `Ic_op')
	     maximize
		 `clopt'
		 `robust'
  	     `mlopts'
		 ;

	tempname b D D1 V;
	matrix `b'  =e(b);
	matrix `V'  =e(V);
	local f = colnumb(`b',"cut11:_cons") -1;
	matrix `D' = I(`f'), J(`f', `NC1_1'+`NC2_1',0);                                   // matrix of derivatives
	matrix `D1' = (J(1,`f',0), 1, J(1,`NC1_1'+`NC2_1'-1, 0));
	matrix `D'  = `D' \ `D1';
	forvalues k = 2/`NC1_1' {;
		matrix `D1'[1,`f'+`k'] = 2*`b'[1,`f'+`k'];
		matrix `D'  = `D' \ `D1';
	   	matrix `b'[1,`f'+`k'] =`b'[1,`f'+`k'-1]+`b'[1,`f'+`k']^2;
	};
	local f = colnumb(`b',"cut21:_cons") -1;
	if (`NC2_1' > 1) matrix `D1' = (J(1,`f', 0), 1, J(1,`NC2_1'-1, 0));
	else             matrix `D1' = (J(1,`f', 0), 1);
   	matrix `D'  = `D' \ `D1';
   	forvalues k = 2/`NC2_1' {;
		matrix `D1'[1,`f'+`k']  = 2*`b'[1,`f'+`k'];
		matrix `D'  = `D' \ `D1';
		matrix `b'[1,`f'+`k'] =`b'[1,`f'+`k'-1]+`b'[1,`f'+`k']^2;
	};
	matrix `V' = `D'*`V'*`D'';
    ereturn repost b=`b' V=`V';

	ereturn scalar ll_0 = `ll_0';               // loglikelihood for non-correlated case
	ereturn scalar k_aux= `NC1_1'+`NC2_1';      // identify ancillary parameters

	ereturn local cmd    "bioprobit";
	ereturn local predict "bioprobit_p";
	ereturn local depvar "`y1' `y2'";
	ereturn local offset1 `offset1';
	ereturn local offset2 `offset2';

	if missing("`robust'") {;
		ereturn local chi2_ct "LR";
		ereturn scalar chi2_c = 2 * (e(ll) - `ll_0');
	};
	else {;
		ereturn local chi2_ct "Wald";
		quietly test [athrho]_b[_cons];
		ereturn scalar chi2_c = r(chi2);
		if $END {;
			ereturn local
			quietly test [athrho]_b[_cons]+[gamma]_b[_cons]==0;
			ereturn scalar chi2_c2 = r(chi2);
			ereturn scalar p_c2    = chiprob(2, e(chi2_c2));
		};
	};
	ereturn scalar p_c = chiprob(1, e(chi2_c));
	global END;
	bioprobit_replay, level(`level');
end; // end program bioprobit_mx

