*! version 3.0.2  21Apr2008 M. Lokshin, Z. Sajaia
#delim ;

program define movestay, sortpreserve;
	version 8;
	if replay() {;
		if ("`e(cmd)'" ~= "movestay") error 301;
		movestay_replay `0';
		};
	else movestay_mx `0';
end; // end program movestay

program define movestay_replay;
	syntax, [Level(int $S_level)];

	local sigma0  diparm(lns0, exp  label("sigma0"));
	local rho0    diparm(r0  , tanh label("rho0"));
	local sigma1  diparm(lns1, exp  label("sigma1"));
	local rho1    diparm(r1  , tanh label("rho1"));

	ml display, level(`level') `sigma0' `sigma1' `rho0' `rho1';

	if ("`e(vcetype)'" != "Robust") local testtyp LR;
	else local testtyp Wald;
	display	in green "`testtyp' test of indep. eqns. :"
		_col(38) "chi2(" in yellow "2" in green ") = " in yellow %8.2f e(chi2_c)
		_col(59) in green  "Prob > chi2 = " in yellow %6.4f e(p_c);
	display	in smcl in green "{hline 78}";
end; // end program movestay_reply

program define movestay_mx, eclass;

	syntax 	anything(id="equations id" equalok) [pweight iweight fweight]
	        [if] [in], SELect(string) [CLuster(varname) Robust *];

	mlopts mlopts, `options';
	local coll `s(collinear)';

	quietly {;

	// read vector of parameters

	gettoken reg0 reg1_ :  anything,  match(parns) bind;
	if ("`parns'"!="(") {;
		local reg0 "`anything'";
		local reg1_ "";
	};
	tokenize "`reg0'", parse("()");
	if ("`1'"!="`reg0'") {;
		display as error "If more than one, equations must be enclosed in brackets";
		exit 198;
	};
	tokenize "`reg0'", parse("=");
	if ("`2'"!="=") {;
		tokenize "`reg0'";
		local y_reg0 `1';
		macro shift;
		local x_reg0 `*';
	};
	else {; // else_63
		if ("`4'"=="=") {;
			display as error "If more than one, equations must be enclosed in brackets";
			exit 198;
		};
	       	local y_reg0 `1'; local x_reg0 `3';
        }; // else_63
	unab x_reg0 : `x_reg0';

	if ("`reg1_'"!="") {; // if_64
	     	gettoken reg1  : reg1_ , match(parns) bind;
		if ("`parns'"!="(") local reg1 "`reg1_'";
	        tokenize "`reg1'", parse("=");
	 	if ("`2'"!="=") {; 	// if_66
	    	tokenize "`reg1'";
			local y_reg1 `1';
			macro shift;
			local x_reg1 `*';
	   	}; // if_66
	    else {;
       		if ("`4'"=="=") {;
			display as error "If more than one, equations must be enclosed in brackets";
			exit 198;
			};
       		local y_reg1 `1'; local x_reg1 `3';
       	}; // else
	}; // if_64
	else {;
		local y_reg1 `y_reg0'; local x_reg1 `x_reg0';
	}; // else
	unab x_reg1 : `x_reg1';

	tokenize "`select'", parse("=");         // define vars for regression equation
	if ("`2'"!="=") {;
		tokenize "`select'";
		local y_prob `1';
		macro shift;
		local x_prob `*';
	};
	else {;
		local y_prob `1';
		local x_prob `3';
	};
	unab x_prob : `x_prob';


 // Drop observations with missing values

	marksample touse;
	markout `touse' `y_prob' `x_prob' `cluster', strok;

	tabulate  `y_prob' if `touse';
	if (r(r) != 2) {;
		noisily display as error "`y_prob' must have exactly two distinct values";
		exit 2000;
	 };
	 else {;
		tempvar sel;
		egen byte `sel' = group(`y_prob');
		replace   `sel' = `sel'-1;            // `sel' has values 0 and 1 now
	};
	tempvar touse0 touse1;
	generate byte `touse0'=`touse';
	generate byte `touse1'=`touse';
	markout `touse0' `y_reg0' `x_reg0';
	markout `touse1' `y_reg1' `x_reg1';
	replace `touse' = cond(`sel'==0,`touse0', `touse1') if `touse';

  	if ("`weight'" == "pweight" | "`cluster'" != "") local robust robust;
  	if ~missing("`cluster'") local clopt cluster(`cluster');
 	if ~missing("`weight'")  local weight "[`weight'`exp']";


	// conditional collinearity
	forvalues i =0/1 {;
		_rmdcoll `y_reg`i'' `x_reg`i'' if `touse' & `sel'==`i', `coll';
		local result "`r(varlist)'";
		local colls`i': list x_reg`i' - result;
		if ~missing("`colls`i''") {;
			noisily display as text "note: `colls`i'' dropped from the `i' equation due to collinearity";
			local x_reg`i' `result';
		};
	}; // i

	 _rmdcoll `y_prob' `x_prob' if `touse', `coll';
		local result "`r(varlist)'";
		local colls: list x_prob - result;
		noisily display _newline;
		if ~missing("`colls'") noisily display as text "note: `colls' dropped from the selection equation due to collinearity";
		local x_prob `result';

	if ("`y_reg0'"=="`y_reg1'") {;
		local eq0 "`y_reg0'0";
		local eq1 "`y_reg0'1";
	};
	else {;
		local eq0 "`y_reg0'";   local eq1 "`y_reg1'";
	};

	// define initial values
	tempname I_reg0 I_reg1 I_prob TMP ll_0;
	tempvar	y_reg0_t y_reg1_t;

	generate double `y_reg0_t'=`y_reg0' if `sel'==0;  // define subsamples for initial value regressions
	generate double `y_reg1_t'=`y_reg1' if `sel'==1;

	scalar `ll_0' = 0;

	noisily display _newline as text "Fitting initial values " _continue;

	probit  `y_prob'   `x_prob' `weight' if `touse', nocoef;
	scalar `ll_0'=`ll_0'+e(ll); noisily display "." _continue;
 	regress `y_reg0_t' `x_reg0' `weight' if `touse', nohead;
	scalar `ll_0'=`ll_0'+e(ll); noisily display "." _continue;
 	regress `y_reg1_t' `x_reg1' `weight' if `touse', nohead;
	scalar `ll_0'=`ll_0'+e(ll); noisily display "." _continue;

	if missing("`coll'") local hopts twostep;
	else 				 local hopts iterate(20);
  	capture heckman `y_reg0_t' `x_reg0' `weight' if `touse', select(`x_prob') `hopts' `coll';
	if _rc exit _rc;
	noisily display "." _continue;
		matrix `TMP'   = e(b);
		matrix `I_reg0'= `TMP'[1,"`y_reg0_t':"];                    		// initial values for regression 1
		matrix colnames `I_reg0'= `x_reg0';
		scalar erho = e(rho);
			if (abs(erho) == 1)  scalar erho = sign(erho)*(1-0.1D-8);  		// check so rho !=1
     	local I_rho0  = atanh(erho); 			         		   			// initial values for rho0
		local I_sig0  = ln(e(sigma));                                       // initial values for sigma0
	capture heckman `y_reg1_t' `x_reg1' `weight' if `touse', select(`x_prob') `hopts' `coll';
	if _rc exit _rc;
	noisily display "." _continue;
		matrix `TMP'   = e(b);
		matrix `I_reg1'= `TMP'[1,"`y_reg1_t':"];                    		// initial values for regression 2
		matrix `I_prob'= `TMP'[1,"select:"];				        		// initial values for probit
		scalar erho = e(rho);
			if (abs(erho) == 1)  scalar erho = sign(erho)*(1-0.1D-8);  		// check so rho !=1
		local I_rho1  = atanh(erho); 			         		   			// initial values for rho1
		local I_sig1  = ln(e(sigma));                                       // initial values for sigma1
	}; // end quietly

	// maximization
	ml model d2 movestay_d2
	      	(`eq0'	   :`y_reg0' = `x_reg0')
	      	(`eq1'	   :`y_reg1' = `x_reg1')
	      	(select    :`sel'    = `x_prob')
	      	/lns0
	      	/lns1
	      	/r0
	      	/r1
	       `weight'
	      	if `touse'
	      	,
	      	title("Endogenous switching regression model")
	        search(on)
			init(`I_reg0' `I_reg1' `I_prob' `I_sig0' `I_sig1' `I_rho0' `I_rho1', copy)
			difficult
			collinear
	      	maximize
			nopreserve
			`clopt'
			`robust'
 	      	`mlopts'
	      	;

	ereturn scalar ll_0 = `ll_0';     	    	    // loglikelihood for non-correlated case
	ereturn scalar k_aux= 4;						// identify ancillary parameters
	ereturn local cmd    "movestay";             	// fill in e(cmd) ;
	ereturn local predict "movestay_p";
	ereturn local depvar "`eq0' `eq1' `y_prob'";

	if missing("`robust'") {;
		ereturn local chi2_ct "LR";
		ereturn scalar chi2_c = 2 * (e(ll) - `ll_0');
	}; // end if
	else {;
		ereturn local chi2_ct "Wald";
		quietly test [r0]_b[_cons] [r1]_b[_cons];
		ereturn scalar chi2_c = r(chi2);
	}; // end else
	ereturn scalar p_c = chiprob(2, e(chi2_c));

	movestay_replay;
end; // end program movestay_mx


