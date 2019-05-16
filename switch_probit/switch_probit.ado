*! version 1.0 14Jan2009 M. Lokshin, Z. Sajaia

# delimit ;

program define switch_probit, properties(ml_score);
	version 8.1;
	if replay() {;
           if ("`e(cmd)'" != "switch_probit") error 301;
           switch_probit_replay `0';
	};
	else switch_probit_mx `0';

end; // program switch_probit

program define switch_probit_replay;
	syntax [,Level(cilevel)];

    local rho1 diparm(athrho1, tanh label("rho1"));
    local rho0 diparm(athrho0, tanh label("rho0"));

 	_coef_table_header;
	display;
   	_coef_table, level(`level') `rho1' `rho0' notest;

	if "`e(vcetype)'" != "Robust" local testtyp LR;
	else 						  local testtyp Wald;
	display as text "`testtyp' test of indep. eqns. (rho1=rho0=0):"
		 			_col(38) "chi2(" as result "2" as text ") = "
		 			as result %8.2f e(chi2_c)
					_col(59) as text "Prob > chi2 = " as result %6.4f e(p_c);
	display in smcl in green "{hline 78}";
	_prefix_footnote;
	exit e(rc);
end; // program switch_probit_replay

program define switch_probit_mx, eclass;

 	syntax 	anything(id="equations id" equalok) [pweight iweight fweight]
	        [if] [in], SELect(string) [CLuster(varname) Robust noSKIP noCONstant noLOg CRITTYPE(passthru)
									   offset_s(varname) offset1(varname) offset0(varname) *];
	mlopts mlopts, `options';
	local coll `s(collinear)';

	// read vector of parameters

	gettoken prob1 prob0_ :  anything,  match(parns) bind;
	if ("`parns'"!="(") {;
           local prob1 "`anything'";
           local prob0_;
	};
	tokenize "`prob1'", parse("()");
	if ("`1'"!="`prob1'") {;
           display as error "If more than one, equations must be enclosed in brackets";
           exit 198;
	};

	tokenize "`prob1'", parse("=");
	if ("`2'"!="=") {;
           tokenize "`prob1'";
           local y_prob1 `1';
           macro shift;
           local x_prob1 `*';
	};
	else {; // else_63
		if ("`4'"=="=") {;
                   display as error "If more than one, equations must be enclosed in brackets";
                   exit 198;
		};
	       	local y_prob1 `1';
			local x_prob1 `3';
        }; // else_63
	unab x_prob1 : `x_prob1', min(0);
	local l_y1 "`y_prob1'";

	if ("`prob0_'"!="") {; // if_64
           gettoken prob0  : prob0_ , match(parns) bind;
           if ("`parns'"!="(") local prob0 "`prob0_'";
           tokenize "`prob0'", parse("=");
           if ("`2'"!="=") {;   // if_66
              tokenize "`prob0'";
              local y_prob0 `1';
              macro shift;
              local x_prob0 `*';
           }; // if_66
           else {;
       		if ("`4'"=="=") {;
                   display as error "If more than one, equations must be enclosed in brackets";
                   exit 198;
                   };
          local y_prob0 `1'; local x_prob0 `3';
       	}; // else
		local l_y0 "`y_prob0'";
	}; // if_64
	else {;
			local y_prob0 `y_prob1';
			local x_prob0 `x_prob1';

			local l_y1 "`y_prob1'_1";
			local l_y0 "`y_prob1'_0";
	}; // else
	unab x_prob0 : `x_prob0', min(0);

	tokenize "`select'", parse("=");         // define vars for regression equation
	if ("`2'"!="=") {;
           tokenize "`select'";
           local y_sel `1';
           macro shift;
           local x_sel `*';
	};
	else {;
           local y_sel `1';
           local x_sel `3';
	};
	unab x_sel : `x_sel', name("selection equation");


  	if ("`weight'" == "pweight" | "`cluster'" != "") {;
		local robust "robust";
		if missing(`"`crittyp'"') local crtype crittype("log pseudolikelihood");
	};
	else local crtype `crittyp';
	if ~missing("`cluster'")  local clopt "cluster(`cluster')";
	if ~missing("`weight'")   local wgt `"[`weight'`exp']"';
	if ~missing("`log'")      local quietly "quietly";
	if ~missing("`offset_s'") local offo_s "offset(`offset_s')";
	if ~missing("`offset1'")  local offo1  "offset(`offset1')";
    if ~missing("`offset0'")  local offo0  "offset(`offset0')";

 // Drop observations with missing values
	marksample touse;
	markout `touse' `y_sel' `x_sel' `cluster' `offset_s', strok;

quietly {;
	tabulate  `y_sel' if `touse';
	if (r(r) != 2) {;
		noisily display as error "`y_sel' must have exactly two distinct values";
		exit 2000;
	};

	tempvar touse1 touse0;
	generate byte `touse1'=`touse';
	generate byte `touse0'=`touse';
	markout `touse1' `y_prob1' `x_prob1' `offset1';
	markout `touse0' `y_prob0' `x_prob0' `offset0';
	replace `touse' = cond(`y_sel'==1,`touse1', `touse0') if `touse';

	// conditional collinearity
	forvalues i =0/1 {;
		tabulate  `y_prob`i'' if `touse`i'';
		if (r(r) != 2) {;
			noisily display as error "`y_prob`i'' must have exactly two distinct values";
			exit 2000;
		};

		_rmdcoll `y_prob`i'' `x_prob`i'' if `touse' & `y_sel'==`i', `coll';
		local result "`r(varlist)'";
		local colls`i': list x_prob`i' - result;
		if ~missing("`colls`i''") {;
			noisily display as text "note: `colls`i'' dropped from the `i' equation due to collinearity";
			local x_prob`i' `result';
		};
	}; // i
	_rmdcoll `y_sel' `x_sel' if `touse', `coll';
	local result "`r(varlist)'";
	local colls: list x_sel - result;
	noisily display _newline;
	if ~missing("`colls'") noisily display as text "note: `colls' dropped from the selection equation due to collinearity";
	local x_sel `result';

}; // quietly
	tempname llc b0 TMP b0sel b00;
	tempvar  nshaz1 nshaz0;

	`quietly' display in green _n "Fitting probit model for `y_sel'=1:";
	`quietly' probit `y_prob1' `x_prob1' `wgt' if `touse' & `y_sel'==1, `constant' nocoef `crtype' `offo1';
	scalar `llc' = e(ll);

	`quietly' display in green _n "Fitting probit model for `y_sel'=0:";
	`quietly' probit `y_prob0' `x_prob0' `wgt' if `touse' & `y_sel'==0, `constant' nocoef `crtype' `offo0';
	scalar `llc' = `llc' + e(ll);

	`quietly' display in green _n "Fitting selection model:";
	`quietly' probit `y_sel' `x_sel' `wgt' if `touse', asis nocoef `crtype' `offo_s';
	matrix `b0sel' = e(b);
	if missing("`robust'") {;
		scalar `llc' = `llc' + e(ll);
		`quietly' display in green _n "Comparison:    log likelihood = " in yellow %10.0g `llc';
	};

	quietly predict  double `nshaz1', xb;
	quietly generate double `nshaz0' = normd(`nshaz1') / (1-normprob(`nshaz1'));  // sel = 0;
	quietly replace         `nshaz1' = normd(`nshaz1') / normprob(`nshaz1'); 	  // sel = 1;

	if missing("`constant'") {;
		tempname one;
		generate byte `one' = 1;
	};

	if ~missing("`skip'") {;
		`quietly' display in green _n "Fitting constant-only starting values:";
		`quietly' probit `y_prob1' `one' `nshaz1' `wgt' if `touse' & `y_sel'==1, noconstant nocoef asis `crtype' `offo1';

		matrix `TMP' = e(b);
		local k = colsof(`TMP');

		tempname rho athrho0 athrho1;

		scalar `rho'     = _b[`nshaz1'];
		scalar `rho'     = max(min(`rho',.85), -.85);
		scalar `athrho1' = atanh(`rho');

		matrix  `b00' = `b0sel' , `TMP'[1,1..`k'-1];

		`quietly' probit `y_prob0' `one' `nshaz0' `wgt' if `touse' & `y_sel'==0, noconstant nocoef asis `crtype' `offo0';

		matrix `TMP' = e(b);
		local k = colsof(`TMP');

		scalar `rho'     = _b[`nshaz0'];
		scalar `rho'     = max(min(`rho',.85), -.85);
		scalar `athrho0' = -atanh(`rho');

		matrix  `b00' = `b00' , `TMP'[1,1..`k'-1], `athrho1', `athrho0';
		local from0 "`b00', copy";
	};

	`quietly' display in green _n "Fitting starting values:";
	`quietly' probit `y_prob1' `x_prob1' `one' `nshaz1' `wgt' if `touse' & `y_sel'==1, noconstant nocoef asis `crtype' `offo1';

	matrix `TMP' = e(b);
	local k = colsof(`TMP');

	tempname rho athrho1 athrho0;

	scalar `rho'     = _b[`nshaz1'];
	scalar `rho'     = max(min(`rho',.85), -.85);
	scalar `athrho1' = atanh(`rho');

	matrix  `b0' = `b0sel' , `TMP'[1,1..`k'-1];

	`quietly' probit `y_prob0' `x_prob0' `one' `nshaz0' `wgt' if `touse' & `y_sel'==0, noconstant nocoef asis `crtype' `offo0';

	matrix `TMP' = e(b);
	local k = colsof(`TMP');

	scalar `rho'     = _b[`nshaz0'];
	scalar `rho'     = max(min(`rho',.85), -.85);
	scalar `athrho0' = -atanh(`rho');

	matrix  `b0' = `b0' , `TMP'[1,1..`k'-1], `athrho1', `athrho0';
  	local from "`b0', copy";

	if ~missing("`skip'") {;
		`quietly' display in green _n "Fitting constant-only model:";

		capture noisily ml model d2 switch_probit_d2
					(`y_sel' : `y_sel'   = `x_sel', `offo_s')
            		(`l_y1'  : `y_prob1' =, `offo1')
	 				(`l_y0'  : `y_prob0' =, `offo0')
            		/athrho1
            		/athrho0
					`wgt' if `touse' , waldtest(0)
					 collinear missing maximize nooutput nopreserve
					 init(`from0') search(off) `log' `mlopts' `crittyp'
					 nocnsnotes;
		if _rc == 1400 & "`from'" == "`b0', copy" {;
			display as text "note:  default initial values infeasible; starting from B=0";

			ml model d2 switch_probit_d2
					(`y_sel' : `y_sel'   = `x_sel', `offo_s')
            		(`l_y1'  : `y_prob1' =, `offo1')
	 				(`l_y0'  : `y_prob0' =, `offo0')
            		/athrho1
            		/athrho0
					`wgt' if `touse' , waldtest(0)
					collinear missing maximize nooutput nopreserve difficult
					init(/athrho1=0 /athrho0=0) search(off) `log' `mlopts' `crittyp' nocnsnotes;
		};
		else if _rc {;
				error _rc;
			};
		local continu "continue";
	};

	global counter=0;

	`quietly' display in green _n "Fitting full model:";

	capture noisily ml model d2 switch_probit_d2
	   		(`y_sel' : `y_sel'    = `x_sel', `offo_s')
            (`l_y1'  : `y_prob1' = `x_prob1', `offo1')
	 		(`l_y0'  : `y_prob0' = `x_prob0', `offo0')
            /athrho1
            /athrho0
		 	if `touse' `wgt',
			collinear missing maximize nooutput nopreserve difficult
			title("Switching probit model")
			init(`from') search(off) `continu' `log' `mlopts' `crittyp' `robust' `clopt';

	if _rc == 1400 & "`from'" == "`b0', copy" {;
		display as text "note:  default initial values infeasible; starting from B=0";

        ml model d2 switch_probit_d2
	   		(`y_sel' : `y_sel'   = `x_sel', `offo_s')
            (`l_y1'  : `y_prob1' = `x_prob1', `offo1')
	 		(`l_y0'  : `y_prob0' = `x_prob0', `offo0')
            /athrho1
            /athrho0
			`wgt'
			if `touse'
			,
			title("Switching probit model")
	 		difficult collinear	missing maximize nooutput nopreserve
        	init(/athrho1=0 /athrho0=0) search(off) `clopt' `continu' `robust' `mlopts' `log';

	};
	else if _rc error _rc;

	ereturn scalar k_aux = 2;

	if missing("`robust'") {; // test of independent equations
		ereturn scalar ll_c = `llc';
		ereturn scalar chi2_c = abs(-2*(e(ll)-e(ll_c)));
		ereturn local chi2_ct "LR";
	};
	else {;
		quietly test [athrho1]_cons = [athrho0]_cons = 0;
		ereturn scalar chi2_c = r(chi2);
		ereturn local chi2_ct "Wald";
	};
	ereturn scalar p_c = chiprob(2, e(chi2_c));

	ereturn scalar rho1 = tanh([athrho1]_b[_cons]);
	ereturn scalar rho0 = tanh([athrho0]_b[_cons]);

	ereturn local marginsok "zb xb1 xb0 p11 p01 p10 p00 psel pcond1 pcond0 te";
	ereturn local marginsnotok "stdp1 stdp0 mte"; // tt tu
	ereturn local depvar "`y_sel' `y_prob1' `y_prob0'";
	ereturn local predict "switch_probit_p";
	ereturn local cmd "switch_probit";

	switch_probit_replay;

end; // program switch_probit_mx

