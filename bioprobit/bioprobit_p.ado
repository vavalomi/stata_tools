*! version 1.0  27Aug2006 Z. Sajaia
#delim ;

program define bioprobit_p;
	version 9;

	quietly {;
	//stolen from bipr_p.ado
	syntax [anything] [if] [in] [, SCores * ];
	if ~missing(`"`scores'"') {;
		ml_score `0';
		exit;
	};

   	local myopts XB1 XB2 STDP1 STDP2 Outcome(string);

	_pred_se "`myopts'" `0';
	if (`s(done)') exit;
	local vtyp  `s(typ)';
	local varn `s(varn)';
	local 0 `"`s(rest)'"';

	syntax [if] [in] [, `myopts' noOFFset];

	local type `xb1'`xb2'`stdp1'`stdp2';

	marksample touse;

	local y1 : word 1 of `e(depvar)';
	local y2 : word 2 of `e(depvar)';

 	if (missing("`type'") & missing("`outcome'")) {;
		local outcome "#1, #1";
		noisily display as text "(option outcome(#1, #1) assumed)";
	};

	if e(df_m) != 0 {;
		tempvar txb1 txb2;
		_predict double `txb1' if `touse', eq(#1) `offset';
		_predict double `txb2' if `touse', eq(#2) `offset';
	};
	else {;
		 if (~missing("`e(offset1)'") & missing("`offset1'")) local txb1 `e(offset1)';
		 else local txb1 0;
		 if (~missing("`e(offset2)'") & missing("`offset1'")) local txb2 `e(offset1)';
		 else local txb2 0;
	}; // if

	if ~missing("`outcome'") {;
		// get cuttoff points
		_parse comma out1 out2: outcome;
		local out2 : subinstr local out2 "," "", all;
		local out2 = trim("`out2'");
		forvalues i = 1/2 {;
			tempvar sy;
			egen `sy' = group(`y`i'');
			capture confirm numeric variable `y`i'';	scalar rc = _rc;
			if substr(`"`out`i''"',1,1)=="#" {;
				local out = substr(`"`out`i''"',2,.);
				Chk confirm integer number `out';
				Chk assert `out' >= 1;
				quietly tabulate `y`i'' if `touse';
				if (`out'>r(r)) {;
					noisily display as error "there is no outcome #`out'" _n
											 "there are only `r(r)' categories in `y`i''";
					exit 111;
				};
				else local c =`out';
			};
			else {;
				if (rc) local out`i' "`"`out`i''"'";
				if missing(`out`i'') local c=-99;
				else {;
					Chk confirm number `out`i'';
					summarize `sy' if `touse' & `y`i''==`out`i'', meanonly; local c = r(mean);
					if missing(`c') {;
					    noisily display as error "`y`i'' never takes value `out`i''";
						exit 111;
					};
				};
 			};
			if (`c' == -99) {;
				local cut`i'   = maxfloat();
				local cut`i'_1 = minfloat();
			};
			else {;
   				quietly summarize `sy' if `touse';
				if (`c'==r(max)) local cut`i' = maxfloat();
				else             local cut`i' = _b[cut`i'`c':_cons];
				local --c;
				if (`c'==0) local cut`i'_1 = minfloat();
				else 		local cut`i'_1 = _b[cut`i'`c':_cons];
			};
		}; // end forvalues

		tempname gamma rr rho;

		capture scalar `gamma' = _b[gamma:_cons];
		if _rc  scalar `gamma' = 0; // it's a SUR model.

		scalar `rr'  = 1/sqrt(1+(`gamma')^2+2*`gamma'*tanh(_b[athrho:_cons]));
		scalar `rho' = (tanh(_b[athrho:_cons])+`gamma')*`rr';

  		generate `vtyp' `varn' =  binorm(`cut1'  -`txb1',`rr'*(`cut2'  -`txb2'-`gamma'*`txb1'),`rho')
							   	- binorm(`cut1_1'-`txb1',`rr'*(`cut2'  -`txb2'-`gamma'*`txb1'),`rho')
							  	- binorm(`cut1'  -`txb1',`rr'*(`cut2_1'-`txb2'-`gamma'*`txb1'),`rho')
							  	+ binorm(`cut1_1'-`txb1',`rr'*(`cut2_1'-`txb2'-`gamma'*`txb1'),`rho');
		exit;
	}; // outcome

	if ("`type'" == "xb1") {;
		generate `vtyp' `varn' = `txb1';
		label variable `varn' "linear prediction of `y1'";
		exit;
	};
	if ("`type'" == "xb2") {;
		generate `vtyp' `varn' = `txb2';
		label variable `varn' "linear prediction of `y2'";
		exit;
	};
	if ("`type'" == "stdp1") {;
		_predict `vtyp' `varn', stdp eq(#1) `offset', if `touse';
		label variable `varn' "S.E. of prediction of `y1'";
		exit;
	};
	if ("`type'" == "stdp2") {;
		_predict `vtyp' `varn', stdp eq(#2) `offset', if `touse';
		label variable `varn' "S.E. of prediction of `dy2'";
		exit;
	};
	error 198;
	}; // quietly
end;

//stolen from ologit_p.ado
program define Chk;
	capture `0';
	if _rc {;
		noisily display as error "outcome() must be pair of either values of `e(depvar)',"
						 _n "or #1, #2, ...";
		exit 111;
	};
end;

