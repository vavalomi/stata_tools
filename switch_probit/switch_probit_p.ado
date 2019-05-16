#delimit ;

program define switch_probit_p;
	version 8.1;

	quietly {;
		//stolen from bipr_p.ado
		syntax [anything] [if] [in] [, SCores * ];
		if ~missing(`"`scores'"') {;
			ml_score `0';
			exit;
		};

		local myopts ZB XB1 XB0 P11 P01 P10 P00 PCOND1 PCOND0 PSEL STDP1 STDP0 TT TE TU MTE;

		_pred_se "`myopts'" `0';
		if (`s(done)') exit;
		local vtyp  `s(typ)';
		local varn `s(varn)';
		local 0 `"`s(rest)'"';

		syntax [if] [in] [, `myopts' noOFFset];

		local type `zb'`xb1'`xb0'`p11'`p10'`p01'`p00'`pcond1'`pcond0'`psel'`stdp1'`stdp0'`tt'`te'`tu' `mte';

		marksample touse;

		local y_sel : word 1 of `e(depvar)';
		local y_1   : word 2 of `e(depvar)';
		local y_0   : word 3 of `e(depvar)';

		if e(df_m) != 0 {;
			tempvar tzb txb1 txb0;
			_predict double `tzb'  if `touse', eq(#1) `offset';
			_predict double `txb1' if `touse', eq(#2) `offset';
			_predict double `txb0' if `touse', eq(#3) `offset';
		};
		else {;
			if (~missing("`e(offset1)'") & missing("`offset1'")) local txb1 `e(offset1)';
			else local txb1 0;
			if (~missing("`e(offset0)'") & missing("`offset0'")) local txb0 `e(offset0)';
			else local txb0 0;
		}; // if
		tempvar psel;
		generate double `psel' = normal(`tzb');
		replace `psel' = epsdouble() if `psel' == 0;
		replace `psel' = 1 - epsdouble() if `psel' == 1;

		tempname rho1 rho0;
		scalar `rho1' = tanh([athrho1]_b[_cons]);
		scalar `rho0' = tanh([athrho0]_b[_cons]);

		if (missing("`type'") | "`type'" == "p11") {;
			if missing("`type'") noisily display as text "(option p11 assumed; Pr(`y_sel'=1, `y_1'=1)";

			generate `vtyp' `varn' = binorm(`tzb',`txb1', `rho1');
			label variable `varn' "Pr(`y_sel'=1,`y_1'=1)";
			exit;
		};
		if "`type'"=="p10" {;
			generate `vtyp' `varn' = binorm(`tzb',-`txb1', -`rho1');
			label variable `varn' "Pr(`y_sel'=1,`y_1'=0)";
			exit;
		};

		if "`type'"=="p01" {;
			generate `vtyp' `varn' = binorm(-`tzb', `txb0', -`rho0');
			label variable `varn' "Pr(`y_sel'=0,`y_0'=1)";
			exit;
		};
		if "`type'"=="p00" {;
			generate `vtyp' `varn' = binorm(-`tzb',-`txb0', `rho0');
			label variable `varn' "Pr(`y_sel'=0,`y_0'=0)";
			exit;
		};
		if "`type'" == "pcond1" {;
			generate `vtyp' `varn' = binorm(`tzb', `txb1', `rho1')/ `psel';
			label variable `varn' "Pr(`y_1'=1 | y_sel'=1)";
			exit;
		};
		if "`type'" == "pcond0" {;
			generate `vtyp' `varn' = binorm(-`tzb', `txb0', -`rho0')/(1-`psel');
			label variable `varn' "Pr(`y_0'=1 | y_sel'=0)";
			exit;
		};
		if "`type'" == "psel" {;
			generate `vtyp' `varn' = `psel';
			label variable `varn' "Pr(`y_sel'=1)";
			exit;
		};
		if ("`type'" == "zb") {;
			generate `vtyp' `varn' = `tzb';
			label variable `varn' "linear prediction of `y_sel'";
			exit;
		};
		if ("`type'" == "xb1") {;
			generate `vtyp' `varn' = `txb1';
			label variable `varn' "linear prediction of `y_1'";
			exit;
		};
		if ("`type'" == "xb0") {;
			generate `vtyp' `varn' = `txb0';
			label variable `varn' "linear prediction of `y_0'";
			exit;
		};
		if ("`type'" == "stdpsel") {;
			_predict `vtyp' `varn', stdp eq(#1) `offset', if `touse';
			label variable `varn' "S.E. of prediction of `y_sel'";
			exit;
		};
		if ("`type'" == "stdp1") {;
			_predict `vtyp' `varn', stdp eq(#2) `offset', if `touse';
			label variable `varn' "S.E. of prediction of `y_1'";
			exit;
		};
		if ("`type'" == "stdp0") {;
			_predict `vtyp' `varn', stdp eq(#3) `offset', if `touse';
			label variable `varn' "S.E. of prediction of `y_0'";
			exit;
		};
		if ("`type'" == "tt") {;
			generate `vtyp' `varn' = (binorm(`tzb',`txb1', `rho1')-binorm(`tzb',`txb0', `rho0'))/`psel' if `touse' & `y_sel' == 1;
			label variable `varn' "tt";
			exit;
		};
		if ("`type'" == "tu") {;
			generate `vtyp' `varn' = (binorm(-`tzb',`txb1',-`rho1')-binorm(-`tzb',`txb0',-`rho0'))/(1-`psel') if `touse' & `y_sel' == 0;
			label variable `varn' "tu";
			exit;
		};
		if ("`type'" == "te") {;
			generate `vtyp' `varn' = normal(`txb1')-normal(`txb0') if `touse';
			label variable `varn' "te";
			exit;
		};
		if ("`type'" == "mte") {;
			local select : word 1 of `=e(depvar)';
			summarize `txb1' if (`y_sel' == 1), meanonly; local mean_x1 = r(mean);
			summarize `txb0' if (`y_sel' == 0), meanonly; local mean_x0 = r(mean);
			generate ___grid = _n/100 in 1/100;

			generate `vtyp' `varn' = normal((`mean_x1'+`rho1'*(-___grid))/sqrt(1-(`rho1')^2) ) -
						normal((`mean_x0'+(`rho0')*(-___grid))/sqrt(1-(`rho0')^2) )
						if `touse' in 1/100;
			label variable `varn' "mte";
			exit;
		};
		display as error "unknown statistic";
		exit 198;
	}; // quietly
end;
