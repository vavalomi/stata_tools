*! version 3.0.2  13Apr2008 M. Lokshin, Z. Sajaia

// vague naming problem of the statistics was corrected
// yc0 gives predicted y (for the whole sample) if they were in regime 0,
// yc1 gives predicted y (for the whole sample) if they were in regime 1,
//---------------------------------------------------------------------------------------------------------
// statistics -yc1_1- -yc2_1- -yc1_2- -yc2_2- are still supported for backward compatibility
// but using new, better named -yc0- -yc1- is recommended

#delim ;
program define movestay_p;
	version 8;
	syntax anything   [if] [in] [, PSel xb0 xb1 yc0 yc1 mills0 mills1 xb2 yc1_1 yc2_1 yc1_2 yc2_2 mills2];

	marksample touse;

	syntax newvarname [if] [in] [, PSel xb0 xb1 yc0 yc1 mills0 mills1 xb2 yc1_1 yc2_1 yc1_2 yc2_2 mills2];

	if ("`e(cmd)'" ~= "movestay") error 301;
	local check "`psel' `xb0' `xb1' `yc0' `yc1' `mills0' `mills1' `xb2' `yc1_1' `yc2_1' `yc1_2' `yc2_2' `mills2'";
	if (wordcount("`check'")>1)  {;
		display as error "Only one statistic is allowed";
		exit 198;
	};
	if wordcount("`check'")==0 {;
		noisily display as text "(Option psel assumed; Pr(`y_prob'))";
		local psel "psel";
	};

	//treat statistics from the old version
	if ~missing("`xb2'")     local xb0    "xb0";
	if ~missing("`mills2'")  local mills0 "mills0";

	quietly {;

	local y_reg0: word 1 of `e(depvar)';
	local y_reg1: word 2 of `e(depvar)';
	local y_prob: word 3 of `e(depvar)';

	tempname b sigma0 sigma1 rho0 rho1 reg0 reg1 prob;
	matrix `b' = e(b);

	scalar `sigma0' =  exp(_b[lns0:_cons]);
  	scalar `sigma1' =  exp(_b[lns1:_cons]);
    scalar `rho0'   = tanh(_b[r0:_cons]);
	scalar `rho1'   = tanh(_b[r1:_cons]);
	matrix `reg0'   = `b'[1,"`y_reg0':"];
	matrix `reg1'   = `b'[1,"`y_reg1':"];
	matrix `prob'   = `b'[1,"select:"];

	tempvar xb_prob x_b0 x_b1 q sel;

	matrix score double `xb_prob' = `prob' if `touse';

	if ~missing("`psel'") {;
		generate `typlist' `varlist' = norm(`xb_prob') if `touse';
		exit;
	};

	if ~missing("`mills0'") {;
		generate `typlist' `varlist' = normden(`xb_prob')/(1-norm(`xb_prob')) if `touse';
		exit;
	};
	if ~missing("`mills1'") {;
		generate `typlist' `varlist' = normden(`xb_prob')/norm(`xb_prob')     if `touse';
		exit;
	};

	matrix score double `x_b0' = `reg0' if `touse';
	matrix score double `x_b1' = `reg1' if `touse';

	egen byte `sel' = group(`y_prob');
	replace   `sel' = `sel'-1;
	generate byte `q'   = 2*`sel'-1;

	if ~missing("`yc0'") {;
		generate `typlist' `varlist' = `x_b0' + `q'*`sigma0'*`rho0'*normden(`xb_prob')/norm(`q'*`xb_prob')  if `touse';
		exit;
	};
	if ~missing("`yc1'") {;
		generate `typlist' `varlist' = `x_b1' + `q'*`sigma1'*`rho1'*normden(`xb_prob')/norm(`q'*`xb_prob') if `touse';
		exit;
	};
	if ~missing("`xb0'") {;
		generate `typlist' `varlist' = `x_b0' if `touse';
		exit;
	};
	if ~missing("`xb1'") {;
		generate `typlist' `varlist' = `x_b1' if `touse';
		exit;
	};

	// old version again
	if ~missing("`yc1_1'") {;
		generate `typlist' `varlist' = `x_b1'+`sigma1'*`rho1'*normden(`xb_prob')/norm(`xb_prob')     if (`sel'==1) & `touse';
		exit;
	};
	if ~missing("`yc2_1'") {;
		generate `typlist' `varlist' = `x_b0'+`sigma0'*`rho0'*normden(`xb_prob')/norm(`xb_prob')     if (`sel'==1) & `touse';
		exit;
	};
	if ~missing("`yc1_2'") {;
		generate `typlist' `varlist' = `x_b1'-`sigma1'*`rho1'*normden(`xb_prob')/(1-norm(`xb_prob')) if (`sel'==0) & `touse';
		exit;
	};
	if ~missing("`yc2_2'") {;
		generate `typlist' `varlist' = `x_b0'-`sigma0'*`rho0'*normden(`xb_prob')/(1-norm(`xb_prob')) if (`sel'==0) & `touse';
		exit;
	};

	}; // end quietly

end; // end mspredict

