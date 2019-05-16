*! version 1.11  17Aug2007
#delimit ;

program define bioprobit_d2, eclass;

	local NC1   = $NC1;
	local NC1_1 = `NC1'-1;
	local NC2   = $NC2;
	local NC2_1 = `NC2'-1;
	local K1    = 4+$END;
	local K2    = `K1'-1+`NC1';
	local N     = `K1'-1+`NC1_1'+`NC2_1';

 	forvalues i = 1/`N' {; local grad "`grad' g`i'"; };

	args todo b lnf g negH `grad';
	quietly {;

	tempvar p1_b p2_b lnf1;
	tempname r rho gamma d dtanh_dr drho_dr drho_dgamma rr drr1_dr d2tanh_dr2;

	mleval `p1_b'  = `b' , eq(1);
	mleval `p2_b'  = `b' , eq(2);
	mleval `r'     = `b' , eq(3) scalar;
	if $END mleval `gamma' = `b' , eq(4) scalar;
	else    scalar `gamma' = 0;
 	forvalues k = 1/`NC1_1' {;
		tempname c`k'1;
		mleval `c`k'1' = `b' , eq(`=`K1'-1+`k'') scalar;
	};
	forvalues k = 1/`NC2_1' {;
		tempname c`k'2;
		mleval `c`k'2' = `b' , eq(`=`K1'-1+`NC1_1'+`k'') scalar;
	};

	scalar `rr'  = 1/sqrt(1+(`gamma')^2+2*`gamma'*tanh(`r'));
	scalar `rho' = (tanh(`r')+`gamma')*`rr';

	tempvar cut1 cut1_1 cut2 cut2_1;
	tempname CUT1 CUT2;

	local neginf = minfloat();
	local posinf = maxfloat();
	if (`NC1_1' > 1) matrix `CUT1' = `neginf' \ `c11' \ J(`NC1_1'-1,1,.) \ `posinf';
	else             matrix `CUT1' = `neginf' \ `c11'                    \ `posinf';
	if (`NC2_1' > 1) matrix `CUT2' = `neginf' \ `c12' \ J(`NC2_1'-1,1,.) \ `posinf';
	else	         matrix `CUT2' = `neginf' \ `c12'                    \ `posinf';
	forvalues k = 2/`NC1_1' {;
		matrix `CUT1'[`k'+1,1] = `CUT1'[`k',1]+`c`k'1'^2;
	};
	forvalues k = 2/`NC2_1' {;
		matrix `CUT2'[`k'+1,1] = `CUT2'[`k',1]+`c`k'2'^2;
	};

  	generate double `cut1'   = `CUT1'[$ML_y1+1,1];
	generate double `cut1_1' = `CUT1'[$ML_y1  ,1];
	generate double `cut2'   = `CUT2'[$ML_y2+1,1];
	generate double `cut2_1' = `CUT2'[$ML_y2  ,1];

	local arg11 "(`cut1'  -`p1_b')"; local arg21 "(`rr'*(`cut2'  -`p2_b'-`gamma'*`p1_b'))"; local cond11 "$ML_y1";   local cond21 "$ML_y2";
	local arg12 "(`cut1_1'-`p1_b')"; local arg22 "(`rr'*(`cut2'  -`p2_b'-`gamma'*`p1_b'))"; local cond12 "$ML_y1-1"; local cond22 "$ML_y2";
	local arg13 "(`cut1'  -`p1_b')"; local arg23 "(`rr'*(`cut2_1'-`p2_b'-`gamma'*`p1_b'))"; local cond13 "$ML_y1";   local cond23 "$ML_y2-1";
	local arg14 "(`cut1_1'-`p1_b')"; local arg24 "(`rr'*(`cut2_1'-`p2_b'-`gamma'*`p1_b'))"; local cond14 "$ML_y1-1"; local cond24 "$ML_y2-1";

	forvalues i = 1/4 {;
		tempvar q`i';
		generate double `q`i''= binorm(`arg1`i'',`arg2`i'',`rho');
	};

	generate double `lnf1' = `q1'-`q2'-`q3'+`q4';
	replace 		`lnf1' = 0.1D-60 if(`lnf1'<=0);
	mlsum `lnf' = log(`lnf1');

  	if (`todo'==0 | `lnf'==.)  exit;                                        // calculate first derivatives

   	scalar `d'           = 1/sqrt(1-(`rho')^2);
	scalar `dtanh_dr'    = 4*exp(2*`r')/((1+exp(2*`r'))^2);
	scalar `drho_dr'     = (`rr'-`gamma'*`rho'*(`rr')^2)*`dtanh_dr';
	scalar `drr1_dr'     = -`gamma'*(`rr')^2*`dtanh_dr';
	scalar `drho_dgamma' = `rr'*(1-(`rho')^2);

	forvalues i = 1/4 {;
		tempvar dF1`i' dF2`i' dF3`i' d1q`i' d3q`i' d`K2'q`i';

		generate double `dF1`i''   = normden(`arg1`i'')*norm(`d'*(`arg2`i'' - `rho'*`arg1`i''));
		generate double `dF2`i''   = normden(`arg2`i'')*norm(`d'*(`arg1`i'' - `rho'*`arg2`i''));
	   	generate double `dF3`i''   = exp(-0.5*((`arg1`i'')^2+(`arg2`i'')^2-2*`rho'*`arg1`i''*`arg2`i'')*`d'^2)/(2*_pi*sqrt(1-(`rho')^2));

		generate double `d1q`i''   = - `dF1`i'' - `gamma'*`rr'*`dF2`i'';
		local 			 d2q`i'    = "(-`dF2`i''*`rr')";
		generate double `d3q`i''   = `dF2`i''*`arg2`i''*`drr1_dr'+`dF3`i''*`drho_dr';
		if $END {;
			tempvar d4q`i';
			generate double `d4q`i'' = -(`rho'*`arg2`i''+`p1_b')*`rr'*`dF2`i''+`dF3`i''*`drho_dgamma';
		};
		local 			 d`K1'q`i'  = "`dF1`i''";
		generate double `d`K2'q`i'' =  `rr'*`dF2`i'';
 		forvalues k = 2/`NC1_1' {;
			local eqn = `K1'-1+`k';
			tempvar d`eqn'q`i';
			generate double `d`eqn'q`i'' = cond(`cond1`i''>=`k',      `dF1`i''*2*`c`k'1', 0);
		}; // k
 		forvalues k = 2/`NC2_1' {;
			local eqn = `K1'-1+`NC1_1'+`k';
			tempvar d`eqn'q`i';
			generate double `d`eqn'q`i'' = cond(`cond2`i''>=`k', `rr'*`dF2`i''*2*`c`k'2', 0);
		}; // k
    }; // i

	capture matrix drop `g';
	forvalues i = 1/`N'{;
		tempvar g`i'i;
		generate double `g`i'i' = `d`i'q1'-`d`i'q2'-`d`i'q3'+`d`i'q4';
		 replace        `g`i''  = `g`i'i'/`lnf1' ;
		mlvecsum `lnf'  `g`i''  = `g`i'', eq(`i');
		matrix `g' = (nullmat(`g'), `g`i'');
	}; // i

	if (`todo'==1 | `lnf'==. )  exit;                            	// calculate second derivative

    scalar `d2tanh_dr2'   = 8*exp(2*`r')*(1-exp(2*`r'))/(1+exp(2*`r'))^3;

	forvalues i = 1/4 {;
		tempvar d1_1q`i' d1_2q`i' d1_3q`i' d1_4q`i' d1_5q`i'
						 d2_2q`i' d2_3q`i' d2_4q`i' d2_5q`i'
								  d3_3q`i' d3_4q`i' d3_5q`i'
										   d4_4q`i' d4_5q`i'
													d5_5q`i';

		tempvar dF1_1`i' dF1_3`i' dF2_2`i' dF2_3`i' dF3_3`i';
		generate double `dF1_1`i''=-(`arg1`i''*`dF1`i''+`rho'*`dF3`i'');
		local  			 dF1_2`i' ="`dF3`i''";
 		generate double `dF1_3`i''=`dF3`i''*(`rho'*`arg2`i''-`arg1`i'')*`d'^2;

		generate double `dF2_2`i''=-(`arg2`i''*`dF2`i''+`rho'*`dF3`i'');
		generate double `dF2_3`i''=`dF3`i''*(`rho'*`arg1`i''-`arg2`i'')*`d'^2;

   		generate double `dF3_3`i''=`dF3`i''*(`arg1`i''*`arg2`i''+`rho'-`rho'*(`arg1`i''^2+`arg2`i''^2-2*`rho'*`arg1`i''*`arg2`i'')*`d'^2)*`d'^2;

		tempvar dA2_dgamma;
		generate double `dA2_dgamma' = -`rr'*(`rho'*`arg2`i''+`p1_b');


		generate double `d1_1q`i''= -(`dF1_1`i''+2*`rr'*`gamma'*`dF1_2`i''+(`gamma'*`rr')^2*`dF2_2`i'');
		generate double `d1_2q`i''= -(`dF1_2`i''+  `rr'*`gamma'*`dF2_2`i'')*`rr';
		generate double `d1_3q`i''=  (`dF1_2`i''+  `rr'*`gamma'*`dF2_2`i'')*`arg2`i''*`drr1_dr' +
									 (`dF1_3`i''+  `rr'*`gamma'*`dF2_3`i'')*`drho_dr'+
								                   `rr'*`gamma'*`dF2`i''*`drr1_dr';
if $END	generate double `d1_4q`i''=  (`dF1_2`i''+  `rr'*`gamma'*`dF2_2`i'')*`dA2_dgamma'+
									 (`dF1_3`i''+  `rr'*`gamma'*`dF2_3`i'')*`drho_dgamma'+
									 (1-`rho'*`rr'*`gamma')*`dF2`i''*`rr';

		generate double `d1_`K1'q`i''=   `dF1_1`i''+  `rr'*`gamma'*`dF1_2`i'';

		generate double `d2_2q`i''= -(`rr')^2*`dF2_2`i'';

		generate double `d2_3q`i''   = (`dF2_2`i''*`arg2`i''*`drr1_dr'+`dF2_3`i''*`drho_dr'+	 `dF2`i''*`drr1_dr')*`rr';
if $END	generate double `d2_4q`i''   = (`dF2_2`i''*`dA2_dgamma'+       `dF2_3`i''*`drho_dgamma'-`dF2`i''*`rr'*`rho')*`rr';
		generate double `d2_`K1'q`i''= `rr'*`dF1_2`i'';

		generate double `d3_3q`i''   =-((`dF2_2`i''*`arg2`i''*`drr1_dr'+`dF2_3`i''*`drho_dr')*`arg2`i''*`drr1_dr'+
									    (`dF2_3`i''*`arg2`i''*`drr1_dr'+`dF3_3`i''*`drho_dr')*`drho_dr'+
									    3*`arg2`i''*(`drr1_dr')^2*`dF2`i''+
									    (`rho'*`drr1_dr'+2*`drho_dr')*`drr1_dr'*`dF3`i'');

if $END	generate double `d3_4q`i''=-((`dF2_2`i''*`arg2`i''*`drr1_dr'+`dF2_3`i''*`drho_dr')*`dA2_dgamma'+
								     (`dF2_3`i''*`arg2`i''*`drr1_dr'+`dF3_3`i''*`drho_dr')*`drho_dgamma'+
									`dF2`i''*(`drr1_dr'*`dA2_dgamma'-`arg2`i''*`rr'*`rho'*`drr1_dr'-`arg2`i''*`rr'*`drho_dr')+
									`dF3`i''*(`drr1_dr'*`drho_dgamma'-2*`rho'*`rr'*`drho_dr'));

		generate double `d3_`K1'q`i''= -(`dF1_2`i''*`arg2`i''*`drr1_dr'+`dF1_3`i''*`drho_dr');

if $END	generate double `d4_4q`i''= -((`dF2_2`i''*`dA2_dgamma'+`dF2_3`i''*`drho_dgamma')*`dA2_dgamma'+
									  (`dF2_3`i''*`dA2_dgamma'+`dF3_3`i''*`drho_dgamma')*`drho_dgamma'-
									  `dF2`i''*`rr'*(2*`rho'*`dA2_dgamma'+`arg2`i''*`drho_dgamma')-
									  `dF3`i''*3*`rr'*`rho'*`drho_dgamma');
if $END	generate double `d4_5q`i''=-(`dF1_2`i''*`dA2_dgamma'+`dF1_3`i''*`drho_dgamma');

		generate double `d`K1'_`K1'q`i''=- `dF1_1`i'';

		forvalues k = 2/`NC1_1' {;
	   		local eqn = `K1'-1+`k';
			forvalues l = 1/`K1' {;
				tempvar          d`l'_`eqn'q`i';
				generate double `d`l'_`eqn'q`i'' = cond(`cond1`i''>=`k', `d`l'_`K1'q`i''*2*`c`k'1', 0);
			};

			forvalues kk = 2/`NC2_1' {;
				local eqn2 = `K1'-1+`NC1_1'+`kk';
				tempvar          d`eqn'_`eqn2'q`i';
				generate double `d`eqn'_`eqn2'q`i'' = cond((`cond1`i''>=`k'   & `cond2`i''>=`kk'), -`d2_`K1'q`i''*4*`c`k'1'*`c`kk'2', 0);
			}; // kk
		}; // k

		forvalues k = 2/`NC2_1' {;
	   		local eqn = `K1'-1+`NC1_1'+`k';
			tempvar d1_`eqn'q`i' d2_`eqn'q`i' d3_`eqn'q`i' d`K1'_`eqn'q`i';
	   		generate double `d1_`eqn'q`i'' = cond(`cond2`i''>=`k', -`d1_2q`i''*2*`c`k'2', 0);
			generate double `d2_`eqn'q`i'' = cond(`cond2`i''>=`k', -`d2_2q`i''*2*`c`k'2', 0);
	   		generate double `d3_`eqn'q`i'' = cond(`cond2`i''>=`k', -`d2_3q`i''*2*`c`k'2', 0);
  if $END {;
  			tempvar 		 d4_`eqn'q`i';
			generate double `d4_`eqn'q`i'' = cond(`cond2`i''>=`k', -`d2_4q`i''*2*`c`k'2', 0);
  };
	   		generate double `d`K1'_`eqn'q`i'' = cond(`cond2`i''>=`k', -`d2_`K1'q`i''*2*`c`k'2', 0);
		}; // k
	}; // i

	forvalues k = 1/`K1' {;
		forvalues l = `k'/`K1' {;
			tempvar g`k'_`l'i;
			generate double `g`k'_`l'i'= `d`k'_`l'q1'-`d`k'_`l'q2'-`d`k'_`l'q3'+`d`k'_`l'q4'+`g`k'i'*`g`l'i'/`lnf1';
		}; // l
	}; //k

			local g1_`K2'i    = "-`g1_2i'";
			local g2_`K2'i    = "-`g2_2i'";
			local g3_`K2'i    = "-`g2_3i'";
	if $END local g4_`K2'i    = "-`g2_4i'";
			local g`K1'_`K2'i = "-`g2_`K1'i'";
			local g`K2'_`K2'i = " `g2_2i'";

	forvalues k = 2/`NC1_1' {;
		local eqn = `K1'-1+`k';
		forvalues l = 1/`K1' {;
			tempvar g`l'_`eqn'i;
			generate double `g`l'_`eqn'i' = `d`l'_`eqn'q1'-`d`l'_`eqn'q2'-`d`l'_`eqn'q3'+`d`l'_`eqn'q4'+`g`l'i'*`g`eqn'i'/`lnf1';
		}; // l
		tempvar g`eqn'_`eqn'i;
		generate double `g`eqn'_`eqn'i' =(`d`K1'_`eqn'q1'-`d`K1'_`eqn'q2'-`d`K1'_`eqn'q3'+`d`K1'_`eqn'q4')*2*`c`k'1'-`g`eqn'i'/(`c`k'1')+`g`eqn'i'^2/`lnf1';
				  local  g`eqn'_`K2'i    = "-`g2_`eqn'i'";
		local k_1 = `k'+1;
		forvalues kk = `k_1'/`NC1_1' {;
			local eqn1 = `K1'-1+`kk';
			tempvar g`eqn'_`eqn1'i;
			generate double `g`eqn'_`eqn1'i' = (`d`K1'_`eqn1'q1'-`d`K1'_`eqn1'q2'-`d`K1'_`eqn1'q3'+`d`K1'_`eqn1'q4')*2*`c`k'1'+`g`eqn'i'*`g`eqn1'i'/`lnf1';
		}; // kk
		forvalues kk = 2/`NC2_1' {;
			local eqn2 = `K1'-1+`NC1_1'+`kk';
			tempvar g`eqn'_`eqn2'i;
			generate double `g`eqn'_`eqn2'i' = `d`eqn'_`eqn2'q1'-`d`eqn'_`eqn2'q2'-`d`eqn'_`eqn2'q3'+`d`eqn'_`eqn2'q4'+ `g`eqn'i'*`g`eqn2'i'/`lnf1';
		}; // kk
	}; // k
	forvalues k = 2/`NC2_1' {;
		local eqn = `K1'-1+`NC1_1'+`k';
		forvalues l = 1/`K1' {;
			tempvar g`l'_`eqn'i;
			generate double `g`l'_`eqn'i' =`d`l'_`eqn'q1'-`d`l'_`eqn'q2'-`d`l'_`eqn'q3'+`d`l'_`eqn'q4'+`g`l'i'*`g`eqn'i'/`lnf1';
		}; // l

		tempvar g`eqn'_`eqn'i;
		generate double `g`eqn'_`eqn'i' = -(`d2_`eqn'q1'-`d2_`eqn'q2'-`d2_`eqn'q3'+`d2_`eqn'q4')*2*`c`k'2'-`g`eqn'i'/(`c`k'2')+`g`eqn'i'^2/`lnf1';
				  local  g`K2'_`eqn'i    = "-`g2_`eqn'i'";
		local k_1 = `k'+1;
		forvalues kk = `k_1'/`NC2_1' {;
			local eqn2 = `K1'-1+`NC1_1'+`kk';
			tempvar g`eqn'_`eqn2'i;
			generate double `g`eqn'_`eqn2'i' =-(`d2_`eqn2'q1'-`d2_`eqn2'q2'-`d2_`eqn2'q3'+`d2_`eqn2'q4')*2*`c`k'2'+`g`eqn'i'*`g`eqn2'i'/`lnf1';
		}; // kk
	}; // k

	replace `g3_3i' = `g3_3i'-(`g3i'/`dtanh_dr')*`d2tanh_dr2';

	capture matrix drop `negH';
	forvalues k = 1/`N' {;
		tempname H;
		tempname g`k'_`k';
		mlmatsum `lnf' `g`k'_`k''  =  `g`k'_`k'i'/`lnf1', eq(`k');
		forvalues l = 1/`N' {;
			if (`l'>`k') {;
				tempname g`k'_`l';
				mlmatsum `lnf' `g`k'_`l''  =  `g`k'_`l'i'/`lnf1', eq(`k',`l');
				matrix `H' =nullmat(`H'),`g`k'_`l'';
			};
			else {;
		        matrix `H' =nullmat(`H'),`g`l'_`k''';
			};
		}; // l
		matrix `negH' = nullmat(`negH') \ `H';
	}; // k
	}; // end queitly
end; // end program bioprobite_d2

