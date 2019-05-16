#delimit ;

program define switch_probit_d2;
	args todo b lnf g negH g1 g2 g3 g4 g5;

	tempname theta1 theta0 zb xb1 xb0;
	tempvar lnf1;

	mleval `zb' 	=`b', eq(1);
	mleval `xb1' 	=`b', eq(2);
	mleval `xb0' 	=`b', eq(3);
	mleval `theta1'	=`b', eq(4) scalar;
	mleval `theta0'	=`b', eq(5) scalar;

	tempname rho1 rho0;

	scalar `rho1' 	= tanh(`theta1');
    scalar `rho0'   = tanh(`theta0');

quietly {;

	tempvar qseli qi xb rhoi;
	generate byte `qseli' = (2*$ML_y1-1);
  	generate byte `qi'    = (2*$ML_y2-1) if $ML_y1==1;
    replace       `qi'    = (2*$ML_y3-1) if $ML_y1==0;

	generate double `xb'  = `xb1' if $ML_y1==1;
    replace         `xb'  = `xb0' if $ML_y1==0;

    generate double `rhoi' = `qseli'*`qi'*`rho1' if $ML_y1==1;
    replace         `rhoi' = `qseli'*`qi'*`rho0' if $ML_y1==0;

    generate double `lnf1' = binorm(`qseli'*`zb', `qi'*`xb', `rhoi');

	mlsum `lnf' = log(`lnf1');

  	if (`todo'==0 | `lnf'==.)  exit;

	tempvar di w1i w2i g1i g2i phi2;

	generate double `di'   = 1/(sqrt(1-`rhoi'^2));
	generate double `w1i'  = `qseli'*`zb';
	generate double `w2i'  = `qi'   *`xb';
	generate double `g1i'  = normden(`w1i')*norm(`di' * (`w2i' - `rhoi'*`w1i'));
	generate double `g2i'  = normden(`w2i')*norm(`di' * (`w1i' - `rhoi'*`w2i'));
	generate double `phi2' = exp(-0.5*(`w1i'^2+`w2i'^2-2*`rhoi'*`w1i'*`w2i')/(1-`rhoi'^2))/(2*_pi*sqrt(1-`rhoi'^2));

	tempname drho1_dtheta1 drho0_dtheta0;
	scalar `drho1_dtheta1' = 4*exp(2*`theta1')/((1+exp(2*`theta1'))^2);
	scalar `drho0_dtheta0' = 4*exp(2*`theta0')/((1+exp(2*`theta0'))^2);

	// scores
	replace `g1' = `qseli'*`g1i'/`lnf1';
    replace `g2' =  `qi'  *`g2i'/`lnf1'                *($ML_y1==1);
    replace `g3' =  `qi'  *`g2i'/`lnf1'                *($ML_y1==0);
	replace `g4' =  `qi' *`phi2'/`lnf1'*`drho1_dtheta1'*($ML_y1==1);
	replace `g5' = -`qi' *`phi2'/`lnf1'*`drho0_dtheta0'*($ML_y1==0);

	tempname d1 d2 d3 d4 d5;

	mlvecsum `lnf' `d1'  = `g1', eq(1);
    mlvecsum `lnf' `d2'  = `g2', eq(2);
    mlvecsum `lnf' `d3'  = `g3', eq(3);
	mlvecsum `lnf' `d4'  = `g4', eq(4);
	mlvecsum `lnf' `d5'  = `g5', eq(5);

	matrix `g' = (`d1', `d2', `d3', `d4', `d5');

	if (`todo'==1 | `lnf'==. )  exit;

	tempname d2rho1_dtheta12 d2rho2_dtheta22;

	scalar `d2rho1_dtheta12' = 8*exp(2*`theta1')*(1-exp(2*`theta1'))/(exp(2*`theta1')+1)^3;
	scalar `d2rho2_dtheta22' = 8*exp(2*`theta0')*(1-exp(2*`theta0'))/(exp(2*`theta0')+1)^3;

	tempname d1d1 d1d2 d1d3 d1d4 d1d5
				  d2d2 d2d3 d2d4 d2d5
					   d3d3 d3d4 d3d5
							d4d4 d4d5
								 d5d5;


	mlmatsum `lnf' `d1d1' = (`w1i'*`g1i'+`rhoi'*`phi2'+`g1i'^2/`lnf1')/`lnf1', eq(1);
	mlmatsum `lnf' `d1d2' = (`qseli'*`qi')/`lnf1'*(`g1i'*`g2i'/`lnf1'-`phi2') if $ML_y1==1, eq(1, 2);
	mlmatsum `lnf' `d1d3' = (`qseli'*`qi')/`lnf1'*(`g1i'*`g2i'/`lnf1'-`phi2') if $ML_y1==0, eq(1, 3);
	mlmatsum `lnf' `d1d4' = `qi'*`phi2'/`lnf1'*`drho1_dtheta1'*
							(-`rhoi'*`di'^2*(`w2i'-`rhoi'*`w1i')+`w1i'+`g1i'/`lnf1') if $ML_y1==1, eq(1, 4);

	mlmatsum `lnf' `d1d5' = `qi'*`phi2'/`lnf1'*`drho0_dtheta0'*
							(-`rhoi'*`di'^2*(`w2i'-`rhoi'*`w1i')+`w1i'+`g1i'/`lnf1') if $ML_y1==0, eq(1, 5);

	mlmatsum `lnf' `d2d2' = (`w2i'*`g2i'+`rhoi'*`phi2'+`g2i'^2/`lnf1')/`lnf1' if $ML_y1==1, eq(2);
	mlmatsum `lnf' `d2d3' = 0, eq(2, 3);
	mlmatsum `lnf' `d2d4' = `qseli'*`phi2'/`lnf1'*`drho1_dtheta1'*
							(-`rhoi'*`di'^2*(`w1i'-`rhoi'*`w2i')+`w2i'+`g2i'/`lnf1') if $ML_y1==1, eq(2, 4);
	mlmatsum `lnf' `d2d5' = 0, eq(2, 5);

	mlmatsum `lnf' `d3d3' = (`w2i'*`g2i'+`rhoi'*`phi2'+`g2i'^2/`lnf1')/`lnf1' if $ML_y1==0, eq(3);
	mlmatsum `lnf' `d3d4' = 0, eq(3, 4);
	mlmatsum `lnf' `d3d5' = `qseli'*`phi2'/`lnf1'*`drho0_dtheta0'*
							(-`rhoi'*`di'^2*(`w1i'-`rhoi'*`w2i')+`w2i'+`g2i'/`lnf1') if $ML_y1==0, eq(3, 5);

	mlmatsum `lnf' `d4d4' = (`phi2'/`lnf1')*(`di'^2*`rhoi'*(-1+`di'^2*(`w1i'^2+`w2i'^2-2*`rhoi'*`w1i'*`w2i'))-
							`di'^2*`w1i'*`w2i'+`phi2'/`lnf1')*((`drho1_dtheta1')^2)-`qseli'*`qi'*`phi2'/`lnf1'*`d2rho1_dtheta12' if $ML_y1==1, eq(4);
	mlmatsum `lnf' `d4d5' = 0, eq(4, 5);

	mlmatsum `lnf' `d5d5' = (`phi2'/`lnf1')*(`di'^2*`rhoi'*(-1+`di'^2*(`w1i'^2+`w2i'^2-2*`rhoi'*`w1i'*`w2i'))-
							`di'^2*`w1i'*`w2i'+`phi2'/`lnf1')*((`drho0_dtheta0')^2)-`qseli'*`qi'*`phi2'/`lnf1'*`d2rho2_dtheta22' if $ML_y1==0, eq(5);

	matrix `negH' = (`d1d1',  `d1d2',  `d1d3',  `d1d4',  `d1d5' \
				     `d1d2'', `d2d2',  `d2d3',  `d2d4',  `d2d5' \
				     `d1d3'', `d2d3'', `d3d3',  `d3d4',  `d3d5' \
				     `d1d4'', `d2d4'', `d3d4'', `d4d4',  `d4d5' \
				     `d1d5'', `d2d5'', `d3d5'', `d4d5'', `d5d5' );
}; // quietly
end;

