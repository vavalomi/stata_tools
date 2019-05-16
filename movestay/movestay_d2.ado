*! version 3.0.1  10Mar2008 M. Lokshin, Z. Sajaia
#delim ;
program define movestay_d2, eclass;
	args todo b lnf g negH g1 g2 g3 g4 g5 g6 g7;

	tempname lns1 lns0 theta0 theta1 xb0 xb1 xb3;

	mleval `xb0' 	=`b', eq(1);
	mleval `xb1' 	=`b', eq(2);
	mleval `xb3' 	=`b', eq(3);
	mleval `lns0'	=`b', eq(4) scalar;
	mleval `lns1'	=`b', eq(5) scalar;
	mleval `theta0'	=`b', eq(6) scalar;
	mleval `theta1'	=`b', eq(7) scalar;

	quietly {;

	tempname sig1 sig0 rho1 rho0 rr1 rr0;

	scalar `sig0' 	= exp(`lns0');
	scalar `sig1' 	= exp(`lns1');
	scalar `rho0' 	= tanh(`theta0');
	scalar `rho1' 	= tanh(`theta1');
	scalar `rr0'	= 1/sqrt(1-`rho0'^2);
	scalar `rr1'	= 1/sqrt(1-`rho1'^2);

	tempvar  eta0 eta1 lf1 eps0 eps1;

	generate double `eta0' = (`xb3' +  ($ML_y1 - `xb0') * `rho0'/`sig0')*`rr0'; replace `eta0'= 37 if (`eta0'> 37);
	generate double `eta1' = (`xb3' +  ($ML_y2 - `xb1') * `rho1'/`sig1')*`rr1'; replace `eta1'=-37 if (`eta1'<-37);
	generate double `eps0' = $ML_y1 - `xb0';
	generate double `eps1' = $ML_y2 - `xb1';

	tempname const0 const1;

	scalar `const0'=0.5*ln(2*_pi*`sig0'^2);
	scalar `const1'=0.5*ln(2*_pi*`sig1'^2);

	mlsum `lnf' = cond($ML_y3==1, ln(norm(  `eta1'))-`const1'-0.5*(`eps1'/`sig1')^2,
				                  ln(norm(- `eta0'))-`const0'-0.5*(`eps0'/`sig0')^2);

	if (`todo'==0 | `lnf'==.)  exit;  // calculate first derivative

	tempname  deta_dx deta_dsig deta_drho M;

	generate double `deta_dx'   = -cond($ML_y3==1,`rho1'/`sig1'*`rr1',`rho0'/`sig0'*`rr0');
	generate double `deta_dsig' = -cond($ML_y3==1,`eps1'*`rho1'/`sig1'^2*`rr1',`eps0'*`rho0'/`sig0'^2*`rr0');

	generate double `deta_drho' =  cond($ML_y3==1,`eps1'*`rr1'/`sig1'+(`rho1'*(`xb3'+`rho1'*`eps1'/`sig1'))*`rr1'^3,
	                                              `eps0'*`rr0'/`sig0'+(`rho0'*(`xb3'+`rho0'*`eps0'/`sig0'))*`rr0'^3);

	generate double `M'         =  cond($ML_y3==1,normden(`eta1')/norm(`eta1'),-normden(`eta0')/norm(-`eta0'));

	tempname drho0_dtheta0 drho1_dtheta1;

	scalar `drho0_dtheta0' 	= 4*exp(2*`theta0')/((1+exp(2*`theta0'))^2);	  // derivative of rho w.r.t theta
  	scalar `drho1_dtheta1' 	= 4*exp(2*`theta1')/((1+exp(2*`theta1'))^2);	  // derivative of rho w.r.t theta

	replace `g1' 	 =  cond($ML_y3==0,`M'*`deta_dx'+`eps0'/(`sig0'^2),0);
	replace `g2' 	 =  cond($ML_y3==1,`M'*`deta_dx'+`eps1'/(`sig1'^2),0);
	replace `g3'	 =  cond($ML_y3==1,`M'*`rr1',`M'*`rr0');
	replace `g4' 	 =  cond($ML_y3==0,(`M'*`deta_dsig'-(1/`sig0')+(`eps0'^2)/(`sig0'^3))*`sig0',0);
	replace `g5' 	 =  cond($ML_y3==1,(`M'*`deta_dsig'-(1/`sig1')+(`eps1'^2)/(`sig1'^3))*`sig1',0);
	replace `g6'	 =  cond($ML_y3==0,`M'*`deta_drho'*`drho0_dtheta0',0);
	replace `g7'	 =  cond($ML_y3==1,`M'*`deta_drho'*`drho1_dtheta1',0);

	tempname d1 d2 d3 d4 d5 d6 d7;

	mlvecsum `lnf' `d1' = `g1' , eq(1);
	mlvecsum `lnf' `d2' = `g2' , eq(2);
	mlvecsum `lnf' `d3' = `g3' , eq(3);
	mlvecsum `lnf' `d4' = `g4' , eq(4);
	mlvecsum `lnf' `d5' = `g5' , eq(5);
	mlvecsum `lnf' `d6' = `g6' , eq(6);
	mlvecsum `lnf' `d7' = `g7' , eq(7);

	matrix `g' = (`d1' , `d2' , `d3' , `d4' , `d5' , `d6', `d7' );

	if (`todo'==1 | `lnf'==. )  exit;   // calculate second derivative

   	tempvar  DM d2eta_drho2 d2rho0_dtheta02 d2rho1_dtheta12;

	generate double `DM'          = -cond($ML_y3==1,`M'*(`eta1'+`M'),`M'*(`eta0'+`M'));

	generate double `d2eta_drho2' =  cond($ML_y3==1,`rr1'^3*(2*`rho1'*`eps1'/`sig1'+(`xb3'+`rho1'*`eps1'/`sig1')*(1+3*(`rho1'*`rr1')^2)),
                                                    `rr0'^3*(2*`rho0'*`eps0'/`sig0'+(`xb3'+`rho0'*`eps0'/`sig0')*(1+3*(`rho0'*`rr0')^2)));

	scalar `d2rho1_dtheta12'      = 8*exp(2*`theta1')*(1-exp(2*`theta1'))/(1+exp(2*`theta1'))^3;	  // 2nd derivative of rho w.r.t theta
	scalar `d2rho0_dtheta02'      = 8*exp(2*`theta0')*(1-exp(2*`theta0'))/(1+exp(2*`theta0'))^3;

	tempvar	 	g11 g12 g13 g14 g15 g16 g17
		    	    g22 g23 g24 g25 g26 g27
	        	        g33 g34 g35 g36 g37
	    		            g44 g45 g46 g47
			                	g55 g56 g57
    			                    g66 g67
					   					g77
				;

	tempvar	 	d11 d12 d13 d14 d15 d16 d17
         		    d22 d23 d24 d25 d26 d27
						d33 d34 d35 d36 d37
		                    d44 d45 d46 d47
								d55 d56 d57
			                    	d66 d67
										d77
				;

	generate double `g22' = - cond($ML_y3==1,`DM'*`deta_dx'^2-1/`sig1'^2,0);
	generate double `g23' = - cond($ML_y3==1,`DM'*`deta_dx'*`rr1',0);
	generate double `g25' = - cond($ML_y3==1,(`DM'*`deta_dsig'-`M'/`sig1')*`deta_dx'-2*`eps1'/`sig1'^3,0);
	generate double `g27' = - cond($ML_y3==1,(`DM'*`deta_drho'+`M'*(1/`rho1'+`rho1'*`rr1'^2))*`deta_dx',0);

	generate double `g33' = - cond($ML_y3==1,`DM'*`rr1'^2,`DM'*`rr0'^2);
	generate double `g35' = - cond($ML_y3==1,`DM'*`deta_dsig'*`rr1',0);
	generate double `g34' = - cond($ML_y3==0,`DM'*`deta_dsig'*`rr0',0);
	generate double `g36' = - cond($ML_y3==0,`DM'*`deta_drho'*`rr0'+`M'*`rho0'*`rr0'^3,0);
	generate double `g37' = - cond($ML_y3==1,`DM'*`deta_drho'*`rr1'+`M'*`rho1'*`rr1'^3,0);

	generate double `g55' = - cond($ML_y3==1,(`DM'*`deta_dsig'^2-2*`M'*`deta_dsig'/`sig1'+1/`sig1'^2-3*`eps1'^2/`sig1'^4)*`sig1',0)-`g5';
	generate double `g57' = - cond($ML_y3==1,`DM'*`deta_drho'*`deta_dsig'+`M'*`deta_dsig'*(1/`rho1'+`rho1'*`rr1'^2),0);

	generate double `g44' = - cond($ML_y3==0,(`DM'*`deta_dsig'^2-2*`M'*`deta_dsig'/`sig0'+1/`sig0'^2-3*`eps0'^2/`sig0'^4)*`sig0',0)-`g4';
	generate double `g46' = - cond($ML_y3==0,`DM'*`deta_drho'*`deta_dsig'+`M'*`deta_dsig'*(1/`rho0'+`rho0'*`rr0'^2),0);

	generate double `g66' = - cond($ML_y3==0,(`DM'*`deta_drho'^2+`M'*`d2eta_drho2')*`drho0_dtheta0'^2,0)-`g6'*`d2rho0_dtheta02';

	generate double `g77' = - cond($ML_y3==1,(`DM'*`deta_drho'^2+`M'*`d2eta_drho2')*`drho1_dtheta1'^2,0)-`g7'*`d2rho1_dtheta12';

	mlmatsum `lnf' `d11'=-(`DM'*`deta_dx'^2-1/`sig0'^2)    											 if ($ML_y3==0), eq(1);
	mlmatsum `lnf' `d12'= 0, 				  		   												 			 	 eq(1,2);
	mlmatsum `lnf' `d13'=-`DM'*`deta_dx'*`rr0' 												 		 if ($ML_y3==0), eq(1,3);
	mlmatsum `lnf' `d14'=-((`DM'*`deta_dsig'-`M'/`sig0')*`deta_dx'-2*`eps0'/`sig0'^3)*`sig0' 		 if ($ML_y3==0), eq(1,4);
	mlmatsum `lnf' `d15'= 0, 				  																		 eq(1,5);
	mlmatsum `lnf' `d16'=-(`DM'*`deta_drho'+`M'*(1/`rho0'+`rho0'*`rr0'^2))*`deta_dx'*`drho0_dtheta0' if ($ML_y3==0), eq(1,6);
	mlmatsum `lnf' `d17'= 0, 																			  			 eq(1,7);

	mlmatsum `lnf' `d22'=`g22', 						eq(2);
	mlmatsum `lnf' `d23'=`g23', 				 		eq(2,3);
	mlmatsum `lnf' `d24'= 0, 				  			eq(2,4);
	mlmatsum `lnf' `d25'=`g25'*`sig1',       			eq(2,5);
	mlmatsum `lnf' `d26'= 0, 				  			eq(2,6);
	mlmatsum `lnf' `d27'=`g27'*`drho1_dtheta1', 		eq(2,7);

	mlmatsum `lnf' `d33'=`g33', 				  		eq(3);
	mlmatsum `lnf' `d34'=`g34'*`sig0',	        		eq(3,4);
	mlmatsum `lnf' `d35'=`g35'*`sig1',       			eq(3,5);
	mlmatsum `lnf' `d36'=`g36'*`drho0_dtheta0',			eq(3,6);
	mlmatsum `lnf' `d37'=`g37'*`drho1_dtheta1', 		eq(3,7);

	mlmatsum `lnf' `d44'=`g44'*`sig0',                  eq(4);
	mlmatsum `lnf' `d45'= 0, 				  			eq(4,5);
	mlmatsum `lnf' `d46'=`g46'*`sig0'*`drho0_dtheta0',	eq(4,6);
	mlmatsum `lnf' `d47'=0,  							eq(4,7);

	mlmatsum `lnf' `d55'=`g55'*`sig1',              	eq(5);
	mlmatsum `lnf' `d56'=0,								eq(5,6);
	mlmatsum `lnf' `d57'=`g57'*`sig1'*`drho1_dtheta1', 	eq(5,7);

	mlmatsum `lnf' `d66'=`g66',	                	  	eq(6);
	mlmatsum `lnf' `d67'= 0, 		                  	eq(6,7);

	mlmatsum `lnf' `d77'=`g77',              		  	eq(7);

	capture matrix drop `negH';
	forvalues k = 1/7 {;
		tempname H;
		forvalues l = 1/7 {;
			if (`l'>`k') matrix `H' =nullmat(`H'),`d`k'`l'';
			else         matrix `H' =nullmat(`H'),`d`l'`k''';
		}; // l
		matrix `negH' = nullmat(`negH') \ `H';
	}; // k

	}; // end quietly
end;

