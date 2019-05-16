*! version 2.31 26Mar2008 M. Lokshin

#delimit ;

program define sedecomposition, rclass;
version 8.0;

syntax using [if] [aweight fweight], sector(varname numeric)
	[pline1(string) pline2(string) var1(varname numeric) var2(string) hc pg pgs] ;

/* possible errors */
if (missing("`poor1'") & (missing("`var1'") | missing("`pline1'")) ) {;
	display as error "You need to specify either pline()--var() pairs, or poor()";
	exit;
};
if ("`poor1'"=="" & "`var1'"=="") {;
	display as error "You need to specify welf() if you do not specify poor()";
	exit;
};
if ("`poor1'"~="" & "`poor2'"=="") {;
	display as error "You need to specify poor2 if you specify poor1";
	exit;
};
if ("`var1'"~="" & "`var2'"=="") {;
	display as error "You need to specify var2 if you specify var1";
	exit;
};

quietly {;

preserve;

tempname POPS1 POOR_S1;
tempvar _fgt1 _fgt2 _poor;

	local fgt=1; local titl1 = "HeadCount";		  	  local titl2="(HC)";	  /* indicator for HC 				*/
if ("`pg'" ~="") {;
	local fgt=2; local titl1 = "Poverty Gap"; 		  local titl2="(PG)"; };  /* indicator for PG		 		*/
if ("`pgs'"~="") {;
	local fgt=3; local titl1 = "Poverty Gap Squared"; local titl2="(PG2)"; }; /* indicator for PG2				*/

if ("`if'"~="") keep `if'; 													/* restrict observation for user-specified if conditions */

generate byte `_poor'=(`var1'<=`pline1');								/* HC        						*/
generate 	  `_fgt1'=`_poor'*((`pline1'-`var1')/`pline1');  		    /* PG        						*/
generate 	  `_fgt2'=`_fgt1'^2;         								/* SP      							*/

tabulate `sector' [`weight' `exp'], matcell(`POPS1'); local n_sec1=r(r); local n_obs1=r(N);
matrix `POPS1' = `POPS1'/`n_obs1';

tabstat `_poor' `_fgt1' `_fgt2' [`weight' `exp'], by(`sector') s(mean) save;
local p1t_1 = el(r(StatTot),1,1); 											/* total poverty level 				*/
local p2t_1 = el(r(StatTot),1,2);                                           /* total poverty gap				*/
local p3t_1 = el(r(StatTot),1,3);                                           /* total poverty gap squared 		*/
local rlbl	= subinstr(r(name1)," ","_",10);
local lbl1 	= r(name1);
matrix `POOR_S1' = r(Stat1);   			                  					/* means and count					*/
forvalues i=2/`n_sec1' {;
	matrix `POOR_S1' = `POOR_S1' \ r(Stat`i');
	local lbl`i' = r(name`i');												/* main labels 						*/
	local rlbl `"`rlbl' `"`lbl`i''"'"';
};
matrix colname `POOR_S1'=P0 P1 P2;
matrix rowname `POOR_S1'=`rlbl';

/* second file */

use `using', clear;
tempname POOR_S2 POPS2 INTER_S INTRA_S SHIFT_S AC_1 AC_2 AC_3 AC PC_1 PC_2 PC_3 PC;

if ("`if'"~="") keep `if';

generate byte 	`_poor' = (`var2'<=`pline2');						 	/* HC       						*/
generate 		`_fgt1'	= `_poor'*((`pline2'-`var2')/`pline2');        	/* PG        						*/
generate 		`_fgt2'	= `_fgt1'^2;         					   		/* SP        						*/

tabulate `sector' [`weight' `exp'], matcell(`POPS2'); local n_sec2=r(r); local n_obs2=r(N);
matrix `POPS2' = `POPS2'/`n_obs2';
if (`n_sec2'~=`n_sec1') {;
	noisily disp _newline in red "Number of sectors is different in year 1 and year 2";
	exit; };

tabstat `_poor' `_fgt1' `_fgt2' [`weight' `exp'], by(`sector') s(mean) save;
local p1t_2 = el(r(StatTot),1,1); 											/* total poverty level 		 	*/
local p2t_2 = el(r(StatTot),1,2);                                           /* total poverty gap				*/
local p3t_2 = el(r(StatTot),1,3);                                           /* total poverty gap squared 		*/
matrix 			`POOR_S2' = r(Stat1);
forvalues i=2/`n_sec2' {;	matrix `POOR_S2' = `POOR_S2' \ r(Stat`i');	};
matrix colname 	`POOR_S2' = P0 P1 P2;


matrix `AC_1'=vecdiag((`POOR_S2'[1...,"P0"]-`POOR_S1'[1...,"P0"])*`POPS1'')'; /* absulute change P0 	  		*/
matrix `AC_2'=vecdiag((`POOR_S2'[1...,"P1"]-`POOR_S1'[1...,"P1"])*`POPS1'')'; /* absulute change P0 	  		*/
matrix `AC_3'=vecdiag((`POOR_S2'[1...,"P2"]-`POOR_S1'[1...,"P2"])*`POPS1'')'; /* absulute change P0 	  		*/
matrix `AC'  =`AC_1',`AC_2',`AC_3';

local p1_change	= `p1t_2'-`p1t_1';
local p2_change	= `p2t_2'-`p2t_1';
local p3_change	= `p3t_2'-`p3t_1';

matrix `PC_1'=`AC_1'*100/`p1_change';   									/* percent change relative to total poverty change */
matrix `PC_2'=`AC_2'*100/`p2_change';                                       /* % change in poverty gap */
matrix `PC_3'=`AC_3'*100/`p3_change';                                       /* % change in poverty gap squared */
matrix `PC'  =`PC_1',`PC_2',`PC_3';

/* Total effect calculation */
matrix `INTRA_S'=(`POOR_S2'[1...,`fgt']-`POOR_S1'[1...,`fgt'])'*`POPS1';				local intra_e=`INTRA_S'[1,1];
matrix `INTER_S'=(`POOR_S2'[1...,`fgt']-`POOR_S1'[1...,`fgt'])'*(`POPS2'-`POPS1');		local inter_e=`INTER_S'[1,1];
matrix `SHIFT_S'=(`POPS2'			   -`POPS1'		  	  )'*`POOR_S1'[1...,`fgt'];	local shift_e=`SHIFT_S'[1,1];
}; // end quietly

	display as text _newline "{hline 70}";
	display as text "   Sectoral Decomposition of a Change in Poverty:  " "`titl1'";
	display as text "{hline 70}";
	display as text "   Poverty in period 1   " "`titl1'" _col(51) in yellow  %6.4f `p`fgt't_1'*100 ;
	display as text "   Poverty in period 2   " "`titl1'" _col(51) in yellow  %6.4f `p`fgt't_2'*100 ;
	display as text "{hline 70}";
	display as text "   Sector                    Population share    Absolute   Percentage";
	display as text "                               in period 1        change       change";
 	display as text "{hline 70}";

 	forvalues i=1(1)`n_sec1'{;
		disp in yellow 	_col(4) "`lbl`i''"
						_col(32) %8.2f `POPS1'[`i',1]*100
						_col(49) %8.4f `AC'[`i',`fgt']*100
						_col(65) %6.2f `PC'[`i',`fgt'];  };
	matrix 		   b_sec = (`POPS1'[1...,1]*100, `AC'[1...,`fgt']*100, `PC'[1...,`fgt']);
	matrix rowname b_sec = `rlbl';
	matrix colname b_sec = "Population share in period 1" "Absolute change" "Percentage change";

	disp in green "{hline 70}";
	disp in green "   Total Intra-sectoral effect"  _col(49)
		 in yellow %8.4f `intra_e'*100 _col(65) %6.2f `intra_e'*100/`p`fgt'_change';
	disp in green "   Population-shift effect"      _col(49)
		 in yellow %8.4f `shift_e'*100 _col(65) %6.2f `shift_e'*100/`p`fgt'_change';
	disp in green "   Interaction effect" 	        _col(49)
	     in yellow %8.4f `inter_e'*100 _col(65) %6.2f `inter_e'*100/`p`fgt'_change';
	disp in green "{hline 70}";
	disp in green "   Change in poverty " "`titl2'" _col(34)
		 in yellow _col(49) %8.4f `p`fgt'_change'*100 _col(65) %6.2f 100 in green _newline "{hline 70}";
	matrix b_tot = (`p`fgt'_change'*100, 100 \
	  				`intra_e'*100      ,`intra_e'*100/`p`fgt'_change' \
					`shift_e'*100      ,`shift_e'*100/`p`fgt'_change' \
					`inter_e'*100      , `inter_e'*100/`p`fgt'_change');
	matrix rowname b_tot ="Change in poverty `titl2'" "Total Intra-sectoral effect" "Population-shift effect" "Interaction effect";
	matrix colname b_tot ="Absolute change" "Percentage change";

	return matrix b_sec = b_sec;
	return matrix b_tot = b_tot;
restore;
end;




