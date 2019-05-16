*! version 2.3  15May2007 M. Lokshin, The World Bank
#delimit ;

program define gidecomposition, rclass;
	version 8.0;
	syntax using [fweight aweight] [if],
		   var1(varname numeric) var2(string) pline1(string) pline2(string) [hc pg pgs];

	if missing("`hc'`pg'`pgs'") | (`: word count `hc' `pg' `pgs' '>1) {;
		local hc "hc";
		local pg;
		local pgs;
	};

	preserve;
	quietly {;
    	tempvar p0 p1 p2;
		tempname M11 M12 M21 M22;

		if (!missing("`if'")) keep `if';

		summarize `var1' [`weight' `exp'], meanonly; local mean_1_1 = r(mean);

		generate  `p0' = (`var1' < `pline1') * 100;  								   // poverty headcount
    	generate  `p1' = (`var1' < `pline1') * ((`pline1' - `var1') / `pline1')* 100;  // poverty gap
    	generate  `p2' = `p1'^2/100;                                                   // poverty severety

		tabstat `p0' `p1' `p2' [`weight' `exp'], save;
		matrix `M11' = r(StatTot)';                          // in version 8 this was the matrix name


		/****************** FILE 2 *****************************************/
		use `using', clear; if ("`if'" ~= "") keep `if';

		summarize `var2' [`weight' `exp'], meanonly;  	local mean_2_2 = r(mean);

		generate  `p0' = (`var2' < `pline2') * 100;									    // poverty headcount 2
    	generate  `p1' = (`var2' < `pline2') * ((`pline2' - `var2') / `pline2') * 100;  // poverty gap 2
    	generate  `p2' = `p1'^ 2/100;                     							    // poverty severity

		tabstat `p0' `p1' `p2' [`weight' `exp'], save;
		matrix `M22' = r(StatTot)';

		replace   `var2' = `var2' * `mean_1_1' / `mean_2_2';  					  // change means in file 2

		replace   `p0' = (`var2' < `pline2') * 100;
    	replace   `p1' = (`var2' < `pline2') * ((`pline2' - `var2') / `pline2') * 100;  // poverty gap 2
    	replace   `p2' = `p1'^ 2/100;                     							  // poverty severity

		tabstat `p0' `p1' `p2' [`weight' `exp'], save;
		matrix `M12' = r(StatTot)';

		restore, preserve;
		/***************** FILE 1 ******************************************/
		if ("`if'" ~= "") keep `if';
		replace   `var1' = `var1' * `mean_2_2' / `mean_1_1';

		generate  `p0' = (`var1' < `pline1') * 100;									  // poverty headcount 2
    	generate  `p1' = (`var1' < `pline1') * ((`pline1' - `var1') / `pline1') * 100;  // poverty gap 2
    	generate  `p2' = `p1'^ 2 / 100;                     							  // poverty severity

		tabstat `p0' `p1' `p2' [`weight' `exp'], save;
		matrix `M21' = r(StatTot)';

		tempname b P CHN GRO RED RES;

		matrix `CHN' = (`M22' - `M11');
		matrix `GRO' = (`M21' - `M11'), (`M22' - `M12');
		matrix `RED' = (`M12' - `M11'), (`M22' - `M21');
		matrix `RES' = (`CHN', `CHN') - `GRO' - `RED';

		matrix `P' = `M11', `M22', `CHN';
		matrix `b' = (`GRO'[1,....], `GRO'[2,....], `GRO'[3,....]) \
				     (`RED'[1,....], `RED'[2,....], `RED'[3,....]) \
				     (`RES'[1,....], `RES'[2,....], `RES'[3,....]) ;

		matrix rownames `b'= growth redistribution interaction;
		matrix colnames `b'=P0_base1 P0_base2 P1_base1 P1_base2 P2_base1 P2_base2;
		matrix rownames `P'=P0 P1 P2;
		matrix colnames `P'=base1 base2 change;
	}; // end quietly

	if (!missing("`hc'"))  {; local index = 1; local title1 "   Poverty rate (P0)";     local title2 "P0"; };
	if (!missing("`pg'"))  {; local index = 2; local title1 "   Poverty gap (P1)";      local title2 "P1"; };
	if (!missing("`pgs'")) {; local index = 3; local title1 "   Poverty severity (P2)"; local title2 "P2"; };

	display as text _newline "{hline 80}";
	display as text "   Growth and Inequality Poverty Decomposition";
	display as text "{hline 80}";
	display as text _col(30) "Base year 1" _col(48) "Base year 2" _col(64) "Average effect";
	display as text "{hline 80}";
	display as text "`title1'"  _col(32) as result  %6.3f (`=`P'[`index', 1]') _col(50) %6.3f (`=`P'[`index', 2]');
	display as text "{hline 80}";
	display as text "   Change in `title2'"
		_col(32) as result %6.3f (`=`P'[`index', 3]')
		_col(50) 		   %6.3f (`=`P'[`index', 3]')
		_col(68) 		   %6.3f (`=`P'[`index', 3]');
	display as text "{hline 80}";
	display as text "   Growth component"
		_col(32) as result %6.3f (`=`GRO'[`index', 1]')
		_col(50) 		   %6.3f (`=`GRO'[`index', 2]')
		_col(68) 		   %6.3f (`=`GRO'[`index', 1]+`GRO'[`index', 2]')/2;
	display as text "   Redistribution component"
		_col(32) as result %6.3f (`=`RED'[`index', 1]')
		_col(50) 		   %6.3f (`=`RED'[`index', 2]')
		_col(68) 		   %6.3f (`=`RED'[`index', 1]+`RED'[`index', 2]')/2;
	display as text "   Interaction component"
	   	_col(32) as result %6.3f (`=`RES'[`index', 1]')
		_col(50) 		   %6.3f (`=`RES'[`index', 1]')
		_col(68) 		   %6.3f (`=`RES'[`index', 1]+`RES'[`index', 2]')/2;
	display as text "{hline 80}";

	return matrix b = `b';
	return matrix P = `P';
	restore;
end;

