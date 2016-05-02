# delimit ;
program define _gfgt;
	syntax newvarname =/exp [if] [in] [, BY(varlist)
										 Weights(varname)
										 Alpha(numlist max=1 >=0)]
										 PLine(string);

	tempname touse poor z nk vk fik;

	capture confirm number `pline';
	if _rc {;
		confirm numeric variable `pline';
		local z `pline';
	};
	else generate `z'=`pline';

	local inc : subinstr local exp "(" "";
	local inc : subinstr local inc ")" "";
	local inc = trim("`inc'");
	confirm variable `inc';

	marksample touse, novarlist;
	markout `touse' `inc' `by';

	if ("`by'" == "") {;
		tempvar by;
		generate byte `by' = 1;
	};
	else local bystr "by(`by')";

	if ("`weights'" == "")  {;
		tempvar wi;
		generate byte `wi' = 1;
	};
	else {;
		tempname wi;
		local wi `weights';
	};
	if ("`alpha'"=="") {;
		local alpha = 0;
	};
	quietly {;{;
		sum `inc' [w=`wi'] if `touse';
		local sumwi = _result(2);

		generate `poor' = (`inc' < `z') if `touse';
		egen `nk'       = sum(`wi') if `touse', by(`by');
		generate `vk'   = `nk'/`sumwi' if `touse';
		generate `fik'  = `wi'/`nk' if `touse';
		egen `varlist'  = sum(`fik'*`poor'*((`z'-`inc')/`z')^(`alpha')) if `touse', by(`by');
		local index = cond(`alpha'==0, "(headcount)",
					  cond(`alpha'==1, "(poverty gap)",
					  cond(`alpha'==2, "(poverty severity)","(fgt(`alpha'))")));
		label variable `varlist' "`index' `inc' `bystr'";
   }; // quietly
end; // program _ggft
