*! version 2.22  17Mar2008 M. Lokshin, Z. Sajaia
# delimit ;

program define lookfor_all;
	version 9.0;

	syntax [anything] [, SUBdir DIRectory(string) Describe CODEbook *];


	if c(changed) {;
		tempfile tmpfile;
		quietly save `tmpfile';
	};
	else local tmpfile `"`c(filename)'"';

	drop _all; // important

	tempname D; matrix `D'=J(1,100,1);
	local i=1;
	local cur_dir `c(pwd)';

	if (!missing(`"`directory'"')) quietly cd `"`directory'"';

	_lookforall `anything', `describe' `options' `codebook';  // search in the root directory

	if (!missing(`"`subdir'"')) {;  // case when subdirectory option specified
		while (1) {;
			local dirs : dir . dirs "*";
			local j=`D'[1,`i'];
			local dirname : word `j' of `dirs';
			if ((!missing("`dirname'")) & _rc!=601) {;
				capture cd "`dirname'";                        // in case we found system folder
				if (_rc!=170) {;
					_lookforall `anything', `describe' `options' `codebook';
					local ++i;
				};
				else matrix `D'[1,`i']=`D'[1,`i']+1;
			};
			else {;
				matrix `D'[1,`i']=1;
				local i=`i'-1;
				quietly cd ..;
				if (`i'==0) continue, break;
				matrix `D'[1,`i']=`D'[1,`i']+1;
			};
		};
	quietly cd "`cur_dir'";
	};

	if ~missing(`"`tmpfile'"') quietly use `"`tmpfile'"';
end;

/*******************************************************************/
/*******************************************************************/

program define _lookforall;
	syntax [anything] [, Describe CODEbook *];

	local all_files : dir . files "*.dta";
	local n_files   : word count `all_files';
	local c_files = 0;
	foreach file of local all_files {;
		if missing("`codebook'") local instr "in 1";
		quietly capture use `"`file'"' `instr';
		if (_rc == 459) display as error  "File "
								as result "`file'"
								as error  " cannot be open in current version of Stata";
		else {;
			local c_files=`c_files'+(_rc==0);
		    local cmd "use "`c(pwd)'`c(dirsep)'`file'"";
			quietly lookfor `anything';
			if (!missing(r(varlist))) {;
				local vlist `r(varlist)';
				if(missing("`codebook'") & missing("`describe'"))
					display `"{ stata `"`cmd'"' : `c(pwd)'`c(dirsep)'`file' }"', as text "`vlist'";
				else {;
					display _newline `"{ stata `"`cmd'"' : `c(pwd)'`c(dirsep)'`file' }"';
					if(!missing("`describe'") & missing("`codebook'")) describe `vlist', `options';
					if(!missing("`codebook'") & missing("`describe'")) codebook `vlist',  compact;
					};
			};
		};
	};
	display _newline as text " Total " as result `c_files'
					 as text " out of " as result `n_files' as text " files checked in " _c;
	display `"{stata `"dir "`c(pwd)'\""' : `c(pwd)'\}"';
	if (`c_files' < `n_files') {;
		display _newline as text "(In order to check all `n_files' files in the folder, " _c;
		display "increase memory allocated to the data area)";
	};
end;

