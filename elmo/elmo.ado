*! version 1.21 15November2023 Z. Sajaia
// version 1.20 6July2010 Z. Sajaia

program elmo, sortpreserve rclass

	version 9.2
	syntax varname(numeric) [aweight fweight] [if/] [in] [, BYgroup(varname numeric) nomax nofail keepzeros]

	local inc "`varlist'"

	tempvar wi
	if missing(`"`weight'"') generate byte `wi' = 1
	else generate `wi' `exp'

	marksample touse
	if ~missing(`"`bygroup'"') markout `touse' `bygroup'


	local cnames "ge0 ge1 ge2"

	quietly {
		if ~missing("`bygroup'") {
			levelsof `bygroup' if `touse', local(lvls)
			if missing("`max'") & (`: word count `lvls'' > 60) {
				noisily display as error "`bygroup' must have no more then 60 categories"
				if missing("`fail'") exit 134
				else exit 0
			}
			foreach lvl of local lvls {
				if int(`lvl') != `lvl' | (`lvl' < 0) {
					noisily display as error "`bygroup' contains non-integer or negative values"
					if missing("`fail'") exit 459
					else exit 0
				}
			}
		}

		count if `inc' < 0 & `touse'
		noisily if r(N) > 0 {
			display " "
			display as txt "Warning: `inc' has `r(N)' values < 0." _c
			display as txt " Not used in calculations"
		}
		replace `touse' = 0 if `inc' < 0

		if "`keepzeros'" == "" {
			count if `inc' == 0 & `touse'
			noisily if r(N) > 0 {
				display " "
				display as txt "Warning: `inc' has `r(N)' values = 0." _c
				display as txt " Not used in calculations"
			}
			replace `touse' = 0 if `inc' == 0
		}

		count if `touse'
		if r(N) == 0 error 2000
		else local N = r(N)

		tempname meany sumwi
		summarize `inc' [w = `wi'] if `touse', meanonly
		scalar `meany' = r(mean)
		scalar `sumwi' = r(sum_w)
		return scalar meany = `meany'

		tempvar fi loginc
		generate double `fi' = `wi' / `sumwi' if `touse'
		generate double `loginc' = log(`inc') if `touse'

		tempname ge0 ge1 ge2
		summarize `loginc' [w = `fi'] if `touse', meanonly
		scalar `ge0' = log(`meany') - r(sum)

		summarize `loginc' [w = `fi' * `inc'] if `touse', meanonly
		scalar `ge1' = r(sum)/`meany' - log(`meany')

		summarize `inc' [w = `fi' * `inc'] if `touse', meanonly
		scalar `ge2' = (r(sum)/`meany'^2 - 1)/2

	 	local rnames "Overall"

		tempname b R F
		matrix `b' = `ge0', `ge1', `ge2'
		matrix `F' = `N', `N', `N'

		matrix colnames `b' = `cnames'
		matrix rownames `b' = `rnames'
		matrix colnames `F' = `cnames'
		matrix rownames `F' = `rnames'

		matrix `R' = `b'

		return matrix overall = `b'
		return matrix overall_F = `F'

		if missing(`"`bygroup'"') {
			noisily matlist `R', format("%10.6f")
			exit
		}

		//========================== SUBGROUP DECOMPOSITIONS

		tempvar fik
		generate double `fik' = .

		tempname mi logmi M F GE0 GE1 GE2 SUMW
		foreach lvl of local lvls {
			summarize `inc' [aw = `wi'] if `touse' & `bygroup' == `lvl', meanonly
			scalar `mi' = r(mean)
			scalar `logmi' = log(`mi')
			matrix `M' = nullmat(`M') \ `mi'
			matrix `SUMW' = nullmat(`SUMW') \ (r(sum_w)/`sumwi')

			replace `fik' = `wi' / r(sum_w) if `touse' & `bygroup' == `lvl'

			summarize `loginc' [aw = `fik'] if `touse' & `bygroup' == `lvl', meanonly
			matrix `GE0' = nullmat(`GE0') \ (`logmi' - r(sum))

			summarize `loginc' [aw = `fik' * `inc'] if `touse' & `bygroup' == `lvl', meanonly
			matrix `GE1' = nullmat(`GE1') \ (r(sum)/`mi' - `logmi')

			summarize `inc' [aw = `fik' * `inc'] if `touse' & `bygroup' == `lvl', meanonly
			matrix `GE2' = nullmat(`GE2') \ ((r(sum) / (`mi'^2) - 1)/2)

			matrix `F' = nullmat(`F') \ r(N)
		}
		return matrix mean = `M', copy
		return matrix sumw = `SUMW', copy

		matrix `b' = `GE0', `GE1', `GE2'
		matrix colnames `b' = `cnames'
		matrix rownames `b' = `lvls'

		matrix `F' = `F', `F', `F'
		matrix colnames `F' = `cnames'
		matrix rownames `F' = `lvls'


		matrix `R' = `R' \ `b'
		local rnames "`rnames' `lvls'"

		return matrix bygroup = `b'
		return matrix bygroup_F = `F'
		return local levels `lvls'

		// ===================== GE index within-group inequalities

		tempname W0 W1 W2

		matrix `W0' = `SUMW'' * `GE0'
		mata : st_matrix("`W1'", st_matrix("`GE1'"):*st_matrix("`M'"))
		matrix `W1' = `SUMW'' * `W1' / `meany'
		mata : st_matrix("`W2'", st_matrix("`GE2'"):*st_matrix("`M'"):*st_matrix("`M'"))
		matrix `W2' = `SUMW'' * `W2' / `meany'^2

		matrix `b' = `W0', `W1', `W2'
		matrix colnames `b' = `cnames'
		matrix rownames `b' = "within"

		matrix `R' = `R' \ `b'
		local rnames "`rnames' within"

		return matrix within = `b'

		// ================== GE index between-group inequalities

		tempvar negtouse
		generate byte `negtouse' = - `touse'

		sort `negtouse' `inc', stable

		capture drop _cummInc _cummPop _sortID
		generate double _cummInc = sum(`inc'*`wi') if `touse'
		generate double _cummPop = sum(`wi') if `touse'

		generate _sortID = _n if `touse'

		tempname logmeany
		scalar `logmeany' = log(`meany')
 		CalcBetweenPart if `touse', meany(`meany') logmeany(`logmeany') sumw(`SUMW') m(`M')

		matrix `b' = r(Result)
		matrix colnames `b' = `cnames'
		matrix rownames `b' = "between"

		matrix `R' = `R' \ `b'
		local rnames "`rnames' between"

		return matrix between = `b'

		// ================== ELMO index between-group inequalities with the original "packing order"

		tempname SORTED_SUMW
		mata : st_matrix("`SORTED_SUMW'", sort((st_matrix("`SUMW'"), st_matrix("`M'")), 2)[., 1])
		CalcBetweenPart if `touse', meany(`meany') logmeany(`logmeany') sumw(`SORTED_SUMW') totw(`sumwi')

		matrix `b' = r(Result)
		matrix colnames `b' = `cnames'
		matrix rownames `b' = "elmo_current"

		matrix `R' = `R' \ `b'
		local rnames "`rnames' elmo_current"
		return matrix elmo_current = `b'

		// =============== ELMO maximum between-group inequality
		if missing("`max'") {
			tempname ELMO_max
			local cmdline "CalcBetweenPart if `touse', meany(`meany') logmeany(`logmeany') sumw(`SORTED_SUMW') totw(`sumwi')"
			mata : mata_ELMOPermutations("`ELMO_max'", "`SUMW'", "`cmdline'", "`SORTED_SUMW'")
			matrix colnames `ELMO_max' = `cnames'
			matrix rownames `ELMO_max' = "elmo_max"

			matrix `R' = `R' \ `ELMO_max'
			local rnames "`rnames' elmo_max"
			return matrix elmo_max = `ELMO_max'
			window manage maintitle reset
		}
	} // quietly

	matrix rownames `R' = `rnames'
	matlist `R'
	return matrix R = `R'
	return scalar N = `N'
	capture drop _cummInc _cummPop _sortID
end

program define CalcBetweenPart, rclass
	syntax [if], meany(string) logmeany(string) sumw(string) [m(string) totw(string)]

	tempname b B0 B1 B2 logM cummPop sumInc

	if missing("`m'") {
		tempname m
		scalar `cummPop' = 0
		scalar  `sumInc' = 0

		forvalues i=1/`=rowsof(`sumw')' {
			scalar `cummPop' = `cummPop' + `sumw'[`i', 1]

			summarize _sortID `if' & _cummPop<=`cummPop'*`totw', meanonly
			local n = r(max)
			tempname inc delta
			if `cummPop' < 0 {
				scalar `delta' = (_cummInc[`=`n'+1'] - _cummInc[`n'])* ///
					(`cummPop'*`totw' - _cummPop[`n']) / (_cummPop[`=`n'+1'] - _cummPop[`n'])
			}
			else scalar `delta' = 0

			scalar `inc' = _cummInc[`n'] + `delta'
			matrix `m' = nullmat(`m') \ (`inc' - `sumInc') / (`sumw'[`i', 1]*`totw')
			scalar `sumInc' = `inc'
		}
	}
	mata: st_matrix("`logM'", log(st_matrix("`m'")))

	mata : st_matrix("`B1'", st_matrix("`m'"):*st_matrix("`logM'"))
	matrix `B1' = `sumw'' * `B1' / `meany' - `logmeany'

	matrix `B0' = `logmeany' - `sumw''*`logM'

	mata : st_matrix("`B2'", st_matrix("`m'"):*st_matrix("`m'"))
	matrix `B2' = (`sumw''*`B2'/(`meany'*`meany') - 1)/2

	matrix `b' = `B0', `B1', `B2'

	return matrix Result = `b'
end // program CalcBetweenPart

mata

void function mata_ELMOPermutations(string scalar matname, string scalar sumw, string scalar cmdline, string scalar sorted_sumw)
	// sumw - matrix name with the original groupweights
	// cmdline will contain something like: CalcBetweenPart `inc' [w] [if], sumw(`sorted_sumw')
	// sorted_sumw - temporary name for the matrix to be created
{
	maxval0 = 0
	maxval1 = 0
	maxval2 = 0

	V1 = st_matrix(sumw)

	info = cvpermutesetup(V1)
	i = 0
	i_max = factorial(rows(st_matrix(sorted_sumw)))

	while ((V1=cvpermute(info)) != J(0,1,.)) {
		st_matrix(sorted_sumw, V1)
		stata(cmdline)
		GE = st_matrix("r(Result)")

		if (maxval0 <= GE[1, 1]) {
			maxval0 = GE[1, 1]
		}
		if (maxval1 <= GE[1, 2]) {
			maxval1 = GE[1, 2]
		}
		if (maxval2 <= GE[1, 3]) {
			maxval2 = GE[1, 3]
		}
		i = i + 1
		p = round(i / i_max * 100, .1)
		stata(`"window manage maintitle ""' + strofreal(i) + " out of " + strofreal(i_max) + " : " + strofreal(p) + `"%""')
	}
	st_matrix(matname, (maxval0, maxval1, maxval2))
}

end
