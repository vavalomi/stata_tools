*! version 1.0  04dec2006 Z. Sajaia

program define bestreg, eclass
	version 9.2

	syntax [anything] [if] [in] [aweight pweight iweight fweight/] [,FIXed(varlist) Graph Vars(numlist integer) *]

	tempname tmp tmp1 tmp2 tmp3
	if ~missing("`anything'") {
		_rmcoll `fixed'
		local fixed "`r(varlist)'"
		local k1 : word count `fixed'
		local ++k1

		_rmcollright `fixed' `anything'

		if `:word count `r(block`k1')''>1 {
			display as error "dependent variable missing"
			exit 197
		}

		local k =r(k)
		local ++k1
		local d = 0
		matrix `tmp' = (1\1\0)
		forvalues i=`k1'/`k' {
			local g :word count `r(block`i')'
			matrix `tmp'=`tmp', (`i'-`k1'+1 \ `d' + 1 \ `d'+`g')
			local d = `d'+ `g'
		}

		local vlist "`r(varblocklist)'"
		local vlist : list vlist - fixed
		local vlist : subinstr local vlist "(" "", all
		local vlist : subinstr local vlist ")" "", all

		local k=`k'-`k1'+2

		marksample touse
		markout `touse' `vlist'

		matrix `tmp1'=`tmp'
		local fn : word count `fixed'

		display _newline

		display as input `k'-1 _c
		if `fn'>0 display as text "(" as input "+`fn'" as text ")" _c
		display as text " distinct variables(blocks) in the regression"
		display as input 2^(`k'-1)-1 as text " possible combinations"

		local vlist "`vlist' `fixed'"

		mata : bestreg_mataleaps("`tmp1'")

		display as input cnt as text " regressions performed"
		window manage maintitle reset

		if ~missing("`weight'") local weight "[`weight'=`exp']"
		gettoken depvar vlist : vlist

		forvalues i=2/`k' {
			matrix `tmp2'=J(`tmp'[3,`i']-`tmp'[2,`i']+1, 1, 1)
			matrix `tmp2'=`tmp2' # `tmp1'[`i'-1,....]
			matrix `tmp3'=nullmat(`tmp3') \ `tmp2'
		}
		if `fn'>0 matrix `tmp3'=`tmp3'\J(`:word count `fixed'', colsof(`tmp3'), 1)
		matrix `tmp3'=`tmp3' \ `tmp1'[`k'...,....]

		ereturn post, esample(`touse') properties ()

		ereturn scalar regs   = cnt
		ereturn scalar bestsub= bestcol
		ereturn scalar blocks = `k'-1
		ereturn scalar fixed  =`fn'

		ereturn local cmd    	"bestreg"
		ereturn local criterium "RSS"
		ereturn local depvar 	`depvar'
		ereturn local vlist  	`vlist'
		ereturn local weight 	"`weight'"

		local k=`k'+`fn'-1

		local cnames
		forvalues i=`fn'/`k' {
			local cnames "`cnames' b`i'"
		}
		matrix roweq    `tmp3'=:
		matrix coleq    `tmp3'=:
		matrix colnames `tmp3'=`cnames'
		matrix rownames `tmp3'=`vlist' RSS R2 AdjR2
		ereturn matrix best = `tmp3', copy
	}
	else {
		if "`e(cmd)'"!="bestreg" {
			display as error "last estimates not found"
			exit 301
		}
		local depvar "`e(depvar)'"
		matrix `tmp3'=e(best)
		local vlist "`e(vlist)'"
		local k  = e(blocks)+e(fixed)
		local weight "`e(weight)'"

		scalar bestcol=e(bestsub)
	}

	if missing("`vars'") local vars =bestcol
	tempname tname
	_estimates hold `tname', copy
	foreach i of numlist `vars' {
		local j=0
		local reglist
		foreach var in `vlist' {
			if `tmp3'[`++j',colnumb(`tmp3',"b`i'")]==1 local reglist "`reglist' `var'"
		}
		display _newline
		display as input "Best `i'-variable(block) regression"
		regress `depvar' `reglist' `weight' if e(sample)
	}
	_estimates unhold `tname'

	if ~missing("`graph'") {
		tempvar N R2
		matrix `tmp3'=`tmp3'["R2".."AdjR2",....]'
		svmat `tmp3', names(`R2')
		quietly generate byte `N'=_n-1 if `R2'2<.
		label variable `R2'1 "R-squared"
		label variable `R2'2 "Adjusted R-squared"
		label variable `N'  "Subset size"
		twoway line `R2'1 `N' || line `R2'2 `N', `options'
	}
end
