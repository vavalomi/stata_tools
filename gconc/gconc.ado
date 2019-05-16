program define gconc, eclass sortpreserve properties(svyj svyb)
 version 9.2

 if !replay() {
   syntax varlist(numeric min=1 max=2) [if] [in] [aweight pweight iweight /], ///
      [nu(integer 2) keepnegatives nozeros svy score(str) over(varlist) subpop(passthru) ///
       vce(passthru) robust CLuster(varlist) level(int $S_level) SPECLabel ///
       VALUEMask(str)]

   tempvar pyk tmp B gu

   if !missing("`vce'") + !missing("`robust'") + !missing("`cluster'") + !missing("`svy'") > 1 {
      di as err "only one of vce, robust, cluster or svy can be specified"
      exit 198
   }

   quietly {

      if missing("`exp'") {
          tempvar wi
          gen byte `wi'=1
      }
      else {
          tempvar wi
          generate `wi' = `exp'
      }

      if !missing("`svy'") local svy svy :
      else local pwgt [pw=`wi']

      unab varlist : `varlist'
      if `: word count `varlist'' == 1 local varlist `varlist' `varlist'
      tokenize `varlist'
      local y `1'
      local x `2'

      tempvar touse0
      marksample touse, novarlist
      generate byte `touse0' = `touse'
      markout `touse' `x' `y' `cluster'

      if missing("`keepnegatives'") replace `touse' = 0 if `y' < 0 // negative incomes are excluded by default
      if ~missing("`zeros'") replace `touse' = 0 if `y' == 0       // zeros are by default included !

      count if `touse'
      if r(N) == 0 exit 2000

      tempvar cdfy negcdfynu xbyy

      * compute the cdf of y
      bysort `touse' `over' (`y' `x') : gen double `cdfy' = sum(`wi'*`touse')

      by `touse' `over' : replace `cdfy' = `cdfy' / `cdfy'[_N]

      gen double `negcdfynu' = -( 1 - `cdfy')^(`nu'-1) if `touse'
      gen double `xbyy' = `x'*`negcdfynu' if `touse'

      if !missing("`cluster'") local cluster cluster(`cluster')

      `svy' mean `x' `negcdfynu' `xbyy' if `touse' `pwgt', over( `over', nolabel ) ///
            `subpop' `vce' `robust' `cluster'
      * store some formatting stuff
      local evce `e(vce)'

      local evcetype `e(vcetype)'
      if !missing("`cluster'") {
         local eclustvar `e(clustvar)'
         local ecluster `e(cluster)'
      }
      if !missing("`score'") {
          gen double `score' = `xbyy'/_b[`x'] - _b[`xbyy']/(_b[`x']*_b[`x'])*`x' - `negcdfynu'
      }

      local e_over `e(over)'
      local e_N_over = e(N_over)
      local e_over_namelist `e(over_namelist)'
      local e_over_labels   `"`e(over_labels)'"'
      if e(N_over) == 1 {
         local first  : word 1 of `: colfullnames e(b)'
         local second : word 2 of `: colfullnames e(b)'
         local third  : word 3 of `: colfullnames e(b)'
         qui nlcom (conc: `nu'*(_b[`third']/_b[`first']-_b[`second']) ), post
      }
      else {
         forvalues k=1/`e_N_over' {
            * generate an estimate for the `k'-th group
            local lab : word `k' of `e(over_namelist)'
            * produce a neat... or ugly?.. label for it
            if !missing("`speclabel'") {
               * parse over_labels
               local cats : word `k' of `e(over_labels)'
               local tolab = strtoname("`cats'")
               if length("`valuemask'")>0 {
                  * tolab is _a_b ...
                  tokenize `tolab', parse("_")
                  local tolab
                  while !missing("`1'") {
                     if "`1'" == "_" local tolab `tolab'_
                     else {
                        * replace the last length("`1'") symbols of `valuemask' with `1'
                        local temp = substr("`valuemask'",1,length("`valuemask'")-length("`1'")) + "`1'"
                        local tolab `tolab'`temp'
                     }
                     mac shift
                  }
               }
            }
            else {
               local tolab _`k'
            }
            local tonlcom `tonlcom' (conc`tolab': `nu'*(_b[`xbyy':`lab']/_b[`x':`lab']-_b[`negcdfynu':`lab']) )
         }
         nlcom `tonlcom', post
      }
   }

   ereturn local sortedby `y'
   ereturn local depvar `x'
   if !missing("`over'") ereturn local over `over'
   ereturn scalar nu = `nu'

   ereturn local predict gconc_p
   ereturn local cmd gconc
   ereturn local scorevar `score'

   ereturn local vce `evce'
   ereturn local vcetype `evcetype'
   if !missing("`cluster'") {
      ereturn local clustvar `eclustvar'
      ereturn local cluster `ecluster'
   }

   if !missing("`e_over'") {
      ereturn local over `"`e_over'"'
      ereturn local over_nolabel nolabel
      ereturn scalar N_over = `e_N_over'
      ereturn local over_namelist `e_over_namelist'
      ereturn local over_labels `e_over_labels'
   }

 }
 else { // replay
   if "`e(cmd)'"!="gconc" error 301
   syntax , [Level(cilevel)]
 }

 di _n "{txt}Generalized concentration coefficient"
 di "{txt}Sorted = {res}`e(sortedby)'{txt}; response variable = {res}`e(depvar)'" _c
 di "{txt}; shape parameter = {res}`e(nu)'" _n

 ereturn display , level(`level')

end
