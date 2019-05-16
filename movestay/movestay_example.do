#delim ;
set more off;

use movestay_example, clear;

local str age age2 edu13 edu4 edu5 reg2 reg3 reg4 ;

movestay lmo_wage `str', select(private =m_s1 job_hold);

