preserve
use switch_probit_example, clear
switch_probit works age age2 wedu_2-wedu_5 hhsize hhsize2 reg_*, select(migrant age age2 wedu_2-wedu_5 hhsize hhsize2 reg_* pmigrants)
predict tt, tt
summarize tt if (migrant == 1)
restore
