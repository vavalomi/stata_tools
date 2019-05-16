program define _gcov
	if strpos("`0'", ",")==0 {
		local t ","
	}
	_ginequal `0' `t' index(cov)
end // program _gcov


