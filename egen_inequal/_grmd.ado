program define _grmd
	if strpos("`0'", ",")==0 {
		local t ","
	}
	_ginequal `0' `t' index(rmd)
end // program _grmd


