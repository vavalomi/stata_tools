program define _ggini
	if strpos("`0'", ",")==0 {
		local t ","
	}
	_ginequal `0' `t' index(gini)
end // program _ggini


