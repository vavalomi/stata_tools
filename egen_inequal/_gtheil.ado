program define _gtheil
	if strpos("`0'", ",")==0 {
		local t ","
	}
	_ginequal `0' `t' index(theil)
end // program _gtheil


