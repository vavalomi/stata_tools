program define _gmld
	if strpos("`0'", ",")==0 {
		local t ","
	}
	_ginequal `0' `t' index(mld)
end // program _gmld


