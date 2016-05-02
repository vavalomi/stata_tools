program define _gentropy
	if strpos("`0'", ",")==0 {
		local t ","
	}
	_ginequal `0' `t' index(entropy)
end // program _gentropy


