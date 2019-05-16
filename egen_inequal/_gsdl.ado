program define _gsdl
	if strpos("`0'", ",")==0 {
		local t ","
	}
	_ginequal `0' `t' index(sdl)
end // program _gsdl


