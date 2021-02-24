pro sanitize_array, array, nonzero=nonzero, width=width, min=min, max=max
	
	if(n_elements(width) eq 0) then width = 5
	if(n_elements(min) eq 0) then min = -1.0e50
	if(keyword_set(nonzero)) then min = 0.0
	if(n_elements(max) eq 0) then max = 1.0e50
	
	
	wbad = where(array lt min or array gt max or ~finite(array),count)
	if(count gt 0) then begin
		array[wbad] = min
		asm = smooth(array,width,/edge_truncate,missing=0.0,/nan)
		array[wbad] = asm[wbad]
	endif
	
end