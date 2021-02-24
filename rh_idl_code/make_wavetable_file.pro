pro make_wavetable_file, lambdamin, lambdamax, nlambda, filename=filename, wavetable=wavetable, lambdas=lambdas

	if(n_elements(lambdas) eq 0) then begin
		if(n_elements(nlambda) eq 0) then nlambda = 1000
		lambdas = lambdamin+(lambdamax-lambdamin)*findgen(nlambda)/(nlambda-1.0)
	endif else begin
		nlambda = n_elements(lambdas)
	endelse
	wavetable = airtovacuum(lambdas)
	
	if(n_elements(filename) eq 0) then filename = 'wavetable'+strjoin(strtrim(string([lambdamin,lambdamax,nlambda]),2),'_')+'.wave'
	
	openw,1,/xdr,filename
	writeu,1,nlambda,wavetable
	close,1
	
end
