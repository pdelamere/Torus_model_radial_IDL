pro read_emission_tables,file,temp,den,emis

ntemp = 31
nden = 41

temp = fltarr(31)
den = fltarr(41)
emis = fltarr(ntemp,nden)

close,1
openr,1,file

readf,1,temp
readf,1,den
readf,1,emis

emis = emis/1.6e-12   ;eV/s

;print,temp
;print,den
;surface,emis

return
end
