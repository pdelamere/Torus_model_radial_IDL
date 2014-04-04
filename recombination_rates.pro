; te = electron temperature in eV

pro recombination_rates,T_e,sp,s2p,s3p,op,o2p,s4p

common recombination_rates_common,rrtemp,rrop,rro2p,rrs2p,rrs3p
forward_function dielectronic_recombination

openr,1,'recombination_O.dat'
openr,2,'recombination_s.dat'

readf,1,a1,a2,a3,a4,a5,a6,a7,a8
rrtemp = a1
rrop = a2
rro2p = a3
rro3p = a4
rrov = a5
rrovi = a6
rrovii = a7
rrov2p = a8

while not(eof(1)) do begin
readf,1,a1,a2,a3,a4,a5,a6,a7,a8
rrtemp = [rrtemp,a1]
rrop = [rrop,a2]
rro2p = [rro2p,a3]
rro3p = [rro3p,a4]
rrov = [rrov,a5]
rrovi = [rrovi,a6]
rrovii = [rrovii,a7]
rrov2p = [rrov2p,a8]
endwhile

readf,2,a1,a2
rrs2p = a1
rrs3p = a2

while not(eof(2)) do begin
readf,2,a1,a2
rrs2p = [rrs2p,a1]
rrs3p = [rrs3p,a2]
endwhile

;if not n_elements(rrtemp) ne 81 then begin
;    readcol,'recombination_O.dat',rrtemp,rrop,rro2p,$
;            rro3p,rrov,rrovi,rrovii,rrov2p,/silent,format='F,F,F,F,F,F,F,F'
    
;    readcol,'recombination_s.dat',rrs2p,rrs3p,/silent,$
;            format='F,F'
;endif

k_eV=8.61739e-5
logtemp=alog10(T_e/k_eV)

; s2p = rate for S(III) + e -> S(II) + [gamma]
sp=rrfit(16,16,T_e/k_eV)+dielectronic_recombination(16,16,T_e)

s2p=interpol(rrs2p,rrtemp,logtemp)
s3p=interpol(rrs3p,rrtemp,logtemp)
if n_params() eq 7 then s4p=rrfit(16,13,T_e/k_eV)+$
  dielectronic_recombination(16,13,T_e)

op=interpol(rrop,rrtemp,logtemp)
o2p=interpol(rro2p,rrtemp,logtemp)

close,1
close,2

end


