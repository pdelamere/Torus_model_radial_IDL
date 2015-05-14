;--------------------------------------------------------------
function func, p, mr, peuv_peuvh
;--------------------------------------------------------------
@common

sp0 = 0.066
sp1 = 0.023
s2p2 = 0.216
s2p3 = 0.146
s3p4 = 0.0318
s3p5 = 0.0608
op6 = 0.224
op7 = 0.201
o2p8 = 0.0321
o2p9 = 0.0903
sp10 = 0.0415
s3p11 = 0.0390
peuv = 1.45e12
dc_peuv = 0.25e12
nelec_cassini = 2200.
tionavg_Voy = 70.0
dtionavg_Voy = 20.0

;cm3_model,mr,peuv
cm3_radial,p,mr,peuv_peuvh,nelec6,tionavg

densflg = 1

;vol = p(5)

chi = sqrt( $
           ( (sp0-mr(0) )/(0.1*sp0))^2 + $
           ( (sp1-mr(1) )/(0.1*sp1))^2 + $
;           ( (s2p2-mr(2) )/(0.1*s2p2))^2 + $
;           ( (s2p3-mr(3) )/(0.1*s2p3))^2 + $
           ( (s3p4-mr(4) )/(0.1*s3p4))^2 + $
           ( (s3p5-mr(5) )/(0.1*s3p5))^2 + $
;           ( (op6-mr(6) )/(0.1*op6))^2 + $
;           ( (op7-mr(7) )/(0.1*op7))^2 + $
;           ( (o2p8-mr(8) )/(0.1*o2p8))^2 + $
;           ( (o2p9-mr(9) )/(0.1*o2p9))^2 + $
           ( (sp10-mr(10) )/(0.1*sp10))^2 + $
           ( (s3p11-mr(11) )/(0.1*s3p11))^2 + $
           ( (nelec_cassini - nelec6 )/(0.1*nelec_cassini))^2 + $
           ( (tionavg - tionavg_Voy)/(dtionavg_Voy))^2 + $
           ( (peuv-(peuv_peuvh)) /(dc_peuv))^2 )

print,'params........'
print,'  fh..........',p(0)
print,'  fh_alpha....',p(1)
print,'  net_source..',p(2)
print,'  ns_alpha....',p(3)
print,'  DLL.........',p(4)
print,'  DLL_alpha...',p(5)
print,'  Teh.........',p(6)
print,'  Teh_alpha...',p(7)
;print,'  otos........',p(8)
print,' '
print,'nelec........',nelec_cassini,nelec6
print,'Tionavg......',tionavg
print,'Peuv.........',peuv_peuvh
print,'sp...........',mr(0),mr(1)
print,'chi_sqrd.....',chi
print,' '

openw,2,'params1.dat',/append
printf,2,p,mr,nelec6,tionavg,peuv_peuvh,chi
close,2

return,chi
end


;-------------------------------------------------------------
;main
;-------------------------------------------------------------
@common

densflg = 0



mr = fltarr(10)
peuv = 0.0
s = fltarr(8)
p=fltarr(8)

fh = 0.0013           ;fh (L = 6.0)
p(0) = fh
s(0) = 0.001
fh_alpha = 5.0          ;fh power law
p(1) = fh_alpha
s(1) = 2.0

net_source = 1e28    ;net source #/s
p(2) = net_source
s(2) = 2.0e28
net_source_alpha = 14.85  ;power law for extended neutral source
p(3) = net_source_alpha
s(3) = 5

DLL_0 = 9e-7         ;DLL (L = 6.0)
p(4) = DLL_0
s(4) = 4.0e-7
DLL_alpha = 4.4        ;DLL power law
p(5) = DLL_alpha
s(5) = 2.0

Teh0 = 30.0
p(6) = Teh0
s(6) = 10.0
Teh0_alpha = 0.0*4.34 
p(7) = Teh0_alpha
s(7) = 1.0

;otos = 1.8
;p(8) = otos
;s(8) = 0.1

;print,func(p,mr,peuv)
;print,mr,peuv

ftol = 0.001
fval = 0.0

rarr = amoeba(ftol,function_value = fval, nmax=150,P0=p,scale=s)

print,p
print,rarr
print,fval

close,2

end
;-------------------------------------------------------------


