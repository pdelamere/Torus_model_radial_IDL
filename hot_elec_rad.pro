;-----------------------------------------------------------------------  
function nu_ee,nelec_c,nelec_h,Tel,Telh
; m is atomic mass
; z is atomic number
; for collision of given ion on thermal electron population
;-----------------------------------------------------------------------  
me = double(9.1e-28)    ;g
mp = double(1.67e-24)   ;g

lambda = 0.0

if (Tel le 10.0) then begin
   lambda = 23.0 - alog(sqrt(nelec_c)*Tel^(-3./2.))
endif

if (Tel gt 10.0) then begin
   lambda = 24.0 - alog(sqrt(nelec_c)/Tel)
endif

nu_ee = 1.8e-19*sqrt(me*me)*nelec_c*nelec_h*lambda/(me*Tel + me*Telh)^(3./2.)

;print,'nu_ei...',nu_ei,lambda

return,nu_ee
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function nu_ie,mu,z,nelec,ni,Tel,Ti
; m is atomic mass
; z is atomic number
; for collision of given ion on thermal electron population
;-----------------------------------------------------------------------  
me = double(9.1e-28)    ;g
mp = double(1.67e-24)   ;g

mi = mu*mp

lambda = 0.0

if ((Ti*me/mi lt Tel) and (Tel le 10*z^2)) then begin
   lambda = 23.0 - alog(sqrt(nelec)*z*Tel^(-3./2.))
endif

if ((Tel gt Ti*me/mi) and (Tel ge 10*Z^2)) then begin
   lambda = 24.0 - alog(sqrt(nelec)/Tel)
endif

if (Tel lt Ti*Z*me/mi) then begin
   lambda = 30.0 - alog(sqrt(ni)*Ti^(-3./2.)*Z^2/mu)
endif

nu_ie = 1.8e-19*sqrt(me*mi)*(Z^2)*nelec*lambda/(mi*Tel + me*Ti)^(3./2.)

;print,'nu_ei...',nu_ie,lambda

return,nu_ie
end
;-----------------------------------------------------------------------  



;-----------------------------------------------------------------------  
function EF_el,n,T,r
;-----------------------------------------------------------------------  
nel = n.fc*n.nel

if (T.Tel le r.emistemp(1)) then T.Tel = r.emistemp(1)
if (T.Tel_sp le r.emistemp(1)) then T.Tel_sp = r.emistemp(1)
if (T.Tel_s2p le r.emistemp(1)) then T.Tel_s2p = r.emistemp(1)
if (T.Tel_s3p le r.emistemp(1)) then T.Tel_s3p = r.emistemp(1)
if (T.Tel_op le r.emistemp(1)) then T.Tel_op = r.emistemp(1)
if (T.Tel_o2p le r.emistemp(1)) then T.Tel_o2p = r.emistemp(1)

;print,'temps...',T.Tel,T.Tel_sp,T.Tel_s2p,T.Tel_s3p,T.Tel_op,T.Tel_o2p

i=where(abs(r.emistemp-T.Tel_sp) eq min(abs(r.emistemp-T.Tel_sp)))
j=where(abs(r.emisden-nel) eq min(abs(r.emisden-nel)))

i = i(0)
j = j(0)

if (T.Tel_sp gt r.emistemp(i)) then ip = i+1
if (T.Tel_sp le r.emistemp(i)) then begin
   ip = i
   i = ip-1
endif

if (nel gt r.emisden(j)) then jp = j+1
if (nel le r.emisden(j)) then begin
   jp = j
   j = jp-1
endif

dT = r.emistemp(ip) - r.emistemp(i)
dn = r.emisden(jp) - r.emisden(j)

A = dT*dn
dT1 = T.Tel_sp - r.emistemp(i)
dT2 = r.emistemp(ip) - T.Tel_sp
dn1 = nel - r.emisden(j)
dn2 = r.emisden(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSII = r.emisSII(i,j)*A4 + r.emisSII(ip,j)*A3 + r.emisSII(i,jp)*A2 + $
          r.emisSII(ip,jp)*A1
L = emisSII*n.nsp
r.psp = emisSII

i=where(abs(r.emistemp-T.Tel_s2p) eq min(abs(r.emistemp-T.Tel_s2p)))
j=where(abs(r.emisden-nel) eq min(abs(r.emisden-nel)))

i = i(0)
j = j(0)

if (T.Tel_s2p gt r.emistemp(i)) then ip = i+1
if (T.Tel_s2p le r.emistemp(i)) then begin
   ip = i
   i = ip-1
endif

if (nel gt r.emisden(j)) then jp = j+1
if (nel le r.emisden(j)) then begin
   jp = j
   j = jp-1
endif

dT = r.emistemp(ip) - r.emistemp(i)
dn = r.emisden(jp) - r.emisden(j)

A = dT*dn
dT1 = T.Tel_s2p - r.emistemp(i)
dT2 = r.emistemp(ip) - T.Tel_s2p
dn1 = nel - r.emisden(j)
dn2 = r.emisden(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSIII = r.emisSIII(i,j)*A4 + r.emisSIII(ip,j)*A3 + r.emisSIII(i,jp)*A2 + $
          r.emisSIII(ip,jp)*A1
L = L+ emisSIII*n.ns2p
r.ps2p = emisSIII

i=where(abs(r.emistemp-T.Tel_s3p) eq min(abs(r.emistemp-T.Tel_s3p)))
j=where(abs(r.emisden-nel) eq min(abs(r.emisden-nel)))

i = i(0)
j = j(0)

if (T.Tel_s3p gt r.emistemp(i)) then ip = i+1
if (T.Tel_s3p le r.emistemp(i)) then begin
   ip = i
   i = ip-1
endif

if (nel gt r.emisden(j)) then jp = j+1
if (nel le r.emisden(j)) then begin
   jp = j
   j = jp-1
endif

dT = r.emistemp(ip) - r.emistemp(i)
dn = r.emisden(jp) - r.emisden(j)

A = dT*dn
dT1 = T.Tel_s3p - r.emistemp(i)
dT2 = r.emistemp(ip) - T.Tel_s3p
dn1 = nel - r.emisden(j)
dn2 = r.emisden(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSIV = r.emisSIV(i,j)*A4 + r.emisSIV(ip,j)*A3 + r.emisSIV(i,jp)*A2 + $
          r.emisSIV(ip,jp)*A1
L = L + emisSIV*n.ns3p
r.ps3p = emisSIV

i=where(abs(r.emistemp-T.Tel_op) eq min(abs(r.emistemp-T.Tel_op)))
j=where(abs(r.emisden-nel) eq min(abs(r.emisden-nel)))

i = i(0)
j = j(0)

if (T.Tel_op gt r.emistemp(i)) then ip = i+1
if (T.Tel_op le r.emistemp(i)) then begin
   ip = i
   i = ip-1
endif

if (nel gt r.emisden(j)) then jp = j+1
if (nel le r.emisden(j)) then begin
   jp = j
   j = jp-1
endif

dT = r.emistemp(ip) - r.emistemp(i)
dn = r.emisden(jp) - r.emisden(j)

A = dT*dn
dT1 = T.Tel_op - r.emistemp(i)
dT2 = r.emistemp(ip) - T.Tel_op
dn1 = nel - r.emisden(j)
dn2 = r.emisden(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisOII = r.emisOII(i,j)*A4 + r.emisOII(ip,j)*A3 + r.emisOII(i,jp)*A2 + $
          r.emisOII(ip,jp)*A1
L = L + emisOII*n.nop
r.pop = emisOII

i=where(abs(r.emistemp-T.Tel_o2p) eq min(abs(r.emistemp-T.Tel_o2p)))
j=where(abs(r.emisden-nel) eq min(abs(r.emisden-nel)))

i = i(0)
j = j(0)

if (T.Tel_o2p gt r.emistemp(i)) then ip = i+1
if (T.Tel_o2p le r.emistemp(i)) then begin
   ip = i
   i = ip-1
endif

if (nel gt r.emisden(j)) then jp = j+1
if (nel le r.emisden(j)) then begin
   jp = j
   j = jp-1
endif

dT = r.emistemp(ip) - r.emistemp(i)
dn = r.emisden(jp) - r.emisden(j)

A = dT*dn
dT1 = T.Tel_o2p - r.emistemp(i)
dT2 = r.emistemp(ip) - T.Tel_o2p
dn1 = nel - r.emisden(j)
dn2 = r.emisden(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisOIII = r.emisOIII(i,j)*A4 + r.emisOIII(ip,j)*A3 + r.emisOIII(i,jp)*A2 + $
          r.emisOIII(ip,jp)*A1
L = L + emisOIII*n.no2p
r.po2p = emisOIII

r.Puv = L

return,L
end
;-----------------------------------------------------------------------  



;main

Te0 = 40.0
Ti0 = 70.0
Teh0 = 500.0
ns0 = 10.0      ;initial S neutral density (cm^-3)
no0 = 50.0     ;initial O neutral density
nsp0= 250.0        ;initial S+ density
ns2p0 = 400.0       ;initial S++ density
ns3p0 = 50.0      ;initial S+++ density
nop0 = 700.0       ;initial O+ density
no2p0 = 50.0      ;initial O++ density 
nel0 = nsp0+ns2p0+ns3p0+nop0+no2p0       ;initial electron density
nelh0 = 0.01*nel0
trans = 1.0/(40.0*8.64e4)
net_source = 3e28
otos = 3.0
s_source = net_source/(1.0+otos)
o_source = net_source*otos/(otos + 1.0)

fc = 0.02
fh = 1.0-fc

n = {density, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
              no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh, $
              nsp_pu: 0.0, ns2p_pu: 0.0, ns3p_pu: 0.0, nop_pu: 0.0, $
              no2p_pu: 0.0}

T = {temp, Tel: Te0, Telh: Teh0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, Top: Ti0, $
           To2p: Ti0, Tsp_pu: 0.0, Ts2p_pu: 0.0, Ts3p_pu: 0.0, $
           Top_pu: 0.0, To2p_pu: 0.0, Tpu_s: 0.0, Tpu_o: 0.0, $
           Tel_sp: Te0, Tel_s2p: Te0, Tel_s3p: Te0, Tel_op: Te0, $
           Tel_o2p: Te0, Tel_el: Te0}


r = {rates, Tel: Te0, is: 0.0, isp: 0.0, is2p: 0.0, rs3p: 0.0, rs2p: 0.0,$
            io: 0.0, iop:0.0, ro2p:0.0, $
            Telh: Teh0,ish: 0.0, isph: 0.0, is2ph: 0.0, ioh: 0.0, $
            ioph: 0.0, $
            S_production: 0.0, O_production: 0.0, net_production: 0.0, $
            Transport: trans, $
            o_to_s: otos, cx_k0: 0.0, cx_k1: 0.0, cx_k2: 0.0, cx_k3: 0.0, $
            cx_k4: 0.0,$
            cx_k5: 0.0, cx_k6: 0.0, cx_k7: 0.0, cx_k8: 0.0, cx_k9: 0.0, $
            cx_k10: 0.0, cx_k11: 0.0, cx_k12: 0.0, cx_k13: 0.0, cx_k14: 0.0,$
            cx_k15: 0.0, cx_k16: 0.0, emisSII: fltarr(31,41), $
            emisSIII: fltarr(31,41), $
            emisSIV: fltarr(31,41), emisOII: fltarr(31,41), $
            emisOIII: fltarr(31,41), $
            emistemp: fltarr(31), emisden: fltarr(41), $
            Puv: 0.0, psp: 0.0, ps2p: 0.0, ps3p: 0.0, $
            pop: 0.0, po2p: 0.0, Sei: 0.0, Secx: 0.0, Oei: 0.0, Oecx:0.0}

read_emission_tables,'tepSII',temp,den,emisSII & r.emisSII = emisSII
read_emission_tables,'tepSIII',temp,den,emisSIII & r.emisSIII = emisSIII
read_emission_tables,'tepSIV',temp,den,emisSIV & r.emisSIV = emisSIV
read_emission_tables,'tepOII',temp,den,emisOII & r.emisOII = emisOII
read_emission_tables,'tepOIII',temp,den,emisOIII & r.emisOIII = emisOIII
r.emistemp=temp
r.emisden=den



dt = 60.0  ;s

nT0 = n.nel*fc*Te0

print,EF_el(n,T,r),nT0

print,(1./nu_ee(1500.,15.0,5.,50.))/3600.

print,(1./nu_ie(32.,1.,1500.,1500.,5.,100.))/3600./24.
stop

for i = 0,5*60. do begin

nT = nT0 - dt*EF_el(n,T,r)
T.Tel = nT/(n.nel*fc)
T.Tel_sp = T.Tel
T.Tel_s2p = T.Tel
T.Tel_s3p = T.Tel
T.Tel_op = T.Tel
T.Tel_o2p = T.Tel

print,T.Tel,EF_el(n,T,r)

nT0 = nT


endfor

end