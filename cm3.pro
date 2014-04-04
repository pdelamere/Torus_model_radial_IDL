;-----------------------------------------------------------------------
pro get_rates,r,T
;-----------------------------------------------------------------------
@common

;Electron impact ionization (Schreier and Eviatar 1998)
r.is = 7.4e-8*sqrt(T.Tel)*exp(-10.36/T.Tel)     ;check this!!! 
r.isp = 1.274e-8*sqrt(T.Tel)*exp(-23.1/T.Tel)  
r.is2p = 5.8e-9*sqrt(T.Tel)*exp(-34.88/T.Tel)  

r.ish = 1.464e-6*alog(T.Telh/27.37)/sqrt(T.Telh)
r.isph = 4.44e-7*alog(T.Telh/27.07)/sqrt(T.Telh)
r.is2ph = 1.76e-7*alog(T.Telh/27.07)/sqrt(T.Telh)

r.io = 1.015e-8*sqrt(T.Tel)*exp(-14.54/T.Tel)
r.iop = 3.78e-9*sqrt(T.Tel)*exp(-33.42/T.Tel)

r.ioh = 5.74e-7*alog(T.Telh/27.07)/sqrt(T.Telh)
r.ioph = 1.76e-7*alog(T.Telh/27.07)/sqrt(T.Telh)

;dielectron recombination (Smith and Strobel 1985)
if (T.Tel le 6.9) then r.rs3p = 5.4e-11*T.Tel^0.53
if (T.Tel gt 6.9) then r.rs3p = 4.3e-12*T.Tel^1.84  

r.rs2p = r.rs3p/3.0   ;massive fudge based on Johnson and Strobel '82

r.ro2p = r.rs3p/7.0   ;another massive fudge base on Johnson and Stobel

;charge exchange

r.cx_k0 = 8.1e-9       ;S+ + S++ -> S++ + S+
r.cx_k1 = 2.4e-8
r.cx_k2 = 3.0e-10
r.cx_k3 = 7.8e-9
r.cx_k4 = k4*1.32e-8
r.cx_k5 = 1.32e-8
r.cx_k6 = 5.2e-10
r.cx_k7 = 5.4e-9
r.cx_k8 = 9.0e-11
r.cx_k9 = 3.0e-9
r.cx_k10 = k10*2.34e-8
r.cx_k11 = k11*1.62e-8
r.cx_k12 = 2.3e-9
r.cx_k13 = k13*1.4e-9   ;Schreier for L=6
r.cx_k14 = k14*1.92e-8
r.cx_k15 = k15*9e-9    ;Schreier for L=6

;print,k4,k10,k11,k13,k14,k15

Rj = 7.14e4*1e5 ;cm
a = 5.0*Rj
b = 8.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2

fs = 1.0/(1.0+otos)
fo = otos/(1.0+otos)
;print,fs,fo

r.S_production = fs*net_source/vol
r.O_production = fo*net_source/vol

return
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function F_s,n,r,S,L 
;-----------------------------------------------------------------------
S = r.S_production
L = r.is*n.ns*n.fc*n.nel + r.ish*n.ns*n.fh*n.nel + $ 
    r.cx_k4*n.ns*n.ns3p + r.cx_k11*n.ns*n.no2p + $
    r.cx_k10*n.ns*n.no2p

F_s = S - L

return,F_s
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_sp,n,r,S,L
;-----------------------------------------------------------------------  
S = r.is*n.ns*n.fc*n.nel + r.ish*n.ns*n.fh*n.nel + r.rs2p*n.ns2p*n.nel + $
    r.cx_k4*n.ns*n.ns3p + r.cx_k10*n.ns*n.no2p + 2.0*r.cx_k2*n.ns*n.ns2p + $
    r.cx_k9*n.ns*n.nop + r.cx_k12*n.no*n.ns2p
L = r.isp*n.nsp*n.fc*n.nel + r.cx_k13*n.no2p*n.nsp + n.nsp*r.Transport + $
    r.isph*n.nsp*n.fh*n.nel + r.cx_k8*n.no*n.nsp

F_sp = S - L

return,F_sp
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_s2p,n,r,S,L
;-----------------------------------------------------------------------  
S = r.isp*n.nsp*n.fc*n.nel + r.isph*n.nsp*n.fh*n.nel + $
    r.rs3p*n.ns3p*n.nel + r.cx_k13*n.no2p*n.nsp + $
    r.cx_k14*n.no*n.ns3p + r.cx_k4*n.ns*n.ns3p + $
    r.cx_k11*n.ns*n.no2p
L = r.is2p*n.ns2p*n.fc*n.nel + r.rs2p*n.ns2p*n.nel + $
    r.cx_k15*n.no2p*n.ns2p + n.ns2p*r.Transport + $
    r.is2ph*n.ns2p*n.fh*n.nel

F_s2p = S - L

return,F_s2p
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_s3p,n,r,S,L
;-----------------------------------------------------------------------  
S = r.is2p*n.ns2p*n.fc*n.nel + r.is2ph*n.ns2p*n.fh*n.nel + $
    r.cx_k15*n.no2p*n.ns2p
L = r.rs3p*n.ns3p*n.nel + n.ns3p*r.Transport + r.cx_k14*n.no*n.ns3p + $
    r.cx_k4*n.ns*n.ns3p

F_s3p = S - L

return,F_s3p
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------
function F_o,n,r,S,L
;-----------------------------------------------------------------------
S = r.O_production
L = r.io*n.no*n.fc*n.nel + r.ioh*n.no*n.fh*n.nel + r.cx_k14*n.no*n.ns3p

F_o = S - L

return,F_o
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_op,n,r,S,L
;-----------------------------------------------------------------------  
S = r.io*n.no*n.fc*n.nel + r.ioh*n.no*n.fh*n.nel + $
    r.ro2p*n.no2p*n.nel + r.cx_k13*n.no2p*n.nsp + $
    r.cx_k15*n.no2p*n.ns2p + r.cx_k14*n.no*n.ns3p + $
    r.cx_k11*n.ns*n.no2p + r.cx_k10*n.ns*n.no2p
L = r.iop*n.nop*n.fc*n.nel + n.nop*r.Transport + r.ioph*n.nop*n.fh*n.nel

F_op = S - L

return,F_op
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_o2p,n,r,S,L
;-----------------------------------------------------------------------  
S = r.iop*n.nop*n.fc*n.nel + r.ioph*n.nop*n.fh*n.nel
L = r.ro2p*n.no2p*n.nel + r.cx_k13*n.no2p*n.nsp + $
    r.cx_k15*n.no2p*n.ns2p + n.no2p*r.Transport + $
    r.cx_k11*n.ns*n.no2p + r.cx_k10*n.ns*n.no2p

F_o2p = S - L

return,F_o2p
end
;-----------------------------------------------------------------------  

;=======================Energy source/loss functions=====================


;-----------------------------------------------------------------------  
function Tpu,m,L
;input mass as atomic number
;-----------------------------------------------------------------------  
mp = 1.67e-27  ;kg

vrel = 12.5*L - 42.0/sqrt(L)   ;km/s
vrel = vrel*1e3  ;m/s

Tpu = 0.5*m*mp*vrel^2   ;J
Tpu = Tpu/1.6e-19       ;eV

return,Tpu
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
   lambda = 23.0 - alog(sqrt(nelec)*Tel^(-3./2.))
endif

if ((Tel gt Ti*me/mi) and (Tel ge 10*Z^2)) then begin
   lambda = 24.0 - alog(sqrt(nelec)/Tel)
endif

if (Tel lt Ti*Z*me/mi) then begin
   lambda = 30.0 - alog(sqrt(ni)*Ti^(-3./2.)*Z^2/mu)
endif

nu_ie = 1.8e-19*sqrt(me*mi)*(Z^2)*nelec*lambda/(mi*Tel + me*Ti)^(3./2.)

;print,'nu_ei...',nu_ei,lambda

return,nu_ie
end
;-----------------------------------------------------------------------  


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
function nu_ii,mu1,mu2,z1,z2,n1,n2,T1,T2
; m is atomic mass
; z is atomic number
; for collision of given ion on ion population
;-----------------------------------------------------------------------  
mp = 1.67e-24   ;g

m1 = double(mu1*mp)
m2 = double(mu2*mp)

lambda = 23.0 - alog((Z1*Z2*(mu1+mu2)/(mu1*T2 + mu2*T1))* $
                     sqrt((n1*Z1^2)/T1 + (n2*Z2^2)/T2))

nu_ii = 1.8e-19*sqrt(m1*m2)*(Z1^2)*(Z2^2)*n2*lambda/(m1*T2 + m2*T1)^(3./2.)

return,nu_ii
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro update_temp,n,nT,T
;-----------------------------------------------------------------------  

T.Tsp = nT.sp/n.nsp
T.Ts2p = nT.s2p/n.ns2p
T.Ts3p = nT.s3p/n.ns3p
T.Top = nT.op/n.nop
T.To2p = nT.o2p/n.no2p
T.Tel = nT.el/n.nel

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function EF_el,n,T,nT,r
;-----------------------------------------------------------------------  
nel = n.fc*n.nel

i=where(abs(r.emistemp-T.Tel) eq min(abs(r.emistemp-T.Tel)))
j=where(abs(r.emisden-nel) eq min(abs(r.emisden-nel)))

i = i(0)
j = j(0)

if (T.Tel gt r.emistemp(i)) then ip = i+1
if (T.Tel le r.emistemp(i)) then begin
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
dT1 = T.Tel - r.emistemp(i)
dT2 = r.emistemp(ip) - T.Tel
dn1 = nel - r.emisden(j)
dn2 = r.emisden(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

;print,'weights...',[A1,A2,A3,A4]/A,total([A1,A2,A3,A4]/A)

;i=where(abs(r.emistemp-T.Tsp) eq min(abs(r.emistemp-T.Tsp)))
;j=where(abs(r.emisden-n.nsp) eq min(abs(r.emisden-n.nsp)))
;L = r.emisSII(i,j)*n.fc*n.nel*n.nsp
emisSII = r.emisSII(i,j)*A4 + r.emisSII(ip,j)*A3 + r.emisSII(i,jp)*A2 + $
          r.emisSII(ip,jp)*A1
L = emisSII*n.nsp
r.psp = emisSII
;L = r.emisSII(i,j)*n.nsp

;i=where(abs(r.emistemp-T.Ts2p) eq min(abs(r.emistemp-T.Ts2p)))
;j=where(abs(r.emisden-n.ns2p) eq min(abs(r.emisden-n.ns2p)))
;L = L + r.emisSIII(i,j)*n.fc*n.nel*n.ns2p
emisSIII = r.emisSIII(i,j)*A4 + r.emisSIII(ip,j)*A3 + r.emisSIII(i,jp)*A2 + $
          r.emisSIII(ip,jp)*A1
L = L+ emisSIII*n.ns2p
r.ps2p = emisSIII
;L = L + r.emisSIII(i,j)*n.ns2p

;i=where(abs(r.emistemp-T.Ts3p) eq min(abs(r.emistemp-T.Ts3p)))
;j=where(abs(r.emisden-n.ns3p) eq min(abs(r.emisden-n.ns3p)))
;L = L + r.emisSIV(i,j)*n.fc*n.nel*n.ns3p
emisSIV = r.emisSIV(i,j)*A4 + r.emisSIV(ip,j)*A3 + r.emisSIV(i,jp)*A2 + $
          r.emisSIV(ip,jp)*A1
L = L + emisSIV*n.ns3p
r.ps3p = emisSIV
;L = L + r.emisSIV(i,j)*n.ns3p

;i=where(abs(r.emistemp-T.Top) eq min(abs(r.emistemp-T.Top)))
;j=where(abs(r.emisden-n.nop) eq min(abs(r.emisden-n.nop)))
;L = L + r.emisOII(i,j)*n.fc*n.nel*n.nop
emisOII = r.emisOII(i,j)*A4 + r.emisOII(ip,j)*A3 + r.emisOII(i,jp)*A2 + $
          r.emisOII(ip,jp)*A1
L = L + emisOII*n.nop
r.pop = emisOII
;L = L + r.emisOII(i,j)*n.nop

;i=where(abs(r.emistemp-T.To2p) eq min(abs(r.emistemp-T.To2p)))
;j=where(abs(r.emisden-n.no2p) eq min(abs(r.emisden-n.no2p)))
;L = L + r.emisOIII(i,j)*n.fc*n.nel*n.no2p
emisOIII = r.emisOIII(i,j)*A4 + r.emisOIII(ip,j)*A3 + r.emisOIII(i,jp)*A2 + $
          r.emisOIII(ip,jp)*A1
L = L + emisOIII*n.no2p
r.po2p = emisOIII
;L = L + r.emisOIII(i,j)*n.no2p

Teq = nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel,T.Tsp)*n.nsp*(T.Tsp-T.Tel) + $
      nu_ie(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel,T.Ts2p)*n.ns2p*(T.Ts2p-T.Tel) + $
      nu_ie(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel,T.Ts3p)*n.ns3p*(T.Ts3p-T.Tel) + $
      nu_ie(16.0,1.0,n.fc*n.nel,n.nop,T.Tel,T.Top)*n.nop*(T.Top-T.Tel) + $
      nu_ie(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel,T.To2p)*n.no2p*(T.To2p-T.Tel) + $
      nu_ee(n.fc*n.nel,n.fh*n.nel,T.Tel,T.Telh)*(T.Telh - T.Tel)

;print,'nu_ie...',nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel,T.Tsp)*n.nsp*(T.Tsp-T.Tel),$
;nu_ie(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel,T.Ts2p)*n.ns2p*(T.Ts2p-T.Tel),$
;nu_ie(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel,T.Ts3p)*n.ns3p*(T.Ts3p-T.Tel), $
;nu_ie(16.0,1.0,n.fc*n.nel,n.nop,T.Tel,T.Top)*n.nop*(T.Top-T.Tel), $
;nu_ie(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel,T.To2p)*n.no2p*(T.To2p-T.Tel), $
;nu_ee(n.fc*n.nel,n.fh*n.nel,T.Tel,T.Telh)*(T.Telh - T.Tel),L

;print,L,Teq,nu_ee(n.fc*n.nel,n.fh*n.nel,T.Tel,T.Telh)*(T.Telh - T.Tel)

r.Puv = L/(n.nsp + n.ns2p + n.ns3p + n.nop + n.no2p)

EF_el = Teq - L/1.0

return,EF_el
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function EF_sp,n,T,nT,r
;-----------------------------------------------------------------------  
S = (r.is*n.ns*n.fc*n.nel + $        ;cold e ionization 
     r.ish*n.ns*n.fh*n.nel + $       ;hot e ionization
     r.cx_k4*n.ns*n.ns3p + $         ;S + S+++ -> S+ + S++
     r.cx_k10*n.ns*n.no2p + $        ;S + O++ -> S+ + O+
     r.cx_k1*n.ns*n.nsp + $          ;S + S+ -> S+ + S
     r.cx_k2*n.ns*n.ns2p + $         ;S + S++ -> S+ + S+ 
     r.cx_k9*n.ns*n.nop)*T.Tpu_s + $ ;S + O+ -> S+ + O
     r.cx_k2*n.ns*n.ns2p*T.Ts2p + $  ;S + S++ -> S+ + S+
     r.cx_k12*n.no*n.ns2p*T.Ts2p + $ ;O + S++ -> O+ + S+
     r.cx_k0*n.nsp*n.ns2p*T.Ts2p     ;S+ + S++ -> S++ + S+

L= r.cx_k1*n.ns*n.nsp*T.Tsp + $     ;S + S+ -> S+ + S
   r.cx_k8*n.no*n.nsp*T.Tsp + $     ;O + S+ -> O+ + S
   r.cx_k13*n.no2p*n.nsp*T.Tsp + $  ;O++ + S+ -> O+ + S++
   r.cx_k0*n.nsp*n.ns2p*T.Tsp + $   ;S+ + S++ -> S++ + S+
   r.Transport*n.nsp*T.Tsp    

Teq= nu_ii(32.0,32.0,1.0,2.0,n.nsp,n.ns2p,T.Tsp,T.Ts2p)*n.nsp*(T.Ts2p-T.Tsp)+$
   nu_ii(32.0,32.0,1.0,3.0,n.nsp,n.ns3p,T.Tsp,T.Ts3p)*n.nsp*(T.Ts3p-T.Tsp)+$
   nu_ii(32.0,16.0,1.0,1.0,n.nsp,n.nop,T.Tsp,T.Top)*n.nsp*(T.Top-T.Tsp)+$
   nu_ii(32.0,16.0,1.0,2.0,n.nsp,n.no2p,T.Tsp,T.To2p)*n.nsp*(T.To2p-T.Tsp)+$
   nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel,T.Tsp)*n.nsp*(T.Tel-T.Tsp) + $
   nu_ie(32.0,1.0,n.fh*n.nel,n.nsp,T.Telh,T.Tsp)*n.nsp*(T.Telh-T.Tsp) 
 
EF_sp = S - L + Teq

return,EF_sp
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function EF_s2p,n,T,nT,r
;-----------------------------------------------------------------------  
S = (r.isp*n.nsp*n.fc*n.nel + $               ;cold e ionization 
     r.isph*n.nsp*n.fh*n.nel)*T.Tsp + $       ;hot e ionization
     r.cx_k3*n.ns*n.ns2p*T.Tpu_s + $          ;S + S++ -> S++ + S
     r.cx_k4*n.ns*n.ns3p*T.Ts3p + $           ;S + S+++ -> S+ + S++
     r.cx_k11*n.ns*n.no2p*T.Tpu_s + $         ;S + O++ -> S++ + O+ + e
     r.cx_k13*n.no2p*n.nsp*T.Tsp + $          ;O++ + S+ -> O+ + S++ 
     r.cx_k14*n.no*n.ns3p*T.Ts3p + $          ;O + S+++ -> O+ + S++
     r.cx_k0*n.nsp*n.ns2p*T.Tsp               ;S+ + S++ -> S++ + S+

L= r.cx_k2*n.ns*n.ns2p*T.Ts2p + $     ;S + S++ -> S+ + S+
   r.cx_k3*n.ns*n.ns2p*T.Ts2p + $     ;S + S++ -> S++ + S
   r.cx_k12*n.no*n.ns2p*T.Ts2p + $    ;O + S++ -> O+ + S+
   r.cx_k15*n.no2p*n.ns2p*T.Ts2p + $  ;O++ + S++ -> O+ + S+++
   r.cx_k0*n.nsp*n.ns2p*T.Ts2p + $    ;S+ + S++ -> S++ + S+
   r.Transport*n.ns2p*T.Ts2p    

Teq= nu_ii(32.0,32.0,2.0,1.0,n.ns2p,n.nsp,T.Ts2p,T.Tsp)*n.ns2p*(T.Tsp-T.Ts2p)+$
   nu_ii(32.0,32.0,2.0,3.0,n.ns2p,n.ns3p,T.Ts2p,T.Ts3p)*n.ns2p*(T.Ts3p-T.Ts2p)+$
   nu_ii(32.0,16.0,2.0,1.0,n.ns2p,n.nop,T.Ts2p,T.Top)*n.ns2p*(T.Top-T.Ts2p)+$
   nu_ii(32.0,16.0,2.0,2.0,n.ns2p,n.no2p,T.Ts2p,T.To2p)*n.ns2p*(T.To2p-T.Ts2p)+$
   nu_ie(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel,T.Ts2p)*n.ns2p*(T.Tel-T.Ts2p) + $
   nu_ie(32.0,2.0,n.fh*n.nel,n.ns2p,T.Telh,T.Ts2p)*n.ns2p*(T.Telh-T.Ts2p) 
 
EF_s2p = S - L + Teq

return,EF_s2p
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function EF_s3p,n,T,nT,r
;-----------------------------------------------------------------------  
S = (r.is2p*n.ns2p*n.fc*n.nel + $               ;cold e ionization 
     r.is2ph*n.ns2p*n.fh*n.nel)*T.Ts2p + $      ;hot e ionization
     r.cx_k15*n.no2p*n.ns2p*T.Ts2p              ;O++ + S++ -> O+ + S+++

L= r.cx_k14*n.no*n.ns3p*T.Ts3p + $      ;O + S+++ -> O+ + S++
   r.Transport*n.ns3p*T.Ts3p    

Teq= nu_ii(32.0,32.0,3.0,1.0,n.ns3p,n.nsp,T.Ts3p,T.Tsp)*n.ns3p*(T.Tsp-T.Ts3p)+$
   nu_ii(32.0,32.0,3.0,2.0,n.ns3p,n.ns2p,T.Ts3p,T.Ts2p)*n.ns3p*(T.Ts2p-T.Ts3p)+$
   nu_ii(32.0,16.0,3.0,1.0,n.ns3p,n.nop,T.Ts3p,T.Top)*n.ns3p*(T.Top-T.Ts3p)+$
   nu_ii(32.0,16.0,3.0,2.0,n.ns3p,n.no2p,T.Ts3p,T.To2p)*n.ns3p*(T.To2p-T.Ts3p)+$
   nu_ie(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel,T.Ts3p)*n.ns3p*(T.Tel-T.Ts3p) + $
   nu_ie(32.0,3.0,n.fh*n.nel,n.ns3p,T.Telh,T.Ts3p)*n.ns3p*(T.Telh-T.Ts3p) 
 
EF_s3p = S - L + Teq

return,EF_s3p
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function EF_op,n,T,nT,r
;-----------------------------------------------------------------------  
S = (r.io*n.no*n.fc*n.nel + $          ;cold e ionization 
     r.ioh*n.no*n.fh*n.nel + $         ;hot e ionization
     r.cx_k5*n.no*n.nop + $            ;O + O+ -> O+ + O
     r.cx_k6*n.no*n.no2p + $           ;O + O++ -> O+ + O+
     r.cx_k8*n.no*n.nsp + $            ;O + S+ -> O+ + S
     r.cx_k12*n.no*n.ns2p + $          ;O + S++ -> O+ + S+ 
     r.cx_k14*n.no*n.ns3p)*T.Tpu_o + $ ;O + S+++ -> O+ + S++
     r.cx_k6*n.no*n.no2p*T.To2p + $    ;O + O++ -> O+ + O+
     r.cx_k10*n.ns*n.no2p*T.To2p + $   ;S + O++ -> S+ + O+
     r.cx_k11*n.ns*n.no2p*T.To2p + $   ;S + O++ -> S+ + O+ + e
     r.cx_k13*n.no2p*n.nsp*T.To2p + $  ;O++ + S+ -> O+ + S++
     r.cx_k15*n.no2p*n.ns2p*T.To2p     ;O++ + S++ -> O+ + S+++


L= r.cx_k5*n.no*n.nop*T.Top + $     ;O + O+ -> O+ + O
   r.cx_k9*n.ns*n.nop*T.Top + $     ;S + O+ -> S+ + O
   r.Transport*n.nop*T.Top    

Teq= nu_ii(16.0,16.0,1.0,2.0,n.nop,n.no2p,T.Top,T.To2p)*n.nop*(T.To2p-T.Top)+$
   nu_ii(16.0,32.0,1.0,1.0,n.nop,n.nsp,T.Top,T.Tsp)*n.nop*(T.Tsp-T.Top)+$
   nu_ii(16.0,32.0,1.0,2.0,n.nop,n.ns2p,T.Top,T.Ts2p)*n.nop*(T.Ts2p-T.Top)+$
   nu_ii(16.0,32.0,1.0,3.0,n.nop,n.ns3p,T.Top,T.Ts3p)*n.nop*(T.Ts3p-T.Top)+$
   nu_ie(16.0,1.0,n.fc*n.nel,n.nop,T.Tel,T.Top)*n.nop*(T.Tel-T.Top) + $
   nu_ie(16.0,1.0,n.fh*n.nel,n.nop,T.Telh,T.Top)*n.nop*(T.Telh-T.Top) 
 
EF_op = S - L + Teq

return,EF_op
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function EF_o2p,n,T,nT,r
;-----------------------------------------------------------------------  
S = (r.iop*n.nop*n.fc*n.nel + $                  ;cold e ionization 
     r.ioph*n.nop*n.fh*n.nel)*T.Top + $          ;hot e ionization
     r.cx_k7*n.no*n.no2p*T.Tpu_o             ;O + O++ -> O++ + O

L= r.cx_k6*n.no*n.no2p*T.To2p + $     ;O + O++ -> O+ + O+
   r.cx_k7*n.no*n.no2p*T.To2p + $     ;O + O++ -> O++ + O
   r.cx_k10*n.ns*n.no2p*T.To2p + $    ;S + O++ -> S+ + O+
   r.cx_k11*n.ns*n.no2p*T.To2p + $    ;S + O++ -> S++ + O+ + e
   r.cx_k13*n.no2p*n.nsp*T.To2p + $   ;O++ + S+ -> O+ + S++
   r.cx_k15*n.no2p*n.ns2p*T.To2p + $  ;O++ + S++ -> O+ + S+++
   r.Transport*n.no2p*T.To2p    

Teq= nu_ii(16.0,16.0,2.0,1.0,n.no2p,n.nop,T.To2p,T.Top)*n.no2p*(T.Top-T.To2p)+$
   nu_ii(16.0,32.0,2.0,1.0,n.no2p,n.nsp,T.To2p,T.Tsp)*n.no2p*(T.Tsp-T.To2p)+$
   nu_ii(16.0,32.0,2.0,2.0,n.no2p,n.ns2p,T.To2p,T.Ts2p)*n.no2p*(T.Ts2p-T.To2p)+$
   nu_ii(16.0,32.0,2.0,3.0,n.no2p,n.ns3p,T.To2p,T.Ts3p)*n.no2p*(T.Ts3p-T.To2p)+$
   nu_ie(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel,T.To2p)*n.no2p*(T.Tel-T.To2p) + $
   nu_ie(16.0,2.0,n.fh*n.nel,n.no2p,T.Telh,T.To2p)*n.no2p*(T.Telh-T.To2p) 
 
EF_o2p = S - L + Teq

return,EF_o2p
end
;-----------------------------------------------------------------------  



;-----------------------------------------------------------------------
pro cm3_model,n,r,T
;main program
;cubic centimeter torus chemistry model
;-----------------------------------------------------------------------
@common

;ns0 = 1.0      ;initial S neutral density (cm^-3)
;no0 = 2.0*ns0     ;initial O neutral density
;nsp0= 1.0        ;initial S+ density
;ns2p0 = 1.0       ;initial S++ density
;ns3p0 = 1.0      ;initial S+++ density
;nop0 = 1.0       ;initial O+ density
;no2p0 = 1.0      ;initial O++ density 
nel0 = nsp0+ns2p0+ns3p0+nop0+no2p0       ;initial electron density
Te0 = float(Te0)         ;initial electron temperature (eV)
dt = dt0          ;time step (seconds)
ntm = fix(runt/dt)

fc = 1.0-fh

n = {density, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
              no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh, $
              nsp_pu: 0.0, ns2p_pu: 0.0, ns3p_pu: 0.0, nop_pu: 0.0, $
              no2p_pu: 0.0}
T = {temp, Tel: Te0, Telh: Teh0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, Top: Ti0, $
           To2p: Ti0, Tsp_pu: 0.0, Ts2p_pu: 0.0, Ts3p_pu: 0.0, $
           Top_pu: 0.0, To2p_pu: 0.0, Tpu_s: 0.0, Tpu_o: 0.0}
nT = {energy, el: nel0*Te0, sp: nsp0*Ti0, s2p: ns2p0*Ti0, $
              s3p: ns3p0*Ti0, op: nop0*Ti0, o2p: no2p0*Ti0, $
              sp_pu: 0.0, s2p_pu: 0.0, s3p_pu: 0.0, op_pu: 0.0,$
              o2p_pu: 0.0} 

;advance density for improved Euler method, half time step
n1 = {density_1, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
                 no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh, $
              nsp_pu: 0.0, ns2p_pu: 0.0, ns3p_pu: 0.0, nop_pu: 0.0, $
              no2p_pu: 0.0}
T1 = {temp_1, Tel: Te0, Telh: Teh0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, Top: Ti0, $
           To2p: Ti0, Tsp_pu: 0.0, Ts2p_pu: 0.0, Ts3p_pu: 0.0, $
           Top_pu: 0.0, To2p_pu: 0.0, Tpu_s: 0.0, Tpu_o: 0.0}
nT1 = {energy_1, el: nel0*Te0, sp: nsp0*Ti0, s2p: ns2p0*Ti0, $
              s3p: ns3p0*Ti0, op: nop0*Ti0, o2p: no2p0*Ti0, $
              sp_pu: 0.0, s2p_pu: 0.0, s3p_pu: 0.0, op_pu: 0.0,$
              o2p_pu: 0.0} 

;density at full time step advance
np = {density_p, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
                 no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh, $
              nsp_pu: 0.0, ns2p_pu: 0.0, ns3p_pu: 0.0, nop_pu: 0.0, $
              no2p_pu: 0.0}
Tp = {temp_p, Tel: Te0, Telh: Teh0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, Top: Ti0, $
           To2p: Ti0, Tsp_pu: 0.0, Ts2p_pu: 0.0, Ts3p_pu: 0.0, $
           Top_pu: 0.0, To2p_pu: 0.0, Tpu_s: 0.0, Tpu_o: 0.0}
nTp = {energy_p, el: nel0*Te0, sp: nsp0*Ti0, s2p: ns2p0*Ti0, $
              s3p: ns3p0*Ti0, op: nop0*Ti0, o2p: no2p0*Ti0, $
              sp_pu: 0.0, s2p_pu: 0.0, s3p_pu: 0.0, op_pu: 0.0,$
              o2p_pu: 0.0} 

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
            cx_k15: 0.0, emisSII: fltarr(31,41), emisSIII: fltarr(31,41), $
            emisSIV: fltarr(31,41), emisOII: fltarr(31,41), $ 
            emisOIII: fltarr(31,41), $
            emistemp: fltarr(31), emisden: fltarr(41), $
            Puv: 0.0, psp: 0.0, ps2p: 0.0, ps3p: 0.0, $
            pop: 0.0, po2p: 0.0}

src = replicate({source_func, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0,$
              o: 0.0, op: 0.0, o2p: 0.0,tm:0.0},ntm+1)
lss = replicate({loss_func, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0,$
              o: 0.0, op: 0.0, o2p: 0.0,tm:0.0},ntm+1)

p = replicate({Power_uv, Puv: 0.0, psp: 0.0, ps2p: 0.0, ps3p: 0.0, $
                         pop: 0.0, po2p: 0.0},ntm+1)

temps = replicate({temperatures, tsp: 0.0, Ts2p: 0.0, Ts3p: 0.0, Top: 0.0, $
                                 To2p: 0.0, Telec: 0.0},ntm+1)

dens = replicate({dens, nsp: 0.0, ns2p: 0.0, ns3p: 0.0, nop: 0.0, $
                  no2p: 0.0, nel: 0.0},ntm+1)

;read emission tables
read_emission_tables,'tepSII',temp,den,emisSII & r.emisSII = emisSII
read_emission_tables,'tepSIII',temp,den,emisSIII & r.emisSIII = emisSIII
read_emission_tables,'tepSIV',temp,den,emisSIV & r.emisSIV = emisSIV
read_emission_tables,'tepOII',temp,den,emisOII & r.emisOII = emisOII
read_emission_tables,'tepOIII',temp,den,emisOIII & r.emisOIII = emisOIII
r.emistemp=temp
r.emisden=den

;set initial pickup temperatures

T.Tpu_s = Tpu(32.0,6.0)
T.Tpu_o = Tpu(16.0,6.0)
T1 = T
Tp = T


if (plot_density eq 1) then begin
plot_io,[0,ntm*dt]/8.64e4,[1.0,no0+ns0],/nodata,$
  title='Density',$
  xtitle='time (days) ',ytitle='n (cm!u-3!n)',yrange=[0.1,10000.0],/ysty,$
  xrange=[0,dt*ntm]/8.64e4,/xsty
endif

if (plot_temperature eq 1) then begin
plot_io,[0,ntm*dt]/8.64e4,[0.1,2000.0],/nodata,$
  title='Temperature',$
  xtitle='time (days) ',ytitle='n (cm!u-3!n)',yrange=[0.1,2000.0],/ysty,$
  xrange=[0,dt*ntm]/8.64e4,/xsty
endif

if (strtrim(!d.name,2) ne 'PS') then device,decompose=0
loadct,27
get_rates,r,T

if (cont eq 'true') then begin
   restore,'restart_cm3.sav'
   fc = 1.0 - fh
   n.fc = fc
   n.fh = fh
   n1.fc = fc
   n1.fh = fh
   np.fc = fc
   np.fh = fh
;   r.transport = trans
;   r.otos = otos
;   T.Tel = Te0
   T.Telh = Teh0
   T1 = T
   Tp = T
endif


for i = 0,ntm do begin
   
   get_rates,r,T

   ;full time step advance
   n1.ns = n.ns + dt*F_s(n,r,S,L)
   src(i).s = S
   lss(i).s = L
   n1.nsp = n.nsp + dt*F_sp(n,r,S,L)
   src(i).sp = S
   lss(i).sp = L
   n1.ns2p = n.ns2p + dt*F_s2p(n,r,S,L)
   src(i).s2p = S
   lss(i).s2p = L
   n1.ns3p = n.ns3p + dt*F_s3p(n,r,S,L)
   src(i).s3p = S
   lss(i).s3p = L
   n1.no = n.no + dt*F_o(n,r,S,L)
   src(i).o = S
   lss(i).o = L
   n1.nop = n.nop + dt*F_op(n,r,S,L)
   src(i).op = S
   lss(i).op = L 
   n1.no2p = n.no2p + dt*F_o2p(n,r,S,L)
   src(i).o2p = S
   lss(i).o2p = L
   src(i).tm = dt*i/8.64e4
   lss(i).tm = dt*i/8.64e4

   n1.nel = 1.0*n1.nsp+2.0*n1.ns2p+3.0*n1.ns3p+1.0*n1.nop+2.0*n1.no2p

   ;full time step advance for energy
   nT1.sp = nT.sp + dt*EF_sp(n,T,nT,r)
   nT1.s2p = nT.s2p + dt*EF_s2p(n,T,nT,r)
   nT1.s3p = nT.s3p + dt*EF_s3p(n,T,nT,r)
   nT1.op = nT.op + dt*EF_op(n,T,nT,r)
   nT1.o2p = nT.o2p + dt*EF_o2p(n,T,nT,r)
   nT1.el = nT.el + dt*EF_el(n,T,nT,r)
   update_temp,n1,nT1,T1
   ;Improved Euler advance

   np.ns = n.ns + dt*0.5*(F_s(n,r)+F_s(n1,r))
   np.nsp = n.nsp + dt*0.5*(F_sp(n,r)+F_sp(n1,r))
   np.ns2p = n.ns2p + dt*0.5*(F_s2p(n,r)+F_s2p(n1,r))
   np.ns3p = n.ns3p + dt*0.5*(F_s3p(n,r)+F_s3p(n1,r))
   np.no = n.no + dt*0.5*(F_o(n,r)+F_o(n1,r))
   np.nop = n.nop + dt*0.5*(F_op(n,r)+F_op(n1,r))
   np.no2p = n.no2p + dt*0.5*(F_o2p(n,r)+F_o2p(n1,r))

   np.nel = 1.0*np.nsp+2.0*np.ns2p+3.0*np.ns3p+1.0*np.nop+2.0*np.no2p

   nTp.sp = nT.sp + dt*0.5*(EF_sp(n1,T1,nT1,r)+ EF_sp(n,T,nT,r))
   nTp.s2p = nT.s2p + dt*0.5*(EF_s2p(n1,T1,nT1,r)+ EF_s2p(n,T,nT,r))
   nTp.s3p = nT.s3p + dt*0.5*(EF_s3p(n1,T1,nT1,r)+ EF_s3p(n,T,nT,r))
   nTp.op = nT.op + dt*0.5*(EF_op(n1,T1,nT1,r)+ EF_op(n,T,nT,r))
   nTp.o2p = nT.o2p + dt*0.5*(EF_o2p(n1,T1,nT1,r)+ EF_o2p(n,T,nT,r))
   nTp.el = nT.el + dt*0.5*(EF_el(n1,T1,nT1,r)+ EF_el(n,T,nT,r))
   update_temp,np,nTp,Tp

   n = np
   nT = nTp  
   T = Tp
   temps(i).Tsp = T.Tsp
   temps(i).Ts2p = T.Ts2p
   temps(i).Ts3p = T.Ts3p
   temps(i).Top = T.Top
   temps(i).To2p = T.To2p
   temps(i).Telec = T.Tel
   dens(i).nsp = n.nsp
   dens(i).ns2p = n.ns2p
   dens(i).ns3p = n.ns3p
   dens(i).nop = n.nop
   dens(i).no2p = n.no2p
   dens(i).nel = n.nel

;   print,T.Tsp,T.Ts2p,T.Ts3p,T.Top,T.To2p,T.Tel

;update Pickup Densities
  
   n.nsp_pu = nT.sp_pu/T.Tsp_pu
   n.ns2p_pu = nT.s2p_pu/T.Ts2p_pu
   n.ns3p_pu = nT.s3p_pu/T.Ts3p_pu
   n.nop_pu = nT.op_pu/T.Top_pu
   n.no2p_pu = nT.o2p_pu/T.To2p_pu

   if (plot_density eq 1) then begin
   ;   plots,i*dt/8.64e4,n.ns,psym=3,color=10
      plots,i*dt/8.64e4,n.nsp,psym=3,color=10,thick=2.0
      plots,i*dt/8.64e4,n.ns2p,psym=3,color=75,thick=2.0
      plots,i*dt/8.64e4,n.ns3p,psym=3,color=50,thick=2.0
      plots,i*dt/8.64e4,n.nel,psym=3,color=100,thick=2.0
   ;   plots,i*dt/8.64e4,n.no,psym=3,color=75,thick=2.0
      plots,i*dt/8.64e4,n.nop,psym=3,color=150,thick=2.0
      plots,i*dt/8.64e4,n.no2p,psym=3,color=200,thick=2.0
   endif

   if (plot_temperature eq 1) then begin
      plots,i*dt/8.64e4,T.Tsp,psym=3,color=10,thick=2.0
      plots,i*dt/8.64e4,T.Ts2p,psym=3,color=75,thick=2.0
      plots,i*dt/8.64e4,T.Ts3p,psym=3,color=50,thick=2.0
      plots,i*dt/8.64e4,T.Tel,psym=3,color=100,thick=2.0
      plots,i*dt/8.64e4,T.Top,psym=3,color=150,thick=2.0
      plots,i*dt/8.64e4,T.To2p,psym=3,color=200,thick=2.0
   endif

   p(i).Puv = r.Puv
   p(i).psp = r.psp
   p(i).ps2p = r.ps2p
   p(i).ps3p = r.ps3p
   p(i).pop = r.pop
   p(i).po2p = r.po2p

endfor

save,filename='src_func.sav',src
save,filename='lss_func.sav',lss
save,filename='temps.sav',temps
save,filename='dens.sav',dens
save,filename='pwr.sav',p
print,'pwr.sav...saved'

if (lonoff eq 1.0) then begin
legend,['S!u+!n','S!u++!n','S!u+++!n','O!u+!n','O!u++!n','e'],$
       linestyle=[0,0,0,0,0,0],colors=[10,75,50,150,200,100],/bottom,$
       /right
endif

upsp = n.nsp
ups2p = n.ns2p
ups3p = n.ns3p
upop = n.nop
upo2p = n.no2p

print,' '
print,'What's up...'
print,'nsp....',n.nsp
print,'ns2p...',n.ns2p
print,'ns3p...',n.ns3p
print,'nop....',n.nop
print,'no2p...',n.no2p
print,' '
print,'S/O:     ',n.ns/n.no
print,'S+/O+:   ',n.nsp/n.nop
print,'S++/S+:  ',n.ns2p/n.nsp
print,'S+++/S+: ',n.ns3p/n.nsp
print,' '
print,'Tsp......',T.Tsp
print,'Ts2p.....',T.Ts2p
print,'Ts3p.....',T.Ts3p
print,'Top......',T.Top
print,'To2p.....',T.To2p
print,'Tel......',T.Tel
print,' '
save,filename='restart_cm3.sav',n,n1,np,T,T1,Tp,nT,nT1,nTp

;cm3_ekappa,n,T

return
end
;-----------------------------------------------------------------------




