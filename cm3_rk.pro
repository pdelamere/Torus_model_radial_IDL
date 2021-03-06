;-----------------------------------------------------------------------
pro get_rates,r,T
;-----------------------------------------------------------------------
@common

;Electron impact ionization (Schreier and Eviatar 1998)
r.is = 7.4e-8*sqrt(T.Tel_el)*exp(-10.36/T.Tel_el)     ;check this!!! 
r.isp = 1.274e-8*sqrt(T.Tel_sp)*exp(-23.1/T.Tel_sp)  
r.is2p = 5.8e-9*sqrt(T.Tel_s2p)*exp(-34.88/T.Tel_s2p)  

;r.is = cfit(16,16,T.Tel_el,c)
;r.isp = cfit(16,15,T.Tel_sp,c)/
;r.is2p = cfit(16,14,T.Tel_s2p,c)

;r.ish = 1.464e-6*alog(T.Telh/27.37)/sqrt(T.Telh)
;r.isph = 4.44e-7*alog(T.Telh/27.07)/sqrt(T.Telh)
;r.is2ph = 1.76e-7*alog(T.Telh/27.07)/sqrt(T.Telh)

;r.ish = cfit(16,16,T.Telh,c)
;r.isph = cfit(16,15,T.Telh,c)
;r.is2ph = cfit(16,14,T.Telh,c)

r.io = 1.015e-8*sqrt(T.Tel_el)*exp(-14.54/T.Tel_el)
r.iop = 3.78e-9*sqrt(T.Tel_op)*exp(-33.42/T.Tel_op)

;r.io = cfit(8,8,T.Tel_el,c)
;r.iop = cfit(8,7,T.Tel_op,c)

;r.ioh = 5.74e-7*alog(T.Telh/27.07)/sqrt(T.Telh)
;r.ioph = 1.76e-7*alog(T.Telh/27.07)/sqrt(T.Telh)

;r.ioh = cfit(8,8,T.Telh,c)
;r.ioph = cfit(8,7,T.Telh,c)

;dielectron recombination (Smith and Strobel 1985)
if (T.Tel_s3p le 6.9) then r.rs3p = 5.4e-11*T.Tel_s3p^0.53
if (T.Tel_s3p gt 6.9) then r.rs3p = 4.3e-12*T.Tel_s3p^1.84  

r.rs2p = r.rs3p/3.0   ;massive fudge based on Johnson and Strobel '82

r.ro2p = r.rs3p/7.0   ;another massive fudge base on Johnson and Stobel

;charge exchange

r.cx_k0 = k0*8.1e-9       ;S+ + S++ -> S++ + S+
r.cx_k1 = k1*2.4e-8
r.cx_k2 = k2*3.0e-10
r.cx_k3 = k3*7.8e-9
r.cx_k4 = k4*1.32e-8
r.cx_k5 = k5*1.32e-8
r.cx_k6 = k6*5.2e-10
r.cx_k7 = k7*5.4e-9
r.cx_k8 = k8*9.0e-11   ;Schreier
;r.cx_k8 = k8*6.0e-11    ;McGrath
r.cx_k9 = k9*3.0e-9
r.cx_k10 = k10*2.34e-8
r.cx_k11 = k11*1.62e-8
r.cx_k12 = k12*2.3e-9  ;new
;r.cx_k12 = k12*7.3e-9   ;old
r.cx_k13 = k13*1.4e-9   ;Schreier for L=6
r.cx_k14 = k14*1.92e-8
;r.cx_k15 = k15*9e-9    ;Schreier for L=6
r.cx_k15 = k15*9e-10    ;McGrath for L=6

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
S = r.S_production ;+ r.cx_k8*n.no*n.nsp
L = r.is*n.ns*n.fc*n.nel + r.ish*n.ns*n.fh*n.nel + $ 
    r.cx_k4*n.ns*n.ns3p + r.cx_k11*n.ns*n.no2p + $
    r.cx_k10*n.ns*n.no2p + $
    r.cx_k2*n.ns*n.ns2p + r.cx_k9*n.ns*n.nop

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

;print,' '
;print,'S+ sources...'
;print,'ionization cold...',r.is*n.ns*n.fc*n.nel
;print,'ionixation hot....',r.ish*n.ns*n.fh*n.nel
;print,'recomination S++..',r.rs2p*n.ns2p*n.nel
;print,'k4................',r.cx_k4*n.ns*n.ns3p
;print,'k10...............',r.cx_k10*n.ns*n.no2p
;print,'k2................',2.0*r.cx_k2*n.ns*n.ns2p
;print,'k9................',r.cx_k9*n.ns*n.nop
;print,'k12...............',r.cx_k12*n.no*n.ns2p 
;print,' '
;print,'S+ losses...'
;print,'ionization s+ cold...',r.isp*n.nsp*n.fc*n.nel
;print,'k13..................',r.cx_k13*n.no2p*n.nsp
;print,'transport............',n.nsp*r.Transport
;print,'ioniation s+ hot.....',r.isph*n.nsp*n.fh*n.nel
;print,'k8...................',r.cx_k8*n.no*n.nsp 
;print,' '
;print,'S-L...',S-L

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
S = r.O_production ;+ r.cx_k9*n.ns*n.nop
L = r.io*n.no*n.fc*n.nel + r.ioh*n.no*n.fh*n.nel + r.cx_k14*n.no*n.ns3p + $
    r.cx_k6*n.no*n.no2p + r.cx_k8*n.no*n.nsp + r.cx_k12*n.no*n.ns2p

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

if (T.Tel le r.emistemp(1)) then T.Tel = r.emistemp(1)
if (T.Tel_sp le r.emistemp(1)) then T.Tel_sp = r.emistemp(1)
if (T.Tel_s2p le r.emistemp(1)) then T.Tel_s2p = r.emistemp(1)
if (T.Tel_s3p le r.emistemp(1)) then T.Tel_s3p = r.emistemp(1)
if (T.Tel_op le r.emistemp(1)) then T.Tel_op = r.emistemp(1)
if (T.Tel_o2p le r.emistemp(1)) then T.Tel_o2p = r.emistemp(1)

;print,T.Tel,T.Tel_sp,T.Tel_s2p,T.Tel_s3p,T.Tel_op,T.Tel_o2p

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

Teq = nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp)*n.nsp*(T.Tsp-T.Tel_sp) + $
      nu_ie(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p)*n.ns2p*(T.Ts2p-T.Tel_s2p) + $
      nu_ie(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p)*n.ns3p*(T.Ts3p-T.Tel_s3p) + $
      nu_ie(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top)*n.nop*(T.Top-T.Tel_op) + $
      nu_ie(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p)*n.no2p*(T.To2p-T.Tel_o2p) + $
      nu_ee(n.fc*n.nel,n.fh*n.nel,T.Tel_el,T.Telh)*(T.Telh - T.Tel_el)

EF_el = Teq - L/1.0

r.Puv = r.psp + r.ps2p + r.ps3p + r.pop + r.po2p

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

r.Sei = (r.is*n.ns*n.fc*n.nel + r.ish*n.ns*n.fh*n.nel)*T.Tpu_s
r.Secx = (r.cx_k4*n.ns*n.ns3p + r.cx_k10*n.ns*n.no2p + $
         r.cx_k1*n.ns*n.nsp + r.cx_k2*n.ns*n.ns2p + $
         r.cx_k9*n.ns*n.nop+r.cx_k3*n.ns*n.ns2p + $
         r.cx_k11*n.ns*n.no2p)*T.Tpu_s

L= r.cx_k1*n.ns*n.nsp*T.Tsp + $     ;S + S+ -> S+ + S
   r.cx_k8*n.no*n.nsp*T.Tsp + $     ;O + S+ -> O+ + S
   r.cx_k13*n.no2p*n.nsp*T.Tsp + $  ;O++ + S+ -> O+ + S++
   r.cx_k0*n.nsp*n.ns2p*T.Tsp + $   ;S+ + S++ -> S++ + S+
   r.Transport*n.nsp*T.Tsp    

Teq= nu_ii(32.0,32.0,1.0,2.0,n.nsp,n.ns2p,T.Tsp,T.Ts2p)*n.nsp*(T.Ts2p-T.Tsp)+$
   nu_ii(32.0,32.0,1.0,3.0,n.nsp,n.ns3p,T.Tsp,T.Ts3p)*n.nsp*(T.Ts3p-T.Tsp)+$
   nu_ii(32.0,16.0,1.0,1.0,n.nsp,n.nop,T.Tsp,T.Top)*n.nsp*(T.Top-T.Tsp)+$
   nu_ii(32.0,16.0,1.0,2.0,n.nsp,n.no2p,T.Tsp,T.To2p)*n.nsp*(T.To2p-T.Tsp)+$
   nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp)*n.nsp*(T.Tel_sp-T.Tsp) + $
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
   nu_ie(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p)*n.ns2p*(T.Tel_s2p-T.Ts2p) + $
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
   nu_ie(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p)*n.ns3p*(T.Tel_s3p-T.Ts3p) + $
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

r.Oei = (r.io*n.no*n.fc*n.nel + r.ioh*n.no*n.fh*n.nel)*T.Tpu_o
r.Oecx = (r.cx_k5*n.no*n.nop + r.cx_k6*n.no*n.no2p + $
         r.cx_k8*n.no*n.nsp + r.cx_k12*n.no*n.ns2p + $
         r.cx_k14*n.no*n.ns3p+r.cx_k7*n.no*n.no2p)*T.Tpu_o

L= r.cx_k5*n.no*n.nop*T.Top + $     ;O + O+ -> O+ + O
   r.cx_k9*n.ns*n.nop*T.Top + $     ;S + O+ -> S+ + O
   r.Transport*n.nop*T.Top    

Teq= nu_ii(16.0,16.0,1.0,2.0,n.nop,n.no2p,T.Top,T.To2p)*n.nop*(T.To2p-T.Top)+$
   nu_ii(16.0,32.0,1.0,1.0,n.nop,n.nsp,T.Top,T.Tsp)*n.nop*(T.Tsp-T.Top)+$
   nu_ii(16.0,32.0,1.0,2.0,n.nop,n.ns2p,T.Top,T.Ts2p)*n.nop*(T.Ts2p-T.Top)+$
   nu_ii(16.0,32.0,1.0,3.0,n.nop,n.ns3p,T.Top,T.Ts3p)*n.nop*(T.Ts3p-T.Top)+$
   nu_ie(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top)*n.nop*(T.Tel_op-T.Top) + $
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
   nu_ie(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p)*n.no2p*(T.Tel_o2p-T.To2p) + $
   nu_ie(16.0,2.0,n.fh*n.nel,n.no2p,T.Telh,T.To2p)*n.no2p*(T.Telh-T.To2p) 
 
EF_o2p = S - L + Teq

return,EF_o2p
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------
pro cm3_model,n,r,T,src,lss,temps,dens
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

kk1 = {rk1, nel:0.0,ns:0.0,nsp:0.0,ns2p:0.0,ns3p:0.0,no:0.0,nop:0.0,no2p:0.0}
kk2 = {rk2, nel:0.0,ns:0.0,nsp:0.0,ns2p:0.0,ns3p:0.0,no:0.0,nop:0.0,no2p:0.0}
kk3 = {rk3, nel:0.0,ns:0.0,nsp:0.0,ns2p:0.0,ns3p:0.0,no:0.0,nop:0.0,no2p:0.0}
kk4 = {rk4, nel:0.0,ns:0.0,nsp:0.0,ns2p:0.0,ns3p:0.0,no:0.0,nop:0.0,no2p:0.0}

Tkk1 = {Trk1, el:0.0,s:0.0,sp:0.0,s2p:0.0,s3p:0.0,o:0.0,op:0.0,o2p:0.0}
Tkk2 = {Trk2, el:0.0,s:0.0,sp:0.0,s2p:0.0,s3p:0.0,o:0.0,op:0.0,o2p:0.0}
Tkk3 = {Trk3, el:0.0,s:0.0,sp:0.0,s2p:0.0,s3p:0.0,o:0.0,op:0.0,o2p:0.0}
Tkk4 = {Trk4, el:0.0,s:0.0,sp:0.0,s2p:0.0,s3p:0.0,o:0.0,op:0.0,o2p:0.0}

T = {temp, Tel: Te0, Telh: Teh0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, Top: Ti0, $
           To2p: Ti0, Tsp_pu: 0.0, Ts2p_pu: 0.0, Ts3p_pu: 0.0, $
           Top_pu: 0.0, To2p_pu: 0.0, Tpu_s: 0.0, Tpu_o: 0.0, $
           Tel_sp: Te0, Tel_s2p: Te0, Tel_s3p: Te0, Tel_op: Te0, $
           Tel_o2p: Te0, Tel_el: Te0}
nT = {energy, el: nel0*Te0, sp: nsp0*Ti0, s2p: ns2p0*Ti0, $
              s3p: ns3p0*Ti0, op: nop0*Ti0, o2p: no2p0*Ti0, $
              sp_pu: 0.0, s2p_pu: 0.0, s3p_pu: 0.0, op_pu: 0.0,$
              o2p_pu: 0.0} 

;advance density for improved Euler method, half time step
n1 = {density_1, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
                 no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh, $
              nsp_pu: 0.0, ns2p_pu: 0.0, ns3p_pu: 0.0, nop_pu: 0.0, $
              no2p_pu: 0.0}
n2 = {density_2, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
                 no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh, $
              nsp_pu: 0.0, ns2p_pu: 0.0, ns3p_pu: 0.0, nop_pu: 0.0, $
              no2p_pu: 0.0}
n3 = {density_3, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
                 no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh, $
              nsp_pu: 0.0, ns2p_pu: 0.0, ns3p_pu: 0.0, nop_pu: 0.0, $
              no2p_pu: 0.0}
T1 = {temp_1, Tel: Te0, Telh: Teh0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, Top: Ti0, $
           To2p: Ti0, Tsp_pu: 0.0, Ts2p_pu: 0.0, Ts3p_pu: 0.0, $
           Top_pu: 0.0, To2p_pu: 0.0, Tpu_s: 0.0, Tpu_o: 0.0, $
           Tel_sp: Te0, Tel_s2p: Te0, Tel_s3p: Te0, Tel_op: Te0, $
           Tel_o2p: Te0, Tel_el: Te0}
T2 = {temp_2, Tel: Te0, Telh: Teh0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, Top: Ti0, $
           To2p: Ti0, Tsp_pu: 0.0, Ts2p_pu: 0.0, Ts3p_pu: 0.0, $
           Top_pu: 0.0, To2p_pu: 0.0, Tpu_s: 0.0, Tpu_o: 0.0, $
           Tel_sp: Te0, Tel_s2p: Te0, Tel_s3p: Te0, Tel_op: Te0, $
           Tel_o2p: Te0, Tel_el: Te0}
T3 = {temp_3, Tel: Te0, Telh: Teh0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, Top: Ti0, $
           To2p: Ti0, Tsp_pu: 0.0, Ts2p_pu: 0.0, Ts3p_pu: 0.0, $
           Top_pu: 0.0, To2p_pu: 0.0, Tpu_s: 0.0, Tpu_o: 0.0, $
           Tel_sp: Te0, Tel_s2p: Te0, Tel_s3p: Te0, Tel_op: Te0, $
           Tel_o2p: Te0, Tel_el: Te0}
nT1 = {energy_1, el: nel0*Te0, sp: nsp0*Ti0, s2p: ns2p0*Ti0, $
              s3p: ns3p0*Ti0, op: nop0*Ti0, o2p: no2p0*Ti0, $
              sp_pu: 0.0, s2p_pu: 0.0, s3p_pu: 0.0, op_pu: 0.0,$
              o2p_pu: 0.0} 
nT2 = {energy_2, el: nel0*Te0, sp: nsp0*Ti0, s2p: ns2p0*Ti0, $
              s3p: ns3p0*Ti0, op: nop0*Ti0, o2p: no2p0*Ti0, $
              sp_pu: 0.0, s2p_pu: 0.0, s3p_pu: 0.0, op_pu: 0.0,$
              o2p_pu: 0.0} 
nT3 = {energy_3, el: nel0*Te0, sp: nsp0*Ti0, s2p: ns2p0*Ti0, $
              s3p: ns3p0*Ti0, op: nop0*Ti0, o2p: no2p0*Ti0, $
              sp_pu: 0.0, s2p_pu: 0.0, s3p_pu: 0.0, op_pu: 0.0,$
              o2p_pu: 0.0} 

;density at full time step advance
np = {density_p, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
                 no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh, $
              nsp_pu: 0.0, ns2p_pu: 0.0, ns3p_pu: 0.0, nop_pu: 0.0, $
              no2p_pu: 0.0}
nave = {density_ave, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
                 no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh, $
              nsp_pu: 0.0, ns2p_pu: 0.0, ns3p_pu: 0.0, nop_pu: 0.0, $
              no2p_pu: 0.0}
Tp = {temp_p, Tel: Te0, Telh: Teh0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, Top: Ti0, $
           To2p: Ti0, Tsp_pu: 0.0, Ts2p_pu: 0.0, Ts3p_pu: 0.0, $
           Top_pu: 0.0, To2p_pu: 0.0, Tpu_s: 0.0, Tpu_o: 0.0, $
           Tel_sp: Te0, Tel_s2p: Te0, Tel_s3p: Te0, Tel_op: Te0, $
           Tel_o2p: Te0, Tel_el: Te0}
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
            pop: 0.0, po2p: 0.0, Sei: 0.0, Secx: 0.0, Oei: 0.0, Oecx:0.0}

src = replicate({source_func, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0,$
              o: 0.0, op: 0.0, o2p: 0.0,tm:0.0},ntm+1)
lss = replicate({loss_func, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0,$
              o: 0.0, op: 0.0, o2p: 0.0,tm:0.0},ntm+1)

temps = replicate({temperatures, tsp: 0.0, Ts2p: 0.0, Ts3p: 0.0, Top: 0.0, $
                                 To2p: 0.0, Telec: 0.0},ntm+1)

dens = replicate({dens, nsp: 0.0, ns2p: 0.0, ns3p: 0.0, nop: 0.0, $
                  no2p: 0.0, nel: 0.0},ntm+1)

p = replicate({Power_uv, Puv: 0.0, psp: 0.0, ps2p: 0.0, ps3p: 0.0, $
                         pop: 0.0, po2p: 0.0},ntm+1)

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

r.ish = cfit(16,16,T.Telh,c)
r.isph = cfit(16,15,T.Telh,c)
r.is2ph = cfit(16,14,T.Telh,c)
r.ioh = cfit(8,8,T.Telh,c)
r.ioph = cfit(8,7,T.Telh,c)

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

one_six = 1./6.
for i = 0,ntm do begin

   get_rates,r,T

   kk1.ns = dt*F_s(n,r,S,L)
   src(i).s = S
   lss(i).s = L
   kk1.nsp = dt*F_sp(n,r,S,L)
   src(i).sp = S
   lss(i).sp = L
   kk1.ns2p = dt*F_s2p(n,r,S,L)
   src(i).s2p = S
   lss(i).s2p = L
   kk1.ns3p =  dt*F_s3p(n,r,S,L)
   src(i).s3p = S
   lss(i).s3p = L
   kk1.no =  dt*F_o(n,r,S,L)
   src(i).o = S
   lss(i).o = L
   kk1.nop = dt*F_op(n,r,S,L)
   src(i).op = S
   lss(i).op = L 
   kk1.no2p = dt*F_o2p(n,r,S,L)
   src(i).o2p = S
   lss(i).o2p = L

   src(i).tm = dt*i/8.64e4
   lss(i).tm = dt*i/8.64e4

   n1.ns = n.ns + 0.5*kk1.ns
   n1.nsp = n.nsp + 0.5*kk1.nsp
   n1.ns2p = n.ns2p + 0.5*kk1.ns2p
   n1.ns3p = n.ns3p + 0.5*kk1.ns3p
   n1.no = n.no + 0.5*kk1.no
   n1.nop = n.nop + 0.5*kk1.nop
   n1.no2p = n.no2p + 0.5*kk1.no2p
   n1.nel = 1.0*n1.nsp+2.0*n1.ns2p+3.0*n1.ns3p+1.0*n1.nop+2.0*n1.no2p

   kk2.ns = dt*F_s(n1,r,S,L)
   kk2.nsp = dt*F_sp(n1,r,S,L)
   kk2.ns2p = dt*F_s2p(n1,r,S,L)
   kk2.ns3p = dt*F_s3p(n1,r,S,L)
   kk2.no = dt*F_o(n1,r,S,L)
   kk2.nop = dt*F_op(n1,r,S,L)
   kk2.no2p = dt*F_o2p(n1,r,S,L)

   n2.ns = n.ns + 0.5*kk2.ns
   n2.nsp = n.nsp + 0.5*kk2.nsp
   n2.ns2p = n.ns2p + 0.5*kk2.ns2p
   n2.ns3p = n.ns3p + 0.5*kk2.ns3p
   n2.no = n.no + 0.5*kk2.no
   n2.nop = n.nop + 0.5*kk2.nop
   n2.no2p = n.no2p + 0.5*kk2.no2p
   n2.nel = 1.0*n2.nsp+2.0*n2.ns2p+3.0*n2.ns3p+1.0*n2.nop+2.0*n2.no2p

   kk3.ns = dt*F_s(n2,r,S,L)
   kk3.nsp = dt*F_sp(n2,r,S,L)
   kk3.ns2p = dt*F_s2p(n2,r,S,L)
   kk3.ns3p = dt*F_s3p(n2,r,S,L)
   kk3.no = dt*F_o(n2,r,S,L)
   kk3.nop = dt*F_op(n2,r,S,L)
   kk3.no2p = dt*F_o2p(n2,r,S,L)
   
   n3.ns = n.ns + kk3.ns
   n3.nsp = n.nsp + kk3.nsp
   n3.ns2p = n.ns2p + kk3.ns2p
   n3.ns3p = n.ns3p + kk3.ns3p
   n3.no = n.no + kk3.no
   n3.nop = n.nop + kk3.nop
   n3.no2p = n.no2p + kk3.no2p
   n3.nel = 1.0*n3.nsp+2.0*n3.ns2p+3.0*n3.ns3p+1.0*n3.nop+2.0*n3.no2p

   kk4.ns = dt*F_s(n3,r,S,L)
   kk4.nsp = dt*F_sp(n3,r,S,L)
   kk4.ns2p = dt*F_s2p(n3,r,S,L)
   kk4.ns3p = dt*F_s3p(n3,r,S,L)
   kk4.no = dt*F_o(n3,r,S,L)
   kk4.nop = dt*F_op(n3,r,S,L)
   kk4.no2p = dt*F_o2p(n3,r,S,L)

   np.ns = n.ns + one_six*(kk1.ns + 2.0*kk2.ns + 2.0*kk3.ns + kk4.ns)
   np.nsp = n.nsp + one_six*(kk1.nsp + 2.0*kk2.nsp + 2.0*kk3.nsp + kk4.nsp)
   np.ns2p = n.ns2p + one_six*(kk1.ns2p + 2.0*kk2.ns2p + 2.0*kk3.ns2p + kk4.ns2p)
   np.ns3p = n.ns3p + one_six*(kk1.ns3p + 2.0*kk2.ns3p + 2.0*kk3.ns3p + kk4.ns3p)
   np.no = n.no + one_six*(kk1.no + 2.0*kk2.no + 2.0*kk3.no + kk4.no)
   np.nop = n.nop + one_six*(kk1.nop + 2.0*kk2.nop + 2.0*kk3.nop + kk4.nop)
   np.no2p = n.no2p + one_six*(kk1.no2p + 2.0*kk2.no2p + 2.0*kk3.no2p + kk4.no2p)
   np.nel = 1.0*np.nsp+2.0*np.ns2p+3.0*np.ns3p+1.0*np.nop+2.0*np.no2p

;   np.ns = n.ns + kk2.ns 
;   np.nsp = n.nsp + kk2.nsp
;   np.ns2p = n.ns2p + kk2.ns2p
;   np.ns3p = n.ns3p + kk2.ns3p
;   np.no = n.no + kk2.no
;   np.nop = n.nop + kk2.nop
;   np.no2p = n.no2p + kk2.no2p
;   np.nel = 1.0*np.nsp+2.0*np.ns2p+3.0*np.ns3p+1.0*np.nop+2.0*np.no2p

   nave.ns = 0.5*(n.ns+np.ns)
   nave.nsp = 0.5*(n.nsp+np.nsp)
   nave.ns2p = 0.5*(n.ns2p+np.ns2p)
   nave.ns3p = 0.5*(n.ns3p+np.ns3p)
   nave.no = 0.5*(n.no+np.no)
   nave.nop = 0.5*(n.nop+np.nop)
   nave.no2p = 0.5*(n.no2p+np.no2p)
   nave.nel = 0.5*(n.nel+np.nel)

   Tkk1.sp = dt*EF_sp(n,T,nT,r)
   Tkk1.s2p = dt*EF_s2p(n,T,nT,r)
   Tkk1.s3p = dt*EF_s3p(n,T,nT,r)
   Tkk1.op = dt*EF_op(n,T,nT,r)
   Tkk1.o2p = dt*EF_o2p(n,T,nT,r)
   Tkk1.el = dt*EF_el(n,T,nT,r)

   nT1.sp = nT.sp + 0.5*Tkk1.sp
   nT1.s2p = nT.s2p + 0.5*Tkk1.s2p
   nT1.s3p = nT.s3p + 0.5*Tkk1.s3p
   nT1.op = nT.op + 0.5*Tkk1.op
   nT1.o2p = nT.o2p + 0.5*Tkk1.o2p
   nT1.el = nT.el + 0.5*Tkk1.el
   update_temp,n1,nT1,T1

   Tkk2.sp = dt*EF_sp(n1,T1,nT1,r)
   Tkk2.s2p = dt*EF_s2p(n1,T1,nT1,r)
   Tkk2.s3p = dt*EF_s3p(n1,T1,nT1,r)
   Tkk2.op = dt*EF_op(n1,T1,nT1,r)
   Tkk2.o2p = dt*EF_o2p(n1,T1,nT1,r)
   Tkk2.el = dt*EF_el(n1,T1,nT1,r)

   nT2.sp = nT.sp + 0.5*Tkk2.sp
   nT2.s2p = nT.s2p + 0.5*Tkk2.s2p
   nT2.s3p = nT.s3p + 0.5*Tkk2.s3p
   nT2.op = nT.op + 0.5*Tkk2.op
   nT2.o2p = nT.o2p + 0.5*Tkk2.o2p
   nT2.el = nT.el + 0.5*Tkk2.el
   update_temp,n2,nT2,T2

   Tkk3.sp = dt*EF_sp(n2,T2,nT2,r)
   Tkk3.s2p = dt*EF_s2p(n2,T2,nT2,r)
   Tkk3.s3p = dt*EF_s3p(n2,T2,nT2,r)
   Tkk3.op = dt*EF_op(n2,T2,nT2,r)
   Tkk3.o2p = dt*EF_o2p(n2,T2,nT2,r)
   Tkk3.el = dt*EF_el(n2,T2,nT2,r)

   nT3.sp = nT.sp + Tkk3.sp
   nT3.s2p = nT.s2p + Tkk3.s2p
   nT3.s3p = nT.s3p + Tkk3.s3p
   nT3.op = nT.op + Tkk3.op
   nT3.o2p = nT.o2p + Tkk3.o2p
   nT3.el = nT.el + Tkk3.el
   update_temp,n3,nT3,T3

   Tkk4.sp = dt*EF_sp(n3,T3,nT3,r)
   Tkk4.s2p = dt*EF_s2p(n3,T3,nT3,r)
   Tkk4.s3p = dt*EF_s3p(n3,T3,nT3,r)
   Tkk4.op = dt*EF_op(n3,T3,nT3,r)
   Tkk4.o2p = dt*EF_o2p(n3,T3,nT3,r)
   Tkk4.el = dt*EF_el(n3,T3,nT3,r)

;   nTp.sp = nT.sp + one_six*(Tkk1.sp + 2.0*Tkk2.sp + 2.0*Tkk3.sp + Tkk4.sp)
;   nTp.s2p = nT.s2p + one_six*(Tkk1.s2p + 2.0*Tkk2.s2p + 2.0*Tkk3.s2p + Tkk4.s2p)
;   nTp.s3p = nT.s3p + one_six*(Tkk1.s3p + 2.0*Tkk2.s3p + 2.0*Tkk3.s3p + Tkk4.s3p)
;   nTp.op = nT.op + one_six*(Tkk1.op + 2.0*Tkk2.op + 2.0*Tkk3.op + Tkk4.op)
;   nTp.o2p = nT.o2p + one_six*(Tkk1.o2p + 2.0*Tkk2.o2p + 2.0*Tkk3.o2p + Tkk4.o2p)
;   nTp.el = nT.el + one_six*(Tkk1.el + 2.0*Tkk2.el + 2.0*Tkk3.el + Tkk4.el)

   nTp.sp = nT.sp + Tkk2.sp
   nTp.s2p = nT.s2p + Tkk2.s2p
   nTp.s3p = nT.s3p + Tkk2.s3p 
   nTp.op = nT.op + Tkk2.op
   nTp.o2p = nT.o2p + Tkk2.o2p
   nTp.el = nT.el + Tkk2.el


;   ;full time step advance for energy
;   nT1.sp = nT.sp + dt*EF_sp(n,T,nT,r)
;   nT1.s2p = nT.s2p + dt*EF_s2p(n,T,nT,r)
;   nT1.s3p = nT.s3p + dt*EF_s3p(n,T,nT,r)
;   nT1.op = nT.op + dt*EF_op(n,T,nT,r)
;   nT1.o2p = nT.o2p + dt*EF_o2p(n,T,nT,r)
;   nT1.el = nT.el + dt*EF_el(n,T,nT,r)
;   update_temp,n1,nT1,T1
;   if (ft_ave eq 1) then cm3_ekappa,n1,T1
;   if (ft_ave eq 0) then begin
;      T1.Tel_sp = T1.Tel
;      T1.Tel_s2p = T1.Tel
;      T1.Tel_s3p = T1.Tel
;      T1.Tel_op = T1.Tel
;      T1.Tel_o2p = T1.Tel
;      T1.Tel_el = T1.Tel
;   endif

;   nTp.sp = nT.sp + dt*0.5*(EF_sp(n1,T1,nT1,r)+ EF_sp(n,T,nT,r))
;   nTp.s2p = nT.s2p + dt*0.5*(EF_s2p(n1,T1,nT1,r)+ EF_s2p(n,T,nT,r))
;   nTp.s3p = nT.s3p + dt*0.5*(EF_s3p(n1,T1,nT1,r)+ EF_s3p(n,T,nT,r))
;   nTp.op = nT.op + dt*0.5*(EF_op(n1,T1,nT1,r)+ EF_op(n,T,nT,r))
;   nTp.o2p = nT.o2p + dt*0.5*(EF_o2p(n1,T1,nT1,r)+ EF_o2p(n,T,nT,r))
;   nTp.el = nT.el + dt*0.5*(EF_el(n1,T1,nT1,r)+ EF_el(n,T,nT,r))
   update_temp,np,nTp,Tp
   if (ft_ave eq 1) then cm3_ekappa,np,Tp
   if (ft_ave eq 0) then begin
      Tp.Tel_sp = Tp.Tel
      Tp.Tel_s2p = Tp.Tel
      Tp.Tel_s3p = Tp.Tel
      Tp.Tel_op = Tp.Tel
      Tp.Tel_o2p = Tp.Tel
      Tp.Tel_el = Tp.Tel
   endif

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
   p(i).Puv = r.Puv
   p(i).psp = r.psp
   p(i).ps2p = r.ps2p
   p(i).ps3p = r.ps3p
   p(i).pop = r.pop
   p(i).po2p = r.po2p

;   print,T.Tsp,T.Ts2p,T.Ts3p,T.Top,T.To2p,T.Tel

;update Pickup Densities
  
;   n.nsp_pu = nT.sp_pu/T.Tsp_pu
;   n.ns2p_pu = nT.s2p_pu/T.Ts2p_pu
;   n.ns3p_pu = nT.s3p_pu/T.Ts3p_pu
;   n.nop_pu = nT.op_pu/T.Top_pu
;   n.no2p_pu = nT.o2p_pu/T.To2p_pu

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

   if (ft_ave eq 1) then cm3_ekappa,n,T
   if (ft_ave eq 0) then begin
      T.Tel_sp = T.Tel
      T.Tel_s2p = T.Tel
      T.Tel_s3p = T.Tel
      T.Tel_op = T.Tel
      T.Tel_o2p = T.Tel
      T.Tel_el = T.Tel
   endif

endfor

save,filename='src_func.sav',src
save,filename='lss_func.sav',lss
save,filename='temps.sav',temps
save,filename='dens.sav',dens
save,filename='pwr.sav',p

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
print,'nsp....',n.nsp/n.nel, 0.12, 0.092
print,'ns2p...',n.ns2p/n.nel, 0.15, 0.17
print,'ns3p...',n.ns3p/n.nel, 0, 0.023
print,'nop....',n.nop/n.nel, 0.36, 0.47
print,'no2p...',n.no2p/n.nel, 0.12, 0.016
print,' '
print,'S/O:     ',n.ns/n.no
print,'S+/O+:   ',n.nsp/n.nop
print,'S++/S+:  ',n.ns2p/n.nsp
print,'S+++/S+: ',n.ns3p/n.nsp
print,'O++/O+:  ',n.no2p/n.nop
print,' '
print,'Tsp......',T.Tsp
print,'Ts2p.....',T.Ts2p
print,'Ts3p.....',T.Ts3p
print,'Top......',T.Top
print,'To2p.....',T.To2p
print,'Tel......',T.Tel
print,' '
print,'Fluxtube averaged electron temp....................',T.Tel_el
print,'Fluxtube averaged electron temp, (O+ weighted).....',T.Tel_op
print,'Fluxtube averaged electron temp, (O++ weighted)....',T.Tel_o2p
print,'Fluxtube averaged electron temp, (S+ weighted).....',T.Tel_sp
print,'Fluxtube averaged electron temp, (S++ weighted)....',T.Tel_s2p
print,'Fluxtube averaged electron temp, (S+++ weighted)...',T.Tel_s3p
print,' '
save,filename='restart_cm3.sav',n,n1,np,T,T1,Tp,nT,nT1,nTp

;cm3_ekappa,n,T

return
end
;-----------------------------------------------------------------------




