; Version 1.3, Based largely on cm3_v1.3, constant hot electron
; density based on feh at equator.


;-----------------------------------------------------------------------
function ftavg,n1,n2,h1,h2,h3
;ions with ions
;-----------------------------------------------------------------------
@common

a = (h1^2 + h2^2)/(h1^2 * h2^2)

dN3 = n1*n2*sqrt(!pi/a)

dn3 = dN3/(sqrt(!pi)*h3)


return,dn3
end
;-----------------------------------------------------------------------

;-----------------------------------------------------------------------
function ftavg_el,n1,n,h1,h,h3
;ions with electrons, sum over each ion species separately and sum for
;total electron contribution
;-----------------------------------------------------------------------
@common

dN3 = 0.0

a = (h1^2 + h.sp^2)/(h1^2 * h.sp^2)
dN3 = dN3 + n1*(1.0*n.nsp)*sqrt(!pi/a)

a = (h1^2 + h.s2p^2)/(h1^2 * h.s2p^2)
dN3 = dN3 + n1*(2.0*n.ns2p)*sqrt(!pi/a)

a = (h1^2 + h.s3p^2)/(h1^2 * h.s3p^2)
dN3 = dN3 + n1*(3.0*n.ns3p)*sqrt(!pi/a)

a = (h1^2 + h.op^2)/(h1^2 * h.op^2)
dN3 = dN3 + n1*(1.0*n.nop)*sqrt(!pi/a)

a = (h1^2 + h.o2p^2)/(h1^2 * h.o2p^2)
dN3 = dN3 + n1*(2.0*n.no2p)*sqrt(!pi/a)

dn3 = dN3/(sqrt(!pi)*h3)


return,dn3
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function ftavg_elh,n1,n,h1,h,h3
;ions with electrons, sum over each ion species separately and sum for
;total electron contribution
;-----------------------------------------------------------------------
@common

dN3 = n1*n.fh*n.nel*sqrt(!pi)*h1
;print,'dN3...',dN3,n1,n.fh*n.nel

;dN3 = 0.0

;a = (h1^2 + h.sp^2)/(h1^2 * h.sp^2)
;dN3 = dN3 + n1*(1.0*n.fh*n.nsp)*sqrt(!pi/a)

;a = (h1^2 + h.s2p^2)/(h1^2 * h.s2p^2)
;dN3 = dN3 + n1*(2.0*n.fh*n.ns2p)*sqrt(!pi/a)

;a = (h1^2 + h.s3p^2)/(h1^2 * h.s3p^2)
;dN3 = dN3 + n1*(3.0*n.fh*n.ns3p)*sqrt(!pi/a)

;a = (h1^2 + h.op^2)/(h1^2 * h.op^2)
;dN3 = dN3 + n1*(1.0*n.fh*n.nop)*sqrt(!pi/a)

;a = (h1^2 + h.o2p^2)/(h1^2 * h.o2p^2)
;dN3 = dN3 + n1*(2.0*n.fh*n.no2p)*sqrt(!pi/a)

;print,'dN3...',dN3

dn3 = dN3/(sqrt(!pi)*h3)


return,dn3
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function ftavg_n,n1,nn,h1,hn,h3
;ions with neutrals, zoff is offset along z in units of h.
;h1 = ion scale height
;hn = neutral scale height
;h3 = originating species scale height, dni/dt
;-----------------------------------------------------------------------
@common

a = (h1^2 + hn^2)/(h1^2 * hn^2)
b = -2*zoff/hn^2
c = zoff^2/hn^2
dN3 = n1*nn*sqrt(!pi/a)*exp((b^2 - 4.0*a*c)/(4.0*a))

dn3 = dN3/(sqrt(!pi)*h3)

return,dn3
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function ftavg_n_el,nn,n,hn,h,h3
;neutrals with electrons, sum over each ion species to get total
;electron contribution.
;-----------------------------------------------------------------------
@common

dN3 = 0.0

a = (hn^2 + h.sp^2)/(hn^2 * h.sp^2)
b = 2*zoff/hn^2
c = zoff^2/hn^2
dN3 = dN3 + nn*(1.0*n.nsp)*sqrt(!pi/a)*exp((b^2 - 4.0*a*c)/(4.0*a)) 

a = (hn^2 + h.s2p^2)/(hn^2 * h.s2p^2)
b = 2*zoff/hn^2
c = zoff^2/hn^2
dN3 = dN3 + nn*(2.0*n.ns2p)*sqrt(!pi/a)*exp((b^2 - 4.0*a*c)/(4.0*a)) 

a = (hn^2 + h.s3p^2)/(hn^2 * h.s3p^2)
b = 2*zoff/hn^2
c = zoff^2/hn^2
dN3 = dN3 + nn*(3.0*n.ns3p)*sqrt(!pi/a)*exp((b^2 - 4.0*a*c)/(4.0*a)) 

a = (hn^2 + h.op^2)/(hn^2 * h.op^2)
b = 2*zoff/hn^2
c = zoff^2/hn^2
dN3 = dN3 + nn*(1.0*n.nop)*sqrt(!pi/a)*exp((b^2 - 4.0*a*c)/(4.0*a)) 

a = (hn^2 + h.o2p^2)/(hn^2 * h.o2p^2)
b = 2*zoff/hn^2
c = zoff^2/hn^2
dN3 = dN3 + nn*(2.0*n.no2p)*sqrt(!pi/a)*exp((b^2 - 4.0*a*c)/(4.0*a)) 

dn3 = dN3/(sqrt(!pi)*h3)


return,dn3
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function ftavg_n_elh,nn,n,hn,h,h3
;neutrals with electrons, sum over each ion species to get total
;electron contribution.
;-----------------------------------------------------------------------
@common

dN3 = nn*n.fh*n.nel*sqrt(!pi)*hn

;a = (hn^2 + h.sp^2)/(hn^2 * h.sp^2)
;b = 2*zoff/hn^2
;c = zoff^2/hn^2
;dN3 = dN3 + nn*(1.0*n.fh*n.nsp)*sqrt(!pi/a)*exp((b^2 - 4.0*a*c)/(4.0*a)) 

;a = (hn^2 + h.s2p^2)/(hn^2 * h.s2p^2)
;b = 2*zoff/hn^2
;c = zoff^2/hn^2
;dN3 = dN3 + nn*(2.0*n.fh*n.ns2p)*sqrt(!pi/a)*exp((b^2 - 4.0*a*c)/(4.0*a)) 

;a = (hn^2 + h.s3p^2)/(hn^2 * h.s3p^2)
;b = 2*zoff/hn^2
;c = zoff^2/hn^2
;dN3 = dN3 + nn*(3.0*n.fh*n.ns3p)*sqrt(!pi/a)*exp((b^2 - 4.0*a*c)/(4.0*a)) 

;a = (hn^2 + h.op^2)/(hn^2 * h.op^2)
;b = 2*zoff/hn^2
;c = zoff^2/hn^2
;dN3 = dN3 + nn*(1.0*n.fh*n.nop)*sqrt(!pi/a)*exp((b^2 - 4.0*a*c)/(4.0*a)) 

;a = (hn^2 + h.o2p^2)/(hn^2 * h.o2p^2)
;b = 2*zoff/hn^2
;c = zoff^2/hn^2
;dN3 = dN3 + nn*(2.0*n.fh*n.no2p)*sqrt(!pi/a)*exp((b^2 - 4.0*a*c)/(4.0*a)) 

dn3 = dN3/(sqrt(!pi)*h3)


return,dn3
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
pro get_scale_heights,h,T,n
;all h in km
;-----------------------------------------------------------------------
@common

omega = 1.76e-4   ;angular frequency of Jupiter's rotation rad/sec
mp = 1.67e-27
ms = 32.*mp 
mo = 16.*mp

;h.s = 6.0*7.14e4*5.0*!dtor
;h.o = 6.0*7.14e4*5.0*!dtor

h.s = 6.0*7.14e4*n_width*!dtor
h.o = 6.0*7.14e4*n_width*!dtor

h.sp = sqrt(2.0*T.Tsp*1.6e-19*(1+T.Tel/T.Tsp)/(3.0*ms*omega^2))/1e3
h.s2p = sqrt(2.0*T.Ts2p*1.6e-19*(1+2.0*T.Tel/T.Ts2p)/(3.0*ms*omega^2))/1e3
h.s3p = sqrt(2.0*T.Ts3p*1.6e-19*(1+3.0*T.Tel/T.Ts3p)/(3.0*ms*omega^2))/1e3
h.op = sqrt(2.0*T.Top*1.6e-19*(1+T.Tel/T.Top)/(3.0*mo*omega^2))/1e3
h.o2p = sqrt(2.0*T.To2p*1.6e-19*(1+2.0*T.Tel/T.To2p)/(3.0*mo*omega^2))/1e3

;h.sp = sqrt(2.0*T.Tsp*1.6e-19/(3.0*ms*omega^2))/1e3
;h.s2p = sqrt(2.0*T.Ts2p*1.6e-19/(3.0*ms*omega^2))/1e3
;h.s3p = sqrt(2.0*T.Ts3p*1.6e-19/(3.0*ms*omega^2))/1e3
;h.op = sqrt(2.0*T.Top*1.6e-19/(3.0*mo*omega^2))/1e3
;h.o2p = sqrt(2.0*T.To2p*1.6e-19/(3.0*mo*omega^2))/1e3

wsp = n.nsp
ws2p = 2.0*n.ns2p 
ws3p = 3.0*n.ns3p  
wop = n.nop   
wo2p = 2.0*n.ns2p
Ntot = wsp+ws2p+ws3p+wop+wo2p  

h.el = (wsp*h.sp + ws2p*h.s2p + ws3p*h.s3p + wop*h.op + wo2p*h.o2p)/(Ntot)

return
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
pro get_rates,r,T,h,Lshl,dLshl
;-----------------------------------------------------------------------
@common

;Electron impact ionization (Schreier and Eviatar 1998)
;r.is = 7.4e-8*sqrt(T.Tel)*exp(-10.36/T.Tel)     ;check this!!! 
;r.isp = 1.274e-8*sqrt(T.Tel)*exp(-23.1/T.Tel)  
r.is2p = 5.8e-9*sqrt(T.Tel)*exp(-34.88/T.Tel)  

r.is = cfit(16,16,T.Tel,c)
r.isp = cfit(16,15,T.Tel,c)
;r.is2p = cfit(16,14,T.Tel,c)

;r.ish = 1.464e-6*alog(T.Telh/27.37)/sqrt(T.Telh)
;r.isph = 4.44e-7*alog(T.Telh/27.07)/sqrt(T.Telh)
;r.is2ph = 1.76e-7*alog(T.Telh/27.07)/sqrt(T.Telh)

;r.ish = cfit(16,16,T.Telh,c)
;r.isph = cfit(16,15,T.Telh,c)
;r.is2ph = cfit(16,14,T.Telh,c)

r.io = 1.015e-8*sqrt(T.Tel)*exp(-14.54/T.Tel)
r.iop = 3.78e-9*sqrt(T.Tel)*exp(-33.42/T.Tel)

;r.io = cfit(8,8,T.Tel,c)
;r.iop = cfit(8,7,T.Tel,c)

;r.ioh = 5.74e-7*alog(T.Telh/27.07)/sqrt(T.Telh)
;r.ioph = 1.76e-7*alog(T.Telh/27.07)/sqrt(T.Telh)

;r.ioh = cfit(8,8,T.Telh,c)
;r.ioph = cfit(8,7,T.Telh,c)

; Total electron recombination rates from the work of Sultana Nahar                                          
;     S II, S III, O III: (ApJS, 101, 423 [1995]; ApJS 106, 213 [1996])                                      
;     O I, O II:          (ApJS, 120, 131 [1999])                                                            
; Recombination rates for S I and S IV come from Mazzotta et al 1998, which                                  
; provides formulae for the dielectronic recombination rate. The                                             
; radiative recombination rate comes from Dima Verners rrfit.f code,                                         
; which uses: Shull & Van Steenberg, (1982, ApJS, 48, 95) for Sulfur                                         
;             Pequignot et al. (1991, A&A, 251, 680) for Oxygen                                              
recombination_rates,T.Tel,rsp,rs2p,rs3p,rop,ro2p
r.rsp  = rsp
r.rs2p = rs2p                                                         
r.rs3p = rs3p                                                                           
r.rop  = rop
r.ro2p = ro2p        
;print,'recombination rates...',r.rs2p,r.rs3p,r.ro2p


;lnt = alog(T.Tel)
;A = 1.0
;y = -23.9+A*cos((lnt*(1.3*!pi)/(alog(10)) + 0.74*!pi))
;r.rs3p = exp(y)

;lnt = alog(T.Tel)
;A = 1.4
;y = -25.2+A*cos((lnt*(1.3*!pi)/(alog(10)) + 0.8*!pi))
;r.rs2p = exp(y)

;lnt = alog(T.Tel)
;A = 1.1
;y = -26.0+A*cos((lnt*(1.2*!pi)/(alog(10)) + 0.74*!pi))
;r.ro2p = exp(y)

;dielectron recombination (Smith and Strobel 1985)
;if (T.Tel le 6.9) then r.rs3p = 5.4e-11*T.Tel^0.53
;if (T.Tel gt 6.9) then r.rs3p = 4.3e-12*T.Tel^1.84  

;r.rs2p = r.rs3p/3.0   ;massive fudge based on Johnson and Strobel '82

;r.ro2p = r.rs3p/7.0   ;another massive fudge base on Johnson and Stobel

;charge exchange

r.cx_k0 = k0*8.1e-9       ;S+ + S++ -> S++ + S+
r.cx_k1 = k1*2.4e-8
r.cx_k2 = k2*3.0e-10
r.cx_k3 = k3*7.8e-9
r.cx_k4 = k4*1.32e-8
r.cx_k5 = k5*1.32e-8
r.cx_k6 = k6*5.2e-10
r.cx_k7 = k7*5.4e-9
;r.cx_k8 = k8*9.0e-11   ;Schreier
r.cx_k8 = k8*6.0e-11    ;McGrath
r.cx_k9 = k9*3.0e-9
r.cx_k10 = k10*2.34e-8
r.cx_k11 = k11*1.62e-8
r.cx_k12 = k12*2.3e-9  ;new
;r.cx_k12 = k12*7.3e-9   ;old
r.cx_k13 = k13*1.4e-9   ;Schreier for L=6
r.cx_k14 = k14*1.92e-8 ;Schreier
;r.cx_k14 = k14*4.3e-9 ;Barbosa94, Johnson82
;r.cx_k15 = k15*9e-9    ;Schreier for L=6
r.cx_k15 = k15*9e-10    ;McGrath for L=6
r.cx_k16 = k16*3.6e-10  ;McGrath

;print,k4,k10,k11,k13,k14,k15

Rj = 7.14e4*1e5 ;cm
;a = 6.0*Rj
;b = 6.5*Rj
;vol = 0.25*!pi^2*(a+b)*(b-a)^2
;Area = !pi*(b^2 - a^2)
dLshl2 = dLshl/2
a = (Lshl-dLshl2)*Rj                                                                         
b = (Lshl+dLshl2)*Rj                                                                 
Hn = h.s*1e5                                                                          

;vol = sqrt(!pi)*Hn*!pi*(b^2 - a^2)   
vol = Hn*!pi*(b^2 - a^2)   
;vol = 2*Hn*!pi*(b^2 - a^2)   
;vol = 2.5e+31
;print,'vol...',vol

fs = 1.0/(1.0+otos)
fo = otos/(1.0+otos)
;print,fs,fo

r.S_production = fs*net_source/vol
r.O_production = fo*net_source/vol

;Calculate the equatorial production rates based on total source rate.
;r.S_production = fs*net_source/(sqrt(!pi)*(h.s*1e5)*Area)
;r.O_production = fo*net_source/(sqrt(!pi)*(h.o*1e5)*Area)

;print,'S source...',r.S_production,Lshl
;print,'O source...',r.O_production,Lshl

;print,'cm3 el vol...',h.el*1e5*Area
;print,'cm3 o vol...',h.op*1e5*Area


return
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
function F_s,n,h,r,S,L 
;-----------------------------------------------------------------------
S = r.S_production ;+ r.cx_k8*ftavg(n.no,n.nsp)
L = r.is*ftavg_n_el(n.ns,n,h.s,h,h.s) + r.ish*ftavg_n_elh(n.ns,n,h.s,h,h.s) + $ 
    r.cx_k4*ftavg_n(n.ns3p,n.ns,h.s3p,h.s,h.s) + r.cx_k11*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.s) + $
    r.cx_k10*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.s) + $
    r.cx_k2*ftavg_n(n.ns2p,n.ns,h.s2p,h.s,h.s) + r.cx_k9*ftavg_n(n.nop,n.ns,h.op,h.s,h.s) + $
    r.cx_k3*ftavg_n(n.ns2p,n.ns,h.s2p,h.s,h.s) + r.cx_k1*ftavg_n(n.nsp,n.ns,h.sp,h.s,h.s) 

F_s = S - L

return,F_s
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_sp,n,h,r,S,L
;-----------------------------------------------------------------------  
S = r.is*ftavg_n_el(n.ns,n,h.s,h,h.sp) + r.ish*ftavg_n_elh(n.ns,n,h.s,h,h.sp) + $
    r.rs2p*ftavg_el(n.ns2p,n,h.s2p,h,h.sp) + $
    r.cx_k4*ftavg_n(n.ns3p,n.ns,h.s3p,h.s,h.sp) + r.cx_k10*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.sp) + $
    2.0*r.cx_k2*ftavg_n(n.ns2p,n.ns,h.s2p,h.s,h.sp) + $
    r.cx_k9*ftavg_n(n.nop,n.ns,h.op,h.s,h.sp) + r.cx_k12*ftavg_n(n.ns2p,n.no,h.s2p,h.o,h.sp)
L = r.isp*ftavg_el(n.nsp,n,h.sp,h,h.sp) + r.cx_k13*ftavg(n.no2p,n.nsp,h.o2p,h.sp,h.sp) + $
    n.nsp*r.Transport + $
    r.isph*ftavg_elh(n.nsp,n,h.sp,h,h.sp) + r.cx_k8*ftavg_n(n.nsp,n.no,h.sp,h.o,h.sp) + $
    r.cx_k16*ftavg(n.ns3p,n.nsp,h.s3p,h.sp,h.sp) + r.rsp*ftavg_el(n.nsp,n,h.sp,h,h.sp)

F_sp = S - L

return,F_sp
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_s2p,n,h,r,S,L
;-----------------------------------------------------------------------  
S = r.isp*ftavg_el(n.nsp,n,h.sp,h,h.s2p) + r.isph*ftavg_elh(n.nsp,n,h.sp,h,h.s2p) + $
    r.rs3p*ftavg_el(n.ns3p,n,h.s3p,h,h.s2p) + r.cx_k13*ftavg(n.no2p,n.nsp,h.o2p,h.sp,h.s2p) + $
    r.cx_k14*ftavg_n(n.ns3p,n.no,h.s3p,h.o,h.s2p) + r.cx_k4*ftavg_n(n.ns3p,n.ns,h.s3p,h.s,h.s2p) + $
    r.cx_k11*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.s2p) + 2.0*r.cx_k16*ftavg(n.ns3p,n.nsp,h.s3p,h.sp,h.s2p)

L = r.is2p*ftavg_el(n.ns2p,n,h.s2p,h,h.s2p) + $
    r.is2ph*ftavg_elh(n.ns2p,n,h.s2p,h,h.s2p) +  $
    r.rs2p*ftavg_el(n.ns2p,n,h.s2p,h,h.s2p) + $
    n.ns2p*r.Transport + $
    r.cx_k2*ftavg_n(n.ns2p,n.ns,h.s2p,h.s,h.s2p) + $        ;S + S++ -> S+ + S+
    r.cx_k12*ftavg_n(n.ns2p,n.no,h.s2p,h.o,h.s2p) + $       ;O + S++ -> O+ + S+ 
    r.cx_k15*ftavg(n.no2p,n.ns2p,h.o2p,h.s2p,h.s2p)         ;O++ + S++ -> O+ + S+++

F_s2p = S - L

return,F_s2p
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_s3p,n,h,r,S,L
;-----------------------------------------------------------------------  
S = r.is2p*ftavg_el(n.ns2p,n,h.s2p,h,h.s3p) + r.is2ph*ftavg_elh(n.ns2p,n,h.s2p,h,h.s3p) + $
    r.cx_k15*ftavg(n.no2p,n.ns2p,h.o2p,h.s2p,h.s3p)
L = r.rs3p*ftavg_el(n.ns3p,n,h.s3p,h,h.s3p) + n.ns3p*r.Transport + $
    r.cx_k14*ftavg_n(n.ns3p,n.no,h.s3p,h.o,h.s3p) + $
    r.cx_k4*ftavg_n(n.ns3p,n.ns,h.s3p,h.s,h.s3p) + r.cx_k16*ftavg(n.ns3p,n.nsp,h.s3p,h.sp,h.s3p)

F_s3p = S - L

return,F_s3p
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------
function F_o,n,h,r,S,L
;-----------------------------------------------------------------------
S = r.O_production ;+ r.cx_k9*n.ns*n.nop
L = r.io*ftavg_n_el(n.no,n,h.o,h,h.o) + r.ioh*ftavg_n_elh(n.no,n,h.o,h,h.o) + $
    r.cx_k14*ftavg_n(n.ns3p,n.no,h.s3p,h.o,h.o) + $
    r.cx_k6*ftavg_n(n.no2p,n.no,h.o2p,h.o,h.o) + r.cx_k8*ftavg_n(n.nsp,n.no,h.sp,h.o,h.o) + $
    r.cx_k12*ftavg_n(n.ns2p,n.no,h.s2p,h.o,h.o) + $
    r.cx_k5*ftavg_n(n.nop,n.no,h.op,h.o,h.o) + r.cx_k7*ftavg_n(n.no2p,n.no,h.o2p,h.o,h.o)

F_o = S - L

return,F_o
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_op,n,h,r,S,L
;-----------------------------------------------------------------------  
S = r.io*ftavg_n_el(n.no,n,h.o,h,h.op) + r.ioh*ftavg_n_elh(n.no,n,h.o,h,h.op) + $
    r.ro2p*ftavg_el(n.no2p,n,h.o2p,h,h.op) + r.cx_k13*ftavg(n.no2p,n.nsp,h.o2p,h.sp,h.op) + $
    r.cx_k15*ftavg(n.no2p,n.ns2p,h.o2p,h.s2p,h.op) + r.cx_k14*ftavg_n(n.ns3p,n.no,h.s3p,h.o,h.op) + $
    r.cx_k11*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.op) + r.cx_k10*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.op) + $
    r.cx_k8*ftavg_n(n.nsp,n.no,h.sp,h.o,h.op) + r.cx_k12*ftavg_n(n.ns2p,n.no,h.s2p,h.o,h.op) + $
    2.0*r.cx_k6*ftavg_n(n.no2p,n.no,h.o2p,h.o,h.op)
L = r.iop*ftavg_el(n.nop,n,h.op,h,h.op) + n.nop*r.Transport + $
    r.ioph*ftavg_elh(n.nop,n,h.op,h,h.op) + $
    r.cx_k9*ftavg_n(n.nop,n.ns,h.op,h.s,h.op) + r.rop*ftavg_el(n.nop,n,h.op,h,h.op)

F_op = S - L

return,F_op
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_o2p,n,h,r,S,L
;-----------------------------------------------------------------------  
S = r.iop*ftavg_el(n.nop,n,h.op,h,h.o2p) + r.ioph*ftavg_elh(n.nop,n,h.op,h,h.o2p)
L = r.ro2p*ftavg_el(n.no2p,n,h.o2p,h,h.o2p) + r.cx_k13*ftavg(n.no2p,n.nsp,h.o2p,h.sp,h.o2p) + $
    r.cx_k15*ftavg(n.no2p,n.ns2p,h.o2p,h.s2p,h.o2p) + n.no2p*r.Transport + $
    r.cx_k11*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.o2p) + r.cx_k10*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.o2p) + $
    r.cx_k6*ftavg_n(n.no2p,n.no,h.o2p,h.o,h.o2p)

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

Tpu = (2./3.)*0.5*m*mp*vrel^2   ;J
Tpu = Tpu/1.6e-19       ;eV

return,Tpu
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function nu_ie,mu,z,nelec,ni,Tel,Ti,hel,hion
; m is atomic mass
; z is atomic number
; for collision of given ion on thermal electron population
;-----------------------------------------------------------------------  
me = double(9.1e-28)    ;g
mp = double(1.67e-24)   ;g

mi = mu*mp

dz = hel/10
nz = 20
zz = findgen(nz)*dz

nel = nelec*exp(-zz^2/hel^2)
nion = ni*exp(-zz^2/hion^2)
nuarr = fltarr(nz)
lambda = fltarr(nz)

if ((Ti*me/mi lt Tel) and (Tel le 10*Z^2)) then begin
   lambda = 23.0 - alog(sqrt(nel)*z*Tel^(-3./2.))
endif

if ((Tel gt Ti*me/mi) and (Tel ge 10*Z^2)) then begin
   lambda = 24.0 - alog(sqrt(nel)/Tel)
endif

if (Tel lt Ti*Z*me/mi) then begin
   lambda = 30.0 - alog(sqrt(nion)*Ti^(-3./2.)*Z^2/mu)
endif

nu_ie = 1.8e-19*sqrt(me*mi)*(Z^2)*nel(0)*nion(0)*lambda(0)/(mi*Tel + me*Ti)^(3./2.)

for k = 0,nz-1 do begin
   nuarr(k) = 1.8e-19*sqrt(me*mi)*(Z^2)*nel(k)*nion(k)*lambda(k)/(mi*Tel + me*Ti)^(3./2.)
endfor

;print,'nu_ie(0)...',nu_ie
nu_ie = total(nuarr*nion*nel)/total(nion*nel)
;print,'<nu_ie>...',nu_ie



return,nu_ie
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function nu_ieh,mu,z,nelec_h,ni,Tel,Ti,hion
; m is atomic mass
; z is atomic number
; for collision of given ion on thermal electron population
;-----------------------------------------------------------------------  
me = double(9.1e-28)    ;g
mp = double(1.67e-24)   ;g

mi = mu*mp

dz = hion/10
nz = 20
zz = findgen(nz)*dz

nion = ni*exp(-zz^2/hion^2)
nuarr = fltarr(nz)
lambda = fltarr(nz)
nel = fltarr(nz)
nel(*) = nelec_h

if ((Ti*me/mi lt Tel) and (Tel le 10*Z^2)) then begin
   lambda = 23.0 - alog(sqrt(nel)*z*Tel^(-3./2.))
endif

if ((Tel gt Ti*me/mi) and (Tel ge 10*Z^2)) then begin
   lambda = 24.0 - alog(sqrt(nel)/Tel)
endif

if (Tel lt Ti*Z*me/mi) then begin
   lambda = 30.0 - alog(sqrt(nion)*Ti^(-3./2.)*Z^2/mu)
endif

nu_ieh = 1.8e-19*sqrt(me*mi)*(Z^2)*nel(0)*nion(0)*lambda(0)/(mi*Tel + me*Ti)^(3./2.)

for k = 0,nz-1 do begin
   nuarr(k) = 1.8e-19*sqrt(me*mi)*(Z^2)*nel(k)*nion(k)*lambda(k)/(mi*Tel + me*Ti)^(3./2.)
endfor

;print,'nu_ieh(0)...',nu_ieh
nu_ieh = total(nuarr*nion)/total(nion)
;print,'<nu_ieh>...',nu_ieh



return,nu_ieh
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function nu_ei,mu,z,nelec,ni,Tel,Ti,hel,hion
;for electron energy equation
; m is atomic mass
; z is atomic number
; for collision of given ion on thermal electron population
;-----------------------------------------------------------------------  
me = double(9.1e-28)    ;g
mp = double(1.67e-24)   ;g

mi = mu*mp

dz = hel/10
nz = 20
zz = findgen(nz)*dz

nel = nelec*exp(-zz^2/hel^2)
nion = ni*exp(-zz^2/hion^2)
nuarr = fltarr(nz)
lambda = fltarr(nz)

if ((Ti*me/mi lt Tel) and (Tel le 10*z^2)) then begin
   lambda = 23.0 - alog(sqrt(nel)*z*Tel^(-3./2.))
endif

if ((Tel gt Ti*me/mi) and (Tel ge 10*Z^2)) then begin
   lambda = 24.0 - alog(sqrt(nel)/Tel)
endif

if (Tel lt Ti*Z*me/mi) then begin
   lambda = 30.0 - alog(sqrt(nion)*Ti^(-3./2.)*Z^2/mu)
endif

nu_ei = 1.8e-19*sqrt(me*mi)*(Z^2)*nel(0)*nion(0)*lambda(0)/(mi*Tel + me*Ti)^(3./2.)

for k = 0,nz-1 do begin
   nuarr(k) = 1.8e-19*sqrt(me*mi)*(Z^2)*nel(k)*nion(k)*lambda(k)/(mi*Tel + me*Ti)^(3./2.)
endfor

;print,'nu_ei(0)...',nu_ei
nu_ei = total(nuarr*nel*nion)/total(nel*nion)
;print,'<nu_ei>...',nu_ei



return,nu_ei
end
;-----------------------------------------------------------------------  




;-----------------------------------------------------------------------  
function nu_ee,n,h,Tel,Telh
; m is atomic mass
; z is atomic number
; for collision of given ion on thermal electron population
;-----------------------------------------------------------------------  
me = double(9.1e-28)    ;g
mp = double(1.67e-24)   ;g

nelec_c = n.nel
nelec_h = n.fh*n.nel

dz = h.el/10
nz = 20
z = findgen(nz)*dz

nel = nelec_c*exp(-z^2/h.el^2)
nuarr = fltarr(nz)
lambda = fltarr(nz)

if (Tel le 10.0) then begin
   lambda = 23.0 - alog(sqrt(nel)*Tel^(-3./2.))
endif

if (Tel gt 10.0) then begin
   lambda = 24.0 - alog(sqrt(nel)/Tel)
endif

nu_ee = 1.8e-19*sqrt(me*me)*nel(0)*nelec_h*lambda(0)/(me*Tel + me*Telh)^(3./2.)

for k = 0,nz-1 do begin
   nuarr(k) = 1.8e-19*sqrt(me*me)*nel(k)*nelec_h*lambda(k)/(me*Tel + me*Telh)^(3./2.)
endfor

;print,'nu_ee(0)...',nu_ee
nu_ee = total(nuarr*nel)/total(nel)
;print,'<nu_ee>...',nu_ee


return,nu_ee
end
;-----------------------------------------------------------------------  



;-----------------------------------------------------------------------  
function nu_ii,mu1,mu2,z1,z2,n1,n2,T1,T2,h1,h2
; m is atomic mass
; z is atomic number
; for collision of given ion on ion population
;-----------------------------------------------------------------------  
mp = 1.67e-24   ;g

m1 = double(mu1*mp)
m2 = double(mu2*mp)


dz = h1/10
nz = 20
zz = findgen(nz)*dz

n1arr = n1*exp(-zz^2/h1^2)
n2arr = n2*exp(-zz^2/h2^2)
nuarr = fltarr(nz)
lambda = fltarr(nz)

lambda = 23.0 - alog((Z1*Z2*(mu1+mu2)/(mu1*T2 + mu2*T1))* $
                     sqrt((n1arr*Z1^2)/T1 + (n2arr*Z2^2)/T2))

;lambda_av = total(lambda_arr*n1arr*n2arr)/total(n1arr*n2arr)

;lambda = 23.0 - alog((Z1*Z2*(mu1+mu2)/(mu1*T2 + mu2*T1))* $
;                     sqrt((n1*Z1^2)/T1 + (n2*Z2^2)/T2))

nu_ii = 1.8e-19*sqrt(m1*m2)*(Z1^2)*(Z2^2)*n1arr*n2arr*lambda/(m1*T2 + m2*T1)^(3./2.)

nu_ii = total(nu_ii*n1arr*n2arr)/total(n1arr*n2arr)

return,nu_ii
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro update_temp,n,nT,T
;-----------------------------------------------------------------------  
@common

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
pro ft_rad,n,T,r,h,Lavg,rad_tot_0
;-----------------------------------------------------------------------  
;nel = n.fc*n.nel
nel = n.nel

dz = h.el/10
nz = 20
z = findgen(nz)*dz

rad_sp = fltarr(nz)
rad_s2p = fltarr(nz)
rad_s3p = fltarr(nz)
rad_op = fltarr(nz)
rad_o2p = fltarr(nz)
rad_tot = fltarr(nz)
deltaT = fltarr(nz)

for k = 0,nz-1 do begin

;print,k,z(k)

nel = n.nel*exp(-z(k)^2 / h.el^2)

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

emisSII = r.emisSII(i,j)*A4 + r.emisSII(ip,j)*A3 + r.emisSII(i,jp)*A2 + $
          r.emisSII(ip,jp)*A1
L = emisSII*n.nsp
rad_sp(k) = emisSII*n.nsp*exp(-z(k)^2/h.sp^2)

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

emisSIII = r.emisSIII(i,j)*A4 + r.emisSIII(ip,j)*A3 + r.emisSIII(i,jp)*A2 + $
          r.emisSIII(ip,jp)*A1
L = L+ emisSIII*n.ns2p
rad_s2p(k) = emisSIII*n.ns2p*exp(-z(k)^2/h.s2p^2)


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

emisSIV = r.emisSIV(i,j)*A4 + r.emisSIV(ip,j)*A3 + r.emisSIV(i,jp)*A2 + $
          r.emisSIV(ip,jp)*A1
L = L + emisSIV*n.ns3p
rad_s3p(k) = emisSIV*n.ns3p*exp(-z(k)^2/h.s3p^2)

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

emisOII = r.emisOII(i,j)*A4 + r.emisOII(ip,j)*A3 + r.emisOII(i,jp)*A2 + $
          r.emisOII(ip,jp)*A1
L = L + emisOII*n.nop
rad_op(k) = emisOII*n.nop*exp(-z(k)^2/h.op^2)

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

emisOIII = r.emisOIII(i,j)*A4 + r.emisOIII(ip,j)*A3 + r.emisOIII(i,jp)*A2 + $
          r.emisOIII(ip,jp)*A1
L = L + emisOIII*n.no2p
rad_o2p(k) = emisOIII*n.no2p*exp(-z(k)^2/h.o2p^2)

;Teq = nu_ei(32.0,1.0,n.nel,n.nsp,T.Tel,T.Tsp)*ftavg(n.nel,n.nsp,h.el,h.sp,h.el)*(T.Tsp-T.Tel) + $
;      nu_ei(32.0,2.0,n.nel,n.ns2p,T.Tel,T.Ts2p)*ftavg(n.nel,n.ns2p,h.el,h.s2p,h.el)*(T.Ts2p-T.Tel) + $
;      nu_ei(32.0,3.0,n.nel,n.ns3p,T.Tel,T.Ts3p)*ftavg(n.nel,n.ns3p,h.el,h.s3p,h.el)*(T.Ts3p-T.Tel) + $
;      nu_ei(16.0,1.0,n.nel,n.nop,T.Tel,T.Top)*ftavg(n.nel,n.nop,h.el,h.op,h.el)*(T.Top-T.Tel) + $
;      nu_ei(16.0,2.0,n.nel,n.no2p,T.Tel,T.To2p)*ftavg(n.nel,n.no2p,h.el,h.o2p,h.el)*(T.To2p-T.Tel) + $
;      nu_ee(n,h,T.Tel,T.Telh)*(T.Telh - T.Tel)

;EF_el = Teq - (2./3.)*L - r.Transport*n.nel*T.Tel

;r.Puv = r.psp*n.nsp + r.ps2p*n.ns2p + r.ps3p*n.ns3p + r.pop*n.nop + $
;        r.po2p*n.no2p

;r.Puv = L


endfor

rad_tot = rad_sp + rad_s2p + rad_s3p + rad_op + rad_o2p

nel = n.nel*exp(-z^2/h.el^2)

Lavg = total(rad_tot*nel)/total(nel)
rad_tot_0 = rad_tot(0)

;plot,z,rad_s2p,linestyle=0
;oplot,z,rad_sp,linestyle=2
;oplot,z,rad_s3p,linestyle=3
;oplot,z,rad_op,linestyle=4
;oplot,z,rad_o2p,linestyle=5

return
end
;-----------------------------------------------------------------------  



;-----------------------------------------------------------------------  
function EF_el,n,T,nT,r,h
;-----------------------------------------------------------------------  
;nel = n.fc*n.nel
nel = n.nel

;i=where(abs(r.emistemp-T.Tel) eq min(abs(r.emistemp-T.Tel)))
;j=where(abs(r.emisden-nel) eq min(abs(r.emisden-nel)))

;i = i(0)
;j = j(0)

;if (T.Tel gt r.emistemp(i)) then ip = i+1
;if (T.Tel le r.emistemp(i)) then begin
;   ip = i
;   i = ip-1
;endif

;if (nel gt r.emisden(j)) then jp = j+1
;if (nel le r.emisden(j)) then begin
;   jp = j
;   j = jp-1
;endif

;dT = r.emistemp(ip) - r.emistemp(i)
;dn = r.emisden(jp) - r.emisden(j)

;A = dT*dn
;dT1 = T.Tel - r.emistemp(i)
;dT2 = r.emistemp(ip) - T.Tel
;dn1 = nel - r.emisden(j)
;dn2 = r.emisden(jp) - nel
;A1 = dT1*dn1/A
;A2 = dT2*dn1/A
;A3 = dT1*dn2/A
;A4 = dT2*dn2/A

;emisSII = r.emisSII(i,j)*A4 + r.emisSII(ip,j)*A3 + r.emisSII(i,jp)*A2 + $
;          r.emisSII(ip,jp)*A1
;L = emisSII*n.nsp
;r.psp = emisSII

;dT = r.emistemp(ip) - r.emistemp(i)
;dn = r.emisden(jp) - r.emisden(j)

;A = dT*dn
;dT1 = T.Tel - r.emistemp(i)
;dT2 = r.emistemp(ip) - T.Tel
;dn1 = nel - r.emisden(j)
;dn2 = r.emisden(jp) - nel
;A1 = dT1*dn1/A
;A2 = dT2*dn1/A
;A3 = dT1*dn2/A
;A4 = dT2*dn2/A

;emisSIII = r.emisSIII(i,j)*A4 + r.emisSIII(ip,j)*A3 + r.emisSIII(i,jp)*A2 + $
;          r.emisSIII(ip,jp)*A1
;L = L+ emisSIII*n.ns2p
;r.ps2p = emisSIII


;dT = r.emistemp(ip) - r.emistemp(i)
;dn = r.emisden(jp) - r.emisden(j)

;A = dT*dn
;dT1 = T.Tel - r.emistemp(i)
;dT2 = r.emistemp(ip) - T.Tel
;dn1 = nel - r.emisden(j)
;dn2 = r.emisden(jp) - nel
;A1 = dT1*dn1/A
;A2 = dT2*dn1/A
;A3 = dT1*dn2/A
;A4 = dT2*dn2/A

;emisSIV = r.emisSIV(i,j)*A4 + r.emisSIV(ip,j)*A3 + r.emisSIV(i,jp)*A2 + $
;          r.emisSIV(ip,jp)*A1
;L = L + emisSIV*n.ns3p
;r.ps3p = emisSIV

;dT = r.emistemp(ip) - r.emistemp(i)
;dn = r.emisden(jp) - r.emisden(j)

;A = dT*dn
;dT1 = T.Tel - r.emistemp(i)
;dT2 = r.emistemp(ip) - T.Tel
;dn1 = nel - r.emisden(j)
;dn2 = r.emisden(jp) - nel
;A1 = dT1*dn1/A
;A2 = dT2*dn1/A
;A3 = dT1*dn2/A
;A4 = dT2*dn2/A

;emisOII = r.emisOII(i,j)*A4 + r.emisOII(ip,j)*A3 + r.emisOII(i,jp)*A2 + $
;          r.emisOII(ip,jp)*A1
;L = L + emisOII*n.nop
;r.pop = emisOII

;dT = r.emistemp(ip) - r.emistemp(i)
;dn = r.emisden(jp) - r.emisden(j)

;A = dT*dn
;dT1 = T.Tel - r.emistemp(i)
;dT2 = r.emistemp(ip) - T.Tel
;dn1 = nel - r.emisden(j)
;dn2 = r.emisden(jp) - nel
;A1 = dT1*dn1/A
;A2 = dT2*dn1/A
;A3 = dT1*dn2/A
;A4 = dT2*dn2/A

;emisOIII = r.emisOIII(i,j)*A4 + r.emisOIII(ip,j)*A3 + r.emisOIII(i,jp)*A2 + $
;          r.emisOIII(ip,jp)*A1
;L = L + emisOIII*n.no2p
;r.po2p = emisOIII

Teq = nu_ei(32.0,1.0,n.nel,n.nsp,T.Tel,T.Tsp,h.el,h.sp)*(T.Tsp-T.Tel) + $
      nu_ei(32.0,2.0,n.nel,n.ns2p,T.Tel,T.Ts2p,h.el,h.s2p)*(T.Ts2p-T.Tel) + $
      nu_ei(32.0,3.0,n.nel,n.ns3p,T.Tel,T.Ts3p,h.el,h.s3p)*(T.Ts3p-T.Tel) + $
      nu_ei(16.0,1.0,n.nel,n.nop,T.Tel,T.Top,h.el,h.op)*(T.Top-T.Tel) + $
      nu_ei(16.0,2.0,n.nel,n.no2p,T.Tel,T.To2p,h.el,h.o2p)*(T.To2p-T.Tel) + $
      nu_ee(n,h,T.Tel,T.Telh)*(T.Telh - T.Tel)


;print,'L(0)....',L
ft_rad,n,T,r,h,L,rad_tot_0
;print,'<L>....',Lavg



EF_el = Teq - (2./3.)*L - r.Transport*n.nel*T.Tel
;EF_el = Teq - (2./3.)*Lavg - r.Transport*n.nel*T.Tel

;print,'Tel...',T.Tel

;r.Puv = r.psp*n.nsp + r.ps2p*n.ns2p + r.ps3p*n.ns3p + r.pop*n.nop + $
;        r.po2p*n.no2p

r.Puv = rad_tot_0
;r.Puv = Lavg

return,EF_el
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function EF_sp,n,T,nT,r,h
;-----------------------------------------------------------------------  
S = (r.is*ftavg_n_el(n.ns,n,h.s,h,h.sp) + $        ;cold e ionization 
     r.ish*ftavg_n_elh(n.ns,n,h.s,h,h.sp) + $       ;hot e ionization
     r.cx_k4*ftavg_n(n.ns3p,n.ns,h.s3p,h.s,h.sp) + $         ;S + S+++ -> S+ + S++
     r.cx_k10*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.sp) + $        ;S + O++ -> S+ + O+
     r.cx_k1*ftavg_n(n.nsp,n.ns,h.sp,h.s,h.sp) + $          ;S + S+ -> S+ + S
     r.cx_k2*ftavg_n(n.ns2p,n.ns,h.s2p,h.s,h.sp) + $         ;S + S++ -> S+ + S+ 
     r.cx_k9*ftavg_n(n.nop,n.ns,h.op,h.s,h.sp))*T.Tpu_s + $ ;S + O+ -> S+ + O
     r.cx_k2*ftavg_n(n.ns2p,n.ns,h.s2p,h.s,h.sp)*T.Ts2p + $  ;S + S++ -> S+ + S+
     r.cx_k12*ftavg_n(n.ns2p,n.no,h.s2p,h.o,h.sp)*T.Ts2p + $ ;O + S++ -> O+ + S+
     r.cx_k0*ftavg(n.nsp,n.ns2p,h.sp,h.s2p,h.sp)*T.Ts2p + $ ;S+ + S++ -> S++ + S+
     r.rs2p*ftavg_el(n.ns2p,n,h.s2p,h,h.sp)*T.Ts2p      ;S++ + e -> S+

r.Sei = (r.is*n.ns*n.nel + r.ish*n.ns*n.nelh)*T.Tpu_s
;r.Seih = r.ish*n.ns*n.nelh*T.Tpu_s
r.Secx = (r.cx_k4*n.ns*n.ns3p + r.cx_k10*n.ns*n.no2p + $
         r.cx_k1*n.ns*n.nsp + r.cx_k2*n.ns*n.ns2p + $
         r.cx_k9*n.ns*n.nop+r.cx_k3*n.ns*n.ns2p + $
         r.cx_k11*n.ns*n.no2p)*T.Tpu_s

L= r.isp*ftavg_el(n.nsp,n,h.sp,h,h.sp)*T.Tsp + $       ;cold e ionization 
   r.isph*ftavg_elh(n.nsp,n,h.sp,h,h.sp)*T.Tsp + $      ;hot e ionization
   r.cx_k1*ftavg_n(n.nsp,n.ns,h.sp,h.s,h.sp)*T.Tsp + $     ;S + S+ -> S+ + S
   r.cx_k8*ftavg_n(n.nsp,n.no,h.sp,h.o,h.sp)*T.Tsp + $     ;O + S+ -> O+ + S
   r.cx_k13*ftavg(n.no2p,n.nsp,h.o2p,h.sp,h.sp)*T.Tsp + $  ;O++ + S+ -> O+ + S++
   r.cx_k0*ftavg(n.nsp,n.ns2p,h.sp,h.s2p,h.sp)*T.Tsp + $   ;S+ + S++ -> S++ + S+
   r.cx_k16*ftavg(n.ns3p,n.nsp,h.s3p,h.sp,h.sp)*T.Tsp + $  ;S+++ + S+ -> S++ + S++
   r.Transport*n.nsp*T.Tsp + $
   r.rsp*ftavg_el(n.nsp,n,h.sp,h,h.sp)*T.Tsp ;S++ + e -> S+     

Teq= nu_ii(32.0,32.0,1.0,2.0,n.nsp,n.ns2p,T.Tsp,T.Ts2p,h.sp,h.s2p)*(T.Ts2p-T.Tsp)+$
   nu_ii(32.0,32.0,1.0,3.0,n.nsp,n.ns3p,T.Tsp,T.Ts3p,h.sp,h.s3p)*(T.Ts3p-T.Tsp)+$
   nu_ii(32.0,16.0,1.0,1.0,n.nsp,n.nop,T.Tsp,T.Top,h.sp,h.op)*(T.Top-T.Tsp)+$
   nu_ii(32.0,16.0,1.0,2.0,n.nsp,n.no2p,T.Tsp,T.To2p,h.sp,h.o2p)*(T.To2p-T.Tsp)+$
   nu_ie(32.0,1.0,n.nel,n.nsp,T.Tel,T.Tsp,h.el,h.sp)*(T.Tel-T.Tsp) + $
   nu_ieh(32.0,1.0,n.nelh,n.nsp,T.Telh,T.Tsp,h.sp)*(T.Telh-T.Tsp) 
  

EF_sp = S - L + Teq

return,EF_sp
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function EF_s2p,n,T,nT,r,h
;-----------------------------------------------------------------------  
S = (r.isp*ftavg_el(n.nsp,n,h.sp,h,h.s2p) + $               ;cold e ionization 
     r.isph*ftavg_elh(n.nsp,n,h.sp,h,h.s2p))*T.Tsp + $       ;hot e ionization
     r.cx_k3*ftavg_n(n.ns2p,n.ns,h.s2p,h.s,h.s2p)*T.Tpu_s + $      ;S + S++ -> S++ + S
     r.cx_k4*ftavg_n(n.ns3p,n.ns,h.s3p,h.s,h.s2p)*T.Ts3p + $       ;S + S+++ -> S+ + S++
     r.cx_k11*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.s2p)*T.Tpu_s + $     ;S + O++ -> S++ + O+ + e
     r.cx_k13*ftavg(n.no2p,n.nsp,h.o2p,h.sp,h.s2p)*T.Tsp + $      ;O++ + S+ -> O+ + S++ 
     r.cx_k14*ftavg_n(n.ns3p,n.no,h.s3p,h.o,h.s2p)*T.Ts3p + $      ;O + S+++ -> O+ + S++
     r.cx_k0*ftavg(n.nsp,n.ns2p,h.sp,h.s2p,h.s2p)*T.Tsp + $       ;S+ + S++ -> S++ + S+
     r.cx_k16*ftavg(n.ns3p,n.nsp,h.s3p,h.sp,h.s2p)*T.Tsp + $      ;S+++ + S+ -> S++ + S++
     r.cx_k16*ftavg(n.ns3p,n.nsp,h.s3p,h.sp,h.s2p)*T.Ts3p + $     ;S+++ + S+ -> S++ + S++
     r.rs3p*ftavg_el(n.ns3p,n,h.s3p,h,h.s2p)*T.Ts3p           ;S+++ + e -> S++     

L= r.is2p*ftavg_el(n.ns2p,n,h.s2p,h,h.s2p)*T.Ts2p + $       ;cold e ionization 
   r.is2ph*ftavg_elh(n.ns2p,n,h.s2p,h,h.s2p)*T.Ts2p + $      ;hot e ionization
   r.cx_k2*ftavg_n(n.ns2p,n.ns,h.s2p,h.s,h.s2p)*T.Ts2p + $     ;S + S++ -> S+ + S+
   r.cx_k3*ftavg_n(n.ns2p,n.ns,h.s2p,h.s,h.s2p)*T.Ts2p + $     ;S + S++ -> S++ + S
   r.cx_k12*ftavg_n(n.ns2p,n.no,h.s2p,h.o,h.s2p)*T.Ts2p + $    ;O + S++ -> O+ + S+
   r.cx_k15*ftavg(n.no2p,n.ns2p,h.o2p,h.s2p,h.s2p)*T.Ts2p + $  ;O++ + S++ -> O+ + S+++
   r.cx_k0*ftavg(n.nsp,n.ns2p,h.sp,h.s2p,h.s2p)*T.Ts2p + $    ;S+ + S++ -> S++ + S+
   r.Transport*n.ns2p*T.Ts2p + $   
   r.rs2p*ftavg_el(n.ns2p,n,h.s2p,h,h.s2p)*T.Ts2p

Teq= nu_ii(32.0,32.0,2.0,1.0,n.ns2p,n.nsp,T.Ts2p,T.Tsp,h.s2p,h.sp)*(T.Tsp-T.Ts2p)+$
   nu_ii(32.0,32.0,2.0,3.0,n.ns2p,n.ns3p,T.Ts2p,T.Ts3p,h.s2p,h.s3p)*(T.Ts3p-T.Ts2p)+$
   nu_ii(32.0,16.0,2.0,1.0,n.ns2p,n.nop,T.Ts2p,T.Top,h.s2p,h.op)*(T.Top-T.Ts2p)+$
   nu_ii(32.0,16.0,2.0,2.0,n.ns2p,n.no2p,T.Ts2p,T.To2p,h.s2p,h.o2p)*(T.To2p-T.Ts2p)+$
   nu_ie(32.0,2.0,n.nel,n.ns2p,T.Tel,T.Ts2p,h.el,h.s2p)*(T.Tel-T.Ts2p) + $
   nu_ieh(32.0,2.0,n.nelh,n.ns2p,T.Telh,T.Ts2p,h.s2p)*(T.Telh-T.Ts2p) 
 
EF_s2p = S - L + Teq

return,EF_s2p
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function EF_s3p,n,T,nT,r,h
;-----------------------------------------------------------------------  
S = (r.is2p*ftavg_el(n.ns2p,n,h.s2p,h,h.s3p) + $               ;cold e ionization 
     r.is2ph*ftavg_elh(n.ns2p,n,h.s2p,h,h.s3p))*T.Ts2p + $     ;hot e ionization
     r.cx_k15*ftavg(n.no2p,n.ns2p,h.o2p,h.s2p,h.s3p)*T.Ts2p         ;O++ + S++ -> O+ + S+++

L= r.cx_k14*ftavg_n(n.ns3p,n.no,h.s3p,h.o,h.s3p)*T.Ts3p + $      ;O + S+++ -> O+ + S++
   r.cx_k16*ftavg(n.ns3p,n.nsp,h.s3p,h.sp,h.s3p)*T.Ts3p + $     ;S+++ + S+ -> S++ + S++
   r.Transport*n.ns3p*T.Ts3p + $
   r.rs3p*ftavg_el(n.ns3p,n,h.s3p,h,h.s3p)*T.Ts3p           ;S+++ + e -> S++   

Teq= nu_ii(32.0,32.0,3.0,1.0,n.ns3p,n.nsp,T.Ts3p,T.Tsp,h.s3p,h.sp)*(T.Tsp-T.Ts3p)+$
   nu_ii(32.0,32.0,3.0,2.0,n.ns3p,n.ns2p,T.Ts3p,T.Ts2p,h.s3p,h.s2p)*(T.Ts2p-T.Ts3p)+$
   nu_ii(32.0,16.0,3.0,1.0,n.ns3p,n.nop,T.Ts3p,T.Top,h.s3p,h.op)*(T.Top-T.Ts3p)+$
   nu_ii(32.0,16.0,3.0,2.0,n.ns3p,n.no2p,T.Ts3p,T.To2p,h.s3p,h.o2p)*(T.To2p-T.Ts3p)+$
   nu_ie(32.0,3.0,n.nel,n.ns3p,T.Tel,T.Ts3p,h.el,h.s3p)*(T.Tel-T.Ts3p) + $
   nu_ieh(32.0,3.0,n.nelh,n.ns3p,T.Telh,T.Ts3p,h.s3p)*(T.Telh-T.Ts3p) 
 
EF_s3p = S - L + Teq

return,EF_s3p
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function EF_op,n,T,nT,r,h
;-----------------------------------------------------------------------  
S = (r.io*ftavg_n_el(n.no,n,h.o,h,h.op) + $          ;cold e ionization 
     r.ioh*ftavg_n_elh(n.no,n,h.o,h,h.op) + $         ;hot e ionization
     r.cx_k5*ftavg_n(n.nop,n.no,h.op,h.o,h.op) + $            ;O + O+ -> O+ + O
     r.cx_k6*ftavg_n(n.no2p,n.no,h.o2p,h.o,h.op) + $           ;O + O++ -> O+ + O+
     r.cx_k8*ftavg_n(n.nsp,n.no,h.sp,h.o,h.op) + $            ;O + S+ -> O+ + S
     r.cx_k12*ftavg_n(n.ns2p,n.no,h.s2p,h.o,h.op) + $          ;O + S++ -> O+ + S+ 
     r.cx_k14*ftavg_n(n.ns3p,n.no,h.s3p,h.o,h.op))*T.Tpu_o + $ ;O + S+++ -> O+ + S++
     r.cx_k6*ftavg_n(n.no2p,n.no,h.o2p,h.o,h.op)*T.To2p + $    ;O + O++ -> O+ + O+
     r.cx_k10*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.op)*T.To2p + $   ;S + O++ -> S+ + O+
     r.cx_k11*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.op)*T.To2p + $   ;S + O++ -> S+ + O+ + e
     r.cx_k13*ftavg(n.no2p,n.nsp,h.o2p,h.sp,h.op)*T.To2p + $  ;O++ + S+ -> O+ + S++
     r.cx_k15*ftavg(n.no2p,n.ns2p,h.o2p,h.s2p,h.op)*T.To2p + $ ;O++ + S++ -> O+ + S+++
     r.ro2p*ftavg_el(n.no2p,n,h.o2p,h,h.op)*T.To2p

r.Oei = (r.io*n.no*n.nel + r.ioh*n.no*n.nelh)*T.Tpu_o
r.Oecx = (r.cx_k5*n.no*n.nop + r.cx_k6*n.no*n.no2p + $
         r.cx_k8*n.no*n.nsp + r.cx_k12*n.no*n.ns2p + $
         r.cx_k14*n.no*n.ns3p+r.cx_k7*n.no*n.no2p)*T.Tpu_o

L= r.iop*ftavg_el(n.nop,n,h.op,h,h.op)*T.Top + $
   r.ioph*ftavg_elh(n.nop,n,h.op,h,h.op)*T.Top + $
   r.cx_k5*ftavg_n(n.nop,n.no,h.op,h.o,h.op)*T.Top + $     ;O + O+ -> O+ + O
   r.cx_k9*ftavg_n(n.nop,n.ns,h.op,h.s,h.op)*T.Top + $     ;S + O+ -> S+ + O
   r.Transport*n.nop*T.Top + $
   r.rop*ftavg_el(n.nop,n,h.op,h,h.op)*T.Top         

Teq= nu_ii(16.0,16.0,1.0,2.0,n.nop,n.no2p,T.Top,T.To2p,h.op,h.o2p)*(T.To2p-T.Top)+$
   nu_ii(16.0,32.0,1.0,1.0,n.nop,n.nsp,T.Top,T.Tsp,h.op,h.sp)*(T.Tsp-T.Top)+$
   nu_ii(16.0,32.0,1.0,2.0,n.nop,n.ns2p,T.Top,T.Ts2p,h.op,h.s2p)*(T.Ts2p-T.Top)+$
   nu_ii(16.0,32.0,1.0,3.0,n.nop,n.ns3p,T.Top,T.Ts3p,h.op,h.s3p)*(T.Ts3p-T.Top)+$
   nu_ie(16.0,1.0,n.nel,n.nop,T.Tel,T.Top,h.el,h.op)*(T.Tel-T.Top) + $
   nu_ieh(16.0,1.0,n.nelh,n.nop,T.Telh,T.Top,h.op)*(T.Telh-T.Top) 
 
EF_op = S - L + Teq

return,EF_op
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function EF_o2p,n,T,nT,r,h
;-----------------------------------------------------------------------  
S = (r.iop*ftavg_el(n.nop,n,h.op,h,h.o2p) + $                  ;cold e ionization 
     r.ioph*ftavg_elh(n.nop,n,h.op,h,h.o2p))*T.Top + $          ;hot e ionization
     r.cx_k7*ftavg_n(n.no2p,n.no,h.o2p,h.o,h.o2p)*T.Tpu_o             ;O + O++ -> O++ + O

L= r.cx_k6*ftavg_n(n.no2p,n.no,h.o2p,h.o,h.o2p)*T.To2p + $     ;O + O++ -> O+ + O+
   r.cx_k7*ftavg_n(n.no2p,n.no,h.o2p,h.o,h.o2p)*T.To2p + $     ;O + O++ -> O++ + O
   r.cx_k10*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.o2p)*T.To2p + $    ;S + O++ -> S+ + O+
   r.cx_k11*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.o2p)*T.To2p + $    ;S + O++ -> S++ + O+ + e
   r.cx_k13*ftavg(n.no2p,n.nsp,h.o2p,h.sp,h.o2p)*T.To2p + $   ;O++ + S+ -> O+ + S++
   r.cx_k15*ftavg(n.no2p,n.ns2p,h.o2p,h.s2p,h.o2p)*T.To2p + $  ;O++ + S++ -> O+ + S+++
   r.Transport*n.no2p*T.To2p + $
   r.ro2p*ftavg_el(n.no2p,n,h.o2p,h,h.o2p)*T.To2p   

Teq= nu_ii(16.0,16.0,2.0,1.0,n.no2p,n.nop,T.To2p,T.Top,h.o2p,h.op)*(T.Top-T.To2p)+$
   nu_ii(16.0,32.0,2.0,1.0,n.no2p,n.nsp,T.To2p,T.Tsp,h.o2p,h.sp)*(T.Tsp-T.To2p)+$
   nu_ii(16.0,32.0,2.0,2.0,n.no2p,n.ns2p,T.To2p,T.Ts2p,h.o2p,h.s2p)*(T.Ts2p-T.To2p)+$
   nu_ii(16.0,32.0,2.0,3.0,n.no2p,n.ns3p,T.To2p,T.Ts3p,h.o2p,h.s3p)*(T.Ts3p-T.To2p)+$
   nu_ie(16.0,2.0,n.nel,n.no2p,T.Tel,T.To2p,h.el,h.o2p)*(T.Tel-T.To2p) + $
   nu_ieh(16.0,2.0,n.nelh,n.no2p,T.Telh,T.To2p,h.o2p)*(T.Telh-T.To2p) 
 
EF_o2p = S - L + Teq

return,EF_o2p
end
;-----------------------------------------------------------------------  

;-----------------------------------------------------------------------  
pro get_neutral_loss_rate,n,r,T,nl,h
;-----------------------------------------------------------------------  

;nl.sei = r.is*n.ns*n.nel
nl.sei = r.is*ftavg_n_el(n.ns,n,h.s,h,h.sp) 
nl.seih =  r.ish*ftavg_n_elh(n.ns,n,h.s,h,h.sp)
;nl.seih = r.ish*n.ns*n.nelh
;nl.oei = r.io*n.no*n.nel
;nl.oeih = r.ioh*n.no*n.nelh
nl.oei = r.io*ftavg_n_el(n.no,n,h.o,h,h.op) 
nl.oeih = r.ioh*ftavg_n_elh(n.no,n,h.o,h,h.op)

;nl.k1 = r.cx_k1*n.ns*n.nsp
nl.k1 = r.cx_k1*ftavg_n(n.nsp,n.ns,h.sp,h.s,h.s)
;nl.k2 = r.cx_k2*n.ns*n.ns2p
nl.k2 = r.cx_k2*ftavg_n(n.ns2p,n.ns,h.s2p,h.s,h.s)
;nl.k3 = r.cx_k3*n.ns*n.ns2p
nl.k3 = r.cx_k3*ftavg_n(n.ns2p,n.ns,h.s2p,h.s,h.s)
;nl.k4 = r.cx_k4*n.ns*n.ns3p
nl.k4 = r.cx_k4*ftavg_n(n.ns3p,n.ns,h.s3p,h.s,h.s)
;nl.k5 = r.cx_k5*n.no*n.nop
nl.k5 = r.cx_k5*ftavg_n(n.nop,n.no,h.op,h.o,h.o)
;nl.k6 = r.cx_k6*n.no*n.no2p
nl.k6 = r.cx_k6*ftavg_n(n.no2p,n.no,h.o2p,h.o,h.o)
;nl.k7 = r.cx_k7*n.no*n.no2p
nl.k7 = r.cx_k7*ftavg_n(n.no2p,n.no,h.o2p,h.o,h.o)
;nl.k8 = r.cx_k8*n.no*n.nsp
nl.k8 = r.cx_k8*ftavg_n(n.nsp,n.no,h.sp,h.o,h.o)
;nl.k9 = r.cx_k9*n.ns*n.nop
nl.k9 = r.cx_k9*ftavg_n(n.nop,n.ns,h.op,h.s,h.s)
;nl.k10 = r.cx_k10*n.ns*n.no2p
nl.k10 = r.cx_k10*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.s)
;nl.k11 = r.cx_k11*n.ns*n.no2p
nl.k11 = r.cx_k11*ftavg_n(n.no2p,n.ns,h.o2p,h.s,h.s)
;nl.k12 = r.cx_k12*n.no*n.ns2p
nl.k12 = r.cx_k12*ftavg_n(n.ns2p,n.no,h.s2p,h.o,h.o)
;nl.k14 = r.cx_k14*n.no*n.ns3p
nl.k14 = r.cx_k14*ftavg_n(n.ns3p,n.no,h.s3p,h.o,h.o)


;nl.fast_S_k1 = r.cx_k1*n.ns*n.nsp*T.Tsp
nl.fast_S_k1 = r.cx_k1*ftavg_n(n.nsp,n.ns,h.sp,h.s,h.s)*T.Tsp
;nl.fast_S_k3 = r.cx_k3*n.ns*n.ns2p*T.Ts2p
nl.fast_S_k3 = r.cx_k3*ftavg_n(n.ns2p,n.ns,h.s2p,h.s,h.s)*T.Ts2p
;nl.fast_O_k5 = r.cx_k5*n.no*n.nop*T.Top
nl.fast_O_k5 = r.cx_k5*ftavg_n(n.nop,n.no,h.op,h.o,h.o)*T.Top
;nl.fast_O_k7 = r.cx_k7*n.no*n.no2p*T.To2p
nl.fast_O_k7 = r.cx_k7*ftavg_n(n.no2p,n.no,h.o2p,h.o,h.o)*T.To2p
;nl.fast_S_k8 = r.cx_k8*n.no*n.nsp*T.Tsp
nl.fast_S_k8 = r.cx_k8*ftavg_n(n.nsp,n.no,h.sp,h.o,h.o)*T.Tsp
;nl.fast_O_k9 = r.cx_k9*n.ns*n.nop*T.Top
nl.fast_O_k9 = r.cx_k9*ftavg_n(n.nop,n.ns,h.op,h.s,h.s)*T.Top

nl.sion_tot = nl.sei+nl.seih
nl.oion_tot = nl.oei+nl.oeih
nl.scx_tot = nl.k1+nl.k2+nl.k3+nl.k4+nl.k9+nl.k10+nl.k11
;nl.scx_tot = nl.k2+nl.k4+nl.k9+nl.k10+nl.k11
;nl.ocx_tot = nl.k6+nl.k8+nl.k12+nl.k14
nl.ocx_tot = nl.k5+nl.k6+nl.k7+nl.k8+nl.k12+nl.k14

nl.s_tot = nl.sion_tot+nl.scx_tot
nl.o_tot = nl.oion_tot+nl.ocx_tot


return
end
;-----------------------------------------------------------------------  


;;-----------------------------------------------------------------------  
;pro get_neutral_loss_rate,n,r,T,nl
;;-----------------------------------------------------------------------  

;nl.sei = r.is*n.ns*n.nel
;nl.seih = r.ish*n.ns*n.nelh
;nl.oei = r.io*n.no*n.nel
;nl.oeih = r.ioh*n.no*n.nelh

;nl.k1 = r.cx_k1*n.ns*n.nsp
;nl.k2 = r.cx_k2*n.ns*n.ns2p
;nl.k3 = r.cx_k3*n.ns*n.ns2p
;nl.k4 = r.cx_k4*n.ns*n.ns3p
;nl.k5 = r.cx_k5*n.no*n.nop
;nl.k6 = r.cx_k6*n.no*n.no2p
;nl.k7 = r.cx_k7*n.no*n.no2p
;nl.k8 = r.cx_k8*n.no*n.nsp
;nl.k9 = r.cx_k9*n.ns*n.nop
;nl.k10 = r.cx_k10*n.ns*n.no2p
;nl.k11 = r.cx_k11*n.ns*n.no2p
;nl.k12 = r.cx_k12*n.no*n.ns2p
;nl.k14 = r.cx_k14*n.no*n.ns3p

;nl.fast_S_k1 = r.cx_k1*n.ns*n.nsp*T.Tsp
;nl.fast_S_k3 = r.cx_k3*n.ns*n.ns2p*T.Ts2p
;nl.fast_O_k5 = r.cx_k5*n.no*n.nop*T.Top
;nl.fast_O_k7 = r.cx_k7*n.no*n.no2p*T.To2p
;nl.fast_S_k8 = r.cx_k8*n.no*n.nsp*T.Tsp
;nl.fast_O_k9 = r.cx_k9*n.ns*n.nop*T.Top

;nl.sion_tot = nl.sei+nl.seih
;nl.oion_tot = nl.oei+nl.oeih
;nl.scx_tot = nl.k1+nl.k2+nl.k3+nl.k4+nl.k9+nl.k10+nl.k11
;;nl.scx_tot = nl.k2+nl.k4+nl.k9+nl.k10+nl.k11
;;nl.ocx_tot = nl.k6+nl.k8+nl.k12+nl.k14
;nl.ocx_tot = nl.k5+nl.k6+nl.k7+nl.k8+nl.k12+nl.k14

;nl.s_tot = nl.sion_tot+nl.scx_tot
;nl.o_tot = nl.oion_tot+nl.ocx_tot

;return
;end
;;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro energy_conservation,n,r,T,nl,Ec,h
;-----------------------------------------------------------------------  

f = 3./2.

Ec.s_ion = f*nl.sei*T.Tpu_s
Ec.s_ion_h = f*nl.seih*T.Tpu_s
Ec.o_ion = f*nl.oei*T.Tpu_o
Ec.o_ion_h = f*nl.oeih*T.Tpu_o
Ec.s_cx = f*nl.scx_tot*T.Tpu_s
Ec.o_cx = f*nl.ocx_tot*T.Tpu_o
Ec.s_tot_in = f*nl.s_tot*T.Tpu_s
Ec.o_tot_in = f*nl.o_tot*T.Tpu_o

;print,' ' 
;print,'Energy conservation (eV cm^-3 s^-1):'
;print,' '
;print,'   S ionization cold........',Ec.s_ion
;print,'   S ionization hot.........',Ec.s_ion_h
;print,'   O ionization cold........',Ec.o_ion
;print,'   O ionization hot.........',Ec.o_ion_h
;print,' '
;print,'   S chex...................',Ec.s_cx
;print,'   O chex...................',Ec.o_cx
;print,' '
;print,'   Total S input............',Ec.s_tot_in
;print,'   Total O input............',Ec.o_tot_in

;Teq_el = nu_ie(32.0,1.0,n.nel,n.nsp,T.Tel,T.Tsp)*n.nsp*(T.Tsp-T.Tel) + $
;      nu_ie(32.0,2.0,n.nel,n.ns2p,T.Tel,T.Ts2p)*n.ns2p*(T.Ts2p-T.Tel) + $
;      nu_ie(32.0,3.0,n.nel,n.ns3p,T.Tel,T.Ts3p)*n.ns3p*(T.Ts3p-T.Tel) + $
;      nu_ie(16.0,1.0,n.nel,n.nop,T.Tel,T.Top)*n.nop*(T.Top-T.Tel) + $
;      nu_ie(16.0,2.0,n.nel,n.no2p,T.Tel,T.To2p)*n.no2p*(T.To2p-T.Tel) + $
;      nu_ee(n.nel,n.nelh,T.Tel,T.Telh)*(T.Telh - T.Tel)

Teq_el = nu_ei(32.0,1.0,n.nel,n.nsp,T.Tel,T.Tsp,h.el,h.sp)*(T.Tsp-T.Tel) + $
      nu_ei(32.0,2.0,n.nel,n.ns2p,T.Tel,T.Ts2p,h.el,h.s2p)*(T.Ts2p-T.Tel) + $
      nu_ei(32.0,3.0,n.nel,n.ns3p,T.Tel,T.Ts3p,h.el,h.s3p)*(T.Ts3p-T.Tel) + $
      nu_ei(16.0,1.0,n.nel,n.nop,T.Tel,T.Top,h.el,h.op)*(T.Top-T.Tel) + $
      nu_ei(16.0,2.0,n.nel,n.no2p,T.Tel,T.To2p,h.el,h.o2p)*(T.To2p-T.Tel) + $
      nu_ee(n,h,T.Tel,T.Telh)*(T.Telh - T.Tel)

Teq_sp = nu_ii(32.0,32.0,1.0,2.0,n.nsp,n.ns2p,T.Tsp,T.Ts2p,h.sp,h.s2p)*(T.Ts2p-T.Tsp)+$
   nu_ii(32.0,32.0,1.0,3.0,n.nsp,n.ns3p,T.Tsp,T.Ts3p,h.sp,h.s3p)*(T.Ts3p-T.Tsp)+$
   nu_ii(32.0,16.0,1.0,1.0,n.nsp,n.nop,T.Tsp,T.Top,h.sp,h.op)*(T.Top-T.Tsp)+$
   nu_ii(32.0,16.0,1.0,2.0,n.nsp,n.no2p,T.Tsp,T.To2p,h.sp,h.o2p)*(T.To2p-T.Tsp)+$
   nu_ie(32.0,1.0,n.nel,n.nsp,T.Tel,T.Tsp,h.el,h.sp)*(T.Tel-T.Tsp) + $
   nu_ieh(32.0,1.0,n.nelh,n.nsp,T.Telh,T.Tsp,h.sp)*(T.Telh-T.Tsp) 

Teq_s2p = nu_ii(32.0,32.0,2.0,1.0,n.ns2p,n.nsp,T.Ts2p,T.Tsp,h.s2p,h.sp)*(T.Tsp-T.Ts2p)+$
   nu_ii(32.0,32.0,2.0,3.0,n.ns2p,n.ns3p,T.Ts2p,T.Ts3p,h.s2p,h.s3p)*(T.Ts3p-T.Ts2p)+$
   nu_ii(32.0,16.0,2.0,1.0,n.ns2p,n.nop,T.Ts2p,T.Top,h.s2p,h.op)*(T.Top-T.Ts2p)+$
   nu_ii(32.0,16.0,2.0,2.0,n.ns2p,n.no2p,T.Ts2p,T.To2p,h.s2p,h.o2p)*(T.To2p-T.Ts2p)+$
   nu_ie(32.0,2.0,n.nel,n.ns2p,T.Tel,T.Ts2p,h.el,h.s2p)*(T.Tel-T.Ts2p) + $
   nu_ieh(32.0,2.0,n.nelh,n.ns2p,T.Telh,T.Ts2p,h.s2p)*(T.Telh-T.Ts2p) 

Teq_s3p = nu_ii(32.0,32.0,3.0,1.0,n.ns3p,n.nsp,T.Ts3p,T.Tsp,h.s3p,h.sp)*(T.Tsp-T.Ts3p)+$
   nu_ii(32.0,32.0,3.0,2.0,n.ns3p,n.ns2p,T.Ts3p,T.Ts2p,h.s3p,h.s2p)*(T.Ts2p-T.Ts3p)+$
   nu_ii(32.0,16.0,3.0,1.0,n.ns3p,n.nop,T.Ts3p,T.Top,h.s3p,h.op)*(T.Top-T.Ts3p)+$
   nu_ii(32.0,16.0,3.0,2.0,n.ns3p,n.no2p,T.Ts3p,T.To2p,h.s3p,h.o2p)*(T.To2p-T.Ts3p)+$
   nu_ie(32.0,3.0,n.nel,n.ns3p,T.Tel,T.Ts3p,h.el,h.s3p)*(T.Tel-T.Ts3p) + $
   nu_ieh(32.0,3.0,n.nelh,n.ns3p,T.Telh,T.Ts3p,h.s3p)*(T.Telh-T.Ts3p) 

Teq_op = nu_ii(16.0,16.0,1.0,2.0,n.nop,n.no2p,T.Top,T.To2p,h.op,h.o2p)*(T.To2p-T.Top)+$
   nu_ii(16.0,32.0,1.0,1.0,n.nop,n.nsp,T.Top,T.Tsp,h.op,h.sp)*(T.Tsp-T.Top)+$
   nu_ii(16.0,32.0,1.0,2.0,n.nop,n.ns2p,T.Top,T.Ts2p,h.op,h.s2p)*(T.Ts2p-T.Top)+$
   nu_ii(16.0,32.0,1.0,3.0,n.nop,n.ns3p,T.Top,T.Ts3p,h.op,h.s3p)*(T.Ts3p-T.Top)+$
   nu_ie(16.0,1.0,n.nel,n.nop,T.Tel,T.Top,h.el,h.op)*(T.Tel-T.Top) + $
   nu_ieh(16.0,1.0,n.nelh,n.nop,T.Telh,T.Top,h.op)*(T.Telh-T.Top) 

Teq_o2p = nu_ii(16.0,16.0,2.0,1.0,n.no2p,n.nop,T.To2p,T.Top,h.o2p,h.op)*(T.Top-T.To2p)+$
   nu_ii(16.0,32.0,2.0,1.0,n.no2p,n.nsp,T.To2p,T.Tsp,h.o2p,h.sp)*(T.Tsp-T.To2p)+$
   nu_ii(16.0,32.0,2.0,2.0,n.no2p,n.ns2p,T.To2p,T.Ts2p,h.o2p,h.s2p)*(T.Ts2p-T.To2p)+$
   nu_ii(16.0,32.0,2.0,3.0,n.no2p,n.ns3p,T.To2p,T.Ts3p,h.o2p,h.s3p)*(T.Ts3p-T.To2p)+$
   nu_ie(16.0,2.0,n.nel,n.no2p,T.Tel,T.To2p,h.el,h.o2p)*(T.Tel-T.To2p) + $
   nu_ieh(16.0,2.0,n.nelh,n.no2p,T.Telh,T.To2p,h.o2p)*(T.Telh-T.To2p) 


;Teq_eh = nu_ee(n.nel,n.nelh,T.Tel,T.Telh)*(T.Telh - T.Tel) + $
; nu_ie(32.0,1.0,n.nelh,n.nsp,T.Telh,T.Tsp)*n.nsp*(T.Telh-T.Tsp) + $
; nu_ie(32.0,2.0,n.nelh,n.ns2p,T.Telh,T.Ts2p)*n.ns2p*(T.Telh-T.Ts2p) + $
; nu_ie(32.0,3.0,n.nelh,n.ns3p,T.Telh,T.Ts3p)*n.ns3p*(T.Telh-T.Ts3p) + $
; nu_ie(16.0,1.0,n.nelh,n.nop,T.Telh,T.Top)*n.nop*(T.Telh-T.Top) + $
; nu_ie(16.0,2.0,n.nelh,n.no2p,T.Telh,T.To2p)*n.no2p*(T.Telh-T.To2p) 

Teq_eh = nu_ee(n,h,T.Tel,T.Telh)*(T.Telh - T.Tel) + $
   nu_ieh(32.0,1.0,n.nelh,n.nsp,T.Telh,T.Tsp,h.sp)*(T.Telh-T.Tsp) + $ 
   nu_ieh(32.0,2.0,n.nelh,n.ns2p,T.Telh,T.Ts2p,h.s2p)*(T.Telh-T.Ts2p) + $ 
   nu_ieh(32.0,3.0,n.nelh,n.ns3p,T.Telh,T.Ts3p,h.s3p)*(T.Telh-T.Ts3p) + $ 
   nu_ieh(16.0,1.0,n.nelh,n.nop,T.Telh,T.Top,h.op)*(T.Telh-T.Top) + $
   nu_ieh(16.0,1.0,n.nelh,n.nop,T.Telh,T.Top,h.op)*(T.Telh-T.Top)  

Ec.P_pu = f*nl.s_tot*T.Tpu_s + f*nl.o_tot*T.Tpu_o

Ec.eh_eq = f*Teq_eh 
;Ec.ion_e_eq = f*(nu_ie(32.0,1.0,n.nel,n.nsp,T.Tel,T.Tsp)*n.nsp*(T.Tsp-T.Tel) + $
;      nu_ie(32.0,2.0,n.nel,n.ns2p,T.Tel,T.Ts2p)*n.ns2p*(T.Ts2p-T.Tel) + $
;      nu_ie(32.0,3.0,n.nel,n.ns3p,T.Tel,T.Ts3p)*n.ns3p*(T.Ts3p-T.Tel) + $
;      nu_ie(16.0,1.0,n.nel,n.nop,T.Tel,T.Top)*n.nop*(T.Top-T.Tel) + $
;      nu_ie(16.0,2.0,n.nel,n.no2p,T.Tel,T.To2p)*n.no2p*(T.To2p-T.Tel))

;Ec.ion_e_eq = f*(nu_ei(32.0,1.0,n.nel,n.nsp,T.Tel,T.Tsp)*n.nel*(T.Tsp-T.Tel) + $
;      nu_ei(32.0,2.0,n.nel,n.ns2p,T.Tel,T.Ts2p)*n.nel*(T.Ts2p-T.Tel) + $
;      nu_ei(32.0,3.0,n.nel,n.ns3p,T.Tel,T.Ts3p)*n.nel*(T.Ts3p-T.Tel) + $
;      nu_ei(16.0,1.0,n.nel,n.nop,T.Tel,T.Top)*n.nel*(T.Top-T.Tel) + $
;      nu_ei(16.0,2.0,n.nel,n.no2p,T.Tel,T.To2p)*n.nel*(T.To2p-T.Tel))

;tot_pu = f*nl.s_tot*T.Tpu_s + f*nl.o_tot*T.Tpu_o
;tot_eh = f*Teq_eh
tot_ee = f*nu_ee(n,h,T.Tel,T.Telh)*(T.Telh - T.Tel)

;print,'   Pickup input.............',Ec.P_pu
;print,'   Hot e equilibration......',Ec.eh_eq
;print,'   Hot e cold e equil.......',tot_ee
;;print,'   Ion e equilibration......',Ec.ion_e_eq
Ec.P_in = Ec.P_pu + Ec.eh_eq
;print,'   P input..................',Ec.P_in
;print,' '
Ec.Puv = r.Puv
;print,'   Ion/elec equil...........',Ec.ion_e_eq
;print,'   Hot elec/elec equil......',tot_ee
;print,'   Total elec input.........',Ec.ion_e_eq+tot_ee
;print,' '
;print,'   Puv,Peuv......................',Ec.Puv,r.Peuv
Ec.Pfast = f*(nl.fast_S_k1 + nl.fast_S_k3 + nl.fast_O_k5 + $
           nl.fast_O_k7 + nl.fast_O_k9 + nl.fast_S_k8) 
Ec.Pfast_O = f*(nl.fast_O_k5 + nl.fast_O_k7 + nl.fast_O_k9) 
Ec.Pfast_S = f*(nl.fast_S_k1 + nl.fast_S_k3 + nl.fast_S_k8) 
;print,' '
;print,'   Pfast....................',Ec.Pfast
;print,'   Pfast O..................',Ec.Pfast_O
;print,'   Pfast S..................',Ec.Pfast_S

Ec.Ptrans = f*r.transport*(n.nsp*T.Tsp + n.ns2p*T.Ts2p + $
         n.ns3p*T.Ts3p + n.nop*T.Top + n.no2p*T.To2p + n.nel*T.Tel) 
Ec.Ptrans_O = f*r.transport*(n.nop*T.Top + n.no2p*T.To2p) 
Ec.Ptrans_S = f*r.transport*(n.nsp*T.Tsp + n.ns2p*T.Ts2p + $
         n.ns3p*T.Ts3p) 
Ec.Ptrans_e = f*r.transport*n.nel*T.Tel
Ec.Ptrans_eh = f*r.transport*n.nel*n.fh*T.Telh
;print,' '
;print,'   Ptrans...................',Ec.Ptrans
;print,'   Ptrans O.................',Ec.Ptrans_O
;print,'   Ptrans S.................',Ec.Ptrans_S
;print,'   Ptrans cold elec.........',Ec.Ptrans_e
;print,'   Ptrans hot elec..........',Ec.Ptrans_eh



Ec.P_out = Ec.Puv + Ec.Pfast + Ec.Ptrans + Ec.Ptrans_eh
;print,'   P output............... .',Ec.P_out

;print,' '
;print,'   Ein/Eout.................',Ec.P_in/Ec.P_out


sinput = Ec.s_tot_in/(f*T.Tpu_s)
oinput = Ec.o_tot_in/(f*T.Tpu_o)

totalinput = sinput+oinput

fastout_s = r.cx_k1*n.ns*n.nsp + r.cx_k3*n.ns*n.ns2p +r.cx_k8*n.no*n.nsp
fastout_o = r.cx_k5*n.no*n.nop + r.cx_k7*n.no*n.no2p + $
            r.cx_k9*n.ns*n.nop

transout_s = r.transport*(n.nsp + n.ns2p + n.ns3p)
transout_o = r.transport*(n.nop + n.no2p) 
;transout_e = f*r.transport*n.nel*T.Tel

;print,' '
;print,'Mass flow'
;print,'S input.............',sinput,sinput/totalinput
;print,'O input.............',oinput,oinput/totalinput
;print,'Fast output.........',fastout_s/totalinput,fastout_o/totalinput,(fastout_s+fastout_o)/totalinput
;print,'Trans output........',transout_s/totalinput,transout_o/totalinput,(transout_s+transout_o)/totalinput

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro get_Peuv,n,T,r,h,Lshl,dLshl,rad_tot_v
;-----------------------------------------------------------------------  
nel = n.fc*n.nel

restore,'euv_emissions.sav'

s2em = s2em/1.6e-12
s3em = s3em/1.6e-12
s4em = s4em/1.6e-12
o2em = o2em/1.6e-12
o3em = o3em/1.6e-12

nel = n.nel

dz = h.el/20
nz = 30
z = findgen(nz)*dz

rad_sp = fltarr(nz)
rad_s2p = fltarr(nz)
rad_s3p = fltarr(nz)
rad_op = fltarr(nz)
rad_o2p = fltarr(nz)
rad_tot = fltarr(nz)
deltaT = fltarr(nz)

for k = 0,nz-1 do begin

;print,k,z(k)

nel = n.nel*exp(-z(k)^2 / h.el^2)

if (nel le 1.0) then nel = 1.01

i=where(abs(temp-T.Tel) eq min(abs(temp-T.Tel)))
j=where(abs(density-nel) eq min(abs(density-nel)))

i = i(0)
j = j(0)

if (T.Tel gt temp(i)) then ip = i+1
if (T.Tel le temp(i)) then begin
   ip = i
   i = ip-1
endif

if (nel gt density(j)) then jp = j+1
if (nel le density(j)) then begin
   jp = j
   j = jp-1
endif

dT = temp(ip) - temp(i)
dn = density(jp) - density(j)

A = dT*dn
dT1 = T.Tel - temp(i)
dT2 = temp(ip) - T.Tel
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSII = s2em(i,j)*A4 + s2em(ip,j)*A3 + s2em(i,jp)*A2 + $
          s2em(ip,jp)*A1
L = emisSII*n.nsp
rad_sp(k) = emisSII*n.nsp*exp(-z(k)^2/h.sp^2)

dT = temp(ip) - temp(i)
dn = density(jp) - density(j)

A = dT*dn
dT1 = T.Tel - temp(i)
dT2 = temp(ip) - T.Tel
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSIII = s3em(i,j)*A4 + s3em(ip,j)*A3 + s3em(i,jp)*A2 + $
          s3em(ip,jp)*A1
L = L+ emisSIII*n.ns2p
rad_s2p(k) = emisSIII*n.ns2p*exp(-z(k)^2/h.s2p^2)

dT = temp(ip) - temp(i)
dn = density(jp) - density(j)

A = dT*dn
dT1 = T.Tel - temp(i)
dT2 = temp(ip) - T.Tel
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSIV = s4em(i,j)*A4 + s4em(ip,j)*A3 + s4em(i,jp)*A2 + $
          s4em(ip,jp)*A1
L = L + emisSIV*n.ns3p
rad_s3p(k) = emisSIV*n.ns3p*exp(-z(k)^2/h.s3p^2)

dT = temp(ip) - temp(i)
dn = density(jp) - density(j)

A = dT*dn
dT1 = T.Tel - temp(i)
dT2 = temp(ip) - T.Tel
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisOII = o2em(i,j)*A4 + o2em(ip,j)*A3 + o2em(i,jp)*A2 + $
          o2em(ip,jp)*A1
L = L + emisOII*n.nop
rad_op(k) = emisOII*n.nop*exp(-z(k)^2/h.op^2)

dT = temp(ip) - temp(i)
dn = density(jp) - density(j)

A = dT*dn
dT1 = T.Tel - temp(i)
dT2 = temp(ip) - T.Tel
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisOIII = o3em(i,j)*A4 + o3em(ip,j)*A3 + o3em(i,jp)*A2 + $
          o3em(ip,jp)*A1
L = L + emisOIII*n.no2p
rad_o2p(k) = emisOIII*n.no2p*exp(-z(k)^2/h.o2p^2)

endfor

rad_tot = rad_sp + rad_s2p + rad_s3p + rad_op + rad_o2p

Rj = 7.14e4*1e5 ;cm
dLshl2 = dLshl/2
a = Lshl-dlshl2
b = Lshl+dlshl2
vol = !pi^(1.5)*((h.el/1.5)*1e5)*(b^2-a^2)*Rj^2
dvol = (dz*1e5)*!pi*(b^2-a^2)*Rj^2

;print,Lshl,vol,h.el

;rad_tot_v = 2.0*total(rad_tot(0:*))*dvol
rad_tot_v = (rad_tot(0) + 2.0*total(rad_tot(1:*)))*dvol

;rad_tot_v = 2.0*total(rad_tot)*vol
;rad_tot_v = rad_tot(0)*vol

nel = n.nel*exp(-z^2/h.el^2)

;plot,exp(-z^2/(h.el/1.5)^2)
;oplot,rad_tot/max(rad_tot),linestyle=1

;stop

Lavg = total(rad_tot*nel)/total(nel)

r.Peuv = Lavg
r.Peuv = rad_tot_v

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro get_Peuv_h,n,T,r,h,Lshl,dLshl,rad_tot_v
;-----------------------------------------------------------------------  
nel = n.fh*n.nel

restore,'euv_emissions.sav'

s2em = s2em/1.6e-12
s3em = s3em/1.6e-12
s4em = s4em/1.6e-12
o2em = o2em/1.6e-12
o3em = o3em/1.6e-12

dz = h.el/20
nz = 30
z = findgen(nz)*dz

rad_sp = fltarr(nz)
rad_s2p = fltarr(nz)
rad_s3p = fltarr(nz)
rad_op = fltarr(nz)
rad_o2p = fltarr(nz)
rad_tot = fltarr(nz)
deltaT = fltarr(nz)

for k = 0,nz-1 do begin

;print,k,z(k)

;nel = n.nel*exp(-z(k)^2 / h.el^2)

if (nel le 1.0) then nel = 1.01

i=where(abs(temp-T.Telh) eq min(abs(temp-T.Telh)))
j=where(abs(density-nel) eq min(abs(density-nel)))

i = i(0)
j = j(0)

if (T.Telh gt temp(i)) then ip = i+1
if (T.Telh le temp(i)) then begin
   ip = i
   i = ip-1
endif

if (nel gt density(j)) then jp = j+1
if (nel le density(j)) then begin
   jp = j
   j = jp-1
endif

dT = temp(ip) - temp(i)
dn = density(jp) - density(j)

A = dT*dn
dT1 = T.Telh - temp(i)
dT2 = temp(ip) - T.Telh
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSII = s2em(i,j)*A4 + s2em(ip,j)*A3 + s2em(i,jp)*A2 + $
          s2em(ip,jp)*A1
L = emisSII*n.nsp
rad_sp(k) = emisSII*n.nsp*exp(-z(k)^2/h.sp^2)

dT = temp(ip) - temp(i)
dn = density(jp) - density(j)

A = dT*dn
dT1 = T.Telh - temp(i)
dT2 = temp(ip) - T.Telh
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSIII = s3em(i,j)*A4 + s3em(ip,j)*A3 + s3em(i,jp)*A2 + $
          s3em(ip,jp)*A1
L = L+ emisSIII*n.ns2p
rad_s2p(k) = emisSIII*n.ns2p*exp(-z(k)^2/h.s2p^2)

dT = temp(ip) - temp(i)
dn = density(jp) - density(j)

A = dT*dn
dT1 = T.Telh - temp(i)
dT2 = temp(ip) - T.Telh
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSIV = s4em(i,j)*A4 + s4em(ip,j)*A3 + s4em(i,jp)*A2 + $
          s4em(ip,jp)*A1
L = L + emisSIV*n.ns3p
rad_s3p(k) = emisSIV*n.ns3p*exp(-z(k)^2/h.s3p^2)

dT = temp(ip) - temp(i)
dn = density(jp) - density(j)

A = dT*dn
dT1 = T.Telh - temp(i)
dT2 = temp(ip) - T.Telh
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisOII = o2em(i,j)*A4 + o2em(ip,j)*A3 + o2em(i,jp)*A2 + $
          o2em(ip,jp)*A1
L = L + emisOII*n.nop
rad_op(k) = emisOII*n.nop*exp(-z(k)^2/h.op^2)

dT = temp(ip) - temp(i)
dn = density(jp) - density(j)

A = dT*dn
dT1 = T.Telh - temp(i)
dT2 = temp(ip) - T.Telh
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisOIII = o3em(i,j)*A4 + o3em(ip,j)*A3 + o3em(i,jp)*A2 + $
          o3em(ip,jp)*A1
L = L + emisOIII*n.no2p
rad_o2p(k) = emisOIII*n.no2p*exp(-z(k)^2/h.o2p^2)

endfor

rad_tot = rad_sp + rad_s2p + rad_s3p + rad_op + rad_o2p

Rj = 7.14e4*1e5 ;cm
dLshl2 = dLshl/2
a = Lshl-dlshl2
b = Lshl+dlshl2
vol = !pi^(1.5)*((h.el/1.5)*1e5)*(b^2-a^2)*Rj^2
dvol = (dz*1e5)*!pi*(b^2-a^2)*Rj^2

;print,Lshl,vol,h.el

rad_tot_v = 2.0*total(rad_tot(0:*))*dvol

;rad_tot_v = 2.0*total(rad_tot)*vol
;rad_tot_v = rad_tot(0)*vol

;nel = n.nel*exp(-z^2/h.el^2)

;plot,exp(-z^2/(h.el/1.5)^2)
;oplot,rad_tot/max(rad_tot),linestyle=1

;stop

Lavg = total(rad_tot*nel)/total(nel)

;r.Peuv = Lavg
r.Peuvh = rad_tot_v

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------
;pro cm3_model,max_theta,n,r,T,src,lss,temps,dens,nl,Ec
;pro cm3_model,n,nT,T
;pro cm3_model,n,T,Lshl,dLshl,mr,emis
pro cm3_model,n,T,nT,n1,T1,nT1,np,Tp,nTp,r,r1,nar,n1ar,h,h1,hp,$
              nl,Ec,src,lss,temps,dens,p,Lshl,dLshl,mr,ii

;main program
;cubic centimeter torus chemistry model
;-----------------------------------------------------------------------
@common

mr = fltarr(5)

;Te0 = float(Te0)         ;initial electron temperature (eV)
;dt = dt0          ;time step (seconds)
;ntm = fix(runt/dt)

;fc = 1.0-fh

n_theta = max_theta*2.0 + 1.0

nar = replicate({density_ar, nel: nel0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
          nop: nop0, no2p: no2p0, nelh: nelh0, ns: ns0, no: no0},n_theta)

;n1ar = replicate({density_ar, nel: nel0, nsp: nsp0, ns2p: ns2p0, $
;       ns3p: ns3p0,nop: nop0, no2p: no2p0, nelh: nelh0, ns: ns0, no: no0},n_theta)


h = {scale_heights, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0, s4p:0.0, $
                              o: 0.0, op: 0.0, o2p: 0.0, el: 0.0}

h1 = {scale_heights_1, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0, s4p:0.0, $
                              o: 0.0, op: 0.0, o2p: 0.0, el: 0.0}

;hp = {scale_heights_p, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0, s4p:0.0, $
;                              o: 0.0, op: 0.0, o2p: 0.0, el: 0.0}

;read emission tables
;read_emission_tables,'tepSII',temp,den,emisSII 
;r.emisSII = emis.emisSII
;read_emission_tables,'tepSIII',temp,den,emisSIII 
;r.emisSIII = emis.emisSIII
;read_emission_tables,'tepSIV',temp,den,emisSIV 
;r.emisSIV = emis.emisSIV
;read_emission_tables,'tepOII',temp,den,emisOII 
;r.emisOII = emis.emisOII
;read_emission_tables,'tepOIII',temp,den,emisOIII 
;r.emisOIII = emis.emisOIII
;r.emistemp=emis.temp
;r.emisden=emis.den

;set initial pickup temperatures

T.Tpu_s = Tpu(32.0,Lshl)
T.Tpu_o = Tpu(16.0,Lshl)
;print,T.Tpu_s,T.Tpu_o

T1 = T
Tp = T

r.ish = cfit(16,16,T.Telh,c)
r.isph = cfit(16,15,T.Telh,c)
r.is2ph = cfit(16,14,T.Telh,c)
r.ioh = cfit(8,8,T.Telh,c)
r.ioph = cfit(8,7,T.Telh,c)

;if (plot_density eq 1) then begin
;!p.thick=2.0
;!x.thick=2.0
;!y.thick=2.0
;!p.charsize=1.2
;plot_io,[0,ntm*dt]/8.64e4,[1.0,no0+ns0],/nodata,$
;  title='Density',$
;  xtitle='time (days) ',ytitle='n (cm!u-3!n)',yrange=[0.1,10000.0],/ysty,$
;  xrange=[0,dt*ntm]/8.64e4,/xsty
;endif

;if (plot_temperature eq 1) then begin
;!p.thick=2.0
;!x.thick=2.0
;!y.thick=2.0
;!p.charsize=1.2
;plot_io,[0,ntm*dt]/8.64e4,[0.1,200.0],/nodata,$
;  title='Temperature',$
;  xtitle='time (days) ',ytitle='n (cm!u-3!n)',yrange=[1.0,200.0],/ysty,$
;  xrange=[0,dt*ntm]/8.64e4,/xsty
;endif

;if (strtrim(!d.name,2) ne 'PS') then device,decompose=0
;loadct,27

get_scale_heights,h,T,n
;get_rates,r,T,h,Lshl,dLshl

;if (cont eq 'true') then begin
;   restore,'restart_cm3.sav'
;   fc = 1.0 - fh
;   n.fc = fc
;   n.fh = fh
;   n1.fc = fc
;   n1.fh = fh
;   np.fc = fc
;   np.fh = fh
;;   r.transport = trans
;;   r.otos = otos
;;   T.Tel = Te0
;   T.Telh = Teh0
;   T1 = T
;   Tp = T
;endif

dz =  6.0*7.14e4*0.1*!dtor
z = 2000*dz/5 -findgen(5000)*dz/5.0

;assume input values represent latitude averages
for i = 0,ntm do begin

    ;convert n to equatorial value and expand to get do averaging.

   n.nelh = n.fh*n.nel

   get_rates,r,T,h,Lshl,dLshl


   ;full time step advance
   n1.ns = n.ns + dt*F_s(n,h,r,S,L)

;   src(i).s = S
;   lss(i).s = L
   n1.nsp = n.nsp + dt*F_sp(n,h,r,S,L)

;   src(i).sp = S
;   lss(i).sp = L
   n1.ns2p = n.ns2p + dt*F_s2p(n,h,r,S,L)

;   src(i).s2p = S
;   lss(i).s2p = L
   n1.ns3p = n.ns3p + dt*F_s3p(n,h,r,S,L)
;   src(i).s3p = S
;   lss(i).s3p = L
   n1.no = n.no + dt*F_o(n,h,r,S,L)
;   src(i).o = S
;   lss(i).o = L
   n1.nop = n.nop + dt*F_op(n,h,r,S,L)
;   if (Lshl lt 6.2) then begin
;      print,Lshl,S,L,S*10.*60.*60.,L*10.*60.*60.,n1.nop 
;      print,' '
;   endif
;   src(i).op = S
;   lss(i).op = L 
   n1.no2p = n.no2p + dt*F_o2p(n,h,r,S,L)
;   src(i).o2p = S
;   lss(i).o2p = L
;   src(i).tm = dt*i/8.64e4
;   lss(i).tm = dt*i/8.64e4

   n1.nel = 1.0*n1.nsp+2.0*n1.ns2p+3.0*n1.ns3p+1.0*n1.nop+2.0*n1.no2p

   ;full time step advance for energy
   nT1.sp = nT.sp + dt*EF_sp(n,T,nT,r,h) 
   nT1.s2p = nT.s2p + dt*EF_s2p(n,T,nT,r,h)
   nT1.s3p = nT.s3p + dt*EF_s3p(n,T,nT,r,h)
   nT1.op = nT.op + dt*EF_op(n,T,nT,r,h)
   nT1.o2p = nT.o2p + dt*EF_o2p(n,T,nT,r,h)
   nT1.el = nT.el + dt*EF_el(n,T,nT,r,h)
   update_temp,n1,nT1,T1
;   T1.Tsp = 70.0 + 150.0*exp(-(i-100)^2/20^2)

   get_scale_heights,h1,T1,n1

   ;Improved Euler advance

   n1.nelh = n1.fh*n1.nel
;   run_iterate_navg,n1,T,h1,max_theta
;   get_lat_dist,n1,T1,nar,max_theta
;   neutral_navg,n1,n1ar,max_theta

   np.ns = n.ns + dt*0.5*(F_s(n,h,r)+F_s(n1,h1,r))
   np.nsp = n.nsp + dt*0.5*(F_sp(n,h,r)+F_sp(n1,h1,r))
   np.ns2p = n.ns2p + dt*0.5*(F_s2p(n,h,r)+F_s2p(n1,h1,r))
   np.ns3p = n.ns3p + dt*0.5*(F_s3p(n,h,r)+F_s3p(n1,h1,r))
   np.no = n.no + dt*0.5*(F_o(n,h,r)+F_o(n1,h1,r))
   np.nop = n.nop + dt*0.5*(F_op(n,h,r)+F_op(n1,h1,r))
   np.no2p = n.no2p + dt*0.5*(F_o2p(n,h,r)+F_o2p(n1,h1,r))

   np.nel = 1.0*np.nsp+2.0*np.ns2p+3.0*np.ns3p+1.0*np.nop+2.0*np.no2p

   nTp.sp = nT.sp + dt*0.5*(EF_sp(n1,T1,nT1,r,h1)+ EF_sp(n,T,nT,r,h))
   nTp.s2p = nT.s2p + dt*0.5*(EF_s2p(n1,T1,nT1,r,h1)+ EF_s2p(n,T,nT,r,h))
   nTp.s3p = nT.s3p + dt*0.5*(EF_s3p(n1,T1,nT1,r,h1)+ EF_s3p(n,T,nT,r,h))
   nTp.op = nT.op + dt*0.5*(EF_op(n1,T1,nT1,r,h1)+ EF_op(n,T,nT,r,h))
   nTp.o2p = nT.o2p + dt*0.5*(EF_o2p(n1,T1,nT1,r,h1)+ EF_o2p(n,T,nT,r,h))
   nTp.el = nT.el + dt*0.5*(EF_el(n1,T1,nT1,r,h1)+ EF_el(n,T,nT,r,h))
   update_temp,np,nTp,Tp

   n = np
   nT = nTp  
   T = Tp

   get_Peuv,n,T,r,h,Lshl,dLshl,rad_tot_v
   get_Peuv_h,n,T,r,h,Lshl,dLshl,rad_tot_v
   
   get_scale_heights,h,T,n

;   ft_rad,n,T,r,h

;   temps(i).Tsp = T.Tsp
;   temps(i).Ts2p = T.Ts2p
;   temps(i).Ts3p = T.Ts3p
;   temps(i).Top = T.Top
;   temps(i).To2p = T.To2p
;   temps(i).Telec = T.Tel
;   dens(i).nsp = n.nsp
;   dens(i).ns2p = n.ns2p
;   dens(i).ns3p = n.ns3p
;   dens(i).nop = n.nop
;   dens(i).no2p = n.no2p
;   dens(i).nel = n.nel
;   p(i).Puv = r.Puv
;   p(i).psp = r.psp
;   p(i).ps2p = r.ps2p
;   p(i).ps3p = r.ps3p
;   p(i).pop = r.pop
;   p(i).po2p = r.po2p

;   print,T.Tsp,T.Ts2p,T.Ts3p,T.Top,T.To2p,T.Tel

;update Pickup Densities
  
;   if (plot_density eq 1) then begin
;   ;   plots,i*dt/8.64e4,n.ns,psym=1,color=10
;      plots,i*dt/8.64e4,n.nsp,psym=1,color=10,thick=1.0
;      plots,i*dt/8.64e4,n.ns2p,psym=1,color=75,thick=1.0
;      plots,i*dt/8.64e4,n.ns3p,psym=1,color=50,thick=1.0
;      plots,i*dt/8.64e4,n.nel,psym=1,color=100,thick=1.0
;   ;   plots,i*dt/8.64e4,n.no,psym=1,color=75,thick=1.0
;      plots,i*dt/8.64e4,n.nop,psym=1,color=150,thick=1.0
;      plots,i*dt/8.64e4,n.no2p,psym=1,color=200,thick=1.0
;   endif

;   if (plot_temperature eq 1) then begin
;      plots,i*dt/8.64e4,T.Tsp,psym=1,color=10,thick=1.0
;      plots,i*dt/8.64e4,T.Ts2p,psym=1,color=75,thick=1.0
;      plots,i*dt/8.64e4,T.Ts3p,psym=1,color=50,thick=1.0
;      plots,i*dt/8.64e4,T.Tel,psym=1,color=100,thick=1.0
;      plots,i*dt/8.64e4,T.Top,psym=1,color=150,thick=1.0
;      plots,i*dt/8.64e4,T.To2p,psym=1,color=200,thick=1.0
;   endif
endfor

get_neutral_loss_rate,n,r,T,nl,h

fastout_s = r.cx_k1*ftavg_n(n.nsp,n.ns,h.sp,h.s,h.s) + r.cx_k3*ftavg_n(n.ns2p,n.ns,h.s2p,h.s,h.s) + $
            r.cx_k8*ftavg_n(n.nsp,n.no,h.sp,h.o,h.o)
fastout_o = r.cx_k5*ftavg_n(n.nop,n.no,h.op,h.o,h.o) + r.cx_k7*ftavg_n(n.no2p,n.no,h.o2p,h.o,h.o) + $
            r.cx_k9*ftavg_n(n.nop,n.ns,h.op,h.s,h.s)


Rj = 7.14e4*1e5 ;cm
;a = 6.0*Rj
;b = 7.5*Rj
;vol = 0.25*!pi^2*(a+b)*(b-a)^2
Area = (sqrt(!pi))*!pi*((Lshl + dLshl/2)^2 - (Lshl - dLshl/2)^2)*Rj^2

;transout_s = r.transport*(n.nsp*h.sp*1e5*Area + n.ns2p*h.s2p*1e5*Area + n.ns3p*h.s3p*1e5*Area)
;transout_o = r.transport*(n.nop*h.op*1e5*Area + n.no2p*h.o2p*1e5*Area) 
;transout_e = f*r.transport*n.nel*T.Tel
totalinput = net_source

;print,' '
;print,'Mass flow'
;print,'S input.............',sinput,sinput/totalinput
;print,'O input.............',oinput,oinput/totalinput
;print,'Fast output.........',fastout_s/totalinput,fastout_o/totalinput,(fastout_s+fastout_o)/totalinput
;print,'Trans output........',transout_s/totalinput,transout_o/totalinput,(transout_s+transout_o)/totalinput
;print,'Scale height...,',h.s
;print,'Fast dot............',ii,(fastout_s*h.s*1e5*Area*32.+fastout_o*h.o*1e5*Area*16)*1.67e-27,net_source;*20*1.67e-27
;print,'Mdot................',(transout_s+transout_o)
;print,'Total mass loss rate...',(fastout_s*h.s*1e5*Area+fastout_o*h.o*1e5*Area) + (transout_s+transout_o)



;energy_conservation,n,r,T,nl,Ec,h

;save,filename='src_func.sav',src
;save,filename='lss_func.sav',lss
;save,filename='temps.sav',temps
;save,filename='dens.sav',dens
;save,filename='pwr.sav',p
;save,filename='density.sav',n

;if (lonoff eq 1.0) then begin
;legend,['S!u+!n','S!u++!n','S!u+++!n','O!u+!n','O!u++!n','e'],$
;       linestyle=[0,0,0,0,0,0],colors=[10,75,50,150,200,100],/bottom,$
;       /right
;endif

;upsp = n.nsp
;ups2p = n.ns2p
;ups3p = n.ns3p
;upop = n.nop
;upo2p = n.no2p

;print,' '
;print,'S h,c ionization rate....',nl.seih/(nl.sion_tot),nl.sei/nl.sion_tot
;print,'O h,c ionization rate....',nl.oeih/(nl.oion_tot),nl.oei/nl.oion_tot
;print,'S/O ionization rate......',nl.sion_tot/nl.oion_tot
;print,'S chex rate k1...........',nl.k1/nl.scx_tot
;print,'S chex rate k2...........',nl.k2/nl.scx_tot
;print,'S chex rate k3...........',nl.k3/nl.scx_tot
;print,'S chex rate k4...........',nl.k4/nl.scx_tot
;print,'S chex rate k9...........',nl.k9/nl.scx_tot
;print,'S chex rate k10..........',nl.k10/nl.scx_tot
;print,'S chex rate k11..........',nl.k11/nl.scx_tot
;print,'O chex rate k5...........',nl.k5/nl.ocx_tot
;print,'O chex rate k6...........',nl.k6/nl.ocx_tot
;print,'O chex rate k7...........',nl.k7/nl.ocx_tot
;print,'O chex rate k8...........',nl.k8/nl.ocx_tot
;print,'O chex rate k12..........',nl.k12/nl.ocx_tot
;print,'O chex rate k14..........',nl.k14/nl.ocx_tot
;print,'S ionization/cx..........',nl.sion_tot/nl.scx_tot
;print,'O ionization/cx..........',nl.oion_tot/nl.ocx_tot
;print,'S/O total loss rate......',nl.s_tot/nl.o_tot
;;print,'total sion,scx...',nl.sion_tot,nl.scx_tot
;;print,'total oion,ocx...',nl.oion_tot,nl.ocx_tot
;print,'ns,no............',n.ns,n.no
;print,'Production of high speed O...',nl.k5
;print,'Production of neutral O......',r.S_production
;print,'High speed O/low speed O.....',nl.k5/r.S_production
;print,' '
nelec = n.nel + 0.1*n.nel
;print,'nsp....',n.nsp/nelec,nar(20).nsp/nar(20).nel     ;0.12, 0.092
;print,'ns2p...',n.ns2p/nelec,nar(20).ns2p/nar(20).nel   ;0.15, 0.17
;print,'ns3p...',n.ns3p/nelec,nar(20).ns3p/nar(20).nel   ;0, 0.023
;print,'nop....',n.nop/nelec,nar(20).nop/nar(20).nel     ;0.36, 0.47
;print,'no2p...',n.no2p/nelec,nar(20).no2p/nar(20).nel   ;0.12, 0.016
;print,' '

dz = h.el/10
nz = 20
zz = findgen(nz)*dz
nel = nelec*exp(-zz^2/h.el^2)
;print,'Flux tube averaged mixing ratios...ion weighted'

nel_tot = n.nsp*sqrt(!pi)*h.sp + 2.0*n.ns2p*sqrt(!pi)*h.s2p + 3.0*n.ns3p*sqrt(!pi)*h.s3p + $
          n.nop*sqrt(!pi)*h.op + 2.0*n.no2p*sqrt(!pi)*h.o2p + 0.1*n.nel*sqrt(!pi)*h.el

;print,'nsp....',total((n.nsp*exp(-zz^2/h.sp^2)/nel)*n.nsp*exp(-zz^2/h.sp^2))/total(n.nsp*exp(-zz^2/h.sp^2)),n.nsp*sqrt(!pi)*h.sp/nel_tot
;print,'ns2p....',total((n.ns2p*exp(-zz^2/h.s2p^2)/nel)*n.ns2p*exp(-zz^2/h.s2p^2))/total(n.ns2p*exp(-zz^2/h.s2p^2)),n.ns2p*sqrt(!pi)*h.s2p/nel_tot
;print,'ns3p....',total((n.ns3p*exp(-zz^2/h.s3p^2)/nel)*n.ns3p*exp(-zz^2/h.s3p^2))/total(n.ns3p*exp(-zz^2/h.s3p^2)),n.ns3p*sqrt(!pi)*h.s3p/nel_tot
;print,'nop....',total((n.nop*exp(-zz^2/h.op^2)/nel)*n.nop*exp(-zz^2/h.op^2))/total(n.nop*exp(-zz^2/h.op^2)),n.nop*sqrt(!pi)*h.op/nel_tot
;print,'no2p....',total((n.no2p*exp(-zz^2/h.o2p^2)/nel)*n.no2p*exp(-zz^2/h.o2p^2))/total(n.no2p*exp(-zz^2/h.o2p^2)),n.no2p*sqrt(!pi)*h.o2p/nel_tot

;print,'Flux tube averaged mixing ratios, electron weighted...'
;print,'nsp....',total((n.nsp*exp(-zz^2/h.sp^2)/nel)*nelec*exp(-zz^2/h.el^2))/total(nelec*exp(-zz^2/h.el^2))
;print,'ns2p....',total((n.ns2p*exp(-zz^2/h.s2p^2)/nel)*nelec*exp(-zz^2/h.el^2))/total(nelec*exp(-zz^2/h.el^2))
;print,'ns3p....',total((n.ns3p*exp(-zz^2/h.s3p^2)/nel)*nelec*exp(-zz^2/h.el^2))/total(nelec*exp(-zz^2/h.el^2))
;print,'nop....',total((n.nop*exp(-zz^2/h.op^2)/nel)*nelec*exp(-zz^2/h.el^2))/total(nelec*exp(-zz^2/h.el^2))
;print,'no2p....',total((n.no2p*exp(-zz^2/h.o2p^2)/nel)*nelec*exp(-zz^2/h.el^2))/total(nelec*exp(-zz^2/h.el^2))

mr(0) = n.nsp*sqrt(!pi)*h.sp/nel_tot
mr(1) = n.ns2p*sqrt(!pi)*h.s2p/nel_tot
mr(2) = n.ns3p*sqrt(!pi)*h.s3p/nel_tot
mr(3) = n.nop*sqrt(!pi)*h.op/nel_tot
mr(4) = n.no2p*sqrt(!pi)*h.o2p/nel_tot
;peuv = r.Peuv

;print,' '
;print,'net_source...',r.O_production+r.S_production
;print,' '
;print,'S/O:     ',n.ns/n.no,nar(20).ns/nar(20).no
;print,'S+/O+:   ',n.nsp/n.nop,nar(20).nsp/nar(20).nop
;print,'S+/S++:  ',n.nsp/n.ns2p,nar(20).nsp/nar(20).ns2p
;print,'S+++/S++: ',n.ns3p/n.ns2p,nar(20).ns3p/nar(20).ns2p
;print,'O++/O+:  ',n.no2p/n.nop,nar(20).no2p/nar(20).nop
;print,' '
;print,'nel.....',n.nel
;print,' '
;neni = n.nel/(n.nsp+n.ns2p+n.ns3p+n.nop+n.no2p)
;sonp = n.nop+n.no2p
;ssnp = n.nsp+n.ns2p+n.ns3p
;print,'ne/ni....',neni
;print,'On+/Sn+..',sonp/ssnp
;print,' '
;print,'Tsp......',T.Tsp
;print,'Ts2p.....',T.Ts2p
;print,'Ts3p.....',T.Ts3p
;print,'Top......',T.Top
;print,'To2p.....',T.To2p
;print,'Tel......',T.Tel
;sti = T.Tsp*n.nsp+T.Ts2p*n.ns2p+T.Ts3p*n.ns3p+T.Top*n.nop+T.To2p*n.no2p
;sni = n.nsp+n.ns2p+n.ns3p+n.nop+n.no2p
;print,'<Ti>.....',sti/sni

;spne1 = n.nsp/n.nel
;s2pne1 = n.ns2p/n.nel
;s3pne1 = n.ns3p/n.nel
;opne1 = n.nop/n.nel
;o2pne1 = n.no2p/n.nel
;neni1 = neni
;onsn1 = sonp/ssnp
;Tiave1 = sti/sni
;Tel1 = T.Tel
;nel1 = n.nel
;save,filename='sensitivity.sav',spne1,s2pne1,s3pne1,opne1,o2pne1,$
;  neni1,onsn1,Tel1,Tiave1,nel1

;save,filename='restart_cm3.sav',n,n1,np,T,T1,Tp,nT,nT1,nTp,Ec,h

;theta = (max_theta-findgen(2.0*max_theta + 1.0))*!dtor      
;set_plot,'ps
;device,filename='neutral.ps'
;!p.thick=2.0
;!x.thick=2.0
;!y.thick=2.0
;!p.font=0
;device,/landscape
;plot_io,theta*!radeg,nar.no,yrange=[1,2000],/ysty,xtitle='lat (deg)',ytitle='density (cm!u-3!n)'
;oplot,theta*!radeg,nar.ns
;oplot,theta*!radeg,nar.nsp,linestyle=1
;oplot,theta*!radeg,nar.ns2p,linestyle=2
;oplot,theta*!radeg,nar.ns3p,linestyle=3                                                         
;oplot,theta*!radeg,nar.nop,linestyle=4                                          
;oplot,theta*!radeg,nar.no2p,linestyle=5                                                     
;device,/close
;set_plot,'x                 
                            

;cm3_ekappa,n,T


return
end
;-----------------------------------------------------------------------




