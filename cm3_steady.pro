;-----------------------------------------------------------------------
pro get_rates,r,T
;-----------------------------------------------------------------------
@common

;Electron impact ionization (Schreier and Eviatar 1998)
;r.is = 7.4e-8*sqrt(T.Tel_el)*exp(-10.36/T.Tel_el)     ;check this!!! 
;r.isp = 1.274e-8*sqrt(T.Tel_sp)*exp(-23.1/T.Tel_sp)  
r.is2p = 5.8e-9*sqrt(T.Tel_s2p)*exp(-34.88/T.Tel_s2p)  

r.is = cfit(16,16,T.Tel_el,c)
r.isp = cfit(16,15,T.Tel_sp,c)
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

lnt = alog(T.Tel_s3p)
A = 1.0
y = -23.9+A*cos((lnt*(1.3*!pi)/(alog(10)) + 0.74*!pi))
r.rs3p = exp(y)

lnt = alog(T.Tel_s2p)
A = 1.4
y = -25.2+A*cos((lnt*(1.3*!pi)/(alog(10)) + 0.8*!pi))
r.rs2p = exp(y)

lnt = alog(T.Tel_o2p)
A = 1.1
y = -26.0+A*cos((lnt*(1.2*!pi)/(alog(10)) + 0.74*!pi))
r.ro2p = exp(y)

;dielectron recombination (Smith and Strobel 1985)
;if (T.Tel_s3p le 6.9) then r.rs3p = 5.4e-11*T.Tel_s3p^0.53
;if (T.Tel_s3p gt 6.9) then r.rs3p = 4.3e-12*T.Tel_s3p^1.84  

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
r.cx_k14 = k14*1.92e-8
;r.cx_k15 = k15*9e-9    ;Schreier for L=6
r.cx_k15 = k15*9e-10    ;McGrath for L=6
r.cx_k16 = k16*3.6e-10  ;McGrath

;print,k4,k10,k11,k13,k14,k15

Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
;vol = 2.5e+31
;print,'vol...',vol

fs = 1.0/(1.0+otos)
fo = otos/(1.0+otos)
;print,fs,fo

r.S_production = fs*net_source/vol
r.O_production = fo*net_source/vol

return
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------
pro get_density_residuals,Rs,n,Rs_arr,r
;-----------------------------------------------------------------------
Rs.s = r.S_production - $
      (r.is*n.ns*n.fc*n.nel + r.ish*n.ns*n.fh*n.nel + $ 
       r.cx_k4*n.ns*n.ns3p + r.cx_k11*n.ns*n.no2p + $
       r.cx_k10*n.ns*n.no2p + $
       r.cx_k2*n.ns*n.ns2p + r.cx_k9*n.ns*n.nop + $
       r.cx_k3*n.ns*n.ns2p + r.cx_k1*n.ns*n.nsp) 

Rs.sp = (r.is*n.ns*n.fc*n.nel + r.ish*n.ns*n.fh*n.nel + $
        r.rs2p*n.ns2p*n.nel + r.cx_k4*n.ns*n.ns3p + $
        r.cx_k10*n.ns*n.no2p + 2.0*r.cx_k2*n.ns*n.ns2p + $
        r.cx_k9*n.ns*n.nop + r.cx_k12*n.no*n.ns2p) - $
        (r.isp*n.nsp*n.fc*n.nel + r.cx_k13*n.no2p*n.nsp + $
         n.nsp*r.Transport + r.isph*n.nsp*n.fh*n.nel + $
         r.cx_k8*n.no*n.nsp + r.cx_k16*n.ns3p*n.nsp)

Rs.s2p = (r.isp*n.nsp*n.fc*n.nel + r.isph*n.nsp*n.fh*n.nel + $
    r.rs3p*n.ns3p*n.nel + r.cx_k13*n.no2p*n.nsp + $
    r.cx_k14*n.no*n.ns3p + r.cx_k4*n.ns*n.ns3p + $
    r.cx_k11*n.ns*n.no2p + 2.0*r.cx_k16*n.ns3p*n.nsp) - $
    (r.is2p*n.ns2p*n.fc*n.nel + $
    r.is2ph*n.ns2p*n.fh*n.nel +  $
    r.rs2p*n.ns2p*n.nel + $
    n.ns2p*r.Transport + $
    r.cx_k2*n.ns*n.ns2p + $        ;S + S++ -> S+ + S+
    r.cx_k12*n.no*n.ns2p + $       ;O + S++ -> O+ + S+ 
    r.cx_k15*n.no2p*n.ns2p)         ;O++ + S++ -> O+ + S+++

Rs.s3p = (r.is2p*n.ns2p*n.fc*n.nel + r.is2ph*n.ns2p*n.fh*n.nel + $
    r.cx_k15*n.no2p*n.ns2p) - $
    (r.rs3p*n.ns3p*n.nel + n.ns3p*r.Transport + r.cx_k14*n.no*n.ns3p + $
    r.cx_k4*n.ns*n.ns3p + r.cx_k16*n.ns3p*n.nsp)

Rs.o = r.O_production - $
    (r.io*n.no*n.fc*n.nel + r.ioh*n.no*n.fh*n.nel + r.cx_k14*n.no*n.ns3p + $
    r.cx_k6*n.no*n.no2p + r.cx_k8*n.no*n.nsp + r.cx_k12*n.no*n.ns2p + $
    r.cx_k5*n.no*n.nop + r.cx_k7*n.no*n.no2p)

Rs.op = (r.io*n.no*n.fc*n.nel + r.ioh*n.no*n.fh*n.nel + $
    r.ro2p*n.no2p*n.nel + r.cx_k13*n.no2p*n.nsp + $
    r.cx_k15*n.no2p*n.ns2p + r.cx_k14*n.no*n.ns3p + $
    r.cx_k11*n.ns*n.no2p + r.cx_k10*n.ns*n.no2p + $
    r.cx_k8*n.no*n.nsp + r.cx_k12*n.no*n.ns2p + 2.0*r.cx_k6*n.no*n.no2p)-$
    (r.iop*n.nop*n.fc*n.nel + n.nop*r.Transport + r.ioph*n.nop*n.fh*n.nel + $
    r.cx_k9*n.ns*n.nop)

Rs.o2p = (r.iop*n.nop*n.fc*n.nel + r.ioph*n.nop*n.fh*n.nel) - $
    (r.ro2p*n.no2p*n.nel + r.cx_k13*n.no2p*n.nsp + $
    r.cx_k15*n.no2p*n.ns2p + n.no2p*r.Transport + $
    r.cx_k11*n.ns*n.no2p + r.cx_k10*n.ns*n.no2p + $
    r.cx_k6*n.no*n.no2p)

Rs_arr(0) = Rs.s
Rs_arr(1) = Rs.sp
Rs_arr(2) = Rs.s2p
Rs_arr(3) = Rs.s3p
Rs_arr(4) = Rs.o
Rs_arr(5) = Rs.op
Rs_arr(6) = Rs.o2p

return
end
;-----------------------------------------------------------------------


;-----------------------------------------------------------------------  
pro get_density_jacobian,Rs,Jc,r,n
;-----------------------------------------------------------------------  

Jc(0,0) = -(r.is*n.fc*n.nel + r.ish*n.fh*n.nel + $ 
       r.cx_k4*n.ns3p + r.cx_k11*n.no2p + $
       r.cx_k10*n.no2p + $
       r.cx_k2*n.ns2p + r.cx_k9*n.nop + $
       r.cx_k3*n.ns2p + r.cx_k1*n.nsp) 
Jc(1,0) = -(r.cx_k1*n.ns) 
Jc(2,0) = -(r.cx_k2*n.ns + r.cx_k3*n.ns) 
Jc(3,0) = -(r.cx_k4*n.ns)
Jc(4,0) = 0.0
Jc(5,0) = -(r.cx_k9*n.ns)
Jc(6,0) = -(r.cx_k11*n.ns + r.cx_k10*n.ns)

Jc(0,1) = (r.is*n.fc*n.nel + r.ish*n.fh*n.nel + $
          r.cx_k4*n.ns3p + $
          r.cx_k10*n.no2p + 2.0*r.cx_k2*n.ns2p + $
          r.cx_k9*n.nop)
Jc(1,1) = -(r.isp*n.fc*n.nel + r.cx_k13*n.no2p + $
            r.Transport + r.isph*n.fh*n.nel + $
            r.cx_k8*n.no + r.cx_k16*n.ns3p)
Jc(2,1) = (r.rs2p*n.nel + 2.0*r.cx_k2*n.ns + r.cx_k12*n.no) 
Jc(3,1) = r.cx_k4*n.ns - (r.cx_k16*n.nsp)
Jc(4,1) = r.cx_k12*n.ns2p -  (r.cx_k8*n.nsp)
Jc(5,1) = r.cx_k9*n.ns
Jc(6,1) = r.cx_k10*n.ns - (r.cx_k13*n.ns)

Jc(0,2) = r.cx_k4*n.ns3p + r.cx_k11*n.no2p - $
          (r.cx_k2*n.ns2p)
Jc(1,2) = r.isp*n.fc*n.nel + r.isph*n.fh*n.nel + $
          r.cx_k13*n.no2p + 2.0*r.cx_k16*n.ns3p
Jc(2,2) = -(r.is2p*n.fc*n.nel + $
    r.is2ph*n.fh*n.nel +  $
    r.rs2p*n.nel + $
    r.Transport + $
    r.cx_k2*n.ns + $        ;S + S++ -> S+ + S+
    r.cx_k12*n.no + $       ;O + S++ -> O+ + S+ 
    r.cx_k15*n.no2p)
Jc(3,2) =  r.rs3p*n.nel + r.cx_k14*n.no + $
           r.cx_k4*n.ns + 2.0*r.cx_k16*n.nsp 
Jc(4,2) = r.cx_k14*n.ns3p - r.cx_k12*n.ns2p
Jc(5,2) = 0.0
Jc(6,2) = r.cx_k13*n.nsp + r.cx_k11*n.ns - (r.cx_k15*n.ns2p)

Jc(0,3) = -r.cx_k4*n.ns3p
Jc(1,3) = -r.cx_k16*n.ns3p
Jc(2,3) = (r.is2p*n.fc*n.nel + r.is2ph*n.fh*n.nel + $
           r.cx_k15*n.no2p) 
Jc(3,3) = -(r.rs3p*n.nel + r.Transport + r.cx_k14*n.no + $
            r.cx_k4*n.ns + r.cx_k16*n.nsp)
Jc(4,3) = -r.cx_k14*n.ns3p
Jc(5,3) = 0.0
Jc(6,3) =  r.cx_k15*n.ns2p

Jc(0,4) = 0.0
Jc(1,4) = -r.cx_k8*n.no
Jc(2,4) = -r.cx_k12*n.no
Jc(3,4) = -r.cx_k14*n.no
Jc(4,4) = -(r.io*n.fc*n.nel + r.ioh*n.fh*n.nel + r.cx_k14*n.ns3p + $
            r.cx_k6*n.no2p + r.cx_k8*n.nsp + r.cx_k12*n.ns2p + $
            r.cx_k5*n.nop + r.cx_k7*n.no2p)
Jc(5,4) = -r.cx_k5*n.no
Jc(6,4) = -(r.cx_k6*n.no + r.cx_k7*n.no)

Jc(0,5) = r.cx_k11*n.no2p + r.cx_k10*n.no2p - r.cx_k9*n.nop
Jc(1,5) = r.cx_k13*n.no2p + r.cx_k8*n.no
Jc(2,5) = r.cx_k15*n.no2p + r.cx_k12*n.no 
Jc(3,5) = r.cx_k14*n.no
Jc(4,5) = (r.io*n.fc*n.nel + r.ioh*n.fh*n.nel + $
  r.cx_k14*n.ns3p + r.cx_k8*n.nsp + r.cx_k12*n.ns2p + 2.0*r.cx_k6*n.no2p)
Jc(5,5) = -(r.iop*n.fc*n.nel + r.Transport + r.ioph*n.fh*n.nel + $
            r.cx_k9*n.ns)
Jc(6,5) = (r.ro2p*n.nel + r.cx_k13*n.nsp + $
           r.cx_k15*n.ns2p + r.cx_k11*n.ns + r.cx_k10*n.ns + $
           2.0*r.cx_k6*n.no)

Jc(0,6) = -(r.cx_k11*n.no2p + r.cx_k10*n.no2p)

Jc(1,6) = -r.cx_k13*n.no2p

Jc(2,6) = -r.cx_k15*n.no2p

Jc(3,6) = 0.0

Jc(4,6) = -r.cx_k6*n.no2p

Jc(5,6) = (r.iop*n.fc*n.nel + r.ioph*n.fh*n.nel) 

Jc(6,6) = -(r.ro2p*n.nel + r.cx_k13*n.nsp + $
    r.cx_k15*n.ns2p + r.Transport + $
    r.cx_k11*n.ns + r.cx_k10*n.ns + $
    r.cx_k6*n.no)

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro get_delta_n,Rs_arr,Jc,delta_n
;-----------------------------------------------------------------------  
;dt = 6e6
;Id = fltarr(7,7)
;for i = 0,6 do begin
;  Id(i,i) = 1.0/dt
;endfor

;A = Id + Jc
;print,' '
;print,Jc
A = Jc
B = -Rs_arr

;Ax = B
;Jc delta_n = - Rs_arr

LUDC, A, INDEX,/double
delta_n = LUSOL(A, INDEX, B,/double)

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro update_density,n,delta_n
;-----------------------------------------------------------------------  

n.ns = n.ns + delta_n(0)
n.nsp = n.nsp + delta_n(1)
n.ns2p = n.ns2p + delta_n(2)
n.ns3p = n.ns3p + delta_n(3)
n.no = n.no + delta_n(4)
n.nop = n.nop + delta_n(5)
n.no2p = n.no2p + delta_n(6)

n.nel = 1.0*n.nsp+2.0*n.ns2p+3.0*n.ns3p+1.0*n.nop+2.0*n.no2p

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
    r.cx_k2*n.ns*n.ns2p + r.cx_k9*n.ns*n.nop + $
    r.cx_k3*n.ns*n.ns2p + r.cx_k1*n.ns*n.nsp 

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
    r.isph*n.nsp*n.fh*n.nel + r.cx_k8*n.no*n.nsp + r.cx_k16*n.ns3p*n.nsp

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
    r.cx_k11*n.ns*n.no2p + 2.0*r.cx_k16*n.ns3p*n.nsp

L = r.is2p*n.ns2p*n.fc*n.nel + $
    r.is2ph*n.ns2p*n.fh*n.nel +  $
    r.rs2p*n.ns2p*n.nel + $
    n.ns2p*r.Transport + $
    r.cx_k2*n.ns*n.ns2p + $        ;S + S++ -> S+ + S+
    r.cx_k12*n.no*n.ns2p + $       ;O + S++ -> O+ + S+ 
    r.cx_k15*n.no2p*n.ns2p         ;O++ + S++ -> O+ + S+++

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
    r.cx_k4*n.ns*n.ns3p + r.cx_k16*n.ns3p*n.nsp

F_s3p = S - L

return,F_s3p
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------
function F_o,n,r,S,L
;-----------------------------------------------------------------------
S = r.O_production ;+ r.cx_k9*n.ns*n.nop
L = r.io*n.no*n.fc*n.nel + r.ioh*n.no*n.fh*n.nel + r.cx_k14*n.no*n.ns3p + $
    r.cx_k6*n.no*n.no2p + r.cx_k8*n.no*n.nsp + r.cx_k12*n.no*n.ns2p + $
    r.cx_k5*n.no*n.nop + r.cx_k7*n.no*n.no2p

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
    r.cx_k11*n.ns*n.no2p + r.cx_k10*n.ns*n.no2p + $
    r.cx_k8*n.no*n.nsp + r.cx_k12*n.no*n.ns2p + 2.0*r.cx_k6*n.no*n.no2p
L = r.iop*n.nop*n.fc*n.nel + n.nop*r.Transport + r.ioph*n.nop*n.fh*n.nel + $
    r.cx_k9*n.ns*n.nop

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
    r.cx_k11*n.ns*n.no2p + r.cx_k10*n.ns*n.no2p + $
    r.cx_k6*n.no*n.no2p

F_o2p = S - L

return,F_o2p
end
;-----------------------------------------------------------------------  

;=======================Energy source/loss functions=====================

;-----------------------------------------------------------------------  
pro get_temp_residuals,RsT,nT,RsT_arr,r
;-----------------------------------------------------------------------  



return
end
;-----------------------------------------------------------------------  

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

nu_ie = 1.8e-19*sqrt(me*mi)*(Z^2)*nelec*ni*lambda/(mi*Tel + me*Ti)^(3./2.)

;print,'nu_ei...',nu_ie,lambda

return,nu_ie
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function d_nu_ie_dTe,mu,z,nelec,ni,Tel,Ti
; m is atomic mass
; z is atomic number
; for collision of given ion on thermal electron population
;-----------------------------------------------------------------------  
me = double(9.1e-28)    ;g
mp = double(1.67e-24)   ;g

mi = mu*mp

lambda = 0.0
dlambda = 0.0

if ((Ti*me/mi lt Tel) and (Tel le 10*z^2)) then begin
   lambda = 23.0 - alog(sqrt(nelec)*z*Tel^(-3./2.))
   dlambda = 1.5/Tel
endif

if ((Tel gt Ti*me/mi) and (Tel ge 10*Z^2)) then begin
   lambda = 24.0 - alog(sqrt(nelec)/Tel)
   dlambda = 1.0/Tel
endif

if (Tel lt Ti*Z*me/mi) then begin
   lambda = 30.0 - alog(sqrt(ni)*Ti^(-3./2.)*Z^2/mu)
   dlambda = 0.0
endif

;nu_ie = 1.8e-19*sqrt(me*mi)*(Z^2)*nelec*lambda/(mi*Tel + me*Ti)^(3./2.)

A = 1.8e-19*sqrt(me*mi)*(Z^2)*nelec*ni

d_nu_ie_dTe = (A*(mi*Tel + me*Ti)^(3./2.)*dlambda - $
         A*lambda*(1.5)*mi*(mi*Tel + me*Ti)^(1./2.))/(mi*Tel + me*Ti)^(3.0)

;print,'nu_ei...',nu_ie,lambda

return,d_nu_ie_dTe
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function d_nu_ie_dTi,mu,z,nelec,ni,Tel,Ti
; m is atomic mass
; z is atomic number
; for collision of given ion on thermal electron population
;-----------------------------------------------------------------------  
me = double(9.1e-28)    ;g
mp = double(1.67e-24)   ;g

mi = mu*mp

lambda = 0.0
dlambda = 0.0

if ((Ti*me/mi lt Tel) and (Tel le 10*z^2)) then begin
   lambda = 23.0 - alog(sqrt(nelec)*z*Tel^(-3./2.))
   dlambda = 0.0
endif

if ((Tel gt Ti*me/mi) and (Tel ge 10*Z^2)) then begin
   lambda = 24.0 - alog(sqrt(nelec)/Tel)
   dlambda = 0.0
endif

if (Tel lt Ti*Z*me/mi) then begin
   lambda = 30.0 - alog(sqrt(ni)*Ti^(-3./2.)*Z^2/mu)
   dlambda = 1.5/Ti
endif

;nu_ie = 1.8e-19*sqrt(me*mi)*(Z^2)*nelec*lambda/(mi*Tel + me*Ti)^(3./2.)

A = 1.8e-19*sqrt(me*mi)*(Z^2)*nelec*ni

d_nu_ie_dTi = (A*(mi*Tel + me*Ti)^(3./2.)*dlambda - $
         A*lambda*(1.5)*me*(mi*Tel + me*Ti)^(1./2.))/(mi*Tel + me*Ti)^(3.0)

;print,'nu_ei...',nu_ie,lambda

return,d_nu_ie_dTi
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
function d_nu_ee_dTe,nelec_c,nelec_h,Tel,Telh
; m is atomic mass
; z is atomic number
; for collision of given ion on thermal electron population
;-----------------------------------------------------------------------  
me = double(9.1e-28)    ;g
mp = double(1.67e-24)   ;g

lambda = 0.0
dlambda = 0.0

if (Tel le 10.0) then begin
   lambda = 23.0 - alog(sqrt(nelec_c)*Tel^(-3./2.))
   dlambda = 1.5/Tel
endif

if (Tel gt 10.0) then begin
   lambda = 24.0 - alog(sqrt(nelec_c)/Tel)
   dlambda = 1.0/(Tel)
endif

A = 1.8e-19*sqrt(me*me)*nelec_c*nelec_h

;nu_ee = 1.8e-19*sqrt(me*me)*nelec_c*nelec_h*lambda/(me*Tel + me*Telh)^(3./2.)

d_nu_ee_dTe = (A*(me*Tel + me*Telh)^(3./2.)*dlambda - $
               A*lambda*1.5*me*(me*Tel + me*Telh)^(1./2.))/(me*Tel + me*Telh)^(3.)

return,d_nu_ee_dTe
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

nu_ii = 1.8e-19*sqrt(m1*m2)*(Z1^2)*(Z2^2)*n1*n2*lambda/(m1*T2 + m2*T1)^(3./2.)


return,nu_ii
end
;-----------------------------------------------------------------------  



;-----------------------------------------------------------------------  
function d_nu_ii_self,mu1,mu2,z1,z2,n1,n2,T1,T2
; m is atomic mass
; z is atomic number
; for collision of given ion on ion population
;-----------------------------------------------------------------------  
mp = 1.67e-24   ;g

m1 = double(mu1*mp)
m2 = double(mu2*mp)

if (T1 eq T2) then begin
   d_nu_ii_self = 0.0
   GOTO,JUMP
endif

lambda = 23.0 - alog((Z1*Z2*(mu1+mu2)/(mu1*T2 + mu2*T1))* $
                     sqrt((n1*Z1^2)/T1 + (n2*Z2^2)/T2))


dlambda = 0.5*((n1*z1^2)/(T1^2))*((n1*(z1^2)/T1) + (n2*(z2^2)/T2))^(-1) + $
                mu2/(mu1*T2 + mu2*T1) 


A = 1.8e-19*sqrt(m1*m2)*(Z1^2)*(Z2^2)*n1*n2

;nu_ii = 1.8e-19*sqrt(m1*m2)*(Z1^2)*(Z2^2)*n2*lambda/(m1*T2 + m2*T1)^(3./2.)

d_nu_ii_self = (A*dlambda*(m1*T2 + m2*T1)^(3./2.) - $
                A*lambda*1.5*m2*(m1*T2 + m2*T1)^(1./2.))/(m1*T2 + m2*T1)^(3.)

JUMP:

return,d_nu_ii_self
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function d_nu_ii_other,mu1,mu2,z1,z2,n1,n2,T1,T2
; m is atomic mass
; z is atomic number
; for collision of given ion on ion population
;-----------------------------------------------------------------------  
mp = 1.67e-24   ;g

m1 = double(mu1*mp)
m2 = double(mu2*mp)

if (T1 eq T2) then begin
   d_nu_ii_other = 0.0
   GOTO,JUMP
endif

lambda = 23.0 - alog((Z1*Z2*(mu1+mu2)/(mu1*T2 + mu2*T1))* $
                     sqrt((n1*Z1^2)/T1 + (n2*Z2^2)/T2))


dlambda = 0.5*((n2*z2^2)/(T2^2))*((n1*(z1^2)/T1) + (n2*(z2^2)/T2))^(-1) + $
                mu1/(mu1*T2 + mu2*T1) 


A = 1.8e-19*sqrt(m1*m2)*(Z1^2)*(Z2^2)*n1*n2

;nu_ii = 1.8e-19*sqrt(m1*m2)*(Z1^2)*(Z2^2)*n2*lambda/(m1*T2 + m2*T1)^(3./2.)

d_nu_ii_other = (A*dlambda*(m1*T2 + m2*T1)^(3./2.) - $
                A*lambda*1.5*m1*(m1*T2 + m2*T1)^(1./2.))/(m1*T2 + m2*T1)^(3.)

JUMP:

return,d_nu_ii_other
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
function get_rad_coeff,n,T,r,delta_Te
;-----------------------------------------------------------------------  
nel = n.fc*n.nel

Tel = T.Tel+delta_Te
Tel_sp = T.Tel_sp+delta_Te
Tel_s2p = T.Tel_s2p+delta_Te
Tel_s3p = T.Tel_s3p+delta_Te
Tel_op = T.Tel_op+delta_Te
Tel_o2p = T.Tel_o2p+delta_Te

if (Tel le r.emistemp(1)) then Tel = r.emistemp(1)
if (Tel_sp le r.emistemp(1)) then Tel_sp = r.emistemp(1)
if (Tel_s2p le r.emistemp(1)) then Tel_s2p = r.emistemp(1)
if (Tel_s3p le r.emistemp(1)) then Tel_s3p = r.emistemp(1)
if (Tel_op le r.emistemp(1)) then Tel_op = r.emistemp(1)
if (Tel_o2p le r.emistemp(1)) then Tel_o2p = r.emistemp(1)

;print,Tel,Tel_sp,Tel_s2p,Tel_s3p,Tel_op,Tel_o2p

i=where(abs(r.emistemp-Tel_sp) eq min(abs(r.emistemp-Tel_sp)))
j=where(abs(r.emisden-nel) eq min(abs(r.emisden-nel)))

i = i(0)
j = j(0)

if (Tel_sp gt r.emistemp(i)) then ip = i+1
if (Tel_sp le r.emistemp(i)) then begin
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
dT1 = Tel_sp - r.emistemp(i)
dT2 = r.emistemp(ip) - Tel_sp
dn1 = nel - r.emisden(j)
dn2 = r.emisden(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSII = r.emisSII(i,j)*A4 + r.emisSII(ip,j)*A3 + r.emisSII(i,jp)*A2 + $
          r.emisSII(ip,jp)*A1
L = emisSII*n.nsp
;r.psp = emisSII

i=where(abs(r.emistemp-Tel_s2p) eq min(abs(r.emistemp-Tel_s2p)))
j=where(abs(r.emisden-nel) eq min(abs(r.emisden-nel)))

i = i(0)
j = j(0)

if (Tel_s2p gt r.emistemp(i)) then ip = i+1
if (Tel_s2p le r.emistemp(i)) then begin
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
dT1 = Tel_s2p - r.emistemp(i)
dT2 = r.emistemp(ip) - Tel_s2p
dn1 = nel - r.emisden(j)
dn2 = r.emisden(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSIII = r.emisSIII(i,j)*A4 + r.emisSIII(ip,j)*A3 + r.emisSIII(i,jp)*A2 + $
          r.emisSIII(ip,jp)*A1
L = L+ emisSIII*n.ns2p
;r.ps2p = emisSIII

i=where(abs(r.emistemp-Tel_s3p) eq min(abs(r.emistemp-Tel_s3p)))
j=where(abs(r.emisden-nel) eq min(abs(r.emisden-nel)))

i = i(0)
j = j(0)

if (Tel_s3p gt r.emistemp(i)) then ip = i+1
if (Tel_s3p le r.emistemp(i)) then begin
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
dT1 = Tel_s3p - r.emistemp(i)
dT2 = r.emistemp(ip) - Tel_s3p
dn1 = nel - r.emisden(j)
dn2 = r.emisden(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSIV = r.emisSIV(i,j)*A4 + r.emisSIV(ip,j)*A3 + r.emisSIV(i,jp)*A2 + $
          r.emisSIV(ip,jp)*A1
L = L + emisSIV*n.ns3p
;r.ps3p = emisSIV

i=where(abs(r.emistemp-Tel_op) eq min(abs(r.emistemp-Tel_op)))
j=where(abs(r.emisden-nel) eq min(abs(r.emisden-nel)))

i = i(0)
j = j(0)

if (Tel_op gt r.emistemp(i)) then ip = i+1
if (Tel_op le r.emistemp(i)) then begin
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
dT1 = Tel_op - r.emistemp(i)
dT2 = r.emistemp(ip) - Tel_op
dn1 = nel - r.emisden(j)
dn2 = r.emisden(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisOII = r.emisOII(i,j)*A4 + r.emisOII(ip,j)*A3 + r.emisOII(i,jp)*A2 + $
          r.emisOII(ip,jp)*A1
L = L + emisOII*n.nop
;r.pop = emisOII

i=where(abs(r.emistemp-Tel_o2p) eq min(abs(r.emistemp-Tel_o2p)))
j=where(abs(r.emisden-nel) eq min(abs(r.emisden-nel)))

i = i(0)
j = j(0)

if (Tel_o2p gt r.emistemp(i)) then ip = i+1
if (Tel_o2p le r.emistemp(i)) then begin
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
dT1 = Tel_o2p - r.emistemp(i)
dT2 = r.emistemp(ip) - Tel_o2p
dn1 = nel - r.emisden(j)
dn2 = r.emisden(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisOIII = r.emisOIII(i,j)*A4 + r.emisOIII(ip,j)*A3 + r.emisOIII(i,jp)*A2 + $
          r.emisOIII(ip,jp)*A1
L = L + emisOIII*n.no2p
;r.po2p = emisOIII

return,L
end
;-----------------------------------------------------------------------  



;----------------------------------------------------------------------- 
pro get_RsT_elec,n,T,nT,r,RsT
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

Teq = nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp)*(T.Tsp-T.Tel_sp) + $
      nu_ie(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p)*(T.Ts2p-T.Tel_s2p) + $
      nu_ie(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p)*(T.Ts3p-T.Tel_s3p) + $
      nu_ie(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top)*(T.Top-T.Tel_op) + $
      nu_ie(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p)*(T.To2p-T.Tel_o2p) + $
      nu_ee(n.fc*n.nel,n.fh*n.nel,T.Tel_el,T.Telh)*(T.Telh - T.Tel_el)

RsT(0) = Teq - (2./3.)*L - r.Transport*n.nel*T.Tel

r.Puv = L

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro get_JcT,n,T,nT,r,RsT,JcT
;-----------------------------------------------------------------------  
delta_Te = 0.1 ;eV


JcT(0,0) = (T.Tsp-T.Tel_sp)*d_nu_ie_dTe(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp) - $
            nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp) + $
           (T.Ts2p-T.Tel_s2p)*d_nu_ie_dTe(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p) - $
            nu_ie(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p) + $
           (T.Ts3p-T.Tel_s3p)*d_nu_ie_dTe(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p) - $
            nu_ie(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p) + $
           (T.Top-T.Tel_op)*d_nu_ie_dTe(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top) - $
            nu_ie(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top) + $
           (T.To2p-T.Tel_o2p)*d_nu_ie_dTe(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p) - $
            nu_ie(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p) + $
           (T.Telh - T.Tel_el)*d_nu_ee_dTe(n.fc*n.nel,n.fh*n.nel,T.Tel_el,T.Telh) - $
            nu_ee(n.fc*n.nel,n.fh*n.nel,T.Tel_el,T.Telh) - $
           (2./3.)*(get_rad_coeff(n,T,r,delta_Te)-get_rad_coeff(n,T,r,-1.0*delta_Te))/(2.0*delta_Te) - $
           n.nel*r.Transport
;print,'JcT(0,0)...',(T.Telh - T.Tel_el)*d_nu_ee_dTe(n.fc*n.nel,n.fh*n.nel,T.Tel_el,T.Telh) - $
;            nu_ee(n.fc*n.nel,n.fh*n.nel,T.Tel_el,T.Telh),(T.Telh - T.Tel_el),$
;            nu_ee(n.fc*n.nel,n.fh*n.nel,T.Tel_el,T.Telh),d_nu_ee_dTe(n.fc*n.nel,n.fh*n.nel,T.Tel_el,T.Telh)


JcT(1,0) = (T.Tsp-T.Tel_sp)*d_nu_ie_dTi(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp) + $
            nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp) 
JcT(2,0) = (T.Ts2p-T.Tel_s2p)*d_nu_ie_dTi(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p) + $
            nu_ie(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p) 
JcT(3,0) = (T.Ts3p-T.Tel_s3p)*d_nu_ie_dTi(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p) + $
            nu_ie(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p)
JcT(4,0) = (T.Top-T.Tel_op)*d_nu_ie_dTi(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top) + $
            nu_ie(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top)
JcT(5,0) = (T.To2p-T.Tel_o2p)*d_nu_ie_dTi(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p) + $
            nu_ie(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p) 


JcT(0,1) =  d_nu_ie_dTe(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp)*(T.Tel_sp-T.Tsp) + $
            nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp)
JcT(1,1) = -(r.isp*n.nsp*n.fc*n.nel + $       ;cold e ionization 
             r.isph*n.nsp*n.fh*n.nel + $      ;hot e ionization
             r.cx_k1*n.ns*n.nsp + $     ;S + S+ -> S+ + S
             r.cx_k8*n.no*n.nsp + $     ;O + S+ -> O+ + S
             r.cx_k13*n.no2p*n.nsp + $  ;O++ + S+ -> O+ + S++
             r.cx_k0*n.nsp*n.ns2p + $   ;S+ + S++ -> S++ + S+
             r.cx_k16*n.ns3p*n.nsp + $  ;S+++ + S+ -> S++ + S++
             r.Transport*n.nsp) + $
            d_nu_ii_self(32.0,32.0,1.0,2.0,n.nsp,n.ns2p,T.Tsp,T.Ts2p)*(T.Ts2p-T.Tsp) - $
              nu_ii(32.0,32.0,1.0,2.0,n.nsp,n.ns2p,T.Tsp,T.Ts2p) + $
            d_nu_ii_self(32.0,32.0,1.0,3.0,n.nsp,n.ns3p,T.Tsp,T.Ts3p)*(T.Ts3p-T.Tsp) - $
              nu_ii(32.0,32.0,1.0,3.0,n.nsp,n.ns3p,T.Tsp,T.Ts3p) + $
            d_nu_ii_self(32.0,16.0,1.0,1.0,n.nsp,n.nop,T.Tsp,T.Top)*(T.Top-T.Tsp) - $
              nu_ii(32.0,16.0,1.0,1.0,n.nsp,n.nop,T.Tsp,T.Top) + $
            d_nu_ii_self(32.0,16.0,1.0,2.0,n.nsp,n.no2p,T.Tsp,T.To2p)*(T.To2p-T.Tsp) - $
              nu_ii(32.0,16.0,1.0,2.0,n.nsp,n.no2p,T.Tsp,T.To2p) + $
            d_nu_ie_dTi(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp)*(T.Tel_sp-T.Tsp) - $
              nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp) + $
            d_nu_ie_dTi(32.0,1.0,n.fh*n.nel,n.nsp,T.Telh,T.Tsp)*(T.Telh-T.Tsp) - $
              nu_ie(32.0,1.0,n.fh*n.nel,n.nsp,T.Telh,T.Tsp)
          
JcT(2,1) = r.cx_k2*n.ns*n.ns2p + $  ;S + S++ -> S+ + S+
           r.cx_k12*n.no*n.ns2p + $ ;O + S++ -> O+ + S+
           r.cx_k0*n.nsp*n.ns2p + $ ;S+ + S++ -> S++ + S+
           r.rs2p*n.ns2p*n.nel + $  ;S++ + e -> S+
           d_nu_ii_other(32.0,32.0,1.0,2.0,n.nsp,n.ns2p,T.Tsp,T.Ts2p)*(T.Ts2p-T.Tsp) + $
           nu_ii(32.0,32.0,1.0,2.0,n.nsp,n.ns2p,T.Tsp,T.Ts2p)
JcT(3,1) = d_nu_ii_other(32.0,32.0,1.0,3.0,n.nsp,n.ns3p,T.Tsp,T.Ts3p)*(T.Ts3p-T.Tsp) + $
           nu_ii(32.0,32.0,1.0,3.0,n.nsp,n.ns3p,T.Tsp,T.Ts3p)
JcT(4,1) = d_nu_ii_other(32.0,16.0,1.0,1.0,n.nsp,n.nop,T.Tsp,T.Top)*(T.Top-T.Tsp) + $
           nu_ii(32.0,16.0,1.0,1.0,n.nsp,n.nop,T.Tsp,T.Top)
JcT(5,1) = d_nu_ii_other(32.0,16.0,1.0,2.0,n.nsp,n.no2p,T.Tsp,T.To2p)*(T.To2p-T.Tsp) + $
           nu_ii(32.0,16.0,1.0,2.0,n.nsp,n.no2p,T.Tsp,T.To2p)


JcT(0,2) = d_nu_ie_dTe(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p)*(T.Tel_s2p-T.Ts2p) + $
           nu_ie(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p)
JcT(1,2) = r.isp*n.nsp*n.fc*n.nel + $               ;cold e ionization 
           r.isph*n.nsp*n.fh*n.nel + $       ;hot e ionization
           r.cx_k13*n.no2p*n.nsp + $          ;O++ + S+ -> O+ + S++ 
           r.cx_k0*n.nsp*n.ns2p + $           ;S+ + S++ -> S++ + S+
           r.cx_k16*n.ns3p*n.nsp + $          ;S+++ + S+ -> S++ + S++
           d_nu_ii_other(32.0,32.0,2.0,1.0,n.ns2p,n.nsp,T.Ts2p,T.Tsp)*(T.Tsp-T.Ts2p) + $ 
           nu_ii(32.0,32.0,2.0,1.0,n.ns2p,n.nsp,T.Ts2p,T.Tsp)
JcT(2,2) =-(r.is2p*n.ns2p*n.fc*n.nel + $       ;cold e ionization 
            r.is2ph*n.ns2p*n.fh*n.nel + $      ;hot e ionization
            r.cx_k2*n.ns*n.ns2p + $     ;S + S++ -> S+ + S+
            r.cx_k3*n.ns*n.ns2p + $     ;S + S++ -> S++ + S
            r.cx_k12*n.no*n.ns2p + $    ;O + S++ -> O+ + S+
            r.cx_k15*n.no2p*n.ns2p + $  ;O++ + S++ -> O+ + S+++
            r.cx_k0*n.nsp*n.ns2p + $    ;S+ + S++ -> S++ + S+
            r.Transport*n.ns2p + $   
            r.rs2p*n.ns2p*n.nel) + $
           d_nu_ii_self(32.0,32.0,2.0,1.0,n.ns2p,n.nsp,T.Ts2p,T.Tsp)*(T.Tsp-T.Ts2p) - $
           nu_ii(32.0,32.0,2.0,1.0,n.ns2p,n.nsp,T.Ts2p,T.Tsp) + $
           d_nu_ii_self(32.0,32.0,2.0,3.0,n.ns2p,n.ns3p,T.Ts2p,T.Ts3p)*(T.Ts3p-T.Ts2p) - $
           nu_ii(32.0,32.0,2.0,3.0,n.ns2p,n.ns3p,T.Ts2p,T.Ts3p) + $
           d_nu_ii_self(32.0,16.0,2.0,1.0,n.ns2p,n.nop,T.Ts2p,T.Top)*(T.Top-T.Ts2p) - $
           nu_ii(32.0,16.0,2.0,1.0,n.ns2p,n.nop,T.Ts2p,T.Top) + $
           d_nu_ii_self(32.0,16.0,2.0,2.0,n.ns2p,n.no2p,T.Ts2p,T.To2p)*(T.To2p-T.Ts2p) - $
           nu_ii(32.0,16.0,2.0,2.0,n.ns2p,n.no2p,T.Ts2p,T.To2p) + $
           d_nu_ie_dTi(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p)*(T.Tel_s2p-T.Ts2p) - $
           nu_ie(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p) + $
           d_nu_ie_dTi(32.0,2.0,n.fh*n.nel,n.ns2p,T.Telh,T.Ts2p)*(T.Telh-T.Ts2p) - $
           nu_ie(32.0,2.0,n.fh*n.nel,n.ns2p,T.Telh,T.Ts2p)
JcT(3,2) = r.cx_k4*n.ns*n.ns3p + $           ;S + S+++ -> S+ + S++
           r.cx_k14*n.no*n.ns3p + $          ;O + S+++ -> O+ + S++
           r.cx_k16*n.ns3p*n.nsp + $         ;S+++ + S+ -> S++ + S++
           r.rs3p*n.ns3p*n.nel + $           ;S+++ + e -> S++     
           d_nu_ii_other(32.0,32.0,2.0,3.0,n.ns2p,n.ns3p,T.Ts2p,T.Ts3p)*(T.Ts3p-T.Ts2p) + $
           nu_ii(32.0,32.0,2.0,3.0,n.ns2p,n.ns3p,T.Ts2p,T.Ts3p)
JcT(4,2) = d_nu_ii_other(32.0,16.0,2.0,1.0,n.ns2p,n.nop,T.Ts2p,T.Top)*(T.Top-T.Ts2p) + $
           nu_ii(32.0,16.0,2.0,1.0,n.ns2p,n.nop,T.Ts2p,T.Top)
JcT(5,2) = d_nu_ii_other(32.0,16.0,2.0,2.0,n.ns2p,n.no2p,T.Ts2p,T.To2p)*(T.To2p-T.Ts2p) + $
           nu_ii(32.0,16.0,2.0,2.0,n.ns2p,n.no2p,T.Ts2p,T.To2p)


JcT(0,3) = d_nu_ie_dTe(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p)*(T.Tel_s3p-T.Ts3p) + $
           nu_ie(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p)
JcT(1,3) = d_nu_ii_other(32.0,32.0,3.0,1.0,n.ns3p,n.nsp,T.Ts3p,T.Tsp)*(T.Tsp-T.Ts3p) + $ 
           nu_ii(32.0,32.0,3.0,1.0,n.ns3p,n.nsp,T.Ts3p,T.Tsp)
JcT(2,3) = r.is2p*n.ns2p*n.fc*n.nel + $            ;cold e ionization 
           r.is2ph*n.ns2p*n.fh*n.nel + $           ;hot e ionization
           r.cx_k15*n.no2p*n.ns2p + $              ;O++ + S++ -> O+ + S+++
           d_nu_ii_other(32.0,32.0,3.0,2.0,n.ns3p,n.ns2p,T.Ts3p,T.Ts2p)*(T.Ts2p-T.Ts3p) + $
           nu_ii(32.0,32.0,3.0,2.0,n.ns3p,n.ns2p,T.Ts3p,T.Ts2p)
JcT(3,3) = r.cx_k14*n.no*n.ns3p + $      ;O + S+++ -> O+ + S++
           r.cx_k16*n.ns3p*n.nsp + $     ;S+++ + S+ -> S++ + S++
           r.Transport*n.ns3p + $
           r.rs3p*n.ns3p*n.nel + $       ;S+++ + e -> S++   
           d_nu_ii_self(32.0,32.0,3.0,1.0,n.ns3p,n.nsp,T.Ts3p,T.Tsp)*(T.Tsp-T.Ts3p) - $
           nu_ii(32.0,32.0,3.0,1.0,n.ns3p,n.nsp,T.Ts3p,T.Tsp) + $
           d_nu_ii_self(32.0,32.0,3.0,2.0,n.ns3p,n.ns2p,T.Ts3p,T.Ts2p)*(T.Ts2p-T.Ts3p) - $
           nu_ii(32.0,32.0,3.0,2.0,n.ns3p,n.ns2p,T.Ts3p,T.Ts2p) + $
           d_nu_ii_self(32.0,16.0,3.0,1.0,n.ns3p,n.nop,T.Ts3p,T.Top)*(T.Top-T.Ts3p) - $
           nu_ii(32.0,16.0,3.0,1.0,n.ns3p,n.nop,T.Ts3p,T.Top) + $
           d_nu_ii_self(32.0,16.0,3.0,2.0,n.ns3p,n.no2p,T.Ts3p,T.To2p)*(T.To2p-T.Ts3p) - $
           nu_ii(32.0,16.0,3.0,2.0,n.ns3p,n.no2p,T.Ts3p,T.To2p) + $
           d_nu_ie_dTi(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p)*(T.Tel_s3p-T.Ts3p) - $
           nu_ie(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p) + $
           d_nu_ie_dTi(32.0,3.0,n.fh*n.nel,n.ns3p,T.Telh,T.Ts3p)*(T.Telh-T.Ts3p) - $
           nu_ie(32.0,3.0,n.fh*n.nel,n.ns3p,T.Telh,T.Ts3p) 
JcT(4,3) = d_nu_ii_other(32.0,16.0,3.0,1.0,n.ns3p,n.nop,T.Ts3p,T.Top)*(T.Top-T.Ts3p) + $
           nu_ii(32.0,16.0,3.0,1.0,n.ns3p,n.nop,T.Ts3p,T.Top)
JcT(5,3) = d_nu_ii_other(32.0,16.0,3.0,2.0,n.ns3p,n.no2p,T.Ts3p,T.To2p)*(T.To2p-T.Ts3p) + $
           nu_ii(32.0,16.0,3.0,2.0,n.ns3p,n.no2p,T.Ts3p,T.To2p)

JcT(0,4) = d_nu_ie_dTe(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top)*(T.Tel_op-T.Top) + $
           nu_ie(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top)
JcT(1,4) = d_nu_ii_other(16.0,32.0,1.0,1.0,n.nop,n.nsp,T.Top,T.Tsp)*(T.Tsp-T.Top) + $
           nu_ii(16.0,32.0,1.0,1.0,n.nop,n.nsp,T.Top,T.Tsp)
JcT(2,4) = d_nu_ii_other(16.0,32.0,1.0,2.0,n.nop,n.ns2p,T.Top,T.Ts2p)*(T.Ts2p-T.Top) + $
           nu_ii(16.0,32.0,1.0,2.0,n.nop,n.ns2p,T.Top,T.Ts2p)
JcT(3,4) = d_nu_ii_other(16.0,32.0,1.0,3.0,n.nop,n.ns3p,T.Top,T.Ts3p)*(T.Ts3p-T.Top) + $
           nu_ii(16.0,32.0,1.0,3.0,n.nop,n.ns3p,T.Top,T.Ts3p)
JcT(4,4) = -(r.iop*n.nop*n.fc*n.nel + $
             r.ioph*n.nop*n.fh*n.nel + $
             r.cx_k5*n.no*n.nop + $     ;O + O+ -> O+ + O
             r.cx_k9*n.ns*n.nop + $     ;S + O+ -> S+ + O
             r.Transport*n.nop) + $
           d_nu_ii_self(16.0,16.0,1.0,2.0,n.nop,n.no2p,T.Top,T.To2p)*(T.To2p-T.Top) - $
           nu_ii(16.0,16.0,1.0,2.0,n.nop,n.no2p,T.Top,T.To2p) + $
           d_nu_ii_self(16.0,32.0,1.0,1.0,n.nop,n.nsp,T.Top,T.Tsp)*(T.Tsp-T.Top) - $
           nu_ii(16.0,32.0,1.0,1.0,n.nop,n.nsp,T.Top,T.Tsp) + $
           d_nu_ii_self(16.0,32.0,1.0,2.0,n.nop,n.ns2p,T.Top,T.Ts2p)*(T.Ts2p-T.Top) - $
           nu_ii(16.0,32.0,1.0,2.0,n.nop,n.ns2p,T.Top,T.Ts2p) + $
           d_nu_ii_self(16.0,32.0,1.0,3.0,n.nop,n.ns3p,T.Top,T.Ts3p)*(T.Ts3p-T.Top) - $
           nu_ii(16.0,32.0,1.0,3.0,n.nop,n.ns3p,T.Top,T.Ts3p) + $
           d_nu_ie_dTi(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top)*(T.Tel_op-T.Top) - $
           nu_ie(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top) + $
           d_nu_ie_dTi(16.0,1.0,n.fh*n.nel,n.nop,T.Telh,T.Top)*(T.Telh-T.Top) - $
           nu_ie(16.0,1.0,n.fh*n.nel,n.nop,T.Telh,T.Top)
JcT(5,4) = r.cx_k6*n.no*n.no2p + $    ;O + O++ -> O+ + O+
           r.cx_k10*n.ns*n.no2p + $   ;S + O++ -> S+ + O+
           r.cx_k11*n.ns*n.no2p + $   ;S + O++ -> S+ + O+ + e
           r.cx_k13*n.no2p*n.nsp + $  ;O++ + S+ -> O+ + S++
           r.cx_k15*n.no2p*n.ns2p + $ ;O++ + S++ -> O+ + S+++
           r.ro2p*n.no2p*n.nel + $
           d_nu_ii_other(16.0,16.0,1.0,2.0,n.nop,n.no2p,T.Top,T.To2p)*(T.To2p-T.Top) + $
           nu_ii(16.0,16.0,1.0,2.0,n.nop,n.no2p,T.Top,T.To2p)

JcT(0,5) = d_nu_ie_dTe(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p)*(T.Tel_o2p-T.To2p) + $
           nu_ie(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p)
JcT(1,5) = d_nu_ii_other(16.0,32.0,2.0,1.0,n.no2p,n.nsp,T.To2p,T.Tsp)*(T.Tsp-T.To2p) + $
           nu_ii(16.0,32.0,2.0,1.0,n.no2p,n.nsp,T.To2p,T.Tsp)
JcT(2,5) = d_nu_ii_other(16.0,32.0,2.0,2.0,n.no2p,n.ns2p,T.To2p,T.Ts2p)*(T.Ts2p-T.To2p) + $
           nu_ii(16.0,32.0,2.0,2.0,n.no2p,n.ns2p,T.To2p,T.Ts2p)
JcT(3,5) = d_nu_ii_other(16.0,32.0,2.0,3.0,n.no2p,n.ns3p,T.To2p,T.Ts3p)*(T.Ts3p-T.To2p) + $
           nu_ii(16.0,32.0,2.0,3.0,n.no2p,n.ns3p,T.To2p,T.Ts3p)
JcT(4,5) = r.iop*n.nop*n.fc*n.nel + $           ;cold e ionization 
           r.ioph*n.nop*n.fh*n.nel + $          ;hot e ionization
           d_nu_ii_other(16.0,16.0,2.0,1.0,n.no2p,n.nop,T.To2p,T.Top)*(T.Top-T.To2p) + $
           nu_ii(16.0,16.0,2.0,1.0,n.no2p,n.nop,T.To2p,T.Top)
JcT(5,5) = -(r.cx_k6*n.no*n.no2p + $     ;O + O++ -> O+ + O+
             r.cx_k7*n.no*n.no2p + $     ;O + O++ -> O++ + O
             r.cx_k10*n.ns*n.no2p + $    ;S + O++ -> S+ + O+
             r.cx_k11*n.ns*n.no2p + $    ;S + O++ -> S++ + O+ + e
             r.cx_k13*n.no2p*n.nsp + $   ;O++ + S+ -> O+ + S++
             r.cx_k15*n.no2p*n.ns2p + $  ;O++ + S++ -> O+ + S+++
             r.Transport*n.no2p + $
             r.ro2p*n.no2p*n.nel) + $   
           d_nu_ii_self(16.0,16.0,2.0,1.0,n.no2p,n.nop,T.To2p,T.Top)*(T.Top-T.To2p) - $
           nu_ii(16.0,16.0,2.0,1.0,n.no2p,n.nop,T.To2p,T.Top) + $
           d_nu_ii_self(16.0,32.0,2.0,1.0,n.no2p,n.nsp,T.To2p,T.Tsp)*(T.Tsp-T.To2p) - $
           nu_ii(16.0,32.0,2.0,1.0,n.no2p,n.nsp,T.To2p,T.Tsp) + $
           d_nu_ii_self(16.0,32.0,2.0,2.0,n.no2p,n.ns2p,T.To2p,T.Ts2p)*(T.Ts2p-T.To2p) - $
           nu_ii(16.0,32.0,2.0,2.0,n.no2p,n.ns2p,T.To2p,T.Ts2p) + $
           d_nu_ii_self(16.0,32.0,2.0,3.0,n.no2p,n.ns3p,T.To2p,T.Ts3p)*(T.Ts3p-T.To2p)+$
           nu_ii(16.0,32.0,2.0,3.0,n.no2p,n.ns3p,T.To2p,T.Ts3p) + $
           d_nu_ie_dTi(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p)*(T.Tel_o2p-T.To2p) - $
           nu_ie(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p) + $
           d_nu_ie_dTi(16.0,2.0,n.fh*n.nel,n.no2p,T.Telh,T.To2p)*(T.Telh-T.To2p) - $
           nu_ie(16.0,2.0,n.fh*n.nel,n.no2p,T.Telh,T.To2p)
 
return
end
;-----------------------------------------------------------------------  

;-----------------------------------------------------------------------  
pro get_delta_nT,RsT,JcT,delta_nT
;-----------------------------------------------------------------------  
;dt = 5e7
;Id = fltarr(6,6)
;for i = 0,5 do begin
;  Id(i,i) = 1.0/dt
;endfor

;A = Id + JcT
;print,' '
;print,JcT

A = JcT
B = -RsT

;Ax = B
;Jc delta_n = - Rs_arr

;print,'a...',a
;print,'b...',b

LUDC, A, INDEX,/double
delta_nT = LUSOL(A, INDEX, B,/double)

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro update_nT,n,T,nT,delta_nT
;-----------------------------------------------------------------------  


T.Tsp = T.Tsp + delta_nT(1)
T.Ts2p = T.Ts2p + delta_nT(2)
T.Ts3p = T.Ts3p + delta_nT(3)
T.Top = T.Top + delta_nT(4)
T.To2p = T.To2p + delta_nT(5)
T.Tel = T.Tel + delta_nT(0)

;nT.el = nT.el + delta_nT(0)
;nT.sp = nT.sp + delta_nT(1)
;nT.s2p = nT.s2p + delta_nT(2)
;nT.s3p = nT.s3p + delta_nT(3)
;nT.op = nT.op + delta_nT(4)
;nT.o2p = nT.o2p + delta_nT(5)

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro get_RsT_sp,n,T,nT,r,RsT
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
     r.cx_k0*n.nsp*n.ns2p*T.Ts2p + $ ;S+ + S++ -> S++ + S+
     r.rs2p*n.ns2p*n.nel*T.Ts2p      ;S++ + e -> S+

r.Sei = (r.is*n.ns*n.fc*n.nel + r.ish*n.ns*n.fh*n.nel)*T.Tpu_s
;r.Seih = r.ish*n.ns*n.fh*n.nel*T.Tpu_s
r.Secx = (r.cx_k4*n.ns*n.ns3p + r.cx_k10*n.ns*n.no2p + $
         r.cx_k1*n.ns*n.nsp + r.cx_k2*n.ns*n.ns2p + $
         r.cx_k9*n.ns*n.nop+r.cx_k3*n.ns*n.ns2p + $
         r.cx_k11*n.ns*n.no2p)*T.Tpu_s

L= r.isp*n.nsp*n.fc*n.nel*T.Tsp + $       ;cold e ionization 
   r.isph*n.nsp*n.fh*n.nel*T.Tsp + $      ;hot e ionization
   r.cx_k1*n.ns*n.nsp*T.Tsp + $     ;S + S+ -> S+ + S
   r.cx_k8*n.no*n.nsp*T.Tsp + $     ;O + S+ -> O+ + S
   r.cx_k13*n.no2p*n.nsp*T.Tsp + $  ;O++ + S+ -> O+ + S++
   r.cx_k0*n.nsp*n.ns2p*T.Tsp + $   ;S+ + S++ -> S++ + S+
   r.cx_k16*n.ns3p*n.nsp*T.Tsp + $  ;S+++ + S+ -> S++ + S++
   r.Transport*n.nsp*T.Tsp    

Teq= nu_ii(32.0,32.0,1.0,2.0,n.nsp,n.ns2p,T.Tsp,T.Ts2p)*(T.Ts2p-T.Tsp)+$
   nu_ii(32.0,32.0,1.0,3.0,n.nsp,n.ns3p,T.Tsp,T.Ts3p)*(T.Ts3p-T.Tsp)+$
   nu_ii(32.0,16.0,1.0,1.0,n.nsp,n.nop,T.Tsp,T.Top)*(T.Top-T.Tsp)+$
   nu_ii(32.0,16.0,1.0,2.0,n.nsp,n.no2p,T.Tsp,T.To2p)*(T.To2p-T.Tsp)+$
   nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp)*(T.Tel_sp-T.Tsp) + $
   nu_ie(32.0,1.0,n.fh*n.nel,n.nsp,T.Telh,T.Tsp)*(T.Telh-T.Tsp) 
 
;print,'nu_ie (S+)...',nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp)

RsT(1) = S - L + Teq

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro get_RsT_s2p,n,T,nT,r,RsT
;-----------------------------------------------------------------------  
S = (r.isp*n.nsp*n.fc*n.nel + $               ;cold e ionization 
     r.isph*n.nsp*n.fh*n.nel)*T.Tsp + $       ;hot e ionization
     r.cx_k3*n.ns*n.ns2p*T.Tpu_s + $          ;S + S++ -> S++ + S
     r.cx_k4*n.ns*n.ns3p*T.Ts3p + $           ;S + S+++ -> S+ + S++
     r.cx_k11*n.ns*n.no2p*T.Tpu_s + $         ;S + O++ -> S++ + O+ + e
     r.cx_k13*n.no2p*n.nsp*T.Tsp + $          ;O++ + S+ -> O+ + S++ 
     r.cx_k14*n.no*n.ns3p*T.Ts3p + $          ;O + S+++ -> O+ + S++
     r.cx_k0*n.nsp*n.ns2p*T.Tsp + $           ;S+ + S++ -> S++ + S+
     r.cx_k16*n.ns3p*n.nsp*T.Tsp + $          ;S+++ + S+ -> S++ + S++
     r.cx_k16*n.ns3p*n.nsp*T.Ts3p + $         ;S+++ + S+ -> S++ + S++
     r.rs3p*n.ns3p*n.nel*T.Ts3p               ;S+++ + e -> S++     

L= r.is2p*n.ns2p*n.fc*n.nel*T.Ts2p + $       ;cold e ionization 
   r.is2ph*n.ns2p*n.fh*n.nel*T.Ts2p + $      ;hot e ionization
   r.cx_k2*n.ns*n.ns2p*T.Ts2p + $     ;S + S++ -> S+ + S+
   r.cx_k3*n.ns*n.ns2p*T.Ts2p + $     ;S + S++ -> S++ + S
   r.cx_k12*n.no*n.ns2p*T.Ts2p + $    ;O + S++ -> O+ + S+
   r.cx_k15*n.no2p*n.ns2p*T.Ts2p + $  ;O++ + S++ -> O+ + S+++
   r.cx_k0*n.nsp*n.ns2p*T.Ts2p + $    ;S+ + S++ -> S++ + S+
   r.Transport*n.ns2p*T.Ts2p + $   
   r.rs2p*n.ns2p*n.nel*T.Ts2p

Teq= nu_ii(32.0,32.0,2.0,1.0,n.ns2p,n.nsp,T.Ts2p,T.Tsp)*(T.Tsp-T.Ts2p)+$
   nu_ii(32.0,32.0,2.0,3.0,n.ns2p,n.ns3p,T.Ts2p,T.Ts3p)*(T.Ts3p-T.Ts2p)+$
   nu_ii(32.0,16.0,2.0,1.0,n.ns2p,n.nop,T.Ts2p,T.Top)*(T.Top-T.Ts2p)+$
   nu_ii(32.0,16.0,2.0,2.0,n.ns2p,n.no2p,T.Ts2p,T.To2p)*(T.To2p-T.Ts2p)+$
   nu_ie(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p)*(T.Tel_s2p-T.Ts2p) + $
   nu_ie(32.0,2.0,n.fh*n.nel,n.ns2p,T.Telh,T.Ts2p)*(T.Telh-T.Ts2p) 
 
RsT(2) = S - L + Teq

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro get_RsT_s3p,n,T,nT,r,RsT
;-----------------------------------------------------------------------  
S = (r.is2p*n.ns2p*n.fc*n.nel + $               ;cold e ionization 
     r.is2ph*n.ns2p*n.fh*n.nel)*T.Ts2p + $      ;hot e ionization
     r.cx_k15*n.no2p*n.ns2p*T.Ts2p              ;O++ + S++ -> O+ + S+++

L= r.cx_k14*n.no*n.ns3p*T.Ts3p + $      ;O + S+++ -> O+ + S++
   r.cx_k16*n.ns3p*n.nsp*T.Ts3p + $     ;S+++ + S+ -> S++ + S++
   r.Transport*n.ns3p*T.Ts3p + $
   r.rs3p*n.ns3p*n.nel*T.Ts3p           ;S+++ + e -> S++   

Teq= nu_ii(32.0,32.0,3.0,1.0,n.ns3p,n.nsp,T.Ts3p,T.Tsp)*(T.Tsp-T.Ts3p)+$
   nu_ii(32.0,32.0,3.0,2.0,n.ns3p,n.ns2p,T.Ts3p,T.Ts2p)*(T.Ts2p-T.Ts3p)+$
   nu_ii(32.0,16.0,3.0,1.0,n.ns3p,n.nop,T.Ts3p,T.Top)*(T.Top-T.Ts3p)+$
   nu_ii(32.0,16.0,3.0,2.0,n.ns3p,n.no2p,T.Ts3p,T.To2p)*(T.To2p-T.Ts3p)+$
   nu_ie(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p)*(T.Tel_s3p-T.Ts3p) + $
   nu_ie(32.0,3.0,n.fh*n.nel,n.ns3p,T.Telh,T.Ts3p)*(T.Telh-T.Ts3p) 
 
RsT(3) = S - L + Teq

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro get_RsT_op,n,T,nT,r,RsT
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
     r.cx_k15*n.no2p*n.ns2p*T.To2p + $ ;O++ + S++ -> O+ + S+++
     r.ro2p*n.no2p*n.nel*T.To2p

r.Oei = (r.io*n.no*n.fc*n.nel + r.ioh*n.no*n.fh*n.nel)*T.Tpu_o
r.Oecx = (r.cx_k5*n.no*n.nop + r.cx_k6*n.no*n.no2p + $
         r.cx_k8*n.no*n.nsp + r.cx_k12*n.no*n.ns2p + $
         r.cx_k14*n.no*n.ns3p+r.cx_k7*n.no*n.no2p)*T.Tpu_o

L= r.iop*n.nop*n.fc*n.nel*T.Top + $
   r.ioph*n.nop*n.fh*n.nel*T.Top + $
   r.cx_k5*n.no*n.nop*T.Top + $     ;O + O+ -> O+ + O
   r.cx_k9*n.ns*n.nop*T.Top + $     ;S + O+ -> S+ + O
   r.Transport*n.nop*T.Top    

Teq= nu_ii(16.0,16.0,1.0,2.0,n.nop,n.no2p,T.Top,T.To2p)*(T.To2p-T.Top)+$
   nu_ii(16.0,32.0,1.0,1.0,n.nop,n.nsp,T.Top,T.Tsp)*(T.Tsp-T.Top)+$
   nu_ii(16.0,32.0,1.0,2.0,n.nop,n.ns2p,T.Top,T.Ts2p)*(T.Ts2p-T.Top)+$
   nu_ii(16.0,32.0,1.0,3.0,n.nop,n.ns3p,T.Top,T.Ts3p)*(T.Ts3p-T.Top)+$
   nu_ie(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top)*(T.Tel_op-T.Top) + $
   nu_ie(16.0,1.0,n.fh*n.nel,n.nop,T.Telh,T.Top)*(T.Telh-T.Top) 
 
RsT(4) = S - L + Teq

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro get_RsT_o2p,n,T,nT,r,RsT
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
   r.Transport*n.no2p*T.To2p + $
   r.ro2p*n.no2p*n.nel*T.To2p   

Teq= nu_ii(16.0,16.0,2.0,1.0,n.no2p,n.nop,T.To2p,T.Top)*(T.Top-T.To2p)+$
   nu_ii(16.0,32.0,2.0,1.0,n.no2p,n.nsp,T.To2p,T.Tsp)*(T.Tsp-T.To2p)+$
   nu_ii(16.0,32.0,2.0,2.0,n.no2p,n.ns2p,T.To2p,T.Ts2p)*(T.Ts2p-T.To2p)+$
   nu_ii(16.0,32.0,2.0,3.0,n.no2p,n.ns3p,T.To2p,T.Ts3p)*(T.Ts3p-T.To2p)+$
   nu_ie(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p)*(T.Tel_o2p-T.To2p) + $
   nu_ie(16.0,2.0,n.fh*n.nel,n.no2p,T.Telh,T.To2p)*(T.Telh-T.To2p) 
 
RsT(5) = S - L + Teq

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro get_neutral_loss_rate,n,r,T,nl
;-----------------------------------------------------------------------  

nl.sei = r.is*n.ns*n.fc*n.nel
nl.seih = r.ish*n.ns*n.fh*n.nel
nl.oei = r.io*n.no*n.fc*n.nel
nl.oeih = r.ioh*n.no*n.fh*n.nel

nl.k1 = r.cx_k1*n.ns*n.nsp
nl.k2 = r.cx_k2*n.ns*n.ns2p
nl.k3 = r.cx_k3*n.ns*n.ns2p
nl.k4 = r.cx_k4*n.ns*n.ns3p
nl.k5 = r.cx_k5*n.no*n.nop
nl.k6 = r.cx_k6*n.no*n.no2p
nl.k7 = r.cx_k7*n.no*n.no2p
nl.k8 = r.cx_k8*n.no*n.nsp
nl.k9 = r.cx_k9*n.ns*n.nop
nl.k10 = r.cx_k10*n.ns*n.no2p
nl.k11 = r.cx_k11*n.ns*n.no2p
nl.k12 = r.cx_k12*n.no*n.ns2p
nl.k14 = r.cx_k14*n.no*n.ns3p

nl.fast_S_k1 = r.cx_k1*n.ns*n.nsp*T.Tsp
nl.fast_S_k3 = r.cx_k3*n.ns*n.ns2p*T.Ts2p
nl.fast_O_k5 = r.cx_k5*n.no*n.nop*T.Top
nl.fast_O_k7 = r.cx_k7*n.no*n.no2p*T.To2p
nl.fast_O_k9 = r.cx_k9*n.ns*n.nop*T.Top

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


;-----------------------------------------------------------------------  
pro energy_conservation,n,r,T,nl,Ec
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

print,' ' 
print,'Energy conservation (eV cm^-3 s^-1):'
print,' '
print,'   S ionization cold........',Ec.s_ion
print,'   S ionization hot.........',Ec.s_ion_h
print,'   O ionization cold........',Ec.o_ion
print,'   O ionization hot.........',Ec.o_ion_h
print,' '
print,'   S chex...................',Ec.s_cx
print,'   O chex...................',Ec.o_cx
print,' '
print,'   Total S input............',Ec.s_tot_in
print,'   Total O input............',Ec.o_tot_in

Teq_el = nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp)*(T.Tsp-T.Tel_sp) + $
      nu_ie(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p)*(T.Ts2p-T.Tel_s2p) + $
      nu_ie(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p)*(T.Ts3p-T.Tel_s3p) + $
      nu_ie(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top)*(T.Top-T.Tel_op) + $
      nu_ie(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p)*(T.To2p-T.Tel_o2p) + $
      nu_ee(n.fc*n.nel,n.fh*n.nel,T.Tel_el,T.Telh)*(T.Telh - T.Tel_el)

Teq_sp= nu_ii(32.0,32.0,1.0,2.0,n.nsp,n.ns2p,T.Tsp,T.Ts2p)*(T.Ts2p-T.Tsp)+$
   nu_ii(32.0,32.0,1.0,3.0,n.nsp,n.ns3p,T.Tsp,T.Ts3p)*(T.Ts3p-T.Tsp)+$
   nu_ii(32.0,16.0,1.0,1.0,n.nsp,n.nop,T.Tsp,T.Top)*(T.Top-T.Tsp)+$
   nu_ii(32.0,16.0,1.0,2.0,n.nsp,n.no2p,T.Tsp,T.To2p)*(T.To2p-T.Tsp)+$
   nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp)*(T.Tel_sp-T.Tsp) + $
   nu_ie(32.0,1.0,n.fh*n.nel,n.nsp,T.Telh,T.Tsp)*(T.Telh-T.Tsp) 
 
Teq_s2p= nu_ii(32.0,32.0,2.0,1.0,n.ns2p,n.nsp,T.Ts2p,T.Tsp)*(T.Tsp-T.Ts2p)+$
   nu_ii(32.0,32.0,2.0,3.0,n.ns2p,n.ns3p,T.Ts2p,T.Ts3p)*(T.Ts3p-T.Ts2p)+$
   nu_ii(32.0,16.0,2.0,1.0,n.ns2p,n.nop,T.Ts2p,T.Top)*(T.Top-T.Ts2p)+$
   nu_ii(32.0,16.0,2.0,2.0,n.ns2p,n.no2p,T.Ts2p,T.To2p)*(T.To2p-T.Ts2p)+$
   nu_ie(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p)*(T.Tel_s2p-T.Ts2p) + $
   nu_ie(32.0,2.0,n.fh*n.nel,n.ns2p,T.Telh,T.Ts2p)*(T.Telh-T.Ts2p) 

Teq_s3p= nu_ii(32.0,32.0,3.0,1.0,n.ns3p,n.nsp,T.Ts3p,T.Tsp)*(T.Tsp-T.Ts3p)+$
   nu_ii(32.0,32.0,3.0,2.0,n.ns3p,n.ns2p,T.Ts3p,T.Ts2p)*(T.Ts2p-T.Ts3p)+$
   nu_ii(32.0,16.0,3.0,1.0,n.ns3p,n.nop,T.Ts3p,T.Top)*(T.Top-T.Ts3p)+$
   nu_ii(32.0,16.0,3.0,2.0,n.ns3p,n.no2p,T.Ts3p,T.To2p)*(T.To2p-T.Ts3p)+$
   nu_ie(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p)*(T.Tel_s3p-T.Ts3p) + $
   nu_ie(32.0,3.0,n.fh*n.nel,n.ns3p,T.Telh,T.Ts3p)*(T.Telh-T.Ts3p) 

Teq_op= nu_ii(16.0,16.0,1.0,2.0,n.nop,n.no2p,T.Top,T.To2p)*(T.To2p-T.Top)+$
   nu_ii(16.0,32.0,1.0,1.0,n.nop,n.nsp,T.Top,T.Tsp)*(T.Tsp-T.Top)+$
   nu_ii(16.0,32.0,1.0,2.0,n.nop,n.ns2p,T.Top,T.Ts2p)*(T.Ts2p-T.Top)+$
   nu_ii(16.0,32.0,1.0,3.0,n.nop,n.ns3p,T.Top,T.Ts3p)*(T.Ts3p-T.Top)+$
   nu_ie(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top)*(T.Tel_op-T.Top) + $
   nu_ie(16.0,1.0,n.fh*n.nel,n.nop,T.Telh,T.Top)*(T.Telh-T.Top) 

Teq_o2p= nu_ii(16.0,16.0,2.0,1.0,n.no2p,n.nop,T.To2p,T.Top)*(T.Top-T.To2p)+$
   nu_ii(16.0,32.0,2.0,1.0,n.no2p,n.nsp,T.To2p,T.Tsp)*(T.Tsp-T.To2p)+$
   nu_ii(16.0,32.0,2.0,2.0,n.no2p,n.ns2p,T.To2p,T.Ts2p)*(T.Ts2p-T.To2p)+$
   nu_ii(16.0,32.0,2.0,3.0,n.no2p,n.ns3p,T.To2p,T.Ts3p)*(T.Ts3p-T.To2p)+$
   nu_ie(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p)*(T.Tel_o2p-T.To2p) + $
   nu_ie(16.0,2.0,n.fh*n.nel,n.no2p,T.Telh,T.To2p)*(T.Telh-T.To2p) 

Teq_eh = nu_ee(n.fc*n.nel,n.fh*n.nel,T.Tel_el,T.Telh)*(T.Telh - T.Tel_el) + $
 nu_ie(32.0,1.0,n.fh*n.nel,n.nsp,T.Telh,T.Tsp)*(T.Telh-T.Tsp) + $
 nu_ie(32.0,2.0,n.fh*n.nel,n.ns2p,T.Telh,T.Ts2p)*(T.Telh-T.Ts2p) + $
 nu_ie(32.0,3.0,n.fh*n.nel,n.ns3p,T.Telh,T.Ts3p)*(T.Telh-T.Ts3p) + $
 nu_ie(16.0,1.0,n.fh*n.nel,n.nop,T.Telh,T.Top)*(T.Telh-T.Top) + $
 nu_ie(16.0,2.0,n.fh*n.nel,n.no2p,T.Telh,T.To2p)*(T.Telh-T.To2p) 

Ec.P_pu = f*nl.s_tot*T.Tpu_s + f*nl.o_tot*T.Tpu_o
Ec.eh_eq = f*Teq_eh 
Ec.ion_e_eq = f*(nu_ie(32.0,1.0,n.fc*n.nel,n.nsp,T.Tel_sp,T.Tsp)*(T.Tsp-T.Tel_sp) + $
      nu_ie(32.0,2.0,n.fc*n.nel,n.ns2p,T.Tel_s2p,T.Ts2p)*(T.Ts2p-T.Tel_s2p) + $
      nu_ie(32.0,3.0,n.fc*n.nel,n.ns3p,T.Tel_s3p,T.Ts3p)*(T.Ts3p-T.Tel_s3p) + $
      nu_ie(16.0,1.0,n.fc*n.nel,n.nop,T.Tel_op,T.Top)*(T.Top-T.Tel_op) + $
      nu_ie(16.0,2.0,n.fc*n.nel,n.no2p,T.Tel_o2p,T.To2p)*(T.To2p-T.Tel_o2p))

;tot_pu = f*nl.s_tot*T.Tpu_s + f*nl.o_tot*T.Tpu_o
;tot_eh = f*Teq_eh
tot_ee = f*nu_ee(n.fc*n.nel,n.fh*n.nel,T.Tel_el,T.Telh)*(T.Telh - T.Tel_el)

print,'   Pickup input.............',Ec.P_pu
print,'   Hot e equilibration......',Ec.eh_eq
print,'   Hot e cold e equil.......',tot_ee
;print,'   Ion e equilibration......',Ec.ion_e_eq
Ec.P_in = Ec.P_pu + Ec.eh_eq
print,'   P input..................',Ec.P_in
print,' '
Ec.Puv = r.Puv
print,'   Ion/elec equil...........',Ec.ion_e_eq
print,'   Hot elec/elec equil......',tot_ee
print,'   Total elec input.........',Ec.ion_e_eq+tot_ee
print,' '
print,'   Puv......................',Ec.Puv
Ec.Pfast = f*(nl.fast_S_k1 + nl.fast_S_k3 + nl.fast_O_k5 + $
           nl.fast_O_k7 + nl.fast_O_k9) 
print,'   Pfast....................',Ec.Pfast
Ec.Ptrans = f*r.transport*(n.nsp*T.Tsp + n.ns2p*T.Ts2p + $
         n.ns3p*T.Ts3p + n.nop*T.Top + n.no2p*T.To2p + n.nel*T.Tel_el) 
Ec.Ptrans_eh = f*r.transport*n.nel*n.fh*T.Telh
print,'   Ptrans...................',Ec.Ptrans
print,'   Ptrans hot elec..........',Ec.Ptrans_eh
Ec.P_out = Ec.Puv + Ec.Pfast + Ec.Ptrans + Ec.Ptrans_eh
print,'   P output............... .',Ec.P_out

print,' '
print,'   Ein/Eout.................',Ec.P_in/Ec.P_out


return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro get_Peuv,n,T,r
;-----------------------------------------------------------------------  
nel = n.fc*n.nel

restore,'euv_emissions.sav

s2em = s2em/1.6e-12
s3em = s3em/1.6e-12
s4em = s4em/1.6e-12
o2em = o2em/1.6e-12
o3em = o3em/1.6e-12

i=where(abs(temp-T.Tel_sp) eq min(abs(temp-T.Tel_sp)))
j=where(abs(density-nel) eq min(abs(density-nel)))

i = i(0)
j = j(0)

if (T.Tel_sp gt temp(i)) then ip = i+1
if (T.Tel_sp le temp(i)) then begin
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
dT1 = T.Tel_sp - temp(i)
dT2 = temp(ip) - T.Tel_sp
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSII = s2em(i,j)*A4 + s2em(ip,j)*A3 + s2em(i,jp)*A2 + $
          s2em(ip,jp)*A1
L = emisSII*n.nsp

i=where(abs(temp-T.Tel_s2p) eq min(abs(temp-T.Tel_s2p)))
j=where(abs(density-nel) eq min(abs(density-nel)))

i = i(0)
j = j(0)

if (T.Tel_s2p gt temp(i)) then ip = i+1
if (T.Tel_s2p le temp(i)) then begin
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
dT1 = T.Tel_s2p - temp(i)
dT2 = temp(ip) - T.Tel_s2p
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSIII = s3em(i,j)*A4 + s3em(ip,j)*A3 + s3em(i,jp)*A2 + $
          s3em(ip,jp)*A1
L = L+ emisSIII*n.ns2p

i=where(abs(temp-T.Tel_s3p) eq min(abs(temp-T.Tel_s3p)))
j=where(abs(density-nel) eq min(abs(density-nel)))

i = i(0)
j = j(0)

if (T.Tel_s3p gt temp(i)) then ip = i+1
if (T.Tel_s3p le temp(i)) then begin
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
dT1 = T.Tel_s3p - temp(i)
dT2 = temp(ip) - T.Tel_s3p
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisSIV = s4em(i,j)*A4 + s4em(ip,j)*A3 + s4em(i,jp)*A2 + $
          s4em(ip,jp)*A1
L = L + emisSIV*n.ns3p

i=where(abs(temp-T.Tel_op) eq min(abs(temp-T.Tel_op)))
j=where(abs(density-nel) eq min(abs(density-nel)))

i = i(0)
j = j(0)

if (T.Tel_op gt temp(i)) then ip = i+1
if (T.Tel_op le temp(i)) then begin
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
dT1 = T.Tel_op - temp(i)
dT2 = temp(ip) - T.Tel_op
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisOII = o2em(i,j)*A4 + o2em(ip,j)*A3 + o2em(i,jp)*A2 + $
          o2em(ip,jp)*A1
L = L + emisOII*n.nop

i=where(abs(temp-T.Tel_o2p) eq min(abs(temp-T.Tel_o2p)))
j=where(abs(density-nel) eq min(abs(density-nel)))

i = i(0)
j = j(0)

if (T.Tel_o2p gt temp(i)) then ip = i+1
if (T.Tel_o2p le temp(i)) then begin
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
dT1 = T.Tel_o2p - temp(i)
dT2 = temp(ip) - T.Tel_o2p
dn1 = nel - density(j)
dn2 = density(jp) - nel
A1 = dT1*dn1/A
A2 = dT2*dn1/A
A3 = dT1*dn2/A
A4 = dT2*dn2/A

emisOIII = o3em(i,j)*A4 + o3em(ip,j)*A3 + o3em(i,jp)*A2 + $
          o3em(ip,jp)*A1
L = L + emisOIII*n.no2p

r.Peuv = L

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------
;pro cm3_model,n,r,T,src,lss,temps,dens,nl,Ec,mr,peuv
pro cm3_model,mr,peuv
;main program
;cubic centimeter torus chemistry model
;-----------------------------------------------------------------------
@common

mr = fltarr(5)
peuv=0.0

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

Rs = {residuals, s: 1.0, sp: 1.0, s2p: 1.0, s3p: 1.0, o: 1.0, op: 1.0, $
                o2p: 1.0}

Jc = fltarr(7,7)
delta_n = fltarr(7)
Rs_arr = fltarr(7)

RsT = fltarr(6)
RsT(0)=1.0
RsT(1)=1.0
RsT(2)=1.0
RsT(3)=1.0
RsT(4)=1.0
RsT(5)=1.0
JcT = fltarr(6,6)
delta_nT = fltarr(6)

n = {density, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
              no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh, $
              nsp_pu: 0.0, ns2p_pu: 0.0, ns3p_pu: 0.0, nop_pu: 0.0, $
              no2p_pu: 0.0}
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
T1 = {temp_1, Tel: Te0, Telh: Teh0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, Top: Ti0, $
           To2p: Ti0, Tsp_pu: 0.0, Ts2p_pu: 0.0, Ts3p_pu: 0.0, $
           Top_pu: 0.0, To2p_pu: 0.0, Tpu_s: 0.0, Tpu_o: 0.0, $
           Tel_sp: Te0, Tel_s2p: Te0, Tel_s3p: Te0, Tel_op: Te0, $
           Tel_o2p: Te0, Tel_el: Te0}
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
            cx_k15: 0.0, cx_k16: 0.0, emisSII: fltarr(31,41), emisSIII: fltarr(31,41), $
            emisSIV: fltarr(31,41), emisOII: fltarr(31,41), $
            emisOIII: fltarr(31,41), $
            emistemp: fltarr(31), emisden: fltarr(41), $
            Puv: 0.0, psp: 0.0, ps2p: 0.0, ps3p: 0.0, Peuv:0.0, $
            pop: 0.0, po2p: 0.0, Sei: 0.0, Secx: 0.0, Oei: 0.0, Oecx:0.0}

r1 = {rates_1, Tel: Te0, is: 0.0, isp: 0.0, is2p: 0.0, rs3p: 0.0, rs2p: 0.0,$
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
            Puv: 0.0, psp: 0.0, ps2p: 0.0, ps3p: 0.0, Peuv:0.0, $
            pop: 0.0, po2p: 0.0, Sei: 0.0, Secx: 0.0, Oei: 0.0, Oecx:0.0}

nl = {neutral_loss, sei: 0.0, seih: 0.0, oei: 0.0, oeih: 0.0, $
       scx: 0.0, ocx: 0.0, k1: 0.0, k2: 0.0, k3: 0.0, k4: 0.0, k5: 0.0, $
       k6: 0.0, k7: 0.0, k8: 0.0, k9: 0.0, k10:0.0, k11: 0.0, k12: 0.0, $
       k14: 0.0, s_tot: 0.0, o_tot: 0.0, scx_tot: 0.0, ocx_tot: 0.0, $
       sion_tot: 0.0, oion_tot: 0.0, fast_O_k5: 0.0, fast_O_k7: 0.0, $
       fast_S_k1: 0.0, fast_S_k3: 0.0, fast_O_k9: 0.0}

Ec = {energy_conserve, s_ion: 0.0, s_ion_h: 0.0, s_cx: 0.0, o_ion: 0.0, o_ion_h: 0.0, $
       o_cx: 0.0, s_tot_in: 0.0, o_tot_in: 0.0, P_pu: 0.0, eh_eq: 0.0, P_in: 0.0, $
       Puv: 0.0, Pfast: 0.0, Ptrans: 0.0, Ptrans_eh: 0.0, P_out: 0.0, $
       in_out: 0.0, ion_e_eq: 0.0, Peuv: 0.0}

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


ep = 1e-7

i = 0
;for i = 0,20 do begin
;while ((abs(Rs.s) gt ep) and (abs(Rs.sp) gt ep) and (abs(Rs.s2p) gt ep) and $
;       (abs(Rs.s3p) gt ep) and (abs(Rs.o) gt ep) and (abs(Rs.op) gt ep) and $
;       (abs(Rs.o2p) gt ep)) do begin 

while (((abs(RsT(0)) gt ep) and (abs(RsT(1)) gt ep) and (abs(RsT(2)) gt ep) and $
       (abs(RsT(3)) gt ep) and (abs(RsT(4)) gt ep) and (abs(RsT(5)) gt ep)) and (i lt 1000)) do begin 

   get_rates,r,T
   
   for m = 0,1 do begin
      get_density_residuals,Rs,n,Rs_arr,r
      get_density_jacobian,Rs,Jc,r,n
      get_delta_n,Rs_arr,Jc,delta_n
      update_density,n,delta_n
   endfor

;   if (n.nsp lt 0.0) then goto,JUMP
;   if (n.ns2p lt 0.0) then goto,JUMP
;   if (n.ns3p lt 0.0) then goto,JUMP
;   if (n.nop lt 0.0) then goto,JUMP
;   if (n.no2p lt 0.0) then goto,JUMP

   get_RsT_elec,n,T,nT,r,RsT
   get_RsT_sp,n,T,nT,r,RsT
   get_RsT_s2p,n,T,nT,r,RsT
   get_RsT_s3p,n,T,nT,r,RsT
   get_RsT_op,n,T,nT,r,RsT
   get_RsT_o2p,n,T,nT,r,RsT

   get_JcT,n,T,nT,r,RsT,JcT

   get_delta_nT,RsT,JcT,delta_nT

   update_nT,n,T,nT,delta_nT
;   print,'T...',T
   if((T.Tel lt 0.0) or (T.Tsp lt 0.0) or (T.Ts2p lt 0.0) or (T.Ts3p lt 0.0) or (T.Top lt 0.0)  or (T.To2p lt 0.0)) then begin
       print,'error...negative temps...'
   endif 
JUMP:

;print,''
;print,'Residuals...',RsT,ep


;   ;full time step advance
;   n1.ns = n.ns + dt*F_s(n,r,S,L)
;   src(i).s = S
;   lss(i).s = L
;   n1.nsp = n.nsp + dt*F_sp(n,r,S,L)
;   src(i).sp = S
;   lss(i).sp = L
;   n1.ns2p = n.ns2p + dt*F_s2p(n,r,S,L)
;   src(i).s2p = S
;   lss(i).s2p = L
;   n1.ns3p = n.ns3p + dt*F_s3p(n,r,S,L)
;   src(i).s3p = S
;   lss(i).s3p = L
;   n1.no = n.no + dt*F_o(n,r,S,L)
;   src(i).o = S
;   lss(i).o = L
;   n1.nop = n.nop + dt*F_op(n,r,S,L)
;   src(i).op = S
;   lss(i).op = L 
;   n1.no2p = n.no2p + dt*F_o2p(n,r,S,L)
;   src(i).o2p = S
;   lss(i).o2p = L
;   src(i).tm = dt*i/8.64e4
;   lss(i).tm = dt*i/8.64e4

;   n1.nel = 1.0*n1.nsp+2.0*n1.ns2p+3.0*n1.ns3p+1.0*n1.nop+2.0*n1.no2p

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

;   ;Improved Euler advance

;   np.ns = n.ns + dt*0.5*(F_s(n,r)+F_s(n1,r))
;   np.nsp = n.nsp + dt*0.5*(F_sp(n,r)+F_sp(n1,r))
;   np.ns2p = n.ns2p + dt*0.5*(F_s2p(n,r)+F_s2p(n1,r))
;   np.ns3p = n.ns3p + dt*0.5*(F_s3p(n,r)+F_s3p(n1,r))
;   np.no = n.no + dt*0.5*(F_o(n,r)+F_o(n1,r))
;   np.nop = n.nop + dt*0.5*(F_op(n,r)+F_op(n1,r))
;   np.no2p = n.no2p + dt*0.5*(F_o2p(n,r)+F_o2p(n1,r))

;   np.nel = 1.0*np.nsp+2.0*np.ns2p+3.0*np.ns3p+1.0*np.nop+2.0*np.no2p


;   nTp.sp = nT.sp + dt*0.5*(EF_sp(n1,T1,nT1,r)+ EF_sp(n,T,nT,r))
;   nTp.s2p = nT.s2p + dt*0.5*(EF_s2p(n1,T1,nT1,r)+ EF_s2p(n,T,nT,r))
;   nTp.s3p = nT.s3p + dt*0.5*(EF_s3p(n1,T1,nT1,r)+ EF_s3p(n,T,nT,r))
;   nTp.op = nT.op + dt*0.5*(EF_op(n1,T1,nT1,r)+ EF_op(n,T,nT,r))
;   nTp.o2p = nT.o2p + dt*0.5*(EF_o2p(n1,T1,nT1,r)+ EF_o2p(n,T,nT,r))
;   nTp.el = nT.el + dt*0.5*(EF_el(n1,T1,nT1,r)+ EF_el(n,T,nT,r))
;   update_temp,np,nTp,Tp
;   if (ft_ave eq 1) then cm3_ekappa,np,Tp
;   if (ft_ave eq 0) then begin
;      Tp.Tel_sp = Tp.Tel
;      Tp.Tel_s2p = Tp.Tel
;      Tp.Tel_s3p = Tp.Tel
;      Tp.Tel_op = Tp.Tel
;      Tp.Tel_o2p = Tp.Tel
;      Tp.Tel_el = Tp.Tel
;   endif

;   n = np
;   nT = nTp  
;   T = Tp
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

;;   print,T.Tsp,T.Ts2p,T.Ts3p,T.Top,T.To2p,T.Tel

;;update Pickup Densities
  
;;   n.nsp_pu = nT.sp_pu/T.Tsp_pu
;;   n.ns2p_pu = nT.s2p_pu/T.Ts2p_pu
;;   n.ns3p_pu = nT.s3p_pu/T.Ts3p_pu
;;   n.nop_pu = nT.op_pu/T.Top_pu
;;   n.no2p_pu = nT.o2p_pu/T.To2p_pu

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

   i = i + 1
;endfor
endwhile

;get_neutral_loss_rate,n,r,T,nl
energy_conservation,n,r,T,nl,Ec
get_Peuv,n,T,r

;save,filename='src_func.sav',src
;save,filename='lss_func.sav',lss
;save,filename='temps.sav',temps
;save,filename='dens.sav',dens
;save,filename='pwr.sav',p

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
print,'ns,no............',n.ns,n.no
;print,'Production of high speed O...',nl.k5
;print,'Production of neutral O......',r.S_production
;print,'High speed O/low speed O.....',nl.k5/r.S_production
nelec = n.nel + 0.1*n.nel
print,' '
print,'nsp....',n.nsp/nelec, 0.12, 0.092
print,'ns2p...',n.ns2p/nelec, 0.15, 0.17
print,'ns3p...',n.ns3p/nelec, 0, 0.023
print,'nop....',n.nop/nelec, 0.36, 0.47
print,'no2p...',n.no2p/nelec, 0.12, 0.016

mr(0) = n.nsp/nelec
mr(1) = n.ns2p/nelec
mr(2) = n.ns3p/nelec
mr(3) = n.nop/nelec
mr(4) = n.no2p/nelec
peuv = r.Peuv

print,' '
print,'S/O:     ',n.ns/n.no
print,'S+/O+:   ',n.nsp/n.nop
print,'S+/S++:  ',n.nsp/n.ns2p
print,'S+++/S++: ',n.ns3p/n.ns2p
print,'O++/O+:  ',n.no2p/n.nop
print,' '
print,'nel.....',n.nel
print,' '
neni = n.nel/(n.nsp+n.ns2p+n.ns3p+n.nop+n.no2p)
sonp = n.nop+n.no2p
ssnp = n.nsp+n.ns2p+n.ns3p
print,'ne/ni....',neni
print,'On+/Sn+..',sonp/ssnp
print,' '
print,'Tsp......',T.Tsp
print,'Ts2p.....',T.Ts2p
print,'Ts3p.....',T.Ts3p
print,'Top......',T.Top
print,'To2p.....',T.To2p
print,'Tel......',T.Tel
sti = T.Tsp*n.nsp+T.Ts2p*n.ns2p+T.Ts3p*n.ns3p+T.Top*n.nop+T.To2p*n.no2p
sni = n.nsp+n.ns2p+n.ns3p+n.nop+n.no2p
print,'<Ti>.....',sti/sni

spne1 = n.nsp/n.nel
s2pne1 = n.ns2p/n.nel
s3pne1 = n.ns3p/n.nel
opne1 = n.nop/n.nel
o2pne1 = n.no2p/n.nel
neni1 = neni
onsn1 = sonp/ssnp
Tiave1 = sti/sni
Tel1 = T.Tel
nel1 = n.nel
save,filename='sensitivity.sav',spne1,s2pne1,s3pne1,opne1,o2pne1,$
  neni1,onsn1,Tel1,Tiave1,nel1

print,' '
print,'Fluxtube averaged electron temp....................',T.Tel_el
print,'Fluxtube averaged electron temp, (O+ weighted).....',T.Tel_op
print,'Fluxtube averaged electron temp, (O++ weighted)....',T.Tel_o2p
print,'Fluxtube averaged electron temp, (S+ weighted).....',T.Tel_sp
print,'Fluxtube averaged electron temp, (S++ weighted)....',T.Tel_s2p
print,'Fluxtube averaged electron temp, (S+++ weighted)...',T.Tel_s3p
print,' '
save,filename='restart_cm3.sav',n,n1,np,T,T1,Tp,nT,nT1,nTp,Ec

;cm3_ekappa,n,T

return
end
;-----------------------------------------------------------------------




