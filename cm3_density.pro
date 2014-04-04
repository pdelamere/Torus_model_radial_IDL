;-----------------------------------------------------------------------
pro get_rates,r
;-----------------------------------------------------------------------
@common

;Electron impact ionization (Schreier and Eviatar 1998)
r.is = 7.4e-8*sqrt(r.Te)*exp(-10.36/r.Te)     ;check this!!! 
r.isp = 1.274e-8*sqrt(r.Te)*exp(-23.1/r.Te)  
r.is2p = 5.8e-9*sqrt(r.Te)*exp(-34.88/r.Te)  

r.ish = 1.464e-6*alog(r.Teh/27.37)/sqrt(r.Teh)
r.isph = 4.44e-7*alog(r.Teh/27.07)/sqrt(r.Teh)
r.is2ph = 1.76e-7*alog(r.Teh/27.07)/sqrt(r.Teh)

r.io = 1.015e-8*sqrt(r.Te)*exp(-14.54/r.Te)
r.iop = 3.78e-9*sqrt(r.Te)*exp(-33.42/r.Te)

r.ioh = 5.74e-7*alog(r.Teh/27.07)/sqrt(r.Teh)
r.ioph = 1.76e-7*alog(r.Teh/27.07)/sqrt(r.Teh)

;dielectron recombination (Smith and Strobel 1985)
if (r.Te le 6.9) then r.rs3p = 5.4e-11*r.Te^0.53
if (r.Te gt 6.9) then r.rs3p = 4.3e-12*r.Te^1.84  

r.rs2p = r.rs3p/3.0   ;massive fudge based on Johnson and Strobel '82

r.ro2p = r.rs3p/7.0   ;another massive fudge base on Johnson and Stobel

;charge exchange

r.cx_s_s3p = k4*1.32e-8
r.cx_s_o2p = k10*2.34e-8
r.cx_s_o2p_e = k11*1.62e-8
r.cx_o2p_sp = k13*1.4e-9   ;Schreier for L=6
r.cx_o_s3p = k14*1.92e-8
r.cx_o2p_s2p = k15*9e-9    ;Schreier for L=6

print,k4,k10,k11,k13,k14,k15

Rj = 7.14e4*1e5 ;cm
a = 5.0*Rj
b = 8.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2

fs = 1.0/(1.0+otos)
fo = otos/(1.0+otos)
print,fs,fo

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
    r.cx_s_s3p*n.ns*n.ns3p + r.cx_s_o2p_e*n.ns*n.no2p + $
    r.cx_s_o2p*n.ns*n.no2p

F_s = S - L

return,F_s
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_sp,n,r,S,L
;-----------------------------------------------------------------------  
S = r.is*n.ns*n.fc*n.nel + r.ish*n.ns*n.fh*n.nel + r.rs2p*n.ns2p*n.nel + $
    r.cx_s_s3p*n.ns*n.ns3p + r.cx_s_o2p*n.ns*n.no2p
L = r.isp*n.nsp*n.fc*n.nel + r.cx_o2p_sp*n.no2p*n.nsp + n.nsp*r.Transport + $
    r.isph*n.nsp*n.fh*n.nel

F_sp = S - L

return,F_sp
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_s2p,n,r,S,L
;-----------------------------------------------------------------------  
S = r.isp*n.nsp*n.fc*n.nel + r.isph*n.nsp*n.fh*n.nel + $
    r.rs3p*n.ns3p*n.nel + r.cx_o2p_sp*n.no2p*n.nsp + $
    r.cx_o_s3p*n.no*n.ns3p + r.cx_s_s3p*n.ns*n.ns3p + $
    r.cx_s_o2p_e*n.ns*n.no2p
L = r.is2p*n.ns2p*n.fc*n.nel + r.rs2p*n.ns2p*n.nel + $
    r.cx_o2p_s2p*n.no2p*n.ns2p + n.ns2p*r.Transport + $
    r.is2ph*n.ns2p*n.fh*n.nel

F_s2p = S - L

return,F_s2p
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_s3p,n,r,S,L
;-----------------------------------------------------------------------  
S = r.is2p*n.ns2p*n.fc*n.nel + r.is2ph*n.ns2p*n.fh*n.nel + $
    r.cx_o2p_s2p*n.no2p*n.ns2p
L = r.rs3p*n.ns3p*n.nel + n.ns3p*r.Transport + r.cx_o_s3p*n.no*n.ns3p + $
    r.cx_s_s3p*n.ns*n.ns3p

F_s3p = S - L

return,F_s3p
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------
function F_o,n,r,S,L
;-----------------------------------------------------------------------
S = r.O_production
L = r.io*n.no*n.fc*n.nel + r.ioh*n.no*n.fh*n.nel + r.cx_o_s3p*n.no*n.ns3p

F_o = S - L

return,F_o
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_op,n,r,S,L
;-----------------------------------------------------------------------  
S = r.io*n.no*n.fc*n.nel + r.ioh*n.no*n.fh*n.nel + $
    r.ro2p*n.no2p*n.nel + r.cx_o2p_sp*n.no2p*n.nsp + $
    r.cx_o2p_s2p*n.no2p*n.ns2p + r.cx_o_s3p*n.no*n.ns3p + $
    r.cx_s_o2p_e*n.ns*n.no2p + r.cx_s_o2p*n.ns*n.no2p
L = r.iop*n.nop*n.fc*n.nel + n.nop*r.Transport + r.ioph*n.nop*n.fh*n.nel

F_op = S - L

return,F_op
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
function F_o2p,n,r,S,L
;-----------------------------------------------------------------------  
S = r.iop*n.nop*n.fc*n.nel + r.ioph*n.nop*n.fh*n.nel
L = r.ro2p*n.no2p*n.nel + r.cx_o2p_sp*n.no2p*n.nsp + $
    r.cx_o2p_s2p*n.no2p*n.ns2p + n.no2p*r.Transport + $
    r.cx_s_o2p_e*n.ns*n.no2p + r.cx_s_o2p*n.ns*n.no2p

F_o2p = S - L

return,F_o2p
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------
pro cm3_model,n
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
;Te0 = 4.0         ;initial electron temperature (eV)
dt = 10000.0          ;time step (seconds)
nt = fix(runt/dt)

fc = 1.0-fh

n = {density, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
              no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh}

;advance density for improved Euler method, half time step
n1 = {density_1, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
                 no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh}
;density at full time step advance
np = {density_p, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
                 no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh}
r = {rates, Te: Te0, is: 0.0, isp: 0.0, is2p: 0.0, rs3p: 0.0, rs2p: 0.0,$
            io: 0.0, iop:0.0, ro2p:0.0, cx_o2p_sp: 0.0, cx_o2p_s2p: 0.0, $
            Teh: Teh0,ish: 0.0, isph: 0.0, is2ph: 0.0, ioh: 0.0, ioph: 0.0, $
            S_production: 0.0, O_production: 0.0, net_production: 0.0, $
            Transport: trans, $
            cx_o_s3p: 0.0, cx_s_s3p: 0.0, cx_s_o2p_e: 0.0, cx_s_o2p: 0.0,$
            o_to_s: otos}

src = replicate({source_func, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0,$
              o: 0.0, op: 0.0, o2p: 0.0,tm:0.0},nt+1)
lss = replicate({loss_func, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0,$
              o: 0.0, op: 0.0, o2p: 0.0,tm:0.0},nt+1)

plot_io,[0,nt*dt]/8.64e4,[1.0,no0+ns0],/nodata,$
  title='T!de!n = '+strtrim(string(Te0),2)+' (eV)',$
  xtitle='time (days) ',ytitle='n (cm!u-3!n)',yrange=[0.1,10000.0],/ysty,$
  xrange=[0,dt*nt]/8.64e4,/xsty

help,n.nsp,n.nop
  
print,!d.name
if (strtrim(!d.name,2) ne 'PS') then device,decompose=0
loadct,27
get_rates,r

if (cont eq 'true') then begin
   restore,'restart_cm3.sav'
   fc = 1.0 - fh
   n.fc = fc
   n.fh = fh
   n1.fc = fc
   n1.fh = fh
   np.fc = fc
   np.fh = fh
endif


for i = 0,nt do begin

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

   ;Improved Euler advance

   np.ns = n.ns + dt*0.5*(F_s(n,r)+F_s(n1,r))
   np.nsp = n.nsp + dt*0.5*(F_sp(n,r)+F_sp(n1,r))
   np.ns2p = n.ns2p + dt*0.5*(F_s2p(n,r)+F_s2p(n1,r))
   np.ns3p = n.ns3p + dt*0.5*(F_s3p(n,r)+F_s3p(n1,r))
   np.no = n.no + dt*0.5*(F_o(n,r)+F_o(n1,r))
   np.nop = n.nop + dt*0.5*(F_op(n,r)+F_op(n1,r))
   np.no2p = n.no2p + dt*0.5*(F_o2p(n,r)+F_o2p(n1,r))

   n = np
  
   n.nel = n.nsp+n.ns2p+n.ns3p+n.nop+n.no2p

;   plots,i*dt/8.64e4,n.ns,psym=3,color=10
   plots,i*dt/8.64e4,n.nsp,psym=3,color=10,thick=2.5
   plots,i*dt/8.64e4,n.ns2p,psym=3,color=75,thick=2.5
   plots,i*dt/8.64e4,n.ns3p,psym=3,color=50,thick=2.5
   plots,i*dt/8.64e4,n.nel,psym=3,color=100,thick=2.5
;   plots,i*dt/8.64e4,n.no,psym=3,color=75,thick=2.5
   plots,i*dt/8.64e4,n.nop,psym=3,color=150,thick=2.5
   plots,i*dt/8.64e4,n.no2p,psym=3,color=200,thick=2.5

endfor

save,filename='src_func.sav',src
save,filename='lss_func.sav',lss

if (lonoff eq 1.0) then begin
legend,['S!u+!n','S!u++!n','S!u+++!n','O!u+!n','O!u++!n','n!de!n'],$
       linestyle=[0,0,0,0,0,0],colors=[10,75,50,150,200,100],/bottom,$
       /right
endif

upsp = n.nsp
ups2p = n.ns2p
ups3p = n.ns3p
upop = n.nop
upo2p = n.no2p

print,'nsp....',n.nsp
print,'ns2p...',n.ns2p
print,'ns3p...',n.ns3p
print,'nop....',n.nop
print,'no2p...',n.no2p
print,'S/O:     ',n.ns/n.no
print,'S+/O+:   ',n.nsp/n.nop
print,'S++/S+:  ',n.ns2p/n.nsp
print,'S+++/S+: ',n.ns3p/n.nsp

save,filename='restart_cm3.sav',n,n1,np

return
end
;-----------------------------------------------------------------------




