;-----------------------------------------------------------------------
pro get_rates,r
;-----------------------------------------------------------------------

;Electron impact ionization (Schreier and Eviatar 1998)
r.is = 7.4e-8*sqrt(r.Te)*exp(-10.0/r.Te)     ;check this!!! 
r.isp = 1.274e-8*sqrt(r.Te)*exp(-23.1/r.Te)  
r.is2p = 5.8e-9*sqrt(r.Te)*exp(-34.88/r.Te)  

r.io = 1.015e-8*sqrt(r.Te)*exp(-14.54/r.Te)
r.iop = 3.78e-9*sqrt(r.Te)*exp(-33.42/r.Te)

;dielectron recombination (Smith and Strobel 1985)
if (r.Te le 6.9) then r.rs3p = 5.4e-11*r.Te^0.53
if (r.Te gt 6.9) then r.rs3p = 4.3e-12*r.Te^1.84  

r.rs2p = r.rs3p/3.0   ;massive fudge based on Johnson and Strobel '82

r.ro2p = r.rs3p/7.0   ;another massive fudge base on Johnson and Stobel

;charge exchange

;r.cx_o2p_sp = 1.4e-9   ;Schreier for L=6
r.cx_o2p_sp = 0.0      ;Schreier for L=6
;r.cx_o2p_s2p = 9e-9    ;Schreier for L=6
r.cx_o2p_s2p = 0.0    ;Schreier for L=6

Rj = 7.14e4*1e5 ;cm
a = 5.0*Rj
b = 8.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2

r.S_production = 6.5e27/vol
r.O_production = 2.0*r.S_production

return
end
;-----------------------------------------------------------------------



;-----------------------------------------------------------------------
pro update_s,n,np,r,dt
;-----------------------------------------------------------------------
S = r.S_production
L = r.is*n.ns*n.nel

np.ns = n.ns + dt*(S - L)

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro update_sp,n,np,r,dt
;-----------------------------------------------------------------------  
S = r.is*n.ns*n.nel + r.rs2p*n.ns2p*n.nel
L = r.isp*n.nsp*n.nel + r.cx_o2p_sp*n.no2p*n.nsp + n.nsp*r.Transport

np.nsp = n.nsp + dt*(S - L)

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro update_s2p,n,np,r,dt
;-----------------------------------------------------------------------  
S = r.isp*n.nsp*n.nel + r.rs3p*n.ns3p*n.nel + r.cx_o2p_sp*n.no2p*n.nsp
L = r.is2p*n.ns2p*n.nel + r.rs2p*n.ns2p*n.nel + $
    r.cx_o2p_s2p*n.no2p*n.ns2p + n.ns2p*r.Transport

np.ns2p = n.ns2p + dt*(S - L)

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro update_s3p,n,np,r,dt
;-----------------------------------------------------------------------  
S = r.is2p*n.ns2p*n.nel + r.cx_o2p_s2p*n.no2p*n.ns2p
L = r.rs3p*n.ns3p*n.nel + n.ns3p*r.Transport

np.ns3p = n.ns3p + dt*(S - L)

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------
pro update_o,n,np,r,dt
;-----------------------------------------------------------------------
S = r.O_production
L = r.io*n.no*n.nel

np.no = n.no + dt*(S - L)

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro update_op,n,np,r,dt
;-----------------------------------------------------------------------  
S = r.io*n.no*n.nel + r.ro2p*n.no2p*n.nel + r.cx_o2p_sp*n.no2p*n.nsp +$
    r.cx_o2p_s2p*n.no2p*n.ns2p
L = r.iop*n.nop*n.nel + n.nop*r.Transport

np.nop = n.nop + dt*(S - L)

return
end
;-----------------------------------------------------------------------  


;-----------------------------------------------------------------------  
pro update_o2p,n,np,r,dt
;-----------------------------------------------------------------------  
S = r.iop*n.nop*n.nel
L = r.ro2p*n.no2p*n.nel + r.cx_o2p_sp*n.no2p*n.nsp + $
    r.cx_o2p_s2p*n.no2p*n.ns2p + n.no2p*r.Transport

np.no2p = n.no2p + dt*(S - L)

return
end
;-----------------------------------------------------------------------  



;-----------------------------------------------------------------------
pro cm3_model
;main program
;cubic centimeter torus chemistry model
;-----------------------------------------------------------------------
@common

ns0 = 1.0      ;initial S neutral density (cm^-3)
no0 = 2.0*ns0     ;initial O neutral density
nsp0= 1.0        ;initial S+ density
ns2p0 = 1.0       ;initial S++ density
ns3p0 = 1.0      ;initial S+++ density
nop0 = 1.0       ;initial O+ density
no2p0 = 1.0      ;initial O++ density 
nel0 = nsp0+ns2p0+ns3p0+nop0+no2p0       ;initial electron density
;Te0 = 4.0         ;initial electron temperature (eV)
dt = 10000.0          ;time step (seconds)
nt = 3200.

n = {density, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
              no: no0, nop: nop0, no2p: no2p0}
np = {density_p, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
                 no: no0, nop: nop0, no2p: no2p0}
r = {rates, Te: Te0, is: 0.0, isp: 0.0, is2p: 0.0, rs3p: 0.0, rs2p: 0.0,$
            io: 0.0, iop:0.0, ro2p:0.0, cx_o2p_sp: 0.0, cx_o2p_s2p: 0.0, $
            S_production: 0.0, O_production: 0.0, Transport: trans}

plot_io,[0,nt*dt]/3600.,[1.0,no0+ns0],/nodata,$
  title='T!de!n = '+strtrim(string(Te0),2)+' (eV)',$
  xtitle='time (hours) ',ytitle='n (cm!u-3!n)',yrange=[1.0,2000.0],/ysty,$
  xrange=[0,dt*nt]/3600.,/xsty
  
device,decompose=0
loadct,27
for i = 0,nt do begin

   get_rates,r
   update_s,n,np,r,dt
   update_sp,n,np,r,dt
   update_s2p,n,np,r,dt
   update_s3p,n,np,r,dt
   update_o,n,np,r,dt
   update_op,n,np,r,dt
   update_o2p,n,np,r,dt
   n = np
   n.nel = n.nsp+n.ns2p+n.ns3p+n.nop+n.no2p
;   plots,i*dt/3600.,n.ns,psym=3,color=10
   plots,i*dt/3600.,n.nsp,psym=3,color=10,thick=2.5
   plots,i*dt/3600.,n.ns2p,psym=3,color=75,thick=2.5
   plots,i*dt/3600.,n.ns3p,psym=3,color=50,thick=2.5
   plots,i*dt/3600.,n.nel,psym=3,color=100,thick=2.5
;   plots,i*dt/3600.,n.no,psym=3,color=75,thick=2.5
   plots,i*dt/3600.,n.nop,psym=3,color=150,thick=2.5
   plots,i*dt/3600.,n.no2p,psym=3,color=200,thick=2.5

endfor

legend,['S!u+!n','S!u++!n','S!u+++!n','O!u+!n','O!u++!n','n!de!n'],$
       linestyle=[0,0,0,0,0,0],colors=[10,75,50,150,200,100],/bottom,$
       /right

return
end
;-----------------------------------------------------------------------



