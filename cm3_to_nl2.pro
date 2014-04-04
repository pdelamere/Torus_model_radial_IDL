;------------------------------------------------------------------
pro get_NL2,nr,tr,nl2,nl2e,nar
;------------------------------------------------------------------
@common

nar = replicate({dens_lat, nsp: 0.0, ns2p: 0.0, ns3p: 0.0, ns4p: 0.0, nop: 0.0, $
                no2p: 0.0, nel: 0.0, nelh: 0.0, ns: 0.0, no: 0.0},61)

theta = max_theta-findgen(2.0*max_theta)
dtheta = 1.0


for ii = 0,nbox-1 do begin

   i = ii+bufcel

   n = nr(i)
   T = tr(i)
   cm3_expand,n,narr,T,nar,max_theta

   nl2(i).nsp = double(8e28)*Lshell(i)^4*total(nar(*).nsp*cos(theta*!dtor)^7*dtheta)
   nl2(i).ns2p = double(8e28)*Lshell(i)^4*total(nar(*).ns2p*cos(theta*!dtor)^7*dtheta)
   nl2(i).ns3p = double(8e28)*Lshell(i)^4*total(nar(*).ns3p*cos(theta*!dtor)^7*dtheta)
   nl2(i).ns4p = double(8e28)*Lshell(i)^4*total(nar(*).ns4p*cos(theta*!dtor)^7*dtheta)
   nl2(i).nop = double(8e28)*Lshell(i)^4*total(nar(*).nop*cos(theta*!dtor)^7*dtheta)
   nl2(i).no2p = double(8e28)*Lshell(i)^4*total(nar(*).no2p*cos(theta*!dtor)^7*dtheta)


;  Schreier
;   nl2e(i).nsp = nl2(i).nsp*T.Tsp*Lshell(i)^(8./3.)
;   nl2e(i).ns2p = nl2(i).ns2p*T.Ts2p*Lshell(i)^(8./3.)
;   nl2e(i).ns3p = nl2(i).ns3p*T.Ts3p*Lshell(i)^(8./3.)
;   nl2e(i).ns4p = nl2(i).ns4p*T.Ts4p*Lshell(i)^(8./3.)
;   nl2e(i).nop = nl2(i).nop*T.Top*Lshell(i)^(8./3.)
;   nl2e(i).no2p = nl2(i).no2p*T.To2p*Lshell(i)^(8./3.)

;  Richardson
   nl2e(i).nsp = nl2(i).nsp*T.Tsp^(4./3.)*Lshell(i)^2
   nl2e(i).ns2p = nl2(i).ns2p*T.Ts2p^(4./3.)*Lshell(i)^2
   nl2e(i).ns3p = nl2(i).ns3p*T.Ts3p^(4./3.)*Lshell(i)^2
   nl2e(i).ns4p = nl2(i).ns4p*T.Ts4p^(4./3.)*Lshell(i)^2
   nl2e(i).nop = nl2(i).nop*T.Top^(4./3.)*Lshell(i)^2
   nl2e(i).no2p = nl2(i).no2p*T.To2p^(4./3.)*Lshell(i)^2

endfor



return
end
;------------------------------------------------------------------


;------------------------------------------------------------------
pro transport_NL2,nl2,nl2e,DLL_0,DLL_alpha               
;pro transport_NL2,nl2,nl2e
;------------------------------------------------------------------
@common
;@params

tbndry_out = 200.0 ;eV
tbndry_in = 70.0
Lo = 6.0
L = [Lshell(0)-dL(0),Lshell,Lshell(nintrvl-1)+dL(nintrvl-1)]
dL2 = (shift(L,-1)-L)/2
Lp = L + dL2

dll = fltarr(n_elements(L))
dll = DLL_0*(Lp/Lo)^DLL_alpha

L_flx = (shift(L,-1) + L)/2.0
L_flx(n_elements(L)-1) = L_flx(n_elements(L)-2) + dL(n_elements(dL)-1)
nL = n_elements(L)
whsm = where(L ge 9.0)

nl2b= [nl2(0).nsp,nl2.nsp,0]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
flux = (dll/Lp^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
nl2.nsp = nl2b(1:nintrvl)

;energy
;solve for boundary temp
;Schreier 
;tspb = tbndry_out*nl2(nintrvl-1).nsp*L(nintrvl)^(8./3.)
;tspb1 = tbndry_in*nl2(0).nsp*L(0)^(8./3)

;Richardson
tspb = tbndry_out^(4./3.)*nl2(nintrvl-1).nsp*L(nintrvl)^(2.)
tspb1 = tbndry_in^(4./3.)*nl2(0).nsp*L(0)^(2.)

nl2b= [tspb1,nl2e.nsp,tspb]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
flux = (dll/Lp^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
nl2e.nsp = nl2b(1:nintrvl)

nl2b = [nl2(0).ns2p,nl2.ns2p,0]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
flux = (dll/Lp^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
nl2.ns2p = nl2b(1:nintrvl)

;energy
;solve for boundary temp
;ts2pb = tbndry_out*nl2(nintrvl-1).ns2p*L(nintrvl)^(8./3.)
;ts2pb1 = tbndry_in*nl2(0).ns2p*L(0)^(8./3)

;Richardson
ts2pb = tbndry_out^(4./3.)*nl2(nintrvl-1).ns2p*L(nintrvl)^(2.)
ts2pb1 = tbndry_in^(4./3.)*nl2(0).ns2p*L(0)^(2.)


nl2b = [ts2pb1,nl2e.ns2p,ts2pb]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
flux = (dll/Lp^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
nl2e.ns2p = nl2b(1:nintrvl)


nl2b = [nl2(0).ns3p,nl2.ns3p,0]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
flux = (dll/Lp^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
nl2.ns3p = nl2b(1:nintrvl)

;energy
;solve for boundary temp
;ts3pb = tbndry_out*nl2(nintrvl-1).ns3p*L(nintrvl)^(8./3.)
;ts3pb1 = tbndry_in*nl2(0).ns3p*L(0)^(8./3)

;Richardson
ts3pb = tbndry_out^(4./3.)*nl2(nintrvl-1).ns3p*L(nintrvl)^(2.)
ts3pb1 = tbndry_in^(4./3.)*nl2(0).ns3p*L(0)^(2.)


nl2b = [ts3pb1,nl2e.ns3p,ts3pb]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
flux = (dll/Lp^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
nl2e.ns3p = nl2b(1:nintrvl)


nl2b = [nl2(0).ns4p,nl2.ns4p,0]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
flux = (dll/Lp^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
nl2.ns4p = nl2b(1:nintrvl)

;energy
;solve for boundary temp
topb = tbndry_out^(4./3.)*nl2(nintrvl-1).nop*L(nintrvl)^2

nl2b = [0,nl2e.ns4p,topb]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
flux = (dll/Lp^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
nl2e.ns4p = nl2b(1:nintrvl)

nl2b = [nl2(0).nop,nl2.nop,0]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
flux = (dll/Lp^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
nl2.nop = nl2b(1:nintrvl)

;energy
;solve for boundary temp
;topb = tbndry_out*nl2(nintrvl-1).nop*L(nintrvl)^(8./3.)
;topb1 = tbndry_in*nl2(0).nop*L(0)^(8./3)

;Richardson
topb = tbndry_out^(4./3.)*nl2(nintrvl-1).nop*L(nintrvl)^(2.)
topb1 = tbndry_in^(4./3.)*nl2(0).nop*L(0)^(2.)


nl2b = [topb1,nl2e.nop,topb]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
flux = (dll/Lp^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
nl2e.nop = nl2b(1:nintrvl)

nl2b = [nl2(0).no2p,nl2.no2p,0]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
flux = (dll/Lp^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
nl2.no2p = nl2b(1:nintrvl)

;energy
;solve for boundary temp
;to2pb = tbndry_out*nl2(nintrvl-1).no2p*L(nintrvl)^(8./3.)
;to2pb1 = tbndry_in*nl2(0).no2p*L(0)^(8./3)

;Richardson
to2pb = tbndry_out^(4./3.)*nl2(nintrvl-1).no2p*L(nintrvl)^(2.)
to2pb1 = tbndry_in^(4./3.)*nl2(0).no2p*L(0)^(2.)

nl2b = [to2pb1,nl2e.no2p,to2pb]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
flux = (dll/Lp^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
nl2e.no2p = nl2b(1:nintrvl)

return
end
;------------------------------------------------------------------


;------------------------------------------------------------------
pro transport_NL2_td,nl2,nl2e,runtime,ft1
;------------------------------------------------------------------
@common
@params

tbndry = 2e2
Lo = 6.0
L = [Lshell(0)-dL(0),Lshell,Lshell(nintrvl-1)+dL(nintrvl-1)]
dll = fltarr(n_elements(L))
whlt = where(L lt 6.0)
whgt = where(L ge 6.0)
dll(whgt) = DLL_0*(L(whgt)/Lo)^DLL_alpha
dll(whlt) = 9.92e-8*(L(whlt)/Lo)^2
;dll(1) = dll(1) + 2.0*dll(1)*ft1
;print,'Dll 0 ...',dll(1),(1./dll(1))/60./60./24.,ft1


L_flx = (shift(L,-1) + L)/2.0
L_flx(n_elements(L)-1) = L_flx(n_elements(L)-2) + dL(n_elements(dL)-1)
nL = n_elements(L)
whsm = where(L ge 9.0)

;nl2b= [0,nl2.nsp,0]
nl2b= [nl2(0).nsp,nl2.nsp,0]
;nl2b = [nl2(0).nsp,nl2.nsp,nl2(nintrvl-1).nsp]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
flux = (dll/L^2)*gradnl2
flxav = (flux + shift(flux,1))/2.0
flxav(0) = flux(0)
flxav(nL-1) = flux(nL-2)
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
;plot,flux(0:nintrvl-2)
;oplot,flxav(0:nintrvl-2)
;wait,0.01
;upwind
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))

;Lax-Wendroff
;dL_flx = shift(L,-1)-shift(L,1)
;nl2b = nl2b + L^2*(dt_trans/dL_flx)*(shift(flxav,-1)-shift(flxav,1)) +  $
;       0.1*(dt_trans/dL_flx)^2*((shift(flxav,-1)-flxav)/(shift(L,-1)-L)-$
;           (flxav - shift(flxav,1))/(L-shift(l,1)))


;nl2b(whsm) = smooth(nl2b(whsm),2)
nl2.nsp = nl2b(1:nintrvl)

;energy
;solve for boundary temp
tspb = tbndry*nl2(nintrvl-1).nsp*L(nintrvl)^(8./3.)*L(nintrvl)^2
;tspb = tbndry*nl2(nintrvl-1).nsp*L(nintrvl)^(8./3.)

nl2b= [0,nl2e.nsp,tspb]
;nl2b= [nl2e(0).nsp,nl2e.nsp,0]
;nl2b = [nl2e(0).nsp,nl2e.nsp,nl2e(nintrvl-1).nsp]
;nl2b = [0,nl2e.nsp,nl2e(nintrvl-1).nsp]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
;gradnl2 = (shift(nl2b,-1) - nl2b)/dL
flux = (dll/L^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
;nl2b = nl2b + L^2*(dt_trans/dL)*(flux-shift(flux,1))
;nl2b = smooth(nl2b,2)
;nl2b(whsm) = smooth(nl2b(whsm),2)
nl2e.nsp = nl2b(1:nintrvl)

;nl2b = [0,nl2.ns2p,0]
nl2b = [nl2(0).ns2p,nl2.ns2p,0]
;nl2b = [nl2(0).ns2p,nl2.ns2p,nl2(nintrvl-1).ns2p]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
;gradnl2 = (shift(nl2b,-1) - nl2b)/dL
flux = (dll/L^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
;nl2b = nl2b + L^2*(dt_trans/dL)*(flux-shift(flux,1))
;nl2b = smooth(nl2b,2)
;nl2b(whsm) = smooth(nl2b(whsm),2)
nl2.ns2p = nl2b(1:nintrvl)

;energy
;solve for boundary temp
ts2pb = tbndry*nl2(nintrvl-1).ns2p*L(nintrvl)^(8./3.)*L(nintrvl)^2
;ts2pb = tbndry*nl2(nintrvl-1).ns2p*L(nintrvl)^(8./3.)

nl2b = [0,nl2e.ns2p,ts2pb]
;nl2b = [nl2e(0).ns2p,nl2e.ns2p,0]
;nl2b = [nl2e(0).ns2p,nl2e.ns2p,nl2e(nintrvl-1).ns2p]
;nl2b = [0,nl2e.ns2p,nl2e(nintrvl-1).ns2p]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
;gradnl2 = (shift(nl2b,-1) - nl2b)/dL
flux = (dll/L^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
;nl2b = nl2b + L^2*(dt_trans/dL)*(flux-shift(flux,1))
;nl2b = smooth(nl2b,2)
;nl2b(whsm) = smooth(nl2b(whsm),2)
nl2e.ns2p = nl2b(1:nintrvl)


;nl2b = [0,nl2.ns3p,0]
nl2b = [nl2(0).ns3p,nl2.ns3p,0]
;nl2b = [nl2(0).ns3p,nl2.ns3p,nl2(nintrvl-1).ns3p]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
;gradnl2 = (shift(nl2b,-1) - nl2b)/dL
flux = (dll/L^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
;nl2b = nl2b + L^2*(dt_trans/dL)*(flux-shift(flux,1))
;nl2b = smooth(nl2b,2)
;nl2b(whsm) = smooth(nl2b(whsm),2)
nl2.ns3p = nl2b(1:nintrvl)

;energy
;solve for boundary temp
ts3pb = tbndry*nl2(nintrvl-1).ns3p*L(nintrvl)^(8./3.)*L(nintrvl)^2
;ts3pb = tbndry*nl2(nintrvl-1).ns3p*L(nintrvl)^(8./3.)

nl2b = [0,nl2e.ns3p,ts3pb]
;nl2b = [nl2e(0).ns3p,nl2e.ns3p,0]
;nl2b = [nl2e(0).ns3p,nl2e.ns3p,nl2e(nintrvl-1).ns3p]
;nl2b = [0,nl2e.ns3p,nl2e(nintrvl-1).ns3p]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
;gradnl2 = (shift(nl2b,-1) - nl2b)/dL
flux = (dll/L^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
;nl2b = nl2b + L^2*(dt_trans/dL)*(flux-shift(flux,1))
;nl2b = smooth(nl2b,2)
;nl2b(whsm) = smooth(nl2b(whsm),2)
nl2e.ns3p = nl2b(1:nintrvl)


;nl2b = [0,nl2.ns4p,0]
nl2b = [nl2(0).ns4p,nl2.ns4p,0]
;nl2b = [nl2(0).ns4p,nl2.ns4p,nl2(nintrvl-1).ns4p]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
;gradnl2 = (shift(nl2b,-1) - nl2b)/dL
flux = (dll/L^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
;nl2b = nl2b + L^2*(dt_trans/dL)*(flux-shift(flux,1))
;nl2b = smooth(nl2b,2)
;nl2b(whsm) = smooth(nl2b(whsm),2)
nl2.ns4p = nl2b(1:nintrvl)

;energy
;solve for boundary temp
topb = 1e3*nl2(nintrvl-1).nop*L(nintrvl)^(8./3.)*L(nintrvl)^2

nl2b = [0,nl2e.ns4p,topb]
;nl2b = [nl2e(0).ns4p,nl2e.ns4p,0]
;nl2b = [nl2e(0).ns4p,nl2e.ns4p,nl2e(nintrvl-1).ns4p]
;nl2b = [0,nl2e.ns4p,nl2e(nintrvl-1).ns4p]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
;gradnl2 = (shift(nl2b,-1) - nl2b)/dL
flux = (dll/L^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
;nl2b = nl2b + L^2*(dt_trans/dL)*(flux-shift(flux,1))
;nl2b = smooth(nl2b,2)
;nl2b(whsm) = smooth(nl2b(whsm),2)
nl2e.ns4p = nl2b(1:nintrvl)



;nl2b = [0,nl2.nop,0]
nl2b = [nl2(0).nop,nl2.nop,0]
;nl2b = [nl2(0).nop,nl2.nop,nl2(nintrvl-1).nop]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
;gradnl2 = (shift(nl2b,-1) - nl2b)/dL
flux = (dll/L^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
;nl2b = nl2b + L^2*(dt_trans/dL)*(flux-shift(flux,1))
;nl2b = smooth(nl2b,2)
;nl2b(whsm) = smooth(nl2b(whsm),2)
nl2.nop = nl2b(1:nintrvl)

;energy
;solve for boundary temp
topb = tbndry*nl2(nintrvl-1).nop*L(nintrvl)^(8./3.)*L(nintrvl)^2
;topb = tbndry*nl2(nintrvl-1).nop*L(nintrvl)^(8./3.)

nl2b = [0,nl2e.nop,topb]
;nl2b = [nl2e(0).nop,nl2e.nop,0]
;nl2b = [nl2e(0).nop,nl2e.nop,nl2e(nintrvl-1).nop]
;nl2b = [0,nl2e.nop,nl2e(nintrvl-1).nop]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
;gradnl2 = (shift(nl2b,-1) - nl2b)/dL
flux = (dll/L^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
;nl2b = nl2b + L^2*(dt_trans/dL)*(flux-shift(flux,1))
;nl2b = smooth(nl2b,2)
;nl2b(whsm) = smooth(nl2b(whsm),2)
nl2e.nop = nl2b(1:nintrvl)

;nl2b = [0,nl2.no2p,0]
nl2b = [nl2(0).no2p,nl2.no2p,0]
;nl2b = [nl2(0).no2p,nl2.no2p,nl2(nintrvl-1).no2p]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
;gradnl2 = (shift(nl2b,-1) - nl2b)/dL
flux = (dll/L^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
;nl2b = nl2b + L^2*(dt_trans/dL)*(flux-shift(flux,1))
;nl2b = smooth(nl2b,2)
;nl2b(whsm) = smooth(nl2b(whsm),2)
nl2.no2p = nl2b(1:nintrvl)

;energy
;solve for boundary temp
to2pb = tbndry*nl2(nintrvl-1).no2p*L(nintrvl)^(8./3.)*L(nintrvl)^2
;to2pb = tbndry*nl2(nintrvl-1).no2p*L(nintrvl)^(8./3.)

nl2b = [0,nl2e.no2p,to2pb]
;nl2b = [nl2e(0).no2p,nl2e.no2p,0]
;nl2b = [nl2e(0).no2p,nl2e.no2p,nl2e(nintrvl-1).no2p]
;nl2b = [0,nl2e.no2p,nl2e(nintrvl-1).no2p]
gradnl2 = (shift(nl2b,-1) - nl2b)/(shift(L,-1)-L)
;gradnl2 = (shift(nl2b,-1) - nl2b)/dL
flux = (dll/L^2)*gradnl2
dL_flx = L_flx - shift(L_flx,1) 
dL_flx(0) = dL_flx(1)
nl2b = nl2b + L^2*(dt_trans/dL_flx)*(flux-shift(flux,1))
;nl2b = nl2b + L^2*(dt_trans/dL)*(flux-shift(flux,1))
;nl2b = smooth(nl2b,2)
;nl2b(whsm) = smooth(nl2b(whsm),2)
nl2e.no2p = nl2b(1:nintrvl)



return
end
;------------------------------------------------------------------


;------------------------------------------------------------------
pro iterate_NL2_to_equator,nl2,nl2e,nr,tr
;------------------------------------------------------------------
@common

nit = 5
ep = 0.00001

nl2_0 = nl2
nl2_p = nl2
nl2e_p = nl2e
nl2_m = nl2
nl2e_m = nl2e
f = nl2
df = nl2

n0 = nr
np = nr
nm = nr

for j = 0,nit do begin

np.nsp = n0.nsp + ep*n0.nsp
np.ns2p = n0.ns2p + ep*n0.ns2p
np.ns3p = n0.ns3p + ep*n0.ns3p
np.ns4p = n0.ns4p + ep*n0.ns4p
np.nop = n0.nop + ep*n0.nop
np.no2p = n0.no2p + ep*n0.no2p

nm.nsp = n0.nsp - ep*n0.nsp
nm.ns2p = n0.ns2p - ep*n0.ns2p
nm.ns3p = n0.ns3p - ep*n0.ns3p
nm.ns4p = n0.ns4p - ep*n0.ns4p
nm.nop = n0.nop - ep*n0.nop
nm.no2p = n0.no2p - ep*n0.no2p

get_NL2,np,tr,nl2_p,nl2e_p
get_NL2,nm,tr,nl2_m,nl2e_m

f.nsp = nl2_p.nsp - nl2_0.nsp
f.ns2p = nl2_p.ns2p - nl2_0.ns2p
f.ns3p = nl2_p.ns3p - nl2_0.ns3p
f.ns4p = nl2_p.ns4p - nl2_0.ns4p
f.nop = nl2_p.nop - nl2_0.nop
f.no2p = nl2_p.no2p - nl2_0.no2p
df.nsp = (nl2_p.nsp - nl2_m.nsp)/(np.nsp - nm.nsp)
df.ns2p = (nl2_p.ns2p - nl2_m.ns2p)/(np.ns2p - nm.ns2p)
df.ns3p = (nl2_p.ns3p - nl2_m.ns3p)/(np.ns3p - nm.ns3p)
df.ns4p = (nl2_p.ns4p - nl2_m.ns4p)/(np.ns4p - nm.ns4p)
df.nop = (nl2_p.nop - nl2_m.nop)/(np.nop - nm.nop)
df.no2p = (nl2_p.no2p - nl2_m.no2p)/(np.no2p - nm.no2p)

whsp = where(abs(f.nsp) gt 0)
whs2p = where(abs(f.ns2p) gt 0)
whs3p = where(abs(f.ns3p) gt 0)
whs4p = where(abs(f.ns4p) gt 0)
whop = where(abs(f.nop) gt 0)
who2p = where(abs(f.no2p) gt 0)

n0(whsp).nsp = n0(whsp).nsp - f(whsp).nsp/df(whsp).nsp
n0(whs2p).ns2p = n0(whs2p).ns2p - f(whs2p).ns2p/df(whs2p).ns2p
n0(whs3p).ns3p = n0(whs3p).ns3p - f(whs3p).ns3p/df(whs3p).ns3p
;n0(whs4p).ns4p = n0(whs4p).ns4p - f(whs4p).ns4p/df(whs4p).ns4p
n0(whop).nop = n0(whop).nop - f(whop).nop/df(whop).nop
n0(who2p).no2p = n0(who2p).no2p - f(who2p).no2p/df(who2p).no2p

endfor

nr = n0

return
end
;------------------------------------------------------------------




