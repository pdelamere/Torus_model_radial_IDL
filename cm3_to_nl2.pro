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
pro Lax_Wend,n,v,x,dt ;cylindrical coords
;------------------------------------------------------------------

f = n*v

n1 = fltarr(n_elements(n))
f1 = fltarr(n_elements(n))
dx = fltarr(n_elements(n))

for i = 0,n_elements(n)-2 do begin

   dx(i) = x(i+1)-x(i)
   x2 = 0.5*(x(i+1)+x(i))
   n1(i) = 0.5*(n(i) + n(i+1)) - (dt/(2*dx(i)*x2))*(x(i+1)*f(i+1) - x(i)*f(i))
   f1(i) = n1(i)*0.5*(v(i+1)+v(i))

endfor


for i = 1,n_elements(n)-1 do begin
   x2 = 0.5*(x(i)+x(i-1))
   n(i) = n(i) - (dt/(dx(i)*x2))*(x(i)*f1(i) - x(i-1)*f1(i-1))
endfor


return
end
;------------------------------------------------------------------



;------------------------------------------------------------------
pro Lax_Wend_E,n,v,x,dt ;cylindrical coords
;------------------------------------------------------------------

f = n*v
gamma = 5./3.

n1 = fltarr(n_elements(n))
f1 = fltarr(n_elements(n))
dx = fltarr(n_elements(n))

for i = 0,n_elements(n)-2 do begin

   dx(i) = x(i+1)-x(i)
   x2 = 0.5*(x(i+1)+x(i))
   n1(i) = 0.5*(n(i) + n(i+1)) - (dt/(2*dx(i)*x2))*(x(i+1)*f(i+1) - x(i)*f(i)) - $
           (gamma-1)*0.5*(n(i+1)+n(i))*(dt/(2*dx(i)*x2))*(x(i+1)*v(i+1) - x(i)*v(i))
   f1(i) = n1(i)*0.5*(v(i+1)+v(i))
          
endfor


for i = 1,n_elements(n)-1 do begin
   x2 = x(i)
;   xp = 0.5*(x(i)+x(i+1))
   xm = 0.5*(x(i)+x(i-1))
   xp = x(i) + (x(i)-xm)
   dx(i) = xp-xm
   n(i) = n(i) - (dt/(dx(i)*x2))*(xp*f1(i) - xm*f1(i-1)) - $
          (gamma-1)*0.5*n(i)*(dt/(2*dx(i)*x2))*(xp*v(i+1) - xm*v(i))          
endfor


return
end
;------------------------------------------------------------------




;------------------------------------------------------------------
pro transport_flux,nr,tr,DLL_0,DLL_alpha,H,mdot               
;pro transport_NL2,nl2,nl2e
;------------------------------------------------------------------
@common
;@params

xsz = n_elements(nr)

;print,H

tbndry_out = 200.0 ;eV
tbndry_in = 70.0
Lo = 6.0
L = [Lshell(0)-dL(0),Lshell,Lshell(nintrvl-1)+dL(nintrvl-1)]
dL2 = L(1)-L(0)
Lp = L + dL2

dll = fltarr(n_elements(L))
dll = DLL_0*(Lp/Lo)^DLL_alpha

dL = shift(L,-1)-L
dL(0) = dL(1)
;dL(nintrvl-2)=dL(nintrvl-3) 
;dL(nintrvl)=dL(nintrvl-1) 
sz = size(dL)
dL(sz(1)-1) = dL(sz(1)-2)
v = dll*(dL)


r = 2*!pi*Lshell*7.14e7
mp = 1.67e-27
mdot = (sqrt(!pi)/2)*r*1e6*mp*v*7.14e7*1000.*(H.sp*nr.nsp*32 + H.s2p*nr.ns2p*32 + H.s3p*nr.ns3p*32 + $
                                H.op*nr.nop*16 + H.o2p*nr.no2p*16)


;v(0) = v(1) 
;v(n_elements(v)-1) = v(n_elements(v)-2)
;print,v*7.14e7
;stop
v = [v(0),v,v(n_elements(v)-1)]
;v = [0,v,v(n_elements(v)-1)]
;print,dL,v*7.14e7
;stop


;nx = nintrvl-1
;vout = 500./(r(nx)*1e6*mp*1000.*(H.sp*nr(nx).nsp*32 + H.s2p*nr(nx).ns2p*32 + H.s3p*nr(nx).ns3p*32 + $
;                                H.op*nr(nx).nop*16 + H.o2p*nr(nx).no2p*16))

;sz = size(v)
;v(sz(1)-1) = vout
;print,v*7.14e7
;stop
;mdot = [mdot(0),mdot,mdot(nintrvl-1)]

;print,total(dL/v)/60./60./24.

;plot,[6,10],[0,100],/nodata
;for i = 1,n_elements(v)-3 do begin
;   plots,L(i),total(dL/v(0:i))/60./60./24.,/data,psym=1
;   print,L(i),total(dL/v(0:i))/60./60./24.
;endfor
;stop



;L_flx = (shift(L,-1) + L)/2.0
;L_flx(n_elements(L)-1) = L_flx(n_elements(L)-2) + dL(n_elements(dL)-1)
;nL = n_elements(L)
;whsm = where(L ge 9.0)

nsp = [nr(0).nsp,nr.nsp,nr(n_elements(nr)-1).nsp]
;flux = nsp*v
Lax_Wend,nsp,v,L,dt_trans
;nsp = nsp + (dt_trans/(shift(L,-1)-shift(L,1)))*(shift(flux,-1)-shift(flux,1))
nr.nsp = nsp(1:nintrvl)

nsp = [nr(0).nsp,nr.nsp,nr(xsz-1).nsp]
;tsp = [tbndry_in,tr.Tsp,tbndry_out]
tsp = [tr(0).Tsp,tr.Tsp,tr(xsz-1).Tsp]
ntsp = tsp*nsp
Lax_Wend_E,ntsp,v,L,dt_trans
tr.Tsp = ntsp(1:nintrvl)/nr.nsp


ns2p = [nr(0).ns2p,nr.ns2p,nr(n_elements(nr)-1).ns2p]
;flux = ns2p*v
;print,ns2p
;ns2p = ns2p + (dt_trans/(shift(L,-1)-shift(L,1)))*(shift(flux,-1)-shift(flux,1))
Lax_Wend,ns2p,v,L,dt_trans
nr.ns2p = ns2p(1:nintrvl)

ns2p = [nr(0).ns2p,nr.ns2p,nr(xsz-1).ns2p]
ts2p = [tr(0).Ts2p,tr.Ts2p,tr(xsz-1).Ts2p]
nts2p = ts2p*ns2p
Lax_Wend_E,nts2p,v,L,dt_trans
tr.Ts2p = nts2p(1:nintrvl)/nr.ns2p

ns3p = [nr(0).ns3p,nr.ns3p,nr(n_elements(nr)-1).ns3p]
;flux = ns3p*v
;print,ns3p
;ns3p = ns3p + (dt_trans/(shift(L,-1)-shift(L,1)))*(shift(flux,-1)-shift(flux,1))
Lax_Wend,ns3p,v,L,dt_trans
nr.ns3p = ns3p(1:nintrvl)

ns3p = [nr(0).ns3p,nr.ns3p,nr(xsz-1).ns3p]
ts3p = [tr(0).Ts3p,tr.Ts3p,tr(xsz-1).Ts3p]
nts3p = ts3p*ns3p
Lax_Wend_E,nts3p,v,L,dt_trans
tr.Ts3p = nts3p(1:nintrvl)/nr.ns3p

nop = [nr(0).nop,nr.nop,nr(n_elements(nr)-1).nop]
;flux = nop*v
;print,nop
;nop = nop + (dt_trans/(shift(L,-1)-shift(L,1)))*(shift(flux,-1)-shift(flux,1))
Lax_Wend,nop,v,L,dt_trans
nr.nop = nop(1:nintrvl)

nop = [nr(0).nop,nr.nop,nr(xsz-1).nop]
top = [tr(0).Top,tr.Top,tr(xsz-1).Top]
ntop = top*nop
Lax_Wend_E,ntop,v,L,dt_trans
tr.Top = ntop(1:nintrvl)/nr.nop

no2p = [nr(0).no2p,nr.no2p,nr(n_elements(nr)-1).no2p]
;flux = no2p*v
;print,no2p
;no2p = no2p + (dt_trans/(shift(L,-1)-shift(L,1)))*(shift(flux,-1)-shift(flux,1))
Lax_Wend,no2p,v,L,dt_trans
nr.no2p = no2p(1:nintrvl)

no2p = [nr(0).no2p,nr.no2p,nr(xsz-1).no2p]
to2p = [tr(0).To2p,tr.To2p,tr(xsz-1).To2p]
nto2p = to2p*no2p
Lax_Wend_E,nto2p,v,L,dt_trans
tr.To2p = nto2p(1:nintrvl)/nr.no2p






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




