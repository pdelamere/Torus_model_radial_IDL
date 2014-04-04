;---------------------------------------------------------------------
pro iterate_navg,navg,T,narr,n3avg,f,nar,max_theta
; navg contains the flux tube averaged densities
; T contains the temps
; narr is returned by cm3_expand containing 61 elements, n(theta) 
;---------------------------------------------------------------------

theta = (max_theta-findgen(2.0*max_theta + 1.0))*!dtor

volarr = (cos(theta)^7)^1
;volarr = 1.0

n0 = navg
m = navg
n1avg = navg
n2avg = navg
n3avg = navg

n1 = navg
n2 = navg
n0 = navg

dn1avg = navg
dn2avg = navg

dn1 = n1
dn2 = n1

f1 = 2.2
f2 = 2.6

ep = 0.01

n1.nel = navg.nel*f.nel-ep
n1.nop = navg.nop*f.nop-ep
n1.no2p = navg.no2p*f.no2p-ep
n1.nsp = navg.nsp*f.nsp-ep
n1.ns2p = navg.ns2p*f.ns2p-ep
n1.ns3p = navg.ns3p*f.ns3p-ep
n1.nelh = navg.nelh*f.nelh-ep

n2.nel = navg.nel*f.nel+ep
n2.nop = navg.nop*f.nop+ep
n2.no2p = navg.no2p*f.no2p+ep
n2.nsp = navg.nsp*f.nsp+ep
n2.ns2p = navg.ns2p*f.ns2p+ep
n2.ns3p = navg.ns3p*f.ns3p+ep
n2.nelh = navg.nelh*f.nelh+ep

cm3_expand,n1,n1arr,T,nar,max_theta
cm3_expand,n2,n2arr,T,nar,max_theta

wf = n1arr(0,*)*volarr
n1avg.nel = total(n1arr(0,*)*wf)/total(wf)
wf = n1arr(1,*)*volarr
n1avg.nop = total(n1arr(1,*)*wf)/total(wf)
wf = n1arr(2,*)*volarr
n1avg.no2p = total(n1arr(2,*)*wf)/total(wf)
wf = n1arr(3,*)*volarr
n1avg.nsp = total(n1arr(3,*)*wf)/total(wf)
wf = n1arr(4,*)*volarr
n1avg.ns2p = total(n1arr(4,*)*wf)/total(wf)
wf = n1arr(5,*)*volarr
n1avg.ns3p = total(n1arr(5,*)*wf)/total(wf)
wf = n1arr(6,*)*volarr
n1avg.nelh = total(n1arr(6,*)*wf)/total(wf)

wf = n2arr(0,*)*volarr
n2avg.nel = total(n2arr(0,*)*wf)/total(wf)
wf = n2arr(1,*)*volarr
n2avg.nop = total(n2arr(1,*)*wf)/total(wf)
wf = n2arr(2,*)*volarr
n2avg.no2p = total(n2arr(2,*)*wf)/total(wf)
wf = n2arr(3,*)*volarr
n2avg.nsp = total(n2arr(3,*)*wf)/total(wf)
wf = n2arr(4,*)*volarr
n2avg.ns2p = total(n2arr(4,*)*wf)/total(wf)
wf = n2arr(5,*)*volarr
n2avg.ns3p = total(n2arr(5,*)*wf)/total(wf)
wf = n2arr(6,*)*volarr
n2avg.nelh = total(n2arr(6,*)*wf)/total(wf)

m.nel = (n2avg.nel - n1avg.nel)/(n2.nel - n1.nel)
m.nop = (n2avg.nop - n1avg.nop)/(n2.nop - n1.nop)
m.no2p = (n2avg.no2p - n1avg.no2p)/(n2.no2p - n1.no2p)
m.nsp = (n2avg.nsp - n1avg.nsp)/(n2.nsp - n1.nsp)
m.ns2p = (n2avg.ns2p - n1avg.ns2p)/(n2.ns2p - n1.ns2p)
m.ns3p = (n2avg.ns3p - n1avg.ns3p)/(n2.ns3p - n1.ns3p)
m.nelh = (n2avg.nelh - n1avg.nelh)/(n2.nelh - n1.nelh)

n0.nel = -(n2avg.nel - navg.nel)/m.nel + n2.nel
n0.nop = -(n2avg.nop - navg.nop)/m.nel + n2.nop
n0.no2p = -(n2avg.no2p - navg.no2p)/m.nel + n2.no2p
n0.nsp = -(n2avg.nsp - navg.nsp)/m.nel + n2.nsp
n0.ns2p = -(n2avg.ns2p - navg.ns2p)/m.nel + n2.ns2p
n0.ns3p = -(n2avg.ns3p - navg.ns3p)/m.nel + n2.ns3p
n0.nelh = -(n2avg.nelh - navg.nelh)/m.nel + n2.nelh

cm3_expand,n0,narr,T,nar,max_theta

n3avg = n1avg

wf = narr(0,*)*volarr
n3avg.nel = total(narr(0,*)*wf)/total(wf)
wf = narr(1,*)*volarr
n3avg.nop = total(narr(1,*)*wf)/total(wf)
wf = narr(2,*)*volarr
n3avg.no2p = total(narr(2,*)*wf)/total(wf)
wf = narr(3,*)*volarr
n3avg.nsp = total(narr(3,*)*wf)/total(wf)
wf = narr(4,*)*volarr
n3avg.ns2p = total(narr(4,*)*wf)/total(wf)
wf = narr(5,*)*volarr
n3avg.ns3p = total(narr(5,*)*wf)/total(wf)
wf = narr(6,*)*volarr
n3avg.nelh = total(narr(6,*)*wf)/total(wf)

return
end
;---------------------------------------------------------------------



;---------------------------------------------------------------------
pro run_iterate_navg,n,T,nar,max_theta
;---------------------------------------------------------------------

;ribbon
;Te0 = 5.0
;Ti0 = 70.0
;Teh0 = 100.0
;nsp0 = 250.
;ns2p0 = 400.0
;ns3p0 = 50.0
;nop0 = 700.0
;no2p0 = 50.0
;nel0 = nsp0 + 2.0*ns2p0 + 3.0*ns3p0 + nop0 + 2.0*no2p0
;nelh0 = 0.01*nel0

f = {frac, nel: 2.0, nsp: 2.0, ns2p: 2.0, ns3p: 2.0,$
              nop: 2.0, no2p: 2.0, nelh: 2.0}
navg = {density, nel: n.nel, nsp:n.nsp, ns2p: n.ns2p, ns3p: n.ns3p,$
              nop: n.nop, no2p: n.no2p, nelh: n.nelh}
;navg = {density, nel: nel0, nsp:nsp0, ns2p: ns2p0, ns3p: ns3p0,$
;              nop: nop0, no2p: no2p0, nelh: nelh0}
;nar = replicate({density, nel: nel0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
;              nop: nop0, no2p: no2p0, nelh: nelh0},61)

Tavg = {temps, Tel: T.Tel, Tsp: T.Tsp, Ts2p: T.Ts2p, Ts3p: T.Ts3p, $
           Top: T.Top, $
           To2p: T.To2p, Telh: T.Telh} 

;Tavg = {temps, Tel: Te0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, $
;           Top: Ti0, $
;           To2p: Ti0, Telh: Teh0} 


;f.nel = 2.0
;f.nop = 2.0
;f.no2p = 2.0
;f.nsp = 2.0
;f.ns2p = 2.0
;f.ns3p = 2.0

ep = 0.00001

for i = 0,5 do begin

   iterate_navg,navg,Tavg,narr,n1avg,f,nar,max_theta
   f.nel = f.nel*navg.nel/n1avg.nel
   f.nop = f.nop*navg.nop/n1avg.nop
   f.no2p = f.no2p*navg.no2p/n1avg.no2p
   f.nsp = f.nsp*navg.nsp/n1avg.nsp
   f.ns2p = f.ns2p*navg.ns2p/n1avg.ns2p
   f.ns3p = f.ns3p*navg.ns3p/n1avg.ns3p
   f.nelh = f.nelh*navg.nelh/n1avg.nelh

;   print,i,navg.nel/n1avg.nel,navg.nop/n1avg.nop,navg.no2p/n1avg.no2p,$
;         navg.nsp/n1avg.nsp,navg.ns2p/n1avg.ns2p,navg.ns3p/n1avg.ns3p

;   print,'f...',f

;      !p.multi=[0,1,1]
;      restore,'torus_profile.sav'

      plot_oi,nar(*).nel,max_theta-findgen(2.0*max_theta),$
         yrange=[-max_theta,max_theta],/ysty,$
         ;title='System III Long = '+strtrim(string(s3*!radeg),2),$
         ytitle = 'Jovigraphic latitude',xtitle='density (cm!u-3!n)',$
         xrange=[1,10000],/xsty
;      oplot,narr(1,*)+2.0*narr(2,*)+narr(3,*)+2.0*narr(4,*)+3.0*narr(5,*),$
;            max_theta-findgen(2.0*max_theta),linestyle=1,thick=3.0
;      oplot,nar(*).nop,max_theta-findgen(2.0*max_theta),linestyle=1
;      oplot,nar(*).no2p,max_theta-findgen(2.0*max_theta),linestyle=2
      oplot,nar(*).nsp,max_theta-findgen(2.0*max_theta),linestyle=3
;      oplot,nar(*).ns2p,max_theta-findgen(2.0*max_theta),linestyle=4
;      oplot,nar(*).ns3p,max_theta-findgen(2.0*max_theta),linestyle=5
;      oplot,nar(*).nelh,max_theta-findgen(2.0*max_theta),linestyle=6

;      plots,[navg.nel,navg.nel],[!y.crange(0),!y.crange(1)],linestyle=0
;      plots,[navg.nop,navg.nop],[!y.crange(0),!y.crange(1)],linestyle=1
;      plots,[navg.no2p,navg.no2p],[!y.crange(0),!y.crange(1)],linestyle=2
      plots,[navg.nsp,navg.nsp],[!y.crange(0),!y.crange(1)],linestyle=3
;      plots,[navg.ns2p,navg.ns2p],[!y.crange(0),!y.crange(1)],linestyle=4
;      plots,[navg.ns3p,navg.ns3p],[!y.crange(0),!y.crange(1)],linestyle=5
;      plots,[navg.nelh,navg.nelh],[!y.crange(0),!y.crange(1)],linestyle=6

;      plots,[n1avg.nel,n1avg.nel],[!y.crange(0),!y.crange(1)],linestyle=0
;      plots,[n1avg.nop,n1avg.nop],[!y.crange(0),!y.crange(1)],linestyle=1
;      plots,[n1avg.no2p,n1avg.no2p],[!y.crange(0),!y.crange(1)],linestyle=2
      plots,[n1avg.nsp,n1avg.nsp],[!y.crange(0),!y.crange(1)],linestyle=3
;      plots,[n1avg.ns2p,n1avg.ns2p],[!y.crange(0),!y.crange(1)],linestyle=4
;      plots,[n1avg.ns3p,n1avg.ns3p],[!y.crange(0),!y.crange(1)],linestyle=5
;      plots,[n1avg.nelh,n1avg.nelh],[!y.crange(0),!y.crange(1)],linestyle=6


;      legend,['e!u-!n','O!u+!n','O!u++!n','S!u+!n','S!u++!n','S!u+++!n',$
;              'ehot'],$
;         linestyle=[0,1,2,3,4,5,6],/right
;;      plot_oi,Tkappa,max_theta-findgen(2.0*max_theta),$
;;         yrange=[-max_theta,max_theta],/ysty,$
;;         title='Electron Temp (Kappa = '+strtrim(string(kappa),2)+')',$
;;         ytitle = 'Jovigraphic latitude',xtitle='Temperature (eV)',$
;;         xrange=[1,100],/xsty

  endfor

print,navg.nel/max(nar(*).nel),navg.nsp/max(nar(*).nsp),navg.ns2p/max(nar(*).ns2p)

JUMP:

return
end
;---------------------------------------------------------------------


;---------------------------------------------------------------------
pro neutral_navg,n,nar,max_theta
;---------------------------------------------------------------------

;assume values in n are averaged values

theta_0 = 0.2/6.0 ;scale height
theta_off = 0.0*!dtor

;theta = (30 - findgen(61))*!dtor
theta = (max_theta-findgen(2.0*max_theta + 1.0))*!dtor
varr = (cos(theta)^7)^1

;n.ns = 10.0
;n.no = 80.0

;print,n.no,n.ns

;ns_0 = n.ns*total(varr)/total(exp(-theta^2/theta_0^2)*varr)
;no_0 = n.no*total(varr)/total(exp(-theta^2/theta_0^2)*varr)

;ns_0 = n.ns
;no_0 = n.no

;ns_0 = n.ns*total(varr)/total(exp(-theta^2/theta_0^2))
;no_0 = n.no*total(varr)/total(exp(-theta^2/theta_0^2))

;nar(*).ns = ns_0*exp(-theta^2/theta_0^2)
;nar(*).no = no_0*exp(-theta^2/theta_0^2)

nar(*).ns = (n.ns/0.71)*exp(-(theta-theta_off)^2/(4.0*theta_0^2))
nar(*).no = (n.no/0.71)*exp(-(theta-theta_off)^2/(4.0*theta_0^2))

;plot_io,theta*!radeg,nar.no,yrange=[1,2000],/ysty
;oplot,theta*!radeg,nar.ns
;oplot,theta*!radeg,nar.nsp,linestyle=1
;oplot,theta*!radeg,nar.ns2p,linestyle=2
;oplot,theta*!radeg,nar.ns3p,linestyle=3
;oplot,theta*!radeg,nar.nop,linestyle=4
;oplot,theta*!radeg,nar.no2p,linestyle=5

;wait,0.05

;print,n.ns,n.no


;print,n.ns,total(nar(30).ns*exp(-theta^2/theta_0^2)*varr)/total(varr)
;print,n.no,total(nar(30).no*exp(-theta^2/theta_0^2)*varr)/total(varr)

return
end
;---------------------------------------------------------------------


