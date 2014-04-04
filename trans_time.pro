
@params

restore,'restart.sav

;set_plot,'ps
;device,/inches,xsize=6.5,xoff=1.0,ysize = 9.0,yoff = 1.0
;device,filename='params.ps
;!p.font=0
;!p.charsize=1.2
;!x.thick=2.0 
!p.charsize=1.2
!p.multi=[0,1,1]

dL = dL0
Lo = 6.0
L = Lo + findgen(21)*dL

dll= DLL_0*(L/Lo)^DLL_alpha 

nl2elec = nl2.nsp+2.0*nl2.ns2p+3.0*nl2.ns3p+nl2.nop+2.0*nl2.no2p
nl2elec = nl2.nsp+nl2.ns2p+nl2.ns3p+nl2.nop+nl2.no2p
N = nl2elec/L^2

plot,nl2elec
;stop

tau = fltarr(n_elements(L))

i = 0
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(i)*dL)
    temp = (t2/t1)/60./60./24.
    tau(0) = temp
    print,tau(i),temp


for i = 1,n_elements(L)-2 do begin
;    t1 = (dll(i)/L(i)^2)*(0.5*(nl2elec(i+1)+nl2elec(i))-0.5*(nl2elec(i)+nl2elec(i-1)))/dL
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    ;print,(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(i)*dL)
    ;print,t2
    temp = (t2/t1)/60./60./24.

    tau(i) = tau(i-1) + temp 
    print,tau(i),temp
endfor

tau_schreier = [63,84,86.5]

L_schreier = [7.0,7.5,8.0]

plot,L,tau,xtitle='L (R!dJ!n)',ytitle='Integrated transport time (days)',$
  xrange=[5.75,9],/xsty,yrange=[0,200],/ysty
oplot,L_schreier,tau_schreier,psym=1

;plot_io,L,tau,xtitle='L (R!dJ!n)',ytitle='Transport time (days)',$
;   linestyle=0,yrange=[0.5,100],/ysty,title='Total = '+string(total(tau))+' (days)

;device,/close

end

