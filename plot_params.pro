
@params


set_plot,'ps
device,/inches,xsize=6.5,xoff=1.0,ysize = 9.0,yoff = 1.0
device,filename='params.ps
!p.font=0
!p.charsize=1.2
!x.thick=2.0 

dL = dL0
Lo = 6.0
L = Lo + findgen(21)*dL

!p.multi=[0,2,3]

print,L

dll= DLL_0*(L/Lo)^DLL_alpha 
tau = (1./dll)/60./60./24.

tau_schreier = (1./[8.6e-7,1.1e-6,1.5e-6])/60./60./24.
L_schreier = [7.0,7.5,8.0]

print,DLL_0,DLL_alpha

print,tau

plot_io,L,tau,xtitle='L (R!dJ!n)',ytitle='Transport time (days)',$
   linestyle=0,yrange=[0.5,100],/ysty,title='Total = '+string(total(tau))+' (days)

oplot,L_schreier,tau_schreier,psym=1

;dll= 1.2e-7*(L/Lo)^7.0 
;tau = (1./dll)/60./60./24.
;oplot,L,tau,linestyle=2



net_source_0 = net_source
net_source = 1.0*net_source_0*(L/Lo)^(-20)
net_source(0) = net_source(0) + 0.5*net_source_0

plot_io,L,net_source,xtitle='L (R!dJ!n)',ytitle='Mass loading (#/s)',$
   title= 'Total = ' +string(total(net_source))+' (#/s)

plot_io,L,net_source*20.*1.67e-27,xtitle='L (R!dJ!n)',ytitle='Mass loading (kg/s)',title='Total = ' +string(total(net_source*20.*1.67e-27))+' (kg/s)
       

fh_0 = fh
fh = fh_0*(L/6.0)^fh_alpha
;wh = where(fh ge fh_max)
;fh(wh) = fh_max     

plot,L,fh,xtitle='L (R!dJ!n)',ytitle='f!deh!n'

;Teh0 = 30.0
Teh = Teh0*(L/Lo)^Teh0_alpha
wh = where(Teh ge Teh0_max) 
Teh(wh) = Teh0_max 

plot,L,Teh,xtitle='L (R!dJ!n)',ytitle='T!deh!n (eV)'

device,/close

end

