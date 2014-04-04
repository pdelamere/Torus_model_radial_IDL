
set_plot,'ps
!p.font=0
!p.thick=2.0
!x.thick=2.0
!y.thick=2.0
device,filename='neutral_cloud.ps

L = findgen(200)/10.+6.0

ns = 10.0/(L-5.0)^3
no = 50.0/(L-5.0)^3

ns(0) = ns(0);+5.0
no(0) = no(0);+25.0

plot_io,L,no,yrange=[0.01,100],/ysty,xrange=[6,10],/xsty,$
    xtitle='L Shell (R!dJ!n)',ytitle='Neutral density (cm!u-3!n)'
oplot,L,ns

xyouts,0.5,0.65,'O',/normal,charsize=1.5
xyouts,0.5,0.35,'S',/normal,charsize=1.5


fh = 0.0025*(L/6.0)^11
wh = where(fh gt 0.2) 
fh(wh) = 0.2

plot_io,L,fh*100.,xrange=[6,10],/xsty,yrange=100*[0.001,0.20],$
    xtitle='L shell (R!dJ!n)',ytitle='f_eh (%)'

dll = 2.0e-7*(L/6.0)^8.0

plot_io,L,(1.0/dll)/60./60./24.,xrange=[6,10],/xsty,ytitle='Transport time (days)',xtitle='L shell (R!dJ!n)'



device,/close

end
