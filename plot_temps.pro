restore,'restart.sav

set_plot,'ps
!p.font = 0
!x.thick=2.0
!y.thick=2.0
!p.thick=2.0
device,/encapsulated
device,filename='iontemp.eps
!p.charsize=1.2
!x.charsize=1.2
!y.charsize=1.2

dL = 0.25
L = 6.0+findgen(20)*dL

!x.range=[6,9]
!y.range=[50,500]

plot_io,L,tr.Tsp,/ysty,linestyle=1,xtitle='L (R!dJ!n)',ytitle='T (eV)'
;oplot,L,tr.Tsp,linestyle=1
oplot,L,tr.Ts2p,linestyle=2
oplot,L,tr.Ts3p,linestyle=3
oplot,L,tr.Top,linestyle=4
oplot,L,tr.To2p,linestyle=5

legend,['S!u+!n','S!u++!n','S!u+++!n','O!u+!n','O!u++!n'],$
    linestyle=[1,2,3,4,5],/bottom,/right

device,/close

end
