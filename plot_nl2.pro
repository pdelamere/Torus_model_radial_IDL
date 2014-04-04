restore,'restart.sav

set_plot,'ps
!p.font = 0
!x.thick=2.0
!y.thick=2.0
!p.thick=2.0
device,/encapsulated
device,filename='nl2.eps
!p.charsize=1.2
!x.charsize=1.2
!y.charsize=1.2

nl2nel = nl2.nsp + nl2.ns2p + nl2.ns3p + nl2.nop + nl2.no2p

dL = 0.25
L = 6.0+findgen(20)*dL

!x.range=[5,9]
!y.range=[1e34,1e37]
plot_io,L,nl2nel,linestyle=0,xtitle='L (R!dJ!n)',ytitle='NL!u2!n'
oplot,L,nl2.nsp,linestyle=1
oplot,L,nl2.ns2p,linestyle=2
oplot,L,nl2.ns3p,linestyle=3
oplot,L,nl2.nop,linestyle=4
oplot,L,nl2.no2p,linestyle=5

oplot,s.L,s.tot33,thick=4
xyouts,0.2,0.8,'Voyager',/normal,charsize=1.0
xyouts,0.5,0.8,'Cassini model',/normal,charsize=1.0


legend,['Total','S!u+!n','S!u++!n','S!u+++!n','O!u+!n','O!u++!n'],$
    linestyle=[0,1,2,3,4,5],/bottom,/left

device,/close


end
