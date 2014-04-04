set_plot,'ps
device,/inches,xsize=6.5,xoff=1.0,ysize = 9.0,yoff = 1.0
device,filename='radial.ps
!p.font=0
!p.charsize=1.2
!x.charsize=1.2
!y.charsize=1.2
!x.thick=2.0
!y.thick=2.0
!p.thick=2.0


nbox = 18
dL0 = 0.25
Lshell = 6.0 + findgen(21)*dL0

restore,'restart.sav'

;window,0
;window,1
;window,2

restore,'janfit.sav
as_r = reverse([6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.15,8.5,8.85])

!p.multi=[0,4,4]
;wset,0

!x.range=[6,9]
!x.title='L (R!dJ!n)'

!y.title='n (cm!u-=3!n)'
plot_io,Lshell,nr.ns,title='S',yrange=[0.1,1000],/ysty
plot_io,Lshell,nr.no,title='O',yrange=[0.1,1000],/ysty
plot_io,Lshell,nr.nsp,title='S+',yrange=[0.01,1000],/ysty
oplot,Lshell,nr.nsph,linestyle=1
plot_io,Lshell,nr.ns2p,title='S++'
plot_io,Lshell,nr.ns3p,title='S+++'
plot_io,Lshell,nr.nop,title='O+',yrange=[0.01,1000],/ysty
oplot,Lshell,nr.noph,linestyle=1
plot_io,Lshell,nr.no2p,title='O++'
plot_io,Lshell,nr.nel,title='Elec'

!y.title='T (eV)'
plot,Lshell,tr.Tsp,title='S+'
plot,Lshell,tr.Ts2p,title='S++'
plot,Lshell,tr.Ts3p,title='S+++'
plot,Lshell,tr.Top,title='O+'
plot,Lshell,tr.To2p,title='O++'
plot,Lshell,tr.Tel,title='Elec'
plot_io,Lshell,tot_peuv*1.6e-19,ytitle='Power (W)',$
     yrange=[1e7,1e12],/ysty,$
     title='Total Peuv...'+string(total(tot_peuv+tot_peuvh)*1.6e-19)+' (W)'
oplot,Lshell,tot_peuvh*1.6e-19,linestyle=1


!p.multi=[0,2,3]

;wset,1
!y.title = 'NL!u2!n'
plot_io,s.l,s.tot33,yrange=[1e34,4e36],/ysty,/xsty,title='Elec'
oplot,Lshell,nl2.nsp+2.0*nl2.ns2p+3.0*nl2.ns3p+nl2.nop+2.0*nl2.no2p,$
   linestyle=1
plot_io,s.l,s.s1,title='S+',yrange=[1e34,4e36],/ysty,/xsty
oplot,Lshell,nl2.nsp,linestyle=1
plot_io,s.l,s.s2,title='S++',yrange=[1e34,4e36],/ysty,/xsty
oplot,Lshell,nl2.ns2p,linestyle=1
plot_io,s.l,s.s3,title='S+++',yrange=[1e34,4e36],/ysty,/xsty
oplot,Lshell,nl2.ns3p,linestyle=1
plot_io,s.l,s.o1,title='O+',yrange=[1e34,4e36],/ysty,/xsty
oplot,Lshell,nl2.nop,linestyle=1
plot_io,s.l,s.o2,title='O++',yrange=[1e34,4e36],/ysty,/xsty
oplot,Lshell,nl2.no2p,linestyle=1

!p.multi=[0,2,3]
;wset,2
!y.title='n/n!de!n'
nelc = nr.nel + 0.1*nr.nel

plot_io,Lshell,mrr.sp,linestyle=5,yrange=[0.01,0.4],/ysty,/xsty,title='S+
oploterr,as_r,s2mix,s2mixerr
plot_io,Lshell,mrr.op,yrange=[0.01,0.4],/ysty,/xsty,$
  linestyle=1,title='O+'
oploterr,as_r,o2mix,o2mixerr
plot_io,Lshell,mrr.s2p,linestyle=2,yrange=[0.01,0.4],/ysty,/xsty,title='S++'
oploterr,as_r,s3mix,s3mixerr

plot_io,Lshell,mrr.o2p,linestyle=4,yrange=[0.01,0.4],/ysty,/xsty,title='O++'
oploterr,as_r,o3mix,o3mixerr

plot_io,Lshell,mrr.s3p,linestyle=3,yrange=[0.01,0.4],/ysty,/xsty,title='S+++'
oploterr,as_r,s4mix,s4mixerr


device,/close

end
;-------------------------------------------------------------------

