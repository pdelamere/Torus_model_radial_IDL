set_plot,'ps
device,/inches,xsize=6.5,xoff=1.0,ysize = 9.0,yoff = 1.0
device,filename='radial.eps
device,/color,bits=8
!p.font=0
!p.charsize=1.6
!x.charsize=1.2
!y.charsize=1.2
!x.thick=2.0
!y.thick=2.0
!p.thick=2.0
loadct,26


nbox = 18
dL0 = 0.25
Lshell = 6.0 + findgen(21)*dL0

restore,'restart.sav'

;window,0
;window,1
;window,2

restore,'janfit.sav
as_r = reverse([6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.15,8.5,8.85])

;!p.multi=[0,4,4]
;wset,0

!p.multi=[0,2,3]
!x.title='L (R!dJ!n)'

!x.range=[6,9]
!y.title='n (cm!u-=3!n)'
;plot_io,Lshell,nr.ns,title='S',/xsty
;plot_io,Lshell,nr.no,title='O',/xsty
;plot_io,Lshell,nr.nsp,title='S+',/xsty
;plot_io,Lshell,nr.ns2p,title='S++',/xsty
;plot_io,Lshell,nr.ns3p,title='S+++',/xsty
;plot_io,Lshell,nr.nop,title='O+',/xsty
;plot_io,Lshell,nr.no2p,title='O++',/xsty

net_src = 0.8*13e27*(Lshell/6.0)^(-20)
net_src(0) = net_src(0) + 0.2*13e27
plot_io,Lshell,net_src,title='Neutral Source',ytitle='S!dn!n (s!u-1!n)'
oplot,Lshell,net_src,color=30,thick=4.0
legend,['Total = 13e27 s!u-1!n ~ 700 kg/s'],charsize = 0.9

trans = 2.4e-7*(Lshell/6.0)^7.0
plot_io,Lshell,trans,title='Radial Transport',ytitle='D!dLL!n (s!u-1!n)'
oplot,Lshell,trans,color=30,thick=4.0
trans = 2.4e-7*(Lshell/6.0)^8.0
oplot,Lshell,trans,color=175,thick=4.0
trans = 2.4e-7*(Lshell/6.0)^6.0
oplot,Lshell,trans,color=65,thick=4.0


plot,Lshell,0.23*(Lshell/6.0)^6.0,title='Percent hot electrons - fixed',$
   xtitle='L (R!dJ!n)',ytitle='%'
oplot,Lshell,0.23*(Lshell/6.0)^6.0,color=30,thick=4.0
oplot,Lshell,0.23*(Lshell/6.0)^7.0,color=175,thick=4.0
oplot,Lshell,0.23*(Lshell/6.0)^5.0,color=65,thick=4.0


plot,Lshell,30*(Lshell/6.0)^6.0,title='Hot electron temperature - fixed',$
   xtitle='L (R!dJ!n)',ytitle='T!deh!n (eV)'
oplot,Lshell,30*(Lshell/6.0)^6.0,color=30,thick=4.0

;plot_io,Lshell,nr.nel,/xsty,title='Electron density',xrange=[6,9],$
;  yrange=[100,4000],/ysty
;oplot,Lshell,nr.nel,color=30,thick=4.0


!y.title='T (eV)'
;plot,Lshell,tr.Tsp,title='S+',/xsty
;plot,Lshell,tr.Ts2p,title='S++',/xsty
;plot,Lshell,tr.Ts3p,title='S+++',/xsty
;plot,Lshell,tr.Top,title='O+',/xsty
;plot,Lshell,tr.To2p,title='O++',/xsty

cas_tel = [5.16,4.69,4.46,4.89,5.94,6.80,7.32,7.37,6.69,9.35,15.0,15.6]
cas_tel_err = [0.0889,0.105,0.101,0.149,0.205,0.248,0.385,0.453,0.677,0.680,1.69,2.87]
cas_L = [6.16,6.38,6.60,6.82,7.03,7.25,7.46,7.66,7.85,8.14,8.51,8.84]


;plot,Lshell,tr.Tel,title='Electron temperature',/xsty,xrange=[6,9],$
;   yrange=[0,20],/ysty
;oplot,Lshell,tr.Tel,color=30,thick=4.0
;oplot,cas_L,cas_tel,linestyle=1
;oploterr,cas_L,cas_tel,cas_tel_err
;legend,['Model','Cassini UVIS'],linestyle=[0,1],charsize=0.8,/top,/left,$
;   color=[30,0]
;plot_io,Lshell,smooth(tot_peuv,2)*1.6e-19,ytitle='Power (W)',$
;     xrange=[6,9],/xsty,$
;     title='Total P!deuv!n = '$
;  +strmid(strtrim(string(total((tot_peuv+tot_peuvh)/1e12)*1.6e-19)),5,5)+' (10!u12!n W)',$
;     yrange=[1e9,4e12],/ysty
;oplot,Lshell,smooth(tot_peuv,2)*1.6e-19,color=30,thick=4.0
;oplot,Lshell,smooth(tot_peuvh,2)*1.6e-19,linestyle=1,color=30,thick=4.0
;xyouts,0.15,0.075,'hot',/normal,charsize=1.0
;xyouts,0.22,0.2,'cold',/normal,charsize=1.0

;legend,['Cassini UVIS: 1.5 +/- 0.5 (10!u12!n W)'],/top,/right,charsize=0.75


device,/close
device,filename='radial1.eps


!p.multi=[0,2,3]
!x.range=[6,9]

;wset,1
!y.title = 'NL!u2!n'
plot_io,s.l,s.tot33,yrange=[1e34,4e36],/ysty,/xsty,title='Elec'
oplot,Lshell,nl2.nsp+2.0*nl2.ns2p+3.0*nl2.ns3p+nl2.nop+2.0*nl2.no2p,$
   linestyle=2,color=30,thick=4.0
legend,['Cassini','Voyager'],linestyle=[2,0],color=[10,0],/bottom,/left,charsize=1.2
plot_io,s.l,s.s1,title='S+',yrange=[1e34,4e36],/ysty,/xsty
oplot,Lshell,nl2.nsp,linestyle=2,color=30,thick=4.0
plot_io,s.l,s.s2,title='S++',yrange=[1e34,4e36],/ysty,/xsty
oplot,Lshell,nl2.ns2p,linestyle=2,color=30,thick=4.0
plot_io,s.l,s.s3,title='S+++',yrange=[1e34,4e36],/ysty,/xsty
oplot,Lshell,nl2.ns3p,linestyle=2,color=30,thick=4.0
plot_io,s.l,s.o1,title='O+',yrange=[1e34,4e36],/ysty,/xsty
oplot,Lshell,nl2.nop,linestyle=2,color=30,thick=4.0
plot_io,s.l,s.o2,title='O++',yrange=[1e34,4e36],/ysty,/xsty
oplot,Lshell,nl2.no2p,linestyle=2,color=30,thick=4.0

device,/close
device,filename='radial2.eps

!p.multi=[0,2,3]
;wset,2
!y.title='n/n!de!n'
nelc = nr.nel + 0.1*nr.nel


mrr_1 = mrr
restore,'restart_trans_6.sav'
mrr_2 = mrr
restore,'restart_trans_8.sav'
mrr_3 = mrr
restore,'restart.sav

!x.range=[6,9]
loadct,26
dL = 0.25
Lo = 6.0
L = Lo + findgen(21)*dL

DLL_0 = 2.5e-7
DLL_alpha = 7.0

dll= DLL_0*(L/Lo)^DLL_alpha 

nl2elec = nl2.nsp+2.0*nl2.ns2p+3.0*nl2.ns3p+nl2.nop+2.0*nl2.no2p
N = nl2elec/L^2

tau = fltarr(n_elements(L))

for i = 0,n_elements(L)-2 do begin
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    ;print,(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(0:i)*dL)
    print,t2
    tau(i) = (t2/t1)/60./60./24.
    ;print,tau(i)
endfor

tau_schreier = [63,84,86.5]

L_schreier = [7.0,7.5,8.0]

plot,L,tau,xtitle='L (R!dJ!n)',ytitle='Integrated transport time (days)',xrange=[6,9],/xsty,yrange=[0,140],/ysty
oplot,L_schreier,tau_schreier,psym=1
oplot,L,tau,color=30,thick=4.0
legend,['Schreier et al., [1998]'],psym=[1],/bottom,/right,charsize=0.8


DLL_0 = 2.5e-7
DLL_alpha = 7.0

dll= DLL_0*(L/Lo)^DLL_alpha 

nl2elec = nl2.nsp+2.0*nl2.ns2p+3.0*nl2.ns3p+nl2.nop+2.0*nl2.no2p
N = nl2elec/L^2

tau = fltarr(n_elements(L))

for i = 0,n_elements(L)-2 do begin
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    ;print,(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(0:i)*dL)
    print,t2
    tau(i) = (t2/t1)/60./60./24.
    ;print,tau(i)
endfor

oplot,L,tau,color=65,linestyle=2,thick=4.0

DLL_0 = 2.5e-7
DLL_alpha = 8.0

dll= DLL_0*(L/Lo)^DLL_alpha 

nl2elec = nl2.nsp+2.0*nl2.ns2p+3.0*nl2.ns3p+nl2.nop+2.0*nl2.no2p
N = nl2elec/L^2

tau = fltarr(n_elements(L))

for i = 0,n_elements(L)-2 do begin
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    ;print,(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(0:i)*dL)
    print,t2
    tau(i) = (t2/t1)/60./60./24.
    ;print,tau(i)
endfor

oplot,L,tau,color=175,linestyle=2,thick=4.0


plot_io,Lshell,mrr_1.sp,linestyle=0,yrange=[0.01,0.4],/ysty,/xsty,title='S+'
oplot,Lshell,mrr_1.sp,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.sp,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.sp,linestyle=1,color=175,thick=4.0
oploterr,as_r,s2mix,s2mixerr
plot_io,Lshell,mrr_1.op,yrange=[0.01,0.4],/ysty,/xsty,$
  linestyle=1,title='O+'
oplot,Lshell,mrr_1.op,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.op,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.op,linestyle=1,color=175,thick=4.0
oploterr,as_r,o2mix,o2mixerr
legend,['D!dLL!n (s!u-1!n): 2.4e-7(L/6)^7','         2.4e-7(L/6)^6','         2.4e-7(L/6)^8'],color=[30,65,175],$
        linestyle=[0,1,1],/bottom,/left,charsize=0.7


plot_io,Lshell,mrr_1.s2p,linestyle=2,yrange=[0.01,0.4],/ysty,/xsty,title='S++'
oplot,Lshell,mrr_1.s2p,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.s2p,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.s2p,linestyle=1,color=175,thick=4.0
oploterr,as_r,s3mix,s3mixerr


plot_io,Lshell,mrr_1.o2p,linestyle=4,yrange=[0.01,0.4],/ysty,/xsty,title='O++'
oplot,Lshell,mrr_1.o2p,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.o2p,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.o2p,linestyle=1,color=175,thick=4.0
oploterr,as_r,o3mix,o3mixerr

plot_io,Lshell,mrr_1.s3p,linestyle=3,yrange=[0.01,0.4],/ysty,/xsty,title='S+++'
oplot,Lshell,mrr_1.s3p,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.s3p,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.s3p,linestyle=1,color=175,thick=4.0
oploterr,as_r,s4mix,s4mixerr



device,/close
device,filename='radial3.eps
restore,'restart.sav'
mrr_1 = mrr
tel1 = tr.Tel
restore,'restart_fh_5.sav'
mrr_2 = mrr
tel2 = tr.Tel
restore,'restart_fh_7.sav'
mrr_3 = mrr
tel3 = tr.Tel

!x.range=[6,9]
loadct,26
plot,Lshell,tel1,title='Electron temperature',/xsty,xrange=[6,9],yrange=[0,20],/ysty,ytitle='T!del!n (eV)'
oplot,cas_L,cas_tel,linestyle=1,color=0,thick=4.0
oploterr,cas_L,cas_tel,cas_tel_err
oplot,Lshell,tel1,color=30,linestyle=0,thick=4.0
oplot,Lshell,tel2,color=65,linestyle=1,thick=4.0
oplot,Lshell,tel3,color=175,linestyle=1,thick=4.0


plot_io,Lshell,mrr_1.sp,linestyle=0,yrange=[0.01,0.4],/ysty,/xsty,title='S+'
oplot,Lshell,mrr_1.sp,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.sp,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.sp,linestyle=1,color=175,thick=4.0
oploterr,as_r,s2mix,s2mixerr
plot_io,Lshell,mrr_1.op,yrange=[0.01,0.4],/ysty,/xsty,$
  linestyle=1,title='O+'
oplot,Lshell,mrr_1.op,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.op,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.op,linestyle=1,color=175,thick=4.0
oploterr,as_r,o2mix,o2mixerr
legend,['f!deh!n (%): 0.23(L/6)^6','         0.23(L/6)^5','         0.23(L/6)^7'],color=[30,65,175],$
        linestyle=[0,1,1],/bottom,/left,charsize=0.8


plot_io,Lshell,mrr_1.s2p,linestyle=2,yrange=[0.01,0.4],/ysty,/xsty,title='S++'
oplot,Lshell,mrr_1.s2p,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.s2p,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.s2p,linestyle=1,color=175,thick=4.0
oploterr,as_r,s3mix,s3mixerr


plot_io,Lshell,mrr_1.o2p,linestyle=4,yrange=[0.01,0.4],/ysty,/xsty,title='O++'
oplot,Lshell,mrr_1.o2p,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.o2p,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.o2p,linestyle=1,color=175,thick=4.0
oploterr,as_r,o3mix,o3mixerr

plot_io,Lshell,mrr_1.s3p,linestyle=3,yrange=[0.01,0.4],/ysty,/xsty,title='S+++'
oplot,Lshell,mrr_1.s3p,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.s3p,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.s3p,linestyle=1,color=175,thick=4.0
oploterr,as_r,s4mix,s4mixerr


device,/close

device,filename='radial4.eps
restore,'restart_source_20.sav'
mrr_1 = mrr
tel1 = tr.Tel
restore,'restart_source_10.sav'
mrr_2 = mrr
tel2 = tr.Tel
restore,'restart_source_30.sav'
mrr_3 = mrr
tel3 = tr.Tel

!x.range=[6,9]
loadct,26

sn1 = 13e27*(Lshell/6.0)^(-20)
sn2 = 13e27*(Lshell/6.0)^(-10)
sn3 = 13e27*(Lshell/6.0)^(-30)

plot_io,Lshell,sn1,title='Neutral Source',/xsty,xrange=[6,9],$
    yrange=[1e24,1e29],$
   /ysty,ytitle='S!dn!n (s!u-1!n)'
oplot,Lshell,sn1,linestyle=0,color=30,thick=4.0
oplot,Lshell,sn2,linestyle=1,color=65,thick=4.0
oplot,Lshell,sn3,linestyle=2,color=175,thick=4.0

plot_io,Lshell,mrr_1.sp,linestyle=0,yrange=[0.01,0.4],/ysty,/xsty,title='S+'
oplot,Lshell,mrr_1.sp,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.sp,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.sp,linestyle=2,color=175,thick=4.0
oploterr,as_r,s2mix,s2mixerr
plot_io,Lshell,mrr_1.op,yrange=[0.01,0.4],/ysty,/xsty,$
  linestyle=1,title='O+'
oplot,Lshell,mrr_1.op,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.op,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.op,linestyle=2,color=175,thick=4.0
oploterr,as_r,o2mix,o2mixerr
legend,['S!dn!n (s!u-1!n): 13e27(L/6)^-20','        13e27(L/6)^-10','         13e27(L/6)^-30'],color=[30,65,175],$
        linestyle=[0,1,1],/bottom,/left,charsize=0.8


plot_io,Lshell,mrr_1.s2p,linestyle=2,yrange=[0.01,0.4],/ysty,/xsty,title='S++'
oplot,Lshell,mrr_1.s2p,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.s2p,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.s2p,linestyle=2,color=175,thick=4.0
oploterr,as_r,s3mix,s3mixerr


plot_io,Lshell,mrr_1.o2p,linestyle=4,yrange=[0.01,0.4],/ysty,/xsty,title='O++'
oplot,Lshell,mrr_1.o2p,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.o2p,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.o2p,linestyle=2,color=175,thick=4.0
oploterr,as_r,o3mix,o3mixerr

plot_io,Lshell,mrr_1.s3p,linestyle=3,yrange=[0.01,0.4],/ysty,/xsty,title='S+++'
oplot,Lshell,mrr_1.s3p,linestyle=0,color=30,thick=4.0
oplot,Lshell,mrr_2.s3p,linestyle=1,color=65,thick=4.0
oplot,Lshell,mrr_3.s3p,linestyle=2,color=175,thick=4.0
oploterr,as_r,s4mix,s4mixerr


device,/close


end
;-------------------------------------------------------------------

