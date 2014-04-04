
@params

restore,'restart.sav

set_plot,'ps
device,/inches,xsize=6.5,xoff=1.0,ysize = 9.0,yoff = 1.0
device,filename='chemistry.ps
!p.font=0
!p.charsize=1.2
!x.thick=2.0 
!y.thick=2.0
!p.thick=2.0
!p.charsize=1.4
!x.range=[5.8,9.2]
!x.style=1
!y.range=[0.01,100]
!y.style=1

;p.multi=[0,1,1]

a = findgen(17)*(!pi*2/16.)
usersym, cos(a),sin(a)


dL = dL0
Lo = 6.0
L = Lo + findgen(21)*dL

dll= DLL_0*(L/Lo)^DLL_alpha 

nl2elec = nl2.nsp+2.0*nl2.ns2p+3.0*nl2.ns3p+nl2.nop+2.0*nl2.no2p
N = nl2elec/L^2

tau = fltarr(n_elements(L))
dtau = fltarr(n_elements(L))

i = 0
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(i)*dL)
    temp = (t2/t1)/60./60./24.
    tau(0) = temp
    dtau(0) = temp
for i = 1,n_elements(L)-2 do begin
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(i)*dL)
    temp = (t2/t1)/60./60./24.
    tau(i) = tau(i-1) + temp
    dtau(i) = temp
endfor

tau_schreier = [63,84,86.5]

L_schreier = [7.0,7.5,8.0]

!p.multi=[0,2,3]
!y.range=[0.4,200]
!x.range=[5.4,9.2]

plot_io,L,dtau,xtitle='L (R!dJ!n)',ytitle='Time scale (days)',$
   /xsty,/ysty,$
   title='S loss'

oplot,L,(1./rr.is)/8.64e4,psym=1
xyouts,5.7,(1./rr(0).is)/8.64e4,'EI',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.ish)/8.64e4,psym=2
xyouts,5.7,(1./rr(0).ish)/8.64e4+2,'HEI',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k1_s)/8.64e4,psym=4
xyouts,5.7,(1./rr(0).k1_s)/8.64e4,'k1',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k2_s)/8.64e4,psym=8
xyouts,5.7,(1./rr(0).k2_s)/8.64e4,'k2',/data,alignment=0.5,charsize=0.8
;oplot,L,(1./rr.k3_s)/8.64e4,psym=8
oplot,L,(1./rr.k4_s)/8.64e4,psym=5
xyouts,5.7,(1./rr(0).k4_s)/8.64e4+2,'k4',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k9_s)/8.64e4,psym=8
xyouts,5.7,(1./rr(0).k9_s)/8.64e4,'k9',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k10_s)/8.64e4,psym=6
xyouts,5.7,(1./rr(0).k10_s)/8.64e4-1,'k10',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k11_s)/8.64e4,psym=7
xyouts,5.7,(1./rr(0).k11_s)/8.64e4,'k11',/data,alignment=0.5,charsize=0.8
;legend,['Ionization','Hot Electon Ionization','k1 (S + S!u+!n)',$
;    'k4 (S + S!u+++!n)','k10 (S + O!u++!n)',$
;    'k11 (S + O!u++!n)','k2,k9'],psym=[1,2,4,5,6,7,8],/bottom,/right,charsize=0.7


plot_io,L,dtau,xtitle='L (R!dJ!n)',ytitle='Time scale (days)',$
   /xsty,/ysty,$
   title='O loss'
oplot,L,(1./rr.io)/8.64e4,psym=1
xyouts,5.7,(1./rr(0).io)/8.64e4,'EI',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.ioh)/8.64e4,psym=2
xyouts,5.7,(1./rr(0).ioh)/8.64e4,'HEI',/data,alignment=0.5,charsize=0.8

;oplot,L,(1./rr.k5_op)/8.64e4,psym=2
;oplot,L,(1./rr.k6_o)/8.64e4,psym=8
;oplot,L,(1./rr.k7_o2p)/8.64e4,psym=8
;oplot,L,(1./rr.k8_o)/8.64e4,psym=4
oplot,L,(1./rr.k12_o)/8.64e4,psym=4
xyouts,5.7,(1./rr(0).k12_o)/8.64e4,'k12',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k14_o)/8.64e4,psym=5
xyouts,5.7,(1./rr(0).k14_o)/8.64e4,'k14',/data,alignment=0.5,charsize=0.8
;legend,['Ionization','Hot Elec Ionization', 'k12 (O + S!u++!n)','k14 (O + S!u+++!n)'],psym=[1,2,4,5],/bottom,/right,charsize=0.7



;Oplot,L,(1./rr.isp)/8.64e4,psym=3
;oplot,L,(1./rr.is2p)/8.64e4,psym=4

nl2elec = nl2.nsp
N = nl2elec/L^2

tau = fltarr(n_elements(L))
dtau = fltarr(n_elements(L))

i = 0
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(i)*dL)
    temp = (t2/t1)/60./60./24.
    tau(0) = temp
    dtau(0) = temp
for i = 1,n_elements(L)-2 do begin
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(i)*dL)
    temp = (t2/t1)/60./60./24.
    tau(i) = tau(i-1) + temp
    dtau(i) = temp
endfor


plot_io,L,dtau,xtitle='L (R!dJ!n)',ytitle='Time scale (days)',$
   /xsty,/ysty,$
   title='S+ Source'
oplot,L,(1./rr.is)/8.64e4,Psym=1
xyouts,5.7,(1./rr(0).is)/8.64e4,'EI',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k2_s)/8.64e4,Psym=2
xyouts,5.7,(1./rr(0).k2_s)/8.64e4,'k2',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k2_s2p)/8.64e4,Psym=3
xyouts,5.7,(1./rr(0).k2_s)/8.64e4,'k2',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k4_s)/8.64e4,Psym=4
xyouts,5.7,(1./rr(0).k4_s)/8.64e4+1,'k4',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k9_s)/8.64e4,Psym=5
xyouts,5.7,(1./rr(0).k9_s)/8.64e4+1.5,'k9',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k10_s)/8.64e4,Psym=6
xyouts,5.7,(1./rr(0).k10_s)/8.64e4,'k10',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k12_s2p)/8.64e4,Psym=7
xyouts,5.7,(1./rr(0).k12_s2p)/8.64e4-1,'k12',/data,alignment=0.5,charsize=0.8
;legend,['Ionization','k2 (S + S!u++!n)', 'k4 (S + S!u+++!n)', 'k9 (S + O!u+!n)', $
;     'k10 (S + O!u++!n)', 'k12 (O + S!u++!n)'],psym=[1,2,4,5,6,7],/bottom,/right,charsize=0.7


plot_io,L,dtau,xtitle='L (R!dJ!n)',ytitle='Time scale (days)',$
   /xsty,/ysty,$
   title='S+ Loss'
oplot,L,(1./rr.isp)/8.64e4,Psym=1
xyouts,5.7,(1./rr(0).isp)/8.64e4,'EI',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.isph)/8.64e4,Psym=2
xyouts,5.7,(1./rr(0).isph)/8.64e4,'HEI',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k8_sp)/8.64e4,Psym=4
xyouts,5.7,(1./rr(0).k8_sp)/8.64e4,'k8',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k13_sp)/8.64e4,Psym=5
xyouts,5.7,(1./rr(0).k13_sp)/8.64e4,'k13',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k16_sp)/8.64e4,Psym=6
xyouts,5.7,(1./rr(0).k16_sp)/8.64e4,'k16',/data,alignment=0.5,charsize=0.8
;legend,['Ionization','Hot Elec Ionization','k8 (O + S!u+!n)', $
;   'k13 (O!u++!n + S!u+!n)', $
;   'k16 (S!u+++!n + S!u+!n)'],psym=[1,2,4,5,6],/bottom,/right,charsize=0.7


a = findgen(17)*(!pi*2/16.)
usersym, cos(a),sin(a)

nl2elec = nl2.nop
N = nl2elec/L^2

tau = fltarr(n_elements(L))
dtau = fltarr(n_elements(L))

i = 0
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(i)*dL)
    temp = (t2/t1)/60./60./24.
    tau(0) = temp
    dtau(0) = temp
for i = 1,n_elements(L)-2 do begin
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(i)*dL)
    temp = (t2/t1)/60./60./24.
    tau(i) = tau(i-1) + temp
    dtau(i) = temp
endfor

plot_io,L,dtau,xtitle='L (R!dJ!n)',ytitle='Time scale (days)',$
   /xsty,/ysty,$
   title='O+ Source'
oplot,L,(1./rr.io)/8.64e4,Psym=1
xyouts,5.7,2.0,'EI',/data,alignment=0.5,charsize=0.8
;;oplot,L,(1./rr.k6_o)/8.64e4,Psym=2
;;oplot,L,(1./rr.k6_o2p)/8.64e4,Psym=3
;;oplot,L,(1./rr.k8_o)/8.64e4,Psym=4
oplot,L,(1./rr.k10_o2p)/8.64e4,Psym=2
xyouts,5.7,(1./rr(0).k10_o2p)/8.64e4,'k10',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k11_o2p)/8.64e4,Psym=8
xyouts,5.7,(1./rr(0).k11_o2p)/8.64e4,'k11',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k12_o)/8.64e4,Psym=4
xyouts,5.7,(1./rr(0).k12_o)/8.64e4,'k12',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k13_o2p)/8.64e4,Psym=5
xyouts,5.7,(1./rr(0).k13_o2p)/8.64e4,'k13',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k14_o)/8.64e4,Psym=6
xyouts,5.7,(1./rr(0).k14_o)/8.64e4,'k14',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k15_o2p)/8.64e4,Psym=7
xyouts,5.7,(1./rr(0).k15_o2p)/8.64e4,'k15',/data,alignment=0.5,charsize=0.8
;legend,['Ionization','k10 (S + O!u++!n)','k11 (S + O!u++!n)','k12 (O + S!u++!n)', $
;     'k13 (O!u++!n + S!u+!n)', 'k14 (O + S!u+++!n)', 'k15 (O!u++!n + S!u++!n)'], $
;      psym=[1,2,8,4,5,6,7],/bottom,/right,charsize=0.7


plot_io,L,dtau,xtitle='L (R!dJ!n)',ytitle='Time scale (days)',$
   /xsty,/ysty,$
   title='O+ Loss'
oplot,L,(1./rr.ioph)/8.64e4,Psym=1
xyouts,5.7,(1./rr(0).ioph)/8.64e4-75,'HEI',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k9_op)/8.64e4,Psym=2
xyouts,5.7,(1./rr(0).k9_op)/8.64e4,'k9',/data,alignment=0.5,charsize=0.8
;oplot,L,(1./rr.k6_o2p)/8.64e4,Psym=3
;oplot,L,(1./rr.k8_o)/8.64e4,Psym=4
;oplot,L,(1./rr.k10_o2p)/8.64e4,Psym=5
;oplot,L,(1./rr.k11_o2p)/8.64e4,Psym=6
;oplot,L,(1./rr.k12_o)/8.64e4,Psym=7
;oplot,L,(1./rr.k13_o2p)/8.64e4,Psym=7
;oplot,L,(1./rr.k14_o)/8.64e4,Psym=7
;oplot,L,(1./rr.k15_o2p)/8.64e4,Psym=7

;legend,['Hot Elec Ionization','k9 (S + O!u+!n)'],psym=[1,2],$
;   /bottom,/right,charsize=0.7

nl2elec = nl2.no2p
N = nl2elec/L^2

tau = fltarr(n_elements(L))
dtau = fltarr(n_elements(L))

i = 0
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(i)*dL)
    temp = (t2/t1)/60./60./24.
    tau(0) = temp
    dtau(0) = temp
for i = 1,n_elements(L)-2 do begin
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(i)*dL)
    temp = (t2/t1)/60./60./24.
    tau(i) = tau(i-1) + temp
    dtau(i) = temp
endfor

dtau = abs(dtau)

plot_io,L,dtau,xtitle='L (R!dJ!n)',ytitle='Time scale (days)',$
   /xsty,/ysty,$
   title='O++ Source'
oplot,L,(1./rr.ioph)/8.64e4,Psym=1
xyouts,5.7,(1./rr(0).ioph)/8.64e4-75,'HEI',/data,alignment=0.5,charsize=0.8
;legend,['Hot Elec Ionization'], $
;      psym=[1],/bottom,/right,charsize=0.7


plot_io,L,dtau,xtitle='L (R!dJ!n)',ytitle='Time scale (days)',$
   /xsty,/ysty,$
   title='O++ Loss'
;oplot,L,(1./rr.k6_o)/8.64e4,Psym=1
oplot,L,(1./rr.k10_o2p)/8.64e4,Psym=1
xyouts,5.7,(1./rr(0).k10_o2p)/8.64e4,'k10',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k11_o2p)/8.64e4,Psym=2
xyouts,5.7,(1./rr(0).k11_o2p)/8.64e4,'k11',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k13_o2p)/8.64e4,Psym=4
xyouts,5.7,(1./rr(0).k13_o2p)/8.64e4,'k13',/data,alignment=0.5,charsize=0.8
;oplot,L,(1./rr.k15_s2p)/8.64e4,Psym=6

;legend,['k10 (S + O!u++!n)','k11 (S + O!u++!n)','k13 (O!u++!n + S!u+!n)'], $
;      psym=[1,2,4],/bottom,/right,charsize=0.7


a = findgen(17)*(!pi*2/16.)
usersym, cos(a),sin(a)

nl2elec = nl2.ns2p
N = nl2elec/L^2

tau = fltarr(n_elements(L))
dtau = fltarr(n_elements(L))

i = 0
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(i)*dL)
    temp = (t2/t1)/60./60./24.
    tau(0) = temp
    dtau(0) = temp
for i = 1,n_elements(L)-2 do begin
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(i)*dL)
    temp = (t2/t1)/60./60./24.
    tau(i) = tau(i-1) + temp
    dtau(i) = temp
endfor

plot_io,L,dtau,xtitle='L (R!dJ!n)',ytitle='Time scale (days)',$
   /xsty,/ysty,$
   title='S++ Source'
oplot,L,(1./rr.isp)/8.64e4,Psym=1
xyouts,5.7,(1./rr(0).isp)/8.64e4-2,'EI',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.isph)/8.64e4,Psym=2
xyouts,5.7,(1./rr(0).isph)/8.64e4,'HEI',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k4_s3p)/8.64e4,Psym=4
xyouts,5.7,(1./rr(0).k4_s3p)/8.64e4,'k4',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k11_s)/8.64e4,Psym=5
xyouts,5.7,(1./rr(0).k11_s)/8.64e4,'k11',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k13_sp)/8.64e4,Psym=6
xyouts,5.7,(1./rr(0).k13_sp)/8.64e4,'k13',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k14_s3p)/8.64e4,Psym=7
xyouts,5.7,(1./rr(0).k14_s3p)/8.64e4,'k14',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k16_s3p)/8.64e4,Psym=8
xyouts,5.7,(1./rr(0).k16_s3p)/8.64e4,'k16',/data,alignment=0.5,charsize=0.8

;legend,['Ionization','Hot Elec Ionization','k4 (S + S!u+++!n)',$
;         'k11 (S + O!u++!n)','k13 (O!u++!n + S!u+!n)','k14 (O + S!u+++!n)',$
;         'k16 (S!u+++!n + S!u+!n)'], $
;      psym=[1,2,4,5,6,7,8],/bottom,/right,charsize=0.7


plot_io,L,dtau,xtitle='L (R!dJ!n)',ytitle='Time scale (days)',$
   /xsty,/ysty,$
   title='S++ Loss'
oplot,L,(1./rr.is2p)/8.64e4,Psym=1
xyouts,5.7,(1./rr(0).is2p)/8.64e4-50,'EI',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.is2ph)/8.64e4,Psym=2
xyouts,5.7,(1./rr(0).is2ph)/8.64e4-75,'HEI',/data,alignment=0.5,charsize=0.8
;oplot,L,(1./rr.k2_s2p)/8.64e4,Psym=4
;xyouts,5.7,(1./rr(0).k2_s2p)/8.64e4,'k2',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k12_s2p)/8.64e4,Psym=5
xyouts,5.7,(1./rr(0).k12_s2p)/8.64e4,'k12',/data,alignment=0.5,charsize=0.8
;oplot,L,(1./rr.k15_s2p)/8.64e4,Psym=6
;xyouts,5.7,(1./rr(0).k15_s2p)/8.64e4,'k15',/data,alignment=0.5,charsize=0.8
;oplot,L,(1./rr.k14_o)/8.64e4,Psym=7
;oplot,L,(1./rr.k16_s3p)/8.64e4,Psym=8

;legend,['Ionization','Hot Elec Ionization','k2 (S + S!u++!n)',$
;         'k12 (O + S!u++!n)','k15 (O!u++!n + S!u++!n)'], $
;      psym=[1,2,4,5,6],/bottom,/right,charsize=0.7

nl2elec = nl2.ns3p
N = nl2elec/L^2

tau = fltarr(n_elements(L))
dtau = fltarr(n_elements(L))

i = 0
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(i)*dL)
    temp = (t2/t1)/60./60./24.
    tau(0) = temp
    dtau(0) = temp
for i = 1,n_elements(L)-2 do begin
    t1 = (dll(i)/L(i)^2)*(nl2elec(i+1)-nl2elec(i))/dL
    t2 = -total(N(i)*dL)
    temp = (t2/t1)/60./60./24.
    tau(i) = tau(i-1) + temp
    dtau(i) = temp
endfor

dtau = abs(dtau)
;dtau = smooth(dtau,2)

plot_io,L,dtau,xtitle='L (R!dJ!n)',ytitle='Time scale (days)',$
   /xsty,/ysty,$
   title='S+++ Source'
;oplot,L,tau_tot
oplot,L,(1./rr.is2p)/8.64e4,Psym=1
xyouts,5.7,(1./rr(0).is2p)/8.64e4-40,'EI',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.is2ph)/8.64e4,Psym=2
xyouts,5.7,(1./rr(0).is2p)/8.64e4-75,'HEI',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k15_s2p)/8.64e4,Psym=4
;oplot,L,(1./rr.k11_o2p)/8.64e4,Psym=5
;oplot,L,(1./rr.k13_o2p)/8.64e4,Psym=6
;oplot,L,(1./rr.k14_o)/8.64e4,Psym=7
;oplot,L,(1./rr.k16_s3p)/8.64e4,Psym=8

;legend,['Ionization','Hot Elec Ionization','k15 (O!u++!n + S!u++!n)'], $
;      psym=[1,2,4],/bottom,/right,charsize=0.7


plot_io,L,dtau,xtitle='L (R!dJ!n)',ytitle='Time scale (days)',$
   /xsty,/ysty,thick=2.0,$
   title='S+++ Loss'
;oplot,L,(1./rr.is2p)/8.64e4,Psym=1
;oplot,L,(1./rr.is2ph)/8.64e4,Psym=2
oplot,L,(1./rr.k4_s3p)/8.64e4,Psym=1
xyouts,5.7,(1./rr(0).k4_s3p)/8.64e4,'k4',/data,alignment=0.5,charsize=0.8
oplot,L,(1./rr.k14_s3p)/8.64e4,Psym=2
xyouts,5.7,(1./rr(0).k14_s3p)/8.64e4,'k14',/data,alignment=0.5,charsize=0.8
;oplot,L,(1./rr.k16_s3p)/8.64e4,Psym=4
;xyouts,5.7,(1./rr(0).k16_s3p)/8.64e4,'k16',/data,alignment=0.5,charsize=0.8
;oplot,L,(1./rr.k14_o)/8.64e4,Psym=7
;oplot,L,(1./rr.k16_s3p)/8.64e4,Psym=8

;legend,['k4 (S + S!u+++!n)','k14 (O + S!u+++!n)','k16 (S!u+++!n + S!u+!n)'], $
;      psym=[1,2,4],/bottom,/right,charsize=0.7






device,/close

end

