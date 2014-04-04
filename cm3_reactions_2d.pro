;-------------------------------------------------------------------
pro sp_production,feh,teh
;-------------------------------------------------------------------
;set_plot,'ps
;@x6x9
;device,/encapsulated
restore,'v1_teh600_feh0023_otos4.0.sav

min_ns = 10e28
max_ns = 50e28
!x.margin=[4,2]
!y.margin=[2,2]
;!y.title='Transport time (days)'
;!x.title='source rate (10!u-4!n cm!u-3!n s!u-1!n)'
!p.charsize=1.0
!p.font=-1
!x.range=[4,20]
!y.range=[15,60]
!y.style = 1
!x.style = 1
Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
xarr = (1e4*xarr/vol)
print,'xarr...',xarr
yarr = (1./yarr)/(3600.*24.)

cx_k0 = 8.1e-9       ;S+ + S++ -> S++ + S+
cx_k1 = 2.4e-8
cx_k2 = 3.0e-10
cx_k3 = 7.8e-9
cx_k4 = 1.32e-8
cx_k5 = 1.32e-8
cx_k6 = 5.2e-10
cx_k7 = 5.4e-9
;r.cx_k8 = k8*9.0e-11   ;Schreier
cx_k8 = 6.0e-11    ;McGrath
cx_k9 = 3.0e-9
cx_k10 = 2.34e-8
cx_k11 = 1.62e-8
cx_k12 = 2.3e-9  ;new
;r.cx_k12 = k12*7.3e-9   ;old
cx_k13 = 1.4e-9   ;Schreier for L=6
cx_k14 = 1.92e-8 ;Schreier
;r.cx_k14 = k14*4.3e-9 ;Barbosa94, Johnson82
;r.cx_k15 = k15*9e-9    ;Schreier for L=6
cx_k15 = 9e-10    ;McGrath for L=6
cx_k16 = 3.6e-10  ;McGrath


nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)

sz = size(narr)
wh = where(nel ge 1.0 )
whx = wh mod sz(1)
why = wh/sz(1)

flr = fltarr(sz(1),sz(2))
flr(whx,why) = 1.0

nel=nel/flr

;device,filename='idl.ps'
;!p.multi=[0,1,1]
;surface,flr
;device,/close

x = findgen(130)+10.
xx = 1e4*x*1e27/2.5e31

;nominal voyager
mdn = (60.0-45.0)/(30.0-50.0)
y0dn = 60.0 - mdn*30.0

;mdn = (25.0-50.0)/(55.0-30.0)
;y0dn = 55.0 - mdn*25.0

ydn = mdn*x + y0dn
whdn = where((ydn ge 45.0) and (ydn le 60))

;device,xsize=6.2,ysize=7.2,/inches,xoffset=1.0,yoffset=3.0
;device,filename='v1_teh600_feh003_otos2.5.eps
;device,filename='v1_test.eps
;device,/color
;!p.multi=[0,4,5]

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)

sz = size(nel)
wh = where(nel le 0.1)
whx = wh mod sz(1)
why = wh/sz(2)
nelec = nel
nelec(whx,why) = 0.000001

r1 = fltarr(sz(1),sz(2)) ;thermal electron impact ionization
r2 = fltarr(sz(1),sz(2)) ;hot electron impact ionization
r3 = fltarr(sz(1),sz(2)) ;k2
r4 = fltarr(sz(1),sz(2)) ;k4
r5 = fltarr(sz(1),sz(2)) ;k9
r6 = fltarr(sz(1),sz(2)) ;k10
r7 = fltarr(sz(1),sz(2)) ;k12
r8 = fltarr(sz(1),sz(2)) ;recombination

r_tot = fltarr(sz(1),sz(2)) ;total reaction rates

for i = 0,sz(1)-1 do begin
   for j = 0,sz(2)-1 do begin
      ris = cfit(16,16,tarr(i,j,0),c)
      rish = cfit(16,16,teh,c)
      r1(i,j) = ris*nel(i,j)*narr(i,j,5)
      r2(i,j) = rish*feh*nel(i,j)*narr(i,j,5)
      r3(i,j) = 2.0*narr(i,j,5)*narr(i,j,1)*cx_k2 
      r4(i,j) = narr(i,j,5)*narr(i,j,2)*cx_k4
      r5(i,j) = narr(i,j,5)*narr(i,j,3)*cx_k9
      r6(i,j) = narr(i,j,5)*narr(i,j,4)*cx_k10
      r7(i,j) = narr(i,j,6)*narr(i,j,1)*cx_k12
      lnt = alog(tarr(i,j,2))
      A = 1.4
      y = -25.2+A*cos((lnt*(1.3*!pi)/(alog(10)) + 0.8*!pi))
      rrs2p = exp(y)
      r8(i,j) = rrs2p*narr(i,j,1)*nel(i,j)
   endfor
endfor

r_tot = r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8

r1 = 100.*r1/r_tot
r2 = 100.*r2/r_tot
r3 = 100.*r3/r_tot
r4 = 100.*r4/r_tot
r5 = 100.*r5/r_tot
r6 = 100.*r6/r_tot
r7 = 100.*r7/r_tot
r8 = 100.*r8/r_tot


;flr(whx,why) = r_tot(whx,why)

!p.multi=[0,2,4]

contour,r1/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,0,1],$
;   levels=[8,10,12,14,16,18,20],$
   title='S + e',c_charsize=1.0,nlev=10
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r2/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S + e(hot)',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r3/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S + S++ -> S+ + S+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r4/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S + S+++ -> S+ + S++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r5/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S + O+ -> S+ + O',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r6/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S + O++ -> S+ + O+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r7/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O + S++ -> O+ + S+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r8/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S++ + e -> S+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1



return
end
;-------------------------------------------------------------------



;-------------------------------------------------------------------
pro sp_loss,feh,teh
;-------------------------------------------------------------------
;set_plot,'ps
;@x6x9
;device,/encapsulated
restore,'v1_teh600_feh0023_otos4.0.sav

min_ns = 10e28
max_ns = 50e28
!x.margin=[4,2]
!y.margin=[2,2]
;!y.title='Transport time (days)'
;!x.title='source rate (10!u-4!n cm!u-3!n s!u-1!n)'
!p.charsize=1.2
!p.font=-1
!x.range=[4,20]
!y.range=[15,60]
!y.style = 1
!x.style = 1
Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
xarr = (1e4*xarr/vol)
print,'xarr...',xarr
yarr = (1./yarr)/(3600.*24.)

cx_k0 = 8.1e-9       ;S+ + S++ -> S++ + S+
cx_k1 = 2.4e-8
cx_k2 = 3.0e-10
cx_k3 = 7.8e-9
cx_k4 = 1.32e-8
cx_k5 = 1.32e-8
cx_k6 = 5.2e-10
cx_k7 = 5.4e-9
;r.cx_k8 = k8*9.0e-11   ;Schreier
cx_k8 = 6.0e-11    ;McGrath
cx_k9 = 3.0e-9
cx_k10 = 2.34e-8
cx_k11 = 1.62e-8
cx_k12 = 2.3e-9  ;new
;r.cx_k12 = k12*7.3e-9   ;old
cx_k13 = 1.4e-9   ;Schreier for L=6
cx_k14 = 1.92e-8 ;Schreier
;r.cx_k14 = k14*4.3e-9 ;Barbosa94, Johnson82
;r.cx_k15 = k15*9e-9    ;Schreier for L=6
cx_k15 = 9e-10    ;McGrath for L=6
cx_k16 = 3.6e-10  ;McGrath

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)

sz = size(narr)
wh = where(nel ge 1.0 )
whx = wh mod sz(1)
why = wh/sz(1)

flr = fltarr(sz(1),sz(2))
flr(whx,why) = 1.0

nel=nel/flr

;device,filename='idl.ps'
;!p.multi=[0,1,1]
;surface,flr
;device,/close

x = findgen(130)+10.
xx = 1e4*x*1e27/2.5e31

;nominal voyager
mdn = (60.0-45.0)/(30.0-50.0)
y0dn = 60.0 - mdn*30.0

;mdn = (25.0-50.0)/(55.0-30.0)
;y0dn = 55.0 - mdn*25.0

ydn = mdn*x + y0dn
whdn = where((ydn ge 45.0) and (ydn le 60))

;device,xsize=6.2,ysize=7.2,/inches,xoffset=1.0,yoffset=3.0
;device,filename='v1_teh600_feh003_otos2.5.eps
;device,filename='v1_test.eps
;device,/color
;!p.multi=[0,4,5]

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)

sz = size(nel)
wh = where(nel le 0.1)
whx = wh mod sz(1)
why = wh/sz(2)
nelec = nel
nelec(whx,why) = 0.000001

r1 = fltarr(sz(1),sz(2)) ;thermal electron impact ionization
r2 = fltarr(sz(1),sz(2)) ;hot electron impact ionization
r3 = fltarr(sz(1),sz(2)) ;k2
r4 = fltarr(sz(1),sz(2)) ;k4
r5 = fltarr(sz(1),sz(2)) ;k9
r6 = fltarr(sz(1),sz(2)) ;k10
r7 = fltarr(sz(1),sz(2)) ;k12

r_tot = fltarr(sz(1),sz(2)) ;total reaction rates

for i = 0,sz(1)-1 do begin
   for j = 0,sz(2)-1 do begin
      risp = cfit(16,15,tarr(i,j,0),c)
      risph = cfit(16,15,teh,c)
      r1(i,j) = risp*nel(i,j)*narr(i,j,0)
      r2(i,j) = risph*feh*nel(i,j)*narr(i,j,0)
      r3(i,j) = narr(i,j,6)*narr(i,j,0)*cx_k8 
      r4(i,j) = narr(i,j,4)*narr(i,j,0)*cx_k13
      r5(i,j) = narr(i,j,2)*narr(i,j,0)*cx_k16
;      r6(i,j) = narr(i,j,5)*narr(i,j,4)*cx_k10
;      r7(i,j) = narr(i,j,6)*narr(i,j,1)*cx_k12
   endfor
endfor

r_tot = r1 + r2 + r3 + r4 + r5

r1 = 100.*r1/r_tot
r2 = 100.*r2/r_tot
r3 = 100.*r3/r_tot
r4 = 100.*r4/r_tot
r5 = 100.*r5/r_tot
;r6 = 100.*r6/r_tot
;r7 = 100.*r7/r_tot


;flr(whx,why) = r_tot(whx,why)

!p.multi=[0,2,4]

contour,r1/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,0,1],$
;   levels=[8,10,12,14,16,18,20],$
   title='S+ + e',c_charsize=1.0,nlev=10
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r2/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S+ + e(hot)',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r3/flr,xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
   level=[2,4,6,8,10,12,14,16],$
   title='O + S+ -> O+ + S',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r4/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O++ + S+ -> O+ + S++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r5/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S+++ + S+ -> S++ + S++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r6/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='S + O++ -> S+ + O+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r7/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O + S++ -> O+ + S+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1



return
end
;-------------------------------------------------------------------



;-------------------------------------------------------------------
pro s2p_production,feh,teh
;-------------------------------------------------------------------
;set_plot,'ps
;@x6x9
;device,/encapsulated
restore,'v1_teh600_feh0023_otos4.0.sav

min_ns = 10e28
max_ns = 50e28
!x.margin=[4,2]
!y.margin=[2,2]
;!y.title='Transport time (days)'
;!x.title='source rate (10!u-4!n cm!u-3!n s!u-1!n)'
!p.charsize=1.2
!p.font=-1
!x.range=[4,20]
!y.range=[15,60]
!y.style = 1
!x.style = 1
Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
xarr = (1e4*xarr/vol)
print,'xarr...',xarr
yarr = (1./yarr)/(3600.*24.)

cx_k0 = 8.1e-9       ;S+ + S++ -> S++ + S+
cx_k1 = 2.4e-8
cx_k2 = 3.0e-10
cx_k3 = 7.8e-9
cx_k4 = 1.32e-8
cx_k5 = 1.32e-8
cx_k6 = 5.2e-10
cx_k7 = 5.4e-9
;r.cx_k8 = k8*9.0e-11   ;Schreier
cx_k8 = 6.0e-11    ;McGrath
cx_k9 = 3.0e-9
cx_k10 = 2.34e-8
cx_k11 = 1.62e-8
cx_k12 = 2.3e-9  ;new
;r.cx_k12 = k12*7.3e-9   ;old
cx_k13 = 1.4e-9   ;Schreier for L=6
cx_k14 = 1.92e-8 ;Schreier
;r.cx_k14 = k14*4.3e-9 ;Barbosa94, Johnson82
;r.cx_k15 = k15*9e-9    ;Schreier for L=6
cx_k15 = 9e-10    ;McGrath for L=6
cx_k16 = 3.6e-10  ;McGrath

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)

sz = size(narr)
wh = where(nel ge 1.0 )
whx = wh mod sz(1)
why = wh/sz(1)

flr = fltarr(sz(1),sz(2))
flr(whx,why) = 1.0

nel=nel/flr

;device,filename='idl.ps'
;!p.multi=[0,1,1]
;surface,flr
;device,/close

x = findgen(130)+10.
xx = 1e4*x*1e27/2.5e31

;nominal voyager
mdn = (60.0-45.0)/(30.0-50.0)
y0dn = 60.0 - mdn*30.0

;mdn = (25.0-50.0)/(55.0-30.0)
;y0dn = 55.0 - mdn*25.0

ydn = mdn*x + y0dn
whdn = where((ydn ge 45.0) and (ydn le 60))

;device,xsize=6.2,ysize=7.2,/inches,xoffset=1.0,yoffset=3.0
;device,filename='v1_teh600_feh003_otos2.5.eps
;device,filename='v1_test.eps
;device,/color
;!p.multi=[0,4,5]

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)

sz = size(nel)
wh = where(nel le 0.1)
whx = wh mod sz(1)
why = wh/sz(2)
nelec = nel
nelec(whx,why) = 0.000001

r1 = fltarr(sz(1),sz(2)) ;thermal electron impact ionization
r2 = fltarr(sz(1),sz(2)) ;hot electron impact ionization
r3 = fltarr(sz(1),sz(2)) ;k2
r4 = fltarr(sz(1),sz(2)) ;k4
r5 = fltarr(sz(1),sz(2)) ;k9
r6 = fltarr(sz(1),sz(2)) ;k10
r7 = fltarr(sz(1),sz(2)) ;k12
r8 = fltarr(sz(1),sz(2)) ;k12

r_tot = fltarr(sz(1),sz(2)) ;total reaction rates

for i = 0,sz(1)-1 do begin
   for j = 0,sz(2)-1 do begin
      risp = cfit(16,15,tarr(i,j,0),c)
      risph = cfit(16,15,teh,c)
      r1(i,j) = risp*nel(i,j)*narr(i,j,0)
      r2(i,j) = risph*feh*nel(i,j)*narr(i,j,0)
      lnt = alog(tarr(i,j,3))
      A = 1.0
      y = -23.9+A*cos((lnt*(1.3*!pi)/(alog(10)) + 0.74*!pi))
      rrs3p = exp(y)
      r3(i,j) = rrs3p*narr(i,j,2)*nel(i,j) 
      r4(i,j) = narr(i,j,5)*narr(i,j,2)*cx_k4
      r5(i,j) = narr(i,j,5)*narr(i,j,4)*cx_k11
      r6(i,j) = narr(i,j,4)*narr(i,j,0)*cx_k13
      r7(i,j) = narr(i,j,6)*narr(i,j,2)*cx_k14
      r8(i,j) = 2.0*narr(i,j,2)*narr(i,j,0)*cx_k16
   endfor
endfor

r_tot = r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8

r1 = 100.*r1/r_tot
r2 = 100.*r2/r_tot
r3 = 100.*r3/r_tot
r4 = 100.*r4/r_tot
r5 = 100.*r5/r_tot
r6 = 100.*r6/r_tot
r7 = 100.*r7/r_tot
r8 = 100.*r8/r_tot


;flr(whx,why) = r_tot(whx,why)

!p.multi=[0,2,4]

contour,r1/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,0,1],$
;   levels=[8,10,12,14,16,18,20],$
   title='S+ + e',c_charsize=1.0,nlev=10
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r2/flr,xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
   levels=[10,12,14,16,18,20,22,24,26,28,30],$
   title='S+ + e(hot)',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r3/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S+++ + e -> S++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r4/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S + S+++ -> S+ + S++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r5/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S + O++ -> S++ + O+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r6/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O++ + S+ -> O+ + S++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r7/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O + S+++ -> O+ + S++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r8/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S+++ + S+ -> S++ + S++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1



return
end
;-------------------------------------------------------------------



;-------------------------------------------------------------------
pro s2p_loss,feh,teh
;-------------------------------------------------------------------
;set_plot,'ps
;@x6x9
;device,/encapsulated
restore,'v1_teh600_feh0023_otos4.0.sav

min_ns = 10e28
max_ns = 50e28
!x.margin=[4,2]
!y.margin=[2,2]
;!y.title='Transport time (days)'
;!x.title='source rate (10!u-4!n cm!u-3!n s!u-1!n)'
!p.charsize=1.2
!p.font=-1
!x.range=[4,20]
!y.range=[15,60]
!y.style = 1
!x.style = 1
Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
xarr = (1e4*xarr/vol)
print,'xarr...',xarr
yarr = (1./yarr)/(3600.*24.)

cx_k0 = 8.1e-9       ;S+ + S++ -> S++ + S+
cx_k1 = 2.4e-8
cx_k2 = 3.0e-10
cx_k3 = 7.8e-9
cx_k4 = 1.32e-8
cx_k5 = 1.32e-8
cx_k6 = 5.2e-10
cx_k7 = 5.4e-9
;r.cx_k8 = k8*9.0e-11   ;Schreier
cx_k8 = 6.0e-11    ;McGrath
cx_k9 = 3.0e-9
cx_k10 = 2.34e-8
cx_k11 = 1.62e-8
cx_k12 = 2.3e-9  ;new
;r.cx_k12 = k12*7.3e-9   ;old
cx_k13 = 1.4e-9   ;Schreier for L=6
cx_k14 = 1.92e-8 ;Schreier
;r.cx_k14 = k14*4.3e-9 ;Barbosa94, Johnson82
;r.cx_k15 = k15*9e-9    ;Schreier for L=6
cx_k15 = 9e-10    ;McGrath for L=6
cx_k16 = 3.6e-10  ;McGrath

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)

sz = size(narr)
wh = where(nel ge 1.0 )
whx = wh mod sz(1)
why = wh/sz(1)

flr = fltarr(sz(1),sz(2))
flr(whx,why) = 1.0

nel=nel/flr

;device,filename='idl.ps'
;!p.multi=[0,1,1]
;surface,flr
;device,/close

x = findgen(130)+10.
xx = 1e4*x*1e27/2.5e31

;nominal voyager
mdn = (60.0-45.0)/(30.0-50.0)
y0dn = 60.0 - mdn*30.0

;mdn = (25.0-50.0)/(55.0-30.0)
;y0dn = 55.0 - mdn*25.0

ydn = mdn*x + y0dn
whdn = where((ydn ge 45.0) and (ydn le 60))

;device,xsize=6.2,ysize=7.2,/inches,xoffset=1.0,yoffset=3.0
;device,filename='v1_teh600_feh003_otos2.5.eps
;device,filename='v1_test.eps
;device,/color
;!p.multi=[0,4,5]

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)

sz = size(nel)
wh = where(nel le 0.1)
whx = wh mod sz(1)
why = wh/sz(2)
nelec = nel
nelec(whx,why) = 0.000001

r1 = fltarr(sz(1),sz(2)) ;thermal electron impact ionization
r2 = fltarr(sz(1),sz(2)) ;hot electron impact ionization
r3 = fltarr(sz(1),sz(2)) ;k2
r4 = fltarr(sz(1),sz(2)) ;k4
r5 = fltarr(sz(1),sz(2)) ;k9
r6 = fltarr(sz(1),sz(2)) ;k10
r7 = fltarr(sz(1),sz(2)) ;k12
r8 = fltarr(sz(1),sz(2)) ;k12

r_tot = fltarr(sz(1),sz(2)) ;total reaction rates

for i = 0,sz(1)-1 do begin
   for j = 0,sz(2)-1 do begin
      ris2p = cfit(16,14,tarr(i,j,0),c)
      ris2ph = cfit(16,14,teh,c)
      r1(i,j) = ris2p*nel(i,j)*narr(i,j,1)
      r2(i,j) = ris2ph*feh*nel(i,j)*narr(i,j,1)
      lnt = alog(tarr(i,j,2))
      A = 1.4
      y = -25.2+A*cos((lnt*(1.3*!pi)/(alog(10)) + 0.8*!pi))
      rrs2p = exp(y)
      r3(i,j) = rrs2p*narr(i,j,1)*nel(i,j) 
      r4(i,j) = narr(i,j,5)*narr(i,j,1)*cx_k2
      r5(i,j) = narr(i,j,6)*narr(i,j,1)*cx_k12
      r6(i,j) = narr(i,j,4)*narr(i,j,1)*cx_k15
;      r7(i,j) = narr(i,j,6)*narr(i,j,2)*cx_k14
;      r8(i,j) = 2.0*narr(i,j,2)*narr(i,j,0)*cx_k16
   endfor
endfor

r_tot = r1 + r2 + r3 + r4 + r5 + r6

r1 = 100.*r1/r_tot
r2 = 100.*r2/r_tot
r3 = 100.*r3/r_tot
r4 = 100.*r4/r_tot
r5 = 100.*r5/r_tot
r6 = 100.*r6/r_tot
r7 = 100.*r7/r_tot
r8 = 100.*r8/r_tot


;flr(whx,why) = r_tot(whx,why)

!p.multi=[0,2,4]

contour,r1/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,0,1],$
;   levels=[8,10,12,14,16,18,20],$
   title='S++ + e',c_charsize=1.0,nlev=10
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r2/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S++ + e(hot)',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r3/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S++ + e -> S+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r4/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S + S++ -> S+ + S+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r5/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O + S++ -> O+ + S+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r6/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O++ + S++ -> O+ + S+++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r7/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O + S+++ -> O+ + S++',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r8/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='S+++ + S+ -> S++ + S++',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1



return
end
;-------------------------------------------------------------------



;-------------------------------------------------------------------
pro s3p_production,feh,teh
;-------------------------------------------------------------------
;set_plot,'ps
;@x6x9
;device,/encapsulated
restore,'v1_teh600_feh0023_otos4.0.sav

min_ns = 10e28
max_ns = 50e28
!x.margin=[4,2]
!y.margin=[2,2]
;!y.title='Transport time (days)'
;!x.title='source rate (10!u-4!n cm!u-3!n s!u-1!n)'
!p.charsize=1.2
!p.font=-1
!x.range=[4,20]
!y.range=[15,60]
!y.style = 1
!x.style = 1
Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
xarr = (1e4*xarr/vol)
print,'xarr...',xarr
yarr = (1./yarr)/(3600.*24.)

cx_k0 = 8.1e-9       ;S+ + S++ -> S++ + S+
cx_k1 = 2.4e-8
cx_k2 = 3.0e-10
cx_k3 = 7.8e-9
cx_k4 = 1.32e-8
cx_k5 = 1.32e-8
cx_k6 = 5.2e-10
cx_k7 = 5.4e-9
;r.cx_k8 = k8*9.0e-11   ;Schreier
cx_k8 = 6.0e-11    ;McGrath
cx_k9 = 3.0e-9
cx_k10 = 2.34e-8
cx_k11 = 1.62e-8
cx_k12 = 2.3e-9  ;new
;r.cx_k12 = k12*7.3e-9   ;old
cx_k13 = 1.4e-9   ;Schreier for L=6
cx_k14 = 1.92e-8 ;Schreier
;r.cx_k14 = k14*4.3e-9 ;Barbosa94, Johnson82
;r.cx_k15 = k15*9e-9    ;Schreier for L=6
cx_k15 = 9e-10    ;McGrath for L=6
cx_k16 = 3.6e-10  ;McGrath

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)

sz = size(narr)
wh = where(nel ge 1.0 )
whx = wh mod sz(1)
why = wh/sz(1)

flr = fltarr(sz(1),sz(2))
flr(whx,why) = 1.0

nel=nel/flr

;device,filename='idl.ps'
;!p.multi=[0,1,1]
;surface,flr
;device,/close

x = findgen(130)+10.
xx = 1e4*x*1e27/2.5e31

;nominal voyager
mdn = (60.0-45.0)/(30.0-50.0)
y0dn = 60.0 - mdn*30.0

;mdn = (25.0-50.0)/(55.0-30.0)
;y0dn = 55.0 - mdn*25.0

ydn = mdn*x + y0dn
whdn = where((ydn ge 45.0) and (ydn le 60))

;device,xsize=6.2,ysize=7.2,/inches,xoffset=1.0,yoffset=3.0
;device,filename='v1_teh600_feh003_otos2.5.eps
;device,filename='v1_test.eps
;device,/color
;!p.multi=[0,4,5]

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)

sz = size(nel)
wh = where(nel le 0.1)
whx = wh mod sz(1)
why = wh/sz(2)
nelec = nel
nelec(whx,why) = 0.000001

r1 = fltarr(sz(1),sz(2)) ;thermal electron impact ionization
r2 = fltarr(sz(1),sz(2)) ;hot electron impact ionization
r3 = fltarr(sz(1),sz(2)) ;k2
r4 = fltarr(sz(1),sz(2)) ;k4
r5 = fltarr(sz(1),sz(2)) ;k9
r6 = fltarr(sz(1),sz(2)) ;k10
r7 = fltarr(sz(1),sz(2)) ;k12
r8 = fltarr(sz(1),sz(2)) ;k12

r_tot = fltarr(sz(1),sz(2)) ;total reaction rates

for i = 0,sz(1)-1 do begin
   for j = 0,sz(2)-1 do begin
      ris2p = cfit(16,14,tarr(i,j,0),c)
      ris2ph = cfit(16,14,teh,c)
      r1(i,j) = ris2p*nel(i,j)*narr(i,j,1)
      r2(i,j) = ris2ph*feh*nel(i,j)*narr(i,j,1)
      r3(i,j) = narr(i,j,4)*narr(i,j,1)*cx_k15 
;      r4(i,j) = narr(i,j,5)*narr(i,j,1)*cx_k2
;      r5(i,j) = narr(i,j,6)*narr(i,j,1)*cx_k12
;      r6(i,j) = narr(i,j,4)*narr(i,j,1)*cx_k15
;      r7(i,j) = narr(i,j,6)*narr(i,j,2)*cx_k14
;      r8(i,j) = 2.0*narr(i,j,2)*narr(i,j,0)*cx_k16
   endfor
endfor

r_tot = r1 + r2 + r3 

r1 = 100.*r1/r_tot
r2 = 100.*r2/r_tot
r3 = 100.*r3/r_tot
r4 = 100.*r4/r_tot
r5 = 100.*r5/r_tot
r6 = 100.*r6/r_tot
r7 = 100.*r7/r_tot
r8 = 100.*r8/r_tot


;flr(whx,why) = r_tot(whx,why)

!p.multi=[0,2,4]

contour,r1/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,0,1],$
;   levels=[8,10,12,14,16,18,20],$
   title='S++ + e',c_charsize=1.0,nlev=10
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r2/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S++ + e(hot)',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r3/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O++ + S++ -> O+ + S+++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r4/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='S + S++ -> S+ + S+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r5/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O + S++ -> O+ + S+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r6/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O++ + S++ -> O+ + S+++',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r7/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O + S+++ -> O+ + S++',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r8/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='S+++ + S+ -> S++ + S++',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1



return
end
;-------------------------------------------------------------------



;-------------------------------------------------------------------
pro s3p_loss,feh,teh
;-------------------------------------------------------------------
;set_plot,'ps
;@x6x9
;device,/encapsulated
restore,'v1_teh600_feh0023_otos4.0.sav

min_ns = 10e28
max_ns = 50e28
!x.margin=[4,2]
!y.margin=[2,2]
;!y.title='Transport time (days)'
;!x.title='source rate (10!u-4!n cm!u-3!n s!u-1!n)'
!p.charsize=1.2
!p.font=-1
!x.range=[4,20]
!y.range=[15,60]
!y.style = 1
!x.style = 1
Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
xarr = (1e4*xarr/vol)
print,'xarr...',xarr
yarr = (1./yarr)/(3600.*24.)

cx_k0 = 8.1e-9       ;S+ + S++ -> S++ + S+
cx_k1 = 2.4e-8
cx_k2 = 3.0e-10
cx_k3 = 7.8e-9
cx_k4 = 1.32e-8
cx_k5 = 1.32e-8
cx_k6 = 5.2e-10
cx_k7 = 5.4e-9
;r.cx_k8 = k8*9.0e-11   ;Schreier
cx_k8 = 6.0e-11    ;McGrath
cx_k9 = 3.0e-9
cx_k10 = 2.34e-8
cx_k11 = 1.62e-8
cx_k12 = 2.3e-9  ;new
;r.cx_k12 = k12*7.3e-9   ;old
cx_k13 = 1.4e-9   ;Schreier for L=6
cx_k14 = 1.92e-8 ;Schreier
;r.cx_k14 = k14*4.3e-9 ;Barbosa94, Johnson82
;r.cx_k15 = k15*9e-9    ;Schreier for L=6
cx_k15 = 9e-10    ;McGrath for L=6
cx_k16 = 3.6e-10  ;McGrath

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)

sz = size(narr)
wh = where(nel ge 1.0 )
whx = wh mod sz(1)
why = wh/sz(1)

flr = fltarr(sz(1),sz(2))
flr(whx,why) = 1.0

nel=nel/flr

;device,filename='idl.ps'
;!p.multi=[0,1,1]
;surface,flr
;device,/close

x = findgen(130)+10.
xx = 1e4*x*1e27/2.5e31

;nominal voyager
mdn = (60.0-45.0)/(30.0-50.0)
y0dn = 60.0 - mdn*30.0

;mdn = (25.0-50.0)/(55.0-30.0)
;y0dn = 55.0 - mdn*25.0

ydn = mdn*x + y0dn
whdn = where((ydn ge 45.0) and (ydn le 60))

;device,xsize=6.2,ysize=7.2,/inches,xoffset=1.0,yoffset=3.0
;device,filename='v1_teh600_feh003_otos2.5.eps
;device,filename='v1_test.eps
;device,/color
;!p.multi=[0,4,5]

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)

sz = size(nel)
wh = where(nel le 0.1)
whx = wh mod sz(1)
why = wh/sz(2)
nelec = nel
nelec(whx,why) = 0.000001

r1 = fltarr(sz(1),sz(2)) ;thermal electron impact ionization
r2 = fltarr(sz(1),sz(2)) ;hot electron impact ionization
r3 = fltarr(sz(1),sz(2)) ;k2
r4 = fltarr(sz(1),sz(2)) ;k4
r5 = fltarr(sz(1),sz(2)) ;k9
r6 = fltarr(sz(1),sz(2)) ;k10
r7 = fltarr(sz(1),sz(2)) ;k12
r8 = fltarr(sz(1),sz(2)) ;k12

r_tot = fltarr(sz(1),sz(2)) ;total reaction rates

for i = 0,sz(1)-1 do begin
   for j = 0,sz(2)-1 do begin
;      ris2p = cfit(16,14,tarr(i,j,0),c)
;      ris2ph = cfit(16,14,teh,c)
      lnt = alog(tarr(i,j,3))
      A = 1.0
      y = -23.9+A*cos((lnt*(1.3*!pi)/(alog(10)) + 0.74*!pi))
      rrs3p = exp(y)
      r1(i,j) = rrs3p*narr(i,j,2)*nel(i,j) 
      r2(i,j) = narr(i,j,5)*narr(i,j,2)*cx_k4
      r3(i,j) = narr(i,j,6)*narr(i,j,2)*cx_k14
      r4(i,j) = narr(i,j,0)*narr(i,j,2)*cx_k16
;      r5(i,j) = narr(i,j,6)*narr(i,j,1)*cx_k12
;      r6(i,j) = narr(i,j,4)*narr(i,j,1)*cx_k15
;      r7(i,j) = narr(i,j,6)*narr(i,j,2)*cx_k14
;      r8(i,j) = 2.0*narr(i,j,2)*narr(i,j,0)*cx_k16
   endfor
endfor

r_tot = r1 + r2 + r3 + r4 

r1 = 100.*r1/r_tot
r2 = 100.*r2/r_tot
r3 = 100.*r3/r_tot
r4 = 100.*r4/r_tot
r5 = 100.*r5/r_tot
r6 = 100.*r6/r_tot
r7 = 100.*r7/r_tot
r8 = 100.*r8/r_tot


;flr(whx,why) = r_tot(whx,why)

!p.multi=[0,2,4]

contour,r1/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,0,1],$
;   levels=[8,10,12,14,16,18,20],$
   title='S+++ + e - > S++',c_charsize=1.0,nlev=10
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r2/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S + S+++ -> S+ + S++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r3/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O + S+++ -> O+ + S++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r4/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S+++ + S+ -> S++ + S++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r5/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O + S++ -> O+ + S+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r6/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O++ + S++ -> O+ + S+++',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r7/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O + S+++ -> O+ + S++',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r8/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='S+++ + S+ -> S++ + S++',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1



return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro op_production,feh,teh
;-------------------------------------------------------------------
;set_plot,'ps
;@x6x9
;device,/encapsulated
restore,'v1_teh600_feh0023_otos4.0.sav

min_ns = 10e28
max_ns = 50e28
!x.margin=[4,2]
!y.margin=[2,2]
;!y.title='Transport time (days)'
;!x.title='source rate (10!u-4!n cm!u-3!n s!u-1!n)'
!p.charsize=1.0
!p.font=-1
!x.range=[4,20]
!y.range=[15,60]
!y.style = 1
!x.style = 1
Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
xarr = (1e4*xarr/vol)
print,'xarr...',xarr
yarr = (1./yarr)/(3600.*24.)

cx_k0 = 8.1e-9       ;S+ + S++ -> S++ + S+
cx_k1 = 2.4e-8
cx_k2 = 3.0e-10
cx_k3 = 7.8e-9
cx_k4 = 1.32e-8
cx_k5 = 1.32e-8
cx_k6 = 5.2e-10
cx_k7 = 5.4e-9
;r.cx_k8 = k8*9.0e-11   ;Schreier
cx_k8 = 6.0e-11    ;McGrath
cx_k9 = 3.0e-9
cx_k10 = 2.34e-8
cx_k11 = 1.62e-8
cx_k12 = 2.3e-9  ;new
;r.cx_k12 = k12*7.3e-9   ;old
cx_k13 = 1.4e-9   ;Schreier for L=6
cx_k14 = 1.92e-8 ;Schreier
;r.cx_k14 = k14*4.3e-9 ;Barbosa94, Johnson82
;r.cx_k15 = k15*9e-9    ;Schreier for L=6
cx_k15 = 9e-10    ;McGrath for L=6
cx_k16 = 3.6e-10  ;McGrath


nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)

sz = size(narr)
wh = where(nel ge 1.0 )
whx = wh mod sz(1)
why = wh/sz(1)

flr = fltarr(sz(1),sz(2))
flr(whx,why) = 1.0

nel=nel/flr

;device,filename='idl.ps'
;!p.multi=[0,1,1]
;surface,flr
;device,/close

x = findgen(130)+10.
xx = 1e4*x*1e27/2.5e31

;nominal voyager
mdn = (60.0-45.0)/(30.0-50.0)
y0dn = 60.0 - mdn*30.0

;mdn = (25.0-50.0)/(55.0-30.0)
;y0dn = 55.0 - mdn*25.0

ydn = mdn*x + y0dn
whdn = where((ydn ge 45.0) and (ydn le 60))

;device,xsize=6.2,ysize=7.2,/inches,xoffset=1.0,yoffset=3.0
;device,filename='v1_teh600_feh003_otos2.5.eps
;device,filename='v1_test.eps
;device,/color
;!p.multi=[0,4,5]

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)

sz = size(nel)
wh = where(nel le 0.1)
whx = wh mod sz(1)
why = wh/sz(2)
nelec = nel
nelec(whx,why) = 0.000001

r1 = fltarr(sz(1),sz(2)) ;thermal electron impact ionization
r2 = fltarr(sz(1),sz(2)) ;hot electron impact ionization
r3 = fltarr(sz(1),sz(2)) ;k2
r4 = fltarr(sz(1),sz(2)) ;k4
r5 = fltarr(sz(1),sz(2)) ;k9
r6 = fltarr(sz(1),sz(2)) ;k10
r7 = fltarr(sz(1),sz(2)) ;k12
r8 = fltarr(sz(1),sz(2)) ;recombination
r9 = fltarr(sz(1),sz(2)) ;recombination
r10 = fltarr(sz(1),sz(2)) ;recombination
r11 = fltarr(sz(1),sz(2)) ;recombination

r_tot = fltarr(sz(1),sz(2)) ;total reaction rates

for i = 0,sz(1)-1 do begin
   for j = 0,sz(2)-1 do begin
      rio = 1.015e-8*sqrt(tarr(i,j,0))*exp(-14.54/tarr(i,j,0))
      rioh = 1.015e-8*sqrt(teh)*exp(-14.54/teh)
      r1(i,j) = rio*nel(i,j)*narr(i,j,6)
      r2(i,j) = rioh*feh*nel(i,j)*narr(i,j,6)
      lnt = alog(tarr(i,j,5))
      A = 1.1
      y = -26.0+A*cos((lnt*(1.2*!pi)/(alog(10)) + 0.74*!pi))
      rro2p = exp(y)
      r3(i,j) = narr(i,j,4)*nel(i,j)*rro2p 
      r4(i,j) = 2.0*narr(i,j,6)*narr(i,j,4)*cx_k6
      r5(i,j) = narr(i,j,6)*narr(i,j,0)*cx_k8
      r6(i,j) = narr(i,j,5)*narr(i,j,4)*cx_k10
      r7(i,j) = narr(i,j,5)*narr(i,j,4)*cx_k11
      r8(i,j) = narr(i,j,6)*narr(i,j,1)*cx_k12
      r9(i,j) = narr(i,j,4)*narr(i,j,0)*cx_k13
      r10(i,j) = narr(i,j,6)*narr(i,j,2)*cx_k14
      r11(i,j) = narr(i,j,4)*narr(i,j,1)*cx_k15
   endfor
endfor

r_tot = r1 + r2 + r3 + r4 + r5 + r6 + r7 + r8 + r9 + r10 + r11

r1 = 100.*r1/r_tot
r2 = 100.*r2/r_tot
r3 = 100.*r3/r_tot
r4 = 100.*r4/r_tot
r5 = 100.*r5/r_tot
r6 = 100.*r6/r_tot
r7 = 100.*r7/r_tot
r8 = 100.*r8/r_tot
r9 = 100.*r9/r_tot
r10 = 100.*r10/r_tot
r11 = 100.*r11/r_tot

;flr(whx,why) = r_tot(whx,why)

!p.multi=[0,3,4]

contour,r1/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,0,1],$
;   levels=[8,10,12,14,16,18,20],$
   title='O + e',c_charsize=1.0,nlev=10
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r2/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O + e(hot)',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r3/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O++ + e -> O+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r4/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O + O++ -> O+ + O+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r5/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O + S+ -> O+ + S',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r6/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S + O++ -> S+ + O+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r7/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S + O++ -> S++ + O+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r8/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O + S++ -> O+ + S+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r9/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O++ + S+ -> O+ + S++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r10/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O + S+++ -> O+ + S++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r11/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O++ + S++ -> O+ + S+++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1



return
end
;-------------------------------------------------------------------



;-------------------------------------------------------------------
pro op_loss,feh,teh
;-------------------------------------------------------------------
;set_plot,'ps
;@x6x9
;device,/encapsulated
restore,'v1_teh600_feh0023_otos4.0.sav

min_ns = 10e28
max_ns = 50e28
!x.margin=[4,2]
!y.margin=[2,2]
;!y.title='Transport time (days)'
;!x.title='source rate (10!u-4!n cm!u-3!n s!u-1!n)'
!p.charsize=1.0
!p.font=-1
!x.range=[4,20]
!y.range=[15,60]
!y.style = 1
!x.style = 1
Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
xarr = (1e4*xarr/vol)
print,'xarr...',xarr
yarr = (1./yarr)/(3600.*24.)

cx_k0 = 8.1e-9       ;S+ + S++ -> S++ + S+
cx_k1 = 2.4e-8
cx_k2 = 3.0e-10
cx_k3 = 7.8e-9
cx_k4 = 1.32e-8
cx_k5 = 1.32e-8
cx_k6 = 5.2e-10
cx_k7 = 5.4e-9
;r.cx_k8 = k8*9.0e-11   ;Schreier
cx_k8 = 6.0e-11    ;McGrath
cx_k9 = 3.0e-9
cx_k10 = 2.34e-8
cx_k11 = 1.62e-8
cx_k12 = 2.3e-9  ;new
;r.cx_k12 = k12*7.3e-9   ;old
cx_k13 = 1.4e-9   ;Schreier for L=6
cx_k14 = 1.92e-8 ;Schreier
;r.cx_k14 = k14*4.3e-9 ;Barbosa94, Johnson82
;r.cx_k15 = k15*9e-9    ;Schreier for L=6
cx_k15 = 9e-10    ;McGrath for L=6
cx_k16 = 3.6e-10  ;McGrath


nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)

sz = size(narr)
wh = where(nel ge 1.0 )
whx = wh mod sz(1)
why = wh/sz(1)

flr = fltarr(sz(1),sz(2))
flr(whx,why) = 1.0

nel=nel/flr

;device,filename='idl.ps'
;!p.multi=[0,1,1]
;surface,flr
;device,/close

x = findgen(130)+10.
xx = 1e4*x*1e27/2.5e31

;nominal voyager
mdn = (60.0-45.0)/(30.0-50.0)
y0dn = 60.0 - mdn*30.0

;mdn = (25.0-50.0)/(55.0-30.0)
;y0dn = 55.0 - mdn*25.0

ydn = mdn*x + y0dn
whdn = where((ydn ge 45.0) and (ydn le 60))

;device,xsize=6.2,ysize=7.2,/inches,xoffset=1.0,yoffset=3.0
;device,filename='v1_teh600_feh003_otos2.5.eps
;device,filename='v1_test.eps
;device,/color
;!p.multi=[0,4,5]

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)

sz = size(nel)
wh = where(nel le 0.1)
whx = wh mod sz(1)
why = wh/sz(2)
nelec = nel
nelec(whx,why) = 0.000001

r1 = fltarr(sz(1),sz(2)) ;thermal electron impact ionization
r2 = fltarr(sz(1),sz(2)) ;hot electron impact ionization
r3 = fltarr(sz(1),sz(2)) ;k2
r4 = fltarr(sz(1),sz(2)) ;k4
r5 = fltarr(sz(1),sz(2)) ;k9
r6 = fltarr(sz(1),sz(2)) ;k10
r7 = fltarr(sz(1),sz(2)) ;k12
r8 = fltarr(sz(1),sz(2)) ;recombination
r9 = fltarr(sz(1),sz(2)) ;recombination
r10 = fltarr(sz(1),sz(2)) ;recombination
r11 = fltarr(sz(1),sz(2)) ;recombination

r_tot = fltarr(sz(1),sz(2)) ;total reaction rates

for i = 0,sz(1)-1 do begin
   for j = 0,sz(2)-1 do begin
      riop = 3.78e-9*sqrt(tarr(i,j,0))*exp(-33.42/tarr(i,j,0))
      rioph = 3.78e-9*sqrt(teh)*exp(-33.42/teh)
      r1(i,j) = riop*nel(i,j)*narr(i,j,3)
      r2(i,j) = rioph*feh*nel(i,j)*narr(i,j,3)
      r3(i,j) = narr(i,j,5)*narr(i,j,3)*cx_k9
;      r6(i,j) = narr(i,j,5)*narr(i,j,4)*cx_k10
;      r7(i,j) = narr(i,j,5)*narr(i,j,4)*cx_k11
;      r8(i,j) = narr(i,j,6)*narr(i,j,1)*cx_k12
;      r9(i,j) = narr(i,j,4)*narr(i,j,0)*cx_k13
;      r10(i,j) = narr(i,j,6)*narr(i,j,2)*cx_k14
;      r11(i,j) = narr(i,j,4)*narr(i,j,1)*cx_k15
   endfor
endfor

r_tot = r1 + r2 + r3 

r1 = 100.*r1/r_tot
r2 = 100.*r2/r_tot
r3 = 100.*r3/r_tot
;r4 = 100.*r4/r_tot
;r5 = 100.*r5/r_tot
;r6 = 100.*r6/r_tot
;r7 = 100.*r7/r_tot
;r8 = 100.*r8/r_tot
;r9 = 100.*r9/r_tot
;r10 = 100.*r10/r_tot
;r11 = 100.*r11/r_tot

;flr(whx,why) = r_tot(whx,why)

!p.multi=[0,3,4]

contour,r1/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,0,1],$
;   levels=[8,10,12,14,16,18,20],$
   title='O+ + e',c_charsize=1.0,nlev=10
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r2/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O+ + e(hot)',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r3/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S + O+ -> S+ + O',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r4/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O + O++ -> O+ + O+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r5/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O + S+ -> O+ + S',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r6/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='S + O++ -> S+ + O+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r7/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='S + O++ -> S++ + O+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r8/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O + S++ -> O+ + S+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1


return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro o2p_production,feh,teh
;-------------------------------------------------------------------
;set_plot,'ps
;@x6x9
;device,/encapsulated
restore,'v1_teh600_feh0023_otos4.0.sav

min_ns = 10e28
max_ns = 50e28
!x.margin=[4,2]
!y.margin=[2,2]
;!y.title='Transport time (days)'
;!x.title='source rate (10!u-4!n cm!u-3!n s!u-1!n)'
!p.charsize=1.0
!p.font=-1
!x.range=[4,20]
!y.range=[15,60]
!y.style = 1
!x.style = 1
Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
xarr = (1e4*xarr/vol)
print,'xarr...',xarr
yarr = (1./yarr)/(3600.*24.)

cx_k0 = 8.1e-9       ;S+ + S++ -> S++ + S+
cx_k1 = 2.4e-8
cx_k2 = 3.0e-10
cx_k3 = 7.8e-9
cx_k4 = 1.32e-8
cx_k5 = 1.32e-8
cx_k6 = 5.2e-10
cx_k7 = 5.4e-9
;r.cx_k8 = k8*9.0e-11   ;Schreier
cx_k8 = 6.0e-11    ;McGrath
cx_k9 = 3.0e-9
cx_k10 = 2.34e-8
cx_k11 = 1.62e-8
cx_k12 = 2.3e-9  ;new
;r.cx_k12 = k12*7.3e-9   ;old
cx_k13 = 1.4e-9   ;Schreier for L=6
cx_k14 = 1.92e-8 ;Schreier
;r.cx_k14 = k14*4.3e-9 ;Barbosa94, Johnson82
;r.cx_k15 = k15*9e-9    ;Schreier for L=6
cx_k15 = 9e-10    ;McGrath for L=6
cx_k16 = 3.6e-10  ;McGrath


nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)

sz = size(narr)
wh = where(nel ge 1.0 )
whx = wh mod sz(1)
why = wh/sz(1)

flr = fltarr(sz(1),sz(2))
flr(whx,why) = 1.0

nel=nel/flr

;device,filename='idl.ps'
;!p.multi=[0,1,1]
;surface,flr
;device,/close

x = findgen(130)+10.
xx = 1e4*x*1e27/2.5e31

;nominal voyager
mdn = (60.0-45.0)/(30.0-50.0)
y0dn = 60.0 - mdn*30.0

;mdn = (25.0-50.0)/(55.0-30.0)
;y0dn = 55.0 - mdn*25.0

ydn = mdn*x + y0dn
whdn = where((ydn ge 45.0) and (ydn le 60))

;device,xsize=6.2,ysize=7.2,/inches,xoffset=1.0,yoffset=3.0
;device,filename='v1_teh600_feh003_otos2.5.eps
;device,filename='v1_test.eps
;device,/color
;!p.multi=[0,4,5]

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)

sz = size(nel)
wh = where(nel le 0.1)
whx = wh mod sz(1)
why = wh/sz(2)
nelec = nel
nelec(whx,why) = 0.000001

r1 = fltarr(sz(1),sz(2)) ;thermal electron impact ionization
r2 = fltarr(sz(1),sz(2)) ;hot electron impact ionization
r3 = fltarr(sz(1),sz(2)) ;k2
r4 = fltarr(sz(1),sz(2)) ;k4
r5 = fltarr(sz(1),sz(2)) ;k9
r6 = fltarr(sz(1),sz(2)) ;k10
r7 = fltarr(sz(1),sz(2)) ;k12
r8 = fltarr(sz(1),sz(2)) ;recombination
r9 = fltarr(sz(1),sz(2)) ;recombination
r10 = fltarr(sz(1),sz(2)) ;recombination
r11 = fltarr(sz(1),sz(2)) ;recombination

r_tot = fltarr(sz(1),sz(2)) ;total reaction rates

for i = 0,sz(1)-1 do begin
   for j = 0,sz(2)-1 do begin
      riop = 3.78e-9*sqrt(tarr(i,j,0))*exp(-33.42/tarr(i,j,0))
      rioph = 3.78e-9*sqrt(teh)*exp(-33.42/teh)
      r1(i,j) = riop*nel(i,j)*narr(i,j,3)
      r2(i,j) = rioph*feh*nel(i,j)*narr(i,j,3)
;      r3(i,j) = narr(i,j,5)*narr(i,j,3)*cx_k9
;      r6(i,j) = narr(i,j,5)*narr(i,j,4)*cx_k10
;      r7(i,j) = narr(i,j,5)*narr(i,j,4)*cx_k11
;      r8(i,j) = narr(i,j,6)*narr(i,j,1)*cx_k12
;      r9(i,j) = narr(i,j,4)*narr(i,j,0)*cx_k13
;      r10(i,j) = narr(i,j,6)*narr(i,j,2)*cx_k14
;      r11(i,j) = narr(i,j,4)*narr(i,j,1)*cx_k15
   endfor
endfor

r_tot = r1 + r2 

r1 = 100.*r1/r_tot
r2 = 100.*r2/r_tot
r3 = 100.*r3/r_tot
;r4 = 100.*r4/r_tot
;r5 = 100.*r5/r_tot
;r6 = 100.*r6/r_tot
;r7 = 100.*r7/r_tot
;r8 = 100.*r8/r_tot
;r9 = 100.*r9/r_tot
;r10 = 100.*r10/r_tot
;r11 = 100.*r11/r_tot

;flr(whx,why) = r_tot(whx,why)

!p.multi=[0,3,4]

contour,r1/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,0,1],$
;   levels=[8,10,12,14,16,18,20],$
   title='O+ + e',c_charsize=1.0,nlev=10
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r2/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O+ + e(hot)',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r3/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='S + O+ -> S+ + O',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r4/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O + O++ -> O+ + O+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r5/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O + S+ -> O+ + S',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r6/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='S + O++ -> S+ + O+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r7/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='S + O++ -> S++ + O+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r8/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O + S++ -> O+ + S+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1


return
end
;-------------------------------------------------------------------



;-------------------------------------------------------------------
pro o2p_loss,feh,teh
;-------------------------------------------------------------------
;set_plot,'ps
;@x6x9
;device,/encapsulated
restore,'v1_teh600_feh0023_otos4.0.sav

min_ns = 10e28
max_ns = 50e28
!x.margin=[4,2]
!y.margin=[2,2]
;!y.title='Transport time (days)'
;!x.title='source rate (10!u-4!n cm!u-3!n s!u-1!n)'
!p.charsize=1.0
!p.font=-1
!x.range=[4,20]
!y.range=[15,60]
!y.style = 1
!x.style = 1
Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
xarr = (1e4*xarr/vol)
print,'xarr...',xarr
yarr = (1./yarr)/(3600.*24.)

cx_k0 = 8.1e-9       ;S+ + S++ -> S++ + S+
cx_k1 = 2.4e-8
cx_k2 = 3.0e-10
cx_k3 = 7.8e-9
cx_k4 = 1.32e-8
cx_k5 = 1.32e-8
cx_k6 = 5.2e-10
cx_k7 = 5.4e-9
;r.cx_k8 = k8*9.0e-11   ;Schreier
cx_k8 = 6.0e-11    ;McGrath
cx_k9 = 3.0e-9
cx_k10 = 2.34e-8
cx_k11 = 1.62e-8
cx_k12 = 2.3e-9  ;new
;r.cx_k12 = k12*7.3e-9   ;old
cx_k13 = 1.4e-9   ;Schreier for L=6
cx_k14 = 1.92e-8 ;Schreier
;r.cx_k14 = k14*4.3e-9 ;Barbosa94, Johnson82
;r.cx_k15 = k15*9e-9    ;Schreier for L=6
cx_k15 = 9e-10    ;McGrath for L=6
cx_k16 = 3.6e-10  ;McGrath


nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)

sz = size(narr)
wh = where(nel ge 1.0 )
whx = wh mod sz(1)
why = wh/sz(1)

flr = fltarr(sz(1),sz(2))
flr(whx,why) = 1.0

nel=nel/flr

;device,filename='idl.ps'
;!p.multi=[0,1,1]
;surface,flr
;device,/close

x = findgen(130)+10.
xx = 1e4*x*1e27/2.5e31

;nominal voyager
mdn = (60.0-45.0)/(30.0-50.0)
y0dn = 60.0 - mdn*30.0

;mdn = (25.0-50.0)/(55.0-30.0)
;y0dn = 55.0 - mdn*25.0

ydn = mdn*x + y0dn
whdn = where((ydn ge 45.0) and (ydn le 60))

;device,xsize=6.2,ysize=7.2,/inches,xoffset=1.0,yoffset=3.0
;device,filename='v1_teh600_feh003_otos2.5.eps
;device,filename='v1_test.eps
;device,/color
;!p.multi=[0,4,5]

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)

sz = size(nel)
wh = where(nel le 0.1)
whx = wh mod sz(1)
why = wh/sz(2)
nelec = nel
nelec(whx,why) = 0.000001

r1 = fltarr(sz(1),sz(2)) ;thermal electron impact ionization
r2 = fltarr(sz(1),sz(2)) ;hot electron impact ionization
r3 = fltarr(sz(1),sz(2)) ;k2
r4 = fltarr(sz(1),sz(2)) ;k4
r5 = fltarr(sz(1),sz(2)) ;k9
r6 = fltarr(sz(1),sz(2)) ;k10
r7 = fltarr(sz(1),sz(2)) ;k12
r8 = fltarr(sz(1),sz(2)) ;recombination
r9 = fltarr(sz(1),sz(2)) ;recombination
r10 = fltarr(sz(1),sz(2)) ;recombination
r11 = fltarr(sz(1),sz(2)) ;recombination

r_tot = fltarr(sz(1),sz(2)) ;total reaction rates

for i = 0,sz(1)-1 do begin
   for j = 0,sz(2)-1 do begin
      lnt = alog(tarr(i,j,0))
      A = 1.1
      y = -26.0+A*cos((lnt*(1.2*!pi)/(alog(10)) + 0.74*!pi))
      rro2p = exp(y)
      r1(i,j) = rro2p*nel(i,j)*narr(i,j,4)
      r2(i,j) = narr(i,j,6)*narr(i,j,4)*cx_k6
      r3(i,j) = narr(i,j,5)*narr(i,j,4)*cx_k10
      r4(i,j) = narr(i,j,5)*narr(i,j,4)*cx_k11
      r5(i,j) = narr(i,j,4)*narr(i,j,0)*cx_k13
      r6(i,j) = narr(i,j,4)*narr(i,j,1)*cx_k15
   endfor
endfor

r_tot = r1 + r2 + r3 + r4 + r5 + r6

r1 = 100.*r1/r_tot
r2 = 100.*r2/r_tot
r3 = 100.*r3/r_tot
r4 = 100.*r4/r_tot
r5 = 100.*r5/r_tot
r6 = 100.*r6/r_tot
;r7 = 100.*r7/r_tot
;r8 = 100.*r8/r_tot
;r9 = 100.*r9/r_tot
;r10 = 100.*r10/r_tot
;r11 = 100.*r11/r_tot

;flr(whx,why) = r_tot(whx,why)

!p.multi=[0,3,4]

contour,r1/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,0,1],$
;   levels=[8,10,12,14,16,18,20],$
   title='O++ + e -> O+',c_charsize=1.0,nlev=10
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r2/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O + O++ -> O+ -> O+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r3/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S + O++ -> S+ + O+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r4/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='S + O++ -> S++ + O+',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r5/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O++ + S+ -> O+ + S++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,r6/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;   levels=[30,40,50,60,70,80,100,120,140,160],$
   title='O++ + S++ -> O+ + S+++',c_charsize=1.0,nlev=10 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r7/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='S + O++ -> S++ + O+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,r8/flr,xarr,yarr,/xsty,$
;;   C_LABELS=[1,1,1,1,1,1,1,1,0,1,1,1,1,1,1],$
;;   levels=[30,40,50,60,70,80,100,120,140,160],$
;   title='O + S++ -> O+ + S+',c_charsize=1.0,nlev=10 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1


return
end
;-------------------------------------------------------------------


set_plot,'ps
@x6x9

SP_PRODUCTION,0.0023,600.
SP_LOSS,0.0023,600.
S2P_PRODUCTION,0.0023,600.
S2P_LOSS,0.0023,600.
S3P_PRODUCTION,0.0023,600.
S3P_LOSS,0.0023,600.
OP_PRODUCTION,0.0023,600.
OP_LOSS,0.0023,600.
O2P_PRODUCTION,0.0023,600.
O2P_LOSS,0.0023,600.

device,/close
end
