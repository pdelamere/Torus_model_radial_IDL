;-------------------------------------------------------------------
pro make_cont_src_tau_cassini
;-------------------------------------------------------------------
set_plot,'ps
@x6x9
device,/encapsulated
!x.margin=[4,2]
!y.margin=[2,2]
restore,'cm3_2d_out.sav
;restore,'cassini_teh40_feh003_otos1.7.sav
;restore,'cm3_2d_teh40_otos3_fehsrc.sav'

min_ns = 1e28
max_ns = 3e28
;!y.title='Transport time (days)'
;!x.title='source rate (10!u-4!n cm!u-3!n s!u-1!n)'
!p.charsize=1.0
!p.font=0
!x.range=[8,24]
!y.type = 1
!y.range=[12,40]
!y.style = 1
!x.style = 1
Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
;print,vol
xarr = (1e4*xarr/vol)
;print,'xarr...',xarr
yarr = (1./yarr)/(3600.*24.)

;narr(*,*,0) = smooth(narr(*,*,0),2)
;narr(*,*,1) = smooth(narr(*,*,1),2)
;narr(*,*,2) = smooth(narr(*,*,2),2)
;narr(*,*,3) = smooth(narr(*,*,3),2)
;narr(*,*,4) = smooth(narr(*,*,4),2)

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)

nel = nel + 0.1*nel

!p.multi=[0,4,4]

sz = size(narr)
wh = where(nel gt 0.0)
whx = wh mod sz(1)
why = wh/sz(1)

flr = fltarr(sz(1),sz(2))
flr(whx,why) = 1.0

y0up = 15
x0up = 140.0
aup = (90.0 - y0up)/(10.0-x0up)^2

x = findgen(130)+10.
yup = aup*(x-x0up)^2 + y0up
xx = 1e4*x*1e27/2.5e31
;print,xx,yup
;stop

mdn = (25.0-50.0)/(55.0-30.0)
y0dn = 55.0 - mdn*25.0

;mdn = (6.0-15.0)/(50.0-10.0)
;y0dn = 15.0 - mdn*10.0


ydn = mdn*x + y0dn
whdn = where((ydn ge 30.0) and (ydn le 60))
;print,xx,ydn
;stop

set_plot,'ps
device,xsize=6.2,ysize=7.2,/inches
;device,filename='cassini_teh40_feh003_otos1.7.eps
!p.multi=[0,4,5]


contour,narr(*,*,5)/flr,xarr,yarr,/xsty,$
  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[2,5,6,7,8,9,10,12,20,25,30,35],$
      title='S',c_charsize=0.5,nlev=10 
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1


contour,narr(*,*,6)/flr,xarr,yarr,/xsty,$
  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[10,20,30,40,50,80,120,160,180,220,260],$
      title='O',c_charsize=0.5,nlev=10 
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1


nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)
nel = nel + 0.1*nel
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)


;contour,nel,xarr,yarr,/xsty,$
;C_LABELS=[1,1],levels=[1650,2000],c_color=[200,255],$
;title='n!de!n',c_charsize=0.5,nlev=10,/fill   
contour,nel/flr,xarr,yarr,/xsty,$
C_LABELS=[1,1,0,1,0,1,0,1,0,1,0],levels=[500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500],$
title='n!de!n',c_charsize=0.5,nlev=10 
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

sz = size(nel)
;wh = where(nel le 0.0)
;whx = wh mod sz(1)
;why = wh/sz(2)
nelec = nel
;nelec(whx,why) = 0.00001
contour,narr(*,*,0)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],$
   levels=[0.095-0.0095,0.095+0.0095],c_color=[200,255],$
   title='S!u+!n/n!de!n',c_charsize=0.5,nlev=10,/fill
contour,narr(*,*,0)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.02,0.04,0.06,0.08,0.1,0.12,0.15,0.20,0.30,0.40],$
   title='S!u+!n/n!de!n',c_charsize=0.5,nlev=10,/overplot 
;plots,xx(whdn),ydn(whdn),/data,linestyle=1


contour,narr(*,*,1)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],levels=[0.228-0.0228,0.228+0.0228],$
   title='S!u++!n/n!de!n',c_charsize=0.5,nlev=10,c_color=[200,255],/fill   
contour,narr(*,*,1)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,0,1,1,1,1,1,1,1,1,1],levels=[0.1,0.14,0.18,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27],$
   title='S!u++!n/n!de!n',c_charsize=0.5,nlev=10,/overplot    
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1



contour,narr(*,*,2)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],levels=[0.0223-0.00223,0.0223+0.00223],$
   title='S!u+++!n/n!de!n',c_charsize=0.5,nlev=10,c_color=[200,255],/fill      
contour,narr(*,*,2)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.001,0.005,0.01,0.015,0.02,0.030,0.04,0.05],$
   title='S!u+++!n/n!de!n',c_charsize=0.5,nlev=10,/overplot    
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1


contour,narr(*,*,3)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],levels=[0.23-0.023,0.23+0.023],$
   title='O!u+!n/n!de!n',c_charsize=0.5,nlev=10,c_color=[200,255],/fill         
contour,narr(*,*,3)/nel,xarr,yarr,/xsty,$   
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],levels=[0.2,0.22,0.24,0.26,0.28,0.30,0.32],$
   title='O!u+!n/n!de!n',c_charsize=0.5,nlev=10,/overplot    
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,narr(*,*,4)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],levels=[0.023-0.0023,0.023+0.0023],$
   title='O!u++!n/n!de!n',c_charsize=0.5,nlev=10,c_color=[200,255],/fill           
contour,narr(*,*,4)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,0,1,0,1,0,1,0,1],$
   levels=[0.005,0.010,0.015,0.02,0.025,0.030,0.035,0.04,0.045,0.050],$
   title='O!u++!n/n!de!n',c_charsize=0.5,nlev=10,/overplot   
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

sz = size(nel)
;wh = where(nel le 0.0)
;whx = wh mod sz(1)
;why = wh/sz(2)
narr1 = narr(*,*,1)
;narr1(whx,why) = 0.00001

contour,(narr(*,*,0)/narr1),xarr,yarr,/xsty,$
    C_LABELS=[1,1],levels=[0.413-0.0413,0.413+0.0413],c_color=[200,255],$
    title='S!u+!n/S!u++!n',c_charsize=0.5,nlev=10,/fill 
contour,(narr(*,*,0)/narr(*,*,1))/flr,xarr,yarr,/xsty,$
    C_LABELS=[1,1,1,1,1,1,1,1],levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0],$
    title='S!u+!n/S!u++!n',c_charsize=0.5,nlev=10,/overplot 
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,(narr(*,*,2)/narr1),xarr,yarr,/xsty,$
    C_LABELS=[1,1],levels=[0.098-0.0098,0.098+0.0098],c_color=[200,255],$
    title='S!u+++!n/S!u++!n',c_charsize=0.5,nlev=10,/fill 
contour,(narr(*,*,2)/narr(*,*,1))/flr,xarr,yarr,/xsty,$
    C_LABELS=[1,1,1,1,1,1,1,1],levels=[0.05,0.1,0.15,0.2,0.25,0.3],$
    title='S!u+++!n/S!u++!n',c_charsize=0.5,nlev=10,/overplot 
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

sonp = narr(*,*,3)+narr(*,*,4)
ssnp = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)

contour,sonp/ssnp,xarr,yarr,/xsty,$
C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1,1],levels=[0.8,0.85,0.9,0.92,0.94,0.96,0.98,1.0],$
title='!9S!3O!un+!n/!9S!3S!un+!n',c_charsize=0.5,nlev=10   
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1
print,sonp/ssnp


contour,nel/sni,xarr,yarr,/xsty,$
C_LABELS=[1,0,1,0,1,0,1,1,1,1],levels=[1.1,1.15,1.2,1.25,1.30,1.35,1.4,1.45,1.5,1.55,1.60],$
title='n!de!n/!9S!3n!di!n',c_charsize=0.5,nlev=10   
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1



contour,tarr(*,*,0)/flr,xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1],levels=[4.0,4.4,4.8,5.2,5.6,6.0,6.4,6.8],$
   title='T!de!n (eV)',c_charsize=0.5,nlev=10   
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

sti = tarr(*,*,1)*narr(*,*,0)+tarr(*,*,2)*narr(*,*,1)+$
      tarr(*,*,3)*narr(*,*,2)+tarr(*,*,4)*narr(*,*,3)+$
      tarr(*,*,5)*narr(*,*,4)
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)

contour,(sti/sni)/flr,xarr,yarr,/xsty,$
    C_LABELS=[1,1,1,1,1,1,1,1],levels=[40,80,120,160,200,240],$
    title='T!di!n (eV)',c_charsize=0.5,nlev=10 
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

ion_cx = (ecarr(*,*,0)+ecarr(*,*,1)+ecarr(*,*,3)+ecarr(*,*,4))/(ecarr(*,*,2)+ecarr(*,*,5))

contour,ion_cx/flr,xarr,yarr,/xsty,$
    C_LABELS=[1,1,1,1,1,1,1,1,1,1],$
    levels=[0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],$
  title='P!dion!n/P!dcx!n', c_charsize=0.5,nlev=10 
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,ecarr(*,*,12)/flr,xarr,yarr,/xsty,$
  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7],$
  title='P!dfast!n ', c_charsize=0.5,nlev=10 
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,(ecarr(*,*,13)+ecarr(*,*,14))/flr,xarr,yarr,/xsty,$
  C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1],levels=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6],$
  title='P!dtrans!n ', c_charsize=0.5,nlev=10 
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,ecarr(*,*,17)/flr,xarr,yarr,/xsty,$
    C_LABELS=[1,1,1,1,0,1,0,1,1,1],levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8],$
  title='P!dion->e!n ', c_charsize=0.5,nlev=10 
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1

puv_obs = 2.0006927e+12/1.6e-19/1.5e31
dpuv_obs = 5.44e11/1.6e-19/1.5e31

contour,ecarr(*,*,11),xarr,yarr,/xsty,$
    C_LABELS=[1,1],levels=[puv_obs-dpuv_obs,puv_obs+dpuv_obs],c_color=[200,255],$
    title='P!duv!n',c_charsize=0.5,nlev=10,/fill 
contour,ecarr(*,*,11)/flr,xarr,yarr,/xsty,$
  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[0.2,0.4,0.6,1.0,1.4,1.8,2.2],$
  title='P!duv!n ', c_charsize=0.5,nlev=10,/overplot
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1


contour,ecarr(*,*,9)/flr,xarr,yarr,/xsty,$
    C_LABELS=[1,1,0,1,1,0,1,1,1,1,1],$
    levels=[0.01,0.1,0.2,0.4,0.6,0.8,1.0,2.0],$
  title='P!deh!n ', c_charsize=0.5,nlev=10 
;plots,xx,yup,/data,linestyle=1
;plots,xx(whdn),ydn(whdn),/data,linestyle=1






xyouts,2500.*0.0,2500.*3.5,'Transport time (days)',$
  orientation=90.0,/device,charsize=0.7
xyouts,2500.*2.5,2500.*(-0.2),'Source rate (10!u-4!n cm!u-3!n s!u-1!n)',$
  /device,charsize=0.7
device,/close

device,filename='cassini_2plots.eps'
;!p.multi=[0,2,1]
contour,narr(*,*,0)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],$
   levels=[0.05,0.07],c_color=[200,255],$
   title='S!u+!n/n!de!n',c_charsize=0.5,nlev=10,/fill
contour,narr(*,*,0)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.02,0.04,0.06,0.08,0.1,0.12,0.15,0.20,0.30,0.40],$
   title='S!u+!n/n!de!n',c_charsize=0.5,nlev=10,/overplot 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,ecarr(*,*,11),xarr,yarr,/xsty,$
    C_LABELS=[1,1],levels=[0.4,0.8],c_color=[200,255],$
    title='P!duv!n',c_charsize=0.5,nlev=10,/fill 
contour,ecarr(*,*,11)/flr,xarr,yarr,/xsty,$
  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[0.2,0.4,0.6,1.0,1.4,1.8,2.2],$
  title='P!duv!n ', c_charsize=0.5,nlev=10,/overplot
;plots,xx,yup,/data,linestyle=1
plots,xx(whdn),ydn(whdn),/data,linestyle=1

device,/close


return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro make_cont_src_tau_6up_oct
;-------------------------------------------------------------------
set_plot,'ps
;@x6x9
device,/encapsulated
device,xsize = 3.3, ysize = 7.0,/inches,xoffset=1.0,yoffset=1.0
!x.margin=[4,2]
!y.margin=[2,2]
restore,'cm3_2d_cassini_oct.sav
;restore,'cm3_2d_out.sav
;restore,'cassini_teh40_feh003_otos1.7.sav
;restore,'cm3_2d_teh40_otos3_fehsrc.sav'

min_ns = 1e28
max_ns = 3e28
;!y.title='Transport time (days)'
;!x.title='source rate (10!u-4!n cm!u-3!n s!u-1!n)'
!p.charsize=1.2
!p.font=0
!x.range=[4,19]
!y.type = 1
!y.range=[20,68]
!y.style = 1
!x.style = 1
Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
;print,vol
xarr = (1e4*xarr/vol)
;print,'xarr...',xarr
yarr = (1./yarr)/(3600.*24.)

;narr(*,*,0) = smooth(narr(*,*,0),2)
;narr(*,*,1) = smooth(narr(*,*,1),2)
;narr(*,*,2) = smooth(narr(*,*,2),2)
;narr(*,*,3) = smooth(narr(*,*,3),2)
;narr(*,*,4) = smooth(narr(*,*,4),2)

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)

nel = nel + 0.1*nel

!p.multi=[0,2,3]

sz = size(narr)
wh = where(nel gt 0.0)
whx = wh mod sz(1)
why = wh/sz(1)

flr = fltarr(sz(1),sz(2))
flr(whx,why) = 1.0

y0up = 15
x0up = 140.0
aup = (90.0 - y0up)/(10.0-x0up)^2

x = findgen(130)+10.
yup = aup*(x-x0up)^2 + y0up
xx = 1e4*x*1e27/2.5e31
;print,xx,yup
;stop

mdn = (40.0-65.0)/(45.0-20.0)
y0dn = 40. - mdn*45.0


;mdn = (6.0-15.0)/(50.0-10.0)
;y0dn = 15.0 - mdn*10.0


ydn = mdn*x + y0dn
whdn = where((ydn ge 38.0) and (ydn le 68))
;print,xx,ydn
;stop

set_plot,'ps
device,xsize=4.8,ysize=7.0,/inches
device,filename='cassini_oct.eps
!p.multi=[0,2,3]
!p.charsize=1.6


;contour,narr(*,*,5)/flr,xarr,yarr,/xsty,$
;  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[2,5,6,7,8,9,10,12,20,25,30,35],$
;      title='S',c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1


;contour,narr(*,*,6)/flr,xarr,yarr,/xsty,$
;  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[10,20,30,40,50,80,120,160,180,220,260],$
;      title='O',c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1


nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)
nel = nel + 0.1*nel
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)


;;contour,nel,xarr,yarr,/xsty,$
;;C_LABELS=[1,1],levels=[1650,2000],c_color=[200,255],$
;;title='n!de!n',c_charsize=0.8,nlev=10,/fill   
;contour,nel/flr,xarr,yarr,/xsty,$
;C_LABELS=[1,1,0,1,0,1,0,1,0,1,0],levels=[500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500],$
;title='n!de!n',c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

sz = size(nel)
;wh = where(nel le 0.0)
;whx = wh mod sz(1)
;why = wh/sz(2)
nelec = nel
;nelec(whx,why) = 0.00001


vol=3.5e31
puv_obs = 2.0006927e+12/1.6e-19/vol
dpuv_obs = 5.44e11/1.6e-19/vol



contour,narr(*,*,0)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],$
   levels=[0.099-0.0099,0.099+0.0099],c_color=[200,255],$
   title='S!u+!n/n!de!n',c_charsize=0.8,nlev=10,/cell_fill
contour,narr(*,*,0)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.02,0.04,0.06,0.08,0.1,0.12,0.15,0.20,0.30,0.40],$
   title='S!u+!n/n!de!n',c_charsize=0.8,nlev=10,/overplot 
plots,xx(whdn),ydn(whdn),/data,linestyle=1
;contour,ecarr(*,*,11),xarr,yarr,/xsty,$
;    C_LABELS=[1,1],levels=[puv_obs-dpuv_obs,puv_obs+dpuv_obs],$
;    c_charsize=0.8,nlev=10,/overplot,c_linestyle=1

tau1 = 27.0
sn1 = 17.2
dtau = 2.0
dsn=0.6
plots,sn1-dsn,tau1-dtau,/data
plots,sn1+dsn,tau1-dtau,/data,/continue
plots,sn1+dsn,tau1+dtau,/data,/continue
plots,sn1-dsn,tau1+dtau,/data,/continue
plots,sn1-dsn,tau1-dtau,/data,/continue



contour,narr(*,*,1)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],levels=[0.216-0.0216,0.216+0.0216],$
   title='S!u++!n/n!de!n',c_charsize=0.8,nlev=10,c_color=[200,255],/cell_fill   
contour,narr(*,*,1)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,0,1,1,1,1,1,1,1,1,1],levels=[0.1,0.14,0.18,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27],$
   title='S!u++!n/n!de!n',c_charsize=0.8,nlev=10,/overplot    
;plots,xx,yup,/data,linestyle=1
plots,xx(whdn),ydn(whdn),/data,linestyle=1
plots,sn1-dsn,tau1-dtau,/data
plots,sn1+dsn,tau1-dtau,/data,/continue
plots,sn1+dsn,tau1+dtau,/data,/continue
plots,sn1-dsn,tau1+dtau,/data,/continue
plots,sn1-dsn,tau1-dtau,/data,/continue

;contour,ecarr(*,*,11),xarr,yarr,/xsty,$
;    C_LABELS=[1,1],levels=[puv_obs-dpuv_obs,puv_obs+dpuv_obs],$
;    c_charsize=0.8,nlev=10,/overplot 



contour,narr(*,*,2)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],levels=[0.0218-0.00218,0.0218+0.00218],$
   title='S!u+++!n/n!de!n',c_charsize=0.8,nlev=10,c_color=[200,255],/cell_fill      
contour,narr(*,*,2)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.001,0.005,0.01,0.015,0.02,0.030,0.04,0.05],$
   title='S!u+++!n/n!de!n',c_charsize=0.8,nlev=10,/overplot    
;plots,xx,yup,/data,linestyle=1
plots,xx(whdn),ydn(whdn),/data,linestyle=1

plots,sn1-dsn,tau1-dtau,/data
plots,sn1+dsn,tau1-dtau,/data,/continue
plots,sn1+dsn,tau1+dtau,/data,/continue
plots,sn1-dsn,tau1+dtau,/data,/continue
plots,sn1-dsn,tau1-dtau,/data,/continue

;contour,ecarr(*,*,11),xarr,yarr,/xsty,$
;    C_LABELS=[1,1],levels=[puv_obs-dpuv_obs,puv_obs+dpuv_obs],$
;    c_charsize=0.8,nlev=10,/overplot 


contour,narr(*,*,3)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],levels=[0.26-0.026,0.26+0.026],$
   title='O!u+!n/n!de!n',c_charsize=0.8,nlev=10,c_color=[200,255],/cell_fill         
contour,narr(*,*,3)/nel,xarr,yarr,/xsty,$   
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],levels=[0.2,0.22,0.24,0.26,0.28,0.30,0.32],$
   title='O!u+!n/n!de!n',c_charsize=0.8,nlev=10,/overplot    
;plots,xx,yup,/data,linestyle=1
plots,xx(whdn),ydn(whdn),/data,linestyle=1

plots,sn1-dsn,tau1-dtau,/data
plots,sn1+dsn,tau1-dtau,/data,/continue
plots,sn1+dsn,tau1+dtau,/data,/continue
plots,sn1-dsn,tau1+dtau,/data,/continue
plots,sn1-dsn,tau1-dtau,/data,/continue

;contour,ecarr(*,*,11),xarr,yarr,/xsty,$
;    C_LABELS=[1,1],levels=[puv_obs-dpuv_obs,puv_obs+dpuv_obs],$
;    c_charsize=0.8,nlev=10,/overplot 

contour,narr(*,*,4)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],levels=[0.017-0.0017,0.017+0.0017],$
   title='O!u++!n/n!de!n',c_charsize=0.8,nlev=10,c_color=[200,255],/cell_fill           
contour,narr(*,*,4)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1],$
   levels=[0.014,0.016,0.018,0.020,0.024,0.028,0.04,0.045,0.050],$
   title='O!u++!n/n!de!n',c_charsize=0.8,/overplot   
;plots,xx,yup,/data,linestyle=1
plots,xx(whdn),ydn(whdn),/data,linestyle=1

plots,sn1-dsn,tau1-dtau,/data
plots,sn1+dsn,tau1-dtau,/data,/continue
plots,sn1+dsn,tau1+dtau,/data,/continue
plots,sn1-dsn,tau1+dtau,/data,/continue
plots,sn1-dsn,tau1-dtau,/data,/continue

;contour,ecarr(*,*,11),xarr,yarr,/xsty,$
;    C_LABELS=[1,1],levels=[puv_obs-dpuv_obs,puv_obs+dpuv_obs],$
;    c_charsize=0.8,nlev=10,/overplot 

sz = size(nel)
;wh = where(nel le 0.0)
;whx = wh mod sz(1)
;why = wh/sz(2)
narr1 = narr(*,*,1)
;narr1(whx,why) = 0.00001

;contour,(narr(*,*,0)/narr1),xarr,yarr,/xsty,$
;    C_LABELS=[1,1],levels=[0.413-0.0413,0.413+0.0413],c_color=[200,255],$
;    title='S!u+!n/S!u++!n',c_charsize=0.8,nlev=10,/fill 
;contour,(narr(*,*,0)/narr(*,*,1))/flr,xarr,yarr,/xsty,$
;    C_LABELS=[1,1,1,1,1,1,1,1],levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0],$
;    title='S!u+!n/S!u++!n',c_charsize=0.8,nlev=10,/overplot 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,(narr(*,*,2)/narr1),xarr,yarr,/xsty,$
;    C_LABELS=[1,1],levels=[0.098-0.0098,0.098+0.0098],c_color=[200,255],$
;    title='S!u+++!n/S!u++!n',c_charsize=0.8,nlev=10,/fill 
;contour,(narr(*,*,2)/narr(*,*,1))/flr,xarr,yarr,/xsty,$
;    C_LABELS=[1,1,1,1,1,1,1,1],levels=[0.05,0.1,0.15,0.2,0.25,0.3],$
;    title='S!u+++!n/S!u++!n',c_charsize=0.8,nlev=10,/overplot 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;sonp = narr(*,*,3)+narr(*,*,4)
;ssnp = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)

;contour,sonp/ssnp,xarr,yarr,/xsty,$
;C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1,1],levels=[0.8,0.85,0.9,0.92,0.94,0.96,0.98,1.0],$
;title='!9S!3O!un+!n/!9S!3S!un+!n',c_charsize=0.8,nlev=10   
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1
;print,sonp/ssnp


;contour,nel/sni,xarr,yarr,/xsty,$
;C_LABELS=[1,0,1,0,1,0,1,1,1,1],levels=[1.1,1.15,1.2,1.25,1.30,1.35,1.4,1.45,1.5,1.55,1.60],$
;title='n!de!n/!9S!3n!di!n',c_charsize=0.8,nlev=10   
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1



;contour,tarr(*,*,0)/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1],levels=[4.0,4.4,4.8,5.2,5.6,6.0,6.4,6.8],$
;   title='T!de!n (eV)',c_charsize=0.8,nlev=10   
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

sti = tarr(*,*,1)*narr(*,*,0)+tarr(*,*,2)*narr(*,*,1)+$
      tarr(*,*,3)*narr(*,*,2)+tarr(*,*,4)*narr(*,*,3)+$
      tarr(*,*,5)*narr(*,*,4)
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)

;contour,(sti/sni)/flr,xarr,yarr,/xsty,$
;    C_LABELS=[1,1,1,1,1,1,1,1],levels=[40,80,120,160,200,240],$
;    title='T!di!n (eV)',c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

ion_cx = (ecarr(*,*,0)+ecarr(*,*,1)+ecarr(*,*,3)+ecarr(*,*,4))/(ecarr(*,*,2)+ecarr(*,*,5))

;contour,ion_cx/flr,xarr,yarr,/xsty,$
;    C_LABELS=[1,1,1,1,1,1,1,1,1,1],$
;    levels=[0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],$
;  title='P!dion!n/P!dcx!n', c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,ecarr(*,*,12)/flr,xarr,yarr,/xsty,$
;  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7],$
;  title='P!dfast!n ', c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,(ecarr(*,*,13)+ecarr(*,*,14))/flr,xarr,yarr,/xsty,$
;  C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1],levels=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6],$
;  title='P!dtrans!n ', c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,ecarr(*,*,17)/flr,xarr,yarr,/xsty,$
;    C_LABELS=[1,1,1,1,0,1,0,1,1,1],levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8],$
;  title='P!dion->e!n ', c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1


puv_obs = 2.0006927e+12/1.6e-19/vol
dpuv_obs = 5.44e11/1.6e-19/vol

contour,ecarr(*,*,18),xarr,yarr,/xsty,$
    C_LABELS=[1,1],levels=[puv_obs-dpuv_obs,puv_obs+dpuv_obs],c_color=[200,255],$
    title='P!dEUV!n',c_charsize=0.8,nlev=10,/cell_fill 
contour,ecarr(*,*,18)/flr,xarr,yarr,/xsty,$
  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[0.1,0.2,0.3,0.4,0.5,1.0,1.8,2.2],$
  title='P!dEUV!n ', c_charsize=0.8,nlev=10,/overplot
;plots,xx,yup,/data,linestyle=1
plots,xx(whdn),ydn(whdn),/data,linestyle=1
plots,sn1-dsn,tau1-dtau,/data
plots,sn1+dsn,tau1-dtau,/data,/continue
plots,sn1+dsn,tau1+dtau,/data,/continue
plots,sn1-dsn,tau1+dtau,/data,/continue
plots,sn1-dsn,tau1-dtau,/data,/continue


;contour,ecarr(*,*,9)/flr,xarr,yarr,/xsty,$
;    C_LABELS=[1,1,0,1,1,0,1,1,1,1,1],$
;    levels=[0.01,0.1,0.2,0.4,0.6,0.8,1.0,2.0],$
;  title='P!deh!n ', c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1






xyouts,2500.*0.0,2500.*2.6,'Transport time (days)',$
  orientation=90.0,/device,charsize=1.2
xyouts,2500.*1.5,2500.*(-0.2),'Source rate (10!u-4!n cm!u-3!n s!u-1!n)',$
  /device,charsize=1.2
device,/close

device,filename='cassini_2plots.eps'
;!p.multi=[0,2,1]
contour,narr(*,*,0)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],$
   levels=[0.05,0.07],c_color=[200,255],$
   title='S!u+!n/n!de!n',c_charsize=0.5,nlev=10,/fill
contour,narr(*,*,0)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.02,0.04,0.06,0.08,0.1,0.12,0.15,0.20,0.30,0.40],$
   title='S!u+!n/n!de!n',c_charsize=0.5,nlev=10,/overplot 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,ecarr(*,*,11),xarr,yarr,/xsty,$
    C_LABELS=[1,1],levels=[0.4,0.8],c_color=[200,255],$
    title='P!duv!n',c_charsize=0.5,nlev=10,/fill 
contour,ecarr(*,*,11)/flr,xarr,yarr,/xsty,$
  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[0.2,0.4,0.6,1.0,1.4,1.8,2.2],$
  title='P!duv!n ', c_charsize=0.5,nlev=10,/overplot
;plots,xx,yup,/data,linestyle=1
plots,xx(whdn),ydn(whdn),/data,linestyle=1

device,/close


return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro make_cont_src_tau_6up_nov
;-------------------------------------------------------------------
set_plot,'ps
@x6x9
device,/encapsulated
!p.thick=2.0
!x.margin=[4,2]
!y.margin=[2,2]
restore,'cm3_2d_cassini_nov.sav
;restore,'cassini_teh40_feh003_otos1.7.sav
;restore,'cm3_2d_teh40_otos3_fehsrc.sav'

min_ns = 1e28
max_ns = 3e28
;!y.title='Transport time (days)'
;!x.title='source rate (10!u-4!n cm!u-3!n s!u-1!n)'
!p.charsize=1.2
!p.font=0
!x.range=[4,19]
!y.type = 1
!y.range=[20,68]
!y.style = 1
!x.style = 1
Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
;print,vol
xarr = (1e4*xarr/vol)
;print,'xarr...',xarr
yarr = (1./yarr)/(3600.*24.)

;narr(*,*,0) = smooth(narr(*,*,0),2)
;narr(*,*,1) = smooth(narr(*,*,1),2)
;narr(*,*,2) = smooth(narr(*,*,2),2)
;narr(*,*,3) = smooth(narr(*,*,3),2)
;narr(*,*,4) = smooth(narr(*,*,4),2)

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)

nel = nel + 0.1*nel

!p.multi=[0,2,3]

sz = size(narr)
wh = where(nel gt 0.0)
whx = wh mod sz(1)
why = wh/sz(1)

flr = fltarr(sz(1),sz(2))
flr(whx,why) = 1.0

y0up = 15
x0up = 140.0
aup = (90.0 - y0up)/(10.0-x0up)^2

x = findgen(130)+10.
yup = aup*(x-x0up)^2 + y0up
xx = 1e4*x*1e27/2.5e31
;print,xx,yup
;stop

;mdn = (25.0-50.0)/(55.0-30.0)
;y0dn = 55.0 - mdn*25.0

mdn = (50.0-60.0)/(45.0-35.0)
y0dn = 50. - mdn*45.0

;mdn = (6.0-15.0)/(50.0-10.0)
;y0dn = 15.0 - mdn*10.0


ydn = mdn*x + y0dn
whdn = where((ydn ge 48.0) and (ydn le 68))
;print,xx,ydn
;stop

set_plot,'ps
device,xsize=4.8,ysize=7.0,/inches
device,filename='cassini_nov.eps
!p.multi=[0,2,3]
!p.charsize=1.6
;device,/color
;loadct,5

;contour,narr(*,*,5)/flr,xarr,yarr,/xsty,$
;  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[2,5,6,7,8,9,10,12,20,25,30,35],$
;      title='S',c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1


;contour,narr(*,*,6)/flr,xarr,yarr,/xsty,$
;  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[10,20,30,40,50,80,120,160,180,220,260],$
;      title='O',c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1


nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)
nel = nel + 0.1*nel
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)


;;contour,nel,xarr,yarr,/xsty,$
;;C_LABELS=[1,1],levels=[1650,2000],c_color=[200,255],$
;;title='n!de!n',c_charsize=0.8,nlev=10,/fill   
;contour,nel/flr,xarr,yarr,/xsty,$
;C_LABELS=[1,1,0,1,0,1,0,1,0,1,0],levels=[500,1000,1500,2000,2500,3000,3500,4000,4500,5000,5500],$
;title='n!de!n',c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1


sz = size(nel)
wh = where(nel le 0.0)
whx = wh mod sz(1)
why = wh/sz(2)
nelec = nel
nelec(whx,why) = 0.00001
contour,narr(*,*,0)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1],$
   levels=[0.0603-0.00603,0.0603+0.00603],c_color=[200,255],$
;   levels=[0.05,0.07],c_color=[200,255],$
   title='S!u+!n/n!de!n',c_charsize=0.5,nlev=10,/cell_fill
contour,narr(*,*,0)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.02,0.04,0.06,0.08,0.1,0.12,0.15,0.20,0.30,0.40],$
   title='S!u+!n/n!de!n',c_charsize=0.8,nlev=10,/overplot 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

tau1 = 64.0
sn1 = 6.4
dtau = 2.0
dsn=0.6
plots,sn1-dsn,tau1-dtau,/data,thick=2.0
plots,sn1+dsn,tau1-dtau,/data,/continue,thick=2.0
plots,sn1+dsn,tau1+dtau,/data,/continue,thick=2.0
plots,sn1-dsn,tau1+dtau,/data,/continue,thick=2.0
plots,sn1-dsn,tau1-dtau,/data,/continue,thick=2.0

contour,narr(*,*,1)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],levels=[0.21-0.021,0.21+0.021],$
   title='S!u++!n/n!de!n',c_charsize=0.8,nlev=10,c_color=[200,200],/cell_fill   
contour,narr(*,*,1)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,0,1,1,1,1,1,1,1,1,1],levels=[0.1,0.14,0.18,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27],$
   title='S!u++!n/n!de!n',c_charsize=0.8,nlev=10,/overplot    
;plots,xx,yup,/data,linestyle=1
plots,xx(whdn),ydn(whdn),/data,linestyle=1

plots,sn1-dsn,tau1-dtau,/data,thick=2.0
plots,sn1+dsn,tau1-dtau,/data,/continue,thick=2.0
plots,sn1+dsn,tau1+dtau,/data,/continue,thick=2.0
plots,sn1-dsn,tau1+dtau,/data,/continue,thick=2.0
plots,sn1-dsn,tau1-dtau,/data,/continue,thick=2.0

contour,narr(*,*,2)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],levels=[0.037-0.0037,0.037+0.0037],$
   title='S!u+++!n/n!de!n',c_charsize=0.8,nlev=10,c_color=[200,255],/cell_fill      
contour,narr(*,*,2)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.001,0.005,0.01,0.015,0.02,0.030,0.04,0.05],$
   title='S!u+++!n/n!de!n',c_charsize=0.8,nlev=10,/overplot    
;plots,xx,yup,/data,linestyle=1
plots,xx(whdn),ydn(whdn),/data,linestyle=1


plots,sn1-dsn,tau1-dtau,/data,thick=2.0
plots,sn1+dsn,tau1-dtau,/data,/continue,thick=2.0
plots,sn1+dsn,tau1+dtau,/data,/continue,thick=2.0
plots,sn1-dsn,tau1+dtau,/data,/continue,thick=2.0
plots,sn1-dsn,tau1-dtau,/data,/continue,thick=2.0


contour,narr(*,*,3)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],levels=[0.27-0.027,0.27+0.027],$
   title='O!u+!n/n!de!n',c_charsize=0.8,nlev=10,c_color=[200,200],/cell_fill         
contour,narr(*,*,3)/nel,xarr,yarr,/xsty,$   
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],levels=[0.2,0.22,0.24,0.26,0.28,0.30,0.32],$
   title='O!u+!n/n!de!n',c_charsize=0.8,nlev=10,/overplot    
;plots,xx,yup,/data,linestyle=1
plots,xx(whdn),ydn(whdn),/data,linestyle=1

plots,sn1-dsn,tau1-dtau,/data,thick=2.0
plots,sn1+dsn,tau1-dtau,/data,/continue,thick=2.0
plots,sn1+dsn,tau1+dtau,/data,/continue,thick=2.0
plots,sn1-dsn,tau1+dtau,/data,/continue,thick=2.0
plots,sn1-dsn,tau1-dtau,/data,/continue,thick=2.0


contour,narr(*,*,4)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],levels=[0.017-0.0017,0.017+0.0017],$
   title='O!u++!n/n!de!n',c_charsize=0.8,nlev=10,c_color=[200,255],/cell_fill           
contour,narr(*,*,4)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1],$
   levels=[0.014,0.016,0.018,0.020,0.025,0.03,0.045,0.050],$
   title='O!u++!n/n!de!n',c_charsize=0.8,/overplot   
;plots,xx,yup,/data,linestyle=1
plots,xx(whdn),ydn(whdn),/data,linestyle=1

plots,sn1-dsn,tau1-dtau,/data,thick=2.0
plots,sn1+dsn,tau1-dtau,/data,/continue,thick=2.0
plots,sn1+dsn,tau1+dtau,/data,/continue,thick=2.0
plots,sn1-dsn,tau1+dtau,/data,/continue,thick=2.0
plots,sn1-dsn,tau1-dtau,/data,/continue,thick=2.0

sz = size(nel)
;wh = where(nel le 0.0)
;whx = wh mod sz(1)
;why = wh/sz(2)
narr1 = narr(*,*,1)
;narr1(whx,why) = 0.00001

;contour,(narr(*,*,0)/narr1),xarr,yarr,/xsty,$
;    C_LABELS=[1,1],levels=[0.413-0.0413,0.413+0.0413],c_color=[200,255],$
;    title='S!u+!n/S!u++!n',c_charsize=0.8,nlev=10,/fill 
;contour,(narr(*,*,0)/narr(*,*,1))/flr,xarr,yarr,/xsty,$
;    C_LABELS=[1,1,1,1,1,1,1,1],levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0],$
;    title='S!u+!n/S!u++!n',c_charsize=0.8,nlev=10,/overplot 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,(narr(*,*,2)/narr1),xarr,yarr,/xsty,$
;    C_LABELS=[1,1],levels=[0.098-0.0098,0.098+0.0098],c_color=[200,255],$
;    title='S!u+++!n/S!u++!n',c_charsize=0.8,nlev=10,/fill 
;contour,(narr(*,*,2)/narr(*,*,1))/flr,xarr,yarr,/xsty,$
;    C_LABELS=[1,1,1,1,1,1,1,1],levels=[0.05,0.1,0.15,0.2,0.25,0.3],$
;    title='S!u+++!n/S!u++!n',c_charsize=0.8,nlev=10,/overplot 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;sonp = narr(*,*,3)+narr(*,*,4)
;ssnp = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)

;contour,sonp/ssnp,xarr,yarr,/xsty,$
;C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1,1],levels=[0.8,0.85,0.9,0.92,0.94,0.96,0.98,1.0],$
;title='!9S!3O!un+!n/!9S!3S!un+!n',c_charsize=0.8,nlev=10   
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1
;print,sonp/ssnp


;contour,nel/sni,xarr,yarr,/xsty,$
;C_LABELS=[1,0,1,0,1,0,1,1,1,1],levels=[1.1,1.15,1.2,1.25,1.30,1.35,1.4,1.45,1.5,1.55,1.60],$
;title='n!de!n/!9S!3n!di!n',c_charsize=0.8,nlev=10   
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1



;contour,tarr(*,*,0)/flr,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1],levels=[4.0,4.4,4.8,5.2,5.6,6.0,6.4,6.8],$
;   title='T!de!n (eV)',c_charsize=0.8,nlev=10   
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

sti = tarr(*,*,1)*narr(*,*,0)+tarr(*,*,2)*narr(*,*,1)+$
      tarr(*,*,3)*narr(*,*,2)+tarr(*,*,4)*narr(*,*,3)+$
      tarr(*,*,5)*narr(*,*,4)
sni = narr(*,*,0)+narr(*,*,1)+narr(*,*,2)+narr(*,*,3)+narr(*,*,4)

;contour,(sti/sni)/flr,xarr,yarr,/xsty,$
;    C_LABELS=[1,1,1,1,1,1,1,1],levels=[40,80,120,160,200,240],$
;    title='T!di!n (eV)',c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

ion_cx = (ecarr(*,*,0)+ecarr(*,*,1)+ecarr(*,*,3)+ecarr(*,*,4))/(ecarr(*,*,2)+ecarr(*,*,5))

;contour,ion_cx/flr,xarr,yarr,/xsty,$
;    C_LABELS=[1,1,1,1,1,1,1,1,1,1],$
;    levels=[0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9],$
;  title='P!dion!n/P!dcx!n', c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,ecarr(*,*,12)/flr,xarr,yarr,/xsty,$
;  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7],$
;  title='P!dfast!n ', c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,(ecarr(*,*,13)+ecarr(*,*,14))/flr,xarr,yarr,/xsty,$
;  C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1],levels=[0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6],$
;  title='P!dtrans!n ', c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

;contour,ecarr(*,*,17)/flr,xarr,yarr,/xsty,$
;    C_LABELS=[1,1,1,1,0,1,0,1,1,1],levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8],$
;  title='P!dion->e!n ', c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1

vol=3.5e31
puv_obs = 1.5006927e+12/1.6e-19/vol
dpuv_obs = 4.2e11/1.6e-19/vol

contour,ecarr(*,*,18),xarr,yarr,/xsty,$
    C_LABELS=[1,1],levels=[puv_obs-dpuv_obs,puv_obs+dpuv_obs],c_color=[200,255],$
    title='P!dEUV!n',c_charsize=0.8,nlev=10,/cell_fill 
contour,ecarr(*,*,18)/flr,xarr,yarr,/xsty,$
  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[0.1,0.2,0.3,0.4,0.5,1.8,2.2],$
  title='P!dEUV!n ', c_charsize=0.8,nlev=10,/overplot
;plots,xx,yup,/data,linestyle=1
plots,xx(whdn),ydn(whdn),/data,linestyle=1

plots,sn1-dsn,tau1-dtau,/data,thick=2.0
plots,sn1+dsn,tau1-dtau,/data,/continue,thick=2.0
plots,sn1+dsn,tau1+dtau,/data,/continue,thick=2.0
plots,sn1-dsn,tau1+dtau,/data,/continue,thick=2.0
plots,sn1-dsn,tau1-dtau,/data,/continue,thick=2.0

;contour,ecarr(*,*,9)/flr,xarr,yarr,/xsty,$
;    C_LABELS=[1,1,0,1,1,0,1,1,1,1,1],$
;    levels=[0.01,0.1,0.2,0.4,0.6,0.8,1.0,2.0],$
;  title='P!deh!n ', c_charsize=0.8,nlev=10 
;;plots,xx,yup,/data,linestyle=1
;;plots,xx(whdn),ydn(whdn),/data,linestyle=1






xyouts,2500.*0.0,2500.*2.6,'Transport time (days)',$
  orientation=90.0,/device,charsize=1.2
xyouts,2500.*1.5,2500.*(-0.2),'Source rate (10!u-4!n cm!u-3!n s!u-1!n)',$
  /device,charsize=1.2
device,/close

device,filename='cassini_2plots.eps'
;!p.multi=[0,2,1]
contour,narr(*,*,0)/nelec,xarr,yarr,/xsty,$
   C_LABELS=[1,1],$
   levels=[0.05,0.07],c_color=[200,255],$
   title='S!u+!n/n!de!n',c_charsize=0.5,nlev=10,/fill
contour,narr(*,*,0)/nel,xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.02,0.04,0.06,0.08,0.1,0.12,0.15,0.20,0.30,0.40],$
   title='S!u+!n/n!de!n',c_charsize=0.5,nlev=10,/overplot 
plots,xx(whdn),ydn(whdn),/data,linestyle=1

contour,ecarr(*,*,11),xarr,yarr,/xsty,$
    C_LABELS=[1,1],levels=[0.4,0.8],c_color=[200,255],$
    title='P!duv!n',c_charsize=0.5,nlev=10,/fill 
contour,ecarr(*,*,11)/flr,xarr,yarr,/xsty,$
  C_LABELS=[1,1,1,1,1,1,1,1,1,1],levels=[0.2,0.4,0.6,1.0,1.4,1.8,2.2],$
  title='P!duv!n ', c_charsize=0.5,nlev=10,/overplot
;plots,xx,yup,/data,linestyle=1
plots,xx(whdn),ydn(whdn),/data,linestyle=1

device,/close


return
end
;-------------------------------------------------------------------





;-------------------------------------------------------------------
pro src_tau_chisqrd_oct
;-------------------------------------------------------------------
set_plot,'ps
@x6x9
device,/encapsulated
!x.margin=[4,2]
!y.margin=[2,2]
restore,'cm3_2d_out.sav
;restore,'cassini_teh40_feh003_otos1.7.sav
;restore,'cm3_2d_teh40_otos3_fehsrc.sav'

min_ns = 1e28
max_ns = 3e28
;!y.title='Transport time (days)'
;!x.title='source rate (10!u-4!n cm!u-3!n s!u-1!n)'
!p.charsize=1.0
!p.font=0
!x.range=[14,24]
!y.type = 1
!y.range=[12,26]
!y.style = 1
!x.style = 1
Rj = 7.14e4*1e5 ;cm
a = 5.5*Rj
b = 7.0*Rj
vol = 0.25*!pi^2*(a+b)*(b-a)^2
xarr = (1e4*xarr/vol)
;print,'xarr...',xarr
yarr = (1./yarr)/(3600.*24.)

;narr(*,*,0) = smooth(narr(*,*,0),2)
;narr(*,*,1) = smooth(narr(*,*,1),2)
;narr(*,*,2) = smooth(narr(*,*,2),2)
;narr(*,*,3) = smooth(narr(*,*,3),2)
;narr(*,*,4) = smooth(narr(*,*,4),2)

nel = narr(*,*,0)+2.0*narr(*,*,1)+3.0*narr(*,*,2)+narr(*,*,3)+ $
     2.0*narr(*,*,4)

nel = nel + 0.1*nel

;flr = fltarr(sz(1),sz(2))
;flr(whx,why) = 1.0

set_plot,'ps
device,xsize=6.2,ysize=7.2,/inches
;device,filename='cassini_teh40_feh003_otos1.7.eps
!p.multi=[0,3,3]

nelec = nel

nsp = narr(*,*,0)/nelec
chisp = (abs(nsp-0.0954456)/0.00963404)

sz = size(chisp)

contour,smooth(rebin(chisp,10.*sz(1),10.*sz(2)),4),rebin(xarr,10.*sz(1)),$
   rebin(yarr,10.*sz(2)),/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.1,0.2,0.5,1.0,2.0,5.0],$
   c_thick=[1,1,1,4,1,1],$
   title='S!u+!n/n!de!n',c_charsize=0.5

ns2p=narr(*,*,1)/nelec
chis2p = (abs(ns2p - 0.227762)/0.0218887)

contour,smooth(rebin(chis2p,10.*sz(1),10.*sz(2)),4),$
   rebin(xarr,10.*sz(1)),rebin(yarr,10.*sz(2)),/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.1,0.2,0.5,1.0,2.0,5.0],$
   c_thick=[1,1,1,4,1,1],$
   title='S!u++!n/n!de!n',c_charsize=0.5   

ns3p = narr(*,*,2)/nelec
chis3p = (abs(ns3p - 0.0223388)/0.00215144)

contour,smooth(rebin(chis3p,10.*sz(1),10.*sz(2)),4),rebin(xarr,10.*sz(1)),$
   rebin(yarr,10.*sz(2)),/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.1,0.2,0.5,1.0,2.0,5.0],$
   c_thick=[1,1,1,4,1,1],$
   title='S!u+++!n/n!de!n',c_charsize=0.5

sps2p = narr(*,*,0)/narr(*,*,1)
chisps2p = (abs(sps2p - 0.419)/0.0419)

contour,smooth(rebin(chisps2p,10.*sz(1),10.*sz(2)),4),rebin(xarr,10.*sz(1)),$
   rebin(yarr,10.*sz(2)),/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.1,0.2,0.5,1.0,2.0,5.0],$
   c_thick=[1,1,1,4,1,1],$
   title='S!u+!n/S!u++!n',c_charsize=0.5

s3ps2p = narr(*,*,2)/narr(*,*,1)
chis3ps2p = (abs(s3ps2p - 0.098408)/0.0098)

contour,smooth(rebin(smooth(chis3ps2p,2),10.*sz(1),10.*sz(2)),4),rebin(xarr,10.*sz(1)),$
   rebin(yarr,10.*sz(2)),/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.1,0.2,0.5,1.0,2.0,5.0],$
   c_thick=[1,1,1,4,1,1],$
   title='S!u+++!n/S!u++!n',c_charsize=0.5
    
nop = narr(*,*,3)/nelec
chiop = (abs(nop - 0.229816)/0.0229876)

contour,smooth(rebin(chiop,10.*sz(1),10.*sz(2)),4),rebin(xarr,10.*sz(1)),$
   rebin(yarr,10.*sz(2)),/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.1,0.2,0.5,1.0,2.0,5.0],$
   c_thick=[1,1,1,4,1,1],$
   title='O!u+!n/n!de!n',c_charsize=0.5

puv_obs = 2.0006927e+12/1.6e-19/1.5e31
dpuv_obs = 5.44e11/1.6e-19/1.5e31

puv = ecarr(*,*,11)
chipuv = (abs(puv - puv_obs)/dpuv_obs)     

contour,smooth(rebin(chipuv,10.*sz(1),10.*sz(2)),4),rebin(xarr,10.*sz(1)),$
   rebin(yarr,10.*sz(2)),/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[0.1,0.2,0.5,1.0,2.0,5.0],$
   c_thick=[1,1,1,4,1,1],$
    title='P!duv!n',c_charsize=0.5

chitot = (chisp+chis2p+chis3p+chiop+6*chipuv+chisps2p+chis3ps2p)/13.0

contour,smooth(rebin(chitot,10.*sz(1),10.*sz(2)),4),$
   rebin(xarr,10.*sz(1)),rebin(yarr,10.*sz(2)),/xsty,$
    C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[1.0,1.1,1.2,1.4,1.6,1.8,2.0,5.0,10.0,20.0,30.0],$
   c_thick=[1,1,1,1,1,4],$
    title='Chi tot',c_charsize=0.5


xyouts,2500.*0.0,2500.*3.5,'Transport time (days)',$
  orientation=90.0,/device,charsize=0.7
xyouts,2500.*2.5,2500.*(-0.2),'Source rate (10!u-4!n cm!u-3!n s!u-1!n)',$
  /device,charsize=0.7
device,/close

return
end
;-------------------------------------------------------------------



;-------------------------------------------------------------------
pro cm3_2d
;-------------------------------------------------------------------
@common

Te0 = 5.0
Ti0 = 70.0
Teh0 = 90.0
fh = 0.003
trans = 1.0/(100.0*8.64e4)
net_source = 2e28
otos = 2.0
s_source = net_source/(1.0+otos)
o_source = net_source*otos/(otos + 1.0)

runt = 50*8.64e4  ;seconds
dt0 = 10000.0      ;seconds

ns0 = 10.0      ;initial S neutral density (cm^-3)
no0 = 2.0*ns0     ;initial O neutral density
nsp0= 10.0        ;initial S+ density
ns2p0 = 10.0       ;initial S++ density
ns3p0 = 10.0      ;initial S+++ density
nop0 = 10.0       ;initial O+ density
no2p0 = 10.0      ;initial O++ density 
nel0 = nsp0+ns2p0+ns3p0+nop0+no2p0       ;initial electron density

k4=1.0
k10=1.0
k11=1.0
k13=1.0
k14=1.0
k15=1.0

rows=1
cols=1

plot_density = 1
plot_temperature = 0

lonoff =1.0
!p.multi=[0,1,1]

upsp = ns0
ups2p = ns2p0
ups3p = ns3p0
upop = nop0
upo2p = no2p0

cont='false'
ft_ave = 0

min_ns = 14e27
max_ns = 26e27
dns = 2e27

nns = round((max_ns-min_ns)/dns)

min_otos = 0.6
max_otos = 3.0
dotos = 0.2

notos = (max_otos - min_otos)/dotos

narr = fltarr(nns+1,notos+1,5)
parr = fltarr(nns+1,notos+1,6)
tarr = fltarr(nns+1,notos+1,6)
xarr = fltarr(nns+1)
yarr = fltarr(notos+1)

;ep = 1e-8
ep = 1e-1
ept = 1e-1

for i = 0,nns do begin
   for j = 0,notos do begin
      net_source = min_ns + i*dns
      otos = min_otos + j*dotos
      print,'net source, O to S ratio...',net_source,otos
      xarr(i) = net_source
      yarr(j) = otos
      JUMP:      
      cm3_model,n,r,T,src,lss,temps,dens
      sz = size(src) & z = sz(1)-1 & zm = sz(1)-2
      if (((dens(z).nsp - dens(zm).nsp) lt ep) and $
          ((dens(z).ns2p - dens(zm).ns2p) lt ep) and $
          ((dens(z).ns3p - dens(zm).ns3p) lt ep) and $
          ((dens(z).nop - dens(zm).nop) lt ep) and $
          ((dens(z).no2p - dens(zm).no2p) lt ep) and $
          ((dens(z).nel - dens(zm).nel) lt ep) and $
          ((temps(z).Tsp - temps(zm).Tsp) lt ept) and $
          ((temps(z).Ts2p - temps(zm).Ts2p) lt ept) and $
          ((temps(z).Ts3p - temps(zm).Ts3p) lt ept) and $
          ((temps(z).Top - temps(zm).Top) lt ept) and $
          ((temps(z).To2p - temps(zm).To2p) lt ept) and $
          ((temps(z).Telec - temps(zm).Telec) lt ept)) then begin
         cont = 'false'
        endif else begin
         cont = 'true'
         goto,JUMP 
	endelse
      print,n
      narr(i,j,0) = n.nsp
      narr(i,j,1) = n.ns2p
      narr(i,j,2) = n.ns3p
      narr(i,j,3) = n.nop
      narr(i,j,4) = n.no2p
      parr(i,j,0) = r.Puv
      parr(i,j,1) = r.psp
      parr(i,j,2) = r.ps2p
      parr(i,j,3) = r.ps3p
      parr(i,j,4) = r.pop
      parr(i,j,5) = r.po2p
      tarr(i,j,0) = T.Tel
      tarr(i,j,1) = T.Tsp
      tarr(i,j,2) = T.Ts2p
      tarr(i,j,3) = T.Ts3p
      tarr(i,j,4) = T.Top
      tarr(i,j,5) = T.To2p
;      cont = 'true'
   endfor
endfor


save,filename='cm3_2d_out.sav',narr,tarr,parr,xarr,yarr

make_cont

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro cm3_2d_fh_teh
;-------------------------------------------------------------------
@common

Te0 = 5.0
Ti0 = 70.0
Teh0 = 34.0
fh = 0.003
trans = 1.0/(80.0*8.64e4)
net_source = 2e28
otos = 2.0
s_source = net_source/(1.0+otos)
o_source = net_source*otos/(otos + 1.0)

runt = 50*8.64e4  ;seconds
dt0 = 5000.0      ;seconds

ns0 = 10.0      ;initial S neutral density (cm^-3)
no0 = 2.0*ns0     ;initial O neutral density
nsp0= 10.0        ;initial S+ density
ns2p0 = 10.0       ;initial S++ density
ns3p0 = 10.0      ;initial S+++ density
nop0 = 10.0       ;initial O+ density
no2p0 = 10.0      ;initial O++ density 
nel0 = nsp0+ns2p0+ns3p0+nop0+no2p0       ;initial electron density

k4=1.0
k10=1.0
k11=1.0
k13=1.0
k14=1.0
k15=1.0

rows=1
cols=1

plot_density = 1
plot_temperature = 0

lonoff =1.0
!p.multi=[0,1,1]

upsp = ns0
ups2p = ns2p0
ups3p = ns3p0
upop = nop0
upo2p = no2p0

cont='false'
ft_ave = 0

min_fh = 0.001
max_fh = 0.005
dfh = 0.001

nfh = round((max_fh-min_fh)/dfh)

min_teh = 40.0
max_teh = 100.0
dteh = 10.0

nteh = (max_teh - min_teh)/dteh

narr = fltarr(nfh+1,nteh+1,5)
parr = fltarr(nfh+1,nteh+1,6)
tarr = fltarr(nfh+1,nteh+1,6)
xarr = fltarr(nfh+1)
yarr = fltarr(nteh+1)

;ep = 1e-8
ep = 1e-1
ept = 1e-1

for i = 0,nfh do begin
   for j = 0,nteh do begin
;      net_source = min_ns + i*dns
;      otos = min_otos + j*dotos
      fh = min_fh + i*dfh
      teh0 = min_teh + j*dteh
      print,'fh, teh0...',fh, teh0
      xarr(i) = fh
      yarr(j) = teh0
      JUMP:      
      cm3_model,n,r,T,src,lss,temps,dens
      sz = size(src) & z = sz(1)-1 & zm = sz(1)-2
      if (((dens(z).nsp - dens(zm).nsp) lt ep) and $
          ((dens(z).ns2p - dens(zm).ns2p) lt ep) and $
          ((dens(z).ns3p - dens(zm).ns3p) lt ep) and $
          ((dens(z).nop - dens(zm).nop) lt ep) and $
          ((dens(z).no2p - dens(zm).no2p) lt ep) and $
          ((dens(z).nel - dens(zm).nel) lt ep) and $
          ((temps(z).Tsp - temps(zm).Tsp) lt ept) and $
          ((temps(z).Ts2p - temps(zm).Ts2p) lt ept) and $
          ((temps(z).Ts3p - temps(zm).Ts3p) lt ept) and $
          ((temps(z).Top - temps(zm).Top) lt ept) and $
          ((temps(z).To2p - temps(zm).To2p) lt ept) and $
          ((temps(z).Telec - temps(zm).Telec) lt ept)) then begin
         cont = 'false'
        endif else begin
         cont = 'true'
         goto,JUMP 
	endelse
      print,n
      narr(i,j,0) = n.nsp
      narr(i,j,1) = n.ns2p
      narr(i,j,2) = n.ns3p
      narr(i,j,3) = n.nop
      narr(i,j,4) = n.no2p
      parr(i,j,0) = r.Puv
      parr(i,j,1) = r.psp
      parr(i,j,2) = r.ps2p
      parr(i,j,3) = r.ps3p
      parr(i,j,4) = r.pop
      parr(i,j,5) = r.po2p
      tarr(i,j,0) = T.Tel
      tarr(i,j,1) = T.Tsp
      tarr(i,j,2) = T.Ts2p
      tarr(i,j,3) = T.Ts3p
      tarr(i,j,4) = T.Top
      tarr(i,j,5) = T.To2p
;      cont = 'true'
   endfor
endfor


save,filename='cm3_2d_out.sav',narr,tarr,parr,xarr,yarr

;make_cont

return
end
;-------------------------------------------------------------------



 
;-------------------------------------------------------------------
pro cm3_2d_otos_tau
;-------------------------------------------------------------------
@common

Te0 = 5.0
Ti0 = 70.0
Teh0 = 50.0
fh = 0.002
trans = 1.0/(100.0*8.64e4)
net_source = 2.0e28
otos = 2.2
s_source = net_source/(1.0+otos)
o_source = net_source*otos/(otos + 1.0)

runt = 50*8.64e4  ;seconds
dt0 = 10000.0      ;seconds

ns0 = 10.0      ;initial S neutral density (cm^-3)
no0 = 2.0*ns0     ;initial O neutral density
nsp0= 200.0        ;initial S+ density
ns2p0 = 400.0       ;initial S++ density
ns3p0 = 50.0      ;initial S+++ density
nop0 = 800.0       ;initial O+ density
no2p0 = 50.0      ;initial O++ density 
nel0 = nsp0+ns2p0+ns3p0+nop0+no2p0       ;initial electron density

k0=1.0
k1=1.0
k2=1.0
k3=1.0
k4=1.0
k5=1.0
k6=1.0
k7=1.0
k8=1.0
k9=1.0
k10=1.0
k11=1.0
k12=1.0
k13=1.0
k14=1.0
k15=1.0
k16=1.0

rows=1
cols=1

plot_density = 1
plot_temperature = 0

lonoff =1.0
!p.multi=[0,1,1]

upsp = ns0
ups2p = ns2p0
ups3p = ns3p0
upop = nop0
upo2p = no2p0

cont='false'
ft_ave = 0

min_ns = 16e27
max_ns = 44e27
dns = 4e27

nns = round((max_ns-min_ns)/dns)

min_otos = 1.0
max_otos = 4.0
dotos = 0.5

notos = (max_otos - min_otos)/dotos

min_tau = 20.0
max_tau = 80.0
dtau = 10

ntau = (max_tau - min_tau)/dtau

narr = fltarr(notos+1,ntau+1,7)
parr = fltarr(notos+1,ntau+1,6)
tarr = fltarr(notos+1,ntau+1,6)
xarr = fltarr(notos+1)
yarr = fltarr(ntau+1)
carr = fltarr(notos+1,ntau+1)  ;convergence
rarr = fltarr(notos+1,ntau+1,4)  ;energy input rates
nlarr = fltarr(notos+1,ntau+1,26) 
ecarr = fltarr(notos+1,ntau+1,18)

;ep = 1e-8
ep = 1e-1
ept = 1e-1

;for i = 0,nns do begin
for i = 0,notos do begin
   for j = 0,ntau do begin
;      net_source = min_ns + i*dns
      otos = min_otos + i*dotos
      tau = min_tau + j*dtau
      trans = 1.0/(tau*8.64e4)
      print,'otos, tau...',otos,tau
      xarr(i) = otos
;      yarr(j) = otos
      yarr(j) = trans
      JUMP:      
      cm3_model,n,r,T,src,lss,temps,dens,nl,Ec
      sz = size(src) & z = sz(1)-1 & zm = sz(1)-2
      carr(i,j) = (src(z).sp + src(z).s2p + src(z).s3p + src(z).op + $
                  src(z).o2p)/(lss(z).sp + lss(z).s2p + lss(z).s3p + $
                  lss(z).op + lss(z).o2p)
      rarr(i,j,0) = r.Sei
      rarr(i,j,1) = r.Secx
      rarr(i,j,2) = r.Oei
      rarr(i,j,3) = r.Oecx

      nlarr(i,j,0) = nl.sei
      nlarr(i,j,1) = nl.seih
      nlarr(i,j,2) = nl.oei
      nlarr(i,j,3) = nl.oeih
      nlarr(i,j,4) = nl.scx_tot
      nlarr(i,j,5) = nl.ocx_tot
      nlarr(i,j,6) = nl.sion_tot
      nlarr(i,j,7) = nl.oion_tot
      nlarr(i,j,8) = nl.k1
      nlarr(i,j,9) = nl.k2
      nlarr(i,j,10) = nl.k3
      nlarr(i,j,11) = nl.k4
      nlarr(i,j,12) = nl.k5
      nlarr(i,j,13) = nl.k6
      nlarr(i,j,14) = nl.k7
      nlarr(i,j,15) = nl.k8
      nlarr(i,j,16) = nl.k9
      nlarr(i,j,17) = nl.k10
      nlarr(i,j,18) = nl.k11
      nlarr(i,j,19) = nl.k12
      nlarr(i,j,20) = nl.k14
      nlarr(i,j,21) = nl.fast_S_k1
      nlarr(i,j,22) = nl.fast_S_k3
      nlarr(i,j,23) = nl.fast_O_k5
      nlarr(i,j,24) = nl.fast_O_k7
      nlarr(i,j,25) = nl.fast_O_k9

      ecarr(i,j,0) = Ec.s_ion
      ecarr(i,j,1) = Ec.s_ion_h
      ecarr(i,j,2) = Ec.s_cx
      ecarr(i,j,3) = Ec.o_ion
      ecarr(i,j,4) = Ec.o_ion_h
      ecarr(i,j,5) = Ec.o_cx
      ecarr(i,j,6) = Ec.s_tot_in
      ecarr(i,j,7) = Ec.o_tot_in
      ecarr(i,j,8) = Ec.P_pu
      ecarr(i,j,9) = Ec.eh_eq
      ecarr(i,j,10) = Ec.P_in
      ecarr(i,j,11) = Ec.Puv
      ecarr(i,j,12) = Ec.Pfast
      ecarr(i,j,13) = Ec.Ptrans
      ecarr(i,j,14) = Ec.Ptrans_eh
      ecarr(i,j,15) = Ec.P_out
      ecarr(i,j,16) = Ec.in_out
      ecarr(i,j,17) = Ec.ion_e_eq

      print,'convergence index...',carr(i,j)
      if (((dens(z).nsp - dens(zm).nsp) lt ep) and $
          ((dens(z).ns2p - dens(zm).ns2p) lt ep) and $
          ((dens(z).ns3p - dens(zm).ns3p) lt ep) and $
          ((dens(z).nop - dens(zm).nop) lt ep) and $
          ((dens(z).no2p - dens(zm).no2p) lt ep) and $
          ((dens(z).nel - dens(zm).nel) lt ep) and $
          ((temps(z).Tsp - temps(zm).Tsp) lt ept) and $
          ((temps(z).Ts2p - temps(zm).Ts2p) lt ept) and $
          ((temps(z).Ts3p - temps(zm).Ts3p) lt ept) and $
          ((temps(z).Top - temps(zm).Top) lt ept) and $
          ((temps(z).To2p - temps(zm).To2p) lt ept) and $
          ((temps(z).Telec - temps(zm).Telec) lt ept) and $
          (carr(i,j) lt 1.01) and (carr(i,j) gt 0.99)) then begin
         cont = 'false'
        endif else begin
         cont = 'true'
         goto,JUMP 
	endelse
      print,n
      carr(i,j) = (src(z).sp + src(z).s2p + src(z).s3p + src(z).op + $
                  src(z).o2p)/(lss(z).sp + lss(z).s2p + lss(z).s3p + $
                  lss(z).op + lss(z).o2p)
      print,'convergence index...',carr(i,j)
      narr(i,j,0) = n.nsp
      narr(i,j,1) = n.ns2p
      narr(i,j,2) = n.ns3p
      narr(i,j,3) = n.nop
      narr(i,j,4) = n.no2p
      narr(i,j,5) = n.ns
      narr(i,j,6) = n.no
      parr(i,j,0) = r.Puv
      parr(i,j,1) = r.psp
      parr(i,j,2) = r.ps2p
      parr(i,j,3) = r.ps3p
      parr(i,j,4) = r.pop
      parr(i,j,5) = r.po2p
      tarr(i,j,0) = T.Tel
      tarr(i,j,1) = T.Tsp
      tarr(i,j,2) = T.Ts2p
      tarr(i,j,3) = T.Ts3p
      tarr(i,j,4) = T.Top
      tarr(i,j,5) = T.To2p
;      cont = 'true'
   endfor
endfor


save,filename='cm3_2d_out.sav',narr,tarr,parr,xarr,yarr,carr,rarr,nlarr,ecarr

make_cont_otos_tau

return
end
;-------------------------------------------------------------------


 
;-------------------------------------------------------------------
pro cm3_2d_src_tau
;-------------------------------------------------------------------
@common

Te0 = 5.0
Ti0 = 70.0
Teh0 = 46.2
fh = 0.00247
trans = 1.0/(100.0*8.64e4)
net_source = 2.0e28
otos = 1.88
s_source = net_source/(1.0+otos)
o_source = net_source*otos/(otos + 1.0)

runt = 50*8.64e4  ;seconds
dt0 = 10000.0      ;seconds

ns0 = 10.0      ;initial S neutral density (cm^-3)
no0 = 2.0*ns0     ;initial O neutral density
nsp0= 200.0        ;initial S+ density
ns2p0 = 400.0       ;initial S++ density
ns3p0 = 50.0      ;initial S+++ density
nop0 = 800.0       ;initial O+ density
no2p0 = 50.0      ;initial O++ density 
nel0 = nsp0+ns2p0+ns3p0+nop0+no2p0       ;initial electron density

k0=1.0
k1=1.0
k2=1.0
k3=1.0
k4=1.0
k5=1.0
k6=1.0
k7=1.0
k8=1.0
k9=1.0
k10=1.0
k11=1.0
k12=1.0
k13=1.0
k14=1.0
k15=1.0
k16=1.0

rows=1
cols=1

plot_density = 1
plot_temperature = 0

lonoff =1.0
!p.multi=[0,1,1]

upsp = ns0
ups2p = ns2p0
ups3p = ns3p0
upop = nop0
upo2p = no2p0

cont='false'
ft_ave = 0

min_ns = 10e27
max_ns = 48e27
dns = 1e27

nns = round((max_ns-min_ns)/dns)

min_otos = 2.0
max_otos = 5.0
dotos = 0.5

notos = (max_otos - min_otos)/dotos

min_tau = 20.0
max_tau = 68.0
dtau = 1
;tau_arr = [6,8,10,12,15,20,30,40,50,60,70]
;ntau = n_elements(tau_arr)

ntau = (max_tau - min_tau)/dtau

y0up = 15
x0up = 140.0
aup = (90.0 - y0up)/(10.0-x0up)^2

;mup = (10.-80.)/(140.-10.)
;y0up = 80.0-mup*10.0

;extended exploration
;mdn = (6.0-15.0)/(50.0-10.0)
;y0dn = 15.0 - mdn*10.0

;cassini
mdn = (50.0-60.0)/(45.0-35.0)
y0dn = 50. - mdn*45.0

;nominal voyager
;mdn = (60.0-45.0)/(30.0-50.0)
;y0dn = 60.0 - mdn*30.0

mtau = (0.001-0.005)/(80.-6.0)
fh0 = 0.005 - mtau*6.0

narr = fltarr(nns+1,ntau+1,7)
parr = fltarr(nns+1,ntau+1,6)
tarr = fltarr(nns+1,ntau+1,6)
xarr = fltarr(nns+1)
yarr = fltarr(ntau+1)
carr = fltarr(nns+1,ntau+1)  ;convergence
rarr = fltarr(nns+1,ntau+1,4)  ;energy input rates
nlarr = fltarr(nns+1,ntau+1,26) 
ecarr = fltarr(nns+1,ntau+1,19)

;ep = 1e-8
ep = 1e-1
ept = 1e-1

for i = 0,nns do begin
;for i = 0,notos do begin
   for j = 0,ntau do begin
      net_source = min_ns + i*dns
;      fh = 0.00067 + 3.33e-5*net_source/1e27
;      otos = min_otos + i*dotos
      tau = min_tau + j*dtau
;;      fh = fh0 + mtau*tau
;;      tau = tau_arr(j)
;      if ((net_source ge 49e27) and (tau ge 59)) then begin
;	print,'Skipping...'
;	goto, NML
;      endif
;;      if ((net_source ge 90e27) and (tau gt 15)) then goto, NML
;;extended exploration
;      if (tau gt aup*(net_source/1e27 - x0up)^2 + y0up) or $
;         (tau lt mdn*net_source/1e27 + y0dn) then begin
;cassini and nominal voyager
      if (tau gt mdn*net_source/1e27 + y0dn) then begin
	 print,'Skipping to next tau...',tau
         narr(i,j,*) = 0.0
         carr(i,j) = 0.0
         rarr(i,j,*) = 0.0
         nlarr(i,j,*) = 0.0
         ecarr(i,j,*) = 0.0
         parr(i,j,*) = 0.0
         tarr(i,j,*) = 0.0
         goto, NML
      endif

      trans = 1.0/(tau*8.64e4)
      print,'src, tau, fh...',net_source,tau,fh
      xarr(i) = net_source
;      yarr(j) = otos
      yarr(j) = trans
      JUMP:      
      cm3_model,n,r,T,src,lss,temps,dens,nl,Ec
      sz = size(src) & z = sz(1)-1 & zm = sz(1)-2
      carr(i,j) = (src(z).sp + src(z).s2p + src(z).s3p + src(z).op + $
                  src(z).o2p)/(lss(z).sp + lss(z).s2p + lss(z).s3p + $
                  lss(z).op + lss(z).o2p)
      rarr(i,j,0) = r.Sei
      rarr(i,j,1) = r.Secx
      rarr(i,j,2) = r.Oei
      rarr(i,j,3) = r.Oecx

      nlarr(i,j,0) = nl.sei
      nlarr(i,j,1) = nl.seih
      nlarr(i,j,2) = nl.oei
      nlarr(i,j,3) = nl.oeih
      nlarr(i,j,4) = nl.scx_tot
      nlarr(i,j,5) = nl.ocx_tot
      nlarr(i,j,6) = nl.sion_tot
      nlarr(i,j,7) = nl.oion_tot
      nlarr(i,j,8) = nl.k1
      nlarr(i,j,9) = nl.k2
      nlarr(i,j,10) = nl.k3
      nlarr(i,j,11) = nl.k4
      nlarr(i,j,12) = nl.k5
      nlarr(i,j,13) = nl.k6
      nlarr(i,j,14) = nl.k7
      nlarr(i,j,15) = nl.k8
      nlarr(i,j,16) = nl.k9
      nlarr(i,j,17) = nl.k10
      nlarr(i,j,18) = nl.k11
      nlarr(i,j,19) = nl.k12
      nlarr(i,j,20) = nl.k14
      nlarr(i,j,21) = nl.fast_S_k1
      nlarr(i,j,22) = nl.fast_S_k3
      nlarr(i,j,23) = nl.fast_O_k5
      nlarr(i,j,24) = nl.fast_O_k7
      nlarr(i,j,25) = nl.fast_O_k9

      ecarr(i,j,0) = Ec.s_ion
      ecarr(i,j,1) = Ec.s_ion_h
      ecarr(i,j,2) = Ec.s_cx
      ecarr(i,j,3) = Ec.o_ion
      ecarr(i,j,4) = Ec.o_ion_h
      ecarr(i,j,5) = Ec.o_cx
      ecarr(i,j,6) = Ec.s_tot_in
      ecarr(i,j,7) = Ec.o_tot_in
      ecarr(i,j,8) = Ec.P_pu
      ecarr(i,j,9) = Ec.eh_eq
      ecarr(i,j,10) = Ec.P_in
      ecarr(i,j,11) = Ec.Puv
      ecarr(i,j,12) = Ec.Pfast
      ecarr(i,j,13) = Ec.Ptrans
      ecarr(i,j,14) = Ec.Ptrans_eh
      ecarr(i,j,15) = Ec.P_out
      ecarr(i,j,16) = Ec.in_out
      ecarr(i,j,17) = Ec.ion_e_eq
      ecarr(i,j,18) = Ec.Peuv

      print,'convergence index...',carr(i,j)
      if (((dens(z).nsp - dens(zm).nsp) lt ep) and $
          ((dens(z).ns2p - dens(zm).ns2p) lt ep) and $
          ((dens(z).ns3p - dens(zm).ns3p) lt ep) and $
          ((dens(z).nop - dens(zm).nop) lt ep) and $
          ((dens(z).no2p - dens(zm).no2p) lt ep) and $
          ((dens(z).nel - dens(zm).nel) lt ep) and $
          ((temps(z).Tsp - temps(zm).Tsp) lt ept) and $
          ((temps(z).Ts2p - temps(zm).Ts2p) lt ept) and $
          ((temps(z).Ts3p - temps(zm).Ts3p) lt ept) and $
          ((temps(z).Top - temps(zm).Top) lt ept) and $
          ((temps(z).To2p - temps(zm).To2p) lt ept) and $
          ((temps(z).Telec - temps(zm).Telec) lt ept) and $
          (carr(i,j) lt 1.01) and (carr(i,j) gt 0.999)) then begin
         cont = 'false'
        endif else begin
         cont = 'true'
         goto,JUMP 
	endelse
      print,n
      carr(i,j) = (src(z).sp + src(z).s2p + src(z).s3p + src(z).op + $
                  src(z).o2p)/(lss(z).sp + lss(z).s2p + lss(z).s3p + $
                  lss(z).op + lss(z).o2p)
      print,'convergence index...',carr(i,j)
      narr(i,j,0) = n.nsp
      narr(i,j,1) = n.ns2p
      narr(i,j,2) = n.ns3p
      narr(i,j,3) = n.nop
      narr(i,j,4) = n.no2p
      narr(i,j,5) = n.ns
      narr(i,j,6) = n.no
      parr(i,j,0) = r.Puv
      parr(i,j,1) = r.psp
      parr(i,j,2) = r.ps2p
      parr(i,j,3) = r.ps3p
      parr(i,j,4) = r.pop
      parr(i,j,5) = r.po2p
      tarr(i,j,0) = T.Tel
      tarr(i,j,1) = T.Tsp
      tarr(i,j,2) = T.Ts2p
      tarr(i,j,3) = T.Ts3p
      tarr(i,j,4) = T.Top
      tarr(i,j,5) = T.To2p
;      cont = 'true'

      ns0 = n.ns      ;initial S neutral density (cm^-3)
      no0 = n.no     ;initial O neutral density
      nsp0= n.nsp        ;initial S+ density
      ns2p0 = n.ns2p       ;initial S++ density
      ns3p0 = n.ns3p    ;initial S+++ density
      nop0 = n.nop     ;initial O+ density
      no2p0 = n.no2p      ;initial O++ density 
      nel0 = n.nel       ;initial electron density

   NML:
   
   endfor
endfor


save,filename='cm3_2d_out.sav',narr,tarr,parr,xarr,yarr,carr,rarr,nlarr,ecarr

;make_cont_src_tau

return
end
;-------------------------------------------------------------------



;-------------------------------------------------------------------
pro cm3_2d_ns_no
;-------------------------------------------------------------------
@common

Te0 = 5.0
Ti0 = 70.0
Teh0 = 40.0
fh = 0.0020
trans = 1.0/(100.0*8.64e4)
net_source = 2e28
otos = 1.6
s_source = net_source/(1.0+otos)
o_source = net_source*otos/(otos + 1.0)

runt = 50*8.64e4  ;seconds
dt0 = 10000.0      ;seconds

;ns0 = 10.0      ;initial S neutral density (cm^-3)
;no0 = 2.0*ns0     ;initial O neutral density
nsp0= 250.0        ;initial S+ density
ns2p0 = 400.0       ;initial S++ density
ns3p0 = 10.0      ;initial S+++ density
nop0 = 900.0       ;initial O+ density
no2p0 = 10.0      ;initial O++ density 
nel0 = nsp0+ns2p0+ns3p0+nop0+no2p0       ;initial electron density

k0=1.0
k1=1.0
k2=1.0
k3=1.0
k4=1.0
k5=1.0
k6=1.0
k7=1.0
k8=1.0
k9=1.0
k10=1.0
k11=1.0
k12=1.0
k13=1.0
k14=1.0
k15=1.0

rows=1
cols=1

plot_density = 1
plot_temperature = 0

lonoff =1.0
!p.multi=[0,1,1]

upsp = ns0
ups2p = ns2p0
ups3p = ns3p0
upop = nop0
upo2p = no2p0

cont='false'
ft_ave = 0

min_ns = 4.0
max_ns = 6.0
dns = 1.0

nns = round((max_ns-min_ns)/dns)

min_no = 40.0
max_no = 60.0

dno = 10.0

nno = round((max_no-min_no)/dno)

narr = fltarr(nns+1,nno+1,5)
parr = fltarr(nns+1,nno+1,6)
tarr = fltarr(nns+1,nno+1,6)
xarr = fltarr(nns+1)
yarr = fltarr(nno+1)
carr = fltarr(nns+1,nno+1)  ;convergence
rarr = fltarr(nns+1,nno+1,4)  ;energy input rates
nlarr = fltarr(nns+1,nno+1,26) 

;ep = 1e-8
ep = 1e-1
ept = 1e-1

for i = 0,nns do begin
;   for j = 0,notos do begin
   for j = 0,nno do begin
      ns0 = min_ns +i *dns      
      no0 = min_no + j*dno
      print,'ns, no....',ns0, no0
      xarr(i) = ns0
      yarr(j) = no0
      JUMP:      
      cm3_model,n,r,T,src,lss,temps,dens,nl,Ec
      sz = size(src) & z = sz(1)-1 & zm = sz(1)-2
      carr(i,j) = (src(z).sp + src(z).s2p + src(z).s3p + src(z).op + $
                  src(z).o2p)/(lss(z).sp + lss(z).s2p + lss(z).s3p + $
                  lss(z).op + lss(z).o2p)
      rarr(i,j,0) = r.Sei
      rarr(i,j,1) = r.Secx
      rarr(i,j,2) = r.Oei
      rarr(i,j,3) = r.Oecx

      nlarr(i,j,0) = nl.sei
      nlarr(i,j,1) = nl.seih
      nlarr(i,j,2) = nl.oei
      nlarr(i,j,3) = nl.oeih
      nlarr(i,j,4) = nl.scx_tot
      nlarr(i,j,5) = nl.ocx_tot
      nlarr(i,j,6) = nl.sion_tot
      nlarr(i,j,7) = nl.oion_tot
      nlarr(i,j,8) = nl.k1
      nlarr(i,j,9) = nl.k2
      nlarr(i,j,10) = nl.k3
      nlarr(i,j,11) = nl.k4
      nlarr(i,j,12) = nl.k5
      nlarr(i,j,13) = nl.k6
      nlarr(i,j,14) = nl.k7
      nlarr(i,j,15) = nl.k8
      nlarr(i,j,16) = nl.k9
      nlarr(i,j,17) = nl.k10
      nlarr(i,j,18) = nl.k11
      nlarr(i,j,19) = nl.k12
      nlarr(i,j,20) = nl.k14


      print,'convergence index...',carr(i,j)
      if (((dens(z).nsp - dens(zm).nsp) lt ep) and $
          ((dens(z).ns2p - dens(zm).ns2p) lt ep) and $
          ((dens(z).ns3p - dens(zm).ns3p) lt ep) and $
          ((dens(z).nop - dens(zm).nop) lt ep) and $
          ((dens(z).no2p - dens(zm).no2p) lt ep) and $
          ((dens(z).nel - dens(zm).nel) lt ep) and $
          ((temps(z).Tsp - temps(zm).Tsp) lt ept) and $
          ((temps(z).Ts2p - temps(zm).Ts2p) lt ept) and $
          ((temps(z).Ts3p - temps(zm).Ts3p) lt ept) and $
          ((temps(z).Top - temps(zm).Top) lt ept) and $
          ((temps(z).To2p - temps(zm).To2p) lt ept) and $
          ((temps(z).Telec - temps(zm).Telec) lt ept) and $
          (carr(i,j) lt 1.001) and (carr(i,j) gt 0.999)) then begin
         cont = 'false'
        endif else begin
         cont = 'true'
         goto,JUMP 
	endelse
      print,n
      carr(i,j) = (src(z).sp + src(z).s2p + src(z).s3p + src(z).op + $
                  src(z).o2p)/(lss(z).sp + lss(z).s2p + lss(z).s3p + $
                  lss(z).op + lss(z).o2p)
      print,'convergence index...',carr(i,j)
      narr(i,j,0) = n.nsp
      narr(i,j,1) = n.ns2p
      narr(i,j,2) = n.ns3p
      narr(i,j,3) = n.nop
      narr(i,j,4) = n.no2p
      parr(i,j,0) = r.Puv
      parr(i,j,1) = r.psp
      parr(i,j,2) = r.ps2p
      parr(i,j,3) = r.ps3p
      parr(i,j,4) = r.pop
      parr(i,j,5) = r.po2p
      tarr(i,j,0) = T.Tel
      tarr(i,j,1) = T.Tsp
      tarr(i,j,2) = T.Ts2p
      tarr(i,j,3) = T.Ts3p
      tarr(i,j,4) = T.Top
      tarr(i,j,5) = T.To2p
;      cont = 'true'
   endfor
endfor


save,filename='cm3_2d_out.sav',narr,tarr,parr,xarr,yarr,carr,rarr,nlarr

make_cont_src_tau

return
end
;-------------------------------------------------------------------
