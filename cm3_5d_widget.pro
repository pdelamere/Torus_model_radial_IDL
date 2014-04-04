;-------------------------------------------------------------------
pro read_np
;-------------------------------------------------------------------
@common_5d

close,1
ntypefile = ''
print,'ntype...',ntype
if (ntype eq 0) then begin
   ntypefile = 'nfall_'
endif 
if (ntype eq 1) then begin
   ntypefile = 'npall_'
endif
if (ntype eq 2) then begin
   ntypefile = 'ddj_'
endif
if (ntype eq 3) then begin
   ntypefile = 't1_'
endif
if (ntype eq 4) then begin
   ntypefile = 't2_'
endif
if (ntype eq 5) then begin
   ntypefile = 't3_'
endif
openr,1,ntypefile+'1.dat',/f77_unformatted

frame=0ll
nt=0ll
nout=0ll
nx=0ll
ny=0ll
nz=0ll
readu,1,nt
readu,1,nout
readu,1,nx
readu,1,nz

close,1

icld = dblarr(nx,nz,/nozero)
ncld = dblarr(nx,nz,/nozero) 

mfile = ((frm-1)/(nfrm/nfiles)) + 1
;print,mfile,frm
frmcount=1
frmcount=1+((frm-1) mod (nfrm/nfiles))
openr,1,ntypefile+strtrim(string(mfile),1)+'.dat',/f77_unformatted
readu,1,nt
readu,1,nout
readu,1,nx
readu,1,nz

cnt=1
while not(eof(1)) do begin
   readu,1,frame
   print, ntypefile+strtrim(string(mfile),1)+' image #.....',frame
   readu,1,icld
   if (cnt eq frmcount) then goto, BAIL
   cnt=cnt+1
endwhile

close,1

BAIL:  

tot_icld = icld
icld = icld(minx:maxx,minz:maxz)

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro read_file,flg
;-------------------------------------------------------------------
@common_5d

close,1

frame=0ll
nt=0ll
nout=0ll
nx=0ll
nz=0ll

mfile = ((frm-1)/(nfrm/nfiles)) + 1
print,mfile,frm
frmcount=1
frmcount=1+((frm-1) mod (nfrm/nfiles))
openr,1,file(0)+strtrim(string(mfile),1)+'.dat',/f77_unformatted
readu,1,nt
print,nt
readu,1,nout
readu,1,nx
readu,1,nz
print,nt,nout,nx,nz

temparr = dblarr(nx,nz,3,/nozero)
cnt=1
while not(eof(1)) do begin
   readu,1,frame
   print,file(0)+strtrim(string(mfile),1)+' image #.....',frame
   readu,1,temparr
   if (cnt eq frmcount) then goto, BAIL
   cnt=cnt+1
endwhile

close,1

BAIL: 

tot_arr=temparr
if (flg eq 0) then begin
    arr = temparr 
endif else begin
   arr = temparr(minx:maxx,minz:maxz,*)
endelse   

end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro make_surface
;-------------------------------------------------------------------
@common_5d
;!p.font=-1

!x.title='x'
!y.title='z'
surface,reform(arr(*,*,vcomp)),charsize=1.0

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro make_contour
;-------------------------------------------------------------------
@common_5d

f_read_coord_2d,'coord.dat',x,z,dzc,dzg,nx,nz
x = x(minx:maxx)
z = z(minz:maxz)

zbin = (max(z)-min(z))/n_elements(z)
xbin = (max(x)-min(x))/n_elements(x)

!x.title='x (km)'
!y.title='z (km)'
ar=reform(arr(*,*,vcomp))
regrid_xz_2d,ar,xx,zx
print,'array sizes...',size(arr)
print,size(ar)
print,size(xx)
print,size(zx)
contour,ar,xx,zx,/c_annotation, $
   charsize=1.0,nlevels=nlev, $
   xrange=[min(x),max(x)],xsty=1,yrange=[min(z),max(z)],ysty=1

return
end
;-------------------------------------------------------------------

;-------------------------------------------------------------------
PRO img_velovec,img,U,V,X,Y, Missing = Missing, Length = length, $ 
                Dots = dots, Color=color, _EXTRA = extra, $
                WINDOW_SCALE = window_scale, ASPECT = aspect, INTERP = interp
;-------------------------------------------------------------------
@common_5d

a = not(bytscl(img))
U = smooth(U,2)
V = smooth(V,2)

sz=size(a)
ab = bytscl(a)
a=rebin(a,4*sz(1),4*sz(2))
sz=size(a)
 
!x.range=[min(X),max(X)]
!y.range=[min(Y),max(Y)]

!x.title = 'x'
!y.title='z'

!y.margin = [4,9]
!x.margin = [6,9]

clrb = findgen(sz(2))*255/max(findgen(sz(2)))
for i=sz(1)-3,sz(1)-1 do begin
   a(i,*) = max(clrb)-clrb
endfor

on_error,2                      ;Return to caller if an error occurs
sz = size(a)			;Size of image
if sz(0) lt 2 then message, 'Parameter not 2D'

	;set window used by contour
contour,[[0,0],[1,1]],/nodata, xstyle=4, ystyle = 4

px = !x.window * !d.x_vsize	;Get size of window in device units
py = !y.window * !d.y_vsize
swx = px(1)-px(0)		;Size in x in device units
swy = py(1)-py(0)		;Size in Y
six = float(sz(1))		;Image sizes
siy = float(sz(2))
aspi = six / siy		;Image aspect ratio
aspw = swx / swy		;Window aspect ratio
f = aspi / aspw			;Ratio of aspect ratios

if (!d.flags and 1) ne 0 then begin	;Scalable pixels?
  if keyword_set(aspect) then begin	;Retain aspect ratio?
				;Adjust window size
	if f ge 1.0 then swy = swy / f else swx = swx * f
	endif

  tv,a*numclr/max(a),px(0),py(0),xsize = swx, ysize = swy, /device

endif else begin	;Not scalable pixels	
   if keyword_set(window_scale) then begin ;Scale window to image?
	tv,a*numclr/max(a),px(0),py(0)	;Output image
	swx = six		;Set window size from image
	swy = siy
    endif else begin		;Scale window
	if keyword_set(aspect) then begin
		if f ge 1.0 then swy = swy / f else swx = swx * f
		endif		;aspect
	tv,poly_2d(a*numclr/max(a),$	;Have to resample image
		[[0,0],[six/swx,0]], [[0,siy/swy],[0,0]],$
		keyword_set(interp),swx,swy), $
		px(0),py(0)
	endelse			;window_scale
  endelse			;scalable pixels

        on_error,2                      ;Return to caller if an error occurs
        s = size(u)
        t = size(v)
        if s(0) ne 2 then begin 
baduv:   message, 'U and V parameters must be 2D and same size.'
                endif
        if total(abs(s(0:2)-t(0:2))) ne 0 then goto,baduv
;
        if n_params(0) lt 4 then x = findgen(s(1)) else $
                if n_elements(x) ne s(1) then begin
badxy:                  message, 'X and Y arrays have incorrect size.'
                        endif
        if n_params(1) lt 5 then y = findgen(s(2)) else $
                if n_elements(y) ne s(2) then goto,badxy
;
        if n_elements(missing) le 0 then missing = 1.0e30
        if n_elements(length) le 0 then length = 1.0

        mag = sqrt(u^2+v^2)             ;magnitude.
                ;Subscripts of good elements
        nbad = 0                        ;# of missing points
        if n_elements(missing) gt 0 then begin
                good = where(mag lt missing) 
                if keyword_set(dots) then bad = where(mag ge missing, nbad)
        endif else begin
                good = lindgen(n_elements(mag))
        endelse

        ugood = u(good)
        vgood = v(good)
        x0 = min(x)                     ;get scaling
        x1 = max(x)
        y0 = min(y)
        y1 = max(y)
	x_step=float(x1-x0)/float(s(1))   ; Convert to float. Integer math
	y_step=float(y1-y0)/float(s(2))   ; could result in divide by 0

	maxmag=max([abs(max(abs(ugood/x_step))),abs(max(abs(vgood/y_step)))])
	sina = length * (ugood/maxmag)
	cosa = length * (vgood/maxmag)
;
        if n_elements(title) le 0 then title = ''
        ;--------------  plot to get axes  ---------------
        if n_elements(color) eq 0 then color = !p.color
        x_b0=x0-x_step
	x_b1=x1+x_step
	y_b0=y0-y_step
	y_b1=y1+y_step
        if n_elements(position) eq 0 then begin
          plot,[x_b0,x_b1],[y_b1,y_b0],/noerase,/nodata,/xst,yst=9, $
            color=color, _EXTRA = extra,$
            pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev
        endif else begin
          plot,[x_b0,x_b1],[y_b1,y_b0],/noerase,/nodata,/xst,yst=9, $
            color=color, _EXTRA = extra,$
            pos = [px(0),py(0), px(0)+swx,py(0)+swy],/dev
        endelse
     
     axis,yaxis=1,ystyle=1,$
;     yticks=2,$
;     ytickv=[min(y),(min(y)+max(y))/2,max(y)], $
;     ytickname=[string(min(a)),string((min(a)+max(a))/2),string(max(a))], $
;     ytitle='Density (10!u6!n ions/cm!u3!n)',$
     yrange = [min(img),max(img)]



;
        r = .3                          ;len of arrow head
        angle = 22.5 * !dtor            ;Angle of arrowhead
        st = r * sin(angle)             ;sin 22.5 degs * length of head
        ct = r * cos(angle)
;
        for i=0,n_elements(good)-1 do begin     ;Each point
                x0 = x(good(i) mod s(1))        ;get coords of start & end
                dx = sina(i)
                x1 = x0 + dx
                y0 = y(good(i) / s(1))
                dy = cosa(i)
                y1 = y0 + dy
		xd=x_step
		yd=y_step
                plots,[x0,x1,x1-(ct*dx/xd-st*dy/yd)*xd, $
			x1,x1-(ct*dx/xd+st*dy/yd)*xd], $
                      [y0,y1,y1-(ct*dy/yd+st*dx/xd)*yd, $
			y1,y1-(ct*dy/yd-st*dx/xd)*yd], $
;                      color=(ab(x0,y0) + vclr ) mod !d.n_colors 
;                       color = abs(vclr-ab(x0,y0))
                       color = 0
                endfor
        if nbad gt 0 then $             ;Dots for missing?
                oplot, x(bad mod s(1)), y(bad / s(1)), psym=3, color=vclr
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
PRO VELOVEC,U,V,X,Y, Missing = Missing, Length = length, Dots = dots,  $
        Color=color, _EXTRA = extra
;-------------------------------------------------------------------

U = smooth(U,2,/edge_truncate)
V = smooth(V,2,/edge_truncate)

        on_error,2                      ;Return to caller if an error occurs
        s = size(u)
        t = size(v)
        if s(0) ne 2 then begin 
baduv:   message, 'U and V parameters must be 2D and same size.'
                endif
        if total(abs(s(0:2)-t(0:2))) ne 0 then goto,baduv
;
        if n_params(0) lt 3 then x = findgen(s(1)) else $
                if n_elements(x) ne s(1) then begin
badxy:                  message, 'X and Y arrays have incorrect size.'
                        endif
        if n_params(1) lt 4 then y = findgen(s(2)) else $
                if n_elements(y) ne s(2) then goto,badxy
;
        if n_elements(missing) le 0 then missing = 1.0e30
        if n_elements(length) le 0 then length = 1.0

        mag = sqrt(u^2+v^2)             ;magnitude.
                ;Subscripts of good elements
        nbad = 0                        ;# of missing points
        if n_elements(missing) gt 0 then begin
                good = where(mag lt missing) 
                if keyword_set(dots) then bad = where(mag ge missing, nbad)
        endif else begin
                good = lindgen(n_elements(mag))
        endelse

        ugood = u(good)
        vgood = v(good)
        x0 = min(x)                     ;get scaling
        x1 = max(x)
        y0 = min(y)
        y1 = max(y)
	x_step=float(x1-x0)/float(s(1))   ; Convert to float. Integer math
	y_step=float(y1-y0)/float(s(2))   ; could result in divide by 0

	maxmag=max([abs(max(abs(ugood/x_step))),abs(max(abs(vgood/y_step)))])
	sina = length * (ugood/maxmag)
	cosa = length * (vgood/maxmag)
;
        if n_elements(title) le 0 then title = ''
        ;--------------  plot to get axes  ---------------
        if n_elements(color) eq 0 then color = !p.color
        x_b0=x0-x_step
	x_b1=x1+x_step
	y_b0=y0-y_step
	y_b1=y1+y_step
        if n_elements(position) eq 0 then begin
          plot,[x_b0,x_b1],[y_b1,y_b0],/nodata,/xst,/yst, $
            color=color, _EXTRA = extra
        endif else begin
          plot,[x_b0,x_b1],[y_b1,y_b0],/nodata,/xst,/yst, $
            color=color, _EXTRA = extra
        endelse
;
        r = .3                          ;len of arrow head
        angle = 22.5 * !dtor            ;Angle of arrowhead
        st = r * sin(angle)             ;sin 22.5 degs * length of head
        ct = r * cos(angle)
;
        for i=0,n_elements(good)-1 do begin     ;Each point
                x0 = x(good(i) mod s(1))        ;get coords of start & end
                dx = sina(i)
                x1 = x0 + dx
                y0 = y(good(i) / s(1))
                dy = cosa(i)
                y1 = y0 + dy
		xd=x_step
		yd=y_step
                plots,[x0,x1,x1-(ct*dx/xd-st*dy/yd)*xd, $
			x1,x1-(ct*dx/xd+st*dy/yd)*xd], $
                      [y0,y1,y1-(ct*dy/yd+st*dx/xd)*yd, $
			y1,y1-(ct*dy/yd-st*dx/xd)*yd], $
                      color=color
                endfor
        if nbad gt 0 then $             ;Dots for missing?
                oplot, x(bad mod s(1)), y(bad / s(1)), psym=3, color=color
end
;-------------------------------------------------------------------


;--------------------------------------------------------------------------
PRO make_vector
;--------------------------------------------------------------------------
@common_5d

vxz = reform(arr(*,*,0))
vzx = reform(arr(*,*,2))

if not((abs(max(vxz)) eq 0) and (abs(max(vzx)) eq 0)) then begin
   !x.title='x'
   !y.title='z'
   !p.title=!p.title + $
            ' (max = '+strtrim(string(max(sqrt(vxz^2 + vzx^2))),2)+')'
   velovec,vxz,vzx,length=1.0,charsize=1.0
   !x.title='x'
   !y.title='y'
   !p.title=''
endif

return
end
;------------------------------------------------------------------------


;--------------------------------------------------------------------------
PRO make_img_vector
;--------------------------------------------------------------------------
@common_5d

f_read_coord_2d,'coord.dat',x,z,dzc,dzg,nx,nz
x = x(minx:maxx)
z = z(minz:maxz)

erase
img = reform(icld(*,*))
vxz = reform(arr(*,*,0))
vzx = reform(arr(*,*,2))

if not((abs(max(vxz)) eq 0) and (abs(max(vzx)) eq 0)) then begin
   !x.title='x'
   !y.title='z'
   ptit = !p.title
   !p.title=''
   img_velovec,img,vxz,vzx,x,z,length=1.0,charsize=1.0,/aspect, $
     title=ptit+' (max = '+strtrim(string(max(sqrt(vxz^2 + vzx^2))),2)+')'
   !x.title='x'
   !y.title='y'
endif

return
end
;------------------------------------------------------------------------


;--------------------------------------------------------------------------
PRO make_profile
;--------------------------------------------------------------------------
@common_5d

erase
img = reform(icld(*,*))
sz = size(img)
img=rebin(img,10*sz(1),10*sz(2))
tvscl,img
profiles,img

return
end
;------------------------------------------------------------------------


;-------------------------------------------------------------------
pro make_btrace
;-------------------------------------------------------------------
@common_5d 

if (file ne 'b1all_') then begin
   file = 'b1all_'
   read_file,1
endif

f_read_coord_2d,'coord.dat',x,z,dzg,dzc,nnx,nnz

erase
img = reform(icld(*,*))
bx = reform(arr(*,*,0))
bz = reform(arr(*,*,1))
nxx = nx
nzz = nz
xx = x
zz = z
!x.title='x'
!y.title='z'

dx = max(xx)/nxx
dz = max(zz)/nzz

print,dx,dz

contour,bx,/nodata,xrange=[0,nxx*dx],yrange=[0,nzz*dz],xstyle=1,ystyle=1,$
        xtitle = 'x (km)', ytitle = 'z (km)'

x0 = 0.0

repeat begin

x0 = x0 + 0.5
x1 = x0
z1 = 0.0
flg = 0

while (flg ne 1) do begin
  i = round(nxx*x1/max(xx))
  k = round(nzz*z1/max(zz))

  if (i gt nxx-1) then goto, BAIL
  if (i lt 0) then goto, BAIL
  if (k gt nzz-1) then goto, BAIL
  if (k lt 0) then goto, BAIL

  bx1 = bx(i,k)
  bz1 = bz(i,k) + 2e-5*1.6e-19/2.3e-25
  z2 = z1+dz
  x2 = x1+dz*bx1/abs(bz1)
  if ((x2-x1) gt dx) then begin
	x2 = x1 + dx
	z2 = z1 + dz*abs(bz1)/bx1
	endif
  plots,[x1,x2],[z1,z2],/data
;  print,x1,x2,z1,z2
  x1 = x2
  z1 = z2

  if (x1 lt 0) then flg = 1
  if (x1 gt max(xx)) then flg = 1

  if (z1 lt 0) then flg = 1
  if (z1 gt max(zz)) then flg = 1



endwhile

BAIL:

endrep until(x1 gt max(xx))

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro get_movie_frames
;-------------------------------------------------------------------
@common_5d

close,1
files = file+strtrim(string(1),1)+'.dat'
openr,1,files,/f77_unformatted

nt=0ll
nout=0ll
nx=0ll
nz=0ll

readu,1,nt
readu,1,nout
readu,1,nx
readu,1,nz

close,1

anmt_arr = dblarr(nx,nz,nfrm)

temparr = dblarr(nx,nz,3,/nozero)

i=0
for m=1,nfiles do begin

   files = file+strtrim(string(m),1)+'.dat'
   openr,1,files,/f77_unformatted
   print,'Reading file.....',files

   readu,1,nt
   readu,1,nout
   readu,1,nx
   readu,1,nz

   while (not(eof(1))) do begin
      readu,1,frame
      readu,1,temparr
      ar = reform(temparr(*,*,vcomp))
      regrid_xz_2d,ar,xx,xz
      anmt_arr(*,*,i) = ar
      i = i + 1
   endwhile
   close,1

endfor

end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro make_animate
;-------------------------------------------------------------------
@common_5d

get_movie_frames
arr1=bytscl(anmt_arr)
arr2=congrid(arr1,100,200,nfrm)
xinteranimate,set=[100,200,nfrm]
for i=0,nfrm-1 do begin
   xinteranimate,image=arr2(*,*,i),frame=i
endfor
xinteranimate

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro get_graphics
;-------------------------------------------------------------------
@common_5d

case graph_type of

   'surf': begin
      make_surface
      end;surf

   'cont': begin
      make_contour
      end;cont

   'anim': begin
      make_animate
      end;anim

   'vect': begin
      make_vector
      end;vect

   'img_vect': begin
      make_img_vector
      end;img_vect

   'profile': begin
      make_profile
      end;profile

   'btrc': begin
      make_btrace
      end;btrc

endcase

return
end
;-------------------------------------------------------------------

;-------------------------------------------------------------------
pro plot_src_tau
;-------------------------------------------------------------------
@common_5d
!p.charsize=1.4
;params = 'feh '+strtrim(string(feharr(ifeh)),2) 
!p.multi=[0,4,5]

!x.title='S!dn!n (10!u-4!n cm!u-3!n s!u-1!n)'
!y.title='Transport time (days)'

!x.range=[5.5,7.5]
!y.range=[60,70]

xarr = nsarr*1e4
yarr = tauarr

nel = narr(*,*,ifeh,iteh,ios,0)+2.0*narr(*,*,ifeh,iteh,ios,1)+$
     3.0*narr(*,*,ifeh,iteh,ios,2)+narr(*,*,ifeh,iteh,ios,3)+ $
     2.0*narr(*,*,ifeh,iteh,ios,4)

nel = nel + 0.1*nel

sni = narr(*,*,ifeh,iteh,ios,0)+narr(*,*,ifeh,iteh,ios,1)+$
      narr(*,*,ifeh,iteh,ios,2)+narr(*,*,ifeh,iteh,ios,3)+$
      narr(*,*,ifeh,iteh,ios,4)

contour,narr(*,*,ifeh,iteh,ios,5),nsarr*1e4,tauarr,/xsty,$
;  C_LABELS=[1,1,1,1,1,0,1,1,1,1],levels=[5,10,15,20,25,30,35],$
      title='S',c_charsize=0.8,nlev=10

contour,narr(*,*,ifeh,iteh,ios,6),nsarr*1e4,tauarr,/xsty,$
;  C_LABELS=[1,1,1,0,1,0,1,1,1,1],levels=[60,80,120,160,180,220,260],$
      title='O',c_charsize=0.8,nlev=10

;contour,nel,xarr,yarr,/xsty,$
;   C_LABELS=[1,1],$
;   levels=[1800,2200],c_color=[200,0],$
;   title='n!de!n',c_charsize=1.0,nlev=10,/fill  
loadct,2
contour,nel,xarr,yarr,levels=[1650,2000],c_color=[200,200],$
        c_labels=[0,0],c_thick=[3,3],title='n!de!n',c_charsize=0.8
contour,nel,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1],$
;   levels=[500,1000,1400,1800,2200,2600,3000,3400,3800,4200],$
   title='n!de!n',c_charsize=0.8,nlev=10,/overplot  



contour,narr(*,*,ifeh,iteh,ios,0)/nel,xarr,yarr,/xsty,$
   C_LABELS=[0,0],c_thick=[3,3],$
;   levels=[0.12,0.13],c_color=[200,0],$ ;voyager
;   levels=[0.088-0.0088,0.088+0.0088],c_color=[200,200],$  ;cassini Oct 5
;   levels=[0.0484-0.00484,0.0484+0.00484],c_color=[200,200],$  ;cassini Nov 2
   levels=[0.05,0.07],c_color=[200,100],/cell_fill,$  ;cassini Jan
   title='S!u+!n/n!de!n',c_charsize=0.8,nlev=10
contour,narr(*,*,ifeh,iteh,ios,0)/nel,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[0.05,0.1,0.12,0.15,0.2,0.25],$
   title='S!u+!n/n!de!n',c_charsize=0.8,nlev=10,/overplot

contour,narr(*,*,ifeh,iteh,ios,1)/nel,xarr,yarr,/xsty,$
   C_LABELS=[0,0],c_thick=[3,3],$
;   levels=[0.15,0.185],c_color=[200,0],$  ;voyager
;   levels=[0.227-0.0227,0.227+0.0227],c_color=[200,200],$   ;cassini Oct
;   levels=[0.21-0.021,0.21+0.021],c_color=[200,200],$   ;cassini Nov
   levels=[0.20,0.23],c_color=[200,100],/cell_fill,$   ;cassini Jan
   title='S!u++!n/n!de!n',c_charsize=0.8,nlev=10   
contour,narr(*,*,ifeh,iteh,ios,1)/nel,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[0.12,0.16,0.18,0.185,0.19,0.2,0.23,0.24],$
   title='S!u++!n/n!de!n',c_charsize=0.8,nlev=10,/overplot   

contour,narr(*,*,ifeh,iteh,ios,2)/nel,xarr,yarr,/xsty,$
   C_LABELS=[0,0],c_thick=[3,3],$
;   levels=[0.02,0.06],c_color=[200,0],$ ;voyager
;   levels=[0.0278-0.00278,0.0278+0.00278],c_color=[200,200],$ ;cassini Oct
;   levels=[0.0455-0.00455,0.0455+0.00455],c_color=[200,200],$ ;cassini Nov
   levels=[0.034,0.040],c_color=[200,100],/cell_fill,$ ;cassini Jan
   title='S!u+++!n/n!de!n',c_charsize=0.8,nlev=10   
contour,narr(*,*,ifeh,iteh,ios,2)/nel,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[0.01,0.02,0.03,0.04,0.05,0.06],$
   title='S!u+++!n/n!de!n',c_charsize=0.8,nlev=10,/overplot   

contour,narr(*,*,ifeh,iteh,ios,3)/nel,xarr,yarr,/xsty,$   
   C_LABELS=[0,0],c_thick=[3,3],$
;   levels=[0.28,0.42],c_color=[200,0],$ ;voyager
;   levels=[0.27-0.027,0.27+0.027],c_color=[200,200],$ ;cassini Oct
;   levels=[0.286-0.0286,0.286+0.0286],c_color=[200,200],$ ;cassini Nov
   levels=[0.26,0.29],c_color=[200,100],/cell_fill,$ ;cassini Jan
   title='O!u+!n/n!de!n',c_charsize=0.8,nlev=10   
contour,narr(*,*,ifeh,iteh,ios,3)/nel,xarr,yarr,/xsty,$   
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[0.32,0.34,0.36,0.38,0.40,0.42,0.44,0.46,0.48,0.5],$
   title='O!u+!n/n!de!n',c_charsize=0.8,nlev=10,/overplot   

;contour,narr(*,*,ifeh,iteh,ios,4)/nel,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[0.007,0.04],c_color=[200,200],$
;   title='O!u++!n/n!de!n',c_charsize=0.8,nlev=10  
contour,narr(*,*,ifeh,iteh,ios,4)/nel,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[0.01,0.02,0.03,0.04,0.05,0.06],$
   levels=[0.016,0.019],c_color=[200,100],/cell_fill,$
   title='O!u++!n/n!de!n',c_charsize=0.8,nlev=10
contour,narr(*,*,ifeh,iteh,ios,4)/nel,xarr,yarr,/xsty,$   
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[0.32,0.34,0.36,0.38,0.40,0.42,0.44,0.46,0.48,0.5],$
   title='O!u++!n/n!de!n',c_charsize=0.8,nlev=10,/overplot   

contour,nel/sni,xarr,yarr,/xsty,$
;  C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1],$
;  levels=[1.2,1.25,1.30,1.36,1.38,1.42,1.46,1.50],$
  title='n!de!n/!9S!3n!di!n',c_charsize=0.8,nlev=10  

sonp = narr(*,*,ifeh,iteh,ios,3)+narr(*,*,ifeh,iteh,ios,4)
ssnp = narr(*,*,ifeh,iteh,ios,0)+narr(*,*,ifeh,iteh,ios,1)+$
       narr(*,*,ifeh,iteh,ios,2)

contour,sonp/ssnp,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[1.0,1.1,1.2,1.3,1.35,1.4,1.5],$
   title='!9S!3O!un+!n/!9S!3S!un+!n',c_charsize=0.8,nlev=10  

;sz = size(nel)
;wh = where(nel le 0.0)
;whx = wh mod sz(1)
;why = wh/sz(2)
;narr1 = narr(*,*,1)
;narr1(whx,why) = 1.0
contour,narr(*,*,ifeh,iteh,ios,0)/narr(*,*,ifeh,iteh,ios,1),xarr,yarr,/xsty,$
  C_LABELS=[0,0],c_thick=[3,3],$
;  levels=[0.67,0.87],c_color=[200,0],$ ;voyager
;  levels=[0.387-0.0387,0.387+0.0387],c_color=[200,200],$ ;cassini Oct
;  levels=[0.232-0.0232,0.232+0.0232],c_color=[200,200],$ ;cassini Nov
  levels=[0.18,0.30],c_color=[200,200],$ ;cassini Jan
  title='S!u+!n/S!u++!n',c_charsize=0.8,nlev=10  
contour,narr(*,*,ifeh,iteh,ios,0)/narr(*,*,ifeh,iteh,ios,1),xarr,yarr,/xsty,$
;  C_LABELS=[1,1,1,1,1,1,1,1,1,1,1],$
;  levels=[0.2,0.3,0.4,0.5,0.6,0.7,1.0,2.0],$
;   levels=[0.387-0.0387,0.387+0.0387],c_color=[200,200],$ ;cassini
;   levels=[0.232-0.0232,0.232+0.0232],c_color=[200,200],$ ;cassini
  title='S!u+!n/S!u++!n',c_charsize=0.8,nlev=10,/overplot  

;print,'narr1...',narr1(where(narr1 le 0.0))
;wh = where(narr1 le 0.0)
;narr1(wh) = 1.0
contour,narr(*,*,ifeh,iteh,ios,2)/narr(*,*,ifeh,iteh,ios,1),xarr,yarr,/xsty,$
  C_LABELS=[0,0],c_thick=[3,3],$
;  levels=[0.11,0.4],c_color=[200,0],$ ;voyager
;  levels=[0.122-0.0122,0.122+0.0122],c_color=[200,200],$ ;cassini Oct
;  levels=[0.218-0.0218,0.218+0.0218],c_color=[200,200],$ ;cassini Nov
  levels=[0.11,0.22],c_color=[200,200],$ ;cassini Jan
  title='S!u+++!n/S!u++!n',c_charsize=0.8,nlev=10  
contour,narr(*,*,ifeh,iteh,ios,2)/narr(*,*,ifeh,iteh,ios,1),xarr,yarr,/xsty,$
;  C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],$
;  levels=[0.02,0.04,0.06,0.08,0.1,0.20,0.3],$
  title='S!u+++!n/S!u++!n',c_charsize=0.8,nlev=10,/overplot  

contour,tarr(*,*,ifeh,iteh,ios,0),xarr,yarr,/xsty,$
   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1],$
   levels=[4.0,4.2,4.4,4.6,4.8,5.0,5.2,5.4,5.6,5.8,6.0],$
   title='T!de!n (eV)',c_charsize=0.8  

sti = tarr(*,*,ifeh,iteh,ios,1)*narr(*,*,ifeh,iteh,ios,0)+$
      tarr(*,*,ifeh,iteh,ios,2)*narr(*,*,ifeh,iteh,ios,1)+$
      tarr(*,*,ifeh,iteh,ios,3)*narr(*,*,ifeh,iteh,ios,2)+$
      tarr(*,*,ifeh,iteh,ios,4)*narr(*,*,ifeh,iteh,ios,3)+$
      tarr(*,*,ifeh,iteh,ios,5)*narr(*,*,ifeh,iteh,ios,4)
sni = narr(*,*,ifeh,iteh,ios,0)+narr(*,*,ifeh,iteh,ios,1)+$
      narr(*,*,ifeh,iteh,ios,2)+narr(*,*,ifeh,iteh,ios,3)+$
      narr(*,*,ifeh,iteh,ios,4)

;contour,sti/sni,xarr,yarr,/xsty,$
;  C_LABELS=[0,0],c_thick=[3,3],$
;  levels=[80.0,100.0],c_color=[200,0],$
;  title='T!di!n (eV)',c_charsize=0.5,nlev=10,/fill  
contour,sti/sni,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[40,60,80,100,120,140,160,180,200,220],$
    title='T!di!n (eV)',c_charsize=0.8,nlev=10

ion_cx = (ecarr(*,*,ifeh,iteh,ios,0)+ecarr(*,*,ifeh,iteh,ios,1)+$
          ecarr(*,*,ifeh,iteh,ios,3)+$
          ecarr(*,*,ifeh,iteh,ios,4))/(ecarr(*,*,ifeh,iteh,ios,2)+$
          ecarr(*,*,ifeh,iteh,ios,5))

contour,ion_cx,xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75],$
  title='P!dion!n/P!dcx!n', c_charsize=0.8,nlev=10

contour,ecarr(*,*,ifeh,iteh,ios,9),xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[0.010,0.05,0.1,0.15,0.2,0.25,0.35,0.45,0.55,0.70,0.90],$
  title='P!deh!n ', c_charsize=0.8,nlev=10

contour,ecarr(*,*,ifeh,iteh,ios,17),xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[0.010,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8],$
  title='P!dion->e!n ', c_charsize=0.8,nlev=10


contour,ecarr(*,*,ifeh,iteh,ios,11),xarr,yarr,/xsty,$
  C_LABELS=[0,0],c_thick=[3,3],$
;  levels=[0.36,0.63],c_color=[200,200],$ ;cassini Oct
;  levels=[0.362,0.57],c_color=[200,200],$ ;cassini Nov
  levels=[0.4,0.8],c_color=[200,100],/cell_fill,$ ;cassini Jan
  title='P!duv!n ',c_charsize=0.8,nlev=10  
contour,ecarr(*,*,ifeh,iteh,ios,11),xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1.0,1.2],$
  title='P!duv!n ', c_charsize=0.8,nlev=10,/overplot

contour,ecarr(*,*,ifeh,iteh,ios,12),xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[0.05,0.1,0.15,0.2,0.25,0.28,0.32,0.36,0.40,0.44,0.48,0.52],$
  title='P!dfast!n ', c_charsize=0.8,nlev=10

contour,(ecarr(*,*,ifeh,iteh,ios,13)+ecarr(*,*,ifeh,iteh,ios,14)),xarr,yarr,/xsty,$
;   C_LABELS=[1,1,1,1,1,1,1,1,1,1,1,1],$
;   levels=[0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18],$
  title='P!dtrans!n ', c_charsize=0.8,nlev=10


return
end
;-------------------------------------------------------------------



;-------------------------------------------------------------------
pro plot_temps
;-------------------------------------------------------------------
@common_5d

      restore,'temps.sav'
      tm = findgen(ntm)*dt/8.64e4
      plot_io,tm,temps.Tsp,/nodata,title='Temperature',$
         xtitle='time (days) ',ytitle='T (eV)',$
         yrange=[0.1,1000.0],/ysty,xrange=[0,dt*ntm]/8.64e4,/xsty
      oplot,tm,temps.Tsp,psym=3,color=10,thick=2.5
      oplot,tm,temps.Ts2p,psym=3,color=75,thick=2.5
      oplot,tm,temps.Ts3p,psym=3,color=50,thick=2.5
      oplot,tm,temps.Top,psym=3,color=150,thick=2.5
      oplot,tm,temps.To2p,psym=3,color=200,thick=2.5
      oplot,tm,smooth(temps.Telec,4),psym=3,color=100,thick=2.5

      if (lonoff eq 1.0) then begin
      legend,['S!u+!n','S!u++!n','S!u+++!n','O!u+!n','O!u++!n','T!de!n'],$
       linestyle=[0,0,0,0,0,0],colors=[10,75,50,150,200,100],/bottom,$
       /right
      endif

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro plot_dens
;-------------------------------------------------------------------
@common_5d

      restore,'dens.sav'
      tm = findgen(ntm)*dt/8.64e4
      plot_io,tm,dens.nel,/nodata,title='Density',$
         xtitle='time (days) ',ytitle='n (cm!u-3!n)',$
         yrange=[0.1,10000.0],/ysty,xrange=[0,dt*ntm]/8.64e4,/xsty
      oplot,tm,dens.nsp,psym=3,color=10,thick=2.5
      oplot,tm,dens.ns2p,psym=3,color=75,thick=2.5
      oplot,tm,dens.ns3p,psym=3,color=50,thick=2.5
      oplot,tm,dens.nop,psym=3,color=150,thick=2.5
      oplot,tm,dens.no2p,psym=3,color=200,thick=2.5
      oplot,tm,dens.nel,psym=3,color=100,thick=2.5

      if (lonoff eq 1.0) then begin
      legend,['S!u+!n','S!u++!n','S!u+++!n','O!u+!n','O!u++!n','n!de!n'],$
       linestyle=[0,0,0,0,0,0],colors=[10,75,50,150,200,100],/bottom,$
       /right
      endif

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro view_event,ev
;-------------------------------------------------------------------
@common_5d

widget_control, ev.id, get_uvalue = uvalue

case uvalue of
 
   'xsn': begin
	fxsn = 0
        fxt = 1
	fxfeh = 1
	fxteh = 1
	fxos = 1
        xflg = 0
	rsn = findgen(nsn)
	rtau = itau
	rfeh = ifeh
	rteh = iteh
	ros = ios
	widget_control,wSn,sensitive=fxsn*fysn
	widget_control,wt,sensitive=fxt*fyt
	widget_control,wfeh,sensitive=fxfeh*fyfeh
        widget_control,wteh,sensitive=fxteh*fyteh
	widget_control,wos,sensitive=fxos*fyos
	end;xsn

   'xt': begin
	fxsn = 1
        fxt = 0
	fxfeh = 1
	fxteh = 1
	fxos = 1
        xflg = 1
	rsn = isn
	rtau = itau
	rfeh = ifeh
	rteh = iteh
	ros = ios
	widget_control,wSn,sensitive=fxsn*fysn
	widget_control,wt,sensitive=fxt*fyt
	widget_control,wfeh,sensitive=fxfeh*fyfeh
        widget_control,wteh,sensitive=fxteh*fyteh
	widget_control,wos,sensitive=fxos*fyos
	end;xsn

   'xfeh': begin
	fxsn = 1
        fxt = 1
	fxfeh = 0
	fxteh = 1
	fxos = 1
	widget_control,wSn,sensitive=fxsn*fysn
	widget_control,wt,sensitive=fxt*fyt
	widget_control,wfeh,sensitive=fxfeh*fyfeh
        widget_control,wteh,sensitive=fxteh*fyteh
	widget_control,wos,sensitive=fxos*fyos
	end;xsn

   'xteh': begin
	fxsn = 1
        fxt = 1
	fxfeh = 1
	fxteh = 0
	fxos = 1
	widget_control,wSn,sensitive=fxsn*fysn
	widget_control,wt,sensitive=fxt*fyt
	widget_control,wfeh,sensitive=fxfeh*fyfeh
        widget_control,wteh,sensitive=fxteh*fyteh
	widget_control,wos,sensitive=fxos*fyos
	end;xsn

   'xos': begin
	fxsn = 1
        fxt = 1
	fxfeh = 1
	fxteh = 1
	fxos = 0
	widget_control,wSn,sensitive=fxsn*fysn
	widget_control,wt,sensitive=fxt*fyt
	widget_control,wfeh,sensitive=fxfeh*fyfeh
        widget_control,wteh,sensitive=fxteh*fyteh
	widget_control,wos,sensitive=fxos*fyos
	end;xsn

   'ysn': begin
	fysn = 0
        fyt = 1
	fyfeh = 1
	fyteh = 1
	fyos = 1
	widget_control,wSn,sensitive=fxsn*fysn
	widget_control,wt,sensitive=fxt*fyt
	widget_control,wfeh,sensitive=fxfeh*fyfeh
        widget_control,wteh,sensitive=fxteh*fyteh
	widget_control,wos,sensitive=fxos*fyos
	end;xsn

   'yt': begin
	fysn = 1
        fyt = 0
	fyfeh = 1
	fyteh = 1
	fyos = 1
	widget_control,wSn,sensitive=fxsn*fysn
	widget_control,wt,sensitive=fxt*fyt
	widget_control,wfeh,sensitive=fxfeh*fyfeh
        widget_control,wteh,sensitive=fxteh*fyteh
	widget_control,wos,sensitive=fxos*fyos
	end;xsn

   'yfeh': begin
	fysn = 1
        fyt = 1
	fyfeh = 0
	fyteh = 1
	fyos = 1
	widget_control,wSn,sensitive=fxsn*fysn
	widget_control,wt,sensitive=fxt*fyt
	widget_control,wfeh,sensitive=fxfeh*fyfeh
        widget_control,wteh,sensitive=fxteh*fyteh
	widget_control,wos,sensitive=fxos*fyos
	end;xsn

   'yteh': begin
	fysn = 1
        fyt = 1
	fyfeh = 1
	fyteh = 0
	fyos = 1
	widget_control,wSn,sensitive=fxsn*fysn
	widget_control,wt,sensitive=fxt*fyt
	widget_control,wfeh,sensitive=fxfeh*fyfeh
        widget_control,wteh,sensitive=fxteh*fyteh
	widget_control,wos,sensitive=fxos*fyos
	end;xsn

   'yos': begin
	fysn = 1
        fyt = 1
	fyfeh = 1
	fyteh = 1
	fyos = 0
	widget_control,wSn,sensitive=fxsn*fysn
	widget_control,wt,sensitive=fxt*fyt
	widget_control,wfeh,sensitive=fxfeh*fyfeh
        widget_control,wteh,sensitive=fxteh*fyteh
	widget_control,wos,sensitive=fxos*fyos
	end;xsn

   'SN': begin
      WIDGET_CONTROL, wSn, GET_VALUE = isn
      WIDGET_CONTROL, wSn_lbl, $
                      SET_VALUE = STRING(nsarr(isn)/1e-4, format='(f6.2)')+'(10^-4)'
      get_graphics
      end;SN

   'TAU': begin
      WIDGET_CONTROL, wt, GET_VALUE = itau
      WIDGET_CONTROL, wt_lbl, $
                      SET_VALUE = STRING(tauarr(itau))
      end;TAU

   'FEH': begin
      WIDGET_CONTROL, wfeh, GET_VALUE = ifeh
      WIDGET_CONTROL, wfeh_lbl, $
                      SET_VALUE = STRING(feharr(ifeh))
      plot_src_tau
      end;TAU

   'TEH': begin
      WIDGET_CONTROL, wteh, GET_VALUE = iteh
      WIDGET_CONTROL, wteh_lbl, $
                      SET_VALUE = STRING(teharr(iteh))
      plot_src_tau
      end;TAU

   'OTOS': begin
      WIDGET_CONTROL, wos, GET_VALUE = ios
      WIDGET_CONTROL, wos_lbl, $
                      SET_VALUE = STRING(osarr(ios))
      plot_src_tau
      end;TAU



   'X': begin
      vcomp=0
      if not(ev.select eq 0) then get_graphics
      end;X

   'Y': begin
      vcomp=1
      if not(ev.select eq 0) then get_graphics
      end;Y

   'Z': begin
      vcomp=2
      if not(ev.select eq 0) then get_graphics
      end;Z

   'FRAME': begin
      widget_control, ev.id, get_value = frame
      frm=frame
      read_file,1
      read_np
      get_graphics
      end;FRAME

   'SURFACE': begin
      widget_control,wclr,sensitive=0
      widget_control,wnlev,sensitive=0
      graph_type = 'surf'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'surf'
      end;SURFACE

   'CONTOUR': begin
      widget_control,wclr,sensitive=0
      widget_control,wnlev,sensitive=1
      graph_type = 'cont'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'cont'
      end;CONTOUR

   'VECTOR': begin
      widget_control,wclr,sensitive=0
      widget_control,wnlev,sensitive=0
      graph_type = 'vect'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'vect'
      end;VECTOR

   'IMG_VECTOR': begin
      widget_control,wclr,sensitive=1
      widget_control,wnlev,sensitive=0
      graph_type = 'img_vect'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'img_vect'
      end;IMG_VECTOR

   'PROFILES': begin
      widget_control,wclr,sensitive=1
      widget_control,wnlev,sensitive=0
      graph_type = 'profile'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'profile'
      end;PROFILES

   'BTRACE': begin
      widget_control,wclr,sensitive=0
      widget_control,wnlev,sensitive=0
      graph_type = 'btrc'
      if not(ev.select eq 0) then get_graphics
      graph_type_prev = 'btrc'
      end;BTRACE


   'ANIMATE': begin 
      graph_type = 'anim'
      get_graphics
      graph_type = graph_type_prev
      end;ANIMATE
  
;   'VIEWXY': partmovie_xy_m,nfrm/nfiles,nfiles
   'VIEWXZ': partmovie_xz_m,nfrm/nfiles,nfiles
;   'VIEWYZ': partmovie_yz_m,nfrm/nfiles,nfiles
   
   'CTABLES':  xloadct

   'VEC_COLOR': begin
      widget_control, ev.id, get_value = vclr
      graph_type = 'img_vect'
      get_graphics
      graph_type_prev = 'img_vect'
      end;VEC_COLOR

   'NLEV': begin
      widget_control, ev.id, get_value = nlev
      graph_type = 'cont'
      get_graphics
      graph_type_prev = 'cont'
      end;NLEV

   'B1': begin
      file = 'b1all_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;B1

   'UF': begin 
      file = 'ufall_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;UF

   'AJ': begin 
      file = 'ajall_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;AJ

   'UGRADU': begin 
      file = 'ugraduall_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;AJ
      
   'E' : begin
      file = 'Eall_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;E      

   'EF' : begin
      file = 'Efall_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;EF      

   'UP' : begin
      file = 'upall_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;UP      

   'UI' : begin
      file = 'uiall_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;UI      

   'UF2' : begin
      file = 'uf2all_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;UF2      

   'UFP2' : begin
      file = 'ufp2all_'
      read_file,1
      if not(ev.select eq 0) then get_graphics
      end;UFP2    
  
   'NP' : begin
      ntype = 1
      read_np
      if not(ev.select eq 0) then get_graphics
      end;NP

   'NF' : begin
      ntype = 0
      read_np
      if not(ev.select eq 0) then get_graphics
      end;NF

   'DDJ' : begin
      ntype = 2
      read_np
      if not(ev.select eq 0) then get_graphics
      end;DDJ

   'T1' : begin
      ntype = 3
      read_np
      if not(ev.select eq 0) then get_graphics
      end;T1

   'T2' : begin
      ntype = 4
      read_np
      if not(ev.select eq 0) then get_graphics
      end;T2

   'T3' : begin
      ntype = 5
      read_np
      if not(ev.select eq 0) then get_graphics
      end;T3

   'PATCH' : begin
      get_sub_arrays
      end;PATCH
   
   'CHOPPER' : begin
        COMMON VOLUME_DATA, A
	A = icld
	SLICER
      end;CHOPPER

   'PS' : begin
;      get_ps_options
      set_plot,'ps
      loadct,0
;      !p.multi=[0,1,2]
      !p.font=0
      !p.charsize=1.0
      device,filename='cassini_plot.ps'
      device,/inches,xsize=7.0,xoffset=0.75,ysize=9.0,yoffset=1.0
      device,/color
      plot_src_tau
;      plot_dens
;      plot_temps
      device,/close
      set_plot,'x
      !p.font=-1
;      !p.multi=[0,cols,rows]
      spawn,'gv cassini_plot.ps'
      loadct,0      
      end;PS

   'ETEMP': begin
      widget_control, ev.id, get_value = etemp
      Te0 = etemp
;      cm3_model
      end;ETEMP

   'ETEMPH': begin
      widget_control, ev.id, get_value = etemph
      Teh0 = etemph
;      cm3_model
      end;ETEMPH

   'EFH': begin
      WIDGET_CONTROL, wefh, GET_VALUE = frac
      frac = frac/10000.0
      fh = frac
;      print,fh
      WIDGET_CONTROL, wefh_lbl, SET_VALUE = STRING(frac*100, format='(f6.3)')
;      goto, set_frac
;      cm3_model
      end;EFH

   'TRANS': begin
      widget_control, ev.id, get_value = tran
      tran = 10^(tran/333.)
      WIDGET_CONTROL, wtrans_lbl, SET_VALUE = $
	     STRING(tran, format='(f7.2)')
      trans = 1.0/(tran*8.64e4)
;      cm3_model
      end;ETEMP

   'NSRC': begin
      widget_control, ev.id, get_value = src
      src = 10^(src/50.0)
      WIDGET_CONTROL, wnsrc_lbl, SET_VALUE = $
	     STRING(src, format='(f7.2)')
      net_source = src*1e27
      print,net_source
      end;SSRC

;   'OTOS': begin
;      WIDGET_CONTROL, wotos, GET_VALUE = rat
;      otos = rat/10.
;      WIDGET_CONTROL, wotos_lbl, $
;                      SET_VALUE = STRING(rat/10.0, format='(f4.2)')
;      end;OTOS

   'DONE': widget_control,/destroy, ev.top

   'RUN': begin
;      plot,sin(findgen(100))/10
      cm3_model
      widget_control,wcont,sensitive=1
      widget_control,wpltden,sensitive=1
      widget_control,wplttmp,sensitive=1
      widget_control,wpltdist,sensitive=1
      widget_control,wpltpwr,sensitive=1
      cont='false'
      end;RUN

   'CONTIN': begin
      cont='true'
      cm3_model
      end

   'RUNSRC': begin
      get_src_func
      end;RUNSRC

   'RUNLSS': begin
      get_lss_func
      end;RUNLSS

   'RUNTIME': begin
      WIDGET_CONTROL, ev.id, GET_VALUE = tm
      tm = 10^(tm/333.)
      WIDGET_CONTROL, wrt_lbl, SET_VALUE = $
	     STRING(tm, format='(f7.2)')
      runt = tm*8.64e4  ;seconds
      end;RUNTIME
   
   'IONDEN': begin
;      WIDGET_CONTROL, ev.id, GET_VALUE = iden
      get_change_ion_density
      widget_control,wcont,sensitive=0
      end;IONDEN

   'CHEM': begin
;      WIDGET_CONTROL, ev.id, GET_VALUE = iden
      get_change_chemistry
      end;CHEM
   
   'ROWS': begin
      widget_control,ev.id,get_value = r
      rows = r
      !p.multi=[0,cols,rows]
      end

   'COLS': begin
      widget_control,ev.id,get_value = c
      cols = c
      !p.multi=[0,cols,rows]
      end

   'LON': begin
      lonoff = 1.0
      end

   'LOFF': begin
      lonoff = 0.0
      end

   'DTUP':  begin
      dt0 = dt0*2.0
      print,'time step (s)...',dt0
      end

   'DTDOWN':  begin
      dt0 = dt0/2.0
      print,'time step (s)...',dt0
      end

   'PLOTDEN': begin
      plot_density = 1
      plot_temperature = 0
      plot_dens
      end

   'PLOTTEMP': begin
      plot_density = 0
      plot_temperature = 1
      plot_temps
      end

   'PLOTDIST': begin
      restore,'restart_cm3.sav
      cm3_ekappa,n,T
      restore,'torus_profile.sav'
      plot_oi,narr(0,*),max_theta-findgen(2.0*max_theta),$
         yrange=[-max_theta,max_theta],/ysty,$
         title='System III Long = '+strtrim(string(s3*!radeg),2),$
         ytitle = 'Jovigraphic latitude',xtitle='density (cm!u-3!n)',$
         xrange=[1,3.0*n0(0)],/xsty
      oplot,narr(1,*),max_theta-findgen(2.0*max_theta),linestyle=1
      oplot,narr(2,*),max_theta-findgen(2.0*max_theta),linestyle=2
      oplot,narr(3,*),max_theta-findgen(2.0*max_theta),linestyle=3
      oplot,narr(4,*),max_theta-findgen(2.0*max_theta),linestyle=4
      oplot,narr(5,*),max_theta-findgen(2.0*max_theta),linestyle=5
      legend,['e!u-!n','O!u+!n','O!u++!n','S!u+!n','S!u++!n','S!u+++!n'],$
         linestyle=[0,1,2,3,4,5],/right
      plot_oi,Tkappa,max_theta-findgen(2.0*max_theta),$
         yrange=[-max_theta,max_theta],/ysty,$
         title='Electron Temp (Kappa = '+strtrim(string(kappa),2)+')',$
         ytitle = 'Jovigraphic latitude',xtitle='Temperature (eV)',$
         xrange=[1,100],/xsty

      set_plot,'ps
      device,filename='ftdist_21.eps
      device,/encapsulated
      !p.font=0
      !p.multi=[0,2,1]
      !p.charsize=1.0
      !x.margin=[4,2]
      !y.margin=[4,2]
      @x6x9
      device,/inches,xsize=10.0,ysize=4.0
      plot_oi,narr(0,*),max_theta-findgen(2.0*max_theta),$
         yrange=[-20,20],/ysty,$
         title='System III Long = '+strtrim(string(s3*!radeg),2),$
         ytitle = 'Jovigraphic latitude',xtitle='density (cm!u-3!n)',$
         xrange=[1,10.0*n0(0)],/xsty
      oplot,narr(1,*),max_theta-findgen(2.0*max_theta),linestyle=1
      oplot,narr(2,*),max_theta-findgen(2.0*max_theta),linestyle=2
      oplot,narr(3,*),max_theta-findgen(2.0*max_theta),linestyle=3
      oplot,narr(4,*),max_theta-findgen(2.0*max_theta),linestyle=4
      oplot,narr(5,*),max_theta-findgen(2.0*max_theta),linestyle=5
      legend,['e!u-!n','O!u+!n','O!u++!n','S!u+!n','S!u++!n','S!u+++!n'],$
         linestyle=[0,1,2,3,4,5],/right,/bottom
      plot,Tkappa,max_theta-findgen(2.0*max_theta),$
         yrange=[-20,20],/ysty,$
         title='Electron Temp (Kappa = '+strtrim(string(kappa),2)+')',$
         ;ytitle = 'Jovigraphic latitude',$
         xtitle='Temperature (eV)',$
         xrange=[0,15],/xsty
      device,/close
      set_plot,'x
      end

   'RADPOWER': begin
      restore,'pwr.sav'
      tm = findgen(ntm)*dt/8.64e4
      plot,tm,p.Puv,xtitle='time (days) ',ytitle='Power Radiated (eV/s)',$
              xrange=[0,dt*ntm]/8.64e4,/xsty
      oplot,tm,p.psp,linestyle=1 
      oplot,tm,p.ps2p,linestyle=2 
      oplot,tm,p.ps3p,linestyle=3 
      oplot,tm,p.pop,linestyle=4 
      oplot,tm,p.po2p,linestyle=5 
      if (lonoff eq 1.0) then begin
      legend,['P!dUV!n','S!u+!n','S!u++!n','S!u+++!n','O!u+!n','O!u++!n'],$
       linestyle=[0,1,2,3,4,5],/top,/left
      endif
      end;RADPOWER

   'FTAVEON': begin
      ft_ave = 1
      end

   'FTAVEOFF': begin
      ft_ave = 0
      end

endcase

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro get_view_widget
;-------------------------------------------------------------------
@common_5d

base = widget_base(title='Torus chemistry model',/row)
lc = widget_base(base,/column,space=10)
;lct = widget_base(lc,/row,space=10)
;lctl = widget_base(lct,/frame,/column,space=10)
;lctr = widget_base(lct,/frame,/column,space=10)
;lcchem = widget_base(lc,/frame,/column,space=5)
;lcm = widget_base(lc,/column,space=5)
;lcleg = widget_base(lc,/column,space=5)
;lcpmulti = widget_base(lc,/frame,/column,space=5)
mc = widget_base(base,/frame,/column)
;mct = widget_base(mc,/row,space=10)
;mctl = widget_base(mct,/column,/frame,space=10)
;mctr = widget_base(mct,/row,space=10)
;rc = widget_base(base,/column,space=10)
;rcelec = widget_base(rc,/frame,/column,space=5)
;rcsrc = widget_base(rc,/frame,/column,space=5)
;rcanim = widget_base(rc,/frame,/column,space=5)
;rcct = widget_base(rc,/column,space=5)


;xmenu,['Sn','t','feh','Teh','O/S'],lctl,/column,/exclusive,$
;   title='x',space=50,uvalue=['xsn','xt','xfeh','xteh','xos']
;xmenu,['Sn','t','feh','Teh','O/S'],lctr,/column,/exclusive,$
;   title='y',space=50,uvalue=['ysn','yt','yfeh','yteh','yos']

;wSn_lbl = WIDGET_LABEL(lc, VALUE = STRING(nsarr(isn)/1e-4,format='(f6.2)')+'(10^-4)')
;wSn = widget_slider(lc, title = 'Sn',uvalue='SN',xsize=150,$
;                   maximum=nsn,minimum=0,/SUPPRESS_VALUE,/DRAG)
;widget_control,wSn,sensitive=fxsn*fysn
;wt_lbl = WIDGET_LABEL(lc, VALUE = STRING(tauarr(itau)))
;wt = widget_slider(lc, title = 't',uvalue='TAU',xsize=150,$
;                   maximum=ntau,minimum=0,/SUPPRESS_VALUE,/DRAG)
;widget_control,wt,sensitive=fxt*fyt
wfeh_lbl = WIDGET_LABEL(lc, VALUE = STRING(feharr(ifeh)))
wfeh = widget_slider(lc, title = 'feh',uvalue='FEH',xsize=150,$
                   maximum=nfeh,minimum=0,/SUPPRESS_VALUE)
;widget_control,wfeh,sensitive=fxfeh*fyfeh
wteh_lbl = WIDGET_LABEL(lc, VALUE = STRING(teharr(iteh)))
wteh = widget_slider(lc, title = 'Teh',uvalue='TEH',xsize=150,$
                   maximum=nteh,minimum=0,/SUPPRESS_VALUE)
;widget_control,wteh,sensitive=fxteh*fyteh
wos_lbl = WIDGET_LABEL(lc, VALUE = STRING(osarr(ios)))
wos = widget_slider(lc, title = 'O/S',uvalue='OTOS',xsize=150,$
                   maximum=nos,minimum=0,/SUPPRESS_VALUE)
;widget_control,wos,sensitive=fxos*fyos

;w1 = widget_button(lct,value = 'start',uvalue='RUN',xsize=100)

;wcont = widget_button(lct,value = 'continue',uvalue='CONTIN',xsize=100)
;widget_control,wcont,sensitive=0


;w1 = widget_button(lct,value = 'dt < 2x',uvalue='DTUP',xsize=50)
;w1 = widget_button(lct,value = 'dt > 2x',uvalue='DTDOWN',xsize=50)

;wrt_lbl = widget_label(lct,value = STRING(runt/8.64e4, $
;                          format='(f7.2)'))
;wrt = widget_slider(lct,title='run time (days)', uvalue='RUNTIME', $
;                    maximum = 1000, minimum = 333,xsize=170,$
;                    /SUPPRESS_VALUE,/DRAG)
;widget_control,wrt,set_value=alog10(runt/8.64e4)*333.

;w1 = widget_button(lcm,value = 'source function',uvalue='RUNSRC',$
;                   xsize=100)
;w1 = widget_button(lcm,value = 'loss function',uvalue='RUNLSS',$
;                   xsize=100)
;w1 = widget_button(lcm,value = 'Initial ion density',uvalue='IONDEN',$
;                   xsize=150)
;w1 = widget_button(lcm,value = 'Chemistry',uvalue='CHEM',$
;                   xsize=150)
;xmenu,['on','off'],lcleg,/row,/exclusive,$
;      space=1,uvalue=['LON','LOFF'],title='Legend'

;w1 = widget_label(lcpmulti,value = '!p.multi')
;w1 = widget_slider(lcpmulti, title = 'rows',uvalue='ROWS',xsize=100,$
;                   maximum=2,minimum=1)
;w1 = widget_slider(lcpmulti, title = 'columns',uvalue='COLS',xsize=100,$
;                   maximum=2,minimum=1)


;w1 = widget_label(rcelec,value = 'Electron Parameters')
;wetemp = widget_slider(rcelec,title='Electron temp (eV)', uvalue='ETEMP', $
;                    maximum = 30, minimum = 1,xsize=170)
;widget_control,wetemp,set_value=Te0
;wetemph = widget_slider(rcelec,title='Hot electron temp (eV)', $
;                        uvalue='ETEMPH', $
;                    maximum = 1000, minimum = 20,xsize=170)
;widget_control,wetemph,set_value=Teh0
;wefh_lbl = WIDGET_LABEL(rcelec, VALUE = STRING(fh*100))
;wefh = widget_slider(rcelec,title='Hot electron frac (%)', uvalue='EFH', $
;                    maximum = 100, minimum = 0,/SUPPRESS_VALUE,/DRAG,$
;                    xsize=170)
;widget_control,wefh,set_value=fh*10000.

;xmenu,['on','off'],rcelec,/row,/exclusive,$
;      space=1,uvalue=['FTAVEON','FTAVEOFF'],title='Flux Tube Ave'

;wtrans_lbl = widget_label(rc,value = STRING((1.0/trans)/8.64e4, $
;                          format='(f7.2)'))
;wtrans = widget_slider(rc,title='Transport time (days)', uvalue='TRANS', $
;                    maximum = 1000, minimum = 250,xsize=170,$
;                    /SUPPRESS_VALUE,/DRAG)
;widget_control,wtrans,set_value=alog10((1.0/trans)/8.64e4)*333.

;w1 = widget_label(rcsrc,value = 'Source Parameters')
;wnsrc_lbl = widget_label(rcsrc,value = STRING(net_source/1e27, $
;                          format='(f7.2)'))
;wnsrc = widget_slider(rcsrc,title='Total source (1e27 1/s)', uvalue='NSRC', $
;                    maximum = 120, minimum = 1,xsize=170,$
;                    /SUPPRESS_VALUE,/DRAG)
;widget_control,wnsrc,set_value=alog10(net_source/1e27)*50.0

;wotos_lbl = WIDGET_LABEL(rcsrc, VALUE = STRING(otos))
;wotos = widget_slider(rcsrc,title='O/S source ratio', uvalue='OTOS', $
;                    maximum = 80, minimum = 0,/SUPPRESS_VALUE,/DRAG, $
;                    xsize=170)
;widget_control,wotos,set_value=otos*10

;wpltden = widget_button(mct,value = ' Plot Density ',uvalue = 'PLOTDEN')
;widget_control,wpltden,sensitive=0
;wplttmp = widget_button(mct,value = ' Plot Temperature ',uvalue = 'PLOTTEMP')
;widget_control,wplttmp,sensitive=0
;wpltdist = widget_button(mct,value = ' Fluxtube Distribution ',$
;   uvalue = 'PLOTDIST')
;widget_control,wpltdist,sensitive=0
;wpltpwr = widget_button(mct,value = ' Radiated Power ',$
;   uvalue = 'RADPOWER')
;widget_control,wpltpwr,sensitive=0

draw = widget_draw(mc,xsize=650,ysize=800)

w1 = widget_button(mc,value = ' PostScript ',uvalue = 'PS')
w1 = widget_button(mc,value = 'Exit',uvalue='DONE')

widget_control,/realize,base,/hourglass

widget_control,get_value = window,draw
wset,window
;im=read_bmp('citepimage.bmp',/rgb)
;tv,im,/true

xmanager,'view',base

end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
PRO src_func_event,ev
;-------------------------------------------------------------------
@common_5d

widget_control, ev.id, get_uvalue = uvalue

case uvalue of

   'DONE': begin 
      widget_control,/destroy, ev.top
      widget_control,get_value = window,draw
      wset,window
      end;DONE

endcase

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro get_src_func
;-------------------------------------------------------------------
@common_5d

base = widget_base(title='Source Function',/column)

draw_src = widget_draw(base,xsize=550,ysize=400)

;w1 = widget_button(mc,value = ' PostScript ',uvalue = 'PS')
w1 = widget_button(base,value = 'Exit',uvalue='DONE')

widget_control,/realize,base,/hourglass

widget_control,get_value = window,draw_src
wset,window
restore,'src_func.sav'
!p.multi=[0,3,2]
!p.charsize=2.0
!x.range=[0,max(src.tm)]
!y.range=[0,max([src.sp,src.s2p,src.s3p,src.op,src.o2p]*1e3)]
plot,src.tm,src.sp*1e3,title='S!u+!n',/xsty,/ysty,color=10,$
     xtitle='time (days)',ytitle='rate (10!u-3!n ions/s)'
plot,src.tm,src.s2p*1e3,title='S!u++!n',/xsty,/ysty,color=75,$
     xtitle='time (days)',ytitle='rate (10!u-3!n ions/s)'
plot,src.tm,src.s3p*1e3,title='S!u+++!n',/xsty,/ysty,color=50,$
     xtitle='time (days)',ytitle='rate (10!u-3!n ions/s)'
plot,src.tm,src.op*1e3,title='O!u+!n',/xsty,/ysty,color=150,$
     xtitle='time (days)',ytitle='rate (10!u-3!n ions/s)'
plot,src.tm,src.o2p*1e3,title='O!u++!n',/xsty,/ysty,color=200,$
     xtitle='time (days)',ytitle='rate (10!u-3!n ions/s)'
!p.multi=[0,rows,cols]
!p.charsize=1.2

xmanager,'src_func',base

end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
PRO lss_func_event,ev
;-------------------------------------------------------------------
@common_5d

widget_control, ev.id, get_uvalue = uvalue

case uvalue of

   'DONE': begin 
      widget_control,/destroy, ev.top
      widget_control,get_value = window,draw
      wset,window
      end;DONE

endcase

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro get_lss_func
;-------------------------------------------------------------------
@common_5d

base = widget_base(title='Loss Function',/column)

draw_lss = widget_draw(base,xsize=550,ysize=400)

;w1 = widget_button(mc,value = ' PostScript ',uvalue = 'PS')
w1 = widget_button(base,value = 'Exit',uvalue='DONE')

widget_control,/realize,base,/hourglass

widget_control,get_value = window,draw_lss
wset,window
restore,'lss_func.sav'
!p.multi=[0,3,2]
!p.charsize=2.0
!x.range=[0,max(lss.tm)]
!y.range=[0,max([lss.sp,lss.s2p,lss.s3p,lss.op,lss.o2p]*1e3)]
plot,lss.tm,lss.sp*1e3,title='S!u+!n',/xsty,/ysty,color=10,$
     xtitle='time (days)',ytitle='rate (10!u-3!n ions/s)'
plot,lss.tm,lss.s2p*1e3,title='S!u++!n',/xsty,/ysty,color=75,$
     xtitle='time (days)',ytitle='rate (10!u-3!n ions/s)'
plot,lss.tm,lss.s3p*1e3,title='S!u+++!n',/xsty,/ysty,color=50,$
     xtitle='time (days)',ytitle='rate (10!u-3!n ions/s)'
plot,lss.tm,lss.op*1e3,title='O!u+!n',/xsty,/ysty,color=150,$
     xtitle='time (days)',ytitle='rate (10!u-3!n ions/s)'
plot,lss.tm,lss.o2p*1e3,title='O!u++!n',/xsty,/ysty,color=200,$
     xtitle='time (days)',ytitle='rate (10!u-3!n ions/s)'
!p.multi=[0,rows,cols]
!p.charsize=1.2

xmanager,'lss_func',base

end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
PRO ionden_event,ev
;-------------------------------------------------------------------
@common_5d

widget_control, ev.id, get_uvalue = uvalue
widget_control,ev.id,get_value = label
label=label(0)

case uvalue of

   'NSP': nsp0 = label      
   'NS2P': ns2p0 = label
   'NS3P': ns3p0 = label
   'NOP' : nop0 = label
   'NO2P': no2p0 = label

   'USEPREV': begin
      widget_control,wsp,set_value=strtrim(string(upsp),2)
      nsp0 = upsp
      widget_control,ws2p,set_value=strtrim(string(ups2p),2)
      ns2p0 = ups2p
      widget_control,ws3p,set_value=strtrim(string(ups3p),2)
      ns3p0 = ups3p
      widget_control,wop,set_value=strtrim(string(upop),2)
      nop0 = upop
      widget_control,wo2p,set_value=strtrim(string(upo2p),2)
      no2p0 = upo2p
      end

   'DONE': begin 
      widget_control,/destroy, ev.top
      end;DONE

endcase

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro get_change_ion_density
;-------------------------------------------------------------------
@common_5d

base = widget_base(title='Initial Ion Density',/column)

;w1 = widget_button(mc,value = ' PostScript ',uvalue = 'PS')

wsp = widget_label(base,value='S+ (cm^-3) ')
wsp = widget_text(base, /editable,/tracking_events,uvalue='NSP')
widget_control,wsp,set_value=strtrim(string(nsp0),2)
ws2p = widget_label(base,value='S++ (cm^-3)')
ws2p = widget_text(base, /editable,/tracking_events,uvalue='NS2P')
widget_control,ws2p,set_value=strtrim(string(ns2p0),2)
ws3p = widget_label(base,value='S+++ (cm^-3)')
ws3p = widget_text(base, /editable,/tracking_events,uvalue='NS3P')
widget_control,ws3p,set_value=strtrim(string(ns3p0),2)
wop = widget_label(base,value='O+ (cm^-3)')
wop = widget_text(base, /editable,/tracking_events,uvalue='NOP')
widget_control,wop,set_value=strtrim(string(nop0),2)
wo2p = widget_label(base,value='O++ (cm^-3)')
wo2p = widget_text(base, /editable,/tracking_events,uvalue='NO2P')
widget_control,wo2p,set_value=strtrim(string(no2p0),2)
w1 = widget_button(base,value = 'Use previous values',uvalue='USEPREV')
w1 = widget_button(base,value = 'Exit',uvalue='DONE')

widget_control,/realize,base,/hourglass

xmanager,'ionden',base

end
;-------------------------------------------------------------------

;-------------------------------------------------------------------
PRO chem_event,ev
;-------------------------------------------------------------------
@common_5d

widget_control, ev.id, get_uvalue = uvalue

case uvalue of
   'CX0ON': begin
      k0 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX0OFF': begin
      k0 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX1ON': begin
      k1 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX1OFF': begin
      k1 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX2ON': begin
      k2 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX2OFF': begin
      k2 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX3ON': begin
      k3 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX3OFF': begin
      k3 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX4ON': begin
      k4 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX4OFF': begin
      k4 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX5ON': begin
      k5 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX5OFF': begin
      k5 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX6ON': begin
      k6 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX6OFF': begin
      k6 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX7ON': begin
      k7 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX7OFF': begin
      k7 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX8ON': begin
      k8 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX8OFF': begin
      k8 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX9ON': begin
      k9 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX9OFF': begin
      k9 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX10ON': begin
      k10 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end
   'CX10OFF': begin
      k10 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX11ON': begin
      k11 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end
   'CX11OFF': begin
      k11 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX12ON': begin
      k12 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end
   'CX12OFF': begin
      k12 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX13ON': begin
      k13 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end
   'CX13OFF': begin
      k13 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX14ON': begin
      k14 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end
   'CX14OFF': begin
      k14 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'CX15ON': begin
      k15 = 1.0
      widget_control,get_value = window,draw
      wset,window
      end
   'CX15OFF': begin
      k15 = 0.0
      widget_control,get_value = window,draw
      wset,window
      end

   'DONE': begin 
      widget_control,/destroy, ev.top
      widget_control,get_value = window,draw
      wset,window
      end;DONE

endcase

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
pro get_change_chemistry
;-------------------------------------------------------------------
@common_5d

base = widget_base(title='Chemistry',/column)
top = widget_base(base,/row)
tlc = widget_base(top,/column)
tmc = widget_base(top,/column)
trc = widget_base(top,/column)

xmenu,['on','off'],tlc,/row,/exclusive,$
      space=1,uvalue=['CX0ON','CX0OFF'],title='k0'
xmenu,['on','off'],tlc,/row,/exclusive,$
      space=1,uvalue=['CX1ON','CX1OFF'],title='k1'
xmenu,['on','off'],tlc,/row,/exclusive,$
      space=1,uvalue=['CX2ON','CX2OFF'],title='k2'
xmenu,['on','off'],tlc,/row,/exclusive,$
      space=1,uvalue=['CX4ON','CX4OFF'],title='k4'
xmenu,['on','off'],tlc,/row,/exclusive,$
      space=1,uvalue=['CX5ON','CX5OFF'],title='k5'
xmenu,['on','off'],tlc,/row,/exclusive,$
      space=1,uvalue=['CX6ON','CX6OFF'],title='k6'
xmenu,['on','off'],tlc,/row,/exclusive,$
      space=1,uvalue=['CX7ON','CX7OFF'],title='k7'
xmenu,['on','off'],trc,/row,/exclusive,$
      space=1,uvalue=['CX8ON','CX8OFF'],title='k8'
xmenu,['on','off'],trc,/row,/exclusive,$
      space=1,uvalue=['CX9ON','CX9OFF'],title='k9'
xmenu,['on','off'],trc,/row,/exclusive,$
      space=1,uvalue=['CX10ON','CX10OFF'],title='k10'
xmenu,['on','off'],trc,/row,/exclusive,$
      space=1,uvalue=['CX11ON','CX11OFF'],title='k11'
xmenu,['on','off'],trc,/row,/exclusive,$
      space=1,uvalue=['CX12ON','CX12OFF'],title='k12'
xmenu,['on','off'],trc,/row,/exclusive,$
      space=1,uvalue=['CX13ON','CX13OFF'],title='k13'
xmenu,['on','off'],trc,/row,/exclusive,$
      space=1,uvalue=['CX14ON','CX14OFF'],title='k14'
xmenu,['on','off'],trc,/row,/exclusive,$
      space=1,uvalue=['CX15ON','CX15OFF'],title='k15'

draw_chem = widget_draw(tmc,xsize=300,ysize=500)

;w1 = widget_button(mc,value = ' PostScript ',uvalue = 'PS')
w1 = widget_button(base,value = 'Exit',uvalue='DONE')

widget_control,/realize,base,/hourglass

widget_control,get_value = window,draw_chem
wset,window

xyouts,0.1,0.95,'k0: S!u+!n + S!u++!n -> S!u++!n + S!u+!n (8.1e-9)',/normal
xyouts,0.1,0.90,'k1:  S + S!u+!n + S!u+!n + S (2.4e-8)',/normal
xyouts,0.1,0.85,'k2:  S + S!u++!n --> S!u+!n + S!u+!n (3e-10)',/normal
xyouts,0.1,0.80,'k3:  S + S!u++!n --> S!u++!n + S (7.8e-9)',/normal
xyouts,0.1,0.75,'k4:  S + S!u3+!n --> S!u+!n + S!u++!n  (1.32e-8)',/normal
xyouts,0.1,0.70,'k5:  O + O!u+!n --> O!u+!n + O (1.32e-8)',/normal
xyouts,0.1,0.65,'k6:  O + O!u++!n --> O!u+!n + O!u+!n (5.2e-10)',/normal
xyouts,0.1,0.60,'k7:  O + O!u++!n --> O!u++!n + O (5.4e-9)',/normal
xyouts,0.1,0.55,'k8:  O + S!u+!n --> O!u+!n + S (9e-11)',/normal
xyouts,0.1,0.50,'k9:  S + O!u+!n --> S!u+!n + O (3e-9)',/normal
xyouts,0.1,0.45,'k10: S + O!u++!n --> S!u+!n + O!u+!n  (2.34e-8)',/normal
xyouts,0.1,0.40,'k11: S + O!u++!n --> S!u++!n + O!u+!n + e!u-!n  (1.62e-8)',/normal
xyouts,0.1,0.35,'k12: O + S!u++!n --> O!u+!n + S!u+!n (2.3e-9)',/normal
xyouts,0.1,0.30,'k13:  O!u++!n + S!u+!n --> O!u+!n + S!u++!n  (1.4e-9)',/normal
xyouts,0.1,0.25,'k14:  O!u!n + S!u3+!n --> O!u+!n + S!u++!n  (1.92e-8)',/norma
xyouts,0.1,0.20,'k15:  O!u++!n + S!u++!n --> O!u+!n + S!u3+!n  (9e-9)',/normal

xmanager,'chem',base

end
;-------------------------------------------------------------------



;-------------------------------------------------------------------
PRO size_event,ev
;-------------------------------------------------------------------
@common_5d

widget_control, ev.id, get_uvalue = uvalue
widget_control,ev.id,get_value = minmax

case uvalue of

   'RESET': begin
      minx=0 & maxx=nx-1
;      miny=0 & maxy=ny-1
      minz=0 & maxz=nz-1
      widget_control,wmnx,/realize,set_slider_min=0
      widget_control,wmxx,/realize,set_slider_max=nx-1
;      widget_control,wmny,/realize,set_slider_min=0
;      widget_control,wmxy,/realize,set_slider_max=ny-1
      widget_control,wmnz,/realize,set_slider_min=0
      widget_control,wmxz,/realize,set_slider_max=nz-1
;      widget_control,wxy,/realize,set_slider_min=0,set_slider_max=nz-1
;      widget_control,wxz,/realize,set_slider_min=0,set_slider_max=ny-1
;      widget_control,wyz,/realize,set_slider_min=0,set_slider_max=nx-1
      icld=tot_icld
      arr=tot_arr
      get_graphics
      end;RESET

   'XMIN': begin
      minx = minmax
;      widget_control,wyz,/realize,set_slider_min=minx,set_slider_max=maxx
;      if (slcyz lt minx) then slcyz = 0
      if (minx ge nx-3) then begin
         widget_control,wmxx,/realize,set_slider_min=nx-2
         minx=nx-3
         icld=tot_icld(minx:maxx,minz:maxz)
         arr=tot_arr(minx:maxx,minz:maxz,*)
         get_graphics
      endif else begin
         widget_control,wmxx,/realize,set_slider_min=minx+1
         print,minx,maxx,minz,maxz
         icld=tot_icld(minx:maxx,minz:maxz)
         arr=tot_arr(minx:maxx,minz:maxz,*)
         get_graphics
      endelse
      end;XMIN
   'XMAX': begin 
      maxx = minmax
;      widget_control,wyz,/realize,set_slider_min=minx,set_slider_max=maxx
;      if (slcyz gt maxx-minx) then slcyz = maxx-minx-1
      widget_control,wmnx,set_slider_max=maxx-1
      icld=tot_icld(minx:maxx,minz:maxz)
      arr=tot_arr(minx:maxx,minz:maxz,*)
      get_graphics
      end;XMAX

;   'YMIN': begin
;      miny = minmax
;      widget_control,wxz,/realize,set_slider_min=miny,set_slider_max=maxy
;      if (slcxz lt miny) then slcxz = 0
;      if (miny ge ny-3) then begin
;         widget_control,wmxy,/realize,set_slider_min=ny-2
;         miny=ny-3
;         icld=tot_icld(minx:maxx,miny:maxy,minz:maxz)
;         arr=tot_arr(minx:maxx,miny:maxy,minz:maxz,*)
;         get_graphics
;      endif else begin
;         widget_control,wmxy,/realize,set_slider_min=miny+1
;         icld=tot_icld(minx:maxx,miny:maxy,minz:maxz)
;         arr=tot_arr(minx:maxx,miny:maxy,minz:maxz,*)
;         get_graphics
;      endelse
;      end;YMIN
;   'YMAX': begin
;      maxy = minmax
;      widget_control,wxz,/realize,set_slider_min=miny,set_slider_max=maxy
;      if (slcxz gt maxy-miny) then slcxz = maxy-miny-1
;      widget_control,wmny,set_slider_max=maxy-1
;      icld=tot_icld(minx:maxx,miny:maxy,minz:maxz)
;      arr=tot_arr(minx:maxx,miny:maxy,minz:maxz,*)
;      get_graphics
;      end;YMAX

   'ZMIN': begin
      minz = minmax
;      widget_control,wxy,/realize,set_slider_min=minz,set_slider_max=maxz
;      if (slcxy lt minz) then slcxy = 0
      if (minz ge nz-3) then begin
         widget_control,wmxz,/realize,set_slider_min=nz-2
         minz=nz-3
         icld=tot_icld(minx:maxx,minz:maxz)
         arr=tot_arr(minx:maxx,minz:maxz,*)
         get_graphics
      endif else begin
         widget_control,wmxz,/realize,set_slider_min=minz+1
         icld=tot_icld(minx:maxx,minz:maxz)
         arr=tot_arr(minx:maxx,minz:maxz,*)
         get_graphics
      endelse
      end;ZMIN
   'ZMAX': begin
      maxz = minmax
;      widget_control,wxy,/realize,set_slider_min=minz,set_slider_max=maxz
;      if (slcxy gt maxz-minz) then slcxy = maxz-minz-1
      widget_control,wmnz,set_slider_max=maxz-1
      icld=tot_icld(minx:maxx,minz:maxz)
      arr=tot_arr(minx:maxx,minz:maxz,*)
      get_graphics
      end;ZMAX

   'DONE': begin 
      widget_control,/destroy, ev.top
;      icld=tot_icld(minx:maxx,miny:maxy,minz:maxz)
;      arr=tot_arr(minx:maxx,miny:maxy,minz:maxz,*)
;      get_graphics
      end;DONE

endcase

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
PRO get_sub_arrays
;-------------------------------------------------------------------
@common_5d

wsizer = widget_base(title='Sizer',/column)

w9 = widget_button(wsizer,value = 'reset',uvalue='RESET')

wmnx = widget_slider(wsizer,title='min x',uvalue='XMIN', $
                         maximum=nx-1, minimum=0)
wmxx = widget_slider(wsizer,title='max x',uvalue='XMAX', $
                         maximum=nx-1, minimum=0)
;wmny = widget_slider(wsizer,title='min y',uvalue='YMIN', $
;                         maximum=ny-1, minimum=0)
;wmxy = widget_slider(wsizer,title='max y',uvalue='YMAX', $
;                         maximum=ny-1, minimum=0)
wmnz = widget_slider(wsizer,title='min z',uvalue='ZMIN', $
                         maximum=nz-1, minimum=0)
wmxz = widget_slider(wsizer,title='max z',uvalue='ZMAX', $
                         maximum=nz-1, minimum=0)
w8 = widget_button(wsizer,value = 'done',uvalue='DONE')

widget_control,/realize,wsizer,/hourglass
xmanager,'size',wsizer

return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
PRO post_event,ev
;-------------------------------------------------------------------
@common_5d

widget_control, ev.id, get_uvalue = uvalue
widget_control,ev.id,get_value = label
label=label(0)
name=uvalue
val=label

case uvalue of

;   'FILE': begin
;       psfile = label 
;       print,label
;       end;FILE
;   'XTIT': !x.title = label
;   'YTIT': !y.title = label
;   'ZTIT': !z.title = label
;   'PTIT': !p.title = label
;   'PORT': begin
;       prt = 1
;       lnd = 0
;       end;PORT 
;   'LAND': begin
;       lnd = 1
;       prt = 0
;       end;LAND
;   'YSIZE': begin
;       ysz = label
;       end;YSIZE
;   'ENCAP': begin
;       encp = 1
;       noencp = 0
;       end;ENCAP
;   'NONENCAP': begin
;       noencp = 1
;       encp = 0
;       end;NONENCAP
;   'COLOR': begin
;       clr = 1
;       nclr = 0
;       end;COLOR
;   'NOCOLOR': begin
;       clr = 0
;       nclr = 1
;       end;COLOR
   'GO': begin 
      numclr = !d.n_colors-1
      set_plot,'ps
      !p.font=0
      !p.charsize=1.0
      device,filename=cassini_plot.ps,/color
;      if (prt eq 1) then device,/portrait
;      if (lnd eq 1) then device,/landscape
;      if (encp eq 1) then device,encapsulated = 1
;      if (noencp eq 1) then device,encapsulated = 0
      device,/inches,xoffset = 1.0,xsize=6.0,yoffset=1.0,ysize=9.0
;      if (clr eq 1) then device,color=1
;      if (nclr eq 1) then device,color=0
      cm3_model
      device,/close
      set_plot,'x'
      !p.font=-1
      !p.title=''
      !x.title='x'
      !y.title='y'
      !z.title='z'

      widget_control,/destroy, ev.top
      return
      end;GO

endcase

WIDGET_CONTROL, ev.top, get_uvalue=out_text
WIDGET_CONTROL, out_text , set_value=name + ': ' + string(val),/append


return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
PRO get_ps_options
;-------------------------------------------------------------------
@common_5d

wps = widget_base(title='PostScript',/row,/frame)
lc = widget_base(wps,/column,/frame,space=10)
rc = widget_base(wps,/column,space=10)
;w1 = widget_label(lc,value='file: ')
;w1 = widget_text(lc, /editable,uvalue='FILE')
;w1 = widget_label(lc,value='plot title: ')
;w1 = widget_text(lc, /editable,uvalue='PTIT')
;w1 = widget_label(lc,value='x title: ')
;w1 = widget_text(lc, /editable,uvalue='XTIT')
;w1 = widget_label(lc,value='y title: ')
;w1 = widget_text(lc, /editable,uvalue='YTIT')
;w1 = widget_label(lc,value='z title: ')
;w1 = widget_text(lc, /editable,uvalue='ZTIT')
;t1 = WIDGET_TEXT(rc, xsize=30, ysize=5, /SCROLL, $
;      value=[ 'Postscript editor'])
;     WIDGET_CONTROL, wps, set_uvalue=t1

;xmenu,['encapsulated','non-encasulated'],rc,/frame,/column,/exclusive,$
;      space=20,uvalue=['ENCAP','NONENCAP']
;xmenu,['portrait','landscape'],rc,/frame,/column,/exclusive,$
;      space=20,uvalue=['PORT','LAND']
;wysz = widget_slider(rc,title='ysize', uvalue='YSIZE', $
;                       maximum = 9, minimum=2)
;xmenu,['color','no color'],rc,/frame,/column,/exclusive,$
;      space=20,uvalue=['COLOR','NOCOLOR']
w8 = widget_button(rc,value = 'GO',uvalue='GO')


widget_control,/realize,wps,/hourglass
xmanager,'post',wps


return
end
;-------------------------------------------------------------------


;-------------------------------------------------------------------
PRO cm3
;-------------------------------------------------------------------
@common_5d


;restore,'cm3_5d_out_110402.sav
;restore,'cm3_5d_out_voyager.sav
restore,'cm3_5d_out.sav


nsarr = nsarr/2.5e31

sz = size(nsarr)
nsn = sz(1)-1
sz = size(tauarr)
ntau = sz(1)-1
sz = size(feharr)
nfeh = sz(1)-1
sz = size(teharr)
nteh = sz(1)-1
sz = size(osarr)
nos = sz(1)-1

fxsn = 0
fxt = 0
fxfeh = 0
fxteh = 0
fxos = 0

fysn = 0
fyt = 0
fyfeh = 0
fyteh = 0
fyos = 0
 
isn = 0
itau = 0
ifeh = 0
iteh = 0
ios = 0



;nfiles=nfil
;nfrm = nfiles*nfr

Te0 = 5.0
Ti0 = 100.0
Teh0 = 40.0
fh = 0.0015
trans = 1.0/(40.0*8.64e4)
net_source = 3e28
otos = 3.0
s_source = net_source/(1.0+otos)
o_source = net_source*otos/(otos + 1.0)

runt = 365*8.64e4  ;seconds
dt0 = 10000.0      ;seconds

ns0 = 10.0      ;initial S neutral density (cm^-3)
no0 = 50.0     ;initial O neutral density
;no0 = 50.0     ;initial O neutral density
;ns0 = 5.0      ;initial S neutral density (cm^-3)
nsp0= 250.0        ;initial S+ density
ns2p0 = 400.0       ;initial S++ density
ns3p0 = 50.0      ;initial S+++ density
nop0 = 700.0       ;initial O+ density
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

;vclr=0
;vcomp=0
;;pln='xy'
;;slcxy=0
;;slcxz=0
;;slcyz=0
;frm=1
;nlev=6
;graph_type = 'surf'
;graph_type_prev = 'surf'
;ntype = 1
;file = 'b1all_'
psfile = 'idl.ps'
prt = 1 & lnd = 0
encp = 0 & noencp = 1
clr = 0 & nclr = 1
numclr=!d.n_colors
;!x.title='x'
;;!y.title='y'
;!z.title='z'
;ysz = 5

;read_file,0
;minx = 0 & maxx = nx-1
;;miny = 0 & maxy = ny-1
;minz = 0 & maxz = nz-1

;read_np

get_view_widget

end
;-------------------------------------------------------------------







