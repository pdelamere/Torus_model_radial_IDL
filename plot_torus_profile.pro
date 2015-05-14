dL = 0.25
r = 6.0+findgen(34)*dL
L = strmid(strtrim(string(r),2),0,5) 
Rj = 7.14e4
max_theta=30

narr_sp = fltarr(n_elements(L),max_theta*2.0+1)
narr_s2p = fltarr(n_elements(L),max_theta*2.0+1)
narr_s3p = fltarr(n_elements(L),max_theta*2.0+1)
narr_op = fltarr(n_elements(L),max_theta*2.0+1)
narr_o2p = fltarr(n_elements(L),max_theta*2.0+1)

!p.multi=[0,1,2]
for j = 0,n_elements(L)-1 do begin
   restore,'torus_profile_1'+L(j)+'.sav'
   nar1 = nar
   restore,'torus_profile'+L(j)+'.sav'

   s = [0,fltarr(max_theta)]
   dtheta=1.0*!pi/180.

   for i = 1,max_theta do begin
      theta = i*!pi/180.
      ds = L(j)*Rj*cos(theta)*sqrt(1+3*sin(theta)^2)*dtheta
      s(i) = s(i-1) + ds
   endfor
   
   s = [-reverse(s),s(1:n_elements(s)-1)]

   plot,nar.nel,s/Rj,/xlog,yrange=[0,4],xrange=[1,3000]
   oplot,nar.nsp,s/Rj,linestyle=1
   oplot,nar.ns2p,s/Rj,linestyle=2
   oplot,nar.ns3p,s/Rj,linestyle=3
   oplot,nar.nop,s/Rj,linestyle=4
   oplot,nar.no2p,s/Rj,linestyle=5

  plot,nar1.nel,s/Rj,/xlog,yrange=[0,4],xrange=[1,3000]
   oplot,nar1.nsp,s/Rj,linestyle=1
   oplot,nar1.ns2p,s/Rj,linestyle=2
   oplot,nar1.ns3p,s/Rj,linestyle=3
   oplot,nar1.nop,s/Rj,linestyle=4
   oplot,nar1.no2p,s/Rj,linestyle=5

   wait,0.5

   narr_sp(j,*) = 100*abs(nar.nsp-nar1.nsp)/nar1.nsp
   narr_s2p(j,*) = 100*abs(nar.ns2p-nar1.ns2p)/nar1.ns2p
   narr_s3p(j,*) = 100*abs(nar.ns3p-nar1.ns3p)/nar1.ns3p
   narr_op(j,*) = 100*abs(nar.nop-nar1.nop)/nar1.nop
   narr_o2p(j,*) = 100*abs(nar.no2p-nar1.no2p)/nar1.no2p

   ;narr_sp(j,*) = (nar.nsp)
   ;narr_s2p(j,*) = (nar.ns2p)
   ;narr_s3p(j,*) = nar.ns3p
   ;narr_op(j,*) = nar.nop
   ;narr_o2p(j,*) = nar.no2p


   ;narr_sp(j,*) = (nar1.nsp)
   ;narr_s2p(j,*) = (nar1.ns2p)
   ;narr_s3p(j,*) = nar1.ns3p
   ;narr_op(j,*) = nar1.nop
   ;narr_o2p(j,*) = nar1.no2p


endfor

w =window(dimensions=[1000,1000])
im = contour(narr_sp,r,s/rj,n_levels=10,layout=[2,3,1],/current,$
             max_value=max(narr_sp),min_value=0.01,$
             /fill,rgb_table=33,xrange=[6,15],aspect_ratio=1,yrange=[-6,6])
cb = colorbar(target=im)

im = contour(narr_s2p,r,s/rj,n_levels=10,layout=[2,3,2],/current,$
             max_value=max(narr_s2p),min_value=0.01,$
             /fill,rgb_table=33,xrange=[6,15],aspect_ratio=1,yrange=[-6,6])
cb = colorbar(target=im)

im = contour(narr_s3p,r,s/rj,n_levels=10,layout=[2,3,3],/current,$
             max_value=max(narr_s3p),min_value=0.01,$
             /fill,rgb_table=33,xrange=[6,15],aspect_ratio=1,yrange=[-6,6])
cb = colorbar(target=im)

im = contour(narr_op,r,s/rj,n_levels=10,layout=[2,3,4],/current,$
             max_value=max(narr_op),min_value=0.01,$
             /fill,rgb_table=33,xrange=[6,15],aspect_ratio=1,yrange=[-6,6])
cb = colorbar(target=im)

im = contour(narr_o2p,r,s/rj,n_levels=10,layout=[2,3,5],/current,$
             max_value=max(narr_o2p),min_value=0.01,$
             /fill,rgb_table=33,xrange=[6,15],aspect_ratio=1,yrange=[-6,6])
cb = colorbar(target=im)




end
