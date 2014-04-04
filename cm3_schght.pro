;---------------------------------------------------------------------
pro sclhgt,s3,L,theta_G_0,theta_G,Z,A,T,anis,n0,n1,fmin
;---------------------------------------------------------------------
;s3....system III longitude
;L.....L-shell
;theta_G_0....Jovigraphic latitude at S_0
;theta_G......Jovigraphic latitude at S
;Z............charge number
;A............atomic mass number
;T............temperature_Perp
;anis.........temperature anisotropy
;n0...........density at S_0
;n1...........density at S

q=1.6e-19       ;electronic charge
mp = 1.67e-27   ;mass of proton
omega = 1.76e-4 ;Jovian angular velocity radians/sec
s30 = 148.0*!dtor
alphat = 9.6*!dtor
alpha = asin(sin(alphat)*sin(s3-s30))
s3t = 292.0*!dtor
d_offset = 0.131   ;R_J
theta_M_0 = theta_G_0 - alpha
theta_M = theta_G - alpha
theta_C_0 = theta_M_0 + alpha/3
theta_C = theta_M + alpha/3

RJ = 7.14e4*1e3  ;meters  
factor = 0.825   ;0.5*mp*omega^2*RJ^2/q

R = L*cos(theta_M)^2
R0 = L*cos(theta_M_0)^2

;print,R,R0

fcent = R^2*cos(theta_G)^2 - R0^2*cos(theta_G_0)^2

f1 = 1 + 3.0*sin(theta_M_0)^2
f2 = 1 + 3.0*sin(theta_M)^2
f3 = (cos(theta_M)/cos(theta_M_0))^6

fmag = alog(sqrt(f1/f2)*f3)

phi = double(0.0)

;p = 50-findgen(100)*10


;plot,[-5,25],[-20,20],/nodata,xsty=4,ysty=4
;axis,xax=0,0,0
;axis,yax=0,0,0
;for i = 0,100 do begin
;p = 50 - i
;e1 = Z*p*anis/T
;e2 = factor*A*fcent*anis/T
;e3 = (anis - 1)*fmag
;ee = e1+e2-e3

;g = total(Z*n0*exp(ee))
;plots,p,g,/data,psym=1
;endfor
;;wait,0.5

it = 10
for i = 0,it-1 do begin

f = 0.0
df = 0.0

e1 = double(Z*phi*anis/T)
e2 = double(factor*A*fcent*anis/T)
e3 = double((anis - 1)*fmag)
ee = e1+e2-e3
n1 = double(n0*exp(ee))
f = double(total(Z*n1))
df = double(total(Z*Z*n1*anis/T))
if (abs(f) lt fmin) then goto, JUMP
print,df
if (df lt 0.1) then begin
   print,'Error...will not converge'
   goto, JUMP       
endif
phi = phi-f/df
print,'i,phi...',i,phi
endfor
if (abs(f) gt fmin) then begin
   print,'Error...check convergence...'
   print,abs(f),fmin
endif
JUMP:

return
end
;---------------------------------------------------------------------


;---------------------------------------------------------------------
pro sclhgt_kappa,s3,L,theta_G_0,theta_G,Z,A,T,kappa,n0,n1,fmin,Tk
;---------------------------------------------------------------------
;s3....system III longitude
;L.....L-shell
;theta_G_0....Jovigraphic latitude at S_0
;theta_G......Jovigraphic latitude at S
;Z............charge number
;A............atomic mass number
;T............temperature_Perp
;kappa........parameter for kappa function
;n0...........density at S_0
;n1...........density at S

gamma = 1 - 1/(kappa - 0.5)
alpha_i = T/n0^(gamma-1)

q=1.6e-19       ;electronic charge
mp = 1.67e-27   ;mass of proton
omega = 1.76e-4 ;Jovian angular velocity radians/sec
s30 = 148.0*!dtor
alphat = 9.6*!dtor
alpha = asin(sin(alphat)*sin(s3-s30))
s3t = 292.0*!dtor
d_offset = 0.131   ;R_J
theta_M_0 = theta_G_0 - alpha
theta_M = theta_G - alpha
theta_C_0 = theta_M_0 + alpha/3
theta_C = theta_M + alpha/3

RJ = 7.14e4*1e3  ;meters  
factor = 0.825   ;0.5*mp*omega^2*RJ^2/q
g_fac = (gamma - 1)/(alpha_i*gamma)
;print,g_fac

R = L*cos(theta_M)^2
R0 = L*cos(theta_M_0)^2

;print,R,R0

fcent = R^2*cos(theta_G)^2 - R0^2*cos(theta_G_0)^2

;f1 = 1 + 3.0*sin(theta_M_0)^2
;f2 = 1 + 3.0*sin(theta_M)^2
;f3 = (cos(theta_M)/cos(theta_M_0))^6

;fmag = alog(sqrt(f1/f2)*f3)

phi = 0.0

;p = 50-findgen(100)*10

;e1 = g_fac(0)*Z(0)*p
;e2 = g_fac(0)*factor*A(0)*fcent
;e3 = n0(0)^(gamma-1)

;g = (e1+e2+e3)^(1/(gamma-1))

;plot,p,g
;print,g
;stop

;plot,[-5,25],[-2000,2000],/nodata,xsty=4,ysty=4
;axis,xax=0,0,0
;axis,yax=0,0,0
;for i = 0,100 do begin
;p = 50 - i
;e1 = g_fac*Z*p
;e2 = g_fac*factor*A*fcent
;e3 = n0^(gamma-1)
;ee = e1+e2+e3

;g = total(Z*ee^(1/(gamma-1)))
;plots,p,g,/data,psym=1
;;print,p,g
;endfor
;wait,0.1


it = 10
for i = 0,it-1 do begin

f = 0.0
df = 0.0

e1 = g_fac*Z*phi
e2 = g_fac*factor*A*fcent
e3 = n0^(gamma-1)
ee = e1+e2+e3
n1 = ee^(1/(gamma-1))
;print,n1,exp(ee)
f = total(Z*n1)
;print,i,f
df = total(Z*Z*ee^((2-gamma)/(gamma-1))/(alpha_i*gamma))
if (abs(f) lt fmin) then goto, JUMP
if (df eq 0.0) then begin
   print,'Error...will not converge'
   goto, JUMP       
endif
phi = phi-f/df
;print,'i,f,phi...',i,f,phi
endfor
if (abs(f) gt fmin) then begin
   print,'Error...check convergence...'
   print,abs(f), fmin
endif
JUMP:

Tk = alpha_i*n1^(gamma-1)
 
return
end
;---------------------------------------------------------------------


;---------------------------------------------------------------------
pro f_df_ions,s3,L,theta_G_0,theta_G,Z,A,T,anis,n0,n1,fions,dfions,phi
;---------------------------------------------------------------------
;s3....system III longitude
;L.....L-shell
;theta_G_0....Jovigraphic latitude at S_0
;theta_G......Jovigraphic latitude at S
;Z............charge number
;A............atomic mass number
;T............temperature_Perp
;anis.........temperature anisotropy
;n0...........density at S_0
;n1...........density at S

q=1.6e-19       ;electronic charge
mp = 1.67e-27   ;mass of proton
omega = 1.76e-4 ;Jovian angular velocity radians/sec
s30 = 148.0*!dtor
s3t = 292.0*!dtor
;alphat = 9.6*!dtor
alphat = 0.0
alpha = asin(sin(alphat)*sin(s3-s3t))
d_offset = 0.131   ;R_J
theta_M_0 = theta_G_0 - alpha
theta_M = theta_G - alpha
theta_C_0 = theta_M_0 + alpha/3
theta_C = theta_M + alpha/3

RJ = 7.14e4*1e3  ;meters  
factor = 0.825   ;0.5*mp*omega^2*RJ^2/q

R = L*cos(theta_M)^2
R0 = L*cos(theta_M_0)^2

;print,R,R0

fcent = R^2*cos(theta_G)^2 - R0^2*cos(theta_G_0)^2

f1 = 1 + 3.0*sin(theta_M_0)^2
f2 = 1 + 3.0*sin(theta_M)^2
f3 = (cos(theta_M)/cos(theta_M_0))^6

fmag = alog(sqrt(f1/f2)*f3)

e1 = double(Z*phi*anis/T)
e2 = double(factor*A*fcent*anis/T)
e3 = double((anis - 1)*fmag)
ee = e1+e2-e3
n1 = double(n0*exp(ee))

fions = (Z*n1)
dfions = (Z*Z*n1*anis/T)

return
end
;---------------------------------------------------------------------

;---------------------------------------------------------------------
pro f_df_elec_kappa,s3,L,theta_G_0,theta_G,Z,A,T,kappa,n0,n1,Tk,felec,dfelec,phi
;---------------------------------------------------------------------
;s3....system III longitude
;L.....L-shell
;theta_G_0....Jovigraphic latitude at S_0
;theta_G......Jovigraphic latitude at S
;Z............charge number
;A............atomic mass number
;T............temperature_Perp
;kappa........parameter for kappa function
;n0...........density at S_0
;n1...........density at S

gamma = 1 - 1/(kappa - 0.5)
alpha_i = T/n0^(gamma-1)

q=1.6e-19       ;electronic charge
mp = 1.67e-27   ;mass of proton
omega = 1.76e-4 ;Jovian angular velocity radians/sec
s30 = 149.0*!dtor
s3t = 292.0*!dtor
;alphat = 9.6*!dtor
alphat = 0.0*!dtor
alpha = asin(sin(alphat)*sin(s3-s3t))
d_offset = 0.131   ;R_J
theta_M_0 = theta_G_0 - alpha
theta_M = theta_G - alpha
theta_C_0 = theta_M_0 + alpha/3
theta_C = theta_M + alpha/3

RJ = 7.14e4*1e3  ;meters  
factor = 0.825   ;0.5*mp*omega^2*RJ^2/q
g_fac = (gamma - 1)/(alpha_i*gamma)

R = L*cos(theta_M)^2
R0 = L*cos(theta_M_0)^2

fcent = R^2*cos(theta_G)^2 - R0^2*cos(theta_G_0)^2
e1 = g_fac*Z*phi
e2 = g_fac*factor*A*fcent
e3 = n0^(gamma-1)
ee = e1+e2+e3
n1 = ee^(1/(gamma-1))

;if (n1 gt n0) then print,'density inconsistency...',n1,n0,theta_G*!radeg

felec = Z*n1
dfelec = (Z*Z*ee^((2-gamma)/(gamma-1))/(alpha_i*gamma))

Tk = alpha_i*n1^(gamma-1)

return
end
;---------------------------------------------------------------------


;---------------------------------------------------------------------
pro cm3_expand,n,narr,T,nar,max_theta
; n contains the equatorial densities
; narr is density array containing [7,61] elements
; nar is a structure containing  61 elements for +/- 30 deg above and 
;     below the equator plane.
; T contains the electron and ion temps.
;---------------------------------------------------------------------
;species...e,O+,O++,S+,S++,S+++

;!p.multi=[0,1,2]
q = 1.6e-19

A = [1/1685.17,16.0,16.0,32.0,32.0,32.0,1/1685.]
;A = [0.0,16.0,16.0,32.0,32.0,32.0]
Z = [-1.0,1.0,2.0,1.0,2.0,3.0,-1.0]
;T = [5.0,70.0,70.0,70.0,70.0,70.0]
Tarr = [T.Tel,T.Top,T.To2p,T.Tsp,T.Ts2p,T.Ts3p,T.Telh]

n0 = [n.nel,n.nop,n.no2p,n.nsp,n.ns2p,n.ns3p,n.nelh]

;n0 = [1000.0,100.0,200.0,500.0,8.0]
;n0 = [total(n0),n0]

;n1_ions = n0(1:*)
n1_ions = n0(*)

n1_elec = n0(0)

s3 = 148.0*!dtor
L = 6.0
theta_G_0 = 0.0
;theta_G = 15.0*!dtor
;theta_G = 0.0*!dtor
anis = 1.0
fmin = 0.001
;max_theta = 30  
kappa = 2.4

narr = fltarr(n_elements(n0),2.0*max_theta + 1.0)
;nkappa = fltarr(n_elements(n0),2.0*max_theta + 1.0)
Tkappa = fltarr(2.0*max_theta + 1.0)
vwght = fltarr(2.0*max_theta + 1.0)
for i = 0,2.0*max_theta do begin
   theta_G = (-max_theta+i)*!dtor
   vwght(i) = cos(theta_G)^7

   phi = 0.0
   it = 20
   for j = 0,it-1 do begin

      f = 0.0
      df = 0.0

;      f_df_ions,s3,L,theta_G_0,theta_G,Z(1:*),A(1:*),Tarr(1:*),anis,n0(1:*),$
;                n1_ions,fions,dfions,phi
      f_df_ions,s3,L,theta_G_0,theta_G,Z(*),A(*),Tarr(*),anis,n0(*),$
                n1_ions,fions,dfions,phi
;      f_df_elec_kappa,s3,L,theta_G_0,theta_G,Z(0),A(0),Tarr(0),kappa,n0(0),$
;                      n1_elec,Tk,felec,dfelec,phi

;      n1 = [n1_elec,n1_ions]
      n1 = [n1_ions]
;      f = total([felec,fions])
;      df = total([dfelec,dfions])
      f = total([fions])
      df = total([dfions])

      if (abs(f) lt fmin) then goto, JUMP
      
      phi = phi-f/df

   endfor
   if (abs(f) gt fmin) then begin
      print,'Error...check convergence...'
      print,abs(f),fmin
   endif
   JUMP:

   narr(*,i) = n1
;   Tkappa(i) = Tk

endfor

;if (n_params(0) eq 4) then begin
nar(*).nel = reform(narr(0,*))
nar(*).nop = reform(narr(1,*))
nar(*).no2p = reform(narr(2,*))
nar(*).nsp = reform(narr(3,*))
nar(*).ns2p = reform(narr(4,*))
nar(*).ns3p = reform(narr(5,*))
nar(*).nelh = reform(narr(6,*))
;endif

;!p.multi=[0,3,3]
;plot,nar.nel,title='nel'
;plot,nar.no,title='no'
;plot,nar.ns,title='ns'
;plot,nar.nsp,title='nsp'
;plot,nar.ns2p,title='ns2p'
;plot,nar.ns3p,title='ns3p'
;plot,nar.nop,title='nop'
;plot,nar.no2p,title='no2p'
;!p.multi=[0,1,1]



;      plot_oi,narr(0,*),max_theta-findgen(2.0*max_theta+1),$
;         yrange=[-max_theta,max_theta],/ysty,$
;         title='System III Long = '+strtrim(string(s3*!radeg),2),$
;         ytitle = 'Jovigraphic latitude',xtitle='density (cm!u-3!n)',$
;         xrange=[1,3.0*n0(0)],/xsty
;      oplot,narr(1,*),max_theta-findgen(2.0*max_theta+1),linestyle=1
;      oplot,narr(2,*),max_theta-findgen(2.0*max_theta+1),linestyle=2
;      oplot,narr(3,*),max_theta-findgen(2.0*max_theta+1),linestyle=3
;      oplot,narr(4,*),max_theta-findgen(2.0*max_theta+1),linestyle=4
;      oplot,narr(5,*),max_theta-findgen(2.0*max_theta+1),linestyle=5
;      legend,['e!u-!n','O!u+!n','O!u++!n','S!u+!n','S!u++!n','S!u+++!n'],$
;         linestyle=[0,1,2,3,4,5],/right

theta = max_theta-findgen(2.0*max_theta)
save,filename='torus_profile.sav',theta,narr,Tkappa,max_theta,S3,n0,kappa

;sum_n = total(narr(0,*)*vwght)
;sum_nT = total(narr(0,*)*Tkappa*vwght)
;T.Tel_el = sum_nT/sum_n
;;print,'Fluxtube averaged electron temp....................',T.Tel_el

;sum_nn = total(narr(1,*)*vwght)
;sum_nnT = total(narr(1,*)*Tkappa*vwght)
;T.Tel_op = sum_nnT/sum_nn
;;print,'Fluxtube averaged electron temp, (O+ weighted).....',T.Tel_op

;sum_nn = total(narr(2,*)*vwght)
;sum_nnT = total(narr(2,*)*Tkappa*vwght)
;T.Tel_o2p = sum_nnT/sum_nn
;;print,'Fluxtube averaged electron temp, (O++ weighted)....',T.Tel_o2p

;sum_nn = total(narr(3,*)*vwght)
;sum_nnT = total(narr(3,*)*Tkappa*vwght)
;T.Tel_sp = sum_nnT/sum_nn
;;print,'Fluxtube averaged electron temp, (S+ weighted).....',T.Tel_sp

;sum_nn = total(narr(4,*)*vwght)
;sum_nnT = total(narr(4,*)*Tkappa*vwght)
;T.Tel_s2p = sum_nnT/sum_nn
;;print,'Fluxtube averaged electron temp, (S++ weighted)....',T.Tel_s2p

;sum_nn = total(narr(5,*)*vwght)
;sum_nnT = total(narr(5,*)*Tkappa*vwght)
;T.Tel_s3p = sum_nnT/sum_nn
;;print,'Fluxtube averaged electron temp, (S+++ weighted)...',T.Tel_s3p

return
end
;---------------------------------------------------------------------


;---------------------------------------------------------------------
pro run_ekappa,n,T
;---------------------------------------------------------------------
;ribbon
Te0 = 5.0
Ti0 = 70.0
nsp0 = 300.0
ns2p0 = 500.0
ns3p0 = 100.0
nop0 = 1400.0
no2p0 = 50.0
nel0 = nsp0 + 2.0*ns2p0 + 3.0*ns3p0 + nop0 + 2.0*no2p0

;east
;Te0 = 4.0
;Ti0 = 15.0
;nsp0 = 400.0
;ns2p0 = 10.0
;ns3p0 = 1.0
;nop0 = 400.0
;no2p0 = 10.0
;nel0 = nsp0 + 2.0*ns2p0 + 3.0*ns3p0 + nop0 + 2.0*no2p0

;west
;Te0 = 5.0
;Ti0 = 70.0
;nsp0 = 200.0
;ns2p0 = 320.0
;ns3p0 = 50.0
;nop0 = 900.0
;no2p0 = 10.0
;nel0 = nsp0 + 2.0*ns2p0 + 3.0*ns3p0 + nop0 + 2.0*no2p0


n = {density, nel: nel0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
              nop: nop0, no2p: no2p0}
T = {temp, Tel: Te0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, Top: Ti0, $
           To2p: Ti0, Tel_sp: Te0, Tel_s2p: Te0, Tel_s3p: Te0, Tel_op: Te0, $
           Tel_o2p: Te0, Tel_el: Te0}

cm3_ekappa,n,T
      !p.multi=[0,1,2]
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


return
end
;---------------------------------------------------------------------



;---------------------------------------------------------------------
pro ftave_ekappa,n,T
;---------------------------------------------------------------------

;ribbon
Te0 = 5.0
Ti0 = 70.0
nsp0 = 300.0
ns2p0 = 500.0
ns3p0 = 100.0
nop0 = 1400.0
no2p0 = 50.0
nel0 = nsp0 + 2.0*ns2p0 + 3.0*ns3p0 + nop0 + 2.0*no2p0

n = {density, nel: nel0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
              nop: nop0, no2p: no2p0}
T = {temp, Tel: Te0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, Top: Ti0, $
           To2p: Ti0, Tel_sp: Te0, Tel_s2p: Te0, Tel_s3p: Te0, Tel_op: Te0, $
           Tel_o2p: Te0, Tel_el: Te0}

Rj = 7.14e4 ;km
L = 6.0

theta = (30-findgen(61))*!dtor

v = 2.0*4*!pi*Rj^3*L^2*total(cos(theta)^7)  ;factor of 2 for up and down

varr = (cos(theta)^7)^1
;print,varr


cm3_ekappa,n,T

restore,'torus_profile.sav'
print,narr(0,*)
print,'nel*nop...',(total(narr(0,*)*narr(1,*)*varr)/total(varr))/$
       (n.nel*n.nop)
print,'nel*no2p...',(total(narr(0,*)*narr(2,*)*varr)/total(varr))/$
       (n.nel*n.no2p)

print,'nel*nsp...',(total(narr(0,*)*narr(3,*)*varr)/total(varr))/$
       (n.nel*n.nsp)
print,'nel*ns2p...',(total(narr(0,*)*narr(4,*)*varr)/total(varr))/$
       (n.nel*n.ns2p)
print,'nel*ns2p...',(total(narr(0,*)*narr(5,*)*varr)/total(varr))/$
       (n.nel*n.ns3p)



      !p.multi=[0,1,2]
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
      plots,[nop0,nop0],[!y.crange(0),!y.crange(1)],linestyle=1

      legend,['e!u-!n','O!u+!n','O!u++!n','S!u+!n','S!u++!n','S!u+++!n'],$
         linestyle=[0,1,2,3,4,5],/right
      plot_oi,Tkappa,max_theta-findgen(2.0*max_theta),$
         yrange=[-max_theta,max_theta],/ysty,$
         title='Electron Temp (Kappa = '+strtrim(string(kappa),2)+')',$
         ytitle = 'Jovigraphic latitude',xtitle='Temperature (eV)',$
         xrange=[1,100],/xsty



return
end
;---------------------------------------------------------------------


;---------------------------------------------------------------------
pro get_lat_dist,n,T,nar,max_theta
; n contains the equatorial densities
; narr is density array containing [7,61] elements
; nar is a structure containing  61 elements for +/- 30 deg above and 
;     below the equator plane.
; T contains the electron and ion temps.
;---------------------------------------------------------------------
;species...e,O+,O++,S+,S++,S+++

;!p.multi=[0,1,2]
q = 1.6e-19

A = [1/1685.17,16.0,16.0,32.0,32.0,32.0,1/1685.]
;A = [0.0,16.0,16.0,32.0,32.0,32.0]
Z = [-1.0,1.0,2.0,1.0,2.0,3.0,-1.0]
;T = [5.0,70.0,70.0,70.0,70.0,70.0]
Tarr = [T.Tel,T.Top,T.To2p,T.Tsp,T.Ts2p,T.Ts3p,T.Telh]

n0 = [n.nel,n.nop,n.no2p,n.nsp,n.ns2p,n.ns3p,n.nelh]

n1_ions = n0(*)

n1_elec = n0(0)

s3 = 148.0*!dtor
L = 6.0
theta_G_0 = 0.0
;theta_G = 15.0*!dtor
;theta_G = 0.0*!dtor
anis = 1.0
fmin = 0.001
;max_theta = 30  
kappa = 2.4

narr = fltarr(n_elements(n0),2.0*max_theta + 1.0)
;nkappa = fltarr(n_elements(n0),2.0*max_theta + 1.0)
Tkappa = fltarr(2.0*max_theta + 1.0)
vwght = fltarr(2.0*max_theta + 1.0)
for i = 0,2.0*max_theta do begin
   theta_G = (-max_theta+i)*!dtor
   vwght(i) = cos(theta_G)^7

   phi = 0.0
   it = 20
   for j = 0,it-1 do begin

      f = 0.0
      df = 0.0

;      f_df_ions,s3,L,theta_G_0,theta_G,Z(1:*),A(1:*),Tarr(1:*),anis,n0(1:*),$
;                n1_ions,fions,dfions,phi
      f_df_ions,s3,L,theta_G_0,theta_G,Z(*),A(*),Tarr(*),anis,n0(*),$
                n1_ions,fions,dfions,phi
;      f_df_elec_kappa,s3,L,theta_G_0,theta_G,Z(0),A(0),Tarr(0),kappa,n0(0),$
;                      n1_elec,Tk,felec,dfelec,phi

;      n1 = [n1_elec,n1_ions]
      n1 = [n1_ions]
;      f = total([felec,fions])
;      df = total([dfelec,dfions])
      f = total([fions])
      df = total([dfions])

      if (abs(f) lt fmin) then goto, JUMP
      
      phi = phi-f/df

   endfor
   if (abs(f) gt fmin) then begin
      print,'Error...check convergence...'
      print,abs(f),fmin
   endif
   JUMP:

   narr(*,i) = n1
;   Tkappa(i) = Tk

endfor

;if (n_params(0) eq 4) then begin
nar(*).nel = reform(narr(0,*))
nar(*).nop = reform(narr(1,*))
nar(*).no2p = reform(narr(2,*))
nar(*).nsp = reform(narr(3,*))
nar(*).ns2p = reform(narr(4,*))
nar(*).ns3p = reform(narr(5,*))
nar(*).nelh = reform(narr(6,*))
;endif

;!p.multi=[0,3,3]
;plot,nar.nel,title='nel'
;plot,nar.no,title='no'
;plot,nar.ns,title='ns'
;plot,nar.nsp,title='nsp'
;plot,nar.ns2p,title='ns2p'
;plot,nar.ns3p,title='ns3p'
;plot,nar.nop,title='nop'
;plot,nar.no2p,title='no2p'
;!p.multi=[0,1,1]

;      plot_oi,narr(0,*),max_theta-findgen(2.0*max_theta+1),$
;         yrange=[-max_theta,max_theta],/ysty,$
;         title='System III Long = '+strtrim(string(s3*!radeg),2),$
;         ytitle = 'Jovigraphic latitude',xtitle='density (cm!u-3!n)',$
;         xrange=[1,3.0*n0(0)],/xsty
;      oplot,[n.nel,n.nel],[!y.crange(0),!y.crange(1)],linestyle=1
;      oplot,narr(1,*),max_theta-findgen(2.0*max_theta+1),linestyle=1
;      oplot,[n.nop,n.nop],[!y.crange(0),!y.crange(1)],linestyle=1
;      oplot,narr(2,*),max_theta-findgen(2.0*max_theta+1),linestyle=2
;      oplot,[n.no2p,n.no2p],[!y.crange(0),!y.crange(1)],linestyle=1
;      oplot,narr(3,*),max_theta-findgen(2.0*max_theta+1),linestyle=3
;      oplot,[n.nsp,n.nsp],[!y.crange(0),!y.crange(1)],linestyle=1
;      oplot,narr(4,*),max_theta-findgen(2.0*max_theta+1),linestyle=4
;      oplot,[n.ns2p,n.ns2p],[!y.crange(0),!y.crange(1)],linestyle=1
;      oplot,narr(5,*),max_theta-findgen(2.0*max_theta+1),linestyle=5
;      oplot,[n.ns3p,n.ns3p],[!y.crange(0),!y.crange(1)],linestyle=1

;      legend,['e!u-!n','O!u+!n','O!u++!n','S!u+!n','S!u++!n','S!u+++!n'],$
;         linestyle=[0,1,2,3,4,5],/right
;wait,0.1



theta = max_theta-findgen(2.0*max_theta)
save,filename='torus_profile.sav',theta,narr,Tkappa,max_theta,S3,n0,kappa

;sum_n = total(narr(0,*)*vwght)
;sum_nT = total(narr(0,*)*Tkappa*vwght)
;T.Tel_el = sum_nT/sum_n
;;print,'Fluxtube averaged electron temp....................',T.Tel_el

;sum_nn = total(narr(1,*)*vwght)
;sum_nnT = total(narr(1,*)*Tkappa*vwght)
;T.Tel_op = sum_nnT/sum_nn
;;print,'Fluxtube averaged electron temp, (O+ weighted).....',T.Tel_op

;sum_nn = total(narr(2,*)*vwght)
;sum_nnT = total(narr(2,*)*Tkappa*vwght)
;T.Tel_o2p = sum_nnT/sum_nn
;;print,'Fluxtube averaged electron temp, (O++ weighted)....',T.Tel_o2p

;sum_nn = total(narr(3,*)*vwght)
;sum_nnT = total(narr(3,*)*Tkappa*vwght)
;T.Tel_sp = sum_nnT/sum_nn
;;print,'Fluxtube averaged electron temp, (S+ weighted).....',T.Tel_sp

;sum_nn = total(narr(4,*)*vwght)
;sum_nnT = total(narr(4,*)*Tkappa*vwght)
;T.Tel_s2p = sum_nnT/sum_nn
;;print,'Fluxtube averaged electron temp, (S++ weighted)....',T.Tel_s2p

;sum_nn = total(narr(5,*)*vwght)
;sum_nnT = total(narr(5,*)*Tkappa*vwght)
;T.Tel_s3p = sum_nnT/sum_nn
;;print,'Fluxtube averaged electron temp, (S+++ weighted)...',T.Tel_s3p

return
end
;---------------------------------------------------------------------
