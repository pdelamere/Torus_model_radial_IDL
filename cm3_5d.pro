;-------------------------------------------------------------------
pro cm3_5d
;-------------------------------------------------------------------
@common

Te0 = 5.0
Ti0 = 70.0
Teh0 = 50.0
fh = 0.0025
trans = 1.0/(100.0*8.64e4)
net_source = 2.0e28
otos = 1.7
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

min_ns = 14e27
max_ns = 18e27
dns = 2e27
nns = round((max_ns-min_ns)/dns)
nsarr = min_ns + findgen(nns+1)*dns


min_otos = 1.9
max_otos = 2.0
dotos = 0.1
nos = (max_otos - min_otos)/dotos
osarr = min_otos + findgen(nos+1)*dotos

min_tau = 60.0
max_tau = 70.0
dtau = 5
ntau = (max_tau - min_tau)/dtau
tauarr = min_tau + findgen(ntau+1)*dtau


min_teh = 40.0
max_teh = 50.0
dteh = 5.0
nteh = (max_teh - min_teh)/dteh
teharr = min_teh + findgen(nteh+1)*dteh


min_feh = 0.0024
max_feh = 0.0027
dfeh = 0.0001
nfeh = (max_feh - min_feh)/dfeh
feharr = min_feh + findgen(nfeh+1)*dfeh

narr = fltarr(nns+1,ntau+1,nfeh+1,nteh+1,nos+1,7)
parr = fltarr(nns+1,ntau+1,nfeh+1,nteh+1,nos+1,6)
tarr = fltarr(nns+1,ntau+1,nfeh+1,nteh+1,nos+1,6)
xarr = fltarr(nns+1)
yarr = fltarr(ntau+1)
;carr = fltarr(nns+1,ntau+1,nfeh+1,nteh+1,nos+1,)  ;convergence
rarr = fltarr(nns+1,ntau+1,nfeh+1,nteh+1,nos+1,4)  ;energy input rates
nlarr = fltarr(nns+1,ntau+1,nfeh+1,nteh+1,nos+1,26) 
ecarr = fltarr(nns+1,ntau+1,nfeh+1,nteh+1,nos+1,19)

;ep = 1e-8
ep = 1e-1
ept = 1e-1
cnt = 0.0
cnt_tot = float(nns)*(ntau+1)*(nfeh+1)*(nteh+1)*(nos+1)

for i = 0,nns do begin
for j = 0,ntau do begin
for k = 0,nfeh do begin
for l = 0,nteh do begin
for m = 0,nos do begin
      net_source = min_ns + i*dns
      tau = min_tau + j*dtau
      fh = min_feh + k*dfeh  
      teh0 = min_teh + l*dteh
      otos = min_otos + m*dotos

      trans = 1.0/(tau*8.64e4)
      print,'src, tau, fh, teh, o/s...',net_source,tau,fh,teh0,otos
      print,'Percent complete.........',100.*cnt/cnt_tot,cnt,cnt_tot
;      xarr(i) = net_source
;      yarr(j) = otos
;      yarr(j) = trans
      JUMP:      
      cm3_model,n,r,T,src,lss,temps,dens,nl,Ec
      sz = size(src) & z = sz(1)-1 & zm = sz(1)-2
      carr = (src(z).sp + src(z).s2p + src(z).s3p + src(z).op + $
                  src(z).o2p)/(lss(z).sp + lss(z).s2p + lss(z).s3p + $
                  lss(z).op + lss(z).o2p)
      rarr(i,j,k,l,m,0) = r.Sei
      rarr(i,j,k,l,m,1) = r.Secx
      rarr(i,j,k,l,m,2) = r.Oei
      rarr(i,j,k,l,m,3) = r.Oecx

      nlarr(i,j,k,l,m,0) = nl.sei
      nlarr(i,j,k,l,m,1) = nl.seih
      nlarr(i,j,k,l,m,2) = nl.oei
      nlarr(i,j,k,l,m,3) = nl.oeih
      nlarr(i,j,k,l,m,4) = nl.scx_tot
      nlarr(i,j,k,l,m,5) = nl.ocx_tot
      nlarr(i,j,k,l,m,6) = nl.sion_tot
      nlarr(i,j,k,l,m,7) = nl.oion_tot
      nlarr(i,j,k,l,m,8) = nl.k1
      nlarr(i,j,k,l,m,9) = nl.k2
      nlarr(i,j,k,l,m,10) = nl.k3
      nlarr(i,j,k,l,m,11) = nl.k4
      nlarr(i,j,k,l,m,12) = nl.k5
      nlarr(i,j,k,l,m,13) = nl.k6
      nlarr(i,j,k,l,m,14) = nl.k7
      nlarr(i,j,k,l,m,15) = nl.k8
      nlarr(i,j,k,l,m,16) = nl.k9
      nlarr(i,j,k,l,m,17) = nl.k10
      nlarr(i,j,k,l,m,18) = nl.k11
      nlarr(i,j,k,l,m,19) = nl.k12
      nlarr(i,j,k,l,m,20) = nl.k14
      nlarr(i,j,k,l,m,21) = nl.fast_S_k1
      nlarr(i,j,k,l,m,22) = nl.fast_S_k3
      nlarr(i,j,k,l,m,23) = nl.fast_O_k5
      nlarr(i,j,k,l,m,24) = nl.fast_O_k7
      nlarr(i,j,k,l,m,25) = nl.fast_O_k9

      ecarr(i,j,k,l,m,0) = Ec.s_ion
      ecarr(i,j,k,l,m,1) = Ec.s_ion_h
      ecarr(i,j,k,l,m,2) = Ec.s_cx
      ecarr(i,j,k,l,m,3) = Ec.o_ion
      ecarr(i,j,k,l,m,4) = Ec.o_ion_h
      ecarr(i,j,k,l,m,5) = Ec.o_cx
      ecarr(i,j,k,l,m,6) = Ec.s_tot_in
      ecarr(i,j,k,l,m,7) = Ec.o_tot_in
      ecarr(i,j,k,l,m,8) = Ec.P_pu
      ecarr(i,j,k,l,m,9) = Ec.eh_eq
      ecarr(i,j,k,l,m,10) = Ec.P_in
      ecarr(i,j,k,l,m,11) = Ec.Puv
      ecarr(i,j,k,l,m,12) = Ec.Pfast
      ecarr(i,j,k,l,m,13) = Ec.Ptrans
      ecarr(i,j,k,l,m,14) = Ec.Ptrans_eh
      ecarr(i,j,k,l,m,15) = Ec.P_out
      ecarr(i,j,k,l,m,16) = Ec.in_out
      ecarr(i,j,k,l,m,17) = Ec.ion_e_eq

;      print,'convergence index...',carr(i,j)
;      if (((dens(z).nsp - dens(zm).nsp) lt ep) and $
;          ((dens(z).ns2p - dens(zm).ns2p) lt ep) and $
;          ((dens(z).ns3p - dens(zm).ns3p) lt ep) and $
;          ((dens(z).nop - dens(zm).nop) lt ep) and $
;          ((dens(z).no2p - dens(zm).no2p) lt ep) and $
;          ((dens(z).nel - dens(zm).nel) lt ep) and $
;          ((temps(z).Tsp - temps(zm).Tsp) lt ept) and $
;          ((temps(z).Ts2p - temps(zm).Ts2p) lt ept) and $
;          ((temps(z).Ts3p - temps(zm).Ts3p) lt ept) and $
;          ((temps(z).Top - temps(zm).Top) lt ept) and $
;          ((temps(z).To2p - temps(zm).To2p) lt ept) and $
;          ((temps(z).Telec - temps(zm).Telec) lt ept) and $
;          (carr lt 1.01) and (carr gt 0.999)) then begin
;         cont = 'false'
;        endif else begin
;         cont = 'true'
;         goto,JUMP 
;	endelse
;      print,n
;      carr = (src(z).sp + src(z).s2p + src(z).s3p + src(z).op + $
;                  src(z).o2p)/(lss(z).sp + lss(z).s2p + lss(z).s3p + $
;                  lss(z).op + lss(z).o2p)
;      print,'convergence index...',carr

      narr(i,j,k,l,m,0) = n.nsp
      narr(i,j,k,l,m,1) = n.ns2p
      narr(i,j,k,l,m,2) = n.ns3p
      narr(i,j,k,l,m,3) = n.nop
      narr(i,j,k,l,m,4) = n.no2p
      narr(i,j,k,l,m,5) = n.ns
      narr(i,j,k,l,m,6) = n.no
      parr(i,j,k,l,m,0) = r.Puv
      parr(i,j,k,l,m,1) = r.psp
      parr(i,j,k,l,m,2) = r.ps2p
      parr(i,j,k,l,m,3) = r.ps3p
      parr(i,j,k,l,m,4) = r.pop
      parr(i,j,k,l,m,5) = r.po2p
      tarr(i,j,k,l,m,0) = T.Tel
      tarr(i,j,k,l,m,1) = T.Tsp
      tarr(i,j,k,l,m,2) = T.Ts2p
      tarr(i,j,k,l,m,3) = T.Ts3p
      tarr(i,j,k,l,m,4) = T.Top
      tarr(i,j,k,l,m,5) = T.To2p
;      cont = 'true'

      ns0 = n.ns      ;initial S neutral density (cm^-3)
      no0 = n.no     ;initial O neutral density
      nsp0= n.nsp        ;initial S+ density
      ns2p0 = n.ns2p       ;initial S++ density
      ns3p0 = n.ns3p    ;initial S+++ density
      nop0 = n.nop     ;initial O+ density
      no2p0 = n.no2p      ;initial O++ density 
      nel0 = n.nel       ;initial electron density

      Te0 = 5.0
      Ti0 = 70.0

      

cnt = cnt+1
save,filename='cm3_5d_out.sav',narr,tarr,parr,xarr,yarr,rarr,nlarr,ecarr,$
  nsarr,osarr,tauarr,feharr,teharr

   NML:
   
   endfor
endfor
endfor
endfor
endfor



;make_cont_src_tau

return
end
;-------------------------------------------------------------------


