;-------------------------------------------------------------------
PRO cm3_radial,prms,mr_r,peuv_peuvh,nelec6,tionavg
;PRO cm3_radial
;-------------------------------------------------------------------
@common
@params

mr_r = fltarr(12)

fh = prms(0)            ;fh (L = 6.0)
fh_alpha = prms(1)          ;fh power law

net_source = prms(2)    ;net source #/s
net_source_alpha = prms(3)  ;power law for extended neutral source

DLL_0 = prms(4)         ;DLL (L = 6.0)
DLL_alpha = prms(5)        ;DLL power law

Teh0 = prms(6)
Teh0_alpha = prms(7)

;otos = prms(8)

window,0,xsize=600,ysize=800
window,1,xsize=600,ysize=800
window,2,xsize=600,ysize=800

read_nl2,s

;===================================================================================
;grid setup
;===================================================================================


bufcel = 0

;nintrvl= nbox;38+bufcel
nintrvl= 38+bufcel
;nbox = 18
;dL0 = 0.25
;dLarr = fltarr(nintrvl)
dLarr = fltarr(nbox+bufcel)
dLarr(*) = dL0
;dL = [dLarr]
;dL = [dLarr];,dL0+1.5*(findgen(nintrvl-nbox)/(nintrvl-nbox))^3]
dL = [dLarr,dL0+1.5*(findgen(nintrvl-nbox)/(nintrvl-nbox))^3]
max_theta = 30

n_theta = max_theta*2.0 + 1.0
;theta_off = 3.5*!dtor
Rj = 7.14e4 ;km
zoff = theta_off*6.0*Rj 
;n_width = 1.0  

;Lshell = 6.0 + dindgen(nintrvl)*dL
Lshell = fltarr(nintrvl)
Lshell(0) = 6.0 - dL0*bufcel
for i = 1,nintrvl-1 do begin
  Lshell(i) = Lshell(i-1) + dL(i)
  ;print,Lshell(i)
endfor

;===================================================================================
;variable declarations
;===================================================================================

Te0 = 5.0
Ti0 = 70.0
Tih0 = 70.0
;Teh0 = 30.0
Teh0_1 = Teh0
;fh = 0.0028
fc = 1.0-fh
fh_0 = fh
fh_1 = fh
trans = 0.0 ;1.0/(64.0*8.64e4)
trans_1 = trans
;net_source = 3.2e27
net_source_0 = net_source
net_source_arr = fltarr(nbox)
;otos = 2.0
otos_1 = otos
s_source = net_source/(1.0+otos)
o_source = net_source*otos/(otos + 1.0)

;dt0 = 10000.0      ;seconds
;runt = 0.1*8.64e4  ;seconds
;runt = 2*dt0  ;seconds
;nit = 400
;dt_trans = dt0/1000.
nit_trans = runt/dt_trans
dt = dt0
ntm = fix(runt/dt)

ns0 = 50.0      ;initial S neutral density (cm^-3)
no0 = 200.0     ;initial O neutral density
nsp0= 150.0        ;initial S+ density
ns2p0 = 600.0       ;initial S++ density
ns3p0 = 50.0      ;initial S+++ density
ns4p0 = 0.0
nop0 = 600.0       ;initial O+ density
no2p0 = 60.0      ;initial O++ density 
nsph0 = 0.01
noph0 = 0.01

;ns0 = 10.      ;initial S neutral density (cm^-3)
;no0 = 10.    ;initial O neutral density
;nsp0= 100.        ;initial S+ density
;ns2p0 = 10.
;ns3p0 = 10.
;ns4p0 = 0.0
;nop0 = 10.
;no2p0 = 10.
;nsph0 = 0.01
;noph0 = 0.01

nel0 = nsp0+2.0*ns2p0+3.0*ns3p0+4.0*ns4p0+nop0+2.0*no2p0+nsph0+noph0   ;initial electron density
nelh0 = 0.01*nel0



tot_peuv = fltarr(nintrvl)
tot_peuvh = fltarr(nintrvl)

nr = replicate({densr, nsp: 0.0d, ns2p: 0.0d, ns3p: 0.0d, ns4p: 0.0d, nop: 0.0d, $
                no2p: 0.0d, nel: 0.0d, nelh: 0.0d, ns: 0.0d, no: 0.0d,$
                nsph: 0.0d, noph: 0.0d},nintrvl)

tr = replicate({tempsr, Tsp: 0.0d, Ts2p: 0.0d, Ts3p: 0.0d, Ts4p: 0.0d, Top: 0.0d, $
                To2p: 0.0d, Tel: 0.0d, Telh: 0.0d,$
                Tsph: 0.0d, Toph: 0.0d},nintrvl)

mrr = replicate({mix_ratio, sp: 0.0d, s2p: 0.0d, s3p: 0.0d, op: 0.0d, $
                o2p: 0.0d},nintrvl)

nl2 = replicate({nlsqrd, nsp: 0.0d, ns2p: 0.0d, ns3p: 0.0d, ns4p: 0.0, nop: 0.0d, $
                no2p: 0.0d, nel: 0.0d},nintrvl)

nl2e = replicate({nlsqrde, nsp: 0.0d, ns2p: 0.0d, ns3p: 0.0d, ns4p: 0.0d, nop: 0.0d, $
                no2p: 0.0d, nel: 0.0d},nintrvl)

sr = replicate({source,neutral:0.0,nsp: 0.0, ns2p: 0.0, ns3p: 0.0, $
                       nop: 0.0, no2p: 0.0, nel: 0.0, trans:0.0},nintrvl)

lr = replicate({loss,neutral:0.0,nsp: 0.0, ns2p: 0.0, ns3p: 0.0, $
                       nop: 0.0, no2p: 0.0, nel: 0.0},nintrvl)

snTr = replicate({e_source,sp: 0.0, s2p: 0.0, s3p: 0.0, $
                       op: 0.0, o2p: 0.0, el: 0.0},nintrvl)

lnTr = replicate({e_loss,sp: 0.0, s2p: 0.0, s3p: 0.0, $
                       op: 0.0, o2p: 0.0, el: 0.0},nintrvl)

rr = replicate({radial_rates,is: 0.0, isp: 0.0, is2p: 0.0, io: 0.0, iop: 0.0, $
                      ish: 0.0, isph: 0.0, is2ph: 0.0, ioh: 0.0, ioph: 0.0, $
                      k0_sp: 0.0, k0_s2p: 0.0, k1_s: 0.0, k1_sp: 0.0, $ 
                      k2_s: 0.0, k2_s2p: 0.0, k3_s: 0.0, k3_s2p: 0.0, $
                      k4_s: 0.0, k4_s3p: 0.0, k5_o: 0.0, k5_op: 0.0, $
                      k6_o: 0.0, k6_o2p: 0.0, k7_o: 0.0, k7_o2p: 0.0, $
                      k8_o: 0.0, k8_sp: 0.0, k9_s: 0.0, k9_op: 0.0, $
                      k10_s: 0.0, k10_o2p: 0.0, k11_s: 0.0, k11_o2p: 0.0, $
                      k12_o: 0.0, k12_s2p: 0.0, k13_o2p: 0.0, k13_sp: 0.0, $
                      k14_o: 0.0, k14_s3p: 0.0, k15_o2p: 0.0, k15_s2p: 0.0, $
                      k16_s3p: 0.0, k16_sp: 0.0, $
                      nu_ee: 0.0, nu_sp_e: 0.0, nu_s2p_e: 0.0, nu_s3p_e: 0.0, $
                      nu_op_e: 0.0, nu_o2p: 0.0},nintrvl)

;cm3 model structures

;nar = replicate({density_ar, nel: nel0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
;          nop: nop0, no2p: no2p0, nelh: nelh0, ns: ns0, no: no0},n_theta)

;n1ar = replicate({density_ar, nel: nel0, nsp: nsp0, ns2p: ns2p0, $
;       ns3p: ns3p0,nop: nop0, no2p: no2p0, nelh: nelh0, ns: ns0, no: no0},n_theta)

n = {density_avg, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
              no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh, $
              nsp_pu: 0.0, ns2p_pu: 0.0, ns3p_pu: 0.0, nop_pu: 0.0, $
              no2p_pu: 0.0, nelh: nelh0, nsph: nsph0, noph: noph0}

T = {temp, Tel: Te0, Telh: Teh0, Tsp: Ti0, Ts2p: Ti0, $
           Ts3p: Ti0, Top: Ti0, $
           To2p: Ti0, Tsp_pu: 0.0, Ts2p_pu: 0.0, Ts3p_pu: 0.0, $
           Top_pu: 0.0, To2p_pu: 0.0, Tpu_s: 0.0, Tpu_o: 0.0, $
           Tel_sp: Te0, Tel_s2p: Te0, Tel_s3p: Te0, Tel_op: Te0, $
           Tel_o2p: Te0, Tel_el: Te0, Tsph: Tih0, Toph: Tih0}

nT = {energy, el: nel0*Te0, sp: nsp0*Ti0, s2p: ns2p0*Ti0, $
              s3p: ns3p0*Ti0, op: nop0*Ti0, o2p: no2p0*Ti0, $
              sp_pu: 0.0, s2p_pu: 0.0, s3p_pu: 0.0, op_pu: 0.0,$
              o2p_pu: 0.0} 

;advance density for improved Euler method, half time step
n1 = {density_1, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
                 no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh, $
              nsp_pu: 0.0, ns2p_pu: 0.0, ns3p_pu: 0.0, nop_pu: 0.0, $
              no2p_pu: 0.0, nelh: nelh0, nsph: nsph0, noph: noph0}
T1 = {temp_1, Tel: Te0, Telh: Teh0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, $
           Top: Ti0, $
           To2p: Ti0, Tsp_pu: 0.0, Ts2p_pu: 0.0, Ts3p_pu: 0.0, $
           Top_pu: 0.0, To2p_pu: 0.0, Tpu_s: 0.0, Tpu_o: 0.0, $
           Tel_sp: Te0, Tel_s2p: Te0, Tel_s3p: Te0, Tel_op: Te0, $
           Tel_o2p: Te0, Tel_el: Te0, Tsph: Tih0, Toph: Tih0}
nT1 = {energy_1, el: nel0*Te0, sp: nsp0*Ti0, s2p: ns2p0*Ti0, $
              s3p: ns3p0*Ti0, op: nop0*Ti0, o2p: no2p0*Ti0, $
              sp_pu: 0.0, s2p_pu: 0.0, s3p_pu: 0.0, op_pu: 0.0,$
              o2p_pu: 0.0} 

;density at full time step advance
np = {density_p, nel: nel0, ns: ns0, nsp: nsp0, ns2p: ns2p0, ns3p: ns3p0,$
                 no: no0, nop: nop0, no2p: no2p0, fc: fc, fh: fh, $
              nsp_pu: 0.0, ns2p_pu: 0.0, ns3p_pu: 0.0, nop_pu: 0.0, $
              no2p_pu: 0.0, nelh: nelh0, nsph: nsph0, noph: noph0}

npar = replicate({density_p_ar, nel: nel0, nsp: nsp0, ns2p: ns2p0, $
      ns3p: ns3p0, nop: nop0, no2p: no2p0, nelh: nelh0, ns: ns0, no: no0},n_theta)

Tp = {temp_p, Tel: Te0, Telh: Teh0, Tsp: Ti0, Ts2p: Ti0, Ts3p: Ti0, Top: Ti0, $
           To2p: Ti0, Tsp_pu: 0.0, Ts2p_pu: 0.0, Ts3p_pu: 0.0, $
           Top_pu: 0.0, To2p_pu: 0.0, Tpu_s: 0.0, Tpu_o: 0.0, $
           Tel_sp: Te0, Tel_s2p: Te0, Tel_s3p: Te0, Tel_op: Te0, $
           Tel_o2p: Te0, Tel_el: Te0, Tsph: Tih0, Toph: Tih0}
nTp = {energy_p, el: nel0*Te0, sp: nsp0*Ti0, s2p: ns2p0*Ti0, $
              s3p: ns3p0*Ti0, op: nop0*Ti0, o2p: no2p0*Ti0, $
              sp_pu: 0.0, s2p_pu: 0.0, s3p_pu: 0.0, op_pu: 0.0,$
              o2p_pu: 0.0} 

r = {rates, Tel: Te0, is: 0.0, isp: 0.0, is2p: 0.0, rs3p: 0.0, $
            rs2p: 0.0, rsp: 0.0, rop: 0.0, $
            io: 0.0, iop:0.0, ro2p:0.0, $
            Telh: Teh0,ish: 0.0, isph: 0.0, is2ph: 0.0, ioh: 0.0, $
            ioph: 0.0, $
            S_production: 0.0, O_production: 0.0, net_production: 0.0, $
            Transport: trans, $
            o_to_s: otos, cx_k0: 0.0, cx_k1: 0.0, cx_k2: 0.0, cx_k3: 0.0, $
            cx_k4: 0.0,$
            cx_k5: 0.0, cx_k6: 0.0, cx_k7: 0.0, cx_k8: 0.0, cx_k9: 0.0, $
            cx_k10: 0.0, cx_k11: 0.0, cx_k12: 0.0, cx_k13: 0.0, cx_k14: 0.0,$
            cx_k15: 0.0, cx_k16: 0.0, emisSII: fltarr(101,101), $
            emisSIII: fltarr(101,101), $
            emisSIV: fltarr(101,101), emisOII: fltarr(101,101), $
            emisOIII: fltarr(101,101), $
            emistemp: fltarr(101), emisden: fltarr(101), $
            Puv: 0.0, psp: 0.0, ps2p: 0.0, ps3p: 0.0, $
            pop: 0.0, po2p: 0.0, Sei: 0.0, Secx: 0.0, Oei: 0.0, Oecx:0.0, Peuv: 0.0, Peuvh: 0.0}

r1 = r

nl = {neutral_loss, sei: 0.0, seih: 0.0, oei: 0.0, oeih: 0.0, $
       scx: 0.0, ocx: 0.0, k1: 0.0, k2: 0.0, k3: 0.0, k4: 0.0, k5: 0.0, $
       k6: 0.0, k7: 0.0, k8: 0.0, k9: 0.0, k10:0.0, k11: 0.0, k12: 0.0, $
       k14: 0.0, s_tot: 0.0, o_tot: 0.0, scx_tot: 0.0, ocx_tot: 0.0, $
       sion_tot: 0.0, oion_tot: 0.0, fast_O_k5: 0.0, fast_O_k7: 0.0, $
       fast_S_k1: 0.0, fast_S_k3: 0.0, fast_O_k9: 0.0, fast_S_k8: 0.0}

Ec = {energy_conserve, s_ion: 0.0, s_ion_h: 0.0, s_cx: 0.0, o_ion: 0.0, o_ion_h: 0.0, $
       o_cx: 0.0, s_tot_in: 0.0, o_tot_in: 0.0, P_pu: 0.0, eh_eq: 0.0, P_in: 0.0, $
       Puv: 0.0, Pfast: 0.0, Pfast_O: 0.0, Pfast_S: 0.0, Ptrans: 0.0, $ 
       Ptrans_O: 0.0, Ptrans_S: 0.0, Ptrans_e: 0.0, Ptrans_eh: 0.0, $
       P_out: 0.0, in_out: 0.0, ion_e_eq: 0.0}

src = replicate({source_func, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0,$
              o: 0.0, op: 0.0, o2p: 0.0,tm:0.0},ntm+1)
lss = replicate({loss_func, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0,$
              o: 0.0, op: 0.0, o2p: 0.0,tm:0.0},ntm+1)

temps = replicate({temperatures, tsp: 0.0, Ts2p: 0.0, Ts3p: 0.0, Top: 0.0, $
                                 To2p: 0.0, Telec: 0.0},ntm+1)

dens = replicate({dens, nsp: 0.0, ns2p: 0.0, ns3p: 0.0, nop: 0.0, $
                  no2p: 0.0, nel: 0.0},ntm+1)

p = replicate({Power_uv, Puv: 0.0, psp: 0.0, ps2p: 0.0, ps3p: 0.0, $
                         pop: 0.0, po2p: 0.0},ntm+1)


h = {scale_heights, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0, s4p:0.0, $
                              o: 0.0, op: 0.0, o2p: 0.0, el: 0.0}

;h1 = {scale_heights_1, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0, s4p:0.0, $
;                              o: 0.0, op: 0.0, o2p: 0.0, el: 0.0}

;hp = {scale_heights_p, s: 0.0, sp: 0.0, s2p: 0.0, s3p: 0.0, s4p:0.0, $
;                              o: 0.0, op: 0.0, o2p: 0.0, el: 0.0}


;===================================================================================
;initialize variables
;===================================================================================



whlt = where(Lshell lt 6.0)
whge = where(Lshell ge 6.0)
wheq = where(Lshell eq 6.0)
print,'net_source...',net_source_0 
net_source_arr = fextend*net_source_0*(Lshell/6.0)^(-net_source_alpha)
net_source_arr(0) = net_source_arr(0) + flocal*net_source_0

mdot_source = (net_source_arr(0)+total(net_source_arr))*20.*1.67e-27
print,total(net_source_arr)+net_source_arr(0),(net_source_arr(0)+total(net_source_arr))*20.*1.67e-27
;print,'intgl..',net_source_arr(0)*6.0^net_source_alpha*((20.0)^(-net_source_alpha + 1) - $
;      6.0^(-net_source_alpha + 1))/(-net_source_alpha + 1)

nr(whge).nsp = nsp0*(Lshell(whge)/6.0)^(-6)
nr(whge).ns2p = ns2p0*(Lshell(whge)/6.0)^(-6)
nr(whge).ns3p = ns3p0*(Lshell(whge)/6.0)^(-6)
nr(whge).ns4p = ns4p0*(Lshell(whge)/6.0)^(-6)
nr(whge).nop = nop0*(Lshell(whge)/6.0)^(-6)
nr(whge).no2p = no2p0*(Lshell(whge)/6.0)^(-6)


nr.nsph = 0.1
nr.noph = 0.1
whhot1 = where((Lshell ge 6.0) and (Lshell le 7.5))
nr(whhot1).nsph = (0.1/1.5)*(Lshell(whhot1)-6.0)
nr(whhot1).noph = (0.1/1.5)*(Lshell(whhot1)-6.0)

;nr.nsph = 0.0
;nr.noph = 0.0

;nr(whlt).nsp = 10.0
;nr(whlt).ns2p = 10.0
;nr(whlt).ns3p = 10.0
;nr(whlt).nop = 10.0
;nr(whlt).no2p = 10.0

;nr(whge).ns = 10.0/(Lshell(whge)-5.0)^3.0
;nr(whge).no = 50.0/(Lshell(whge)-5.0)^3.0 + 10.0*exp(-(Lshell(whge)-9.0)^2/4.0)

nr(whge).ns = 1e-10
nr(whge).no = 1e-10

;plot_io,Lshell,nr.no
;stop
;nr(wheq).ns = nr(wheq).ns + 5.0
;nr(wheq).no = nr(wheq).no + 20.0 

plot_io,Lshell,nr.no,yrange=[0.1,100],/ysty
oplot,Lshell,nr.ns


nr(*).nel =  nr(*).nsp+2.0*nr(*).ns2p+3.0*nr(*).ns3p+4.0*nr(*).ns4p+nr(*).nop+$
             2.0*nr(*).no2p + nr(*).nsph + nr(*).noph
nr(*).nelh = fh*nel0

tr(whge).Tsp = Ti0 ;+ (1000.0-Ti0)/(Lshell(nintrvl-1) - 6.0)*(Lshell(whge)-6.0)
tr(whge).Ts2p = Ti0 ;+ (1000.0-Ti0)/(Lshell(nintrvl-1) - 6.0)*(Lshell(whge)-6.0)
tr(whge).Ts3p = Ti0 ;+ (1000.0-Ti0)/(Lshell(nintrvl-1) - 6.0)*(Lshell(whge)-6.0)
tr(whge).Ts4p = Ti0 ;+ (1000.0-Ti0)/(Lshell(nintrvl-1) - 6.0)*(Lshell(whge)-6.0)
tr(whge).Top = Ti0 ;+ (1000.0-Ti0)/(Lshell(nintrvl-1) - 6.0)*(Lshell(whge)-6.0)
tr(whge).To2p = Ti0 ;+ (1000.0-Ti0)/(Lshell(nintrvl-1) - 6.0)*(Lshell(whge)-6.0)
tr(whge).Tel = Te0 ;+ (20.0-Te0)/(Lshell(nintrvl-1) - 6.0)*(Lshell(whge)-6.0)
tr(whge).Telh = Teh0_1 ;+ (1e3-Teh0_1)/(Lshell(nintrvl-1) - 6.0)*(Lshell(whge)-6.0)

Tr.Tsph = 0.0 ;10e3
Tr.Toph = 0.0 ;10e3
;whhot1 = where((Lshell ge 6.0) and (Lshell le 7.5))
;Tr(whhot1).Tsph = 1.0+(50e3/1.5)*(Lshell(whhot1)-6.0)
;Tr(whhot1).Toph = 1.0+(50e3/1.5)*(Lshell(whhot1)-6.0)

;tr(whlt).Tsp = 20.0
;tr(whlt).Ts2p = 20.0
;tr(whlt).Ts3p = 20.0
;tr(whlt).Top = 20.0
;tr(whlt).To2p = 20.0
;tr(whlt).Tel = Te0
;tr(whlt).Telh = Teh0

;===================================================================================
;initialize cm3 model stuff
;===================================================================================


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

plot_density = 0
plot_temperature = 0

lonoff =0.0
!p.multi=[0,1,1] 

upsp = ns0
ups2p = ns2p0
ups3p = ns3p0
upop = nop0
upo2p = no2p0

cont='false'
ft_ave = 0

psfile = 'idl.ps'
prt = 1 & lnd = 0
encp = 0 & noencp = 1
clr = 0 & nclr = 1
numclr=!d.n_colors

wh = where(Lshell ge 6.0) 
sr(wh).neutral = 23e27/(1+findgen(n_elements(wh)))^2

;restore,'restart4.sav'

restore,'janfit.sav'
as_r = reverse([6.0,6.2,6.4,6.6,6.8,7.0,7.2,7.4,7.6,7.8,8.15,8.5,8.85])


read_emission_tables,'tepSII',temp,den,emisSII & r.emisSII = emisSII
read_emission_tables,'tepSIII',temp,den,emisSIII & r.emisSIII = emisSIII
read_emission_tables,'tepSIV',temp,den,emisSIV & r.emisSIV = emisSIV
read_emission_tables,'tepOII',temp,den,emisOII & r.emisOII = emisOII
read_emission_tables,'tepOIII',temp,den,emisOIII & r.emisOIII = emisOIII
r.emisden = den
r.emistemp = temp

;nr.ns = nr.ns
;nr.no = nr.no + 10.0*exp(-(Lshell(whge)-9.0)^2/1.0)



;===================================================================================
;Start time and nbox loops
;===================================================================================



for j = 0,nit do begin

    print,'time...',j*runt/60./60./24.

for ii = 0,nbox-1 do begin
    
    i = ii+bufcel

;    print,'box...',i+1

    Rj = 7.14e4*1e5             ;cm
    a = Lshell(i)*Rj
    b = (Lshell(i)+dL(i))*Rj
    Lvol = 0.25*!pi^2*(a+b)*(b-a)^2

    fh = fh_0*(Lshell(i)/6.0)^fh_alpha
    if (fh ge fh_max) then fh = fh_max

    Teh0 = Teh0_1*(Lshell(i)/6.0)^Teh0_alpha
    if (T.Telh ge Teh0_max) then T.Telh = Teh0_max

    n.ns = nr(i).ns 
    n.no = nr(i).no 
    n.nsp= nr(i).nsp 
    n.ns2p = nr(i).ns2p 
    n.ns3p = nr(i).ns3p 
;ns4p0 = nr(i).ns4p 
    n.nop = nr(i).nop 
    n.no2p = nr(i).no2p 
    n.nsph = nr(i).nsph
    n.noph = nr(i).noph
    n.nel = nr(i).nsp+2.0*nr(i).ns2p+3.0*nr(i).ns3p+4.0*nr(i).ns4p+nr(i).nop+2.0*nr(i).no2p
;n.nelh = 0.01*nel0
    n.fh = fh
    n.fc = 1.0 - fh

    n1 = n
    np = n

    T.Tsp = tr(i).Tsp 
    T.Ts2p = tr(i).Ts2p 
    T.Ts3p = tr(i).Ts3p 
;Ts4p0 = tr(i).Ts4p 
    T.Top = tr(i).Top 
    T.To2p = tr(i).To2p 
    T.Tel = tr(i).Tel 
    T.Tsph = tr(i).Tsph 
    T.Toph = tr(i).Toph
    T.Telh = Teh0

    T1 = T
    Tp = T

    nT.el = n.nel*T.Tel
    nT.sp = n.nsp*T.Tsp
    nT.s2p = n.ns2p*T.Ts2p
    nT.s3p = n.ns3p*T.Ts3p
    nT.op = n.nop*T.Top
    nT.o2p = n.no2p*T.To2p

    nT1 = nT
    nTp = nTp

;    ns0 = nr(i).ns 
;    no0 = nr(i).no 
;    nsp0= nr(i).nsp             ;+ sr(i).nsp - lr(i).nsp
;    ns2p0 = nr(i).ns2p          ;+ sr(i).ns2p - lr(i).ns2p
;    ns3p0 = nr(i).ns3p          ;+ sr(i).ns3p - lr(i).ns3p
;;ns4p0 = nr(i).ns4p 
;    nop0 = nr(i).nop            ;+ sr(i).nop - lr(i).nop
;    no2p0 = nr(i).no2p          ;+ sr(i).no2p - lr(i).no2p
;    nsph0 = nr(i).nsph
;    noph0 = nr(i).noph
;    nel0 = nr(i).nsp+2.0*nr(i).ns2p+3.0*nr(i).ns3p+4.0*nr(i).ns4p+nr(i).nop+$
;           2.0*nr(i).no2p+nr(i).nsph+nr(i).noph

;    Tsp0 = tr(i).Tsp 
;    Ts2p0 = tr(i).Ts2p 
;    Ts3p0 = tr(i).Ts3p 
;;Ts4p0 = tr(i).Ts4p 
;    Top0 = tr(i).Top 
;    To2p0 = tr(i).To2p
;    Tsph0 = tr(i).Tsph
;    Toph0 = tr(i).Toph
;    Tel0 = tr(i).Tel 

    trans = 0.0

    net_source = net_source_arr(i)

;    print,'params...',fh,net_source,DLL_0
    cm3_model,n,T,nT,n1,T1,nT1,np,Tp,nTp,r,r1,nar,n1ar,h,h1,hp,$
              nl,Ec,src,lss,temps,dens,p,Lshell(i),dL(i),mr,ii

    tot_peuv(i) = r.peuv 
    tot_peuvh(i) = r.peuvh
;    print,'Peuv..',r.peuv*1.6e-19,r.peuvh*1.6e-19,100.*r.peuvh/r.peuv

    nr(i).ns = n.ns
    nr(i).no = n.no
    nr(i).nsp = n.nsp
    nr(i).ns2p = n.ns2p
    nr(i).ns3p = n.ns3p
;nr(i).ns4p = n.ns4p
    nr(i).nop = n.nop
    nr(i).no2p = n.no2p
    nr(i).nel = nr(i).nsp+2.0*nr(i).ns2p+3.0*nr(i).ns3p+4.0*nr(i).ns4p+ $
                nr(i).nop+2.0*nr(i).no2p

    nelec6 = nr(0).nel + 0.1*nr(0).nel

    tr(i).Tsp = T.Tsp
    tr(i).Ts2p = T.Ts2p
    tr(i).Ts3p = T.Ts3p
;tr(i).Ts4p = T.Ts4p
    tr(i).Top = T.Top
    tr(i).To2p = T.To2p
    tr(i).Tel = T.Tel
    tr(i).Telh = T.Telh

    get_reaction_rates,i,r,n,T,rr,nr(i).nel
    
    mrr(i).sp = mr(0)
    mrr(i).s2p = mr(1)
    mrr(i).s3p = mr(2)
    mrr(i).op = mr(3)
    mrr(i).o2p = mr(4)
    
    
endfor

sti = tr(0).Tsp*nr(0).nsp+tr(0).Ts2p*nr(0).ns2p+tr(0).Ts3p*nr(0).ns3p+tr(0).Top*nr(0).nop+ $
      tr(0).To2p*nr(0).no2p                       
sni = nr(0).nsp+nr(0).ns2p+nr(0).ns3p+nr(0).nop+nr(0).no2p 
tionavg=sti/sni

get_NL2,nr,tr,nl2,nl2e

!p.multi=[0,1,1]

for k = 0,nit_trans-1 do begin
   transport_NL2,nl2,nl2e,DLL_0,DLL_alpha
endfor


;for k = 0,nit_trans-1 do begin
;   transport_flux,nr,tr,DLL_0,DLL_alpha,H,mdot
;endfor



;Schreier
;tr.Tsp = nl2e.nsp/(nl2.nsp*Lshell^(8./3.))
;tr.Ts2p = nl2e.ns2p/(nl2.ns2p*Lshell^(8./3.))
;tr.Ts3p = nl2e.ns3p/(nl2.ns3p*Lshell^(8./3.))
;tr.Top = nl2e.nop/(nl2.nop*Lshell^(8./3.))
;tr.To2p = nl2e.no2p/(nl2.no2p*Lshell^(8./3.))

;Richardson
tr.Tsp = (nl2e.nsp/(nl2.nsp*Lshell^(2)))^(3./4.)
tr.Ts2p = (nl2e.ns2p/(nl2.ns2p*Lshell^(2)))^(3./4.)
tr.Ts3p = (nl2e.ns3p/(nl2.ns3p*Lshell^(2)))^(3./4.)
tr.Top = (nl2e.nop/(nl2.nop*Lshell^(2)))^(3./4.)
tr.To2p = (nl2e.no2p/(nl2.no2p*Lshell^(2)))^(3./4.)

tr(nintrvl-1).Tsp = tr(nintrvl-2).Tsp 
tr(nintrvl-1).Ts2p = tr(nintrvl-2).Ts2p 
tr(nintrvl-1).Ts3p = tr(nintrvl-2).Ts3p 
tr(nintrvl-1).Top = tr(nintrvl-2).Top 
tr(nintrvl-1).To2p = tr(nintrvl-2).To2p 

iterate_NL2_to_equator,nl2,nl2e,nr,tr

;tr.Tsp = 100.
;tr.Ts2p = 100.
;tr.Ts3p = 100.
;tr.Top = 100.
;tr.To2p = 100.


get_NL2,nr,tr,nl2,nl2e

;erase
!p.multi=[0,3,5]

!p.charsize=1.8

wset,0

plot_io,Lshell,nr.ns,title='S',xrange=[6,6+(nbox-1)*dL0],/xsty,yrange=[0.01,1000],/ysty
plot_io,Lshell,nr.no,title='O',xrange=[6,6+(nbox-1)*dL0],/xsty,yrange=[0.01,1000],/ysty
plot_io,Lshell,nr.nsp,title='S+',xrange=[6,6+(nbox-1)*dL0],/xsty,yrange=[0.01,1000],/ysty
oplot,Lshell,nr.nsph,linestyle=1
plot_io,Lshell,nr.ns2p,title='S++',xrange=[6,6+(nbox-1)*dL0],/xsty
plot_io,Lshell,nr.ns3p,title='S+++',xrange=[6,6+(nbox-1)*dL0],/xsty
plot_io,Lshell,nr.nop,title='O+',xrange=[6,6+(nbox-1)*dL0],/xsty,yrange=[0.01,1000],/ysty
oplot,Lshell,nr.noph,linestyle=1
plot_io,Lshell,nr.no2p,title='O++',xrange=[6,6+(nbox-1)*dL0],/xsty

peuv_peuvh = (tot_peuv(0) + tot_peuvh(0)+total(tot_peuv+tot_peuvh))*1.6e-19
;peuv_peuvh = total(tot_peuv+tot_peuvh)*1.6e-19

plot,Lshell,tr.Tsp,xrange=[6,10],title='S+',/xsty
plot,Lshell,tr.Ts2p,xrange=[6,6+(nbox-1)*dL0],title='S++',/xsty
plot,Lshell,tr.Ts3p,xrange=[6,6+(nbox-1)*dL0],title='S+++',/xsty
plot,Lshell,tr.Top,xrange=[6,6+(nbox-1)*dL0],title='O+',/xsty
plot,Lshell,tr.To2p,xrange=[6,6+(nbox-1)*dL0],title='O++',/xsty
plot,Lshell,tr.Tel,xrange=[6,6+(nbox-1)*dL0],title='elec',/xsty
plot,Lshell,nr.nel,xrange=[6,6+(nbox-1)*dL0],/xsty
plot_io,Lshell,tot_peuv*1.6e-19,ytitle='Power (W)',$
     xrange=[6,6+(nbox-1)*dL0],/xsty,yrange=[1e8,1e12],/ysty,$
     title='Total Peuv...'+string(peuv_peuvh)+' (W)'
oplot,Lshell,tot_peuvh*1.6e-19,linestyle=1


!p.multi=[0,2,3]

wset,1
plot_io,s.l,s.tot33,yrange=[1e34,8e36],/ysty,xrange=[5,12],/xsty
oplot,Lshell,nl2.nsp+2.0*nl2.ns2p+3.0*nl2.ns3p+nl2.nop+2.0*nl2.no2p,$
   linestyle=1
plot_io,s.l,s.s1,title='S+',yrange=[1e34,4e36],/ysty,xrange=[5,12],/xsty
oplot,Lshell,nl2.nsp,linestyle=1
plot_io,s.l,s.s2,title='S++',yrange=[1e34,4e36],/ysty,xrange=[5,12],/xsty
oplot,Lshell,nl2.ns2p,linestyle=1
plot_io,s.l,s.s3,title='S+++',yrange=[1e34,4e36],/ysty,xrange=[5,12],/xsty
oplot,Lshell,nl2.ns3p,linestyle=1
plot_io,s.l,s.o1,title='O+',yrange=[1e34,4e36],/ysty,xrange=[5,12],/xsty
oplot,Lshell,nl2.nop,linestyle=1
plot_io,s.l,s.o2,title='O++',yrange=[1e34,4e36],/ysty,xrange=[5,12],/xsty
oplot,Lshell,nl2.no2p,linestyle=1

!p.multi=[0,2,3]
wset,2
nelc = nr.nel + 0.1*nr.nel

plot_io,Lshell,mrr.sp,linestyle=5,yrange=[0.01,0.4],/ysty,$
   xrange=[6,Lshell(nbox-1)],/xsty,title='S+
oploterr,as_r,s2mix,s2mixerr
plot_io,Lshell,mrr.op,yrange=[0.01,0.4],/ysty,xrange=[6,Lshell(nbox-1)],/xsty,$
  linestyle=1,title='O+'
oploterr,as_r,o2mix,o2mixerr
plot_io,Lshell,mrr.s2p,linestyle=2,yrange=[0.01,0.4],/ysty,$
  xrange=[6,Lshell(nbox-1)],/xsty,title='S++'
oploterr,as_r,s3mix,s3mixerr

plot_io,Lshell,mrr.o2p,linestyle=4,yrange=[0.01,0.4],/ysty,$
   xrange=[6,Lshell(nbox-1)],/xsty,title='O++'
oploterr,as_r,o3mix,o3mixerr

plot_io,Lshell,mrr.s3p,linestyle=3,yrange=[0.01,0.4],/ysty,$
   xrange=[6,Lshell(nbox-1)],/xsty,title='S+++'
oploterr,as_r,s4mix,s4mixerr

;plot_io,Lshell,mdot,title='Mdot'
;oplot,[!x.crange(0),!x.crange(1)],[mdot_source,mdot_source],linestyle=1


mr_r(0) = mrr(0).sp
mr_r(1) = mrr(12).sp
mr_r(2) = mrr(0).s2p
mr_r(3) = mrr(12).s2p
mr_r(4) = mrr(0).s3p
mr_r(5) = mrr(12).s3p
mr_r(6) = mrr(0).op
mr_r(7) = mrr(12).op
mr_r(8) = mrr(0).o2p
mr_r(9) = mrr(12).o2p
mr_r(10) = mrr(4).sp
mr_r(11) = mrr(4).s3p


save,filename='restart.sav',nr,tr,nl2,nl2e,mrr,s,tot_peuv,tot_peuvh,rr

endfor

end
;-------------------------------------------------------------------

