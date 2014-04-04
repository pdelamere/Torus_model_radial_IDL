;-------------------------------------------------------------------
pro get_reaction_rates,i,r,n,T,rr,nel
;-------------------------------------------------------------------

rr(i).is = r.is*nel
rr(i).isp = r.isp*nel
rr(i).is2p = r.is2p*nel
rr(i).io = r.io*nel
rr(i).iop = r.iop*nel

rr(i).ish = r.ish*n.fh*nel
rr(i).isph = r.isph*n.fh*nel
rr(i).is2ph = r.is2ph*n.fh*nel
rr(i).ioh = r.ioh*n.fh*nel
rr(i).ioph = r.ioph*n.fh*nel

rr(i).k1_s = r.cx_k1*n.nsp
rr(i).k1_sp = r.cx_k1*n.ns
rr(i).k2_s = r.cx_k2*n.ns2p
rr(i).k2_s2p = r.cx_k2*n.ns
rr(i).k3_s = r.cx_k3*n.ns2p
rr(i).k3_s2p = r.cx_k3*n.ns
rr(i).k4_s = r.cx_k4*n.ns3p
rr(i).k4_s3p = r.cx_k4*n.ns
rr(i).k5_o = r.cx_k5*n.nop
rr(i).k5_op = r.cx_k5*n.no
rr(i).k6_o = r.cx_k6*n.no2p
rr(i).k6_o2p = r.cx_k6*n.no
rr(i).k7_o = r.cx_k7*n.no2p
rr(i).k7_o2p = r.cx_k7*n.no
rr(i).k8_o = r.cx_k8*n.nsp
rr(i).k8_sp = r.cx_k8*n.no
rr(i).k9_s = r.cx_k9*n.nop
rr(i).k9_op = r.cx_k9*n.ns
rr(i).k10_s = r.cx_k10*n.no2p
rr(i).k10_o2p = r.cx_k10*n.ns
rr(i).k11_s = r.cx_k11*n.no2p
rr(i).k11_o2p = r.cx_k11*n.ns
rr(i).k12_o = r.cx_k12*n.ns2p
rr(i).k12_s2p = r.cx_k12*n.no
rr(i).k13_o2p = r.cx_k13*n.nsp
rr(i).k13_sp = r.cx_k13*n.no2p
rr(i).k14_o = r.cx_k14*n.ns3p
rr(i).k14_s3p = r.cx_k14*n.no
rr(i).k15_o2p = r.cx_k15*n.ns2p
rr(i).k15_s2p = r.cx_k15*n.no2p
rr(i).k16_s3p = r.cx_k16*n.nsp
rr(i).k16_sp = r.cx_k16*n.ns3p


return
end
;-------------------------------------------------------------------

