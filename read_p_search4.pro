close,1
openr,1,'params4.dat
junk=' '

fh = 0.0
fha = 0.0
ns = 0.0
nsa =0.0
dll = 0.0
dlla = 0.0
teh = 0.0
teha = 0.0
peuv = 0.0

s = {params,fh: 0.0, fha: 0.0, ns: 0.0, nsa:0.0, dll:0.0, dlla: 0.0, $
            mr: fltarr(12), nelec: 0.0, peuv:0.0, chi: 0.0}

readf,1,s
p = s
print,s


while (not(eof(1))) do begin

   readf,1,s
   print,s
   p = [p,s]

endwhile

!p.multi=[0,2,5]
nrun = indgen(n_elements(p.chi))
!p.charsize=1.4

plot,nrun,p.fh,title='fh',yrange=[0.002,0.003]
plot,nrun,p.fha,title='fh_alpha',yrange=[3,7]
plot,nrun,p.ns,title='Sn
plot,nrun,p.nsa,title='Sn_alpha
plot,nrun,p.dll,title='DLL
plot,nrun,p.dlla,title='DLL_alpha
;plot,nrun,p.teh,title='Teh
;plot,nrun,p.teha,title='Teh_alpha
plot,nrun,p.peuv,title='PEUV'
plots,[!x.crange(0),!x.crange(1)],[1.2e12,1.2e12],linestyle=1
plots,[!x.crange(0),!x.crange(1)],[1.7e12,1.7e12],linestyle=1
plot,nrun,p.chi,title='CHI^2


;endwhile

wh = where((p.dlla ge 3.0) and (p.dlla le 4.0))

print,wh

stop
end
