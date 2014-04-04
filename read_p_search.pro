set_plot,'ps
@x6x9


close,1
openr,1,'p_search.dat
junk=' '

fh = 0.0
fha = 0.0
ns = 0.0
nsa =0.0
dll = 0.0
dlla = 0.0
peuv = 0.0

readf,1,junk
readf,1,junk
readf,1,junk
readf,1,junk

readf,1,format="(a16,f10.5)",junk,a
fh = a
readf,1,format="(a16,f10.5)",junk,a
fha = a
readf,1,format="(a16,f20.5)",junk,a
ns = a
readf,1,format="(a16,f10.5)",junk,a
nsa = a
readf,1,format="(a16,f20.5)",junk,a
dll = a
readf,1,format="(a16,f10.5)",junk,a
dlla = a

readf,1,junk
readf,1,junk
readf,1,format="(a15,f20.5)",junk,a
peuv = a
readf,1,junk
readf,1,format="(a15,f10.5)",junk,a
chi = a


readf,1,junk



print,fh,fha,ns,nsa,dll,dlla,peuv,chi


while (not(eof(1))) do begin
readf,1,junk
readf,1,junk
readf,1,junk
readf,1,junk

readf,1,format="(a16,f10.5)",junk,a
fh = [fh,a]
readf,1,format="(a16,f10.5)",junk,a
fha = [fha,a]
readf,1,format="(a16,f20.5)",junk,a
ns = [ns,a]
readf,1,format="(a16,f10.5)",junk,a
nsa = [nsa,a]
readf,1,format="(a16,f20.5)",junk,a
dll = [dll,a]
readf,1,format="(a16,f10.5)",junk,a
dlla = [dlla,a]

readf,1,junk
readf,1,junk
readf,1,format="(a15,f20.5)",junk,a
peuv = [peuv,a]
readf,1,junk
readf,1,format="(a15,f10.5)",junk,a
chi = [chi,a]
readf,1,junk

print,chi

!p.multi=[0,2,4]
nrun = indgen(n_elements(chi))
!p.charsize=1.4

plot,nrun,fh,title='fh',yrange=[0.002,0.003]
plot,nrun,fha,title='fh_alpha',yrange=[3,7]
plot,nrun,ns,title='Sn
plot,nrun,nsa,title='Sn_alpha
plot,nrun,dll,title='DLL
plot,nrun,dlla,title='DLL_alpha
plot,nrun,peuv,title='PEUV'
plots,[!x.crange(0),!x.crange(1)],[1.2e12,1.2e12],linestyle=1
plots,[!x.crange(0),!x.crange(1)],[1.7e12,1.7e12],linestyle=1
plot,nrun,chi,title='CHI^2


endwhile







stop
end
