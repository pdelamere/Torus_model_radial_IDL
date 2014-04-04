set_plot,'ps
device,/encapsulated
device,/inches,xsize=6.5,ysize=6.0
device,filename='compo.eps
!p.font=0
!p.thick=2.0


restore,'dens.sav

oct = [0.0994,0.216,0.0218,0.263,0.0166]
jan = [0.0603,0.207,0.0372,0.269,0.0169]
nov = [0.047,0.195,0.0445,0.284,0.0179]

nelec = dens.nel + 0.1*dens.nel

dt = 10000.0 
runt = 105*8.64e4
ntm = fix(runt/dt)


mint = 275
tm = mint+findgen(ntm)*dt/8.64e4

plot_io,tm,dens.nop/nelec,xrange=[mint-5,max(tm)+5],/xsty,$
   yrange=[0.01,0.35],/ysty,xtitle='t (days)',ytitle='n!di!n /n!de!n',$
   charsize=1.2
oplot,tm,dens.nsp/nelec,linestyle=1
oplot,tm,dens.ns2p/nelec,linestyle=2
oplot,tm,dens.ns3p/nelec,linestyle=4
oplot,tm,dens.no2p/nelec,linestyle=5

plots,mint,oct(0),psym=1,symsize=1.6
plots,mint,oct(1),psym=2,symsize=1.6
plots,mint,oct(2),psym=4,symsize=1.6
plots,mint,oct(3),psym=5,symsize=1.6
plots,mint,oct(4),psym=6,symsize=1.6

plots,mint+45,nov(0),psym=1,symsize=1.6
plots,mint+45,nov(1),psym=2,symsize=1.6
plots,mint+45,nov(2),psym=4,symsize=1.6
plots,mint+45,nov(3),psym=5,symsize=1.6
plots,mint+45,nov(4),psym=6,symsize=1.6

plots,mint+105,jan(0),psym=1,symsize=1.6
plots,mint+105,jan(1),psym=2,symsize=1.6
plots,mint+105,jan(2),psym=4,symsize=1.6
plots,mint+105,jan(3),psym=5,symsize=1.6
plots,mint+105,jan(4),psym=6,symsize=1.6


xyouts,0.23,0.89,'O!u+!n',/normal,charsize=1.5
xyouts,0.23,0.8,'S!u++!n',/normal,charsize=1.5
xyouts,0.23,0.64,'S!u+!n',/normal,charsize=1.5
xyouts,0.23,0.38,'S!u+++!n',/normal,charsize=1.5
xyouts,0.23,0.18,'O!u++!n',/normal,charsize=1.5

device,/close


end
