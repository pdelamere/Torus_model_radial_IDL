set_plot,'ps
device,/encapsulated
device,filename='power_euv.eps
device,/inches,xsize=6.5,ysize=6.0,xoffset=1.0,yoffset=1.0
!p.font=0
!p.thick=1.0
!p.charsize = 1.2



dt = 10000.0 
runt = 60*8.64e4
ntm = fix(runt/dt)

mint = 275
;restore,'torus_power.dat
      restore,'pwr.sav'
      tm = mint+findgen(ntm)*dt/8.64e4
;      plot,tm,p.ps2p,xtitle='time (days) ',$
;           ytitle='Power Radiated (eV cm!u-3!n s!u-1!n)',$
;           xrange=[0,(dt*ntm + dt)]/8.64e4,/xsty,linestyle=3,$
;           yrange=[0,max(p.ps2p)*1.1],/ysty
;      oplot,tm,p.psp,linestyle=1 
;;      oplot,tm,p.ps2p,linestyle=2 
;      oplot,tm,p.ps3p,linestyle=2 
;      oplot,tm,p.pop,linestyle=4 
;;      oplot,tm,p.po2p,linestyle=5 
;      legend,['S!u+!n','S!u++!n','S!u+++!n','O!u+!n'],$
;       linestyle=[1,3,2,4],/top,/right

      vol = 2.6e31
      print,p.peuv
      plot,tm,p.peuv*vol*1.6e-19/1e12,xtitle='time (days) ',$
              ytitle='Power Radiated (Terawatts)',$
              xrange=mint+[0,dt*ntm]/8.64e4,/xsty,thick=3,$
              yrange=[0,0.05+max(p.lps2p*vol*1.6e-19/1e12)],/ysty
;      legend,['S!u+!n','S!u++!n','S!u+++!n','O!u+!n'],$
;       linestyle=[5,0,2,4],psym=[6,1,4,5],/top,/right
;      endif

device,/close

end
