set_plot,'ps
device,/encapsulated
device,filename='power.eps
device,/inches,xsize=6.5,ysize=6.0,xoffset=1.0,yoffset=1.0
!p.font=0
!p.thick=1.0
!p.charsize = 1.2



dt = 10000.0 
runt = 60*8.64e4
ntm = fix(runt/dt)

mint = 275
restore,'torus_power.dat
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

      vol = 3.1e31
      plot,tm,p.lps2p*vol*1.6e-19/1e12,xtitle='time (days) ',$
              ytitle='Power Radiated (Terawatts)',$
              xrange=mint+[0,dt*ntm]/8.64e4,/xsty,thick=3,$
              yrange=[0,0.05+max(p.lps2p*vol*1.6e-19/1e12)],/ysty
      oplot,congrid(time,50),congrid(smooth(siii,10),50)/1e12,psym=1,$
         symsize=0.8
      oplot,tm,p.lpsp*vol*1.6e-19/1e12,linestyle=4,thick=3 
      oplot,congrid(time,50),congrid(smooth(sii,10),50)/1e12,psym=6,$
        symsize=0.8
      oplot,tm,p.lps3p*vol*1.6e-19/1e12,linestyle=2,thick=3
      oplot,congrid(time,50),congrid(smooth(siv,10),50)/1e12,psym=4,$
        symsize=0.8
      oplot,tm,p.lpop*vol*1.6e-19/1e12,linestyle=5,thick=3
      oplot,congrid(time,50),congrid(smooth(oii,10),50)/1e12,psym=5,$
          symsize=0.8
;      oplot,tm,p.po2p,linestyle=5 
;      if (lonoff eq 1.0) then begin
      legend,['S!u+!n','S!u++!n','S!u+++!n','O!u+!n'],$
       linestyle=[4,0,2,5],psym=-[6,1,4,5],/top,/right
;      endif

device,/close

end
