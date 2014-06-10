L = 6+findgen(100)/10
mperm = 11*(L/6.0)^(-5)
mdot = 500.
dL = L(1)-L(0)

plot,L,mperm,/ylog,/xlog,xrange=[6,15],/xsty

v = (mdot/mperm)/7.14e7

plot,L,v

plot,[6,15],[0,100],/nodata
for i = 0,n_elements(v)-1 do begin
   plots,L(i),total(dL/v(0:i))/60./60./24.,psym=1
   print,total(dL/v(0:i))/60./60./24.
endfor

f = 60e-7
alpha = 6.0

v = dL*(f*(L/6.0)^(alpha))

for i = 0,n_elements(v)-1 do begin
   plots,L(i),total(dL/v(0:i))/60./60./24.,psym=5
   print,total(dL/v(0:i))/60./60./24.
endfor



end



