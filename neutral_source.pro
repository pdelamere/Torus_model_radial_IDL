
alpha = 30
Lo = 6.0
lmax = 9.0
ndot0 = 1e28
x0 = 1.0
x1 = 9.0/Lo

ndot_tot = 2e28

aintgl = Lo*(x1^(-alpha + 1) - x0^(-alpha + 1))/(-alpha + 1)

ndot0 = ndot_tot/aintgl

;print,ndot0
;tot = ndot0*Lo*(x1^(-alpha + 1) - x0^(-alpha + 1))/(-alpha + 1)


L = Lo + findgen(10)*0.25
dL = L(1) - L(0)
wh = where(L le 9.0)

nintgl = total((L(wh)/Lo)^(-alpha))

print,ndot_tot/nintgl,ndot_tot*20.*1.67e-27


end