L = 12.0
H = 1.0

t = H/L

intgl = (sin(t)*cos(t)^6)/7 + (6./7.)*( (sin(t)*cos(t)^4)/5 + (4./5.)*(sin(t) - sin(t)^3/3))

print,intgl


theta = (findgen(10))*!dtor
dtheta = theta(1)-theta(0)

print,total(cos(theta)^7)*dtheta

end

