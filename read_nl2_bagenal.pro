pro read_nl2,s


ss = {input,i:0,l: 0.0, tot33: 0.0, s1: 0.0, s2: 0.0, s3: 0.0, $
               o1: 0.0, o2: 0.0, hot: 0.0, oth: 0.0}
s = replicate({nl2,i:0,l: 0.0, tot33: 0.0, s1: 0.0, s2: 0.0, s3: 0.0, $
               o1: 0.0, o2: 0.0, hot: 0.0, oth: 0.0},46)


close,1
openr,1,'nl2_bagenal.txt'

d1=0
d2=0.
d3=0.
d4=0.

for i=0,45 do begin
  ;readf,1,d1,d2,d3,d4,d5,d6,d7,d8,d9
  readf,1,ss
  ;readf,1,d1,d2,d3,d4
  ;print,d1,d2,d3,d4
;  print,ss
  s(i).i = ss.i
  s(i).l = ss.l
  s(i).tot33 = ss.tot33
  s(i).s1 = ss.s1
  s(i).s2 = ss.s2
  s(i).s3 = ss.s3
  s(i).o1 = ss.o1
  s(i).o2 = ss.o2
  s(i).hot = ss.hot
  s(i).oth = ss.oth
endfor

close,1

end
