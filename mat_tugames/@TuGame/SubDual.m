function sdv=SubDual(clv,sS)
% Usage: sdv=SubDual(clv,sS)

v=clv.tuvalues;
bd=length(sS)
if bd>2
for k=1:bd-1;
 CN=sS(bd)-sS(k);
 cv=v(CN);
 cv(bd)=0;
 sdv=(v(bd)-cv);
end
 else
   cv=0;
   sdv=v-cv;
end

