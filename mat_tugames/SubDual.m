function sdv=SubDual(v,sS)
% Usage: sdv=SubDual(v,sS)

bd=length(sS)
if bd>2
for k=1:bd-1;
% N=sS(end);
 CN=sS(bd)-sS(k);
 cv=v(CN);
 cv(bd)=0;
 sdv=(v(bd)-cv);
end
 else
   cv=0;
   sdv=v-cv;
end

