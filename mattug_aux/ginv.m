function ginv = ginv(X)
if isempty(X)
% quick return
ginv = zeros(size(X'),class(X));
return
end
[n,m]=size(X);
if n > m
   C = X'*X ;
   ginv = C\X';
else
   C = X*X';
   G = C\X;
   ginv = G';
end
end
