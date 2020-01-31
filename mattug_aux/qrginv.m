function A = qrginv(B) 
[N,M] = size(B);
if issparse(B) 
  [Q,R,P] = spqr(B);
else
  [Q,R,P] = qr(B);
end
r=sum(any(abs(R)>1e-7,2)); 
R1 = R(1:r,:); 
R2 = ginv(R1); 
R3 = [R2 zeros(M,N-r)]; 
A = P*R3*Q'; 
%qrginv = A;
