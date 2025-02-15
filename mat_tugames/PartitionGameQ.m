function PGQ=PartitionGameQ(th,w_vec,hrQ,tol)
% PartitionGameQ verifies if (th,w_vec) induces a partition game.
%
% Source: Isbell (1956) A Class Of Majority Games.   
%         Sudhoelter (1996), Star-shapedness of the kernel for homogeneous games.
%
% Usage: PGQ=PartitionGameQ(th,w_vec);
%
% Define field variables:
%  output:
%  pQ          -- Returns one (true) whenever the simple game is a partition game.
%  MR.mrQ      -- Returns one (true) whenever a minimal homogeneous representation was found, 
%                 otherwise zero (false).
%  MR.hrQ      -- Returns one whenever the representation is homogeneous.
%  MR.mwgs     -- Vector of minimal weights.
%  MR.th       -- Quorum of the minimal homogeneous representation.
%  MR.eQ       -- Returns one whenever the both weighed majority games 
%                 i.e., (th0,w_vec0) vs. (th,mwgs), are equal.
%  MR.smpl     -- Players with character sum.
%  MR.stpl     -- Players with character step.
%  MR.npl      -- Players with character null-player.
%  MR.th0      -- Original quorum.
%  MR.wvec0    -- Original weights.
%
%  input:
%  th       -- Threshold/quorum to pass a bill (positive number).
%  w_vec    -- Vector of weights (descend ordering).
%  hrQ      -- Set one (true) for a homogeneous representation, 
%              otherwise empty set '', whereas a zero interrupts 
%              the evaluation Default is one, and suppresses the 
%              costly check of homogeneity.
%  tol      -- Numerical tolerance.
%


%    
%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   05/17/2022        1.9.1           hme    
%    
%
if nargin<3
   hrQ=true; 
   tol=10^8*eps;
end   
n=length(w_vec);
MR=min_homogrep(th,w_vec,hrQ,tol);
s1=size(MR.M);
if n>s1
   pQ=false;
elseif n==s1
   pQ=CheckWeightsQ(MR.mwgs,n);
else
   pQ=false;  
end    
    
PGQ.pQ=pQ;
PGQ.MR=MR;


%%%%%%%%%%%%%%%%%%%%%%%%%%
function cwQ=CheckWeightsQ(wgs,n)

acw=false(1,n);    
acw(n-1:n)=wgs(n-1)==wgs(n)==1;
acw(2:3)=wgs(2)==wgs(3);
acw(1)=wgs(1)==1+sum(wgs(4:n));

n2=n-2;
for kk=n2:-1:3
    c4a=wgs(kk)>=wgs(kk+1);
    c4b=wgs(kk)<=1+sum(wgs(kk+2:n));
    acw(kk)= c4a & c4b;
end 
cwQ=all(acw);    
    
