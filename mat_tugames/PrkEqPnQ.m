function PKSVQ=PrkEqPnQ(v,x,tol)
% PRKEQPNQ checks whether the pre-kernel of a game v is a singleton.
% 
%  Usage: PKSVQ=PrkEqPnQ(v)
%
%  Source: Meinhardt 2022, On the Replication of the Pre-Kernel and Related Solutions, Thm 3.4.
%
%
% Define variables:
%  output: Fields
%  pksvQ    -- Returns one (true) whenever the pre-kernel is a singleton, otherwise zero (false).
%              Composte by xzeqQ and rdQ.
%  x        -- Returns the original pre-kernel element.
%  z        -- Returns the computed pre-kernel element, must be a replication of x, otherwise the 
%              the pre-kernel element cannot be a singleton.
%  xzeqQ    -- Returns true if the payoffs x and z are equal, otherwise false.
%  pkQ      -- Returns true if the x is replicated as a pre-kernel element by z, otherwise false. 
%  dmrd     -- Returns the list of Davis-Maschler reduced games w.r.t. the partition set and x.
%  dmrdcv   -- Returns the list of Davis-Maschler reduced games w.r.t. the partition set and x that
%              must be given by v_Sx(T)=max{v(T),v(T \cup S^c)-x(S^c)} for T subset in S.
%  rdQ      -- Returns true if both types of reduced games are equal for all R in the partition. 
%  ptnQ     -- Returns true if D(x) contains a partition.
%  ptn      -- Returns the set of partitions (or anti-partition) which are included in D(x).
%  setD     -- Returns the collection of coalitions with largest nontrivial excess, i.e., D(x).
%
%  input:
%  v        -- A Tu-Game v of length 2^n-1. 
%  x        -- Payoff vector of size(1,n), pre-kernel vector.
%  tol      -- Tolerance value. By default, it is set to 10^6*eps.
%

%  Author:        Holger I. Meinhardt (hme)
%  E-Mail:        Holger.Meinhardt@wiwi.uni-karlsruhe.de
%  Institution:   University of Karlsruhe (KIT)  
%
%  Record of revisions:
%   Date              Version         Programmer
%   ====================================================
%   01/28/2022        1.9.1           hme
%

if nargin<2
  tol=10^6*eps;
  x=PreKernel(v); %% Taking an arbitrary pre-kernel element.
elseif nargin<3
  tol=10^6*eps;	
end	

N=length(v);
[~, n]=log2(N);
pl=1:n;

SATIS=satisfaction(v,x);
z=zeros(1,n);
if SATIS.ptnQ==0
   PKSVQ.pksvQ=false;
   PKSVQ.x=x;
   PKSVQ.z=z;
   PKSVQ.xzeQ=false;
   PKSVQ.pkQ=false;
   PKSVQ.dmrd='none';
   PKSVQ.dmrdcv='none';
   PKSVQ.rdQ=false;
   PKSVQ.ptnQ=false;
   PKSVQ.ptn=SATIS.ptn;
   PKSVQ.setD=SATIS.setD;
   return;
else	
  ptn=SATIS.ptn';
  aptn=SATIS.aptn;
  [sa1,sa2]=size(aptn);
  if sa1>1 & sa2>=1; 
     aptn=SATIS.aptn';
  end
  setD=SATIS.setD;
  ptn=sort(ptn(1,:),'ascend'); % Taking the first partition.
  aptn=sort(aptn(1,:),'ascend'); % Taking the first anti-partition.
  if nnz(ismember(setD,ptn))<=1;
     ptn=aptn;   % Then setD must contain an anti-partition instead of a partition.
  end	  
  sz=size(ptn);
  v_t=cell(1,sz(2));
  v_tcv=cell(1,sz(2));
  rdeQ=false(1,sz(2));
  for k=1:sz(2)
      sS=ptn(1,k);
      inS=bitget(sS,pl)==1;
      v_t{1,k}=RedGame(v,x,sS);
      v_tcv{1,k}=CvRedGame(v,x,sS);
      rdeQ(k)=all(abs(v_t{1,k}-v_tcv{1,k})<tol);
      if rdeQ(k)==1
         pk_vx=PreKernel(v_t{1,k});
      else 
         pk_vx=PreKernel(v_tcv{1,k});	      
      end	   
      z(inS)=pk_vx;   
  end
  pkQ=PrekernelQ(v,z,tol);
  xzeqQ=all(abs(x-z)<tol);
  rdQ=all(rdeQ);
end
PKSVQ.pksvQ=rdQ & xzeqQ;
PKSVQ.x=x;
PKSVQ.z=z;
PKSVQ.pkQ=pkQ;
PKSVQ.xzeqQ=xzeqQ;
PKSVQ.dmrd=v_t;
PKSVQ.dmrdcv=v_tcv;
PKSVQ.rdQ=all(rdeQ);
PKSVQ.ptnQ=SATIS.ptnQ;
PKSVQ.ptn=ptn;
PKSVQ.setD=setD;


